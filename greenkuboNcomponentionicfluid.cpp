/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
  *
**/

#include "greenkuboNcomponentionicfluid.h"
#include "chargefluxts.h"
#include "heatfluxts.h"
#include <string>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include "cronometro.h"
#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif

GreenKuboNComponentIonicFluid::GreenKuboNComponentIonicFluid(ReadLog *traiettoria,
                                                             std::string log,
                                                             unsigned int skip,
                                                             std::vector<std::string> headers,
                                                             bool dump,
                                                             unsigned int lunghezza_funzione_max,
                                                             unsigned int nthreads,
                                                             bool subtract_mean,
                                                             unsigned int start_mean,
                                                             unsigned int n_seg,
                                                             bool do_bench,
                                                             unsigned int n_seg_start,
                                                             unsigned int n_seg_stop) : OperazioniSuLista<GreenKuboNComponentIonicFluid>(),
    traiettoria (traiettoria), log(log), ntimesteps(0),skip(skip), scrivi_file(dump),
    lmax(lunghezza_funzione_max),nthread(nthreads),subtract_mean(subtract_mean),
    start_mean(start_mean),n_seg(n_seg),bench(false),
    n_seg_start(n_seg_start), n_seg_stop(n_seg_stop)
{


    if (n_seg<1){
        std::cerr << "Attenzione: n_seg < 1 . Imposto a 1.\n";
        n_seg=1;
    }

    if(!do_bench)
        benchmarked=true;

    std::pair<unsigned int ,bool> res;


    //Trova gli indici degli header delle correnti da utilizzare per il calcolo
    for (unsigned int j=0;j<headers.size();j++){
        res=traiettoria->get_index_of(headers.at(j));
        if (res.second)
            idx_j.push_back(res.first);
        else {
            std::cerr << "Errore: header '" << headers.at(j) << "' non trovato!\n";
            abort();
        }
    }

#ifdef HALF_CORR
    N_corr=idx_j.size()*(idx_j.size()+1)/2;
#else
    N_corr=idx_j.size()*idx_j.size();
#endif
    narr=3*N_corr+2;
}

GreenKuboNComponentIonicFluid::~GreenKuboNComponentIonicFluid(){
#ifdef DEBUG2
    //    std::cerr << "Called delete GreenKuboNComponentIonicFluid, je="<<je<<", jz="<<jz<<"\n";
#endif
}

GreenKuboNComponentIonicFluid & GreenKuboNComponentIonicFluid::operator =(const GreenKuboNComponentIonicFluid & destra) {
#ifdef DEBUG2
    std::cerr << "Chiamato GreenKuboNComponentIonicFluid::operator =\n";
#endif
    OperazioniSuLista<GreenKuboNComponentIonicFluid>::operator =( destra);
    return *this;
}

unsigned int GreenKuboNComponentIonicFluid::numeroTimestepsOltreFineBlocco(unsigned int n_b) {
    return (traiettoria->n_timestep()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->n_timestep()/(n_b+1)+1 : lmax;
}

void GreenKuboNComponentIonicFluid::reset(unsigned int numeroTimestepsPerBlocco) {
    leff=(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    lunghezza_lista=(leff)*narr;
    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new double [lunghezza_lista];
    for (unsigned int i=0;i<lunghezza_lista;i++){
        lista[i]=0.0;
    }
}

double * GreenKuboNComponentIonicFluid::jN(unsigned int N,unsigned int ts){
    return&traiettoria->line(ts)[idx_j[N]];
}



unsigned int GreenKuboNComponentIonicFluid::get_narr(){
    return narr;
}

unsigned int GreenKuboNComponentIonicFluid::get_indexOfKappa(){
    return 3*N_corr;
}

void GreenKuboNComponentIonicFluid::calcola(unsigned int primo) {


    if(!benchmarked)
        n_seg_bench();


    if(nthread<1)
        nthread=1;

    //nthread > 1

    /*dividi il lavoro in gruppi
         itimestep  --->  [0,leff[
         npassi     --->  leff

         ciascun gruppo avrà npassith=npassi/nthread passi
         l'ultimo deve finire alla fine

         gli integrali li calcolo alla fine, è uguale.
        Nel caso di 1 thread tengo l'implementazione precedente
            (anche per poter confrontare i risultati)
        */


    double *matr=new double [idx_j.size()*idx_j.size()];
    double *intJJ=new double[N_corr];
    double *int_ein_JJ=new double[N_corr];
    double *JJm_T=new double[N_corr];
    unsigned int cont_JJ=0;

    for (unsigned int j=0;j<N_corr;j++){
        intJJ[j]=0.0;
        int_ein_JJ[j]=0.0;
    }
    cont_JJ=0;

    unsigned int npassith=leff/nthread;
    std::vector<std::thread> threads;
    unsigned int timesteps_seg=ntimesteps/n_seg;
    for (unsigned int iseg=0;iseg<n_seg;iseg++){
        ++cont_JJ;
        for (unsigned int ith=0;ith<nthread;ith++){
            threads.push_back(std::thread([&,ith](){
                double *JJ=new double[N_corr];
//                unsigned int ultimo= (ith != nthread-1 )?npassith*(ith+1):leff;
//                for (unsigned int itimestep=npassith*ith;itimestep<ultimo;itimestep++) {
                  //fa fare ai threads diversi conti con i dati vicini (per ottimizzare l'uso della cache in comune fra i vari core)
                    for (unsigned int itimestep=ith;itimestep<leff;itimestep+=nthread) {
                    //fa la media sulla traiettoria dei vari prodotti,
                    //con una differenza di timesteps fissata "itimestep"

                    for (unsigned int j=0;j<N_corr;j++){
                        JJ[j]=0.0;
                    }


                    unsigned int cont=0;
                    // cambio questo: la media viene fatta a pezzettini

                    unsigned int allinea=(timesteps_seg*iseg)+(timesteps_seg*iseg)%skip;
                    for (unsigned int jmedia=allinea;jmedia<timesteps_seg*(iseg+1);jmedia+=skip) {
                        //prodotto JzJz
                        cont++;
                        unsigned int idxj=0;
                        for (unsigned int j1=0;j1<idx_j.size();j1++)
#ifdef HALF_CORR
                            for (unsigned int j2=j1;j2<idx_j.size();j2++)
#else
                            for (unsigned int j2=0;j2<idx_j.size();j2++)
#endif
                            {
                                double delta=(jN(j1,primo+jmedia)[0]*jN(j2,primo+jmedia+itimestep)[0]+
                                        jN(j1,primo+jmedia)[1]*jN(j2,primo+jmedia+itimestep)[1]+
                                        jN(j1,primo+jmedia)[2]*jN(j2,primo+jmedia+itimestep)[2]
                                        )/3.0 - JJ[idxj];
                                JJ[idxj]+=delta/(cont);
                                idxj++;
                            }

                    }

                    for (unsigned int j1=0;j1<N_corr;j1++)
                    {
                        //questa diventa la media della media nei pezzettini (devo aggiungere un contatore e usare la formula della media)
                        //possibile perdita di precisione!
                        double delta_JJ=JJ[j1]-lista[(itimestep)*narr+j1];
                        lista[(itimestep)*narr+j1]+=delta_JJ/cont_JJ;
                    }

                    //N_corr (funzioni di correlazione), N_corr (integrali,integrali di einstein), 1 (kappa), 1 (kappa_einstein)
                    // totale 3*N_corr+2
                }
                delete [] JJ;

            }));
        }
        for (unsigned int ithread=0;ithread<nthread;ithread++){
            threads[ithread].join();
        }
        threads.clear();
    }
    if (bench){
        delete [] intJJ;
        delete [] int_ein_JJ;
        delete [] JJm_T;
        delete [] matr;
        return;
    }
    //calcola la media delle funzioni di correlazione

    if (subtract_mean) {
        unsigned int cont_JJm_T=0;
        for (unsigned int itimestep=start_mean;itimestep<leff;itimestep++){
            cont_JJm_T++;
            for (unsigned int j1=0;j1<N_corr;j1++)  {
                double delta=lista[(itimestep)*narr+j1]-JJm_T[j1];
                JJm_T[j1]+=delta/cont_JJm_T;
            }
        }
    }

    // calcola gli integrali
    if (true) {
        unsigned int istart=0;
        if (subtract_mean) { //toglie la media a tutte le funzioni di correlazione prima di fare gli integrali
            for (unsigned int itimestep=0;itimestep<leff;itimestep++)
                for (unsigned int j=0;j<N_corr;j++)
                    lista[(itimestep)*narr+j]-=JJm_T[j];
        }

        for (unsigned int itimestep=istart;itimestep<leff;itimestep++) {

            //calcola tutti gli integrali
            for (unsigned int j=0;j<N_corr;j++){
                if (itimestep>0){
                    intJJ[j]+=lista[(itimestep-1)*narr+j];
                    int_ein_JJ[j]+=lista[(itimestep-1)*narr+j]*itimestep;
                }

                lista[(itimestep)*narr+N_corr+2*j]=intJJ[j]+lista[(itimestep)*narr+j]/2.0;
                lista[(itimestep)*narr+N_corr+2*j+1]=int_ein_JJ[j]+lista[(itimestep)*narr+j]*itimestep/2.0;

            }
            //calcola il coefficiente di conducibilità come 1/(inversa della matrice(0,0))
            for (unsigned int j=0;j<idx_j.size()*idx_j.size();j++){
                matr[j]=0.0;
            }
            unsigned int idxj=0;
            for (unsigned int j1=0;j1<idx_j.size();j1++)
#ifdef HALF_CORR
                for (unsigned int j2=j1;j2<idx_j.size();j2++) {
                    matr[j2*idx_j.size()+j1]=lista[(itimestep)*narr+N_corr+2*idxj];
                    matr[j1*idx_j.size()+j2]=lista[(itimestep)*narr+N_corr+2*idxj];
                    idxj++;
                }
#else
                for (unsigned int j2=0;j2<idx_j.size();j2++) {
                    matr[j2*idx_j.size()+j1]+=lista[(itimestep)*narr+N_corr+2*idxj]/2.0;
                    matr[j1*idx_j.size()+j2]+=lista[(itimestep)*narr+N_corr+2*idxj]/2.0;
                    idxj++;
                }
#endif

            Eigen::Map<Eigen::MatrixXd> coeff(matr,idx_j.size(),idx_j.size());

            //calcola il complemento di schur di (0,0)  -- questo è equivalente alla componente (0,0)^-1 della matrice inversa:

            double k;
            if (idx_j.size()>1)
                k= (coeff.block(0,0,1,1) - coeff.block(0,1,1,idx_j.size()-1)*coeff.block(1,1,idx_j.size()-1,idx_j.size()-1).inverse()*coeff.block(1,0,idx_j.size()-1,1))(0,0);
            else
                k=coeff(0,0);
            lista[(itimestep)*narr+3*N_corr+0]=k;

            //stessa cosa con la formula di einstein
            for (unsigned int j=0;j<idx_j.size()*idx_j.size();j++){
                matr[j]=0.0;
            }
            idxj=0;
            for (unsigned int j1=0;j1<idx_j.size();j1++)
#ifdef HALF_CORR
                for (unsigned int j2=j1;j2<idx_j.size();j2++) {
                    matr[j2*idx_j.size()+j1]=lista[(itimestep)*narr+N_corr+2*idxj]-lista[(itimestep)*narr+N_corr+2*idxj+1]/itimestep;
                    matr[j1*idx_j.size()+j2]=lista[(itimestep)*narr+N_corr+2*idxj]-lista[(itimestep)*narr+N_corr+2*idxj+1]/itimestep;
                    idxj++;
                }
#else
                for (unsigned int j2=0;j2<idx_j.size();j2++) {
                    matr[j2*idx_j.size()+j1]+=(lista[(itimestep)*narr+N_corr+2*idxj]-lista[(itimestep)*narr+N_corr+2*idxj+1]/itimestep)/2.0;
                    matr[j1*idx_j.size()+j2]+=(lista[(itimestep)*narr+N_corr+2*idxj]-lista[(itimestep)*narr+N_corr+2*idxj+1]/itimestep)/2.0;
                    idxj++;
                }
#endif

            //calcola il complemento di schur di (0,0)  -- questo è equivalente alla componente (0,0)^-1 della matrice inversa:

            if (idx_j.size()>1)
                k= (coeff.block(0,0,1,1) - coeff.block(0,1,1,idx_j.size()-1)*coeff.block(1,1,idx_j.size()-1,idx_j.size()-1).inverse()*coeff.block(1,0,idx_j.size()-1,1))(0,0);
            else
                k=coeff(0,0);
            lista[(itimestep)*narr+3*N_corr+1]=k;


        }
    }
    //divide per itimestep tutti gli integrali einsteniani
    for (unsigned int itimestep=0;itimestep<leff;itimestep++) {
        for (unsigned int ieinst=0;ieinst<N_corr;ieinst++){
            lista[(itimestep)*narr+N_corr+2*ieinst+1]/=itimestep;
        }
    }


    if (scrivi_file) {
        std::ofstream outfile(log+".greekdump",std::ios::app);
        for (unsigned int itimestep=0;itimestep<leff;itimestep++) {
            for (unsigned int j=0;j<narr;j++){
                outfile << lista[(itimestep)*narr+j] << " ";
            }
            outfile << "\n";
        }
        outfile << "\n\n";
    }

    delete [] intJJ;
    delete [] int_ein_JJ;
    delete [] JJm_T;
    delete [] matr;


}

bool GreenKuboNComponentIonicFluid::benchmarked=false;

unsigned int GreenKuboNComponentIonicFluid::n_seg_bench(){
    //salva la vecchia dimensione totale e il vecchio numero n_seg
    unsigned int orig_n_seg=n_seg,orig_ntimesteps=ntimesteps,ris=0;
    double min=std::numeric_limits<double>::max() ;
    benchmarked=true;
    bench=true;
    n_seg=1;
    std::cerr << "# k    cputime\n";
    for (unsigned int i=n_seg_start;i<n_seg_stop;i+=5){
        ntimesteps=orig_ntimesteps/i;
        cronometro cron;

        cron.start();
        calcola(0);
        cron.stop();
        std::cerr << i<<" "<<cron.time()*i<<"\n";
        if (cron.time()*i<min){
            min=cron.time()*i;
            ris=i;
        }
    }


    bench=false;
    n_seg=orig_n_seg;
    ntimesteps=orig_ntimesteps;
    return ris;
}
