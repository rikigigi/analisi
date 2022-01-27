/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/

#include "greenkuboNcomponentionicfluid.h"
#include "chargefluxts.h"
#include "heatfluxts.h"
#include <string>
#include <sstream>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include "cronometro.h"
#include "config.h"

#include "eigen_include.h"

#ifdef USE_MPI
#include "mp.h"
#endif


template<class READLOG, class TFLOAT, class TFLOAT_READ> GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::GreenKuboNComponentIonicFluid(READLOG *traiettoria,
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
                                                             unsigned int n_seg_stop) : OperazioniSuLista<GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ>,TFLOAT>(),
    traiettoria (traiettoria), log(log), ntimesteps(0),skip(skip), scrivi_file(dump),
    lmax(lunghezza_funzione_max),nthread(nthreads),subtract_mean(subtract_mean),
    start_mean(start_mean),n_seg(n_seg),bench(false),
    n_seg_start(n_seg_start), n_seg_stop(n_seg_stop)
{


    if (n_seg<1){
        std::cerr << "Warning: n_seg < 1 . I am setting it to 1.\n";
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
            throw std::runtime_error("Errore: header '" + headers.at(j) + "' not found!\n");
        }
    }

#ifdef HALF_CORR
    N_corr=idx_j.size()*(idx_j.size()+1)/2;
#else
    N_corr=idx_j.size()*idx_j.size();
#endif
    narr=3*N_corr+2;


    unsigned int c=0;
    std::stringstream descr;
    for (unsigned int i=0;i<idx_j.size();i++) {
#ifdef HALF_CORR
        for (unsigned int j=i;j<idx_j.size();j++) {
#else
        for (unsigned int j=0;j<idx_j.size();j++) {
#endif
            descr<<"#"<<++c<<": c("<<i<<", "<<j<<")\n#"<<++c<<": var[c("<<i<<", "<<j<<")]\n";
        }
    }
    for (unsigned int i=0;i<idx_j.size();i++) {
#ifdef HALF_CORR
        for (unsigned int j=i;j<idx_j.size();j++) {
#else
        for (unsigned int j=0;j<idx_j.size();j++) {
#endif
            descr<<"#"<<++c<<": L("<<i<<", "<<j<<")\n#"<<++c<<": var[L("<<i<<", "<<j<<")]\n";
            descr<<"#"<<++c<<": L*("<<i<<", "<<j<<")\n#"<<++c<<": var[L*("<<i<<", "<<j<<")]\n";
        }
    }
    descr<<"#"<<++c<<": k_gk\n#"<<++c<<": var[k_gk]\n";
    descr<<"#"<<++c<<": k_einst\n#"<<++c<<": var[k_einst]\n";
    c_descr=descr.str();
}

template<class READLOG, class TFLOAT, class TFLOAT_READ> GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::~GreenKuboNComponentIonicFluid(){
#ifdef DEBUG2
#endif
}

template<class READLOG, class TFLOAT, class TFLOAT_READ> GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ> &GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::operator =(const GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ> &destra) {
#ifdef DEBUG2
    std::cerr << "Chiamato GreenKuboNComponentIonicFluid<TFLOAT>::operator =\n";
#endif
    OperazioniSuLista<GreenKuboNComponentIonicFluid<READLOG, TFLOAT,TFLOAT_READ>,TFLOAT >::operator =( destra);
    return *this;
}

template<class READLOG, class TFLOAT, class TFLOAT_READ> unsigned int GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::numeroTimestepsOltreFineBlocco(unsigned int n_b) {
    return (traiettoria->n_timestep()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->n_timestep()/(n_b+1)+1 : lmax;
}

template<class READLOG, class TFLOAT, class TFLOAT_READ> void GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::reset(unsigned int numeroTimestepsPerBlocco) {
    leff=(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    lunghezza_lista=(leff)*narr;
    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new TFLOAT [lunghezza_lista];
    for (unsigned int i=0;i<lunghezza_lista;i++){
        lista[i]=0.0;
    }
}
template<class READLOG, class TFLOAT, class TFLOAT_READ>
std::vector<ssize_t> GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ>::get_shape(){
    return {leff,narr};
}
template<class READLOG, class TFLOAT, class TFLOAT_READ>
std::vector<ssize_t> GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ>::get_stride(){
    return {static_cast<long>(narr*sizeof(TFLOAT)),sizeof(TFLOAT)};
}


template<class READLOG, class TFLOAT, class TFLOAT_READ> TFLOAT_READ * GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::jN(unsigned int N,unsigned int ts){
    return&traiettoria->line(ts)[idx_j[N]];
}



template<class READLOG, class TFLOAT, class TFLOAT_READ> unsigned int GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::get_narr(){
    return narr;
}

template<class READLOG, class TFLOAT, class TFLOAT_READ> unsigned int GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::get_indexOfKappa(){
    return 3*N_corr;
}

template<class READLOG, class TFLOAT, class TFLOAT_READ> void GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::calcola(unsigned int primo) {


    if(!benchmarked)
        n_seg_bench();


    if(nthread<1)
        nthread=1;

    //nthread > 1

    if (leff+primo+ntimesteps > traiettoria->n_timestep())
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the average or the lenght of the correlation function");


    TFLOAT *matr=new TFLOAT [idx_j.size()*idx_j.size()];
    TFLOAT *intJJ=new TFLOAT[N_corr];
    TFLOAT *int_ein_JJ=new TFLOAT[N_corr];
    TFLOAT *JJm_T=new TFLOAT[N_corr];
    unsigned int cont_JJ=0;

    for (unsigned int j=0;j<N_corr;j++){
        intJJ[j]=0.0;
        int_ein_JJ[j]=0.0;
        JJm_T[j]=0.0;
    }
    cont_JJ=0;

    std::vector<std::thread> threads;
    unsigned int timesteps_seg=ntimesteps/n_seg;
    for (unsigned int iseg=0;iseg<n_seg;iseg++){
        ++cont_JJ;
        for (unsigned int ith=0;ith<nthread;ith++){
            threads.push_back(std::thread([&,ith](){

                TFLOAT *JJ=new TFLOAT[N_corr];
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
                                TFLOAT delta=(jN(j1,primo+jmedia)[0]*jN(j2,primo+jmedia+itimestep)[0]+
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
                        //possibile perdita di precisione(?)
                        TFLOAT delta_JJ=JJ[j1]-lista[(itimestep)*narr+j1];
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
                TFLOAT delta=lista[(itimestep)*narr+j1]-JJm_T[j1];
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

        for (unsigned int itimestep=1;itimestep<leff;itimestep++) {

            //calcola tutti gli integrali (metodo dei trapezi)
        // I[i] = I[i-1]+ f[i-1]/2.0 + f[i]/2.0
            for (unsigned int j=0;j<N_corr;j++){

            lista[(itimestep)*narr+N_corr+2*j]  = lista[(itimestep-1)*narr+N_corr+2*j]    +lista[(itimestep-1)*narr+j]/2.0 +          lista[(itimestep)*narr+j]/2.0;
                lista[(itimestep)*narr+N_corr+2*j+1]= lista[(itimestep-1)*narr+N_corr+2*j+1]  +lista[(itimestep-1)*narr+j]*itimestep/2.0 +lista[(itimestep)*narr+j]*itimestep/2.0  ;

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

            Eigen::Map<Eigen::Matrix<TFLOAT, Eigen::Dynamic, Eigen::Dynamic> > coeff(matr,idx_j.size(),idx_j.size());

            //calcola il complemento di schur di (0,0)  -- questo è equivalente alla componente (0,0)^-1 della matrice inversa:

            TFLOAT k;
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
    for (unsigned int itimestep=1;itimestep<leff;itimestep++) {
        for (unsigned int ieinst=0;ieinst<N_corr;ieinst++){
            lista[(itimestep)*narr+N_corr+2*ieinst+1]/=itimestep;
        }
    }


    if (scrivi_file) {
#ifndef USE_MPI
        std::ofstream outfile(log+".greekdump",std::ios::app);
#else
        std::ofstream outfile(Mp::mpi().outname(log+".gkdump"),std::ios::app);
#endif
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

template<class READLOG, class TFLOAT, class TFLOAT_READ> bool GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::benchmarked=false;

template<class READLOG, class TFLOAT, class TFLOAT_READ> unsigned int GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::n_seg_bench(){
    //salva la vecchia dimensione totale e il vecchio numero n_seg
    unsigned int orig_n_seg=n_seg,orig_ntimesteps=ntimesteps,ris=0;
    TFLOAT min=std::numeric_limits<TFLOAT>::max() ;
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

template<class READLOG, class TFLOAT, class TFLOAT_READ> std::string GreenKuboNComponentIonicFluid<READLOG, TFLOAT, TFLOAT_READ>::get_columns_description() {
    return c_descr;
}



#ifdef BUILD_MMAP
#include "readlog.h"
template class GreenKuboNComponentIonicFluid<ReadLog<double>, double,double>;
template class GreenKuboNComponentIonicFluid<ReadLog<double>, long double,double>;
template class GreenKuboNComponentIonicFluid<ReadLog<long double>, long double,long double>;
#endif

#ifdef PYTHON_SUPPORT
#include "readlog_numpy.h"
template class GreenKuboNComponentIonicFluid<ReadLog_numpy<double>,double,double>;
#endif
