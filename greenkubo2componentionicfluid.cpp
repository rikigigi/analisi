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

#include "greenkubo2componentionicfluid.h"
#include "chargefluxts.h"
#include "heatfluxts.h"
#include <string>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>

GreenKubo2ComponentIonicFluid::GreenKubo2ComponentIonicFluid(ReadLog *traiettoria, std::string log, double * cariche, unsigned int skip, bool dump, unsigned int lunghezza_funzione_max, unsigned int nthreads,unsigned int n_ris) : OperazioniSuLista<GreenKubo2ComponentIonicFluid>(),
    traiettoria (traiettoria), log(log), ntimesteps(0),skip(skip), scrivi_file(dump),lmax(lunghezza_funzione_max),nthread(nthreads),n_ris(n_ris)
{
    carica[0]=cariche[0];
    carica[1]=cariche[1];

    std::pair<unsigned int ,bool> res;
    res=traiettoria->get_index_of("c_flusso[1]");
    if(res.second)
        idx_je=res.first;
    else
    {std::cerr << "Errore: header di una corrente non trovato!";assert(0);abort();}
    res=traiettoria->get_index_of("c_vcm[1][1]");
    if(res.second)
        idx_j0=res.first;
    else
    {std::cerr << "Errore: header di una corrente non trovato!";assert(0);abort();}
    res=traiettoria->get_index_of("c_vcm[2][1]");
    if(res.second)
        idx_j1=res.first;
    else
    {std::cerr << "Errore: header di una corrente non trovato!";assert(0);abort();}

}

GreenKubo2ComponentIonicFluid::~GreenKubo2ComponentIonicFluid(){
#ifdef DEBUG2
//    std::cerr << "Called delete GreenKubo2ComponentIonicFluid, je="<<je<<", jz="<<jz<<"\n";
#endif
}

GreenKubo2ComponentIonicFluid & GreenKubo2ComponentIonicFluid::operator =(const GreenKubo2ComponentIonicFluid & destra) {
#ifdef DEBUG2
    std::cerr << "Chiamato GreenKubo2ComponentIonicFluid::operator =\n";
#endif
    OperazioniSuLista<GreenKubo2ComponentIonicFluid>::operator =( destra);
    return *this;
}

unsigned int GreenKubo2ComponentIonicFluid::numeroTimestepsOltreFineBlocco(unsigned int n_b) {
    return (traiettoria->n_timestep()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->n_timestep()/(n_b+1)+1 : lmax;
}

void GreenKubo2ComponentIonicFluid::reset(unsigned int numeroTimestepsPerBlocco) {
    leff=(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    lunghezza_lista=(leff)*9; // Jee,Jzz,Jez,Jintee,Jintzz,Jintez,lambda
    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new double [lunghezza_lista];
}

std::array<double,3> GreenKubo2ComponentIonicFluid::jz(unsigned int ts){
    double *j0=&traiettoria->line(ts)[idx_j0];
    double *j1=&traiettoria->line(ts)[idx_j1];
    return std::array<double,3> {{
            j0[0]*carica[0]+j1[0]*carica[1],
            j0[1]*carica[0]+j1[1]*carica[1],
            j0[2]*carica[0]+j1[2]*carica[1],
        }};
}

double* GreenKubo2ComponentIonicFluid::je(unsigned int ts){
    return&traiettoria->line(ts)[idx_je];
}

void GreenKubo2ComponentIonicFluid::calcola(unsigned int primo) {


    unsigned int allinea=0;//primo%skip;

    double intzz=0.0;
    double intee=0.0;
    double intez=0.0;
    double intze=0.0;
    double jeeo=0.0,jzzo=0.0,jezo=0.0,jzeo=0.0;


    if (nthread<1)
    for (unsigned int itimestep=0;itimestep<leff;itimestep++) {
        //fa la media sulla traiettoria dei vari prodotti,
        //con una differenza di timesteps fissata "itimestep"
        double jee=0.0,jzz=0.0,jez=0.0,jze=0.0;
        unsigned int cont=0;
        for (unsigned int jmedia=allinea;jmedia<ntimesteps;jmedia+=skip) {
            //prodotto JzJz
            double deltazz=(jz(primo+jmedia)[0]*jz(primo+jmedia+itimestep)[0]+
                            jz(primo+jmedia)[1]*jz(primo+jmedia+itimestep)[1]+
                            jz(primo+jmedia)[2]*jz(primo+jmedia+itimestep)[2])/3.0
                                - jzz;
            double deltaez=(je(primo+jmedia)[0]*jz(primo+jmedia+itimestep)[0]+
                            je(primo+jmedia)[1]*jz(primo+jmedia+itimestep)[1]+
                            je(primo+jmedia)[2]*jz(primo+jmedia+itimestep)[2])/3.0
                                - jez;
            double deltaze=(je(primo+jmedia+itimestep)[0]*jz(primo+jmedia)[0]+
                            je(primo+jmedia+itimestep)[1]*jz(primo+jmedia)[1]+
                            je(primo+jmedia+itimestep)[2]*jz(primo+jmedia)[2])/3.0
                                - jze;
            double deltaee=(je(primo+jmedia)[0]*je(primo+jmedia+itimestep)[0]+
                            je(primo+jmedia)[1]*je(primo+jmedia+itimestep)[1]+
                            je(primo+jmedia)[2]*je(primo+jmedia+itimestep)[2])/3.0
                                - jee;
            jzz+=deltazz/(++cont);
            jee+=deltaee/(cont);
            jez+=deltaez/(cont);
            jze+=deltaze/(cont);
        }
        lista[(itimestep)*9+0]=jee;
        lista[(itimestep)*9+1]=jzz;
        lista[(itimestep)*9+2]=jez;
        lista[(itimestep)*9+7]=jze;
        //integrale con il metodo dei trapezi

            intee+=jeeo;
            intez+=jezo;
            intze+=jzeo;
            intzz+=jzzo;

        lista[(itimestep)*9+3]=intee+jee/2.0;
        lista[(itimestep)*9+4]=intzz+jzz/2.0;
        lista[(itimestep)*9+5]=intez+jez/2.0;
        lista[(itimestep)*9+6]=intee+jee/2.0-(intez+jez/2.0)*(intez+jez/2.0)/(intzz+jzz/2.0);
        lista[(itimestep)*9+8]=intze+jze/2.0;

        jeeo=jee;
        jezo=jez;
        jzeo=jze;
        jzzo=jzz;
    }
    else { //nthread > 1

        /*dividi il lavoro in gruppi
         itimestep  --->  [0,leff[
         npassi     --->  leff

         ciascun gruppo avrà npassith=npassi/nthread passi
         l'ultimo deve finire alla fine

         gli integrali li calcolo alla fine, è uguale.
        Nel caso di 1 thread tengo l'implementazione precedente
            (anche per poter confrontare i risultati)
        */

        unsigned int npassith=leff/nthread;
        std::vector<std::thread> threads;
        for (unsigned int ith=0;ith<nthread;ith++){
            threads.push_back(std::thread([&,ith](){
                unsigned int ultimo= (ith != nthread-1 )?npassith*(ith+1):leff;
                for (unsigned int itimestep=npassith*ith;itimestep<ultimo;itimestep++) {
                    //fa la media sulla traiettoria dei vari prodotti,
                    //con una differenza di timesteps fissata "itimestep"
                    double jee=0.0,jzz=0.0,jez=0.0,jze=0.0;
                    unsigned int cont=0;
                    for (unsigned int jmedia=allinea;jmedia<ntimesteps;jmedia+=skip) {
                        //prodotto JzJz
                        double deltazz=(jz(primo+jmedia)[0]*jz(primo+jmedia+itimestep)[0]+
                                        jz(primo+jmedia)[1]*jz(primo+jmedia+itimestep)[1]+
                                        jz(primo+jmedia)[2]*jz(primo+jmedia+itimestep)[2])/3.0
                                            - jzz;
                        double deltaez=(je(primo+jmedia)[0]*jz(primo+jmedia+itimestep)[0]+
                                        je(primo+jmedia)[1]*jz(primo+jmedia+itimestep)[1]+
                                        je(primo+jmedia)[2]*jz(primo+jmedia+itimestep)[2])/3.0
                                            - jez;
                        double deltaze=(je(primo+jmedia+itimestep)[0]*jz(primo+jmedia)[0]+
                                        je(primo+jmedia+itimestep)[1]*jz(primo+jmedia)[1]+
                                        je(primo+jmedia+itimestep)[2]*jz(primo+jmedia)[2])/3.0
                                            - jze;
                        double deltaee=(je(primo+jmedia)[0]*je(primo+jmedia+itimestep)[0]+
                                        je(primo+jmedia)[1]*je(primo+jmedia+itimestep)[1]+
                                        je(primo+jmedia)[2]*je(primo+jmedia+itimestep)[2])/3.0
                                            - jee;
                        jzz+=deltazz/(++cont);
                        jee+=deltaee/(cont);
                        jez+=deltaez/(cont);
                        jze+=deltaze/(cont);
                    }
                    lista[(itimestep)*9+0]=jee;
                    lista[(itimestep)*9+1]=jzz;
                    lista[(itimestep)*9+2]=jez;
                    lista[(itimestep)*9+7]=jze;
                    //integrale con il metodo dei trapezi, solo per il primo

                    if (ith==0){
                        intee+=jeeo;
                        intez+=jezo;
                        intze+=jzeo;
                        intzz+=jzzo;

                        lista[(itimestep)*9+3]=intee+jee/2.0;
                        lista[(itimestep)*9+4]=intzz+jzz/2.0;
                        lista[(itimestep)*9+5]=intez+jez/2.0;
                        lista[(itimestep)*9+6]=intee+jee/2.0-(intez+jez/2.0)*(intez+jez/2.0)/(intzz+jzz/2.0);
                        lista[(itimestep)*9+8]=intze+jze/2.0;

                        jeeo=jee;
                        jezo=jez;
                        jzeo=jze;
                        jzzo=jzz;
                    }
                }
            }));
        }
        for (unsigned int ithread=0;ithread<nthread;ithread++)
            threads[ithread].join();
        threads.clear();

        // calcola gli integrali a partire da dove ci eravamo fermati
        if (nthread>1) {
            for (unsigned int itimestep=npassith;itimestep<leff;itimestep++) {

                intee+=lista[((itimestep-1))*9+0];
                intzz+=lista[((itimestep-1))*9+1];
                intez+=lista[((itimestep-1))*9+2];
                intze+=lista[((itimestep-1))*9+7];

                lista[(itimestep)*9+3]=intee+lista[(itimestep)*9+0]/2.0;
                lista[(itimestep)*9+4]=intzz+lista[(itimestep)*9+1]/2.0;
                lista[(itimestep)*9+5]=intez+lista[(itimestep)*9+2]/2.0;
                lista[(itimestep)*9+6]=intee+lista[(itimestep)*9+0]/2.0
                        -(intez+lista[(itimestep)*9+2]/2.0)*(intez+lista[(itimestep)*9+2]/2.0)/(intzz+lista[(itimestep)*9+1]/2.0);
                lista[(itimestep)*9+8]=intze+lista[(itimestep)*9+7]/2.0;

            }
        }
    }

    if (scrivi_file) {
        std::ofstream outfile(log+".greekdump",std::ios::app);
        std::ofstream outfile_s11(log+".greekdump_s11",std::ios::app);
        std::ofstream outfile_s12(log+".greekdump_s12",std::ios::app);
        std::ofstream outfile_s21(log+".greekdump_s21",std::ios::app);
        std::ofstream outfile_s22(log+".greekdump_s22",std::ios::app);
        for (unsigned int itimestep=0;itimestep<leff;itimestep++) {
            for (unsigned int j=0;j<9;j++){
                outfile << lista[(itimestep)*9+j] << " ";
            }
            outfile << "\n";
        }
        outfile_s11 << " " << lista[n_ris*9+0];
        outfile_s22 << " " << lista[n_ris*9+1];
        outfile_s12 << " " << lista[n_ris*9+2];
        outfile_s21 << " " << lista[n_ris*9+7];
        outfile << "\n\n";
    }

}
