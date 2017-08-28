#include "msd.h"
#include <cmath>
#include<thread>
#include <vector>

MSD::MSD(Traiettoria *t, unsigned int skip, unsigned int tmax, unsigned int nthreads) :
    traiettoria(t), lmax(tmax),skip(skip),nthread(nthreads)
{
}

unsigned int MSD::numeroTimestepsOltreFineBlocco(unsigned int n_b){
    return (traiettoria->get_ntimesteps()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->get_ntimesteps()/(n_b+1)+1 : lmax;
}

void MSD::reset(const unsigned int numeroTimestepsPerBlocco) {

    leff=(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    lunghezza_lista=leff*traiettoria->get_ntypes();
    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new double [lunghezza_lista];
}


void MSD::calcola(unsigned int primo) {



    if (nthread<1) {
        unsigned int *cont=new unsigned int [traiettoria->get_ntypes()];

        for (unsigned int t=0;t<leff;t++) {
            for (unsigned int i=0;i<traiettoria->get_ntypes();i++){
                lista[traiettoria->get_ntypes()*t+i]=0.0;
                cont[i]=0;
            }
            for (unsigned int imedia=0;imedia<ntimesteps;imedia+=skip){
                for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                    double delta=(pow(traiettoria->posizioni(primo+imedia,iatom)[0]-traiettoria->posizioni(primo+imedia+t,iatom)[0],2)+
                            pow(traiettoria->posizioni(primo+imedia,iatom)[1]-traiettoria->posizioni(primo+imedia+t,iatom)[1],2)+
                            pow(traiettoria->posizioni(primo+imedia,iatom)[2]-traiettoria->posizioni(primo+imedia+t,iatom)[2],2))/3.0
                            -lista[traiettoria->get_ntypes()*t+traiettoria->get_type(iatom)];
                    lista[traiettoria->get_ntypes()*t+traiettoria->get_type(iatom)]+=delta/(++cont[traiettoria->get_type(iatom)]);

                }
            }
        }
        delete [] cont;
    } else { //nthread >1 (il codice precedente viene tenuto per verificare il corretto funzionamento di quello parallelo

        /*dividi il lavoro in gruppi
         t  --->  [0,leff[
         npassi     --->  leff

         ciascun gruppo avrà npassith=npassi/nthread passi
         l'ultimo deve finire alla fine
        */
        unsigned int npassith=leff/nthread;
        std::vector<std::thread> threads;

        for (unsigned int ith=0;ith<nthread;ith++) {
            threads.push_back(std::thread([&,ith](){
                unsigned int *cont=new unsigned int [traiettoria->get_ntypes()];
                unsigned int ultimo= (ith != nthread-1 )?npassith*(ith+1):leff;

                for (unsigned int t=npassith*ith;t<ultimo;t++) {
                    for (unsigned int i=0;i<traiettoria->get_ntypes();i++){
                        lista[traiettoria->get_ntypes()*t+i]=0.0;
                        cont[i]=0;
                    }
                    for (unsigned int imedia=0;imedia<ntimesteps;imedia+=skip){
                        for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                            double delta=(pow(traiettoria->posizioni(primo+imedia,iatom)[0]-traiettoria->posizioni(primo+imedia+t,iatom)[0],2)+
                                    pow(traiettoria->posizioni(primo+imedia,iatom)[1]-traiettoria->posizioni(primo+imedia+t,iatom)[1],2)+
                                    pow(traiettoria->posizioni(primo+imedia,iatom)[2]-traiettoria->posizioni(primo+imedia+t,iatom)[2],2))/3.0
                                    -lista[traiettoria->get_ntypes()*t+traiettoria->get_type(iatom)];
                            lista[traiettoria->get_ntypes()*t+traiettoria->get_type(iatom)]+=delta/(++cont[traiettoria->get_type(iatom)]);

                        }
                    }
                }
                delete [] cont;
            }));
        }

        for (unsigned int ith=0;ith<nthread;ith++)
            threads[ith].join();

        threads.clear();



    }
}

MSD & MSD::operator=(const MSD & destra) {
    OperazioniSuLista<MSD>::operator =( destra);
    return *this;
}
