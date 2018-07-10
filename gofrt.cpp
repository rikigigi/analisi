#include "gofrt.h"
#include <cmath>
#include<thread>
#include <vector>
#include <fstream>

Gofrt::Gofrt(Traiettoria *t, double rmin, double rmax, unsigned int nbin, unsigned int tmax, unsigned int nthreads, unsigned int skip, bool debug) :
    traiettoria(t),rmin(rmin),rmax(rmax),nbin(nbin), skip(skip), lmax(tmax), nthreads(nthreads), debug(debug)
{

}

Gofrt::~Gofrt() {

}

unsigned int Gofrt::numeroTimestepsOltreFineBlocco(unsigned int n_b){
    return (traiettoria->get_ntimesteps()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->get_ntimesteps()/(n_b+1)+1 : lmax;
}

void Gofrt::reset(const unsigned int numeroTimestepsPerBlocco) {

    //lunghezza in timestep
    leff =(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    //numero di timestep su cui fare la media
    ntimesteps=numeroTimestepsPerBlocco;
    lunghezza_lista=leff*traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2*nbin;

    delete [] lista;
    lista=new double [lunghezza_lista];


}

double * Gofrt::gofr(unsigned int ts,unsigned int itype) {
    return &lista[ts*traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2*nbin];
}

void Gofrt::calcola(unsigned int primo) {

    if (nthreads<=1){
        std::cerr << "Attenzione: sto usando un solo thread.\n";
        nthreads=1;
    }

    std::vector<std::thread> threads;

    for (unsigned int  ith=0;ith<nthreads;ith++) {
        threads.push_back(std::thread([&,ith](){
///TODO: fare le cose.
        }));
    }
    for (unsigned int  ith=0;ith<nthreads;ith++)
        threads[ith].join();
    threads.clear();

    if (debug) {
        std::ofstream out ("analisi_gofrt.debug",std::fstream::app);
        for (unsigned int ts=0;ts<leff;ts++) {
            for (unsigned int r=0;r<nbin;r++){
                out << ts<< " "<< r;
                for (unsigned int itype=0;itype<traiettoria->get_ntypes()*(itype<traiettoria->get_ntypes()+1)/2;itype++){
                    out << " "<< gofr(ts,itype)[r];
                }
                out << "\n";
            }
        }
        out << "\n\n";
    }

}

Gofrt & Gofrt::operator =(const Gofrt &destra) {
    OperazioniSuLista<Gofrt>::operator = (destra);
    return *this;
}
