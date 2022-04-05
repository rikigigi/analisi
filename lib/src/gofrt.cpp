/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "gofrt.h"
#include <cmath>
#include<thread>
#include <vector>
#include <fstream>
#include <sstream>
#include "config.h"

#ifdef USE_MPI
#include "mp.h"
#endif

template <class TFLOAT,class T> Gofrt<TFLOAT,T>::Gofrt(T *t, TFLOAT rmin, TFLOAT rmax, unsigned int nbin, unsigned int tmax, unsigned int nthreads, unsigned int skip,unsigned int every, bool debug) :
    CalcolaMultiThread_T {nthreads, skip, t->get_natoms(), every},
    traiettoria(t),rmin(rmin),rmax(rmax),nbin(nbin), lmax(tmax), debug(debug)
{

    dr=(rmax-rmin)/nbin;
    rmax2=rmax*rmax;
    rmin2=rmin*rmin;

}

template <class TFLOAT, class T> Gofrt<TFLOAT,T>::~Gofrt() {

}

template <class TFLOAT, class T> unsigned int Gofrt<TFLOAT,T>::numeroTimestepsOltreFineBlocco(unsigned int n_b){
    return (traiettoria->get_ntimesteps()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->get_ntimesteps()/(n_b+1)+1 : lmax;
}

template <class TFLOAT, class T> void Gofrt<TFLOAT,T>::reset(const unsigned int numeroTimestepsPerBlocco) {

    std::stringstream descr;
    descr << "# every column is followed by the variance. Then you have the following: "<<std::endl;
    for (unsigned int t1=0;t1<traiettoria->get_ntypes();t1++) {
        for (unsigned int t2=t1;t2<traiettoria->get_ntypes();t2++){
            descr << "#g("<<t1<<", "<<t2<<"), different atom index: "<<get_itype(t1,t2)*2<<std::endl;
            descr << "#g("<<t1<<", "<<t2<<"), same atom index: "<<(get_itype(t1,t2)+traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2)*2<<std::endl;
        }
    }
    c_descr=descr.str();

    //lunghezza in timestep
    leff =(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    //numero di timestep su cui fare la media
    ntimesteps=numeroTimestepsPerBlocco;
    lunghezza_lista=leff*traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)*nbin;

    delete [] lista;
    lista=new TFLOAT [lunghezza_lista];
}

template <class TFLOAT, class T> std::vector<ssize_t> Gofrt<TFLOAT,T>::get_shape(){
    return {leff,traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1),nbin};
}
template <class TFLOAT, class T> std::vector<ssize_t> Gofrt<TFLOAT,T>::get_stride(){
    return {static_cast<long>(traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)*nbin*sizeof(TFLOAT)),
             static_cast<long>(nbin*sizeof(TFLOAT)), sizeof(TFLOAT)};
}

template <class TFLOAT, class T>
void Gofrt<TFLOAT,T>::calc_init(int primo) {
    if (nthreads<=1){
        std::cerr << "WARNING: using a single thread.\n";
        nthreads=1;
    }


    if (leff+ntimesteps+primo > traiettoria->get_ntimesteps() + 1){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the average or the lenght of the time lag");
    }

    //init
    th_data = new TFLOAT[lunghezza_lista*(nthreads-1)];
    azzera();
    for (unsigned int i=0;i<lunghezza_lista*(nthreads-1);++i){
        th_data[i]=0;
    }
    if (ntimesteps/skip>0) incr=1.0/int(ntimesteps/skip);
    else                   incr=1;
}

template <class TFLOAT, class T>
void Gofrt<TFLOAT,T>::calc_single_th(int t, int imedia, int atom_start, int atom_stop,int primo, int ith) {
    TFLOAT * th_data_ = th_data + (ith-1)*lunghezza_lista;
    if (ith==0) {
        th_data_ = lista;
    }
    for (unsigned int iatom=atom_start;iatom<atom_stop;iatom++) {
        for (unsigned int jatom=0;jatom<traiettoria->get_natoms();jatom++) {
            double d=traiettoria->d2_minImage(iatom,jatom,primo+imedia,primo+imedia+t);
            if (d>rmax2 || d<rmin2) continue;
            unsigned int type1=traiettoria->get_type(iatom);
            unsigned int type2=traiettoria->get_type(jatom);
            unsigned int itype=get_itype(type1,type2);
            if (iatom==jatom){
                itype=itype+traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2;
            }

            //calcola il quadrato della distanza della minima immagine
            //aggiorna l'istogramma
            d=sqrt(d);
            int idx=(int)floorf((d-rmin)/dr);

            if (idx<nbin && idx >= 0)
                th_data_[gofr_idx(t,itype,idx)]+=incr;


        }
    }

}



template <class TFLOAT, class T>
void Gofrt<TFLOAT,T>::calc_end() {
    for (int ith=0;ith<nthreads-1;ith++) {
        for (int i=0;i<lunghezza_lista;++i) {
            lista[i]+=th_data[ith*lunghezza_lista+i];
        }
    }

    delete [] th_data;
    th_data = nullptr;

    if (debug) {
#ifndef USE_MPI
        std::ofstream out("gofrt.dump",std::ios::app);
#else
        std::ofstream out(Mp::mpi().outname("gofrt.dump"));
#endif
        for (unsigned int ts=0;ts<leff;ts++) {
            for (unsigned int r=0;r<nbin;r++){
                out << ts<< " "<< r;
                for (unsigned int itype=0;itype<traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1);itype++){
                    out << " "<< *gofr(ts,itype,r);
                }
                out << "\n";
            }
        }
        out << "\n\n";
    }

}

template <class TFLOAT, class T> Gofrt<TFLOAT,T> & Gofrt<TFLOAT,T>::operator =(const Gofrt<TFLOAT,T> &destra) {
    OperazioniSuLista<Gofrt<TFLOAT,T>,TFLOAT >::operator = (destra);
    return *this;
}

#ifdef BUILD_MMAP
template class Gofrt<double,Traiettoria>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class Gofrt<double, Traiettoria_numpy>;
#endif
