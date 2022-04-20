/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "Skt.h"
#include <cmath>
#include<thread>
#include <vector>
#include <fstream>
#include <sstream>
#include "config.h"

#ifdef USE_MPI
#include "mp.h"
#endif

template <class TFLOAT,class T> Skt<TFLOAT,T>::Skt(T *t, TFLOAT rmin, TFLOAT rmax ,TFLOAT kmin, TFLOAT kmax, unsigned int nk , unsigned int tmax, unsigned int nthreads, unsigned int skip,unsigned int every, bool debug) :
    CalcolaMultiThread_T {nthreads, skip, t->get_natoms(), every},
    traiettoria(t),rmin(rmin),rmax(rmax),kmin(kmin),kmax(kmax),nk(nk), lmax(tmax), debug(debug)
{
    // totn = nkx*nky*nkz;
    dk = (kmax-kmin)/nk
    dr=(nky-nkx)/nbin;
    nky2=nky*nky;
    nkx2=nkx*nkx;

}

template <class TFLOAT, class T> Skt<TFLOAT,T>::~Skt() {

}

template <class TFLOAT, class T> unsigned int Skt<TFLOAT,T>::numeroTimestepsOltreFineBlocco(unsigned int n_b){
    return (traiettoria->get_ntimesteps()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->get_ntimesteps()/(n_b+1)+1 : lmax;
}

template <class TFLOAT, class T> void Skt<TFLOAT,T>::reset(const unsigned int numeroTimestepsPerBlocco) {

    std::stringstream descr;
    descr << "# The first column is the time difference in timesteps, then you have the bin index. Every column after is followed by the variance. Then you have the following: "<<std::endl;
    for (unsigned int t1=0;t1<traiettoria->get_ntypes();t1++) {
        for (unsigned int t2=t1;t2<traiettoria->get_ntypes();t2++){
            descr << "#g("<<t1<<", "<<t2<<"), different atom index: "<<get_itype(t1,t2)*2+3<<std::endl;
            descr << "#g("<<t1<<", "<<t2<<"), same atom index: "<<(get_itype(t1,t2)+traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2)*2+3<<std::endl;
        }
    }
    descr << "# same atom index means that the atom is tracked around and the average self-spread is shown with larger time differences."<<std::endl;
    descr << "# different atom index is something that for t=0 is the traditional g(r) "<<std::endl;
    c_descr=descr.str();

    //lunghezza in timestep
    leff =(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    //numero di timestep su cui fare la media
    ntimesteps=numeroTimestepsPerBlocco;
    lunghezza_lista=leff*traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)*nbin;

    delete [] lista;
    lista=new TFLOAT [lunghezza_lista];
}

template <class TFLOAT, class T> std::vector<ssize_t> Skt<TFLOAT,T>::get_shape(){
    return {leff,traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1),nbin};
}
template <class TFLOAT, class T> std::vector<ssize_t> Skt<TFLOAT,T>::get_stride(){
    return {static_cast<long>(traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)*nbin*sizeof(TFLOAT)),
             static_cast<long>(nbin*sizeof(TFLOAT)), sizeof(TFLOAT)};
}

template <class TFLOAT, class T>
void Skt<TFLOAT,T>::calc_init(int primo) {
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
void Skt<TFLOAT,T>::calc_single_th(int t, int imedia, int atom_start, int atom_stop,int primo, int ith) {
    TFLOAT * th_data_ = th_data + (ith-1)*lunghezza_lista;
    if (ith==0) {
        th_data_ = lista;
    }
    for (unsigned int iatom=atom_start;iatom<atom_stop;iatom++) {
        for (unsigned int jatom=0;jatom<traiettoria->get_natoms();jatom++) {
            double x[3];
            double d2=traiettoria->d2_minImage(iatom,jatom,primo+imedia,primo+imedia+t,x);
            double d = sqrt(d2) ; 
            if (d2>rmax2 || d2<rmin2) continue;
            unsigned int type1=traiettoria->get_type(iatom);
            unsigned int type2=traiettoria->get_type(jatom);
            unsigned int itype=get_itype(type1,type2);
            if (iatom==jatom){
                itype=itype+traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2;
            }

            //calcola il quadrato della distanza della minima immagine
            //aggiorna l'istogramma
            for (size_t k = 0; k < nk; k++)
            {
                int k_ = k*dk+kmin ;
                double kd = k_ * d;
                th_data_[sk_idx(t,itype,k)]+=sin(kd)/kd
            }

        }
    }

}



template <class TFLOAT, class T>
void Skt<TFLOAT,T>::calc_end() {
    for (int ith=0;ith<nthreads-1;ith++) {
        for (int i=0;i<lunghezza_lista;++i) {
            lista[i]+=th_data[ith*lunghezza_lista+i];
        }
    }

    delete [] th_data;
    th_data = nullptr;

    if (debug) {
#ifndef USE_MPI
        std::ofstream out("Skt.dump",std::ios::app);
#else
        std::ofstream out(Mp::mpi().outname("Skt.dump"));
#endif
        for (unsigned int ts=0;ts<leff;ts++) {
            for (unsigned int r=0;r<nbin;r++){
                out << ts<< " "<< r;
                for (unsigned int itype=0;itype<traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1);itype++){
                    out << " "<< *sk(ts,itype,r);
                }
                out << "\n";
            }
        }
        out << "\n\n";
    }

}

template <class TFLOAT, class T> Skt<TFLOAT,T> & Skt<TFLOAT,T>::operator =(const Skt<TFLOAT,T> &destra) {
    OperazioniSuLista<Skt<TFLOAT,T>,TFLOAT >::operator = (destra);
    return *this;
}

#ifdef BUILD_MMAP
template class Skt<double,Traiettoria>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class Skt<double, Traiettoria_numpy>;
#endif
