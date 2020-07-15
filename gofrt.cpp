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

template <class TFLOAT,class T> Gofrt<TFLOAT,T>::Gofrt(T *t, TFLOAT rmin, TFLOAT rmax, unsigned int nbin, unsigned int tmax, unsigned int nthreads, unsigned int skip, bool debug) :
    traiettoria(t),rmin(rmin),rmax(rmax),nbin(nbin), skip(skip), lmax(tmax), nthreads(nthreads), debug(debug)
{

    dr=(rmax-rmin)/nbin;


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

template <class TFLOAT, class T> TFLOAT * Gofrt<TFLOAT,T>::gofr(unsigned int ts,unsigned int itype,unsigned int r) {
    unsigned int idx= ts   * traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)*nbin
                     +nbin * itype
                     +r;
    if (idx >= lunghezza_lista) {
        std::cerr << "Errore: indice fuori dal range!\n";
        abort();
    }
    return &lista[idx];
}

template <class TFLOAT, class T> unsigned int Gofrt<TFLOAT,T>::get_itype(unsigned int & type1,unsigned int & type2) const {
     /*
     * xxxxx  ntypes*(ntypes+1)/2 - (m+2)*(m+1)/2 +
     * xxxxo  + altra coordinata (che deve essere la più grande)
     * xxxoo  = indice della coppia nella memoria
     * xxooo
     * xoooo
     *
     * xx
     * x
    */
    if (type2<type1) {
        unsigned int tmp=type2;
        type2=type1;
        type1=tmp;
    }
    return traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2
            -(type2+1)*(type2+2)/2 +type1;
}

template <class TFLOAT, class T> void Gofrt<TFLOAT,T>::calcola(unsigned int primo) {

    if (nthreads<=1){
        std::cerr << "Attenzione: sto usando un solo thread.\n";
        nthreads=1;
    }


    if (leff+ntimesteps+primo > traiettoria->get_ntimesteps()){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the average or the lenght of the time lag");
    }

    unsigned int npassith=leff/nthreads;
    std::vector<std::thread> threads;

    for (unsigned int  ith=0;ith<nthreads;ith++) {
        threads.push_back(std::thread([&,ith](){
            unsigned int ultimo= (ith != nthreads-1 )?npassith*(ith+1):leff;
            unsigned int itimestep = primo;
            double l[3]={traiettoria->scatola(itimestep)[1]-traiettoria->scatola(itimestep)[0],
                         traiettoria->scatola(itimestep)[3]-traiettoria->scatola(itimestep)[2],
                         traiettoria->scatola(itimestep)[5]-traiettoria->scatola(itimestep)[4]};

            for (unsigned int t=npassith*ith;t<ultimo;t++){

                //azzera g(r) per tutti i tipi
                for (unsigned int itype=0;itype<traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1);itype++) {
                    for (unsigned int r=0;r<nbin;r++){
                        *gofr(t,itype,r)=0.0;
                    }
                }
                TFLOAT incr=1.0/int(ntimesteps/skip);
                for (unsigned int imedia=0;imedia<ntimesteps;imedia+=skip){
                    for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                        for (unsigned int jatom=0;jatom<traiettoria->get_natoms();jatom++) {
                            unsigned int type1=traiettoria->get_type(iatom);
                            unsigned int type2=traiettoria->get_type(jatom);
                            unsigned int itype=get_itype(type1,type2);
                            if (iatom==jatom){
                                itype=itype+traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2;
                            }

                            //calcola il quadrato della distanza della minima immagine
                            double d=sqrt(traiettoria->d2_minImage(iatom,jatom,primo+imedia,primo+imedia+t,l));
                            //aggiorna l'istogramma
                            int idx=(int)floorf((d-rmin)/dr);

                            if (idx<nbin && idx >= 0)
                                (*gofr(t,itype,idx))+=incr;


                        }
                    }

                }

            }



        }));
    }
    for (unsigned int  ith=0;ith<nthreads;ith++)
        threads[ith].join();
    threads.clear();

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

template class Gofrt<double,Traiettoria>;
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class Gofrt<double, Traiettoria_numpy>;
#endif
