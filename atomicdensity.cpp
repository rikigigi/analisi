#include "atomicdensity.h"
#include "config.h"

template <class T, class Hist>
AtomicDensity<T,Hist>::AtomicDensity(T *t, std::array<size_t, 3> nbin, unsigned int nthreads, unsigned int skip) : CalcolaMultiThread<This> {nthreads, skip}, nbin{nbin},t{t}, ntypes{t->get_ntypes()}
{
    lunghezza_lista=nbin[0]*nbin[1]*nbin[2]*ntypes;
    lista=new Hist [lunghezza_lista];
    hist=new Hist [ lunghezza_lista*(nthreads-1)];
    azzera();
}

template <class T, class Hist>
AtomicDensity<T,Hist>::~AtomicDensity() {
    delete [] hist; // lista deleted in OperazioniSuLista
}

template <class T, class Hist>
void AtomicDensity<T, Hist>::calc_single_th(const unsigned int &start, const unsigned int & stop, const unsigned int & primo, const unsigned int & ith)   {
    if (ntimesteps+primo > t->get_ntimesteps()){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the block");
    }
    Hist * hist_=0;
    if (ith>0) hist_= hist+lunghezza_lista*(ith-1);
    else       hist_= lista;
    unsigned int natoms=t->get_natoms();
    for (unsigned int i=start;i<stop; i+=skip) {
        for (unsigned int iatom=0; iatom<natoms;++iatom) {
            int itype=t->get_type(iatom);
            auto * pos = t->posizioni(i,iatom);
            auto * l = t->scatola(i);
            bool in_range=false;
            auto ih=idx_(pos,l,in_range,itype);
            if (in_range){
                hist_[ih]+=1;
            }
        }
    }

}

template <class T, class Hist>
void AtomicDensity<T,Hist>::join_data() {
    for (unsigned int ith=0;ith<nthreads-1;++ith) {
        for (size_t i=0;i<lunghezza_lista;++i){
            lista[i]+=hist[i+ith*lunghezza_lista];
        }
    }
}

template <class T, class Hist>
unsigned int  AtomicDensity<T,Hist>::numeroTimestepsOltreFineBlocco(unsigned int n_b)  {return 0;}
template <class T, class Hist>
void  AtomicDensity<T,Hist>::reset(const unsigned int numeroTimestepsPerBlocco) { ntimesteps=numeroTimestepsPerBlocco;                                                                                   azzera(); }
template <class T, class Hist>
std::vector<ssize_t>  AtomicDensity<T,Hist>::get_shape() const{
    return {static_cast<long>(ntypes),static_cast<long>(nbin[0]),static_cast<long>(nbin[1]),static_cast<long>(nbin[2])};
}
template <class T, class Hist>
std::vector<ssize_t>  AtomicDensity<T,Hist>::get_stride() const {
    return {sizeof (Hist)*static_cast<long>(nbin[0]*nbin[1]*nbin[2]),static_cast<long>(nbin[1]*nbin[2]*sizeof (Hist)),static_cast<long>(nbin[2]*sizeof (Hist)),static_cast<long>(sizeof (Hist))};
}

#include "traiettoria.h"
template class AtomicDensity<Traiettoria,long>;

#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class AtomicDensity<Traiettoria_numpy,long>;

#endif

