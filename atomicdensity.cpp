#include "atomicdensity.h"
#include "config.h"

template <class T, class Hist>
AtomicDensity<T,Hist>::AtomicDensity(T *t, std::array<size_t, 3> nbin, unsigned int nthreads, unsigned int skip) : CalcolaMultiThread<This> {nthreads, skip}, nbin{nbin},t{t}
{
    hist=new Hist [ nbin[0]*nbin[1]*nbin[2]*nthreads];
    lista=hist;
    lunghezza_lista=nbin[0]*nbin[1]*nbin[2]*nthreads;
    azzera();
}

template <class T, class Hist>
AtomicDensity<T,Hist>::~AtomicDensity() {
    //delete [] hist; // deleted in OperazioniSuLista
}

template <class T, class Hist>
void AtomicDensity<T, Hist>::calc_single_th(const unsigned int &start, const unsigned int & stop, const unsigned int & primo, const unsigned int & ith)   {
    if (ntimesteps+primo > t->get_ntimesteps()){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the block");
    }
    size_t offset = ith*nbin[0]*nbin[1]*nbin[2];
    unsigned int natoms=t->get_natoms();
    for (unsigned int i=start;i<stop; i+=skip) {
        for (unsigned int iatom=0; iatom<natoms;++iatom) {
            auto * pos = t->posizioni(i,iatom);
            auto * l = t->scatola(i);
            bool in_range=false;
            auto ih=idx_(pos,l,in_range);
            if (in_range){
                hist[offset+ih]++;
            }
        }
    }

}

template <class T, class Hist>
void AtomicDensity<T,Hist>::join_data() {
    for (unsigned int ith=1;ith<nthreads;++ith) {
        for (size_t i=0;i<nbin[0]*nbin[1]*nbin[2];++i){
            hist[i]+=hist[i+ith*nbin[0]*nbin[1]*nbin[2]];
        }
    }
}

template <class T, class Hist>
unsigned int  AtomicDensity<T,Hist>::numeroTimestepsOltreFineBlocco(unsigned int n_b)  {return 0;}
template <class T, class Hist>
void  AtomicDensity<T,Hist>::reset(const unsigned int numeroTimestepsPerBlocco) { ntimesteps=numeroTimestepsPerBlocco;                                                                                   azzera(); }
template <class T, class Hist>
std::vector<ssize_t>  AtomicDensity<T,Hist>::get_shape() const{
    return {static_cast<long>(nbin[0]),static_cast<long>(nbin[1]),static_cast<long>(nbin[2])};
}
template <class T, class Hist>
std::vector<ssize_t>  AtomicDensity<T,Hist>::get_stride() const {
    return {sizeof (Hist),static_cast<long>(nbin[0]*sizeof (Hist)),static_cast<long>(nbin[0]*nbin[1]*sizeof (Hist))};
}

#include "traiettoria.h"
template class AtomicDensity<Traiettoria,long>;

#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class AtomicDensity<Traiettoria_numpy,long>;

#endif

