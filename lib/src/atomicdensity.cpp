#include "atomicdensity.h"
#include "config.h"

template <class T, class Hist>
AtomicDensity<T,Hist>::AtomicDensity(T *t, std::array<ssize_t, 3> nbin, unsigned int nthreads, unsigned int skip) :
    CalcolaMultiThread<This> {nthreads, skip}, nbin{nbin},t{t}, ntypes{t->get_ntypes()}
{
    data_length=nbin[0]*nbin[1]*nbin[2]*ntypes;
    vdata=new Hist [data_length];
    hist=new Hist [ data_length*(nthreads-1)];
    azzera();
}

template <class T, class Hist>
AtomicDensity<T,Hist>::~AtomicDensity() {
    delete [] hist; // vdata deleted in VectorOp
}

template <class T, class Hist>
void AtomicDensity<T, Hist>::calc_single_th(const unsigned int &start, const unsigned int & stop, const unsigned int & primo, const unsigned int & ith)   {
    if (ntimesteps+primo > t->get_ntimesteps()){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the block");
    }
    Hist * hist_=0;
    if (ith>0) hist_= hist+data_length*(ith-1);
    else       hist_= vdata;
    unsigned int natoms=t->get_natoms();
    for (unsigned int i=start;i<stop; i+=skip) {
        for (unsigned int iatom=0; iatom<natoms;++iatom) {
            int itype=t->get_type(iatom);
            auto * pos = t->posizioni(i+primo,iatom);
            auto * l = t->scatola(i+primo);
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
        for (size_t i=0;i<data_length;++i){
            vdata[i]+=hist[i+ith*data_length];
        }
    }
}

template <class T, class Hist>
unsigned int  AtomicDensity<T,Hist>::nExtraTimesteps(unsigned int n_b)  {return 0;}
template <class T, class Hist>
void  AtomicDensity<T,Hist>::reset(const unsigned int numeroTimestepsPerBlocco) {
    ntimesteps=numeroTimestepsPerBlocco;
    azzera();
    for (int i=0;i<data_length*(nthreads-1);++i)
        hist[i]=0;
}
template <class T, class Hist>
std::vector<ssize_t>  AtomicDensity<T,Hist>::get_shape() const{
    return {static_cast<long>(ntypes),static_cast<long>(nbin[0]),static_cast<long>(nbin[1]),static_cast<long>(nbin[2])};
}
template <class T, class Hist>
std::vector<ssize_t>  AtomicDensity<T,Hist>::get_stride() const {
    return {static_cast<long>(sizeof (Hist))*static_cast<long>(nbin[0]*nbin[1]*nbin[2]),static_cast<long>(nbin[1]*nbin[2]*sizeof (Hist)),static_cast<long>(nbin[2]*sizeof (Hist)),static_cast<long>(sizeof (Hist))};
}

#ifdef BUILD_MMAP
#include "traiettoria.h"
template class AtomicDensity<Traiettoria,long>;
#endif

#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class AtomicDensity<Traiettoria_numpy,long>;

#endif

