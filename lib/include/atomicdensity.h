#ifndef ATOMICDENSITY_H
#define ATOMICDENSITY_H

#include "calculatemultithread.h"
#include "operazionisulista.h"
#include <array>
#include <cmath>

template <class T, class Hist >
class AtomicDensity : public CalculateMultiThread<AtomicDensity<T,Hist>>, public VectorOp<AtomicDensity<T,Hist>,Hist>
{
public:
    AtomicDensity(T *t,
                  std::array<ssize_t,3> nbin,
    unsigned int nthreads=0,
    unsigned int skip=1
            );
    ~AtomicDensity();
    void calc_single_th(const unsigned int &start, const unsigned int & stop, const unsigned int & primo, const unsigned int & ith);
    void join_data();
    unsigned int nExtraTimesteps(unsigned int n_b);
    void reset(const unsigned int numeroTimestepsPerBlocco) ;
    std::vector<ssize_t> get_shape() const;
    std::vector<ssize_t> get_stride() const;
    using This = AtomicDensity<T,Hist>;
    using VectorOp<This,Hist>::azzera;

private:

    std::array<ssize_t,3> nbin;
    Hist * hist;
    T * t;
    size_t ntypes;
    size_t idx(const size_t x, const size_t y, const size_t z, const size_t itype) const {
        return x + nbin[0]*(y + nbin[1]*(z + itype*nbin[2]));
    }
    size_t idx_(const double * pos, const double * l, bool & ok,int itype) const {
        size_t idxs[3];
        for (unsigned int icoord=0;icoord<3;++icoord){
            double b=l[3+icoord]*2;
            double l_=pos[icoord]/b;//-std::floor(pos[icoord]/b);
            if (l_<0 || l_>=1) {
                ok=false;
                std::cerr << l_ << " b="<<b<<" coord="<<pos[icoord]<<std::endl;
                throw std::runtime_error("Error: wrong pbc wrapping");
                return 0;
            } else {
                ok=true;
            }
            idxs[icoord] = std::floor(l_*nbin[icoord]);
        }
        return idx(idxs[0],idxs[1],idxs[2],itype);
    }

    using CalculateMultiThread<This>::ntimesteps;
    using CalculateMultiThread<This>::skip;
    using CalculateMultiThread<This>::nthreads;
    using VectorOp<This,Hist>::data_length;
    using VectorOp<This,Hist>::vdata;

};

#endif // ATOMICDENSITY_H
