#ifndef ATOMICDENSITY_H
#define ATOMICDENSITY_H

#include "calcolamultithread.h"
#include "operazionisulista.h"
#include <array>
#include <cmath>

template <class T, class Hist >
class AtomicDensity : public CalcolaMultiThread<AtomicDensity<T,Hist>>, public OperazioniSuLista<AtomicDensity<T,Hist>,Hist>
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
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void reset(const unsigned int numeroTimestepsPerBlocco) ;
    std::vector<ssize_t> get_shape() const;
    std::vector<ssize_t> get_stride() const;
    using This = AtomicDensity<T,Hist>;
    using OperazioniSuLista<This,Hist>::azzera;

private:

    std::array<ssize_t,3> nbin;
    Hist * hist;
    T * t;
    int ntypes;
    size_t idx(const size_t x, const size_t y, const size_t z, const size_t itype) const {
        return x + nbin[2]*(y + nbin[1]*(z + itype*nbin[0]));
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
        return idx(idxs[2],idxs[1],idxs[0],itype);
    }

    using CalcolaMultiThread<This>::ntimesteps;
    using CalcolaMultiThread<This>::skip;
    using CalcolaMultiThread<This>::nthreads;
    using OperazioniSuLista<This,Hist>::lunghezza_lista;
    using OperazioniSuLista<This,Hist>::lista;

};

#endif // ATOMICDENSITY_H
