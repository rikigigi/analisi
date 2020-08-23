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
                  std::array<size_t,3> nbin,
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

    std::array<size_t,3> nbin;
    Hist * hist;
    T * t;
    size_t idx(const size_t x, const size_t y, const size_t z) const {
        return x + nbin[0]*(y + nbin[1]*z);
    }
    size_t idx_(const double * pos, const double * l, bool & ok) const {
        size_t idxs[3];
        for (unsigned int icoord=0;icoord<3;++icoord){
            double l_=(pos[icoord]-l[2*icoord])/(l[2*icoord+1]-l[2*icoord]);
            if (l_<0 || l_>=1) {
                ok=false;
                return 0;
            } else {
                ok=true;
            }
            idxs[icoord] = std::floor(l_*nbin[icoord]);
        }
        return idx(idxs[0],idxs[1],idxs[2]);
    }

    using CalcolaMultiThread<This>::ntimesteps;
    using CalcolaMultiThread<This>::skip;
    using CalcolaMultiThread<This>::nthreads;
    using OperazioniSuLista<This,Hist>::lunghezza_lista;
    using OperazioniSuLista<This,Hist>::lista;

};

#endif // ATOMICDENSITY_H
