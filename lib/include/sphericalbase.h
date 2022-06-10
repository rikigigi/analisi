#ifndef SPHERICALBASE_H
#define SPHERICALBASE_H
#include <vector>
#include "neighbour.h"

template <int l, class TFLOAT, class T>
class SphericalBase
{
public:
    using Rminmax_t = std::vector<std::pair<TFLOAT,TFLOAT> >;
    using Neighbours_T = Neighbours<T,double>;
    SphericalBase(T*t,
                  const size_t nbin,
                  const Rminmax_t rminmax);
    void check_rminmax_size();
    void calc(int timestep,
              TFLOAT *result,// get_result_size() size
              TFLOAT *workspace,//workspace array TFLOAT workspace[(lmax+1)*(lmax+1)];
              TFLOAT * cheby,//workspace array, TFLOAT cheby[2*(lmax+1)];
              int * counter=nullptr, // counter of #of atoms in the bin, for every atom
              Neighbours_T * nnl = nullptr /// auto nnl = Neighbours<T,double>{&t,detail.neighList}; in the caller, if needed
              ) const;

    size_t get_single_atom_size() const {
        return (l+1)*(l+1)*nbin*ntypes;
    }
    size_t get_result_size()const{
        return get_single_atom_size()*natoms;
    }

private:
    T & t;
    const size_t natoms;
    const size_t ntypes;
    const size_t nbin;
    const struct _detail {
        _detail(const Rminmax_t rminmax,const TFLOAT nbin) :
            rminmax{rminmax}
        {
            dr.reserve(rminmax.size());
            for (auto & range : rminmax) {
              dr.push_back((range.second-range.first)/nbin);
            }
        }
        Rminmax_t rminmax;
        std::vector<TFLOAT> dr;
    } detail;
protected:

    inline int index_wrk_counter(const int iatom,const int jtype,const int ibin=0) const noexcept {
        return (nbin*(ntypes*iatom+jtype)+ibin);
    }

    inline int index_wrk(const int iatom,const int jtype,const int ibin=0) const noexcept {
        return (l+1)*(l+1)*index_wrk_counter(iatom,jtype,ibin);
    }
};

#endif // SPHERICALBASE_H
