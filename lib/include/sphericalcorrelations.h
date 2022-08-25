#ifndef SPHERICALCORRELATIONS_H
#define SPHERICALCORRELATIONS_H

#include "operazionisulista.h"
#include "neighbour.h"
#include "sphericalbase.h"

template <int l,class TFLOAT, class T>
class SphericalCorrelations : public VectorOp<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>,
        public SphericalBase<l,TFLOAT,T>
{
public:
    using rminmax_t = std::vector<std::pair<TFLOAT,TFLOAT> >;
    using NeighListSpec = typename Neighbours<T,double>::ListSpec;
    using SPB = SphericalBase<l,TFLOAT,T>;
    using Neighbours_T = typename SPB::Neighbours_T;
    SphericalCorrelations(T *t,
                          const rminmax_t rminmax,
                          size_t nbin,
                          size_t tmax=0,
                          size_t nthreads=0,
                          size_t skip=1,
                          size_t buffer_size=10,
                          bool debug=false,
                          const NeighListSpec neighList={}
                         );
    ~SphericalCorrelations();
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int);
    unsigned int nExtraTimesteps(unsigned int n_b);
    SphericalCorrelations<l,TFLOAT,T> & operator =(const SphericalCorrelations<l,TFLOAT,T> & destra){
        VectorOp<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::operator = (destra);
        return *this;
    }
    const std::vector<ssize_t> get_shape()const { return {(ssize_t)leff,(ssize_t)ntypes,(ssize_t)ntypes,(ssize_t)nbin,(ssize_t)(l+1)};}
    const std::vector<ssize_t> get_stride()const {
        auto s=get_shape();
        auto s_old=s[s.size()-1];
        s[s.size()-1]=sizeof (TFLOAT);
        for (int i=s.size()-2;i>=0;--i) {
            auto t=s_old*s[i+1];
            s_old=s[i];
            s[i]=t;
        }
        return s;
    }

    std::string get_columns_description() {return c_descr;}
    inline int index(const int t, const int type1, const int type2,const int ibin=0) const noexcept {
        return (l+1)*(nbin*(ntypes*(ntypes*t + type1) + type2)+ibin);
    }
    using VectorOp<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::azzera;
    size_t get_single_type_size() const {
        return SPB::get_single_atom_size();
    }
    size_t get_snap_size()const{
        return SPB::get_result_size();
    }
    size_t get_final_snap_size() const {
        return get_single_type_size()*ntypes/(l+1);
    }
    void corr_sh_calc(const TFLOAT * sh1, const TFLOAT *sh2, TFLOAT * aveTypes, TFLOAT * aveWork1, int sh_snap_size , int sh_final_size, int *avecont) const noexcept;

    using SPB::calc;
protected:
    using VectorOp<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::lista;
    using VectorOp<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::lunghezza_lista;
    T & t;

    size_t nbin, tmax,nthreads,skip,leff,ntimesteps,buffer_size;
    const size_t ntypes,natoms;
    bool debug;
    std::string c_descr;

private:
    const NeighListSpec neighList;
    using SPB::index_wrk;
    using SPB::index_wrk_counter;
};

#endif // SPHERICALCORRELATIONS_H
