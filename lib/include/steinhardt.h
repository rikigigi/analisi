#ifndef STEINHARDT_H
#define STEINHARDT_H

#include "sphericalcorrelations.h"
#include "calcolamultithread.h"

namespace Steinhardt_Flags {

constexpr int FLAGS =
        CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE |
        CalcolaMultiThread_Flags::CALL_INNER_JOIN_DATA |
        CalcolaMultiThread_Flags::CALL_CALC_INIT;

}

template <int l, class TFLOAT, class T>
class Steinhardt : public SphericalCorrelations<l,TFLOAT,T>, public CalcolaMultiThread< Steinhardt<l,TFLOAT,T>, Steinhardt_Flags::FLAGS > {
public:
    using SPHC = SphericalCorrelations<l,TFLOAT,T>;
    using CMT = CalcolaMultiThread<Steinhardt<l,TFLOAT,T>, Steinhardt_Flags::FLAGS >;
    using rminmax_t = typename SPHC::rminmax_t;
    using CMT::calcola;
    Steinhardt(T *t,
               const rminmax_t rminmax,
               unsigned int nbin,
               unsigned int steinhardt_histogram_size,
               std::vector<unsigned int> steinhardt_l_histogram, /// list of l to use to make an histogram. If the list is not emplty, use the tmax argument as the size of each dimension of the histogram
               unsigned int nthreads=2,
               unsigned int skip=1,
               bool debug=false
            );

    unsigned int numeroTimestepsOltreFineBlocco(unsigned int nb) {return 0;}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calc_init(int primo);
    void calc_single_th(int,//average index, begin
                        int,//average index, end
                        int,//first index of the block
                        int//thread index
                        );
    void join_data();

    std::vector<size_t> compute_stride() {
        stride = {steinhardt_histogram_size*ntypes*ntypes,steinhardt_histogram_size*ntypes,steinhardt_histogram_size};
        size_t strideh=steinhardt_histogram_size;
        for (const auto & i : steinhardt_l_histogram) {
            strideh/=nbin_steinhardt;
            stride.push_back(strideh);
        }
        return stride;
    }

    //for the python interface
    std::vector<ssize_t> get_stride() const {
        std::vector<ssize_t> nstride;
        for (const auto & s : stride){
            nstride.push_back(s*sizeof(TFLOAT));
        }
        return nstride;
    }

    std::vector<ssize_t> get_shape(){
        std::vector<ssize_t> res {nbin,ntypes,ntypes};
        for (int i=0;i<steinhardt_l_histogram.size();++i){
            res.push_back(nbin_steinhardt);
        }
        return res;
    }

    size_t get_index(const int ibin, const int type1, const int type2) const {
        return stride[0]*ibin+stride[1]*type1+stride[2]*type2;
    }

    size_t get_index(const int ibin, const int type1, const int type2, const int * hists ) const {
        size_t idx = stride[0]*ibin+stride[1]*type1+stride[2]*type2;
        for (unsigned i=0;i<steinhardt_l_histogram.size();++i) {
            idx+=stride[3+i]*hists[i];
        }
        return idx;
    }


private:
    using CMT::ntimesteps; //note that due to mess, also SPHC has a variable with the same name
    using SPHC::natoms;
    using SPHC::ntypes;
    using SPHC::c_descr;
    using SPHC::t;
    using SPHC::lunghezza_lista;
    using SPHC::nbin;
    using SPHC::lista;
    const std::vector<unsigned int> steinhardt_l_histogram;
    size_t nbin_steinhardt,steinhardt_histogram_size;
    TFLOAT * threadResults;
    TFLOAT incr;
    std::vector<size_t> stride;
};

#endif // STEINHARDT_H
