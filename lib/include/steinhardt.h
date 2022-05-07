#ifndef STEINHARDT_H
#define STEINHARDT_H

#include "sphericalbase.h"
#include "calcolamultithread.h"
#include "operazionisulista.h"

namespace Steinhardt_Flags {

constexpr int FLAGS =
        CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE |
        CalcolaMultiThread_Flags::CALL_INNER_JOIN_DATA |
        CalcolaMultiThread_Flags::CALL_CALC_INIT;

}

template <int l, class TFLOAT, class T>
class Steinhardt : public SphericalBase<l,TFLOAT,T>,
        public CalcolaMultiThread< Steinhardt<l,TFLOAT,T>, Steinhardt_Flags::FLAGS >,
        public OperazioniSuLista<Steinhardt<l,TFLOAT,T>,TFLOAT>
{
public:
    using LISTA = OperazioniSuLista<Steinhardt<l,TFLOAT,T>,TFLOAT>;
    using SPB = SphericalBase<l,TFLOAT,T>;
    using CMT = CalcolaMultiThread<Steinhardt<l,TFLOAT,T>, Steinhardt_Flags::FLAGS >;
    using CMT::calcola;
    using SPB::calc;
    using typename SPB::Rminmax_t;
    using typename SPB::Neighbours_T;
    using NeighListSpec = typename Neighbours_T::ListSpec;
    Steinhardt(T *t,
               const Rminmax_t rminmax,
               unsigned int nbin,
               unsigned int steinhardt_histogram_size,
               std::vector<unsigned int> steinhardt_l_histogram, /// list of l to use to make an histogram.
               unsigned int nthreads=2,
               unsigned int skip=1,
               bool do_histogram=true,
               const NeighListSpec nls={}
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
        if (do_histogram){
            stride = {steinhardt_histogram_size*ntypes*ntypes,steinhardt_histogram_size*ntypes,steinhardt_histogram_size};
            size_t strideh=steinhardt_histogram_size;
            for (const auto & i : steinhardt_l_histogram) {
                strideh/=nbin_steinhardt;
                stride.push_back(strideh);
            }
        } else {
            stride = {ntimesteps*natoms*ntypes*l,natoms*ntypes*l, ntypes*l,l,1};
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
        if (do_histogram){
            std::vector<ssize_t> res {(ssize_t)nbin,(ssize_t)ntypes,(ssize_t)ntypes};
            for (int i=0;i<steinhardt_l_histogram.size();++i){
                res.push_back(nbin_steinhardt);
            }
            return res;
        } else {
            return std::vector<ssize_t> {(ssize_t)nbin,(ssize_t)ntimesteps,(ssize_t)natoms,(ssize_t)ntypes,(ssize_t)l};
        }
    }

    size_t get_index(const size_t ibin, const size_t itimestep, const size_t iatom, const size_t itype, const size_t il) const {
        return stride[0]*ibin + stride[1]*itimestep + stride[2]*iatom + stride[3]*itype + stride[4]*il;
    }

    size_t get_index(const size_t ibin, const size_t type1, const size_t type2) const {
        return stride[0]*ibin+stride[1]*type1+stride[2]*type2;
    }

    size_t get_index(const size_t ibin, const size_t type1, const size_t type2, const size_t * hists ) const {
        size_t idx = stride[0]*ibin+stride[1]*type1+stride[2]*type2;
        for (unsigned i=0;i<steinhardt_l_histogram.size();++i) {
            idx+=stride[3+i]*hists[i];
        }
        return idx;
    }


private:
    using CMT::ntimesteps;
    const bool do_histogram;
    const size_t natoms;
    const size_t ntypes;
    const size_t nbin;
    std::string c_descr;
    T & t;
    using LISTA::lunghezza_lista;
    using LISTA::lista;
    const std::vector<unsigned int> steinhardt_l_histogram;
    size_t nbin_steinhardt,steinhardt_histogram_size;
    TFLOAT * threadResults;
    TFLOAT incr;
    std::vector<size_t> stride;
    const NeighListSpec neighListSpec;
};

#endif // STEINHARDT_H
