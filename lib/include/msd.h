/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef MSD_H
#define MSD_H

#include "operazionisulista.h"
#include "calcolamultithread.h"
#include <vector>

namespace MSD_Flags {
constexpr int FLAGS = CalcolaMultiThread_Flags::PARALLEL_SPLIT_TIME |
         CalcolaMultiThread_Flags::CALL_DEBUG_ROUTINE |
         CalcolaMultiThread_Flags::CALL_CALC_INIT;
}

template <class T,bool FPE=false>
class MSD : public VectorOp<MSD<T,FPE> >, public CalcolaMultiThread<MSD<T,FPE>, MSD_Flags::FLAGS >
{
public:
    MSD(T *t,
        unsigned int skip=1,
        unsigned int tmax=0,
        unsigned int nthreads=0,
        bool calcola_msd_centro_di_massa=false,
        bool calcola_msd_nel_sistema_del_centro_di_massa=false,
        bool debug=false
            );
    std::vector<ssize_t> get_shape() const { return {static_cast<ssize_t> (leff),static_cast<ssize_t>(f_cm),static_cast<ssize_t>(ntypes)} ; }
    std::vector<ssize_t> get_stride() const { return { static_cast<ssize_t> (ntypes*f_cm*sizeof(double)),static_cast<ssize_t>(ntypes*sizeof (double)),static_cast<ssize_t>(sizeof(double))};}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    //calcolaMultiThread interface
    void calc_init(size_t primo);
    void calc_single_th(size_t begin, size_t end, size_t primo, unsigned int ith );
    void calc_end();

    MSD<T,FPE> & operator =(const MSD<T,FPE> & destra);
    unsigned int nExtraTimesteps(unsigned int n_b);
    using This = MSD<T, FPE>;
    using CMT = CalcolaMultiThread<This,MSD_Flags::FLAGS>;
private:
    using VectorOp<MSD<T,FPE> >::lista;
    using VectorOp<MSD<T,FPE> >::lunghezza_lista;
    T * traiettoria;
    size_t lmax,f_cm,ntypes;
    using CMT::nthreads;
    using CMT::skip;
    using CMT::ntimesteps;
    using CMT::leff;
    bool cm_msd,cm_self,debug;
};

#endif // MSD_H
