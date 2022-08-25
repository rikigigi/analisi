/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/


#ifndef GREENKUBONCOMPONENTIONICFLUID_H
#define GREENKUBONCOMPONENTIONICFLUID_H



#include "operazionisulista.h"
#include "mediablocchi.h"
#include <array>
#include <vector>
#include <string>

template <class READLOG,class TFLOAT=double, class TFLOAT_READ=double >
class GreenKuboNComponentIonicFluid : public VectorOp<GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ>,TFLOAT>
{
public:
    GreenKuboNComponentIonicFluid(READLOG * traiettoria,
                                  std::string log,
                                  unsigned int skip,
                                  std::vector<std::string> headers,
                                  bool dump=false,
                                  unsigned int lunghezza_funzione_max=0,
                                  unsigned int nthreads=0,
                                  bool subtract_mean=false,
                                  unsigned int start_mean=100,
                                  unsigned int n_seg=1,
                                  bool do_bench=false,
                                  unsigned int n_seg_start=10,
                                  unsigned int n_seg_stop=100
            );
    ~GreenKuboNComponentIonicFluid();
    unsigned int nExtraTimesteps(unsigned int n_b);
    void reset(unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ> & operator =(const GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ> &);
    unsigned int get_narr();
    unsigned int get_indexOfKappa();
    unsigned int n_seg_bench();

    std::string get_columns_description();
    std::vector<ssize_t> get_shape();
    std::vector<ssize_t> get_stride();

    private:
    using VectorOp<GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ>,TFLOAT>::lista;
    using VectorOp<GreenKuboNComponentIonicFluid<READLOG,TFLOAT,TFLOAT_READ>,TFLOAT>::lunghezza_lista;
    static bool benchmarked;
    unsigned int narr,N_corr,start_mean,n_seg,n_seg_start,n_seg_stop;
    bool scrivi_file,subtract_mean,bench;
    unsigned int idx_je,idx_j0,idx_j1;
    std::vector<unsigned int> idx_j;
    std::string log,c_descr;
    READLOG *traiettoria;
    unsigned int ntimesteps,lmax,leff,nthread,skip;
    TFLOAT_READ*  jN(unsigned int N,unsigned int ts);
};



#endif // GREENKUBONCOMPONENTIONICFLUID_H
