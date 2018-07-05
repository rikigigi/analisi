/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
  *
**/


#ifndef GREENKUBONCOMPONENTIONICFLUID_H
#define GREENKUBONCOMPONENTIONICFLUID_H


#include "readlog.h"

#include "operazionisulista.h"
#include "mediablocchi.h"
#include <array>
#include <vector>
#include <string>

template <class TFLOAT=double,class TFLOAT_READ=double>
class GreenKuboNComponentIonicFluid : public OperazioniSuLista<GreenKuboNComponentIonicFluid<TFLOAT,TFLOAT_READ>,TFLOAT>
{
public:
    GreenKuboNComponentIonicFluid(ReadLog<TFLOAT_READ> * traiettoria,
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
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void reset(unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    GreenKuboNComponentIonicFluid<TFLOAT,TFLOAT_READ> & operator =(const GreenKuboNComponentIonicFluid<TFLOAT,TFLOAT_READ> &);
    unsigned int get_narr();
    unsigned int get_indexOfKappa();
    unsigned int n_seg_bench();
    private:
    using OperazioniSuLista<GreenKuboNComponentIonicFluid<TFLOAT,TFLOAT_READ>,TFLOAT>::lista;
    using OperazioniSuLista<GreenKuboNComponentIonicFluid<TFLOAT,TFLOAT_READ>,TFLOAT>::lunghezza_lista;
    static bool benchmarked;
    unsigned int narr,N_corr,start_mean,n_seg,n_seg_start,n_seg_stop;
    bool scrivi_file,subtract_mean,bench;
    unsigned int idx_je,idx_j0,idx_j1;
    std::vector<unsigned int> idx_j;
    std::string log;
    ReadLog<TFLOAT_READ> *traiettoria;
    unsigned int ntimesteps,lmax,leff,nthread,skip;
    TFLOAT_READ*  jN(unsigned int N,unsigned int ts);
};



#endif // GREENKUBONCOMPONENTIONICFLUID_H
