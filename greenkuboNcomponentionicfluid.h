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


class GreenKuboNComponentIonicFluid : public OperazioniSuLista<GreenKuboNComponentIonicFluid>
{
public:
    GreenKuboNComponentIonicFluid(ReadLog * traiettoria,
                                  std::string log,
                                  double *cariche,
                                  unsigned int skip,
				  std::vector<string> headers,
                                  bool dump=false,
                                  unsigned int lunghezza_funzione_max=0,
                                  unsigned int nthreads=0,
				  unsigned int n_ris=100
            );
    ~GreenKuboNComponentIonicFluid();
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void reset(unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    GreenKuboNComponentIonicFluid & operator =(const GreenKuboNComponentIonicFluid &);
    static const unsigned int narr;
private:
    bool scrivi_file;
    unsigned int idx_je,idx_j0,idx_j1,n_ris;
    std:vector<unsigned int> idx_j;
    std::string log;
    ReadLog *traiettoria;
    unsigned int ntimesteps,lmax,leff,nthread,skip;
    double carica[2];
    double*  jN(unsigned int N,unsigned int ts);
    std::array<double,3>  jz(unsigned int ts);
    double* je(unsigned int ts);
};



#endif // GREENKUBONCOMPONENTIONICFLUID_H