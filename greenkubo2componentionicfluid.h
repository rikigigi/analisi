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

#ifndef GREENKUBO2COMPONENTIONICFLUID_H
#define GREENKUBO2COMPONENTIONICFLUID_H

#include "readlog.h"

#include "operazionisulista.h"
#include "mediablocchi.h"
#include <array>


class GreenKubo2ComponentIonicFluid : public OperazioniSuLista<GreenKubo2ComponentIonicFluid>
{
public:
    GreenKubo2ComponentIonicFluid(ReadLog * traiettoria,
                                  std::string log,
                                  double *cariche,
                                  unsigned int skip,
                                  bool dump=false,
                                  unsigned int lunghezza_funzione_max=0,
                                  unsigned int nthreads=0,
				  unsigned int n_ris=100
            );
    ~GreenKubo2ComponentIonicFluid();
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void reset(unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    GreenKubo2ComponentIonicFluid & operator =(const GreenKubo2ComponentIonicFluid &);
    static const unsigned int narr;
private:
    bool scrivi_file;
    unsigned int idx_je,idx_j0,idx_j1,n_ris;
    std::string log;
    ReadLog *traiettoria;
    unsigned int ntimesteps,lmax,leff,nthread,skip;
    double carica[2];
    std::array<double,3>  jz(unsigned int ts);
    double* je(unsigned int ts);
};

#endif // GREENKUBO2COMPONENTIONICFLUID_H
