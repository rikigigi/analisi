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

#include "traiettoria.h"
#include "operazionisulista.h"
#include "mediablocchi.h"

class HeatFluxTs;
class ChargeFluxTs;

class GreenKubo2ComponentIonicFluid : public OperazioniSuLista<GreenKubo2ComponentIonicFluid>
{
public:
    GreenKubo2ComponentIonicFluid(Traiettoria * t,
                                  std::string log,
                                  double *cariche,
                                  unsigned int skip=1,
                                  bool dump=false,
                                  unsigned int lunghezza_funzione_max=0,
                                  unsigned int nthreads=0
            );
    ~GreenKubo2ComponentIonicFluid();
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void reset(unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    HeatFluxTs * puntatoreHeatFluxTs(){return je;}
    GreenKubo2ComponentIonicFluid & operator =(const GreenKubo2ComponentIonicFluid &);
private:
    bool scrivi_file;
    HeatFluxTs * je;
    ChargeFluxTs * jz;
    std::string log;
    Traiettoria *traiettoria;
    unsigned int ntimesteps,skip,lmax,leff,nthread;
};

#endif // GREENKUBO2COMPONENTIONICFLUID_H
