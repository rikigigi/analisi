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


#ifndef CHARGEFLUXTS_H
#define CHARGEFLUXTS_H

#include "traiettoria.h"

class ChargeFluxTs
{
public:
    ChargeFluxTs(Traiettoria * t);
    double *J_z(const unsigned int & timestep);
    void calcola(unsigned int primo);
    void reset(unsigned int numeroTimestepPerBlocco);
private:
    Traiettoria * traiettoria;
    double *J;
    unsigned int timestep_finestra,timestep_corrente;
    bool calcolato;
};

#endif // CHARGEFLUXTS_H
