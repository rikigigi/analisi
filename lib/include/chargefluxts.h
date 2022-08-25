/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/


#ifndef CHARGEFLUXTS_H
#define CHARGEFLUXTS_H

#include "trajectory.h"

class ChargeFluxTs
{
public:
    ChargeFluxTs(Trajectory * t);
    ~ChargeFluxTs();
    double *J_z(const unsigned int & timestep);
    void calculate(unsigned int primo);
    void reset(unsigned int numeroTimestepPerBlocco);
private:
    Trajectory * traiettoria;
    double *J;
    unsigned int timestep_finestra,timestep_corrente;
    bool calcolato;
};

#endif // CHARGEFLUXTS_H
