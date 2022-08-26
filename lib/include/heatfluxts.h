/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef HEATFLUXTS_H
#define HEATFLUXTS_H

#include <string>
#include "trajectory.h"

class HeatFluxTs
{
public:
    HeatFluxTs(std::string filename,Trajectory * t,unsigned int skip=1);
    ~HeatFluxTs();
    double * flux(unsigned int ts);
    double * temp(unsigned int ts);
    unsigned int get_skip() {return skip;}
    double get_L(){return L;}
private:
    Trajectory * traiettoria;
    unsigned int skip;
    double * heatflux;
    double * T,L;
    unsigned int size;
};

#endif // HEATFLUXTS_H
