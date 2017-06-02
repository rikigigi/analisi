#ifndef HEATFLUXTS_H
#define HEATFLUXTS_H

#include <string>
#include "traiettoria.h"

class HeatFluxTs
{
public:
    HeatFluxTs(std::string filename,Traiettoria * t,unsigned int skip=1);
    ~HeatFluxTs();
    double * flux(unsigned int ts);
    unsigned int get_skip() {return skip;}
private:
    Traiettoria * traiettoria;
    unsigned int skip;
    double * heatflux;
    unsigned int size;
};

#endif // HEATFLUXTS_H
