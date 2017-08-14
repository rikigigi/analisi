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
    double * temp(unsigned int ts);
    unsigned int get_skip() {return skip;}
    double get_L(){return L;}
private:
    Traiettoria * traiettoria;
    unsigned int skip;
    double * heatflux;
    double * T,L;
    unsigned int size;
};

#endif // HEATFLUXTS_H
