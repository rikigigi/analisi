/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef ISTOGRAMMAVELOCITA_H
#define ISTOGRAMMAVELOCITA_H

#include "operazionisulista.h"
#include "traiettoria.h"

class IstogrammaVelocita : public VectorOp<IstogrammaVelocita>
{
public:
    IstogrammaVelocita(Traiettoria *t,unsigned int nbins,double vminmax_);
    unsigned int nExtraTimesteps(unsigned int n_b);
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    IstogrammaVelocita & operator =(const IstogrammaVelocita & destra);
private:
    Traiettoria *traiettoria;
    unsigned int ntimestep,bins;
    double vminmax;

};

#endif // ISTOGRAMMAVELOCITA_H
