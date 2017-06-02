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



#ifndef ISTOGRAMMAVELOCITA_H
#define ISTOGRAMMAVELOCITA_H

#include "operazionisulista.h"
#include "traiettoria.h"

class IstogrammaVelocita : public OperazioniSuLista<IstogrammaVelocita>
{
public:
    IstogrammaVelocita(Traiettoria *t,unsigned int nbins,double vminmax_);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
private:
    Traiettoria *traiettoria;
    unsigned int ntimestep,bins;
    double vminmax;

};

#endif // ISTOGRAMMAVELOCITA_H
