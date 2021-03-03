/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/

#ifndef HEATC_H
#define HEATC_H

#include "operazionisulista.h"
#include <vector>

class Traiettoria;

class HeatC: public OperazioniSuLista<HeatC, double>
{
public:
    HeatC(Traiettoria *t,
          double sigma,
          unsigned int nthreads=0,
          unsigned int skip=1);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) {return 1;}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    HeatC & operator = (const HeatC &);
    ~HeatC();
    std::vector<ssize_t> get_shape() const;
    std::vector<ssize_t> get_stride() const;

private:
    double sigma;
    unsigned int nthreads;
    unsigned int skip;
    Traiettoria *t;
    unsigned int ntimesteps;
    unsigned int tipi_atomi;
};

#endif // HEATC_H
