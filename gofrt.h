#ifndef GOFRT_H
#define GOFRT_H

#include "operazionisulista.h"

class Traiettoria;

class Gofrt : public OperazioniSuLista<Gofrt>
{
public:
    Gofrt(Traiettoria *t,
          double rmin,
          double rmax,
          unsigned int nbin,
          unsigned int tmax=0,
          unsigned int nthreads=0,
          unsigned int skip=1,
          bool debug=false);
    ~Gofrt();
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    Gofrt & operator =(const Gofrt & destra);

private:
    double rmin,rmax;
    bool debug;
    Traiettoria * traiettoria;
    unsigned int nbin,ntimesteps,lmax,leff,skip,nthreads;
    double * gofr(unsigned int ts,unsigned int itype);

};

#endif // GOFRT_H
