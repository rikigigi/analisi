/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef GOFRT_H
#define GOFRT_H

#include "operazionisulista.h"

class Traiettoria;

template <class TFLOAT>
class Gofrt : public OperazioniSuLista<Gofrt<TFLOAT>,TFLOAT>
{
public:
    Gofrt(Traiettoria *t,
          TFLOAT rmin,
          TFLOAT rmax,
          unsigned int nbin,
          unsigned int tmax=0,
          unsigned int nthreads=0,
          unsigned int skip=1,
          bool debug=false);
    ~Gofrt();
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    Gofrt<TFLOAT> & operator =(const Gofrt<TFLOAT> & destra);
    std::vector<ssize_t> get_shape();
    std::vector<ssize_t> get_stride();
    std::string get_columns_description() {return c_descr;};

private:
    using OperazioniSuLista<Gofrt<TFLOAT>,TFLOAT>::lista;
    using OperazioniSuLista<Gofrt<TFLOAT>,TFLOAT>::lunghezza_lista;
    TFLOAT rmin,rmax,dr;
    bool debug;
    Traiettoria * traiettoria;
    unsigned int nbin,ntimesteps,lmax,leff,skip,nthreads;
    TFLOAT * gofr(unsigned int ts, unsigned int itype=0, unsigned int r=0);
    std::string c_descr;
    unsigned int get_itype(unsigned int &t1, unsigned int &t2) const;

};

#endif // GOFRT_H
