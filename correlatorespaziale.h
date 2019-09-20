/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef CORRELATORESPAZIALE_H
#define CORRELATORESPAZIALE_H

#include "operazionisulista.h"

#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif

class Traiettoria;


class CorrelatoreSpaziale : public OperazioniSuLista<CorrelatoreSpaziale, double>
{
public:
    CorrelatoreSpaziale(Traiettoria *t,
                        unsigned int n,
                        double sigma2,
                        unsigned int nthreads=0,
                        unsigned int skip=1,
                        bool debug=false
            );
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) {return 0;}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    void s_fac_k(double  k[3], unsigned int i_t,fftw_complex * out );
    double corr(unsigned int rx, unsigned int ry, unsigned int rz, unsigned int itype, unsigned int idim=0);
    CorrelatoreSpaziale & operator = (const CorrelatoreSpaziale &);
    ~CorrelatoreSpaziale();
private:
    using OperazioniSuLista<CorrelatoreSpaziale, double>::lista;
    using OperazioniSuLista<CorrelatoreSpaziale, double>::lunghezza_lista;

    fftw_complex * sfac;
    double sigma2;
    unsigned int size;
    unsigned int nthreads,nk,skip,tipi_atomi,ntimesteps;
    bool debug;
    Traiettoria *t;
    static fftw_plan fftw3;
};

#endif // CORRELATORESPAZIALE_H
