#ifndef CORRELATORESPAZIALE_H
#define CORRELATORESPAZIALE_H

#include "operazionisulista.h"

#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif

class Traiettoria;


template <class TFLOAT>
class CorrelatoreSpaziale : public OperazioniSuLista<CorrelatoreSpaziale<TFLOAT>, TFLOAT>
{
public:
    CorrelatoreSpaziale(Traiettoria *t,
                        unsigned int n,
                        TFLOAT sigma2,
                        unsigned int nthreads=0,
                        unsigned int skip=1,
                        bool debug=false
            );
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) {return 0;}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    void s_fac_k(TFLOAT  k[3], unsigned int i_t,fftw_complex * out );
    CorrelatoreSpaziale<TFLOAT> & operator = (const CorrelatoreSpaziale<TFLOAT> &);
    ~CorrelatoreSpaziale();
private:
    using OperazioniSuLista<CorrelatoreSpaziale<TFLOAT>, TFLOAT>::lista;
    using OperazioniSuLista<CorrelatoreSpaziale<TFLOAT>, TFLOAT>::lunghezza_lista;

    fftw_complex * sfac;
    TFLOAT sigma2;
    unsigned int nthreads,nk,skip,tipi_atomi,ntimesteps;
    bool debug;
    Traiettoria *t;
    static fftw_plan fftw3;
};

#endif // CORRELATORESPAZIALE_H
