#ifndef CORRELATORESPAZIALE_H
#define CORRELATORESPAZIALE_H

#include "operazionisulista.h"

class Traiettoria;
typename fftw_complex;

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
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) {return 0;};
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    fftw_complex[3] s_fac_k(TFLOAT  k[3], unsigned int i_t );
    void s_fac_k(TFLOAT  k[3], unsigned int i_t,fftw_complex * out );
    CorrelatoreSpaziale<TFLOAT> & operator = (const CorrelatoreSpaziale<TFLOAT> &);
    ~CorrelatoreSpaziale();
private:
    fftw_complex * sfac;
    TFLOAT sigma2;
    unsigned int nthreads,nk,skip,tipi_atomi,ntimesteps;
    bool debug;
    Traiettoria *t;
};

#endif // CORRELATORESPAZIALE_H
