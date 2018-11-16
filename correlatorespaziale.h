#ifndef CORRELATORESPAZIALE_H
#define CORRELATORESPAZIALE_H

#include "operazionisulista.h"

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
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) {return 0;};
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    CorrelatoreSpaziale<TFLOAT> & operator = (const CorrelatoreSpaziale<TFLOAT> &);
    ~CorrelatoreSpaziale();
private:
    TFLOAT * e_of_r;
    TFLOAT sigma2;
    unsigned int nthreads,n,skip    ;
    bool debug;
};

#endif // CORRELATORESPAZIALE_H
