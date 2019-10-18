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
#include <vector>

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
                        std::vector< std::array<double,3> >  k,
                        double sigma2,
                        unsigned int nthreads=0,
                        unsigned int skip=1,
                        bool debug=false
            );
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) {return 1;}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    void s_fac_k(double  k[3], unsigned int i_t,double * out );
    int get_sfac_size(){return size_sfac;}
    void print(std::ostream & out);
    CorrelatoreSpaziale & operator = (const CorrelatoreSpaziale &);
    ~CorrelatoreSpaziale();
    std::vector<ssize_t> get_shape();
    std::vector<ssize_t> get_stride();

private:
    using OperazioniSuLista<CorrelatoreSpaziale, double>::lista;
    using OperazioniSuLista<CorrelatoreSpaziale, double>::lunghezza_lista;

    double * sfac;
    std::vector< std::array<double,3> > klist;
    double sigma2;
    unsigned int size_sfac,size_k;
    unsigned int nthreads,nk,skip,tipi_atomi,ntimesteps;
    bool debug;
    Traiettoria *t;
};

#endif // CORRELATORESPAZIALE_H
