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
#include "calcolamultithread.h"
#include "operazionisulista.h"
#include <vector>


class Traiettoria;


class CorrelatoreSpaziale : public CalcolaMultiThread<CorrelatoreSpaziale>, public OperazioniSuLista<CorrelatoreSpaziale,double>
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
    void calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int & ith) noexcept;
    void s_fac_k(const double  k[3], const unsigned int i_t,double * out ) const;
    int get_sfac_size()const {return size_sfac;}
    void print(std::ostream & out);
    using CalcolaMultiThread::operator=;
    ~CorrelatoreSpaziale();
    std::vector<ssize_t> get_shape() const;
    std::vector<ssize_t> get_stride() const;
    void join_data(){}

private:
    using OperazioniSuLista<CorrelatoreSpaziale,double>::lista;
    using OperazioniSuLista<CorrelatoreSpaziale,double>::lunghezza_lista;
    Traiettoria *t;
    double * sfac;
    std::vector< std::array<double,3> > klist;
    double sigma2;
    unsigned int size_sfac,size_k;
    unsigned int nk,tipi_atomi,ntimesteps;
    bool debug;
};

#endif // CORRELATORESPAZIALE_H
