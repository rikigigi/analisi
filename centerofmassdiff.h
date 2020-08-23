#ifndef CENTEROFMASSDIFF_H
#define CENTEROFMASSDIFF_H

#include "calcolamultithread.h"
#include "operazionisulista.h"
#include <valarray>

class CenterOfMassDiff : public CalcolaMultiThread<CenterOfMassDiff>, public OperazioniSuLista<CenterOfMassDiff,double>
{
public:
    CenterOfMassDiff(Traiettoria *t, unsigned int nthreads=0, unsigned int skip=1);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) {return 0;}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int & ith) noexcept;
    std::vector<ssize_t> get_shape() const {return {ntimesteps/skip,ntype,3}; }
    std::vector<ssize_t> get_stride() const {return { static_cast<long> (sizeof (double)*3*ntype),sizeof (double)*3,sizeof (double)};}
    void set_starting_center(const std::valarray<double> & s) {if (s.size()!=starting_center.size()) throw std::runtime_error("wrong size of provided starting centers"); starting_center=s;}
    void set_zero_threshold(const double & z ){zero=z;}
    bool get_reset_error() {bool e=error; error=false; return e;}
    void join_data(){}
    ~CenterOfMassDiff();
private:
    Traiettoria *t;
    using OperazioniSuLista<CenterOfMassDiff,double>::lista;
    using OperazioniSuLista<CenterOfMassDiff,double>::lunghezza_lista;
    unsigned int lista_alloc;
    int ntype,maxiter;
    std::valarray<double> starting_center;
    double zero;
    bool error;

};

#endif // CENTEROFMASSDIFF_H
