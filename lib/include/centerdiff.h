#ifndef CENTERDIFF_H
#define CENTERDIFF_H

#include "calcolamultithread.h"
#include "operazionisulista.h"
#include <array>

class CenterDiff : public CalcolaMultiThread<CenterDiff>, public VectorOp<CenterDiff,double>
{
public:
    CenterDiff(Traiettoria *t, unsigned int nthreads=0, unsigned int skip=1, unsigned int nit=1,bool sum_first_two_and_ignore_vz=false,bool sum_1=false);
    unsigned int nExtraTimesteps(unsigned int n_b) {return 0;}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int & ith) noexcept;
    std::vector<ssize_t> get_shape() const {return {data_length/3/nit/3,nit,3,3}; }
    std::vector<ssize_t> get_stride() const {return { static_cast<long> (sizeof (double)*3*3*nit),sizeof (double)*3*3,sizeof (double)*3,sizeof (double)};}
    void set_starting_center(const std::array<double,3*3> & s) {starting_center=s;}
    void join_data(){}
    ~CenterDiff();
private:
    Traiettoria *t;
    using VectorOp<CenterDiff,double>::vdata;
    using VectorOp<CenterDiff,double>::data_length;
    unsigned int nit,lista_alloc;
    bool sum_first_two_and_ignore_vz,sum_1;
    std::array<double,3*3> starting_center;

};

#endif // CENTERDIFF_H
