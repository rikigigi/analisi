/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef SPETTROVIBRAZIONALE_H
#define SPETTROVIBRAZIONALE_H

#include "config.h"

#include "operazionisulista.h"
#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif
#include <vector>

class Traiettoria;


template <class T>
class SpettroVibrazionale : public VectorOp<SpettroVibrazionale<T> >
{
public:
    unsigned int nExtraTimesteps(unsigned int n_b);
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calculate(unsigned int primo);
    SpettroVibrazionale(T* t,bool dump=false);
    ~SpettroVibrazionale();
    std::vector<ssize_t> get_shape() const { return {static_cast<ssize_t> (tipi_atomi),static_cast<ssize_t>(size/2+1),static_cast<ssize_t>(3)} ; }
    std::vector<ssize_t> get_stride() const { return { static_cast<ssize_t> (3*(size/2+1)*sizeof(double)),static_cast<ssize_t>(3*sizeof (double)),static_cast<ssize_t>(sizeof(double))};}
    static void deallocate_plan(); // da chiamare alla fine del programma!
    double spettro(unsigned int frequenza, unsigned int dim, unsigned int tipo_atomo);
    SpettroVibrazionale<T> & operator = (const SpettroVibrazionale<T> &);


private:
    T * traiettoria;
    unsigned int size;
    int tipi_atomi;
    using VectorOp<SpettroVibrazionale<T> >::vdata;
    using VectorOp<SpettroVibrazionale<T> >::data_length;
    unsigned int trasformata_size;
    fftw_complex * trasformata;
    static fftw_plan fftw3;
    static unsigned int fplan_natoms,fplan_size;
    bool dump;


};

//per fare anche le varie medie a blocchi
//template class BlockAverage<SpettroVibrazionale>;

#endif // SPETTROVIBRAZIONALE_H
