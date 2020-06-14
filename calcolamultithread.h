#ifndef CALCOLAMULTITHREAD_H
#define CALCOLAMULTITHREAD_H

#include <vector>
#include <sys/types.h>

/**
  * virtual class that makes easy to implement a multithreaded calculation. The calculation is divided in blocks of contiguous timesteps, and then a function that calculates the quantity for every block is called.
  * The user must be careful about multithreading safety of the calc_single_th function (that should write on different memory address for different threads)
**/

class CalcolaMultiThread
{
public:
    CalcolaMultiThread(unsigned int nthreads=0, unsigned int skip=0);
    void calcola(unsigned int primo);
    /*
     *  for array access in python
    */

    /*
     * number of elements for every dimension of the array
    */
    virtual std::vector<ssize_t> get_shape() const=0;
    /*
     * difference, in bytes, between two elements with the corresponding index incremented by one
    */
    virtual std::vector<ssize_t> get_stride() const =0;
    /*
     * This is called by many threads at the same time in "calcola" function
    */
    virtual void calc_single_th(const unsigned int &start, const unsigned int & stop, const unsigned int & primo, const unsigned int & ith) noexcept =0;
    /*
     * those must be implemented to interface with block averages stuff. Anyway, there are useful
    */
    virtual unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b)=0;
    /*
     * This function must allocate all the memory:  The argument is the number of timesteps to be calculated.
     * remember to set ntimesteps here, so the division of the work can be done correctly
    */
    virtual void reset(const unsigned int numeroTimestepsPerBlocco)=0;

protected:
    unsigned int nthreads,skip,ntimesteps;

private:

};

#endif // CALCOLAMULTITHREAD_H
