#ifndef CALCOLAMULTITHREAD_H
#define CALCOLAMULTITHREAD_H

#include <vector>
#include <sys/types.h>

/**
  * CRTP class that makes easy to implement a multithreaded calculation. The calculation is divided in blocks of contiguous timesteps, and then a function that calculates the quantity for every block is called.
  * The user must be careful about multithreading safety of the calc_single_th function (that should write on different memory address for different threads)
**/


#include <thread>

template <class T>
class CalcolaMultiThread
{
public:
    CalcolaMultiThread(unsigned int nthreads=0, unsigned int skip=0) : nthreads{nthreads},skip{skip},ntimesteps{0}
{
    if (nthreads==0) nthreads=1;
    if (skip==0) skip=1;
}

    void calcola(unsigned int primo){
        unsigned int npassith=ntimesteps/skip/nthreads;
        std::vector<std::thread> threads;
        for (unsigned int ith=0;ith<nthreads;++ith){
            threads.push_back(std::thread([&,ith](){
                unsigned int start=primo+npassith*ith*skip;
                unsigned int stop=start+npassith*skip;
                if (ith==nthreads-1)
                    stop=primo+ntimesteps;
                // calculate given start and stop timestep. note that & captures everything, user must take care of multithread safety of calc_single_th function
                static_cast<T*>(this)->calc_single_th(start,stop,primo,ith);

            }));
        }
        for (auto & t : threads){
            t.join();
        }
        threads.clear();
        static_cast<T*>(this)->join_data();

    }
    /*
     *  for array access in python
    */

    /*
     * number of elements for every dimension of the array
    */
    //std::vector<ssize_t> get_shape() const;
    /*
     * difference, in bytes, between two elements with the corresponding index incremented by one
    */
    //std::vector<ssize_t> get_stride() const ;
    /*
     * This is called by many threads at the same time in "calcola" function
    */
    // void calc_single_th(const unsigned int &start, const unsigned int & stop, const unsigned int & primo, const unsigned int & ith) ;
    /*This function is called when all threads finished their work*/
    // void join_data(){}
    /*
     * those must be implemented to interface with block averages stuff. Anyway, there are useful
    */
    //unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b)=0;
    /*
     * This function must allocate all the memory:  The argument is the number of timesteps to be calculated.
     * remember to set ntimesteps here, so the division of the work can be done correctly
    */
    // void reset(const unsigned int numeroTimestepsPerBlocco)=0;

protected:
    unsigned int nthreads,skip,ntimesteps;

private:

};

#endif // CALCOLAMULTITHREAD_H
