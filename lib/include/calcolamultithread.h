#ifndef CALCOLAMULTITHREAD_H
#define CALCOLAMULTITHREAD_H

#include <vector>
#include <sys/types.h>

/**
  * CRTP class that makes easy to implement a multithreaded calculation. The calculation is divided in blocks of contiguous timesteps, and then a function that calculates the quantity for every block is called.
  * The user must be careful about multithreading safety of the calc_single_th function (that should write on different memory address for different threads)
**/


#include <thread>

namespace CalcolaMultiThread_Flags {
constexpr int PARALLEL_SPLIT_AVERAGE = 0b00000001;
constexpr int PARALLEL_SPLIT_TIME =    0b00000010;
constexpr int PARALLEL_SPLIT_ATOM =    0b00000100;
constexpr int PARALLEL_LOOP_AVERAGE =  0b00010000;
constexpr int PARALLEL_LOOP_TIME =     0b00100000;
}


template <class T, int FLAGS = CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE>
class CalcolaMultiThread
{
public:
    CalcolaMultiThread(unsigned int nthreads=0, unsigned int skip=0,unsigned int natoms=0) : nthreads{nthreads},skip{skip},ntimesteps{0},natoms{natoms}
{
    if (nthreads==0) nthreads=1;
    if (skip==0) skip=1;

}

    std::pair<int,int> splitter(int ith, int nthreads) const {
        std::pair<int,int> res;
        res.first=ith*npassith;
        if (ith==nthreads-1) {
            res.second=end;
        } else {
            res.second=(ith+1)*npassith;
        }
        return res;
    }

    void init_split() {
        if constexpr (FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE){
            npassith=ntimesteps/skip/nthreads;
            end=ntimesteps;
        } else if (FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_TIME) {
            npassith=leff/nthreads;
            end=leff;
        } else if (FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_ATOM) {
            npassith=natoms/nthreads;
            end=natoms;
        } else {
            static_assert (! (FLAGS & (CalcolaMultiThread_Flags::PARALLEL_SPLIT_ATOM | CalcolaMultiThread_Flags::PARALLEL_SPLIT_TIME | CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE)),
                    "YOU MUST SPECIFY ONE OF CalcolaMultiThread_Flags" );
        }
        t0=0;t1=1;
        i0=0;i1=1;
        if constexpr (FLAGS & CalcolaMultiThread_Flags::PARALLEL_LOOP_AVERAGE) {
            i0=0;
            i1=ntimesteps;
        }
        if constexpr (FLAGS & CalcolaMultiThread_Flags::PARALLEL_LOOP_TIME) {
            t0=0;
            t1=leff;
        }
    }

    void calcola(unsigned int primo){
        init_split();
        std::vector<std::thread> threads;
        for (int t=t0;t<t1;t++){ // loop over time lags. Can be a loop over a single value if disabled
            for (int i=i0;i<i1;i+=skip){ //loop over trajectory. Can be a loop over a single value if disabled
                for (unsigned int ith=0;ith<nthreads;++ith){
                    threads.push_back(std::thread([&,ith](){
                        auto range=splitter(ith,nthreads);
                        // calculate given start and stop timestep. note that & captures everything, user must take care of multithread safety of calc_single_th function
                        // select a different signature of the calculation function depending on where the loops (if present) are parallelized
                        if constexpr (FLAGS & CalcolaMultiThread_Flags::PARALLEL_LOOP_AVERAGE & CalcolaMultiThread_Flags::PARALLEL_LOOP_TIME) {
                            static_cast<T*>(this)->calc_single_th(t,i,range.first+primo,range.second+primo,primo,ith);
                        } else if constexpr(FLAGS & CalcolaMultiThread_Flags::PARALLEL_LOOP_AVERAGE) {
                            static_cast<T*>(this)->calc_single_th(i,range.first+primo,range.second+primo,primo,ith);
                        } else if constexpr(FLAGS & CalcolaMultiThread_Flags::PARALLEL_LOOP_TIME) {
                            static_cast<T*>(this)->calc_single_th(t,range.first+primo,range.second+primo,primo,ith);
                        } else {
                            static_cast<T*>(this)->calc_single_th(range.first+primo,range.second+primo,primo,ith);
                        }

                    }));
                }
                for (auto & t : threads){
                    t.join();
                }
                threads.clear();
                static_cast<T*>(this)->join_data();
            }
        }

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
    unsigned int nthreads,skip,ntimesteps, leff;

private:
    unsigned int npassith,end, natoms;
    int t0,t1,i0,i1;
};

#endif // CALCOLAMULTITHREAD_H
