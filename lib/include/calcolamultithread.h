#ifndef CALCOLAMULTITHREAD_H
#define CALCOLAMULTITHREAD_H

#include <vector>
#include <sys/types.h>
#include "compiler_types.h"

#include <iostream>
#include <exception>
/**
  * CRTP class that makes easy to implement a multithreaded calculation. The calculation is divided in blocks of contiguous timesteps, and then a function that calculates the quantity for every block is called.
  * The user must be careful about multithreading safety of the calc_single_th function (that should write on different memory address for different threads)
**/


#include <thread>

namespace CalcolaMultiThread_Flags {
//mutually exclusive SPLITS:
constexpr int PARALLEL_SPLIT_AVERAGE =  0b00000001; // assign to the single thread worker different averages range
constexpr int PARALLEL_SPLIT_TIME =     0b00000010; // assigne to the single thread worker different time differences ranges
constexpr int PARALLEL_SPLIT_ATOM =     0b00000100; // assign to the single thread worker different atom ranges
//loop to perform in serial: this will use different signatures of the single thread worker
constexpr int SERIAL_LOOP_AVERAGE =   0b00010000; // do the loop for computing the average outside the parallelized region; pass time average index in first position
constexpr int SERIAL_LOOP_TIME =      0b00100000; // do the loop over the times outside the parallelized region; pass time index in first position
//if both AVERAGE and TIME are performed in serial the code will pass first the time and then the average index
//if no SERIAL_LOOP_* flag is present the code will call the worker with the argument:
// rangeA, rangeB, first index of the block, thread id, where rangeA and rangeB are computed to split the work across the thread according to the PARALLEL_SPLIT_* option choosen

//options to call additional routines in some part of the computation
constexpr int CALL_INNER_JOIN_DATA =    0b01000000;
constexpr int CALL_DEBUG_ROUTINE =      0b10000000;
constexpr int CALL_CALC_INIT    =      0b100000000;
}


template <class T, int FLAGS_T = CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE | CalcolaMultiThread_Flags::CALL_INNER_JOIN_DATA>
class CalcolaMultiThread
{
public:
    CalcolaMultiThread(const ssize_t nthreads=0, const ssize_t skip=0, const size_t natoms=0, const ssize_t every=0) : nthreads{nthreads},skip{skip},ntimesteps{0},natoms{natoms},every{every}
{
    if (nthreads==0) CalcolaMultiThread::nthreads=1;
    if (skip==0) CalcolaMultiThread::skip=1;
    if (every==0) CalcolaMultiThread::every=1 ; 

}
    static constexpr int FLAGS=FLAGS_T;

    std::pair<size_t,size_t> splitter(size_t ith, size_t primo) const {
        std::pair<size_t,size_t> res;
        res.first=ith*npassith;
        if (ith==nthreads-1) {
            res.second=end;
            if constexpr ( !!(FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE)){
                //align the last timestep so it is a multiple of skip and the total number of points is ntimestep/skip
                res.second = end - end%skip;
            }
            if constexpr ( !!(FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_TIME)) {
                //align the last timestep so it is a multiple of every and the total number of points is leff/every
                res.second= end - end%every;
            }
        } else {
            res.second=(ith+1)*npassith;
        }
        if constexpr ( !!(FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE)) {
            if (skip > 1 && res.first % skip > 0) { 
                //align the first timestep so it is a multiple of skip
                res.first = res.first - res.first % skip + skip;
            }
        }
        if constexpr ( !!(FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_TIME)) {
            if (every > 1 && res.first % every > 0 ) {
                res.first = res.first - res.first % every + every ;
            }
        }
        return res;
    }

    void init_split() {
        if constexpr ( !! (FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE)){
            npassith=ntimesteps/nthreads;
            end=ntimesteps;
        } else if (FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_TIME) {
            npassith=leff/nthreads;
            end=leff;
        } else if (FLAGS & CalcolaMultiThread_Flags::PARALLEL_SPLIT_ATOM) {
            npassith=natoms/nthreads;
            end=natoms;
        } else {
            static_assert (FLAGS & (CalcolaMultiThread_Flags::PARALLEL_SPLIT_ATOM | CalcolaMultiThread_Flags::PARALLEL_SPLIT_TIME | CalcolaMultiThread_Flags::PARALLEL_SPLIT_AVERAGE),
                    "YOU MUST SPECIFY ONE OF CalcolaMultiThread_Flags" );
        }
        t0=0;t1=1;
        i0=0;i1=1;
        if constexpr ( !!(FLAGS & CalcolaMultiThread_Flags::SERIAL_LOOP_AVERAGE)) {
            i0=0;
            i1=ntimesteps;
        }
        if constexpr ( !!(FLAGS & CalcolaMultiThread_Flags::SERIAL_LOOP_TIME)) {
            t0=0;
            t1=leff;
        }
    }

    void calcola(size_t primo){
        init_split();
        if constexpr (!!(FLAGS & CalcolaMultiThread_Flags::CALL_CALC_INIT)) {
            static_cast<T*>(this)->calc_init(primo);
        }
        std::exception_ptr * thread_exception = new std::exception_ptr[nthreads];
        for (size_t i=0;i<nthreads;++i) thread_exception[i]=nullptr;
        std::vector<std::thread> threads;
        for (ssize_t t=t0;t<t1;t+=every){ // loop over time lags. Can be a loop over a single value if disabled
            for (ssize_t i=i0;i<i1;i+=skip){ //loop over trajectory. Can be a loop over a single value if disabled
                for (unsigned int ith=0;ith<nthreads;++ith){ //this is a loop on TIME/AVERAGE INDEX/ATOM INDEX depending on what PARALLEL_SPLIT_* you choose
                    threads.push_back(std::thread([&,ith,t,i](){
                        try {
                            auto range=splitter(ith,primo);
                            // calculate given start and stop timestep. note that & captures everything, user must take care of multithread safety of calc_single_th function
                            // select a different signature of the calculation function depending on where the loops (if present) are parallelized
                            // rangeA and rangeB are passed to each thread and their meaning depends on what PARALLEL_SPLIT_* option you choose
                            if constexpr ( !! (FLAGS & CalcolaMultiThread_Flags::SERIAL_LOOP_AVERAGE && (FLAGS & CalcolaMultiThread_Flags::SERIAL_LOOP_TIME))) {
                                static_cast<T*>(this)->calc_single_th(t,i,range.first,range.second,primo,ith); // time, average, rangeA, rangeB, primo, thread id
                            } else if constexpr( !! (FLAGS & CalcolaMultiThread_Flags::SERIAL_LOOP_AVERAGE) ) {
                                static_cast<T*>(this)->calc_single_th(i,range.first,range.second,primo,ith); // average, rangeA, rangeB, primo, thread id
                            } else if constexpr( !! (FLAGS & CalcolaMultiThread_Flags::SERIAL_LOOP_TIME) ) {
                                static_cast<T*>(this)->calc_single_th(t,range.first,range.second,primo,ith); // time, rangeA, rangeB, primo, thread id
                            } else {
                                static_cast<T*>(this)->calc_single_th(range.first,range.second,primo,ith); // rangeA, rangeB, primo, thread id
                            }
                        } catch (...) {
                            thread_exception[ith] = std::current_exception();
                        }
                    }));
//                    std::cerr << "thread " << & threads.back() << " started" <<std::endl;
                }
                for (auto & t : threads){
                    t.join();
                }
                for (size_t i=0;i<nthreads;++i){
                    if (thread_exception[i]) {
                         std::cerr << "thread "<< i <<" exited with an exception"<<std::endl;
                    }
                }
                for (size_t i=0;i<nthreads;++i){
                    if (thread_exception[i]) {
                         std::rethrow_exception(thread_exception[i]);
                    }
                }
                threads.clear();
                if constexpr (!!(FLAGS & CalcolaMultiThread_Flags::CALL_INNER_JOIN_DATA))
                    static_cast<T*>(this)->join_data();
            }
        }

        if constexpr (!!(FLAGS & CalcolaMultiThread_Flags::CALL_DEBUG_ROUTINE)) {
            static_cast<T*>(this)->calc_end();
        }
        delete [] thread_exception;

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
    ssize_t nthreads,skip,ntimesteps, every , leff;

private:
    size_t npassith,end, natoms;
    ssize_t t0,t1,i0,i1;
};

#endif // CALCOLAMULTITHREAD_H
