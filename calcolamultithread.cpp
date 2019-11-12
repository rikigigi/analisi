#include "calcolamultithread.h"
#include <thread>

CalcolaMultiThread::CalcolaMultiThread(unsigned int nthreads, unsigned int skip) : nthreads{nthreads},skip{skip},ntimesteps{0}
{
    if (nthreads==0) nthreads=1;
    if (skip==0) skip=1;
}


void CalcolaMultiThread::calcola(unsigned int primo) {
    unsigned int npassith=ntimesteps/skip/nthreads;
    std::vector<std::thread> threads;
    for (unsigned int ith=0;ith<nthreads;++ith){
        threads.push_back(std::thread([&,ith](){
            unsigned int start=primo+npassith*ith*skip;
            unsigned int stop=start+npassith*skip;
            if (ith==nthreads-1)
                stop=primo+ntimesteps;
            // calculate given start and stop timestep. note that & captures everything, user must take care of multithread safety of calc_single_th function
            calc_single_th(start,stop,primo,ith);

        }));
    }
    for (auto & t : threads){
        t.join();
    }
    threads.clear();

}
