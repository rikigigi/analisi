/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "istogrammaatomiraggio.h"
#include "traiettoria.h"
#include<thread>

IstogrammaAtomiRaggio::IstogrammaAtomiRaggio(Traiettoria *t, double r, unsigned int skip,unsigned int nthreads) : r2(r*r), skip(skip),traiettoria(t),nthreads(nthreads),hist(0)
{
    if (IstogrammaAtomiRaggio::skip<1) IstogrammaAtomiRaggio::skip=1;
    if (IstogrammaAtomiRaggio::nthreads<1) IstogrammaAtomiRaggio::nthreads=1;
}

void IstogrammaAtomiRaggio::reset(const unsigned int numeroTimestepsPerBlocco) {
    numeroTimestepsBlocco=numeroTimestepsPerBlocco;
    ntypes=traiettoria->get_ntypes();
    natoms=traiettoria->get_natoms();
    delete [] hist;
    hist= new std::map<unsigned int, unsigned int> [ntypes];
}

void IstogrammaAtomiRaggio::calculate(unsigned int tstart) {
    unsigned int natomith=natoms/nthreads;
    std::vector<std::thread> threads;

    unsigned int tstop=tstart+numeroTimestepsBlocco;

    std::map<unsigned int, unsigned int> **histth= new std::map<unsigned int, unsigned int>* [nthreads];
    for (unsigned int i=0;i<nthreads;i++)
        histth[i]=new std::map<unsigned int, unsigned int> [ntypes];

    for (unsigned int itimestep=tstart;itimestep<tstop;itimestep+=skip){
        for (unsigned int ith=0;ith<nthreads;ith++) {
            threads.push_back(std::thread([&,ith](){
                unsigned int primo= natomith*ith;
                unsigned int ultimo= (ith != nthreads-1 )?natomith*(ith+1):natoms;

                unsigned int *cont = new unsigned int [ntypes];


                for (unsigned int iatom=primo;iatom<ultimo;iatom++) {
                    for (unsigned int j=0;j<ntypes;j++)
                        cont[j]=0;
                    for (unsigned int jatom=0;jatom<natoms;jatom++) {
                        if (   traiettoria->d2_minImage(iatom,jatom,itimestep,itimestep)    <   r2)
                            cont[traiettoria->get_type(jatom)]++;
                    }
                    for (unsigned int j=0;j<ntypes;j++)
                        ((histth[ith])[j])[cont[j]]++;
                }

            }));
        }
        for (unsigned int ith=0;ith<nthreads;ith++)
            threads[ith].join();
        threads.clear();
    }



    //somma i vari risultati dei threads
    for (unsigned int i=0;i<nthreads;i++){
        for (unsigned int j=0;j<ntypes;j++)
            for (auto it=histth[i][j].begin();it!=histth[i][j].end();++it){
                (hist[j])[it->first]+=it->second;
            }
    }

    //libera la memoria

    for (unsigned int i=0;i<nthreads;i++){
            delete [] histth[i];
    }
    delete [] histth;

}


