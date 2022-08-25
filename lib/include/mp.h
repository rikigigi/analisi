/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/




#ifndef MPI_H
#define MPI_H

#include "config.h"

#ifdef USE_MPI

#include <string>
#include <mpi.h>

template <class T,class TFLOAT> class VectorOp;

class Mp
{
public:
    static Mp & mpi(int *argc=0, char ***argv=0);
    int me();
    int size();
    bool ionode();
    std::string outname(std::string s);

     template <class T,class T2> void send_to_root(VectorOp<T,T2> * l) {
        MPI_Send( l->accesso_lista(),l->lunghezza(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }

     template <class T,class T2> void recv_root(VectorOp<T,T2> * l,int source){
    MPI_Recv(l->accesso_lista(),l->lunghezza(),MPI_DOUBLE,source,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

    ~Mp();

private:
    Mp(Mp const &) = delete;
    void operator=(Mp const&) = delete;
    Mp(int *argc, char ***argv);
    static bool initialized;
    int world_rank,world_size;

};

#endif
#endif // MPI_H
