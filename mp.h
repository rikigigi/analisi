#ifndef MPI_H
#define MPI_H

//#include "operazionisulista.h"

template <class T,class TFLOAT> class OperazioniSuLista;
#define OMPI_SKIP_MPICXX
#include <mpi.h>

class Mpi
{
public:
    static Mpi & mpi(int *argc=0, char ***argv=0);
    int me();
    int size();
    bool ionode();

     template <class T,class T2> void send_to_root(OperazioniSuLista<T,T2> * l) {
        MPI_Send( l->accesso_lista(),l->lunghezza(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }

     template <class T,class T2> void recv_root(OperazioniSuLista<T,T2> * l,int source){
    MPI_Recv(l->accesso_lista(),l->lunghezza(),MPI_DOUBLE,source,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

    ~Mpi();

private:
    Mpi(Mpi const &) = delete;
    void operator=(Mpi const&) = delete;
    Mpi(int *argc, char ***argv);
    static bool initialized;
    int world_rank,world_size;

};

#endif // MPI_H
