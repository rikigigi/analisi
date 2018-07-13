
#include "Mpi.h"

#include <mpi.h>
#include "operazionisulista.h"


bool Mpi::initialized=false;

Mpi & Mpi::mpi(int *argc, char ***argv) {

    static Mpi istanza(argc,argv);
    return istanza;

}

Mpi::Mpi(int *argc, char ***argv)
{

    int provided_thread=0;

    if (initialized) {
        std::cerr << "Errore: non posso inizializzare due volte l'ambiente MPI!\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    } else {
        MPI_Init_thread(argc,argv,MPI_THREAD_SERIALIZED,&provided_thread);
        MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
        MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    }

}


template <> void Mpi::send_to_root(OperazioniSuLista<double,double> * l) {
    MPI_Send( l->accesso_lista(),l->lunghezza(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
}

template <> void Mpi::recv_root(OperazioniSuLista<double,double> *l) {
    MPI_Recv(l->accesso_lista(),l->lunghezza(),MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

int Mpi::me(){
    return world_rank;
}
int Mpi::size(){
    return world_size;
}
bool Mpi::ionode() {
    if (world_rank==0)
        return true;
    return false;
}
