
#include "mp.h"

#include <mpi.h>
#include "operazionisulista.h"
#include <iostream>

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
        if (!ionode()) { //suppress all output from other nodes
            std::cout.setstate(std::ios_base::failbit);
            std::cerr.setstate(std::ios_base::failbit);
            std::clog.setstate(std::ios_base::failbit);
        }
    }

}

Mpi::~Mpi() {
    MPI_Finalize();
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
