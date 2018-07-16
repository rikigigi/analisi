
#include "mp.h"

#include <mpi.h>
#include "operazionisulista.h"
#include <iostream>

bool Mp::initialized=false;

Mp & Mp::mpi(int *argc, char ***argv) {

    static Mp istanza(argc,argv);
    return istanza;

}

Mp::Mp(int *argc, char ***argv)
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

Mp::~Mp() {
    MPI_Finalize();
}



int Mp::me(){
    return world_rank;
}
int Mp::size(){
    return world_size;
}
bool Mp::ionode() {
    if (world_rank==0)
        return true;
    return false;
}
