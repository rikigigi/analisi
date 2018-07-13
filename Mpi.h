#ifndef MPI_H
#define MPI_H

template <class T,class TFLOAT=double> class OperazioniSuLista;

class Mpi
{
public:
    static Mpi & mpi(int *argc=0, char ***argv=0);
    int me();
    int size();
    bool ionode();

    template <class T,class TFLOAT> void send_to_root(OperazioniSuLista<T,TFLOAT> * l);
    template <class T,class TFLOAT> void recv_root(OperazioniSuLista<T,TFLOAT> * l);


private:
    Mpi(Mpi const &) = delete;
    void operator=(Mpi const&) = delete;
    Mpi(int *argc, char ***argv);
    static bool initialized;
    int world_rank,world_size;

};

#endif // MPI_H
