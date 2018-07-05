#ifndef LAMMPS_STRUCT_H
#define LAMMPS_STRUCT_H

typedef int tagint;
typedef int64_t bigint;

//poi basterà convertire il puntatore char a un puntatore Intestazione_timestep*

//pragma pack(1) dice al compilatore (dalla documentazione funziona con gcc e il compilatore microsoft) di non inserire byte di allineamento
#pragma pack(push,1)
struct Intestazione_timestep {
    bigint timestep;
    bigint natoms;
    int triclinic;
    int condizioni_al_contorno[3][2];
    double scatola[6]; //xlo,xhi,ylo,yhi,zlo,zhi
    int dimensioni_riga_output;
    int nchunk;
};


// se si usa questo va sostituito ovunque al posto di Intestazione_timestep
struct Intestazione_timestep_triclinic {
    bigint timestep;
    bigint natoms;
    int triclinic;
    int condizioni_al_contorno[3][2];
    double scatola[6];
    double xy_xz_yz[3];
    int dimensioni_riga_output;
    int nchunk;
};

//dipende dal formato imposto alla simulazione LAMMPS!
//poi basterà convertire il puntatore char+ ad un puntatore Atomo*
//qui vanno messi in ordine i dati come li scrive LAMMPS
struct Atomo {
    double id;
    double tipo;
    double posizione[3];
    double velocita[3];
};
//questo va impostato al numero di dati per atomo
#define NDOUBLE_ATOMO 8
#pragma pack(pop)

//ce ne sono nchunk
//questa non basta traslarla
struct Chunk {
    int n_atomi;
    Atomo * atomi; // questa cella di memoria andrà impostata

};


#endif // LAMMPS_STRUCT_H
