/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/

#include <stdexcept>
#include <algorithm>
#include "stdint.h"

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
    int condizioni_al_contorno[6];
    double scatola[6]; //xlo,xhi,ylo,yhi,zlo,zhi
    int dimensioni_riga_output;
    int nchunk;
};


// se si usa questo va sostituito ovunque al posto di Intestazione_timestep
struct Intestazione_timestep_triclinic {
    bigint timestep;
    bigint natoms;
    int triclinic;
    int condizioni_al_contorno[6];
    double scatola[6];
    double xy_xz_yz[3];
    int dimensioni_riga_output;
    int nchunk;
};


//questo va impostato al numero di dati per atomo
#define NDOUBLE_ATOMO 8
//dipende dal formato imposto alla simulazione LAMMPS!
//poi basterà convertire il puntatore char+ ad un puntatore Atomo*
//qui vanno messi in ordine i dati come li scrive LAMMPS
struct Atomo {
    double id;
    double tipo;
    double posizione[3];
    double velocita[3];
    static_assert (NDOUBLE_ATOMO==8, "You have to modify class Atomo if you change NDOUBLE_ATOMO" );
};
#pragma pack(pop)

//questo non potrà essere puntato direttamente ad una posizione nel file...
struct Intestazione_timestep_new : public Intestazione_timestep_triclinic {
    bigint magic_string_len;
    char * magic_string;
    char * columns = nullptr;
    char * unit_style = nullptr;
    int endian;
    int revision;
    int unit_style_len;
    char time_flag;
    double time;
    int columns_len;

    Intestazione_timestep_new(): magic_string{nullptr}, columns{nullptr},unit_style{nullptr} {
    }
    ~Intestazione_timestep_new() {
        delete [] magic_string;
        delete [] columns;
        delete [] unit_style;
    }
};

//#include <iostream>
class TimestepManager {
    enum class Type {
        Old, Old_triclinic, F2020
    };
    template <class T> char * read_and_advance(char * start, T * destination, int n =1) const {
        //destination = * (T*) start;
        std::copy_n((T*) start, n, destination);
        //std::cerr<<sizeof (T) <<std::endl;
        return start + sizeof (T)*n;
    }

public:

    TimestepManager() : intestazione_old{nullptr}, intestazione_old_tri{nullptr}, intestazione_new{nullptr} {}

    size_t read(char * begin, char * file_end){
        char * current_ptr=begin;
        if (begin+sizeof (Intestazione_timestep)>=file_end) {
            throw std::runtime_error("Error: end of file reached");
        }
        intestazione_old = (Intestazione_timestep *) begin;
        if (intestazione_old->timestep <0) {
            type=Type::F2020;
            delete intestazione_new;
            intestazione_new = new Intestazione_timestep_new;
            //read intestazione
            intestazione_new->magic_string_len=-intestazione_old->timestep;
            current_ptr+=sizeof(bigint);
            intestazione_new->magic_string = new char[intestazione_new->magic_string_len];
            current_ptr=read_and_advance(current_ptr,intestazione_new->magic_string,intestazione_new->magic_string_len);
#define READ_ADV(what) \
    current_ptr=read_and_advance(current_ptr,&(intestazione_new->what));
#define READ_ADV_A(what,n) \
    current_ptr=read_and_advance(current_ptr,intestazione_new->what,n);

            READ_ADV(endian);
            READ_ADV(revision);
            READ_ADV(timestep);
            READ_ADV(natoms);
            READ_ADV(triclinic);
            READ_ADV_A(condizioni_al_contorno,6);
            READ_ADV_A(scatola,6);
            if(intestazione_new->triclinic) {
                READ_ADV_A(xy_xz_yz,6);
            }
            READ_ADV(dimensioni_riga_output);
            if (intestazione_new->revision > 0x0001) {
                READ_ADV(unit_style_len);
                if (intestazione_new->unit_style_len>0){
                    intestazione_new->unit_style = new char[intestazione_new->unit_style_len];
                    READ_ADV_A(unit_style,intestazione_new->unit_style_len);
                }
                READ_ADV(time_flag);
                if (intestazione_new->time_flag){
                    READ_ADV(time);
                }
                READ_ADV(columns_len);
                intestazione_new->columns = new char[intestazione_new->columns_len];
                READ_ADV_A(columns,intestazione_new->columns_len);
            }
            READ_ADV(nchunk);

        } else {
            if (intestazione_old->triclinic) {
                type=Type::Old_triclinic;
                intestazione_old_tri = (Intestazione_timestep_triclinic*) begin;
                intestazione_old=nullptr;
                current_ptr+=sizeof (Intestazione_timestep_triclinic);
            } else {
                type= Type::Old;
                current_ptr+=sizeof (Intestazione_timestep);
            }
        }
        if (current_ptr>=file_end) {
            throw std::runtime_error("Error: end of file reached");
        }
        if (dimensioni_riga_output() != NDOUBLE_ATOMO) {
            throw std::runtime_error("ERROR: the binary format does not have "+std::to_string(NDOUBLE_ATOMO)+" numbers per atom");
        }
        return current_ptr-begin;
    }

#define select_var(what, kind)\
    kind what() const {\
        switch (type) {\
        case Type::Old:\
            return intestazione_old->what;\
        case Type::Old_triclinic:\
            return intestazione_old_tri->what;\
        case Type::F2020:\
            return intestazione_new->what;\
        }\
        throw std::runtime_error("Undefined binary type");\
    }
    select_var(nchunk,int)
    select_var(natoms,int)
    select_var(triclinic,int)
    select_var(timestep,bigint)
    select_var(scatola, double *)
    select_var(dimensioni_riga_output,int)

    ~TimestepManager () {
        delete intestazione_new;
    }


private:
    Type type;
    Intestazione_timestep * intestazione_old;
    Intestazione_timestep_triclinic * intestazione_old_tri;
    Intestazione_timestep_new * intestazione_new;
};

//ce ne sono nchunk
//questa non basta traslarla
struct Chunk {
    int n_atomi;
    Atomo * atomi; // questa cella di memoria andrà impostata

};


#endif // LAMMPS_STRUCT_H
