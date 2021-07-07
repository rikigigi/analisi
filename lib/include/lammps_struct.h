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
#include <cstring>

#ifndef LAMMPS_STRUCT_H
#define LAMMPS_STRUCT_H

typedef int tagint;
typedef int64_t bigint;



template <class T> char * read_and_advance(char * start, T * destination, int n =1) {
    std::copy_n((T*) start, n, destination);
    return start + sizeof (T)*n;
}


struct Intestazione_timestep {
    bigint timestep;
    bigint natoms;
    int triclinic;
    int condizioni_al_contorno[6];
    double scatola[6]; //xlo,xhi,ylo,yhi,zlo,zhi
    int dimensioni_riga_output;
    int nchunk;
    char * read(char * file) {
        file=read_and_advance(file, &timestep);
        file=read_and_advance(file, &natoms);
        file=read_and_advance(file,&triclinic);
        file=read_and_advance(file, condizioni_al_contorno,6);
        file=read_and_advance(file,scatola,6);
        file=read_and_advance(file,&dimensioni_riga_output);
        file=read_and_advance(file,&nchunk);
        return file;
    }
    static int get_triclinic(char * file) {
        int triclinic;
        read_and_advance(file+sizeof (timestep)+sizeof (natoms),&triclinic);
        return triclinic;
    }
};


struct Intestazione_timestep_triclinic {
    bigint timestep;
    bigint natoms;
    int triclinic;
    int condizioni_al_contorno[6];
    double scatola[6];
    double xy_xz_yz[3];
    int dimensioni_riga_output;
    int nchunk;
    char * read(char * file) {
        file=read_and_advance(file, &timestep);
        file=read_and_advance(file, &natoms);
        file=read_and_advance(file,&triclinic);
        file=read_and_advance(file, condizioni_al_contorno,6);
        file=read_and_advance(file,scatola,6);
        file=read_and_advance(file,xy_xz_yz,3);
        file=read_and_advance(file,&dimensioni_riga_output);
        file=read_and_advance(file,&nchunk);
        return file;
    }
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
    static char * read_id_tipo(char * begin, double & id_, double & tipo_) {
        begin = read_and_advance(begin,&id_);
        begin = read_and_advance(begin,&tipo_);
        return begin;
    }
    static char * read_pos_vel(char * begin, double * pos, double * vel) {
        begin = read_and_advance(begin,pos,3);
        begin = read_and_advance(begin,vel,3);
        return begin;
    }
};
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

    char * read(char * current_ptr) {
        current_ptr=read_and_advance(current_ptr, &magic_string_len);
        magic_string_len=-magic_string_len;
        delete [] magic_string;
        magic_string = new char[magic_string_len];

        current_ptr=read_and_advance(current_ptr,magic_string,magic_string_len);
#define READ_ADV(what) \
current_ptr=read_and_advance(current_ptr,&what);
#define READ_ADV_A(what,n) \
current_ptr=read_and_advance(current_ptr,what,n);

        READ_ADV(endian);
        READ_ADV(revision);
        READ_ADV(timestep);
        READ_ADV(natoms);
        READ_ADV(triclinic);
        READ_ADV_A(condizioni_al_contorno,6);
        READ_ADV_A(scatola,6);
        if(triclinic) {
            READ_ADV_A(xy_xz_yz,6);
        }
        READ_ADV(dimensioni_riga_output);
        if (revision > 0x0001) {
            READ_ADV(unit_style_len);
            if (unit_style_len>0){
                delete [] unit_style;
                unit_style = new char[unit_style_len];
                READ_ADV_A(unit_style,unit_style_len);
            }
            READ_ADV(time_flag);
            if (time_flag){
                READ_ADV(time);
            }
            READ_ADV(columns_len);
            delete columns;
            columns = new char[columns_len];
            READ_ADV_A(columns,columns_len);
        }
        READ_ADV(nchunk);
        return current_ptr;
    }
};

//#include <iostream>
class TimestepManager {
    enum class Type {
        Old, Old_triclinic, F2020
    };

public:

    TimestepManager() {}

    size_t read(char * begin, char * file_end){
        char * current_ptr=begin;
        if (begin+sizeof (Intestazione_timestep)>=file_end) {
            throw std::runtime_error("Error: end of file reached");
        }
        bigint timestep;
        read_and_advance(begin, &timestep);
        if (timestep <0) {
            type=Type::F2020;
            current_ptr=intestazione_new.read(begin);
        } else {
            if (Intestazione_timestep::get_triclinic(begin)) {
                type=Type::Old_triclinic;
                current_ptr = intestazione_old_tri.read(begin);

            } else {
                type= Type::Old;
                current_ptr=intestazione_old.read(begin);
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
    kind what() {\
        switch (type) {\
        case Type::Old:\
            return intestazione_old.what;\
        case Type::Old_triclinic:\
            return intestazione_old_tri.what;\
        case Type::F2020:\
            return intestazione_new.what;\
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
    }


private:
    Type type;
    Intestazione_timestep intestazione_old;
    Intestazione_timestep_triclinic intestazione_old_tri;
    Intestazione_timestep_new intestazione_new;
};

struct Chunk {
    int n_atomi;
    char * atomi; // questa cella di memoria andrà impostata

};


#endif // LAMMPS_STRUCT_H
