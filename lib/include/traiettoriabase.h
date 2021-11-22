#ifndef TRAIETTORIABASE_H
#define TRAIETTORIABASE_H

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <map>
#include <algorithm>
#include "compiler_types.h"

#define DECL_ARG(a,b) a b
#define NAME_ARG(a,b) b
#define DECL_CALL_BASE_0(ret, fname) ret fname () { return static_cast<T*>(this)-> fname ( ); }
#define DECL_CALL_BASE_1(ret, fname, arg1) ret fname ( DECL_ARG arg1) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1 ); }
#define DECL_CALL_BASE_2(ret, fname, arg1, arg2) ret fname ( DECL_ARG arg1, DECL_ARG arg2) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2 ); }
#define DECL_CALL_BASE_3(ret, fname, arg1, arg2, arg3) ret fname ( DECL_ARG arg1, DECL_ARG arg2, DECL_ARG arg3) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2, NAME_ARG arg3 ); }
#define DECL_CALL_BASE_4(ret, fname, arg1, arg2, arg3, arg4) ret fname ( DECL_ARG arg1, DECL_ARG arg2, DECL_ARG arg3, DECL_ARG arg4) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2, NAME_ARG arg3, NAME_ARG arg4 ); }
#define DECL_CALL_BASE_5(ret, fname, arg1, arg2, arg3, arg4, arg5) ret fname ( DECL_ARG arg1, DECL_ARG arg2, DECL_ARG arg3, DECL_ARG arg4, DECL_ARG arg5) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2, NAME_ARG arg3, NAME_ARG arg4, NAME_ARG arg5 ); }

template <class T>
class TraiettoriaBase {

public:

    TraiettoriaBase() : ntypes{0},n_timesteps{0},natoms{0},min_type{0},max_type{0},
    wrap_pbc{true}, buffer_posizioni{nullptr}, buffer_velocita{nullptr},
    buffer_scatola{nullptr}, buffer_posizioni_cm{nullptr},
    buffer_velocita_cm{nullptr}, masse{nullptr}, cariche{nullptr},
    buffer_tipi{nullptr},buffer_tipi_id{nullptr}, serve_pos{true} {}

    //double * posizioni (const int & timestep, const int & atomo) { return static_cast<T*>(this)->posizioni(timestep,atomo);}
    DECL_CALL_BASE_2(double *, posizioni, (const int &, timestep), (const int &, atomo))
    DECL_CALL_BASE_2(double *, velocita,(const int &, timestep), (const int &, atomo) )
    DECL_CALL_BASE_1(double *, scatola, (const int &, timestep))
    DECL_CALL_BASE_2(double *, posizioni_cm,(const int &, timestep), (const int &, tipo))
    DECL_CALL_BASE_2(double *, velocita_cm,(const int &, timestep), (const int &, tipo))
    DECL_CALL_BASE_0(double *,scatola_last)
    DECL_CALL_BASE_0(int, ntimesteps_loaded)
    std::vector<unsigned int> get_types(){
        get_ntypes();
        return types;

    }
    unsigned int get_type(const unsigned int &atomo){
        if (atomo < natoms) {
            return buffer_tipi_id[atomo];
        } else {
            std::cerr << "Errore: tipo atomo fuori dal range! ("<< atomo <<" su un massimo di "<<natoms<<")\n";
            abort();
            return 0;
        }
    }
    enum Errori {non_inizializzato=0,oltre_fine_file=2,Ok=1};
    Errori imposta_dimensione_finestra_accesso(const size_t & timesteps){std::cerr << "Warning: doing nothing (not reading in blocks)"<<std::endl; return Errori::Ok;}
    Errori imposta_inizio_accesso(const size_t & timesteps){std::cerr << "Warning: doing nothing (not reading in blocks)"<<std::endl;return Errori::Ok;}
    //void index_all();

    //this can produce multiple bugs: what happens if I call it in the wrong moment?
    void set_pbc_wrap(bool p){
        wrap_pbc=p;
    }


    int get_ntypes (){
        if (ntypes==0) {
            types.clear();
            min_type=buffer_tipi[0];
            max_type=buffer_tipi[0];
            bool *duplicati = new bool[natoms];
            for (size_t i=0;i<natoms;i++)
                duplicati[i]=false;
            for (size_t i=0;i<natoms;i++) {
                if (!duplicati[i]) {
                    if (buffer_tipi[i]>max_type)
                        max_type=buffer_tipi[i];
                    if (buffer_tipi[i]<min_type)
                        min_type=buffer_tipi[i];
                    for (size_t j=i+1;j<natoms;j++){
                        if (buffer_tipi[j]==buffer_tipi[i]){
                            duplicati[j]=true;
                        }
                    }
                    types.push_back(buffer_tipi[i]);
                    ntypes++;
                }
            }
            std::sort(types.begin(),types.end());
            type_map.clear();
            for (unsigned int i=0;i<types.size(); i++){
                type_map[types[i]]=i;
            }
            delete [] duplicati;
            masse = new double [ntypes];
            cariche = new double [ntypes];

            for (size_t i=0;i<natoms;i++) {
                buffer_tipi_id[i]=type_map.at(buffer_tipi[i]);
            }
        }
        return ntypes;
    }
    double * posizioni_inizio(){return buffer_posizioni;}
    double * velocita_inizio(){return buffer_velocita;}
    int get_type_min() {return min_type;}
    int get_type_max() {return max_type;}
    size_t get_natoms()const {return natoms;}
    size_t get_ntimesteps() const{return n_timesteps;}
    double get_mass(unsigned int i) {if (i<get_ntypes()) return masse[i]; throw std::runtime_error("Cannot get mass for a type that does not exists!\n");}
    void set_mass(unsigned int i,double m) {if (i<get_ntypes()) masse[i]=m;}
    void set_charge(unsigned int i, double c){if (i<get_ntypes()) cariche[i]=c;}
    double get_charge(unsigned int  i){if (i<get_ntypes()) return cariche[i]; std::cerr<< "Cannot get charge for a type that does not exists!\n";abort(); return 0.0;}
    double d2_minImage(size_t i,size_t j, size_t itimestep,double *l){
        return d2_minImage(i,j,itimestep,itimestep,l);
    }
    double d2_minImage(size_t i,size_t j, size_t itimestep, size_t jtimestep, double *l){
        double x[3];
        return d2_minImage(i,j,itimestep,jtimestep,l,x);
    }
    double d2_minImage(size_t i,size_t j, size_t itimestep, size_t jtimestep,double *l,double *x){
        double d2=0.0;
        double *xi=posizioni(itimestep,i);
        double *xj=posizioni(jtimestep,j);
        for (unsigned int idim=0;idim<3;idim++) {
            x[idim]=xi[idim]-xj[idim];
            int k= (x[idim]/l[idim] + ((x[idim]>=0.0) ? 0.5 : -0.5));
            x[idim]-=l[idim]*k;
            /*
            if (x[idim] >   l[idim] * 0.5) x[idim] = x[idim] - l[idim];
            if (x[idim] <= -l[idim] * 0.5) x[idim] = x[idim] + l[idim];*/
            d2+=x[idim]*x[idim];
        }
        return d2;
    }

    size_t get_size_posizioni() const { return buffer_posizioni_size; }
    size_t get_size_cm() const {return buffer_cm_size;}
    std::vector<ssize_t> get_shape() {
        return {static_cast<ssize_t>(loaded_timesteps),
                static_cast<ssize_t>(natoms),
                3};
    }
    std::vector<ssize_t> get_stride() {
        return { static_cast<ssize_t> (natoms*3*sizeof(double)),
                 static_cast<ssize_t>(3*sizeof (double)),
                 static_cast<ssize_t>(sizeof(double))};
    }

    void toggle_pos_vel() {
        serve_pos=!serve_pos;
    }

    bool serving_pos() {
        return serve_pos;
    }

    size_t get_nloaded_timesteps() const {
        return loaded_timesteps;
    }

protected:

    double * buffer_posizioni; //velocita' e posizioni copiate dal file caricato con mmap, in ordine (nela traiettoria di LAMMPS sono disordinate)
    double * buffer_velocita;
    double * masse;
    double * cariche;
    double * buffer_scatola; //dimensioni della simulazione ad ogni timestep
    double * buffer_posizioni_cm; // posizioni del centro di massa
    double * buffer_velocita_cm; // velocità del centro di massa
    size_t buffer_posizioni_size, buffer_cm_size; //sizes of allocated buffers

    int * buffer_tipi,*buffer_tipi_id;
    size_t natoms,ntypes,min_type,max_type,n_timesteps, loaded_timesteps;
    bool wrap_pbc;


    std::vector<unsigned int> types;
    std::map<int,unsigned int>type_map;

private:
    bool serve_pos;
};

#endif // TRAIETTORIABASE_H
