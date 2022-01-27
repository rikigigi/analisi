#ifndef TRAIETTORIABASE_H
#define TRAIETTORIABASE_H

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <map>
#include <algorithm>
#include "compiler_types.h"

#include "restrict.h"
#include <cmath>

#define DECL_ARG(a,b) a b
#define NAME_ARG(a,b) b
#define DECL_CALL_BASE_0(ret, fname) ret fname () { return static_cast<T*>(this)-> fname ( ); }
#define DECL_CALL_BASE_1(ret, fname, arg1) ret fname ( DECL_ARG arg1) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1 ); }
#define DECL_CALL_BASE_2(ret, fname, arg1, arg2) ret fname ( DECL_ARG arg1, DECL_ARG arg2) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2 ); }
#define DECL_CALL_BASE_3(ret, fname, arg1, arg2, arg3) ret fname ( DECL_ARG arg1, DECL_ARG arg2, DECL_ARG arg3) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2, NAME_ARG arg3 ); }
#define DECL_CALL_BASE_4(ret, fname, arg1, arg2, arg3, arg4) ret fname ( DECL_ARG arg1, DECL_ARG arg2, DECL_ARG arg3, DECL_ARG arg4) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2, NAME_ARG arg3, NAME_ARG arg4 ); }
#define DECL_CALL_BASE_5(ret, fname, arg1, arg2, arg3, arg4, arg5) ret fname ( DECL_ARG arg1, DECL_ARG arg2, DECL_ARG arg3, DECL_ARG arg4, DECL_ARG arg5) { return static_cast<T*>(this)-> fname ( NAME_ARG arg1, NAME_ARG arg2, NAME_ARG arg3, NAME_ARG arg4, NAME_ARG arg5 ); }


struct Cell_ {
    const double xlo,ylo,zlo,lxhalf,lyhalf,lzhalf, xy,xz,yz;
    void middle(double * c) const {
        c[0]=lxhalf;
        c[1]=lyhalf;
        c[2]=lzhalf;
    }
};
template <bool TRICLINIC>
struct Cell;

template<>
struct Cell<true> : public Cell_ {
    Cell(double *c) :
        Cell_{c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8]}
    {}
};
template<>
struct Cell<false> : public Cell_ {
    Cell(double *c) :
        Cell_{c[0],c[1],c[2],c[3],c[4],c[5],0,0,0}
    {}
};

template <class T>
class TraiettoriaBase {


public:
    enum class BoxFormat{Lammps_ortho,Lammps_triclinic,Cell_vectors, Invalid};

    TraiettoriaBase() : ntypes{0},n_timesteps{0},natoms{0},min_type{0},max_type{0},
    wrap_pbc{true}, buffer_posizioni{nullptr}, buffer_velocita{nullptr},
    buffer_scatola{nullptr}, buffer_posizioni_cm{nullptr},
    buffer_velocita_cm{nullptr}, masse{nullptr}, cariche{nullptr},
    buffer_tipi{nullptr},buffer_tipi_id{nullptr}, serve_pos{true},
    box_format{BoxFormat::Invalid}{}

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
            throw std::runtime_error("Atom index out of range\n");
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

    static void lammps_to_internal(double * c) {
        // xlo,  xhi,   ylo,   yhi,     zlo,     zhi
        // becomes:
        // xlo,  ylo,   zlo,   (xhi-xlo) / 2, (yhi-ylo) / 2, (zhi-zlo) / 2
        //first low coords and then high-low coordinates
        double t=c[4];
        c[5]=(c[5]-c[4])/2;
        c[4]=(c[3]-c[2])/2;
        c[3]=(c[1]-c[0])/2;
        c[1]=c[2];
        c[2]=t; //c[4]
    }
    static void internal_to_lammps(double * c) {
        //first low coords and then high-low coordinates
        // xlo,  ylo,   zlo,   (xhi-xlo) / 2, (yhi-ylo) / 2, (zhi-zlo) / 2
        // becomes:
        // xlo,  xhi,   ylo,   yhi,     zlo,     zhi
        double t=c[3];
        c[3]=c[4]*2+c[1];
        c[5]=c[5]*2+c[2];
        c[4]=c[2];
        c[2]=c[1];
        c[1]=c[0]+t*2;
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
    double get_charge(unsigned int  i){if (i<get_ntypes()) return cariche[i]; throw std::runtime_error("Cannot get charge for a type that does not exists!\n");}


    /**
      * this functions wraps the positions around the center of the orthorombic cell. Positions if the cell is not orthorombic will be outside the cell,
      * but they are more compact making the min image algorithm more efficient
    **/
    template<bool TRICLINIC>
    void pbc_wrap(ssize_t idx) {
        double * c = buffer_scatola+buffer_scatola_stride*idx;
        Cell<TRICLINIC> s(c);
        double x0[3];
        s.middle(x0);
        for (size_t iatom=0;iatom<natoms;++iatom) {
            double *xa=buffer_posizioni+idx*natoms*3+iatom*3;
            for (int icoord=0;icoord<3;++icoord){
                xa[icoord]=xa[icoord]-x0[icoord];
            }
            minImage_triclinic<TRICLINIC>(xa,c+3,c+6);
            for (int icoord=0;icoord<3;++icoord){
                xa[icoord]=xa[icoord]+x0[icoord];
            }
        }
    }


    double d2_minImage(size_t i,
                       size_t j,
                       size_t itimestep,
                       size_t jtimestep
                       ) {
        double x[3];
        return d2_minImage(i,j,itimestep,jtimestep,x);
    }

    double d2_minImage(size_t i,
                       size_t j,
                       size_t itimestep,
                       size_t jtimestep,
                       double *x
                       ) {
        double d2=0.0;
        double *xi=posizioni(itimestep,i);
        double *xj=posizioni(jtimestep,j);
        const double *l=scatola(itimestep)+3;
        const double *xy_xz_yz = scatola(itimestep)+6;
        d2=d2_minImage_triclinic(xi,xj,l,x,xy_xz_yz);
        return d2;
    }
    double d2_minImage_triclinic(double * xi, ///position of first atom
                                 double *xj, ///position of second atom
                                 const double *l, /// half orthorombic edges
                                 double *x, /// min image distance
                                 const double * xy_xz_yz ///offsets of triclinic wrt orthorombic
                                 ) const{
        double d2=0.0;
        for (unsigned int idim=0;idim<3;idim++) {
            x[idim]=xi[idim]-xj[idim];
        }
        if (triclinic){
            minImage_triclinic<true>(x,l,xy_xz_yz);
        } else {
            minImage_triclinic<false>(x,l,xy_xz_yz);
        }
        for (unsigned int idim=0;idim<3;idim++) {
            d2+=x[idim]*x[idim];
        }
        return d2;
    }
    template <bool TRICLINIC>
    static void minImage_triclinic(double * __restrict delta,
                                   const double * __restrict l_half,
                                   const double * __restrict xy_xz_yz) {
        double xy;
        double xz;
        double yz;
        if constexpr (TRICLINIC) {
            xy=xy_xz_yz[0];
            xz=xy_xz_yz[1];
            yz=xy_xz_yz[2];
        }
        while (fabs(delta[2]) > l_half[2]) {
            if (delta[2] < 0.0) {
                delta[2] += l_half[2]*2;
                if constexpr(TRICLINIC){
                    delta[1] += yz;
                    delta[0] += xz;
                }
            } else {
                delta[2] -= l_half[2]*2;
                if constexpr(TRICLINIC){
                    delta[1] -= yz;
                    delta[0] -= xz;
                }
            }
        }
        while (fabs(delta[1]) > l_half[1]) {
            if (delta[1] < 0.0) {
                delta[1] += l_half[1]*2;
                if constexpr(TRICLINIC){
                    delta[0] += xy;
                }
            } else {
                delta[1] -= l_half[1]*2;
                if constexpr(TRICLINIC){
                    delta[0] -= xy;
                }
            }
        }
        while (fabs(delta[0]) > l_half[0]) {
            if (delta[0] < 0.0) delta[0] += l_half[0]*2;
            else delta[0] -= l_half[0]*2;
        }
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

    bool is_triclinic() const {
        return triclinic;
    }

    size_t get_box_stride() const {return buffer_scatola_stride;}

    //functions that are compiled in cpp file
    void dump_lammps_bin_traj(const std::string &fname, int start_ts, int stop_ts);
    int get_ntypes ();



protected:

    double * buffer_posizioni; //velocita' e posizioni copiate dal file caricato con mmap, in ordine (nela traiettoria di LAMMPS sono disordinate)
    double * buffer_velocita;
    double * masse;
    double * cariche;
    double * buffer_scatola; //dimensioni della simulazione ad ogni timestep
    size_t buffer_scatola_stride; // 6 for orthorombic, 9 for triclinic
    BoxFormat box_format; //format used to store box information
    double * buffer_posizioni_cm; // posizioni del centro di massa
    double * buffer_velocita_cm; // velocità del centro di massa
    size_t buffer_posizioni_size, buffer_cm_size; //sizes of allocated buffers

    int * buffer_tipi,*buffer_tipi_id;
    ssize_t natoms,ntypes,min_type,max_type,n_timesteps, loaded_timesteps;
    bool wrap_pbc, triclinic;


    std::vector<unsigned int> types;
    std::map<int,unsigned int>type_map;

private:
    bool serve_pos;
};

#endif // TRAIETTORIABASE_H
