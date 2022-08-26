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
#define DECL_CALL_BASE_1(ret, fname, arg1) template<bool SAFE=true> ret fname ( DECL_ARG arg1) { return static_cast<T*>(this)->template fname<SAFE> ( NAME_ARG arg1 ); }
#define DECL_CALL_BASE_2(ret, fname, arg1, arg2) template<bool SAFE=true> ret fname ( DECL_ARG arg1, DECL_ARG arg2) { return static_cast<T*>(this)->template fname<SAFE> ( NAME_ARG arg1, NAME_ARG arg2 ); }


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
class BaseTrajectory {


public:
    enum class BoxFormat{Lammps_ortho,Lammps_triclinic,Cell_vectors, Invalid};

    BaseTrajectory() : ntypes{0},n_timesteps{0},natoms{0},min_type{0},max_type{0},current_timestep{0},
    wrap_pbc{true}, buffer_positions{nullptr}, buffer_velocity{nullptr},
    buffer_boxes{nullptr}, buffer_positions_cm{nullptr},
    buffer_velocity_cm{nullptr}, mass{nullptr}, charge{nullptr},
    buffer_type{nullptr},buffer_type_id{nullptr}, serve_pos{true},
    box_format{BoxFormat::Invalid}{}

    //double * positions (const int & timestep, const int & atomo) { return static_cast<T*>(this)->positions(timestep,atomo);}
    DECL_CALL_BASE_2(double *, positions, (const int &, timestep), (const int &, atomo))
    DECL_CALL_BASE_2(double *, velocity,(const int &, timestep), (const int &, atomo) )
    DECL_CALL_BASE_1(double *, box, (const int &, timestep))
    DECL_CALL_BASE_2(double *, positions_cm,(const int &, timestep), (const int &, tipo))
    DECL_CALL_BASE_2(double *, velocity_cm,(const int &, timestep), (const int &, tipo))
    DECL_CALL_BASE_0(double *,box_last)
    std::vector<unsigned int> get_types(){
        get_ntypes();
        return types;

    }
    unsigned int get_type(const unsigned int &atomo){
        if (atomo < natoms) {
            return buffer_type_id[atomo];
        } else {
            throw std::runtime_error("Atom index out of range\n");
        }
    }
    enum Errori {non_inizializzato=0,oltre_fine_file=2,Ok=1};
    Errori set_data_access_block_size(const size_t & timesteps){std::cerr << "Warning: doing nothing (not reading in blocks)"<<std::endl; return Errori::Ok;}
    Errori set_access_at(const size_t & timesteps){std::cerr << "Warning: doing nothing (not reading in blocks)"<<std::endl;return Errori::Ok;}

    //I have to set this before loading the trajectory
    void set_pbc_wrap(bool p){
        wrap_pbc=p;
    }

    bool get_pbc_wrap() const {
        return wrap_pbc;
    }

    /**
     * Convert lammps format to internal format
     * **/
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
    /**
     * Convert internal format to lammps format
     * **/
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

    /**
     * return a pointer to the beginning of the stored positions
     * **/
    double * positions_data(){return buffer_positions;}
    /**
     * return a pointer to the beginning of the stored velocities
     * **/
    double * velocity_data(){return buffer_velocity;}
    int get_type_min() {return min_type;}
    int get_type_max() {return max_type;}
    size_t get_natoms()const {return natoms;}
    size_t get_ntimesteps() const{return n_timesteps;}
    ssize_t get_current_timestep() const {return current_timestep;}
    double get_mass(unsigned int i) {if (i<get_ntypes()) return mass[i]; throw std::runtime_error("Cannot get mass for a type that does not exists!\n");}
    void set_mass(unsigned int i,double m) {if (i<get_ntypes()) mass[i]=m;}
    void set_charge(unsigned int i, double c){if (i<get_ntypes()) charge[i]=c;}
    double get_charge(unsigned int  i){if (i<get_ntypes()) return charge[i]; throw std::runtime_error("Cannot get charge for a type that does not exists!\n");}


    /**
      * this functions wraps the positions around the center of the orthorombic cell. Positions if the cell is not orthorombic will be outside the cell,
      * but they are more compact making the min image algorithm more efficient
    **/
    template<bool TRICLINIC>
    void pbc_wrap(ssize_t idx) {
        double * c = buffer_boxes+buffer_boxes_stride*idx;
        Cell<TRICLINIC> s(c);
        double x0[3];
        s.middle(x0);
        for (size_t iatom=0;iatom<natoms;++iatom) {
            double *xa=buffer_positions+idx*natoms*3+iatom*3;
            for (int icoord=0;icoord<3;++icoord){
                xa[icoord]=xa[icoord]-x0[icoord];
            }
            minImage_triclinic<TRICLINIC>(xa,c+3,c+6);
            for (int icoord=0;icoord<3;++icoord){
                xa[icoord]=xa[icoord]+x0[icoord];
            }
        }
    }


    /**
     * This returns the square of the minimum image distance between the specified timesteps and atom
     * using the box data of timestep itimestep
     * **/
    double d2_minImage(size_t i,
                       size_t j,
                       size_t itimestep,
                       size_t jtimestep
                       ) {
        double x[3];
        return d2_minImage(i,j,itimestep,jtimestep,x);
    }

    /**
     * This returns the square of the minimum image distance between the specified timesteps and atom
     * using the box data of timestep itimestep. Sets x(i) - x(j) in memory located at x
     * **/
    double d2_minImage(size_t i,
                       size_t j,
                       size_t itimestep,
                       size_t jtimestep,
                       double *x
                       ) {
        double d2=0.0;
        double *xi=positions<false>(itimestep,i);
        double *xj=positions<false>(jtimestep,j);
        const double *l=box<false>(itimestep)+3;
        const double *xy_xz_yz = box<false>(itimestep)+6;
        d2=d2_minImage_triclinic(xi,xj,l,x,xy_xz_yz);
        return d2;
    }
    /**
     * This returns the square of the minimum image distance between the specified timesteps and atom
     * using the box data specified in the arguments (if not triclinc the last argument is not used)
     * Sets x(i) - x(j) in memory located at x
     * **/
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
    /**
     * apply the minimum image convention in the triclinic and not triclinc case.
     * Must provide box data as in the internal format used by the code.
     * **/
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


    size_t get_positions_size() const { return buffer_positions_size; }
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

    size_t get_box_stride() const {return buffer_boxes_stride;}

    //functions that are compiled in cpp file
    void dump_lammps_bin_traj(const std::string &fname, int start_ts, int stop_ts);
    size_t get_ntypes();



protected:

    double * buffer_positions; //positions and velocities may need a copy because in the lammps format they are not ordered
    double * buffer_velocity;
    double * mass; //used almost nowhere
    double * charge;
    double * buffer_boxes; //box data in specified format
    size_t buffer_boxes_stride; // 6 for orthorombic, 9 for triclinic
    BoxFormat box_format; //format used to store box information
    double * buffer_positions_cm; // center of mass position for each atomic type
    double * buffer_velocity_cm; // center of mass velocity for each atomic type
    size_t buffer_positions_size, buffer_cm_size; //sizes of allocated buffers

    int * buffer_type,*buffer_type_id;
    ssize_t natoms,ntypes,min_type,max_type,n_timesteps, loaded_timesteps,current_timestep;
    bool wrap_pbc, triclinic;


    std::vector<unsigned int> types;
    std::map<int,unsigned int>type_map;

private:
    bool serve_pos;
};

#endif // TRAIETTORIABASE_H
