/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
  *
**/



#include <stdint.h>
#include <stddef.h>
#include "config.h"
#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif
#include <string>
#include "stdio.h"
#include "string.h"
#include "stdint.h"
#ifdef DEBUG
#include <mcheck.h>
#endif

#define FIT_LIN_NPAR 160
/**Pigreco*/
#ifndef PI
#define PI 3.141592653589793238462643383279502884197169399
#endif
using namespace Eigen;


template <typename T,typename T2> class History{
public:
        History (unsigned int, T * ini=NULL);
        History ();
        ~History();
        bool is_full();
        void add(T);
        void reset_size(unsigned int, T * ini=NULL);
        void reset();
        unsigned int get_n();
        T & operator [] (unsigned int);
        void media_var(T & media,T2 & varianza,unsigned int ini=0,unsigned int fine=0); // media e var devono essere già zero
    void media_var_calc_block(T & media,T2 & varianza,const unsigned int & block_size); // media e var vengono azzerati
        void stat_ineff_calc(unsigned int N_punti=100);
    T2 stat_ineff_media_ultimi(unsigned int N_punti);
        unsigned int get_s_size(){return s_xy_size;}
        T2 get_s_y(unsigned int i) {return s_y[i];}
        unsigned int get_s_x(unsigned int i) {return s_x[i];}
private:
        T *history,media;
        void media_var_calc(T & media,T2 & varianza,unsigned int ini=0,unsigned int fine=0,bool calc_var=true);
        T2 s,* s_y,var; // statistical inefficency
        unsigned int * s_x, size, cur, media_ini,media_fine,s_xy_size;
        bool full,media_var_fatto;
};


class Position_history {
protected:
    double ** pos, **vel, cell[3],*J;
        double * t;
    double *d,*d_var,*n_at,*n_at_var,*corr_v,*corr_v_int,*corr_v_int_dif,*corr_v_int_ecc,*corr_v_var,*corr_v_int_var,m,q,em,eq,chi2,D,D_err,*A,*A_err,tau,tau_err;
    double *corr_J, *corr_J_var,*corr_J_int,*corr_J_int_var,DJ,DJ_err,*fitted_f;
    unsigned int cur_n,max_d,max_d_sc,n_autoc_v,n_autoc_J;
        unsigned int size;
        unsigned int natoms;
    unsigned int medie_cont,media_traiettoria_skip,fit1,fit2;
    unsigned int **idx_A;
    int *types;
    bool full,media_sezioni_scatole;
        unsigned int find_atom(unsigned int,double *);
        double rimasti(unsigned int * , unsigned int * );
    unsigned int destra(unsigned int *);
public:
    unsigned int max_(unsigned int);
        Position_history();
    ~Position_history();
        void allocate_J();
        void allocate_J(unsigned int);
    void hist_natoms_dx(unsigned int *&ist, double &media, double &var, unsigned int &min, unsigned int &max, unsigned int skip_tstep=100);
    void set_media_traiettoria_skip(unsigned int n_);
    void set_scatola_media_sezioni(bool n_);
    bool is_full() {return full;}
        unsigned int get_n(){return max_d;}
    unsigned int get_n_velocity_autocorr(){return n_autoc_v;}
    unsigned int get_n_J_autocorr(){return n_autoc_J;}
        unsigned int get_n_scatola(){return max_d_sc;}
        double operator [] (unsigned int);
    double get_d_var(unsigned int i) {return d_var[i];}
    double scatola(unsigned int i) {return n_at[i];}
    double scatola_fit(unsigned int i) {return fitted_f[i];}
    double scatola_var(unsigned int i) {return n_at_var[i];}
        double * element (unsigned int);
        void add(double, double *,double *,double * J_=NULL);
        void add(double,double * J_);
        void reset_types(int * types_);
        void reset_natoms_size_cell_types(unsigned int n,unsigned int s,double cell_[3],int * types_);
        void reset_natoms_size_cell(unsigned int,unsigned int,double cell_[3]);
        void get_m_q_em_eq_chi2(double &,double &,double &,double &,double &);
        void reset();
        void media_temporale_pos2(unsigned int d_max=0);
        void media_temporale_scatole(unsigned int d_max=0);
    void media_temporale_scatole_set_fit_range(unsigned int f1,unsigned int f2);
    void calc_v_autocorrelation_old(unsigned int d_max=500);
    void calc_v_autocorrelation(unsigned int d_max=500,unsigned int err_last=10);
    void calc_J_autocorrelation(unsigned int d_max=500,unsigned int err_last=10);
        void scatola_exp_fit();
    void get_param_scatola_exp_fit(double *, double &, double *A__, double &tau__);
    double get_D_v_autocorrelation(){return D;}
    double get_D_err_v_autocorrelation(){return D_err;}
    double get_D_J_autocorrelation(){return DJ;}
    double get_D_err_J_autocorrelation(){return DJ_err;}
    double velocity_correlation(unsigned int i) {if(i<n_autoc_v) return corr_v[i]; else return 0.0;}
    double velocity_correlation_var(unsigned int i) {if(i<n_autoc_v) return corr_v_var[i]; else return 0.0;}
    double velocity_correlation_int(unsigned int i) {if(i<n_autoc_v) return corr_v_int[i]; else return 0.0;}
    double velocity_correlation_int_var(unsigned int i) {if(i<n_autoc_v) return corr_v_int_var[i]; else return 0.0;}
    double velocity_correlation_int_dif(unsigned int i) {if(i<n_autoc_v) return corr_v_int_dif[i]; else return 0.0;}
    double velocity_correlation_int_ecc(unsigned int i) {if(i<n_autoc_v) return corr_v_int_ecc[i]; else return 0.0;}
    double J_correlation(unsigned int i) {if(i<n_autoc_J) return corr_J[i]; else return 0.0;}
    double J_correlation_var(unsigned int i) {if(i<n_autoc_J) return corr_J_var[i]; else return 0.0;}
    double J_correlation_int(unsigned int i) {if(i<n_autoc_J) return corr_J_int[i]; else return 0.0;}
    double J_correlation_int_var(unsigned int i) {if(i<n_autoc_J) return corr_J_int_var[i]; else return 0.0;}
    void media_temporale_pos2_(unsigned int d_max);
    void media_temporale_scatole2(unsigned int d_max);
};


//usati da lammps (in binary2txt.cpp)
typedef int tagint;
typedef int64_t bigint;

#ifndef PRId64
#define PRId64 "ld"
#endif
#define BIGINT_FORMAT "%" PRId64
#define INT_DUMP double
#define POS_DUMP double

class Analisi{
public:
    Analisi(std::string traiettoria_fname, std::string energia_fname, unsigned int size);
    ~Analisi();
    Position_history traiettoria;
    History<double,double> T_history, time_history,pressure_history,kin_energy_history,pot_energy_history,T_t,t_t;
    bigint natoms;
    double cell[3];
private:
    bool heat_autocorr;
    double cell_[3];
    bigint ntimestep,natoms_;
    int size_one,nchunk,triclinic_;
    double xlo_,xhi_,ylo_,yhi_,zlo_,zhi_,xy,xz,yz;
    int boundary_[3][2];
    char boundstr[9];
};
