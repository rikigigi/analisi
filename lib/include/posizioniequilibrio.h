/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef POSIZIONIEQUILIBRIO_H
#define POSIZIONIEQUILIBRIO_H

#include "operazionisulista.h"
#include "config.h"
#include <map>
#include <array>
#include <vector>
#include "eigen_include.h"
#include "traiettoria.h"

class PosizioniEquilibrio : public VectorOp<PosizioniEquilibrio>
{
public:
    PosizioniEquilibrio(Traiettoria *, unsigned int timesteps_sottoblocco=0);
    ~PosizioniEquilibrio();
    unsigned int nExtraTimesteps(unsigned int n_b);
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calculate(unsigned int primo);
    void fit_nacl();
    double d2_reticolo_spostamento_medio(double *min, double *max, double *spostamento);
    bool zona_brillouin(double * k_test);
    bool zona_brillouin_simple_cubic(double *k_test);
    std::map<std::array<int,4>,std::pair<Eigen::Vector3d,unsigned int> > & get_reticolo();
    unsigned int get_atom_index_origin_cell(unsigned int icell=0);
    double * get_atom_position_origin_cell (unsigned int icell=0);
    double * get_atom_position(unsigned int iatom);
    double * get_fitted_pos(unsigned int iatom);
    void get_displacement(unsigned int iatom,unsigned int tstep,double * disp);
    void get_brillouin_limit(double *kmax,double *kmin);
    unsigned int get_atoms_cell() {return atomi_per_cella;}
    unsigned int get_number_cells();
    unsigned int get_atom_base_index(unsigned int iatom);
    unsigned int get_type_base(unsigned int ibase);
    double get_simulation_size(unsigned int icoord);
    void calculate_new_pos();
    std::vector<std::array<double,3> > & get_atom_nearest_image_translation(unsigned int iatom);
private:
    Eigen::Matrix3d reticolo,reticolo_inv,reticolo_reciproco;
    Eigen::Matrix3Xd base_reticolo;
    std::vector< std::array<double,3> > *traslation;
    unsigned int *base_tipi,atomi_per_cella;
    unsigned int *origin_cell_idx;
    double origine_reticolo[3];
    double simulation_size[3];
    double *posizioni_fittate_reticolo;
    std::map<std::array<int,4>,std::pair<Eigen::Vector3d,unsigned int> > reticolo_xyz;
    std::vector<std::array<int,4>> coordinate_reticolo_atomi;
    void reticolo_inizializza(double * base_r, double * base, unsigned int *base_type, const unsigned int & nbase);
    void coord_reticolo(double *xyz,double *uvw_min,double *uvw_max);
    bool calcolato;
    Traiettoria *traiettoria;
    unsigned int lunghezza_media,timestepSottoblocco;
};

#endif // POSIZIONIEQUILIBRIO_H
