#ifndef TRAIETTORIA_NUMPY_H
#define TRAIETTORIA_NUMPY_H

#include "traiettoriabase.h"
#include "pybind11/pybind11.h"

class Traiettoria_numpy : public TraiettoriaBase<Traiettoria_numpy>
{
public:
    Traiettoria_numpy(pybind11::buffer &&buffer_pos, pybind11::buffer &&buffer_vel, pybind11::buffer &&buffer_types, pybind11::buffer &&buffer_box,bool lammps_box=true,bool pbc_wrap=false);
    ~Traiettoria_numpy();
    double * posizioni (const int & timestep, const int & atomo) {return buffer_posizioni+natoms*3*timestep+atomo*3;}
    double * velocita (const int & timestep, const int & atomo){return buffer_velocita+natoms*3*timestep+atomo*3;}
    double * scatola (const int & timestep){return buffer_scatola+timestep*6;}
    double * posizioni_cm(const int & timestep, const int & tipo){return buffer_posizioni_cm+timestep*ntypes*3 + tipo*3;}
    double * velocita_cm(const int & timestep, const int & tipo){return buffer_velocita_cm+timestep*ntypes*3 + tipo*3;}
    double *scatola_last(){return buffer_scatola + (n_timesteps-1)*6; }
    void dump_lammps_bin_traj(const std::string &fname, int start_ts, int stop_ts);

private:
    pybind11::buffer buffer_pos,buffer_vel,buffer_types,buffer_box;
    bool lammps_box,posizioni_allocated;
    void calc_cm_pos_vel(double * a, double *&cm);
};

#endif // TRAIETTORIA_NUMPY_H
