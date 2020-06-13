#ifndef TRAIETTORIA_NUMPY_H
#define TRAIETTORIA_NUMPY_H

#include "traiettoriabase.h"
#include "pybind11/pybind11.h"

class Traiettoria_numpy : public TraiettoriaBase<Traiettoria_numpy>
{
public:
    Traiettoria_numpy(pybind11::buffer buffer_pos, pybind11::buffer buffer_vel, pybind11::buffer buffer_types, pybind11::buffer buffer_box,bool lammps_box=true,bool pbc_wrap=false);
    ~Traiettoria_numpy();
    double * posizioni (const int & timestep, const int & atomo) {return buffer_posizioni+natoms*3*timestep+atomo*3;}
    double * velocita (const int & timestep, const int & atomo){return buffer_velocita+natoms*3*timestep+atomo*3;}
    double * scatola (const int & timestep){return buffer_scatola+timestep*6;}
    double * posizioni_cm(const int & timestep, const int & tipo){throw std::runtime_error("Center of mass positions not implemented");}
    double * velocita_cm(const int & timestep, const int & tipo){throw std::runtime_error("Center of mass velocities not implemented");}
    double *scatola_last(){return buffer_scatola + (n_timesteps-1)*6; }
    using TraiettoriaBase<Traiettoria_numpy>::get_ntimesteps;

private:
    bool lammps_box,posizioni_allocated;
};

#endif // TRAIETTORIA_NUMPY_H
