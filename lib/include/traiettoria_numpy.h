#ifndef TRAIETTORIA_NUMPY_H
#define TRAIETTORIA_NUMPY_H

#include "traiettoriabase.h"
#include "pybind11/pybind11.h"

class Traiettoria_numpy : public TraiettoriaBase<Traiettoria_numpy>
{
public:
    Traiettoria_numpy(pybind11::buffer buffer_pos,
                      pybind11::buffer buffer_vel,
                      pybind11::buffer buffer_types,
                      pybind11::buffer buffer_box, ///data of the boxes, in the format specified by matrix_box
                      TraiettoriaBase::BoxFormat matrix_box=TraiettoriaBase::BoxFormat::Cell_vectors,
                      bool pbc_wrap=false,
                      bool save_rotation_matrix=false);
    ~Traiettoria_numpy();
    template <bool SAFE=true>
    double * posizioni (const int & timestep, const int & atomo) {return buffer_posizioni+natoms*3*timestep+atomo*3;}
    template <bool SAFE=true>
    double * velocita (const int & timestep, const int & atomo){return buffer_velocita+natoms*3*timestep+atomo*3;}
    template <bool SAFE=true>
    double * scatola (const int & timestep){return buffer_scatola+timestep*buffer_scatola_stride;}
    template <bool SAFE=true>
    double * posizioni_cm(const int & timestep, const int & tipo){return buffer_posizioni_cm+timestep*ntypes*3 + tipo*3;}
    template <bool SAFE=true>
    double * velocita_cm(const int & timestep, const int & tipo){return buffer_velocita_cm+timestep*ntypes*3 + tipo*3;}
    double *scatola_last(){return buffer_scatola + (n_timesteps-1)*buffer_scatola_stride; }
    double * get_rotation_matrix(size_t t){return rotation_matrix+9*t;}
    using TraiettoriaBase::BoxFormat;
private:
    pybind11::buffer buffer_pos,buffer_vel,buffer_types,buffer_box;
    double * rotation_matrix;
    bool posizioni_allocated, velocities_allocated, box_allocated;
    void calc_cm_pos_vel(double * a, double *&cm);
};

#endif // TRAIETTORIA_NUMPY_H
