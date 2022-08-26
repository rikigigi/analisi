#ifndef TRAIETTORIA_NUMPY_H
#define TRAIETTORIA_NUMPY_H

#include "basetrajectory.h"
#include "pybind11/pybind11.h"

class Trajectory_numpy : public BaseTrajectory<Trajectory_numpy>
{
public:
    Trajectory_numpy(pybind11::buffer buffer_pos,
                      pybind11::buffer buffer_vel,
                      pybind11::buffer buffer_types,
                      pybind11::buffer buffer_box, ///data of the boxes, in the format specified by matrix_box
                      BaseTrajectory::BoxFormat matrix_box=BaseTrajectory::BoxFormat::Cell_vectors,
                      bool pbc_wrap=false,
                      bool save_rotation_matrix=false);
    ~Trajectory_numpy();
    template <bool SAFE=true>
    double * positions (const int & timestep, const int & atomo) {return buffer_positions+natoms*3*timestep+atomo*3;}
    template <bool SAFE=true>
    double * velocity (const int & timestep, const int & atomo){return buffer_velocity+natoms*3*timestep+atomo*3;}
    template <bool SAFE=true>
    double * box (const int & timestep){return buffer_boxes+timestep*buffer_boxes_stride;}
    template <bool SAFE=true>
    double * positions_cm(const int & timestep, const int & tipo){return buffer_positions_cm+timestep*ntypes*3 + tipo*3;}
    template <bool SAFE=true>
    double * velocity_cm(const int & timestep, const int & tipo){return buffer_velocity_cm+timestep*ntypes*3 + tipo*3;}
    double *box_last(){return buffer_boxes + (n_timesteps-1)*buffer_boxes_stride; }
    double * get_rotation_matrix(size_t t){return rotation_matrix+9*t;}
    using BaseTrajectory::BoxFormat;
private:
    pybind11::buffer buffer_pos,buffer_vel,buffer_types,buffer_box;
    double * rotation_matrix;
    bool posizioni_allocated, velocities_allocated, box_allocated;
    void calc_cm_pos_vel(double * a, double *&cm);
};

#endif // TRAIETTORIA_NUMPY_H
