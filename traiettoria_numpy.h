#ifndef TRAIETTORIA_NUMPY_H
#define TRAIETTORIA_NUMPY_H

#include "traiettoriabase.h"
#include "pybind11/pybind11.h"

class Traiettoria_numpy : public TraiettoriaBase<Traiettoria_numpy>
{
public:
    Traiettoria_numpy(pybind11::buffer buffer_pos, pybind11::buffer buffer_vel, pybind11::buffer buffer_types, pybind11::buffer buffer_box,bool lammps_box=true,bool pbc_wrap=false);
    ~Traiettoria_numpy();

private:
    bool lammps_box,posizioni_allocated;
};

#endif // TRAIETTORIA_NUMPY_H
