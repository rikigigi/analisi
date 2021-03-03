/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef CONVERTIBINARIO_H
#define CONVERTIBINARIO_H

#include <string>

class ConvertiBinario
{
public:
    enum Type{lammps,natoms_box_xyz_vxvyvz,gromax_trr};
    ConvertiBinario(const std::string filein, const std::string fileout, Type tipo=natoms_box_xyz_vxvyvz, const std::string typefile="");
};

#endif // CONVERTIBINARIO_H
