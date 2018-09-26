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
