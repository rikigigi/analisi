#ifndef CONVERTIBINARIO_H
#define CONVERTIBINARIO_H

#include <string>

class ConvertiBinario
{
public:
    enum Type{lammps,natoms_box_xyz_vxvyvz};
    ConvertiBinario(std::string filein,std::string fileout,Type tipo=natoms_box_xyz_vxvyvz);
};

#endif // CONVERTIBINARIO_H
