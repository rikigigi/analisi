#include "convertibinario.h"
#include <fstream>
#include <iostream>
#include "lammps_struct.h"

ConvertiBinario::ConvertiBinario(std::string filein, std::string fileout, Type tipo)
{

    std::ofstream out(fileout,std::ofstream::binary);
    std::ifstream in(filein);

    if (tipo!=natoms_box_xyz_vxvyvz) {
        std::cerr << "Non implementato!\n";
        abort();
    }

    bigint itim=0;
    while (in.good()) {
        //per ogni timestep, legge tutti i dati necessari nella struttura dati
        Intestazione_timestep head;
        in >> head.natoms >> head.scatola[0] >> head.scatola[1] >> head.scatola[2] >> head.scatola[3] >> head.scatola[4] >> head.scatola[5];
        head.timestep=++itim;
        if (!in.good()){
            std::cerr << "Errore nella lettura del file (timestep "<<itim<<")!\n";
            break;
        }

        head.triclinic=false;
        head.condizioni_al_contorno[0][0]=0;
        head.condizioni_al_contorno[1][0]=0;
        head.condizioni_al_contorno[2][0]=0;
        head.condizioni_al_contorno[0][1]=0;
        head.condizioni_al_contorno[1][1]=0;
        head.condizioni_al_contorno[2][1]=0;
        head.dimensioni_riga_output=8;
        head.nchunk=1;
        int n_data=head.natoms*8;
        //qui ho impostato l'intestazione, la scrivo
        out.write((char*) &head,sizeof(Intestazione_timestep));
        out.write((char*) &n_data,sizeof(int));
        //adesso devo scrivere tutte le righe, una dopo l'altra
        for (unsigned int i=0;i<head.natoms;++i) {
            double data[8];
            in >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5] >> data[6] >> data[7];
            if (!in.good()){
                std::cerr << "Errore nella lettura del file (timestep "<<itim<<")!\n";
                break;
            }

            out.write((char*) data,sizeof(double)*8);
        }



    }



}
