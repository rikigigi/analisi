/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "testtraiettoria.h"
#include <iostream>
#include <fstream>
TestTraiettoria::TestTraiettoria(std::string filename) : Trajectory(filename)
{

    std::ofstream out("analisi_test.debug");

    int n_timesteps=get_ntimesteps();
    double *coord=new double [n_timesteps];
    double *coord_cm=new double [n_timesteps];
    for (unsigned int i=0;i<n_timesteps;i++){
        coord[i]=0.0;
        coord_cm[i]=0.0;
    }

    unsigned int n_b=31;
    const unsigned int iatom=47,icoord=3;

    unsigned int ts_b=n_timesteps/(n_b+1);

    set_data_access_block_size(ts_b*2);
    for (unsigned int ib=0;ib<n_b;ib++) {
        std::cout << "Blocco "<<ib<<"\n";
        set_access_at(ts_b*ib);
        for (unsigned int ts=ts_b*ib;ts<ts_b*(ib+2);ts++) {
            out << ts <<" " <<positions(ts,iatom)[icoord] << " " << coord[ts] << " "<<positions_cm(ts,0)[icoord]<< " " << coord_cm[ts] <<"\n";
            coord[ts]=positions(ts,iatom)[icoord];
            coord_cm[ts]=positions_cm(ts,0)[icoord];
        }


        std::cout << "\n";

    }

    std::cout << "#--------\n";

    delete [] coord;
    delete [] coord_cm;
}
