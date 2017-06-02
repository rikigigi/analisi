/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
  *
**/


#include "chargefluxts.h"


ChargeFluxTs::ChargeFluxTs(Traiettoria *t) :
    traiettoria(t),J(0),ntimesteps(0)
{

}

void ChargeFluxTs::reset(unsigned int numeroTimestepPerBlocco) {
    if (ntimesteps!=numeroTimestepPerBlocco) {
        delete [] J;
        ntimesteps=numeroTimestepPerBlocco;
        J=new double [ntimesteps*3];
    }
}

void ChargeFluxTs::calcola(unsigned int primo_) {
    primo=primo_;
    for (unsigned int istep=0;istep<ntimesteps;istep++){
        unsigned int timestep=primo+istep;
        for (unsigned int i=0;i<3;i++)
            J[istep*3+i]=0.0;

        for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
            for (unsigned int i=0;i<3;i++){
                J[istep*3+i]+=traiettoria->velocita(timestep,iatom)[i]*traiettoria->get_charge(traiettoria->get_type(iatom));
            }
        }
        for (unsigned int i=0;i<3;i++)
            J[istep*3+i]/=double(traiettoria->get_natoms());
    }
}

double * ChargeFluxTs::J_z(const unsigned int & timestep) {

    if (timestep >= primo && timestep <primo+ntimesteps) {
        return & J[(timestep-primo)*3];
    } else {
        std::cerr << "Errore: richiesta la corrente di carica per un timestep fuori dall'intervallo per cui la corrente è stata calcolata!\n";
        abort();
        return 0;
    }

}
