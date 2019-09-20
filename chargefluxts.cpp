/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/


#include "chargefluxts.h"


ChargeFluxTs::ChargeFluxTs(Traiettoria *t) :
    traiettoria(t),J(0),timestep_finestra(0),calcolato(false)
{

}

ChargeFluxTs::~ChargeFluxTs() {
    delete [] J;
}

void ChargeFluxTs::reset(unsigned int numeroTimestepPerBlocco) {
    if (timestep_finestra!=numeroTimestepPerBlocco) {
        delete [] J;
        timestep_finestra=numeroTimestepPerBlocco;
        J=new double [timestep_finestra*3];
        calcolato=false;
    }
}

void ChargeFluxTs::calcola(unsigned int timestep) {



    //se ho già calcolato i dati per certi timesteps, li ricopia semplicemente senza calcolare tutto di nuovo.

    int     timestep_copy_tstart=0,timestep_copy_tend=0,
            timestep_read_start=0, timestep_read_end=timestep_finestra,
            finestra_differenza=timestep_corrente-timestep;
    if (calcolato) {
        if (abs(finestra_differenza)<timestep_finestra) { // si sovrappongono
            if (finestra_differenza>0){
                timestep_copy_tstart=0;
                timestep_copy_tend=timestep_finestra-finestra_differenza;

                timestep_read_end=finestra_differenza;
                timestep_read_start=0;
            } else {
                timestep_copy_tstart=-finestra_differenza;
                timestep_copy_tend=timestep_finestra;

                timestep_read_start=finestra_differenza+timestep_finestra;
                timestep_read_end=timestep_finestra;
            }
        }
    }

    //copia i dati già calcolati

    for (int i=timestep_copy_tstart;i<timestep_copy_tend;i++) {
        for (int idata=0;idata<3;idata++) {
            J[(finestra_differenza+i)*3+idata]=J[i*3+idata];
        }
    }




    for (unsigned int istep=timestep_read_start;istep<timestep_read_end;istep++){
        unsigned int t=timestep+istep;
        for (unsigned int i=0;i<3;i++)
            J[istep*3+i]=0.0;

        for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
            for (unsigned int i=0;i<3;i++){
                J[istep*3+i]+=traiettoria->velocita(t,iatom)[i]
                        *traiettoria->get_charge(traiettoria->get_type(iatom));
            }
        }
        for (unsigned int i=0;i<3;i++)
            J[istep*3+i]/=double(traiettoria->get_natoms());
    }
    timestep_corrente=timestep;
    calcolato=true;
}

double * ChargeFluxTs::J_z(const unsigned int & timestep) {

    if (timestep >= timestep_corrente && timestep <timestep_corrente+timestep_finestra) {
        return & J[(timestep-timestep_corrente)*3];
    } else {
        std::cerr << "Errore: richiesta la corrente di carica per un timestep fuori dall'intervallo per cui la corrente è stata calcolata!\n";
        abort();
        return 0;
    }

}
