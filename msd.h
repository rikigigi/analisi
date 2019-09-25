/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef MSD_H
#define MSD_H

#include "operazionisulista.h"
#include "traiettoria.h"

class MSD : public OperazioniSuLista<MSD>
{
public:
    MSD(Traiettoria *t,
        unsigned int skip=1,
        unsigned int tmax=0,
        unsigned int nthreads=0,
        bool calcola_msd_centro_di_massa=false,
        bool calcola_msd_nel_sistema_del_centro_di_massa=false,
        bool debug=false
            );

    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    MSD & operator =(const MSD & destra);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
private:
    Traiettoria * traiettoria;
    unsigned int ntimesteps,skip,lmax,leff,nthread,f_cm;
    bool cm_msd,cm_self,debug;
};

#endif // MSD_H
