#ifndef MSD_H
#define MSD_H

#include "operazionisulista.h"
#include "traiettoria.h"

class MSD : public OperazioniSuLista<MSD>
{
public:
    MSD(Traiettoria *t,
        unsigned int skip=1,
        unsigned int tmax=0
            );

    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    MSD & operator =(const MSD & destra);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
private:
    Traiettoria * traiettoria;
    unsigned int ntimesteps,skip,lmax,leff;
};

#endif // MSD_H
