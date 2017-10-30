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



#ifndef SPETTROVIBRAZIONALE_H
#define SPETTROVIBRAZIONALE_H

#include "mediablocchi.h"
#include "operazionisulista.h"
#include <fftw3.h>
#include "mediablocchi.h"
//#include <memory>

class Traiettoria;


class SpettroVibrazionale : public OperazioniSuLista<SpettroVibrazionale>
{
public:
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    SpettroVibrazionale(Traiettoria *);
    ~SpettroVibrazionale();
    static void deallocate_plan(); // da chiamare alla fine del programma!
    double spettro(unsigned int frequenza, unsigned int dim, unsigned int tipo_atomo);
    SpettroVibrazionale & operator = (const SpettroVibrazionale &);


private:
    Traiettoria * traiettoria;
    unsigned int size;
    int tipi_atomi;
    unsigned int trasformata_size;
    fftw_complex * trasformata;
    static fftw_plan fftw3;
    static unsigned int fplan_natoms,fplan_size;



};

//per fare anche le varie medie a blocchi
//template class MediaBlocchi<SpettroVibrazionale>;

#endif // SPETTROVIBRAZIONALE_H
