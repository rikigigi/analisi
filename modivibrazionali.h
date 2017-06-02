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



#ifndef MODIVIBRAZIONALI_H
#define MODIVIBRAZIONALI_H

#include "operazionisulista.h"
#include "posizioniequilibrio.h"
#include <string>
#include <complex>

class ModiVibrazionali : public OperazioniSuLista<ModiVibrazionali>
{
public:
    ModiVibrazionali(Traiettoria * tr, std::string ifcfile, std::string fononefile, unsigned int n_threads, unsigned int timestep_blocco=0);
    void read_force_file(std::string f);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    void azzera();
    void reset(unsigned int s);
    void calcola(unsigned int primo);

    /*
    ModiVibrazionali &operator =(const ModiVibrazionali &);
    ModiVibrazionali & operator+= (const ModiVibrazionali&);
    ModiVibrazionali & operator-= (const ModiVibrazionali&);
    ModiVibrazionali & operator*= (const ModiVibrazionali&);
    ModiVibrazionali & operator/= (const ModiVibrazionali&);
    */
private:
    Eigen::MatrixXcd autovettori;
    Eigen::VectorXd autovalori;
//    Eigen
    Eigen::Matrix3Xd vettori_onda;
    std::string fileFononi;
    unsigned int numero_timesteps,numero_threads,timestepBlocco;
    Traiettoria * traiettoria;
    PosizioniEquilibrio * posizioni_equilibrio;
    Eigen::MatrixXd ifc;
};

#endif // MODIVIBRAZIONALI_H
