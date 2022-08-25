/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef MODIVIBRAZIONALI_H
#define MODIVIBRAZIONALI_H

#include "operazionisulista.h"
#include <string>
#include <complex>
#include "traiettoria.h"
#include "posizioniequilibrio.h"

class ModiVibrazionali : public VectorOp<ModiVibrazionali>
{
public:
    ModiVibrazionali(Traiettoria * tr, std::string ifcfile, std::string fononefile, unsigned int n_threads, unsigned int timestep_blocco=0);
    void read_force_file(std::string f);
    unsigned int nExtraTimesteps(unsigned int n_b);
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
