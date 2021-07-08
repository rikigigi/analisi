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
#include <vector>

template <class T,bool FPE=false>
class MSD : public OperazioniSuLista<MSD<T,FPE> >
{
public:
    MSD(T *t,
        unsigned int skip=1,
        unsigned int tmax=0,
        unsigned int nthreads=0,
        bool calcola_msd_centro_di_massa=false,
        bool calcola_msd_nel_sistema_del_centro_di_massa=false,
        bool debug=false
            );
    std::vector<ssize_t> get_shape() const { return {static_cast<ssize_t> (leff),static_cast<ssize_t>(f_cm),static_cast<ssize_t>(ntypes)} ; }
    std::vector<ssize_t> get_stride() const { return { static_cast<ssize_t> (ntypes*f_cm*sizeof(double)),static_cast<ssize_t>(ntypes*sizeof (double)),static_cast<ssize_t>(sizeof(double))};}
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int primo);
    MSD<T,FPE> & operator =(const MSD<T,FPE> & destra);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);

private:
    using OperazioniSuLista<MSD<T,FPE> >::lista;
    using OperazioniSuLista<MSD<T,FPE> >::lunghezza_lista;
    T * traiettoria;
    size_t ntimesteps,skip,lmax,leff,nthread,f_cm,ntypes;
    bool cm_msd,cm_self,debug;
};

#endif // MSD_H
