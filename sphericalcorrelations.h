#ifndef SPHERICALCORRELATIONS_H
#define SPHERICALCORRELATIONS_H

#include "operazionisulista.h"



template <int l,class TFLOAT, class T>
class SphericalCorrelations : public OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>
{
public:
    SphericalCorrelations(T *t,
                          TFLOAT rmin,
                          TFLOAT rmax,
                          unsigned int nbin,
                          unsigned int tmax=0,
                          unsigned int nthreads=0,
                          unsigned int skip=1,
                          bool debug=false);
    ~SphericalCorrelations();
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    SphericalCorrelations<l,TFLOAT,T> & operator =(const SphericalCorrelations<l,TFLOAT,T> & destra){
        OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::operator = (destra);
        return *this;
    }
    std::vector<ssize_t> get_shape();
    std::vector<ssize_t> get_stride();
    std::string get_columns_description() {return c_descr;}
    inline int index(const int t, const int type1, const int type2,const int ibin=0) const noexcept {
        return (l+1)*(l+1)*(nbin*(ntypes*(ntypes*t + type1) + type2)+ibin);
    }
    using OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::azzera;
private:
    using OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::lista;
    using OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::lunghezza_lista;
    T & t;

    inline int index_wrk(const int iatom,const int jtype,const int ibin=0) const noexcept {
        return (l+1)*(l+1)*(nbin*(ntypes*iatom+jtype)+ibin);
    }

    inline void sh_snapshot(int timestep, TFLOAT * workspace, TFLOAT * cheby, TFLOAT * result, double *) ;

    TFLOAT rmin, rmax,dr;
    unsigned int nbin, tmax,nthreads,skip,leff,ntimesteps,ntypes,natoms;
    bool debug;
    std::string c_descr;

};

#endif // SPHERICALCORRELATIONS_H
