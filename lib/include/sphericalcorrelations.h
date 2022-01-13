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
                          unsigned int buffer_size=10,
                          bool debug=false);
    ~SphericalCorrelations();
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calcola(unsigned int);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    SphericalCorrelations<l,TFLOAT,T> & operator =(const SphericalCorrelations<l,TFLOAT,T> & destra){
        OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::operator = (destra);
        return *this;
    }
    const std::vector<ssize_t> get_shape()const { return {leff,ntypes,ntypes,nbin,(l+1)};}
    const std::vector<ssize_t> get_stride()const {
        auto s=get_shape();
        auto s_old=s[s.size()-1];
        s[s.size()-1]=sizeof (TFLOAT);
        for (int i=s.size()-2;i>=0;--i) {
            auto t=s_old*s[i+1];
            s_old=s[i];
            s[i]=t;
        }
        return s;
    }

    std::string get_columns_description() {return c_descr;}
    inline int index(const int t, const int type1, const int type2,const int ibin=0) const noexcept {
        return (l+1)*(nbin*(ntypes*(ntypes*t + type1) + type2)+ibin);
    }
    using OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::azzera;
    int get_single_type_size() const {
        return (l+1)*(l+1)*nbin*ntypes;
    }
    int get_snap_size()const{
        return get_single_type_size()*natoms;
    }
    int get_final_snap_size() const {
        return get_single_type_size()*ntypes/(l+1);
    }
    inline void calc(int timestep, TFLOAT * result, TFLOAT * workspace, TFLOAT * cheby) const;
    void corr_sh_calc(const TFLOAT * sh1, const TFLOAT *sh2, TFLOAT * aveTypes, TFLOAT * aveWork1, int sh_snap_size , int sh_final_size, int *avecont) const noexcept;

private:
    using OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::lista;
    using OperazioniSuLista<SphericalCorrelations<l,TFLOAT,T>,TFLOAT>::lunghezza_lista;
    T & t;

    inline int index_wrk(const int iatom,const int jtype,const int ibin=0) const noexcept {
        return (l+1)*(l+1)*(nbin*(ntypes*iatom+jtype)+ibin);
    }





    TFLOAT rmin, rmax,dr;
    unsigned int nbin, tmax,nthreads,skip,leff,ntimesteps,ntypes,natoms,buffer_size;
    bool debug;
    std::string c_descr;

};

#endif // SPHERICALCORRELATIONS_H
