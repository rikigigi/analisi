/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef GOFRT_H
#define GOFRT_H

#include "operazionisulista.h"
#include "calcolamultithread.h"


namespace Gofrt_Flags {
constexpr int FLAGS = CalcolaMultiThread_Flags::PARALLEL_SPLIT_ATOM |
        CalcolaMultiThread_Flags::PARALLEL_LOOP_TIME |
        CalcolaMultiThread_Flags::PARALLEL_LOOP_AVERAGE |
        CalcolaMultiThread_Flags::CALL_DEBUG_ROUTINE |
        CalcolaMultiThread_Flags::CALL_CALC_INIT;
}

template <class TFLOAT, class T>
class Gofrt : public OperazioniSuLista<Gofrt<TFLOAT,T>,TFLOAT>, public CalcolaMultiThread<Gofrt<TFLOAT,T>, Gofrt_Flags::FLAGS  >
{
public:
    using This = Gofrt<TFLOAT,T>;
    using CalcolaMultiThread_T = CalcolaMultiThread<This, Gofrt_Flags::FLAGS>;
    using CalcolaMultiThread_T::FLAGS;
    using OperazioniSuLista_T = OperazioniSuLista<This,TFLOAT>;

    Gofrt(T *t,
          TFLOAT rmin,
          TFLOAT rmax,
          unsigned int nbin,
          unsigned int tmax=0,
          unsigned int nthreads=0,
          unsigned int skip=1,
	  unsigned int every=1,
          bool debug=false);
    ~Gofrt();
    void reset(const unsigned int numeroTimestepsPerBlocco);
    unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b);
    This & operator =(const This & destra);
    std::vector<ssize_t> get_shape();
    std::vector<ssize_t> get_stride();
    std::string get_columns_description() {return c_descr;}
    using OperazioniSuLista_T::azzera;

    void calc_init(int);
    void calc_single_th(int,int,int,int,int,int);
    void calc_end();

private:
    using OperazioniSuLista_T::lista;
    using OperazioniSuLista_T::lunghezza_lista;
    TFLOAT * th_data;
    TFLOAT rmin,rmax,rmax2,rmin2,dr,incr;
    bool debug;
    T * traiettoria;
    unsigned int nbin,lmax;
    using CalcolaMultiThread_T::ntimesteps;
    using CalcolaMultiThread_T::skip;
    using CalcolaMultiThread_T::nthreads;
    using CalcolaMultiThread_T::leff;
    int gofr_idx(unsigned int ts, unsigned int itype=0, unsigned int r=0) {
        unsigned int idx= ts   * traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)*nbin
                         +nbin * itype
                         +r;
        return idx;
    }
    TFLOAT * gofr(unsigned int ts, unsigned int itype=0, unsigned int r=0){
    unsigned int idx= gofr_idx(ts,itype,r);
    if (idx >= lunghezza_lista) {
        std::cerr << "Errore: indice fuori dal range!\n";
        abort();
    }
    return &lista[idx];
}
    std::string c_descr;
    unsigned int get_itype(unsigned int & type1,unsigned int & type2) const {
        /*
        * xxxxx  ntypes*(ntypes+1)/2 - (m+2)*(m+1)/2 +
        * xxxxo  + altra coordinata (che deve essere la più grande)
        * xxxoo  = indice della coppia nella memoria
        * xxooo
        * xoooo
        *
        * xx
        * x
       */
       if (type2<type1) {
           unsigned int tmp=type2;
           type2=type1;
           type1=tmp;
       }
       return traiettoria->get_ntypes()*(traiettoria->get_ntypes()+1)/2
               -(type2+1)*(type2+2)/2 +type1;
   }


};

#endif // GOFRT_H
