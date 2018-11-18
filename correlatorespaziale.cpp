#include "correlatorespaziale.h"
#include "traiettoria.h"
#include "operazionisulista.h"

#include "config.h"
#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif
CorrelatoreSpaziale::CorrelatoreSpaziale(Traiettoria *t,
                                          unsigned int nk,
                                          TFLOAT sigma2,
                                          unsigned int nthreads=0,
                                          unsigned int skip=1,
                                          bool debug=false) :
    sfac(0), nk(nk), sigma2(sigma2),nthreads(nthreads),debug(debug),skip(skip),tipi_atomi(0),t(t)
{


}

template  <class TFLOAT>
void CorrelatoreSpaziale<TFLOAT>::reset(const unsigned int numeroTimestepsPerBlocco) {
     tipi_atomi=traiettoria->get_ntypes();
    if (sfac==0){
        sfac = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2);
        lista = new double[nk*nk*nk];
    }
    ntimesteps=numeroTimestepsPerBlocco;
}

template <class TFLOAT>
fftw_complex[3] CorrelatoreSpaziale<TFLOAT>::s_fac_k(TFLOAT k[3],unsigned int i_t) {
    fftw_complex res[3]={0.0};
    for (unsigned int i_at=0;i_at<t->get_natoms();i_at++){
        double * pos = t->posizioni(i_t,i_at);
        double * vel = t->velocita(i_t,i_at);
        TFLOAT arg=k[0]*pos[0]+k[1]*pos[1]+k[2]*pos[2];
        TFLOAT s=sin(arg),c=cos(arg);
        for (unsigned int icoord=0;icoord<3;icoord++){
            res[icoord][0]+=c*vel[icoord];
            res[icoord][1]+=s*vel[icoord];
        }
    }
    return res;
}

template <class TFLOAT>
void CorrelatoreSpaziale<TFLOAT>::s_fac_k(TFLOAT  k[3], unsigned int i_t,fftw_complex * out ) {
    for (unsigned int i_at=0;i_at<t->get_natoms();i_at++){
        double * pos = t->posizioni(i_t,i_at);
        double * vel = t->velocita(i_t,i_at);
        unsigned int type=t->get_type(i_at);
        TFLOAT arg=k[0]*pos[0]+k[1]*pos[1]+k[2]*pos[2];
        TFLOAT s=sin(arg),c=cos(arg);
        for (unsigned int icoord=0;icoord<3;icoord++){
            out[3*type+icoord][0]+=c*vel[icoord];
            out[3*type+icoord][1]+=s*vel[icoord];
        }
    }
}

template <class TFLOAT>
void CorrelatoreSpaziale<TFLOAT>::calcola(unsigned int primo) {
   
   for (unsigned int itimestep=0;itimestep<ntimesteps;itimestep++){
       //ciclo su kx,ky e kz. kz va da 0 a kz/2+1
       TFLOAT k[3]={0.0}
       for (unsigned int kx=0;kx<nk;kx++){
          k[0]=
          for (unsigned int ky=0;kx
       }
   }

}



template <class TFLOAT>
CorrelatoreSpaziale::~CorrelatoreSpaziale(){
    fftw_free(sfac);
}
