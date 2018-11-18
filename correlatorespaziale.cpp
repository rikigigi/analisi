#include "correlatorespaziale.h"
#include "traiettoria.h"
#include "operazionisulista.h"

#include "config.h"
#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628

template  <class TFLOAT>
CorrelatoreSpaziale<TFLOAT>::CorrelatoreSpaziale(Traiettoria *t,
                                                 unsigned int nk,
                                                 TFLOAT sigma2,
                                                 unsigned int nthreads,
                                                 unsigned int skip,
                                                 bool debug) :
    sfac(0), nk(nk), sigma2(sigma2),nthreads(nthreads),debug(debug),skip(skip),tipi_atomi(0),t(t)
{


}

template  <class TFLOAT>
fftw_plan CorrelatoreSpaziale<TFLOAT>::fftw3;

template  <class TFLOAT>
void CorrelatoreSpaziale<TFLOAT>::reset(const unsigned int numeroTimestepsPerBlocco) {
    tipi_atomi=t->get_ntypes();
    if (sfac==0){
        sfac = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2);
        lista = new TFLOAT[nk*nk*nk*tipi_atomi*(tipi_atomi+1)/2];
    }
    ntimesteps=numeroTimestepsPerBlocco;
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

    for (unsigned int itimestep=primo;itimestep<ntimesteps+primo;itimestep++){
        double * box=t->scatola(itimestep);
        for (unsigned int i=0;i<nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2;i++) {
            sfac[i][0]=0.0;
            sfac[i][1]=0.0;
        }
        //ciclo su kx,ky e kz. kz va da 0 a kz/2+1
        TFLOAT k[3]={0.0};
        for (int kx=0;kx<nk;kx++){
            unsigned int ikx=kx*nk*(nk/2+1)*3*(tipi_atomi+1)*tipi_atomi/2;
            k[0]=PI/(box[1]-box[0])*kx;
            for (int ky=0;ky<nk;ky++) {
                unsigned int iky=ky*(nk/2+1)*3*(tipi_atomi+1)*tipi_atomi/2;
                k[1]=PI/(box[1]-box[0])*ky;
                for (int kz=0;kz<nk/2+1;kz++) {
                    unsigned int ik=ikx+iky+kz*3*(tipi_atomi+1)*tipi_atomi/2;
                    k[2]=PI/(box[1]-box[0])*ky;
                    s_fac_k(k,itimestep,&sfac[ik]);
                }

            }
        }
        //fare trasformata di fourier discreta (copia da spettrovibrazionale.cpp:84)

        //calcolare i vari moduli quadri

        //media online
        //cosa mettere in lista[nk*nk*nk*tipi_atomi*(tipi_atomi+1)/2]???
    }

}



template <class TFLOAT>
CorrelatoreSpaziale<TFLOAT>::~CorrelatoreSpaziale(){
    fftw_free(sfac);
}
