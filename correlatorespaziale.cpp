#include "correlatorespaziale.h"
#include "traiettoria.h"
#include "operazionisulista.h"
#include <cmath>

#include "config.h"
#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628

CorrelatoreSpaziale::CorrelatoreSpaziale(Traiettoria *t,
                                                 unsigned int nk,
                                                 double sigma2,
                                                 unsigned int nthreads,
                                                 unsigned int skip,
                                                 bool debug) :
    sfac(0), nk(nk), sigma2(sigma2),nthreads(nthreads),debug(debug),skip(skip),tipi_atomi(0),t(t)
{


}


fftw_plan CorrelatoreSpaziale::fftw3;

void CorrelatoreSpaziale::reset(const unsigned int numeroTimestepsPerBlocco) {
    tipi_atomi=t->get_ntypes();
    if (sfac==0){
        sfac = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2);
        lista = new double[nk*nk*nk*tipi_atomi*(tipi_atomi+1)/2];
    }
    ntimesteps=numeroTimestepsPerBlocco;
}


void CorrelatoreSpaziale::s_fac_k(double  k[3], unsigned int i_t,fftw_complex * out ) {
    for (unsigned int i_at=0;i_at<t->get_natoms();i_at++){
        double * pos = t->posizioni(i_t,i_at);
        double * vel = t->velocita(i_t,i_at);
        unsigned int type=t->get_type(i_at);
        double arg=k[0]*pos[0]+k[1]*pos[1]+k[2]*pos[2];
        double s=sin(arg),c=cos(arg);
        for (unsigned int icoord=0;icoord<3;icoord++){
            out[3*type+icoord][0]+=c*vel[icoord];
            out[3*type+icoord][1]+=s*vel[icoord];
        }
    }
}

void CorrelatoreSpaziale::calcola(unsigned int primo) {

    size = nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2;
    int itime = 0;

    for (unsigned int i=0;i<nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2;i++) {
        sfac[i][0]=0.0;
        sfac[i][1]=0.0;
    }
    fftw_complex * sfac_t=new fftw_complex[tipi_atomi];
    for (unsigned int itimestep=primo;itimestep<ntimesteps+primo;itimestep++){
        itime += 1;
        for (unsigned int i=0;i<tipi_atomi;i++) {
            sfac_t[i][0]=0.0;
            sfac_t[i][1]=0.0;
        }
        double * box=t->scatola(itimestep);
        //ciclo su kx,ky e kz. kz va da 0 a kz/2+1
        double k[3]={0.0};
        for (int kx=0;kx<nk;kx++){
            unsigned int ikx=kx*nk*(nk/2+1)*3*(tipi_atomi+1)*tipi_atomi/2;
            k[0]=PI/(box[1]-box[0])*kx;
            for (int ky=0;ky<nk;ky++) {
                unsigned int iky=ky*(nk/2+1)*3*(tipi_atomi+1)*tipi_atomi/2;
                k[1]=PI/(box[1]-box[0])*ky;
                for (int kz=0;kz<nk/2+1;kz++) {
                    unsigned int ik=ikx+iky+kz*3*(tipi_atomi+1)*tipi_atomi/2;
                    k[2]=PI/(box[1]-box[0])*ky;
                    s_fac_k(k,itimestep,sfac_t);
                    for (unsigned int i=0;i<tipi_atomi;i++)
                        for (unsigned int j=0;j<=i;j++){
                            unsigned int m = (tipi_atomi+1)*tipi_atomi/2 - (j+1)*(j+2)/2 + i;

                            sfac[ik+m][0]+=(sfac_t[i][0]*sfac_t[j][0]-sfac[ik+m][0])/(itime);
                            sfac[ik+m][1]+=(-sfac_t[i][0]*sfac_t[j][1]-sfac[ik+m][1])/(itime);
                        }
                }

            }
        }



    }


    //fare trasformata di fourier discreta (copia da spettrovibrazionale.cpp:84)
    int n[]={nk,nk,nk};
    fftw3 = fftw_plan_many_dft_c2r(3, // rango della trasformata (3D)
                                   n, // lunghezza di ciascuna trasformata
                                    tipi_atomi*(tipi_atomi+1)/2, // numero di trasformate
                                   /*  ****** input ******  */
                                   sfac,//traiettoria->velocita_inizio(), // puntatore ai dati
                                   NULL, // i dati sono tutti compatti, non fanno parte di array più grandi
                                   (tipi_atomi+1)*tipi_atomi/2, // la distanza fra due dati successivi
                                   1, // la distanza fra due serie di dati adiacenti
                                   /*  ****** output ******  */
                                   lista, // puntatore all'array di output
                                   NULL,
                                   (tipi_atomi+1)*tipi_atomi/2, // la distanza fra due dati successivi
                                   1, // la distanza fra due serie di dati adiacenti
                                   FFTW_PRESERVE_INPUT
                                   );



    delete [] sfac_t;

}



CorrelatoreSpaziale::~CorrelatoreSpaziale(){
    fftw_free(sfac);
}

CorrelatoreSpaziale & CorrelatoreSpaziale::operator =(const CorrelatoreSpaziale & destra) {
    OperazioniSuLista<CorrelatoreSpaziale>::operator =(destra);

}


double CorrelatoreSpaziale::corr(unsigned int rx, unsigned int ry, unsigned int rz, unsigned int itype) {
    if (itype >=tipi_atomi*(tipi_atomi+1)/2) {
            std::cerr << "Errore: indice dei tipi fuori dal range in "__FILE__"\n";
            abort();
    }
    if (rx>=nk || ry>=nk || rz>=nk) {
            std::cerr << "Errore: indice spaziale fuori dal range in "__FILE__"\n";
            abort();
    }
    return lista[ rx*nk*nk*tipi_atomi*(tipi_atomi+1)/2+
            ry*nk*tipi_atomi*(tipi_atomi+1)/2+
            rz*tipi_atomi*(tipi_atomi+1)/2+
            itype
            ];
}
