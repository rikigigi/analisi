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
        sfac = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*3*nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2);
        //lista = new double[nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2*3];
        lunghezza_lista=nk*nk*nk*tipi_atomi*(tipi_atomi+1)/2*3;
        lista = new double[lunghezza_lista];
       // lista = (double *) fftw_malloc(sizeof(double)*nk*nk*nk*tipi_atomi*(tipi_atomi+1)/2*3);
    }
    ntimesteps=numeroTimestepsPerBlocco;
    std::cerr<< "Sto usando il CorrSpaziale Davide"<<std::endl;
}


void CorrelatoreSpaziale::s_fac_k(double  k[3], unsigned int i_t, fftw_complex * out1 ) {
    for (unsigned int i_at=0;i_at<t->get_natoms();i_at++){
        double * pos = t->posizioni(i_t,i_at);
        double * vel = t->velocita(i_t,i_at);
        unsigned int type=t->get_type(i_at);
        //std::cerr << type << "\n";
        double arg=k[0]*pos[0]+k[1]*pos[1]+k[2]*pos[2];
        double s=sin(arg),c=cos(arg);
        for (unsigned int icoord=0;icoord<3;icoord++){
            out1[3*type+icoord][0]+=c*vel[icoord];
            out1[3*type+icoord][1]+=s*vel[icoord];
        }
    }
}

void CorrelatoreSpaziale::calcola(unsigned int primo) {

    //size = nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2;
    size = nk*nk*nk*tipi_atomi*(tipi_atomi+1)/2;
    //int size_half=nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2;
    int size_half=nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2;
    int itime = 0;

    for (unsigned int i=0;i<3*nk*nk*(nk/2+1)*tipi_atomi*(tipi_atomi+1)/2;i++) {
        sfac[i][0]=0.0;
        sfac[i][1]=0.0;
    }
    fftw_complex * sfac_t=(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * tipi_atomi*3);
    int n[]={nk,nk,nk};
    //int n[]={1,1,1};
    //int n[]={nk,nk,nk/2+1};
    fftw3 = fftw_plan_dft_c2r_3d ( n[0],n[1],n[2], // lunghezza di ciascuna trasformata
                                   /*  ****** input ******  */
                                   sfac,//traiettoria->velocita_inizio(), // puntatore ai dati
                                   /*  ****** output ******  */
                                   lista, // puntatore all'array di output
                                   FFTW_DESTROY_INPUT
                                   );

 //   fftw3 = fftw_plan_many_dft_c2r(3, // rango della trasformata (3D)
 //                                  n, // lunghezza di ciascuna trasformata
 //                                  tipi_atomi*(tipi_atomi+1)/2, // numero di trasformate
 //                                  /*  ****** input ******  */
 //                                  sfac,//traiettoria->velocita_inizio(), // puntatore ai dati
 //                                  NULL, // i dati sono tutti compatti, non fanno parte di array più grandi
 //                                  1, // la distanza fra due dati successivi
 //                                  nk*nk*nk, // la distanza fra due serie di dati adiacenti
 //                                  /*  ****** output ******  */
 //                                  lista, // puntatore all'array di output
 //                                  NULL,
 //                                  1, // la distanza fra due dati successivi
 //                                  nk*nk*nk, // la distanza fra due serie di dati adiacenti
 //                                  FFTW_DESTROY_INPUT
 //                                  );

    for (unsigned int itimestep=primo;itimestep<ntimesteps+primo;itimestep+=skip){
        itime += 1;
        double * box=t->scatola(itimestep);
        //ciclo su kx,ky e kz. kz va da 0 a kz/2+1
        double k[3]={0.0};
        unsigned int ii_tipi=0;
        for (unsigned int i=0;i<tipi_atomi;i++){
            for (unsigned int j=0;j<=i;j++){
                int ii_tipi_matrix=ii_tipi*nk*nk*(nk/2+1);
                for (int kx=0;kx<nk;kx++){
                    unsigned int ikx=ii_tipi_matrix+kx*nk*(nk/2+1);  
                    k[0]=PI/(box[1]-box[0])*kx;
                    for (int ky=0;ky<nk;ky++) {
                        unsigned int iky=ikx+ky*(nk/2+1);  
                        k[1]=PI/(box[1]-box[0])*ky;
                        for (int kz=0;kz<nk/2+1;kz++) {
                            for (unsigned int i=0;i<tipi_atomi*3;i++) {
                                sfac_t[i][0]=0.0;
                                sfac_t[i][1]=0.0;
                            }
                            unsigned int ik=iky+kz;
                            //k[2]=PI/(box[1]-box[0])*ky; ///// MERDA IL COPIA/INCOLLA
                            k[2]=PI/(box[1]-box[0])*kz;
                          //  sfac_t[0][0]=0.1 ;/////// TOGLIERE
                            s_fac_k(k,itimestep,sfac_t);

                            sfac[ik+ii_tipi_matrix][0]+=       ( sfac_t[i*3][0]*  sfac_t[j*3][0]  + sfac_t[i*3][1]*  sfac_t[j*3][1] -sfac[ik+ii_tipi_matrix][0])/(itime);
                            sfac[ik+ii_tipi_matrix][1]+=       (-sfac_t[i*3][0]*  sfac_t[j*3][1] + sfac_t[i*3][1]*  sfac_t[j*3][0] -sfac[ik+ii_tipi_matrix][1])/(itime);
                            sfac[size_half  +ik+ii_tipi_matrix][0]+=  ( sfac_t[i*3+1][0]*sfac_t[j*3+1][0] + sfac_t[i*3+1][1]*sfac_t[j*3+1][1]  - sfac[  size_half+ik+ii_tipi_matrix][0])/(itime);
                            sfac[size_half  +ik+ii_tipi_matrix][1]+=  (-sfac_t[i*3+1][0]*sfac_t[j*3+1][1] + sfac_t[i*3+1][1]*sfac_t[j*3+1][0] -  sfac[  size_half+ik+ii_tipi_matrix][1])/(itime);
                            sfac[2*size_half+ik+ii_tipi_matrix][0]+=( sfac_t[i*3+2][0]*sfac_t[j*3+2][0] + sfac_t[i*3+2][1]*sfac_t[j*3+2][1]-    sfac[2*size_half+ik+ii_tipi_matrix][0])/(itime);
                            sfac[2*size_half+ik+ii_tipi_matrix][1]+=(-sfac_t[i*3+2][0]*sfac_t[j*3+2][1] + sfac_t[i*3+2][1]*sfac_t[j*3+2][0]-     sfac[2*size_half+ik+ii_tipi_matrix][1])/(itime);
     //                         sfac[ik+ii_tipi_matrix][0]+=               cos(2*PI/(nk)*ik);
     //                         sfac[ik+ii_tipi_matrix][1]+=               sin(2*PI/(nk)*ik); 
     //                         sfac[size_half  +ik+ii_tipi_matrix][0]+=   cos(2*PI/(nk)*ik); 
     //                         sfac[size_half  +ik+ii_tipi_matrix][1]+=   sin(2*PI/(nk)*ik);
     //                         sfac[2*size_half+ik+ii_tipi_matrix][0]+=   cos(2*PI/(nk)*ik); 
     //                         sfac[2*size_half+ik+ii_tipi_matrix][1]+=   sin(2*PI/(nk)*ik);

                        }
                    }
                }
            ii_tipi += 1;  
            }

        }
     }
    //fare trasformata di fourier discreta (copia da spettrovibrazionale.cpp:84)



    fftw_free( sfac_t);
   for (unsigned int i=0;i<tipi_atomi*(tipi_atomi+1)/2;i++){  
    fftw_execute_dft_c2r(fftw3,&sfac[              i*nk*nk*(nk/2+1)], &lista[i*nk*nk*nk]);
    fftw_execute_dft_c2r(fftw3,&sfac[  size_half + i*nk*nk*(nk/2+1)], &lista[i*nk*nk*nk + size]);
    fftw_execute_dft_c2r(fftw3,&sfac[2*size_half + i*nk*nk*(nk/2+1)], &lista[i*nk*nk*nk + 2*size]);
   }

}



CorrelatoreSpaziale::~CorrelatoreSpaziale(){
    fftw_free(sfac);
}

CorrelatoreSpaziale & CorrelatoreSpaziale::operator =(const CorrelatoreSpaziale & destra) {
    OperazioniSuLista<CorrelatoreSpaziale>::operator =(destra);

}


double CorrelatoreSpaziale::corr(unsigned int rx, unsigned int ry, unsigned int rz, unsigned int itype,unsigned int idim) {
    if (itype >=tipi_atomi*(tipi_atomi+1)/2) {
            std::cerr << " Errore: indice dei tipi fuori dal range in " __FILE__ "\n";
            abort();
    }
    if (rx>=nk || ry>=nk || rz>=nk) {
            std::cerr << " Errore: indice spaziale fuori dal range in " __FILE__ "\n";
            abort();
    }
    if (idim >= 3) {
        std::cerr << " Errore: indice della coordinata fuori dal range in " __FILE__ "\n";
    }
    return lista[idim*nk*nk*nk*tipi_atomi*(tipi_atomi+1)/2+ rx*nk*nk*tipi_atomi*(tipi_atomi+1)/2+
            ry*nk*tipi_atomi*(tipi_atomi+1)/2+
            rz*tipi_atomi*(tipi_atomi+1)/2+
            itype
            ];
}
