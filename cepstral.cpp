#include "cepstral.h"
#include "fftw3.h"
#include "convolution.h"

Cepstral::Cepstral(std::vector<Serie> in_) : correnti(in_),fplan(0),size_in(0)
{

    if (correnti.size()!=0)
        size_in=correnti.first().length;

    for (unsigned int i=0;i<correnti.size();i++){
        if(correnti[i].length!=size_in){
            std::cerr << "Non implementato: lunghezza delle correnti di input diverse! ("<<size_in<<" e "<< correnti[i].length <<")";
            abort();
        }
        //controlla eventuali duplicati e alloca gli array nuovi, se necessario
        correnti[i].fft_idx[0]=-1;
        correnti[i].fft_idx[1]=-1;
        for (unsigned int j=0;j<i;j++){
            for (unsigned int k1=0;k1<2;k1++)
                for (unsigned int k2=0;k2<2;k2++)
                    if (correnti[i].serie[k2]<0 && correnti[j].serie[k1]==correnti[i].serie[k2]){
                        correnti[i].fft_idx[k2]=correnti[j].fft_idx[k1];
                    }
        }

        for (unsigned int k1=0;k1<2;k1++){
            if (correnti[i].fft_idx[k1]<0){
                assert(in.size()==fft.size());
                correnti[i].fft_idx[k1]=in.size();
                in_src.append(correnti[i].serie[k1]);
                fft.append(fftw_alloc_complex(correnti[i].length/2+1));
                fft_conv.append(fftw_alloc_complex(correnti[i].length/2+1));
            }
        }
    }

    in=fftw_alloc_real(size_in);

    // inizializza la trasformata di fourier veloce
    fplan=fftw_plan_dft_r2c_1d(size_in,in,fft.first(),FFTW_MEASURE);
}


void Cepstral::calcola(){

    //copia i dati in ingresso negli array e calcola le trasformate di fourier

    for (unsigned int i=0;i<in_src.size();i++){
        for (unsigned int j=0;j<size_in;j++){
            in[j]=in_src[i][j];
        }
        fftw_execute_dft_r2c(fplan,in,fft.at(i));
        //fai la convoluzione

    }






}
