#include "cepstral.h"
#include "convolution.h"
#include <iostream>

Cepstral::Cepstral(std::vector<Serie> in_) : correnti(in_),fplan(0),size_in(0)
{

    if (correnti.size()!=0)
        size_in=correnti.front().length;

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
                correnti[i].fft_idx[k1]=fft.size();
                in_src.push_back(correnti[i].serie[k1]);
                fft.push_back((std::complex<double> *)fftw_alloc_complex(correnti[i].length/2+1));
                fft_conv.push_back((std::complex<double> *)fftw_alloc_complex(correnti[i].length/2+1));
            }
        }
    }

    in=fftw_alloc_real(size_in);

    // inizializza la trasformata di fourier veloce
    fplan=fftw_plan_dft_r2c_1d(size_in,in,reinterpret_cast<fftw_complex*>(fft.front()),FFTW_MEASURE);
}


void Cepstral::calcola(unsigned int conv_n){

    //copia i dati in ingresso negli array e calcola le trasformate di fourier convolute

    Convolution<std::complex<double> > convoluzione(std::function< std::complex<double> (const std::complex<double>  &)> ([&conv_n](const std::complex<double> & x)->std::complex<double>{
        return exp(-x*x/(2.0*conv_n*conv_n));
    }), (conv_n*6+1),-3*conv_n,3*conv_n,3*conv_n);
    for (unsigned int i=0;i<in_src.size();i++){
        for (unsigned int j=0;j<size_in;j++){
            in[j]=in_src[i][j];
        }
        fftw_execute_dft_r2c(fplan,in,reinterpret_cast<fftw_complex*>(fft.at(i)));
        //fai la convoluzione
        convoluzione.calcola(fft.at(i),fft_conv.at(i),size_in/2+1);

    }

    // esegue il logaritmo complesso di J1^* J2






}
