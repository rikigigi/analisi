/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "convolution.h"
#include <iostream>

template <typename T> Convolution<T>::Convolution(std::function<double (const double  &)> f,int n,double start,double stop,int center)
{
    if (n<=0 || center <0 || center >= n) {
        std::cerr<< "Errore: richiesta una lunghezza "<<n<<" ed un centro "<<center<<" per la convoluzione!\n";
        abort();
    }
    centro=center;
    double norm=0;
    n_fc=n;
    fc=new T[n];
    for (int i=0;i<n;i++){
        fc[i]=f(start+(stop-start)*double(i)/double(n-1));
        norm+=fc[i];
    }

    for (int i=0;i<n;i++){
        fc[i]/=norm;
#ifdef DEBUG
    std::cerr << fc[i]<<"\n";
#endif
    }

}

template <typename T> void Convolution<T>::calcola(T *in, T*out, int n, int skip) {
    for (unsigned int in_index=0;in_index<n;in_index++) {
        out[in_index]=0;
        for (unsigned int i=0;i<n_fc;i++ ) {
            int idx=in_index- centro+i;
            if (idx <0)
                idx=-idx;
            if (idx>=n)
                idx=n-2 - idx%n;
            out[in_index]+=in[idx*skip]*fc[i];
        }
    }
}

template <typename T> void Convolution<T>::calcola(T *in, T*out, int n) {
    for (unsigned int in_index=0;in_index<n;in_index++) {
        out[in_index]=0;
        for (unsigned int i=0;i<n_fc;i++ ) {
            int idx=in_index- centro+i;
            if (idx <0)
                idx=-idx;
            if (idx>=n)
                idx=n-2 - idx%n;
            out[in_index]+=in[idx]*fc[i];
        }
    }
}

template <typename T> Convolution<T>::~Convolution(){
    delete [] fc;
}


template class Convolution<double>;
