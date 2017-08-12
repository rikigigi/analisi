#include "convolution.h"
#include <iostream>

template <typename T> Convolution<T>::Convolution(T (*f)(const T &),int n,T start,T stop,int center)
{
    if (n<=0 || center <0 || center >= n) {
        std::cerr<< "Errore: richiesta una lunghezza "<<n<<" ed un centro "<<center<<" per la convoluzione!\n";
        abort();
    }
    n_fc=n;
    fc=new T[n];
    for (int i=0;i<n;i++){
        fc[i]=f(start+(stop-start)*double(i)/double(n-1));
    }

}

template <typename T> Convolution<T>::calcola(T *in,T*out,int n) {
    for (unsigned int in_index=0;in_index<n;in_index++) {
        out[in_index]=0;
        for (unsigned int i=0;i<n_fc;i++ ) {
            int idx=in_index- centro+i;
            if (idx <0)
                idx=0;
            else if (idx>=n)
                idx=n-1;
            out[in_index]+=in[idx]*fc[i];
        }
    }
}
