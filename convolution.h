#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <functional>
template < typename T >class Convolution
{
public:
    Convolution(std::function<T (const T  &)> f, int n, T start, T stop, int center);
    ~Convolution();
    void calcola(T* in , T *out,int n,int skip);
    void calcola(T* in , T *out,int n);
private:
    T * fc;
    unsigned int n_fc,centro;


};

#endif // CONVOLUTION_H
