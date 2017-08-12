#ifndef CONVOLUTION_H
#define CONVOLUTION_H


template < typename T >class Convolution
{
public:
    Convolution(T (*f)(const T &), int n, T start, T stop, int center);
    ~Convolution();
    void calcola(T* in , T *out,int n);
private:
    T * fc;
    unsigned int n_fc,centro;


};

#endif // CONVOLUTION_H
