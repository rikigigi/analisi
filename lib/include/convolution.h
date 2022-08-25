/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <functional>
template < typename T >class Convolution
{
public:
    /**
      * convoluzione della funzione f con un array di tipo T. Per il tipo T devono essere defiie le operazioni di moltiplicazione con il tipo double.
    **/
    Convolution(std::function<double(const double &)> f, int n, double start, double stop, int center);
    ~Convolution();
    void calculate(T* in , T *out,int n,int skip);
    void calculate(T* in , T *out,int n);
private:
    T * fc;
    unsigned int n_fc,centro;


};

#endif // CONVOLUTION_H
