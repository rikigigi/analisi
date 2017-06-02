/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
  *
**/



#ifndef FIT_H
#define FIT_H
#define N__lin Dynamic

#include "config.h"
#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif
#include <fstream>

#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939
#endif
using namespace Eigen;

template <unsigned int N_lin> class fit{
private:
    static double genX(const unsigned int &i, const unsigned int &j,const double &tau,const double *x);
    static void ols_multidim(const Matrix<double,Dynamic,N__lin> &X , const Matrix<double,Dynamic,1> &Y, Matrix<double, N__lin, 1> &B);
    static double aus_func(const double &t, const Matrix<double,N__lin,1> &B, const double &tau);
    static double chi2(const unsigned int &N, const double *x, const double *y, const Matrix<double,N__lin,1> &B,const double &tau);
    static double chi2_derivat (const unsigned int &N,const double *x,const double *y, const Matrix<double,N__lin,1> &B,const double &tau);
    static bool output;
    static char * fout;
    static unsigned int fitcont;
    static unsigned int fitskip;
    static double fit_func(const double &t,const Matrix<double,N__lin,1> &B, const double &tau);
public:
    static void set_fitskip(unsigned int);
    static void set_output(char * );
    static void fit_lin_exp(const unsigned int &,const double *,const double*, double *,double&,bool teor_lin=false,unsigned int natoms=1024);
    static double get_fit_func(const unsigned int & N, const double *x, double * y, const double * b, const double & tau);
};

#endif
