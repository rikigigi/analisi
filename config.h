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



#ifndef CONFIG_H
#define CONFIG_H


#define HAVEeigen3EigenDense 1
/* #undef HAVEEigenDense */
#define FFTW_OMP
/* #undef FFTW_TH */
#if defined(FFTW_OMP) || defined(FFTW_TH)
#define FFTW3_THREADS
#endif
#define CMAKE_CXX_COMPILER "/usr/bin/g++"
#define CMAKE_CXX_FLAGS " -std=c++11 -fopenmp"
#define CMAKE_SYSTEM "Linux-4.9.16-gentoo"
#define CMAKE_SYSTEM_PROCESSOR "x86_64"








#endif
