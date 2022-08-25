/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef ISTOGRAMMAATOMIRAGGIO_H
#define ISTOGRAMMAATOMIRAGGIO_H

#include <map>

class Trajectory;

class IstogrammaAtomiRaggio
{
public:
    IstogrammaAtomiRaggio(Trajectory * t,double r,unsigned int skip=1,unsigned int nthreads=0);
    void reset(const unsigned int numeroTimestepsPerBlocco);
    void calculate(unsigned int);
    unsigned int nExtraTimesteps(unsigned int n_b) {return 0;}
    std::map<unsigned int, unsigned int> * get_hist() {return hist;};
private:
    double r2;
    unsigned int skip,nthreads,ntypes,natoms,numeroTimestepsBlocco;
    Trajectory * traiettoria;
    std::map<unsigned int, unsigned int> * hist;
};

#endif // ISTOGRAMMAATOMIRAGGIO_H
