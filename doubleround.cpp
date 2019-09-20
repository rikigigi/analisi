/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "doubleround.h"
#include <stdint.h>

DoubleRound::DoubleRound()
{

}


void DoubleRound::round(double *d) {
    uint64_t t;
    //                                    |      ||      ||      ||      ||      ||      ||      ||      |
    t = (  * (uint64_t*) d &  (uint64_t)0b1111111111111111111111111111111111111111111111111111111100000000) ;
    *d= * (double*) &t;
}
