/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef TESTTRAIETTORIA_H
#define TESTTRAIETTORIA_H

#include "trajectory.h"

class Trajectory;

class TestTraiettoria : public Trajectory
{
public:
    TestTraiettoria(std::string filename);

};

#endif // TESTTRAIETTORIA_H
