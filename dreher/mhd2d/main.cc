// -*- C++ -*-
// $Id: main.cc,v 1.7 2011/07/21 09:08:56 jd Exp $

#include <cstdio>
#include "global.h"

using namespace mhd2d;

// program entry point
int main(int argc, char** argv)
{
  printf("Hello, MHD 2d!\ndx=%g dy=%g dt=%g nStep=%d\n", dx, dy, dt, nStep);

  // setup
  allocateArrays();
  initFields();

  // write initial data
  writeVTKFile(0);

  // integration loop
  for(int step=1; step<=nStep; step++){
    doStep(); // see timeStep.cc
    
    if(0 == (step%outputStep)){ // % is modulo-operator in C
      writeVTKFile(step);
    }
  }
}
