// -*- C++ -*-
// $Id: main.cc,v 1.6 2011/06/28 07:44:29 jd Exp $

#include <cstdio>
#include "global.h"

using namespace mhd2d;

int main(int argc, char** argv)
{
  printf("Hello, MHD 2d!\ndx=%g dy=%g dt=%g nStep=%d\n", dx, dy, dt, nStep);
  allocateArrays();
  initFields();

  writeVTKFile(0);
  for(int step=1; step<=nStep; step++){
    doStep();
    
    if(0 == (step%outputStep)){
      writeVTKFile(step);
    }
  }
}
