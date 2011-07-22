// -*- C++ -*-
// $Id: boundary.cc,v 1.4 2011/07/21 09:08:56 jd Exp $

#include "global.h"

// apply boundary conditions
namespace mhd2d
{
  
  // periodic b.c. for 2 ghost cell layers in x
  static void boundaryX(Array& f)
  {
    for(int j=0; j<ny; ++j){
      f(-2, j) = f(nx-2, j); f(-1, j) = f(nx-1, j);
      f(nx, j) = f(0, j); f(nx+1, j) = f(1, j);
    }
  }
  // and y
  static void boundaryY(Array& f)
  {
    for(int i=-2; i<nx+2; ++i){
      f(i, -2) = f(i, ny-2); f(i, -1) = f(i, ny-1);
      f(i, ny) = f(i, 0); f(i, ny+1) = f(i, 1);
    }
  }

  // this is called from outside
  void boundaryCond()
  {
    for(int iField=0; iField<N_FIELDS; iField++){
      // to implement fixed boundaries in x and/or y, comment out the 
      // appropriate lines
      boundaryX(fields[iField]);
      boundaryY(fields[iField]);
    }
  }
  
}


