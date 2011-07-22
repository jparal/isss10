// -*- C++ -*-
// $Id: boundary.cc,v 1.3 2011/06/28 07:44:29 jd Exp $

#include "global.h"

namespace mhd2d
{

  static void boundaryX(Array& f)
  {
    for(int j=0; j<ny; ++j){
      f(-2, j) = f(nx-2, j); f(-1, j) = f(nx-1, j);
      f(nx, j) = f(0, j); f(nx+1, j) = f(1, j);
    }
  }

  static void boundaryY(Array& f)
  {
    for(int i=-2; i<nx+2; ++i){
      f(i, -2) = f(i, ny-2); f(i, -1) = f(i, ny-1);
      f(i, ny) = f(i, 0); f(i, ny+1) = f(i, 1);
    }
  }

  void boundaryCond()
  {
    for(int iField=0; iField<N_FIELDS; iField++){
      boundaryX(fields[iField]);
      boundaryY(fields[iField]);
    }
  }
  
}


