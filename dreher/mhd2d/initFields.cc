// -*- C++ -*-
// $Id: initFields.cc,v 1.7 2011/06/28 07:44:29 jd Exp $

#include <cmath>
#include "global.h"

namespace mhd2d
{

  static real xcenter(int i) { return xmin + (i+0.5)*dx; }
  static real ycenter(int j) { return ymin + (j+0.5)*dy; }

  void initFields()
  {
    for(int j=-2; j<ny+2; j++){
      const real y = ycenter(j);
      for(int i=-2; i<nx+2; i++){
      const real x = xcenter(i);

      rho(i, j) = ((x > 0) && (y<2-x) && (2*y>x-2)) ? 1. : 0.5;

      ux(i, j) = 0.*rho(i, j);
      uy(i, j) = 0.*rho(i, j);

//       bx(i, j) = 1.;
//       by(i, j) = 0.;

      w(i, j) = 0.5*(ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)) / rho(i, j) + 1. * rho(i, j);
      
      }
    }

    boundaryCond();
  }
  
}


