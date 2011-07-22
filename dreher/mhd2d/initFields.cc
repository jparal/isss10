// -*- C++ -*-
// $Id: initFields.cc,v 1.11 2011/07/21 14:49:48 jd Exp $

#include <cmath>
#include "global.h"

namespace mhd2d
{

  // position of cell center computed from cell index
  static real xcenter(int i) { return xmin + (i+0.5)*dx; }
  static real ycenter(int j) { return ymin + (j+0.5)*dy; }

  // initial conditions
  void initFields()
  {
    // loop over all cells, incl. ghost cells (2 layers)
    for(int j=-2; j<ny+2; j++){
      const real y = ycenter(j);

      for(int i=-2; i<nx+2; i++){
        // x-center of current cells
        const real x = xcenter(i);
        
        // Alfven wave parameters
        const real A = 1e-4;
        const real k = 2*M_PI / (xmax-xmin);
        const real rho0 = 1.;

        // set the field variables
        rho(i, j) = rho0;
        ux(i, j) = 0.;
        uy(i, j) = A * rho0 * sin(k*x);
        bx(i, j) = 1.;
        by(i, j) = A * sqrt(rho0) * sin(k*x);
        
        const real p = 1.;

        // internal energy is computed from p according to its definition
        e(i, j) = 
          + p / (gamma-1.)
          + 0.5*( sqr(ux(i, j)) + sqr(uy(i, j)) ) / rho(i, j)
          + 0.5*( sqr(bx(i, j)) + sqr(by(i, j)) );
      }
    }
    
    // finally, apply b.c. to ensure periodicity
    boundaryCond();
  }
  
}


