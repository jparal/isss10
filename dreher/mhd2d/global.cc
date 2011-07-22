// -*- C++ -*-
// $Id: global.cc,v 1.15 2011/07/21 14:50:03 jd Exp $

#include <cmath>
#include "global.h"

// Modify parameter values here and recompile with "make"
namespace mhd2d
{
  // domain ranges
  const real xmin = 0.;
  const real ymin = 0.;
  const real xmax = 2.*M_PI; // two pi
  const real ymax = 2.*M_PI;
  
  // No. of grid cells
  const int nx = 100;
  const int ny = 10;

  // No. of time steps and outputs
  const int nStep = 300;
  const int outputStep = 20;
  const char* outFileStem = "data/MHD";

  // do not change dx/dy!
  const real dx = (xmax-xmin) / nx;
  const real dy = (ymax-ymin) / ny;
  
  // time step size.
  // Can also be set explicitly, e.g. 
  //  const real dt = 1e-2 ;
  const real dt = .5 * fmin(dx, dy);

  // adiabatic index
  const real gamma = 5./3.;

  // do not change these lines
  const real dxi = 1. / dx;
  const real dyi = 1. / dy;
  Array fields[N_FIELDS];
  Array fields_old[N_FIELDS];
  Array fluxX[N_FIELDS];
  Array fluxY[N_FIELDS];

  // set references to working fields
  Array& rho = fields[Rho];
  Array& ux = fields[Ux];
  Array& uy = fields[Uy];
  Array& bx = fields[Bx];
  Array& by = fields[By];
  Array& e = fields[E];
  
  // do not change this
  void allocateArrays()
  {
    for(int iField=0; iField<N_FIELDS; ++iField){
      fields[iField].allocate(-2, -2, nx+1, ny+1);
      fields_old[iField].allocate(-2, -2, nx+1, ny+1);
      fluxX[iField].allocate(0, 0, nx, ny-1);
      fluxY[iField].allocate(0, 0, nx-1, ny);
    }
  }
}

