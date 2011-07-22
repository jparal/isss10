// -*- C++ -*-
// $Id: global.cc,v 1.10 2011/06/28 07:44:29 jd Exp $

#include <cmath>
#include "global.h"

namespace mhd2d
{

  const int nx = 150;
  const int ny = 120;
    
  const real xmin = -M_PI;
  const real ymin = -M_PI;
    
  const real xmax = M_PI;
  const real ymax = M_PI;
  
  const int nStep = 200;
  const int outputStep = 20;
  const char* outFileStem = "data/MHD";

  const real dx = (xmax-xmin) / nx;
  const real dy = (ymax-ymin) / ny;
  
  const real dt = .4*fmin(dx, dy);

  const real dxi = 1. / dx;
  const real dyi = 1. / dy;
    
  Array fields[N_FIELDS];

  Array fields_old[N_FIELDS];

  Array fluxX[N_FIELDS];
  Array fluxY[N_FIELDS];

  Array& rho = fields[Rho];
  
  Array& ux = fields[Ux];
  Array& uy = fields[Uy];
  
//   Array& bx = fields[Bx];
//   Array& by = fields[By];
  
  Array& w = fields[W];

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

