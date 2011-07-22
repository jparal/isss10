// -*- C++ -*-
// $Id: timeStep.cc,v 1.4 2011/07/21 09:08:56 jd Exp $

#include "global.h"

// Do not change this file

namespace mhd2d
{

  // one increment in Runge-Kutta sequence. Fluxes are computed, now
  // f  <-  frac * f^n + (1-frac) * ( f - dt * div(H) )
  static void increment(real frac_old)
  {
    for(int iField=0; iField<N_FIELDS; ++iField){
      Array& f = fields[iField];
      Array& f_old = fields_old[iField];
      Array& Hx = fluxX[iField];
      Array& Hy = fluxY[iField];
      
      for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
          f(i, j) = frac_old*f_old(i, j) + (1.-frac_old)*
            ( f(i, j) - dt* ( dxi*(Hx(i+1, j) - Hx(i, j)) + dyi*(Hy(i, j+1) - Hy(i, j)) ) );
        }
      }
    }
  }

  // this is called from outside
  void doStep()
  {
    // low-storage Runge-Kutta scheme: save old field data
    for(int iField=0; iField<N_FIELDS; ++iField){
      Array& f = fields[iField];
      Array& f_old = fields_old[iField];
      for(int j=-2; j<ny+2; j++){
        for(int i=-2; i<nx+2; i++){
          f_old(i, j) = f(i, j);
        }
      }
    }

    // 3rd order Runge-Kutta time step (Shu & Osher, 1988)
    real frac_old[] = { 0., 3./4., 1./3. };

    // 3 substeps: compute flux, do increment according to R.-K. formula. Apply b.c.
    for(int subStep=0; subStep<3; ++subStep){
      calcFlux();
      increment(frac_old[subStep]);
      boundaryCond();
    }
  }
  
}


