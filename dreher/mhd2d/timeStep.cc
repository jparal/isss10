// -*- C++ -*-
// $Id: timeStep.cc,v 1.3 2011/06/28 07:44:29 jd Exp $

#include "global.h"

namespace mhd2d
{

  static void increment(real frac_old)
  {
    for(int iField=0; iField<N_FIELDS; ++iField){
      Array& f = fields[iField];
      Array& f_old = fields_old[iField];
      Array& Fx = fluxX[iField];
      Array& Fy = fluxY[iField];
      
      for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
          f(i, j) = frac_old*f_old(i, j) + (1.-frac_old)*
            ( f(i, j) - dt* ( dxi*(Fx(i+1, j) - Fx(i, j)) + dyi*(Fy(i, j+1) - Fy(i, j)) ) );
        }
      }
    }
  }

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

    for(int subStep=0; subStep<3; ++subStep){
      calcFlux();
      increment(frac_old[subStep]);
      boundaryCond();
    }
  }
  
}


