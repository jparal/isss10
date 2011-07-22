// -*- C++ -*-
// $Id: flux.cc,v 1.5 2011/06/28 07:44:29 jd Exp $

#include <cmath>
#include "global.h"

namespace mhd2d
{
  static const real GAMMA=5./3.;

  static real speed(real* fields)
  {
    return sqrt( fields[Ux]*fields[Ux] + fields[Uy]*fields[Uy]);
  }

  static void fluxFunctionX(real* fluxes, real* fields)
  {
    real rho = fields[Rho];
    real ux = fields[Ux];
    real uy = fields[Uy];
    real w = fields[W];
    real p = (GAMMA-1.)*(w - 0.5*(ux*ux+uy*uy));

    fluxes[Rho] = ux;
    fluxes[Ux] = ux*ux/rho + p;
    fluxes[Uy] = ux*uy/rho;
    fluxes[W] = ux*(w+p);
  }

  static void fluxFunctionY(real* fluxes, real* fields)
  {
    real rho = fields[Rho];
    real ux = fields[Ux];
    real uy = fields[Uy];
    real w = fields[W];
    real p = (GAMMA-1.)*(w - 0.5*(ux*ux+uy*uy));

    fluxes[Rho] = uy;
    fluxes[Ux] = uy*ux/rho;
    fluxes[Uy] = uy*uy/rho + p;
    fluxes[W] = uy*(w+p);
  }

  void calcFlux()
  {
    real leftReco[N_FIELDS];
    real rightReco[N_FIELDS];

    real fluxL[N_FIELDS];
    real fluxR[N_FIELDS];

    real speedL;
    real speedR;

    // x-component of fluxes: loop over the cell boundaries, compute flux at *left* cell boundary
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx+1; i++){

        for(int iField=0; iField<N_FIELDS; ++iField){
          Array& f=fields[iField];
          leftReco[iField] = recoPlus(f(i-2, j), f(i-1, j), f(i, j));
          rightReco[iField] = recoMinus(f(i-1, j), f(i, j), f(i+1, j));
        }

        fluxFunctionX(fluxL, leftReco);
        fluxFunctionX(fluxR, rightReco);

        speedL=speed(leftReco);
        speedR=speed(rightReco);
        real maxSpeed=fmax(speedL, speedR);
        
        for(int iField=0; iField<N_FIELDS; ++iField){
          fluxX[iField](i, j) = 0.5*( (fluxL[iField] + fluxR[iField]) - maxSpeed*(rightReco[iField]-leftReco[iField]) );
        }        
      }
    }
  
    // same for y-component of fluxes:
    for(int j=0; j<ny+1; j++){
      for(int i=0; i<nx; i++){

        for(int iField=0; iField<N_FIELDS; ++iField){
          Array& f=fields[iField];
          leftReco[iField] = recoPlus(f(i, j-2), f(i, j-1), f(i, j));
          rightReco[iField] = recoMinus(f(i, j-1), f(i, j), f(i, j+1));
        }

        fluxFunctionY(fluxL, leftReco);
        fluxFunctionY(fluxR, rightReco);

        speedL=speed(leftReco);
        speedR=speed(rightReco);
        real maxSpeed=fmax(speedL, speedR);
        
        for(int iField=0; iField<N_FIELDS; ++iField){
          fluxY[iField](i, j) = 0.5*( (fluxL[iField] + fluxR[iField]) - maxSpeed*(rightReco[iField]-leftReco[iField]) );
        }        
      }
    }
  }
  
}


