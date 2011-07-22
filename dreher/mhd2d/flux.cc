// -*- C++ -*-
// $Id: flux.cc,v 1.7 2011/07/21 09:08:56 jd Exp $

#include <cmath>
#include "global.h"

/* 
   Computation of the numerical fluxes H(u) from conserved variables at cell interfaces.

   speedX/Y() estimate the wave speed from reconstructed field. They enter the diffusive
   flux. Change if you change the equations to solve.

   fluxFunctionX/Y() compute the physical fluxes from reconstructed fields at cell interfaces.
   They define the equations that are integrated (e.g. MHD). Change this to integrate other
   conservation laws than MHD.

   calcFlux() at the bottom does the entire flux computation and calls the other functions.
   Do not change.

*/
namespace mhd2d
{
  static real speedX(real* fields)
  // fields contain the reconstructed values at cell interface. Compute max. phase speed.
  {
    //
    real rhoInv = 1./fields[Rho];
    real p_mag = 0.5*(sqr(fields[Bx]) + sqr(fields[By]));
    real p = (gamma-1.)*(fields[E] - 0.5*rhoInv*(sqr(fields[Ux])+sqr(fields[Uy])) - p_mag);
    return sqrt( (gamma*p + 2.*p_mag)*rhoInv ) + rhoInv*fabs(fields[Ux]); // MHD fast speed + convection |ux|/rho
  }

  static real speedY(real* fields)
  {
    real rhoInv = 1./fields[Rho];
    real p_mag = 0.5*(sqr(fields[Bx]) + sqr(fields[By]));
    real p = (gamma-1.)*(fields[E] - 0.5*rhoInv*(sqr(fields[Ux])+sqr(fields[Uy])) - p_mag);
    return sqrt( (gamma*p + 2.*p_mag)*rhoInv ) + rhoInv*fabs(fields[Uy]); // MHD fast speed + convection |uy|/rho
  }

  // The MHD flux function, x-components:
  static void fluxFunctionX(real* fluxes, real* fields)
  {
    real rhoInv = 1./fields[Rho];
    real ux = fields[Ux];
    real uy = fields[Uy];
    real bx = fields[Bx];
    real by = fields[By];
    real e = fields[E];

    real p_mag = 0.5*(sqr(bx) + sqr(by));
    real p = (gamma-1.)*(e - 0.5*rhoInv*(sqr(ux)+sqr(uy)) - p_mag);

    fluxes[Rho] = ux;
    fluxes[Ux] = ux*ux*rhoInv + p + p_mag - bx*bx;
    fluxes[Uy] = ux*uy*rhoInv - bx*by;
    fluxes[Bx] = 0.;
    fluxes[By] = rhoInv*(ux*by - bx*uy);
    fluxes[E] = rhoInv*( ux*(e + p + p_mag) - bx*(ux*bx + uy*by) );
  }

  // The MHD flux function, y-components:
  static void fluxFunctionY(real* fluxes, real* fields)
  {
    real rhoInv = 1./fields[Rho];
    real ux = fields[Ux];
    real uy = fields[Uy];
    real bx = fields[Bx];
    real by = fields[By];
    real e = fields[E];

    real p_mag = 0.5*(sqr(bx)+sqr(by));
    real p = (gamma-1.)*(e - 0.5*rhoInv*(sqr(ux)+sqr(uy)) - p_mag);

    fluxes[Rho] = uy;
    fluxes[Ux] = uy*ux*rhoInv - by*bx;
    fluxes[Uy] = uy*uy*rhoInv + p + p_mag - by*by;
    fluxes[Bx] = rhoInv*(uy*bx - by*ux);
    fluxes[By] = 0.;
    fluxes[E] = rhoInv*( uy*(e + p + p_mag) - by*(ux*bx + uy*by) );
  }

  // Entire flux computation is done here. Result is written in global arrays fluxX and fluxY.
  // It is then used in the Runge-Kutta increment (file timeStep.cc).
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

        // for each conserved variable...
        for(int iField=0; iField<N_FIELDS; ++iField){
          Array& f=fields[iField];
          // compute left and right reconstruction at left interface
          leftReco[iField] = recoPlus(f(i-2, j), f(i-1, j), f(i, j));
          rightReco[iField] = recoMinus(f(i-1, j), f(i, j), f(i+1, j));
        }

        // compute physical fluxes from left- and right-reconstruced values
        fluxFunctionX(fluxL, leftReco);
        fluxFunctionX(fluxR, rightReco);

        // estimate the maximum speed from the reconstructed values
        speedL=speedX(leftReco);
        speedR=speedX(rightReco);
        real maxSpeed=fmax(speedL, speedR);
        
        // mix the physical fluxes and the diffusive fluxes into arrays fluxX
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

        speedL=speedY(leftReco);
        speedR=speedY(rightReco);
        real maxSpeed=fmax(speedL, speedR);
        
        for(int iField=0; iField<N_FIELDS; ++iField){
          fluxY[iField](i, j) = 0.5*( (fluxL[iField] + fluxR[iField]) - maxSpeed*(rightReco[iField]-leftReco[iField]) );
        }        
      }
    }
  }
  
}


