// -*- C++ -*-
// $Id: global.h,v 1.10 2011/07/21 09:08:56 jd Exp $

#ifndef GLOBAL_H
#define GLOBAL_H

#include "types.h"

// Global declarations of variables and functions (subroutines)
namespace mhd2d
{
  // No. of grid cells
  extern const int nx;
  extern const int ny;
  
  // domain ranges. real=double
  extern const real xmin;
  extern const real ymin;
  extern const real xmax;
  extern const real ymax;

  // grid spacings and inverse grid spacings
  extern const real dx;
  extern const real dy;
  extern const real dxi;
  extern const real dyi;

  // time step size
  extern const real dt;

  // no. of integration steps for the simulation
  extern const int nStep;
  // data output every ... step, directory & name of output files
  extern const int outputStep;
  extern const char* outFileStem;

  // definition of the dependent variables (MHD fields)
  enum { Rho=0, Ux, Uy, Bx, By, E, N_FIELDS};
  // adiabtic index
  extern const real gamma;

  // references to individual arrays for easy access:
  extern Array& rho;
  extern Array& ux;
  extern Array& uy;
  extern Array& bx;
  extern Array& by;
  extern Array& e;

  // -- do not change lines below --

  // working fields, helper fields for Runge-Kutta (old step), 
  // fluxes for finite-volume update
  extern Array fields[N_FIELDS];
  extern Array fields_old[N_FIELDS];
  extern Array fluxX[N_FIELDS];
  extern Array fluxY[N_FIELDS];

  // global functions
  extern void allocateArrays();
  extern void initFields();
  extern void doStep();
  extern void calcFlux();
  extern void boundaryCond();
  extern void writeVTKFile(int step);

  // reconstruction in central scheme
  extern real recoPlus(real l, real c, real r);
  extern real recoMinus(real l, real c, real r);
}

#endif
