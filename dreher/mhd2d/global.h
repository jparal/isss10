// -*- C++ -*-
// $Id: global.h,v 1.8 2011/06/28 07:44:29 jd Exp $

#ifndef GLOBAL_H
#define GLOBAL_H

#include "types.h"

namespace mhd2d
{

  extern const int nx;
  extern const int ny;
  
  extern const real xmin;
  extern const real ymin;

  extern const real xmax;
  extern const real ymax;

  extern const real dx;
  extern const real dy;

  extern const real dxi;
  extern const real dyi;

  extern const real dt;

  extern const char* outFileStem;

  enum { Rho=0, Ux, Uy, W, N_FIELDS};

  extern Array fields[N_FIELDS];

  extern Array fields_old[N_FIELDS];

  extern Array fluxX[N_FIELDS];
  extern Array fluxY[N_FIELDS];

  extern Array& rho;

  extern Array& ux;
  extern Array& uy;

  extern Array& bx;
  extern Array& by;

  extern Array& w;


  extern void allocateArrays();
  extern void initFields();
  extern void doStep();
  extern void calcFlux();
  extern void boundaryCond();
  extern void writeVTKFile(int step);

  extern real recoPlus(real l, real c, real r);
  extern real recoMinus(real l, real c, real r);

  extern const int nStep;
  extern const int outputStep;
}

#endif
