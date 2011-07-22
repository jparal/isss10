// -*- C++ -*-
// $Id: reconstruct.cc,v 1.5 2011/07/21 14:49:38 jd Exp $

#include <cmath>
#include "global.h"

/*
  Functions to reconstruct variable values at cell interfaces
  from their cell avergages.
  Two methods are implemented: 
  a) a piecewise linear minMod-based reconstruction and
  b) a piecewise parabolic reconstruction with CWENO limiter
     (see e.g.  Kurganov & Levy, J. Sci. Comput., 22, p. 1461 (2000).

     To select either, (un-)comment the appropriate lines in function 
     recoMinus() below. Do not change the other functions.
*/

namespace mhd2d
{
  static real minMod(real d1, real d2)
  {
    // minMod-function
    return d1*d2>0. ? ( d1>0. ? fmin(d1, d2) : fmax(d1, d2) ) : 0.;
  }

  static real recoMinusLinear(real l, real c, real r)
  {
    // linear interpolation to lower cell interface
    return c - 0.5*minMod(r-c, c-l);
  }

  // these are smoothness indicators for the CWENO reconstruction. Do not change.
  const real EPS = 1e-6;
  static inline real alpha_l(real is_l) { return 0.25/ ( (EPS+is_l)*(EPS+is_l) ); }
  static inline real alpha_c(real is_c) { return 0.5 / ( (EPS+is_c)*(EPS+is_c) ); }

  static real recoMinusCWENO(real l, real c, real r)
  // CWENO reconstruction. See paper by Kurganov & Levy. Do not change.
  {
    real w_l = l-c;
    real w_r = c-r;
    
    w_l = alpha_l(w_l*w_l);
    w_r = alpha_l(w_r*w_r);
    
    real delta_c = 0.5*(r-l);
    real lap_c = l + r - 2.*c;
    real w_c = alpha_c(13./3.*lap_c*lap_c + delta_c*delta_c);

    real wNorm = 1. / (w_l + w_r + w_c);
    w_l *= wNorm;
    w_c *= wNorm;
    w_r *= wNorm;
     
    // calulate the parabola coefficients
    const real c1 = w_c*lap_c;
    const real a = c - c1/12.;
    const real b = ( (w_r + 0.5*w_c) * r + (-w_r + w_l) * c - (w_l + 0.5*w_c) * l );
  
    return a - 0.5*b + 0.25*c1;
  }

  // lower reconstruction.
  real recoMinus(real l, real c, real r)
  {
    /* Call either method here. (Un-)comment the appropriate lines. */
    //    return recoMinusLinear(l, c, r);
    return recoMinusCWENO(l, c, r);
  }

  // upper reconstruction = lower reconstruction mirrored. Do not change.
  real recoPlus(real l, real c, real r)
  {
    return recoMinus(r, c, l);
  }

}
