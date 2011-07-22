// -*- C++ -*-
// $Id: reconstruct.cc,v 1.2 2011/06/28 07:44:29 jd Exp $

#include <cmath>
#include "global.h"

namespace mhd2d
{
  static real minMod(real d1, real d2)
  {
    return d1*d2>0. ? ( d1>0. ? fmin(d1, d2) : fmax(d1, d2) ) : 0.;
  }

  static real recoMinusLinear(real l, real c, real r)
  {
    return c - 0.5*minMod(r-c, c-l);
  }

  const real EPS = 1e-6;
  static inline real alpha_l(real is_l) { return 0.25/ ( (EPS+is_l)*(EPS+is_l) ); }
  static inline real alpha_c(real is_c) { return 0.5 / ( (EPS+is_c)*(EPS+is_c) ); }

  static real recoMinusCWENO(real l, real c, real r)
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

  real recoMinus(real l, real c, real r)
  {
    //return recoMinusLinear(l, c, r);
    return recoMinusCWENO(l, c, r);
  }

  real recoPlus(real l, real c, real r)
  {
    return recoMinus(r, c, l);
  }

}
