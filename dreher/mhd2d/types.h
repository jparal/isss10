// -*- C++ -*-
// $Id: types.h,v 1.7 2011/07/21 09:08:57 jd Exp $

#ifndef TYPES_H
#define TYPES_H

// definition of arrays. do not change.

namespace mhd2d
{
  typedef double real;

  inline real sqr(real x) { return x*x; }
  
  struct Array
  {
    Array();

    void allocate(int i_min, int j_min, int i_max, int j_max);
    void deAllocate();

    real operator() (int i, int j) const { return zero_[index(i, j)]; }
    real& operator() (int i, int j) { return zero_[index(i, j)]; }
    ~Array();
    
  private:
    Array(const Array& a);
    real* data_;
    real* zero_;
    int s_;
    int index(int i, int j) const { return i + j*s_; }
  };

}

#endif
