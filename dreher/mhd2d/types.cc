// -*- C++ -*-
// $Id: types.cc,v 1.5 2011/06/28 07:44:29 jd Exp $

#include <cassert>
#include "types.h"

namespace mhd2d
{

  Array::Array() 
    : data_(0), zero_(0), s_(0)
  {
  }

  void Array::allocate(int i_min, int j_min, int i_max, int j_max)
  {
    deAllocate();

    assert(i_max-i_min+1 > 0);
    assert(j_max-j_min+1 > 0);

    s_ = i_max-i_min+1;
    data_ = new real[ s_ * (j_max-j_min+1) ];
    assert(0!=data_);
    zero_ = data_ - index(i_min, j_min);
  }
  
  Array::~Array()
  {
    deAllocate();
  }

  void Array::deAllocate()
  {
    if(0!=data_){
      delete[] data_;
    }
    data_ = zero_ = 0;
    s_ = 0;
  }
  
}
