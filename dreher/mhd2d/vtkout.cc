// -*- C++ -*-
// $Id: vtkout.cc,v 1.8 2011/07/21 09:08:57 jd Exp $

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "global.h"

// output to vtk file. do not change.

namespace mhd2d
{
  // big <-> little endian swap.
  static float swap(float f)
  {
    float res;
    char* in=(char*)(&f);
    char* out=(char*)(&res);
    out[0]=in[3];
    out[1]=in[2];
    out[2]=in[1];
    out[3]=in[0];
    return res;
  }

  static void writeArray(FILE* fp, const Array& f, const char* fieldName)
  {
    float line[nx];

    fprintf(fp, "SCALARS %s float\nLOOKUP_TABLE default\n", fieldName);
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        line[i] = swap(f(i, j));
      }
      fwrite(line, sizeof(float), nx, fp);
    }
  }
  
  void writeVTKFile(int step)
  {
    std::cout << "Output step " << step << " time=" << step*dt << std::endl;

    char filename[500];
    sprintf(filename, "%s_%04d.vtk", outFileStem, step);

    FILE* fp = fopen(filename, "wb");
    if(0==fp){
      std::cout << "Can't create output file " << filename << std::endl;
      exit(1);
    }

    fprintf(fp, "# vtk DataFile Version 2.0\nMHD2D\nBINARY\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n", nx, ny, 1);
    fprintf(fp, "ORIGIN %g %g %g\nSPACING %g %g %g\n", xmin, ymin, 0., dx, dy, 1.);
    fprintf(fp, "POINT_DATA %d\n", nx*ny);

    writeArray(fp, rho, "rho");
    writeArray(fp, ux, "ux");
    writeArray(fp, uy, "uy");
    writeArray(fp, bx, "bx");
    writeArray(fp, by, "by");
    writeArray(fp, e, "E");

    fclose(fp);
  }

}
