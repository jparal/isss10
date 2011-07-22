/* 1d PIC library for pushing particles with darwin electric and magnetic */
/* fields and depositing current and derivative of current */

#include "dpush1lib_f.h"

void gmjpost1_(double *part, double *amu, double *qm, int *nop,
               int *idimp, int *nxv);

void gsmjpost1_(double *part, double *amu, double *qm, int *nop,
                int *idimp, int *nxv);

void gdcjpost1_(double *part, double *fxyz, double *byz, double *cu,
                double *dcu, double *amu, double *omx, double *qm,
                double *qbm, double *dt, int *idimp, int *nop,
                int *nxv);

void gsdcjpost1_(double *part, double *fxyz, double *byz, double *cu,
                 double *dcu, double *amu, double *omx, double *qm,
                 double *qbm, double *dt, int *idimp, int *nop,
                 int *nxv);

void gmjpost1l_(double *part, double *amu, double *qm, int *nop,
                int *idimp, int *nxv);

void gsmjpost1l_(double *part, double *amu, double *qm, int *nop,
                 int *idimp, int *nxv);

void gdcjpost1l_(double *part, double *fxyz, double *byz, double *cu,
                 double *dcu, double *amu, double *omx, double *qm,
                 double *qbm, double *dt, int *idimp, int *nop,
                 int *nxv);

void gsdcjpost1l_(double *part, double *fxyz, double *byz, double *cu,
                  double *dcu, double *amu, double *omx, double *qm,
                  double *qbm, double *dt, int *idimp, int *nop,
                  int *nxv);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void gmjpost1(double part[], double amu[], double qm, int nop,
              int idimp, int nxv) {
   gmjpost1_(part,amu,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gsmjpost1(double part[], double amu[], double qm, int nop,
               int idimp, int nxv) {
   gsmjpost1_(part,amu,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gdcjpost1(double part[], double fxyz[], double byz[], double cu[],
               double dcu[], double amu[], double omx, double qm,
               double qbm, double dt, int idimp, int nop, int nxv) {
   gdcjpost1_(part,fxyz,byz,cu,dcu,amu,&omx,&qm,&qbm,&dt,&idimp,&nop,
              &nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gsdcjpost1(double part[], double fxyz[], double byz[], double cu[],
                double dcu[], double amu[], double omx, double qm,
                double qbm, double dt, int idimp, int nop, int nxv) {
   gsdcjpost1_(part,fxyz,byz,cu,dcu,amu,&omx,&qm,&qbm,&dt,&idimp,&nop,
               &nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gmjpost1l(double part[], double amu[], double qm, int nop,
               int idimp, int nxv) {
   gmjpost1l_(part,amu,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gsmjpost1l(double part[], double amu[], double qm, int nop,
                int idimp, int nxv) {
   gsmjpost1l_(part,amu,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gdcjpost1l(double part[], double fxyz[], double byz[], double cu[],
                double dcu[], double amu[], double omx, double qm,
                double qbm, double dt, int idimp, int nop, int nxv) {
   gdcjpost1l_(part,fxyz,byz,cu,dcu,amu,&omx,&qm,&qbm,&dt,&idimp,&nop,
               &nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gsdcjpost1l(double part[], double fxyz[], double byz[],
                 double cu[], double dcu[], double amu[], double omx,
                 double qm, double qbm, double dt, int idimp, int nop,
                 int nxv) {
   gsdcjpost1l_(part,fxyz,byz,cu,dcu,amu,&omx,&qm,&qbm,&dt,&idimp,&nop,
                &nxv);
   return;
}
