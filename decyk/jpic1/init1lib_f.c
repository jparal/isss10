/* 1d PIC library for initialization */
/* Wrappers for calling the Fortran routines from a C main program */

#include "init1lib_f.h"

void distr1_(double *part, double *vtx, double *vdx, int *npx,
             int *idimp, int *nop, int *nx, int *ipbc);

void distr1h_(double *part, double *vtx, double *vty, double *vtz,
              double *vdx, double *vdy, double *vdz, int *npx,
              int *idimp, int *nop, int *nx, int *ipbc);

double fnx_(double *x, double *argx1, double *argx2, double *argx3,
            int *intg);

void fdistr1_(double *part, double (*fnx_)(), double *argx1,
              double *argx2, double *argx3, int *npx, int *idimp,
              int *nop, int *nx, int *ipbc, int *ierr);

void vdistr1_(double *part, double *vtx, double *vdx, int *idimp,
              int *nop);

void vdistr1h_(double *part, double *vtx, double *vty, double *vtz,
               double *vdx, double *vdy, double *vdz, int *idimp,
               int *nop);

void gbdistr1l_(double *part, double *byz, double *qbm, int *idimp,
                int *nop, int *nx, int *nxv, int *ipbc);

double ranorm_();

double randum_();

double erfn_();

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void distr1(double part[], double vtx, double vdx, int npx, int idimp,
            int nop, int nx, int ipbc) {
   distr1_(part,&vtx,&vdx,&npx,&idimp,&nop,&nx,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void distr1h(double part[], double vtx, double vty, double vtz,
             double vdx, double vdy, double vdz, int npx, int idimp,
             int nop, int nx, int ipbc) {
   distr1h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&idimp,&nop,&nx,
            &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void fdistr1(double part[], double (*fnx_)(), double argx1,
             double argx2, double argx3, int npx, int idimp, int nop,
             int nx, int ipbc, int *ierr) {
   fdistr1_(part,fnx_,&argx1,&argx2,&argx3,&npx,&idimp,&nop,&nx,&ipbc,
            ierr);
   return;
}

/*--------------------------------------------------------------------*/
void vdistr1(double part[], double vtx, double vdx, int idimp, int nop) {
   vdistr1_(part,&vtx,&vdx,&idimp,&nop);
   return;
}

/*--------------------------------------------------------------------*/
void vdistr1h(double part[], double vtx, double vty, double vtz,
              double vdx, double vdy, double vdz, int idimp, int nop) {
   vdistr1h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&idimp,&nop);
   return;
}

/*--------------------------------------------------------------------*/
void gbdistr1l(double part[], double byz[], double qbm, int idimp,
               int nop, int nx, int nxv, int ipbc) {
   gbdistr1l_(part,byz,&qbm,&idimp,&nop,&nx,&nxv,&ipbc);
   return;
}

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

double erfn() {
   return erfn_();
}
