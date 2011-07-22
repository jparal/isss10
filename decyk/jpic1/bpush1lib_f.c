/* 1d PIC library for pushing particles with magnetic field */
/* and depositing current */

#include "bpush1lib_f.h"

void gjpost1_(double *part, double *cu, double *qm, double *dt,
              int *nop, int *idimp, int *nx, int *nxv, int *ipbc);

void gsjpost1_(double *part, double *cu, double *qm, double *dt,
               int *nop, int *idimp, int *nx, int *nxv, int *ipbc);

void gsjpost1x_(double *part, double *cu, double *qm, double *dt,
                int *nop, int *idimp, int *nx, int *nxv, int *ipbc);

void gjpost1l_(double *part, double *cu, double *qm, double *dt,
               int *nop, int *idimp, int *nx, int *nxv, int *ipbc);

void gsjpost1l_(double *part, double *cu, double *qm, double *dt,
                int *nop, int *idimp, int *nx, int *nxv, int *ipbc);

void gsjpost1xl_(double *part, double *cu, double *qm, double *dt,
                 int *nop, int *idimp, int *nx, int *nxv, int *ipbc);

void gbpush13_(double *part, double *fxyz, double *byz, double *omx,
               double *qbm, double *dt, double *dtc, double *ek,
               int *idimp, int *nop, int *nx, int *nxv, int *ipbc);

void gsbpush13_(double *part, double *fxyz, double *byz, double *omx,
                double *qbm, double *dt, double *dtc, double *ek,
                int *idimp, int *nop, int *nx, int *nxv, int *ipbc);

void gbpush13l_(double *part, double *fxyz, double *byz, double *omx,
                double *qbm, double *dt, double *dtc, double *ek,
                int *idimp, int *nop, int *nx, int *nxv, int *ipbc);

void gsbpush13l_(double *part, double *fxyz, double *byz, double *omx,
                 double *qbm, double *dt, double *dtc, double *ek,
                 int *idimp, int *nop, int *nx, int *nxv, int *ipbc);

void retard1_(double *part, double *dtc, int *idimp, int *nop, int *nx,
              int *ipbc);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void gjpost1(double part[], double cu[], double qm, double dt,
             int nop, int idimp, int nx, int nxv, int ipbc) {
   gjpost1_(part,cu,&qm,&dt,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gsjpost1(double part[], double cu[], double qm, double dt,
              int nop, int idimp, int nx, int nxv, int ipbc) {
   gsjpost1_(part,cu,&qm,&dt,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gsjpost1x(double part[], double cu[], double qm, double dt,
               int nop, int idimp, int nx, int nxv, int ipbc) {
   gsjpost1x_(part,cu,&qm,&dt,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gjpost1l(double part[], double cu[], double qm, double dt,
              int nop, int idimp, int nx, int nxv, int ipbc) {
   gjpost1l_(part,cu,&qm,&dt,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gsjpost1l(double part[], double cu[], double qm, double dt,
               int nop, int idimp, int nx, int nxv, int ipbc) {
   gsjpost1l_(part,cu,&qm,&dt,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gsjpost1xl(double part[], double cu[], double qm, double dt,
                int nop, int idimp, int nx, int nxv, int ipbc) {
   gsjpost1xl_(part,cu,&qm,&dt,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gbpush13(double part[], double fxyz[], double byz[], double omx,
              double qbm, double dt, double dtc, double *ek, int idimp,
              int nop, int nx, int nxv, int ipbc) {
   gbpush13_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&nxv,
             &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gsbpush13(double part[], double fxyz[], double byz[], double omx,
               double qbm, double dt, double dtc, double *ek, int idimp,
               int nop, int nx, int nxv, int ipbc) {
   gsbpush13_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&nxv,
              &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gbpush13l(double part[], double fxyz[], double byz[], double omx,
               double qbm, double dt, double dtc, double *ek, int idimp,
               int nop, int nx, int nxv, int ipbc) {
   gbpush13l_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&nxv,
              &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gsbpush13l(double part[], double fxyz[], double byz[], double omx,
                double qbm, double dt, double dtc, double *ek,
                int idimp, int nop, int nx, int nxv, int ipbc) {
   gsbpush13l_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&nxv,
               &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void retard1(double part[], double dtc, int idimp, int nop, int nx,
             int ipbc) {
   retard1_(part,&dtc,&idimp,&nop,&nx,&ipbc);
   return;
}
