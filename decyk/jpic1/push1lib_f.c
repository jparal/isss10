/* 1d C PIC library for pushing particles and depositing charge */
/* Wrappers for calling the Fortran routines from a C main program */

#include "push1lib_f.h"

void gpost1_(double *part, double *q, double *qm, int *nop, int *idimp,
             int *nxv);

void gspost1_(double *part, double *q, double *qm, int *nop,
              int *idimp, int *nxv);

void gspost1x_(double *part, double *q, double *qm, int *nop,
               int *idimp, int *nxv);

void gpost1l_(double *part, double *q, double *qm, int *nop,
              int *idimp, int *nxv);

void gspost1l_(double *part, double *q, double *qm, int *nop,
               int *idimp, int *nxv);

void gspost1xl_(double *part, double *q, double *qm, int *nop,
                int *idimp, int *nxv);

void gpush1_(double *part, double *fx, double *qbm, double *dt,
             double *ek, int *idimp, int *nop, int *nx, int *nxv,
             int *ipbc);

void gspush1_(double *part, double *fx, double *qbm, double *dt,
              double *ek, int *idimp, int *nop, int *nx, int *nxv,
              int *ipbc);

void gpush1l_(double *part, double *fx, double *qbm, double *dt,
              double *ek, int *idimp, int *nop, int *nx, int *nxv,
              int *ipbc);

void gspush1l_(double *part, double *fx, double *qbm, double *dt,
               double *ek, int *idimp, int *nop, int *nx, int *nxv,
               int *ipbc);

void sortp1x_(double *part, double *pt, int *ip, int *npic, int *idimp,
              int *nop, int *nx1);

void sortp1xl_(double *part, double *pt, int *ip, int *npic, int *idimp,
               int *nop, int *nx1);

void dsortp1x_(double *parta, double *partb, int *npic, int *idimp,
               int *nop, int *nx1);

void dsortp1xl_(double *parta, double *partb, int *npic, int *idimp,
                int *nop, int *nx1);

void rmove1_(double *part, int *ihole, int *nx, int *idimp, int *nop,
             int *ntmax, int *ipbc);

void dpost1gl_(double *part, double complex *q, double complex *sctx,
               double *qm, int *nop, int *idimp, int *nx, int *nxh,
               int *nxvh);

void push1gl_(double *part, double complex *fx, double complex *sctx,
              double *qbm, double *dt, double *ek, int *idimp, int *nop,
              int *nx, int *nxh, int *nxvh);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void gpost1(double part[], double q[], double qm, int nop, int idimp,
            int nxv) {
   gpost1_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gspost1(double part[], double q[], double qm, int nop, int idimp,
             int nxv) {
   gspost1_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gspost1x(double part[], double q[], double qm, int nop, int idimp,
              int nxv) {
   gspost1x_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gpost1l(double part[], double q[], double qm, int nop, int idimp,
             int nxv) {
   gpost1l_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gspost1l(double part[], double q[], double qm, int nop, int idimp,
              int nxv) {
   gspost1l_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gspost1xl(double part[], double q[], double qm, int nop, int idimp,
               int nxv) {
   gspost1xl_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void gpush1(double part[], double fx[], double qbm, double dt,
            double *ek, int idimp, int nop, int nx, int nxv, int ipbc) {
   gpush1_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gspush1(double part[], double fx[], double qbm, double dt,
             double *ek, int idimp, int nop, int nx, int nxv, int ipbc) {
   gspush1_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gpush1l(double part[], double fx[], double qbm, double dt,
             double *ek, int idimp, int nop, int nx, int nxv, int ipbc) {
   gpush1l_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void gspush1l(double part[], double fx[], double qbm, double dt,
              double *ek, int idimp, int nop, int nx, int nxv, int ipbc) {
   gspush1l_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void sortp1x(double part[], double pt[], int ip[], int npic[],
             int idimp, int nop, int nx1) {
   sortp1x_(part,pt,ip,npic,&idimp,&nop,&nx1);
   return;
}

/*--------------------------------------------------------------------*/
void sortp1xl(double part[], double pt[], int ip[], int npic[],
              int idimp, int nop, int nx1) {
   sortp1xl_(part,pt,ip,npic,&idimp,&nop,&nx1);
   return;
}

/*--------------------------------------------------------------------*/
void dsortp1x(double parta[], double partb[], int npic[], int idimp,
              int nop, int nx1) {
   dsortp1x_(parta,partb,npic,&idimp,&nop,&nx1);
   return;
}

/*--------------------------------------------------------------------*/
void dsortp1xl(double parta[], double partb[], int npic[], int idimp,
                int nop, int nx1) {
   dsortp1xl_(parta,partb,npic,&idimp,&nop,&nx1);
   return;
}

/*--------------------------------------------------------------------*/
void rmove1(double part[], int ihole[], int nx, int idimp, int nop,
            int ntmax, int ipbc) {
   rmove1_(part,ihole,&nx,&idimp,&nop,&ntmax,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void dpost1gl(double part[], double complex q[], double complex sctx[],
              double qm, int nop, int idimp, int nx, int nxh,
              int nxvh) {
   dpost1gl_(part,q,sctx,&qm,&nop,&idimp,&nx,&nxh,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void push1gl(double part[], double complex fx[], double complex sctx[],
              double qbm, double dt, double *ek, int idimp, int nop,
              int nx, int nxh, int nxvh) {
   push1gl_(part,fx,sctx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxh,&nxvh);
   return;
}

