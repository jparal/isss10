/* header file for push1lib_f.c */

#include <complex.h>

void gpost1(double part[], double q[], double qm, int nop, int idimp,
            int nxv);

void gspost1(double part[], double q[], double qm, int nop, int idimp,
             int nxv);

void gspost1x(double part[], double q[], double qm, int nop, int idimp,
              int nxv);

void gpost1l(double part[], double q[], double qm, int nop, int idimp,
             int nxv);

void gspost1l(double part[], double q[], double qm, int nop, int idimp,
              int nxv);

void gspost1xl(double part[], double q[], double qm, int nop, int idimp,
               int nxv);

void gpush1(double part[], double fx[], double qbm, double dt,
            double *ek, int idimp, int nop, int nx, int nxv, int ipbc);

void gspush1(double part[], double fx[], double qbm, double dt,
             double *ek, int idimp, int nop, int nx, int nxv, int ipbc);

void gpush1l(double part[], double fx[], double qbm, double dt,
             double *ek, int idimp, int nop, int nx, int nxv, int ipbc);

void gspush1l(double part[], double fx[], double qbm, double dt,
              double *ek, int idimp, int nop, int nx, int nxv, int ipbc);

void sortp1x(double part[], double pt[], int ip[], int npic[],
             int idimp, int nop, int nx1) ;

void sortp1xl(double part[], double pt[], int ip[], int npic[],
              int idimp, int nop, int nx1);

void dsortp1x(double parta[], double partb[], int npic[], int idimp,
              int nop, int nx1);

void dsortp1xl(double parta[], double partb[], int npic[], int idimp,
                int nop, int nx1);

void rmove1(double part[], int ihole[], int nx, int idimp, int nop,
            int ntmax, int ipbc);

void dpost1gl(double part[], double complex q[], double complex sctx[],
              double qm, int nop, int idimp, int nx, int nxh,
              int nxvh);

void push1gl(double part[], double complex fx[], double complex sctx[],
              double qbm, double dt, double *ek, int idimp, int nop,
              int nx, int nxh, int nxvh);

