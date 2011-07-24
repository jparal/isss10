/* header file for mpush1lib_f.c */

#include <complex.h>

void mgpost1(double part[], double q[], double qm, int nop, int idimp,
             int nxv, double qp[], int idtask[], int nmt, int *ierr);

void mgspost1(double part[], double q[], double qm, int nop, int idimp,
              int nxv, double qp[], int idtask[], int nmt, int *ierr);

void mgspost1x(double part[], double q[], double qm, int nop, int idimp,
               int nxv, double qp[], int idtask[], int nmt, int *ierr);

void mgpost1l(double part[], double q[], double qm, int nop, int idimp,
              int nxv, double qp[], int idtask[], int nmt, int *ierr);

void mgspost1l(double part[], double q[], double qm, int nop, int idimp,
               int nxv, double qp[], int idtask[], int nmt, int *ierr);

void mgspost1xl(double part[], double q[], double qm, int nop,
                int idimp, int nxv, double qp[], int idtask[], int nmt, 
                int *ierr);

void mgpush1(double part[], double fx[], double qbm, double dt,
             double *ek, int idimp, int nop, int nx, int nxv,
             double ekp[], int idtask[], int nmt, int *ierr);

void mgspush1(double part[], double fx[], double qbm, double dt,
              double *ek, int idimp, int nop, int nx, int nxv,
              double ekp[], int idtask[], int nmt, int *ierr);

void mgpush1l(double part[], double fx[], double qbm, double dt,
              double *ek, int idimp, int nop, int nx, int nxv,
              double ekp[], int idtask[], int nmt, int *ierr);

void mgspush1l(double part[], double fx[], double qbm, double dt,
               double *ek, int idimp, int nop, int nx, int nxv,
               double ekp[], int idtask[], int nmt, int *ierr);

void msortp1x(double part[], double pt[], int ip[], int npic[],
              int idimp, int nop, int nx1, int npicp[], int idtask[],
              int nmt, int *ierr);

void msortp1xl(double part[], double pt[], int ip[], int npic[],
               int idimp, int nop, int nx1, int npicp[], int idtask[],
               int nmt, int *ierr);

void mdsortp1x(double parta[], double partb[], int npic[], int idimp,
               int nop, int nx1, int npicp[], int idtask[], int nmt,
               int *ierr);

void mdsortp1xl(double parta[], double partb[], int npic[], int idimp,
                int nop, int nx1, int npicp[], int idtask[], int nmt,
                int *ierr);

void mdpost1gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int nxh,
               int nxvh, double complex qp[], double complex sctxp[],
               int idtask[], int nmt, int *ierr);

void mpush1gl(double part[], double complex fx[], double complex sctx[],
              double qbm, double dt, double *ek, int idimp, int nop,
              int nx, int nxh, int nxvh, double complex sctxp[],
              double ekp[], int idtask[], int nmt, int *ierr);
