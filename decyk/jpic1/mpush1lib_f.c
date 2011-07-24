/* 1d PIC multi-tasking library for pushing particles and depositing */
/* charge                                                            */
/* Wrappers for calling the Fortran routines from a C main program */

#include "mpush1lib_f.h"

void mgpost1_(double *part, double *q, double *qm, int *nop, int *idimp,
              int *nxv, double *qp, int *idtask, int *nmt, int *ierr);

void mgspost1_(double *part, double *q, double *qm, int *nop,
               int *idimp, int *nxv, double *qp, int *idtask, int *nmt,
               int *ierr);

void mgspost1x_(double *part, double *q, double *qm, int *nop,
                int *idimp, int *nxv, double *qp, int *idtask, int *nmt,
                int *ierr);

void mgpost1l_(double *part, double *q, double *qm, int *nop,
               int *idimp, int *nxv, double *qp, int *idtask, int *nmt,
               int *ierr);

void mgspost1l_(double *part, double *q, double *qm, int *nop,
                int *idimp, int *nxv, double *qp, int *idtask, int *nmt,
                int *ierr);

void mgspost1xl_(double *part, double *q, double *qm, int *nop,
                 int *idimp, int *nxv, double *qp, int *idtask,
                 int *nmt, int *ierr);

void mgpush1_(double *part, double *fx, double *qbm, double *dt,
              double *ek, int *idimp, int *nop, int *nx, int *nxv,
              double *ekp, int *idtask, int *nmt, int *ierr);

void mgspush1_(double *part, double *fx, double *qbm, double *dt,
               double *ek, int *idimp, int *nop, int *nx, int *nxv,
               double *ekp, int *idtask, int *nmt, int *ierr);

void mgpush1l_(double *part, double *fx, double *qbm, double *dt,
               double *ek, int *idimp, int *nop, int *nx, int *nxv,
               double *ekp, int *idtask, int *nmt, int *ierr);

void mgspush1l_(double *part, double *fx, double *qbm, double *dt,
                double *ek, int *idimp, int *nop, int *nx, int *nxv,
                double *ekp, int *idtask, int *nmt, int *ierr);

void msortp1x_(double *part, double *pt, int *ip, int *npic, int *idimp,
               int *nop, int *nx1, int *npicp, int *idtask, int *nmt,
               int *ierr);

void msortp1xl_(double *part, double *pt, int *ip, int *npic,
                int *idimp, int *nop, int *nx1, int *npicp, int *idtask,
                int *nmt, int *ierr);

void mdsortp1x_(double *parta, double *partb, int *npic, int *idimp,
                int *nop, int *nx1, int *npicp, int *idtask, int *nmt,
                int *ierr);

void mdsortp1xl_(double *parta, double *partb, int *npic, int *idimp,
                 int *nop, int *nx1, int *npicp, int *idtask, int *nmt,
                 int *ierr);

void mdpost1gl_(double *part, double complex *q, double complex *sctx,
                double *qm, int *nop, int *idimp, int *nx, int *nxh,
                int *nxvh, double complex *qp, double complex *sctxp,
                int *idtask, int *nmt, int *ierr);

void mpush1gl_(double *part, double complex *fx, double complex *sctx,
               double *qbm, double *dt, double *ek, int *idimp,
               int *nop, int *nx, int *nxh, int *nxvh,
               double complex *sctxp, double *ekp, int *idtask,
               int *nmt, int *ierr);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void mgpost1(double part[], double q[], double qm, int nop, int idimp,
             int nxv, double qp[], int idtask[], int nmt, int *ierr) {
   mgpost1_(part,q,&qm,&nop,&idimp,&nxv,qp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgspost1(double part[], double q[], double qm, int nop, int idimp,
              int nxv, double qp[], int idtask[], int nmt, int *ierr) {
   mgspost1_(part,q,&qm,&nop,&idimp,&nxv,qp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgspost1x(double part[], double q[], double qm, int nop, int idimp,
               int nxv, double qp[], int idtask[], int nmt, int *ierr) {
   mgspost1x_(part,q,&qm,&nop,&idimp,&nxv,qp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgpost1l(double part[], double q[], double qm, int nop, int idimp,
              int nxv, double qp[], int idtask[], int nmt, int *ierr) {
   mgspost1x_(part,q,&qm,&nop,&idimp,&nxv,qp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgspost1l(double part[], double q[], double qm, int nop, int idimp,
               int nxv, double qp[], int idtask[], int nmt, int *ierr) {
   mgspost1l_(part,q,&qm,&nop,&idimp,&nxv,qp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgspost1xl(double part[], double q[], double qm, int nop,
                int idimp, int nxv, double qp[], int idtask[], int nmt, 
                int *ierr) {
   mgspost1xl_(part,q,&qm,&nop,&idimp,&nxv,qp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgpush1(double part[], double fx[], double qbm, double dt,
             double *ek, int idimp, int nop, int nx, int nxv,
             double ekp[], int idtask[], int nmt, int *ierr) {
   mgpush1_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,ekp, idtask,&nmt,
            ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgspush1(double part[], double fx[], double qbm, double dt,
              double *ek, int idimp, int nop, int nx, int nxv,
              double ekp[], int idtask[], int nmt, int *ierr) {
   mgspush1_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,ekp,idtask,&nmt,
             ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgpush1l(double part[], double fx[], double qbm, double dt,
              double *ek, int idimp, int nop, int nx, int nxv,
              double ekp[], int idtask[], int nmt, int *ierr) {
   mgspush1_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,ekp,idtask,&nmt,
             ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mgspush1l(double part[], double fx[], double qbm, double dt,
               double *ek, int idimp, int nop, int nx, int nxv,
               double ekp[], int idtask[], int nmt, int *ierr) {
   mgspush1l_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,ekp,idtask,&nmt,
              ierr);
   return;
}

/*--------------------------------------------------------------------*/
void msortp1x(double part[], double pt[], int ip[], int npic[],
              int idimp, int nop, int nx1, int npicp[], int idtask[],
              int nmt, int *ierr) {
   msortp1x_(part,pt,ip,npic,&idimp,&nop,&nx1,npicp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void msortp1xl(double part[], double pt[], int ip[], int npic[],
               int idimp, int nop, int nx1, int npicp[], int idtask[],
               int nmt, int *ierr) {
   msortp1xl_(part,pt,ip,npic,&idimp,&nop,&nx1,npicp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mdsortp1x(double parta[], double partb[], int npic[], int idimp,
               int nop, int nx1, int npicp[], int idtask[], int nmt,
               int *ierr) {
   mdsortp1x_(parta,partb,npic,&idimp,&nop,&nx1,npicp,idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mdsortp1xl(double parta[], double partb[], int npic[], int idimp,
                int nop, int nx1, int npicp[], int idtask[], int nmt,
                int *ierr) {
   mdsortp1xl_(parta,partb,npic,&idimp,&nop,&nx1,npicp,idtask,&nmt,
               ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mdpost1gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int nxh,
               int nxvh, double complex qp[], double complex sctxp[],
               int idtask[], int nmt, int *ierr) {
   mdpost1gl_(part,q,sctx,&qm,&nop,&idimp,&nx,&nxh,&nxvh,qp,sctxp,
              idtask,&nmt,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void mpush1gl(double part[], double complex fx[], double complex sctx[],
              double qbm, double dt, double *ek, int idimp, int nop,
              int nx, int nxh, int nxvh, double complex sctxp[],
              double ekp[], int idtask[], int nmt, int *ierr) {
   mpush1gl_(part,fx,sctx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxh,&nxvh,sctxp,
             ekp,idtask,&nmt,ierr);
   return;
}
