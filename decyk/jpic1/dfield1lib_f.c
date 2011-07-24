/* 1d PIC library for solving field equations with dirichlet boundary */
/* conditions                                                         */
/* Wrappers for calling the Fortran routines from a C main program */

#include "dfield1lib_f.h"

void lcguard1_(double *byz, int *nx, int *nxe);

void lbguard1_(double *fxyz, int *nx, int *nxe);

void ldguard1_(double *fx, int *nx, int *nxe);

void lscguard1_(double *cu, double *yj0, double *zj0, int *nx,
                int *ngx, int *nxe);

void lsguard1_(double *q, double *qi0, int *nx, int *ngx, int *nxe);

void lacguard1_(double *cu, int *nx, int *nxe);

void laguard1_(double *q, int *nx, int *nxe);

void lscguard1l_(double *cu, double *yj0, double *zj0, int *nx,
                 int *ngx, int *nxe);

void lsguard1l_(double *q, double *qi0, int *nx, int *ngx, int *nxe);

void dblsin1a_(double *cu, double *cu2, int *nx, int *nxv, int *nx2v);

void dblsin1d_(double *q, double *q2, int *nx, int *nxv, int *nx2v);

void dblcos1a_(double *cu, double *cu2, int *nx, int *nxv, int *nx2v);

void dblcos1d_(double *q, double *q2, int *nx, int *nxv, int *nx2v);

void hafdbl1c_(double *byz, double *byz2, int *nx, int *nxe, int *nx2v);

void hafdbl1b_(double *fxyz, double *fxyz2, int *nx, int *nxe,
               int *nx2v);

void hafdbl1d_(double *q, double *q2, int *nx, int *nxe, int *nx2v);

void poisdx1_(double *q, double *fx, int *isign, double *ffd,
              double *ax, double *affp, double *we, int *nx);

void poisd1_(double *q, double *fx, int *isign, double *ffd, double *ax,
             double *affp, double *we, int *nx, int *nxe, int *nx2v);

void bpoisdx13_(double complex *cu, double complex *byz, int *isign,
                double complex *ffd, double *ax, double *affp,
                double *ci, double *wm, int *nx, int *nxv, int *nxd);

void bpoisd13_(double *cu, double *byz, int *isign, double complex *ffd,
               double *ax, double *affp, double *ci, double *wm,
               int *nx, int *nxe, int *nxv);

void ibpoisdx13_(double complex *cu, double complex *byz,
                 double complex *ffd, double *ci, double *wm, int *nx,
                 int *nxv, int *nxd);

void ibpoisd13_(double *cu, double *byz, double complex *ffd,
                double *ci, double *wm, int *nx, int *nxe, int *nxv);

void maxweldx1_(double complex *eyz, double complex *byz,
                double complex *cu, double complex *ffd, double *ci,
                double *dt, double *wf, double *wm, int *nx, int *nxv,
                int *nxd);

void maxweld1_(double *eyz, double *byz, double *cu,
               double complex *ffd, double *ci, double *dt, double *wf,
               double *wm, int *nx, int *nxe, int *nxv);

void dmfieldd1_(double complex *q2, double *q, int *nx, int *nxv,
                int *nxe);

void cmfieldd1_(double complex *cu2, double *cu, int *nx, int *nxv,
                int *nxe);

void amfieldd1_(double complex *amu2, double *amu, int *nx, int *nxv,
                int *nxe);

void emfieldd1_(double complex *fxyz, double complex *fx, double *eyz,
                double complex *ffd, int *nx, int *nxv, int *nxe,
                int *nxd);

void bmfieldd1_(double complex *fxyz, double *eyz, double complex *ffd,
                int *nx, int *nxv, int *nxe, int *nxd);

void avpotdx13_(double complex *byz, double complex *ayz, int *nx,
                int *nxv);

void avpotd13_(double *byz, double *ayz, int *nx, int *nxe);

void avrpotdx13_(double complex *ayz, double complex *byz,
                 double complex *ffd, double *ci, int *nx, int *nxv,
                 int *nxd);

void avrpotd13_(double *ayz, double *byz, double complex *ffd,
                double *ci, int *nx, int *nxe, int *nxd);

void gtsmodes1_(double *pot, double *pott, int *nx, int *it,
                int *modesx, int *nxe, int *nt, int *modesxd);

void ptsmodes1_(double *pot, double *pott, int *nx, int *it,
                int *modesx, int *nxe, int *nt, int *modesxd);

void gtvsmodes1_(double *vpot, double *vpott, int *nx, int *it,
                 int *modesx, int *ndim, int *nxv, int *nt,
                 int *modesxd);

void ptvsmodes1_(double *vpot, double *vpott, int *nx, int *it,
                 int *modesx, int *ndim, int *nxv, int *nt,
                 int *modesxd);

void lscfguard1_(double *cus, double *cu, double *q2m0, int *nx,
                 int *nxe);

void lscfguard1l_(double *cus, double *cu, double *q2m0, int *nx,
                  int *nxe);

void dcuperpdx13_(double complex *dcu, double complex *amu, int *nx,
                  int *nxv);

void dcuperpd13_(double *dcu, double *amu, int *nx, int *nxe);

void adcuperpdx13_(double complex *dcu, double complex *amu, int *nx,
                   int *nxv);

void adcuperpd13_(double *dcu, double *amu, int *nx, int *nxe);

void epoisdx13_(double complex *dcu, double complex *eyz, int *isign,
                double complex *fff, double *ax, double *affp,
                double *wp0, double *ci, double *wf, int *nx, int *nxv,
                int *nxd);

void epoisd13_(double *dcu, double *eyz, int *isign,
               double complex *fff, double *ax, double *affp,
               double *wp0, double *ci, double *wf, int *nx, int *nxe,
               int *nxv);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void lcguard1(double byz[], int nx, int nxe) {
   lcguard1_(byz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void lbguard1(double fxyz[], int nx, int nxe) {
   lbguard1_(fxyz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void ldguard1(double fx[], int nx, int nxe) {
   ldguard1_(fx,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void lscguard1(double cu[], double yj0, double zj0, int nx, int ngx,
               int nxe) {
   lscguard1_(cu,&yj0,&zj0,&nx,&ngx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void lsguard1(double q[], double qi0, int nx, int ngx, int nxe) {
   lsguard1_(q,&qi0,&nx,&ngx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void lacguard1(double cu[], int nx, int nxe) {
   lacguard1_(cu,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void laguard1(double q[], int nx, int nxe) {
   laguard1_(q,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void lscguard1l(double cu[], double yj0, double zj0, int nx, int ngx,
                int nxe) {
   lscguard1l_(cu,&yj0,&zj0,&nx,&ngx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void lsguard1l(double q[], double qi0, int nx, int ngx, int nxe) {
   lsguard1l_(q,&qi0,&nx,&ngx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void dblsin1a(double cu[], double cu2[], int nx, int nxv, int nx2v) {
   dblsin1a_(cu,cu2,&nx,&nxv,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void dblsin1d(double q[], double q2[], int nx, int nxv, int nx2v) {
   dblsin1d_(q,q2,&nx,&nxv,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void dblcos1a(double cu[], double cu2[], int nx, int nxv, int nx2v) {
   dblcos1a_(cu,cu2,&nx,&nxv,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void dblcos1d(double q[], double q2[], int nx, int nxv, int nx2v) {
   dblcos1d_(q,q2,&nx,&nxv,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void hafdbl1c(double byz[], double byz2[], int nx, int nxe, int nx2v) {
   hafdbl1c_(byz,byz2,&nx,&nxe,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void hafdbl1b(double fxyz[], double fxyz2[], int nx, int nxe, int nx2v) {
   hafdbl1b_(fxyz,fxyz2,&nx,&nxe,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void hafdbl1d(double q[], double q2[], int nx, int nxe, int nx2v) {
   hafdbl1d_(q,q2,&nx,&nxe,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void poisdx1(double q[], double fx[], int isign, double ffd[],
             double ax, double affp, double *we, int nx) {
   poisdx1_(q,fx,&isign,ffd,&ax,&affp,we,&nx);
   return;
}

/*--------------------------------------------------------------------*/
void poisd1(double q[], double fx[], int isign, double ffd[], double ax,
            double affp, double *we, int nx, int nxe, int nx2v) {
   poisd1_(q,fx,&isign,ffd,&ax,&affp,we,&nx,&nxe,&nx2v);
   return;
}

/*--------------------------------------------------------------------*/
void bpoisdx13(double complex cu[], double complex byz[], int isign,
               double complex ffd[], double ax, double affp, double ci,
               double *wm, int nx, int nxv, int nxd) {
   bpoisdx13_(cu,byz,&isign,ffd,&ax,&affp,&ci,wm,&nx,&nxv,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void bpoisd13(double cu[], double byz[], int isign,
              double complex ffd[], double ax, double affp, double ci,
              double *wm, int nx, int nxe, int nxv) {
   bpoisd13_(cu,byz,&isign,ffd,&ax,&affp,&ci,wm,&nx,&nxe,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void ibpoisdx13(double complex cu[], double complex byz[],
                double complex ffd[], double ci, double *wm, int nx,
                int nxv, int nxd) {
   ibpoisdx13_(cu,byz,ffd,&ci,wm,&nx,&nxv,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void ibpoisd13(double cu[], double byz[], double complex ffd[],
                double ci, double *wm, int nx, int nxe, int nxv) {
   ibpoisd13_(cu,byz,ffd,&ci,wm,&nx,&nxe,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void maxweldx1(double complex eyz[], double complex byz[],
               double complex cu[], double complex ffd[], double ci,
               double dt, double *wf, double *wm, int nx, int nxv,
               int nxd) {
   maxweldx1_(eyz,byz,cu,ffd,&ci,&dt,wf,wm,&nx,&nxv,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void maxweld1(double eyz[], double byz[], double cu[],
              double complex ffd[], double ci, double dt, double *wf,
              double *wm, int nx, int nxe, int nxv) {
   maxweld1_(eyz,byz,cu,ffd,&ci,&dt,wf,wm,&nx,&nxe,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void dmfieldd1(double complex q2[], double q[], int nx, int nxv,
               int nxe) {
   dmfieldd1_(q2,q,&nx,&nxv,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cmfieldd1(double complex cu2[], double cu[], int nx, int nxv,
               int nxe) {
   cmfieldd1_(cu2,cu,&nx,&nxv,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void amfieldd1(double complex amu2[], double amu[], int nx, int nxv,
               int nxe) {
   amfieldd1_(amu2,amu,&nx,&nxv,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void emfieldd1(double complex fxyz[], double complex fx[], double eyz[],
               double complex ffd[], int nx, int nxv, int nxe, int nxd) {
   emfieldd1_(fxyz,fx,eyz,ffd,&nx,&nxv,&nxe,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void bmfieldd1(double complex fxyz[], double eyz[],
               double complex ffd[], int nx, int nxv, int nxe, int nxd) {
   bmfieldd1_(fxyz,eyz,ffd,&nx,&nxv,&nxe,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void avpotdx13(double complex byz[], double complex ayz[], int nx,
               int nxv) {
   avpotdx13_(byz,ayz,&nx,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void avpotd13(double byz[], double ayz[], int nx, int nxe) {
   avpotd13_(byz,ayz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void avrpotdx13(double complex ayz[], double complex byz[],
                double complex ffd[], double ci, int nx, int nxv,
                int nxd) {
   avrpotdx13_(ayz,byz,ffd,&ci,&nx,&nxv,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void avrpotd13(double ayz[], double byz[], double complex ffd[],
               double ci, int nx, int nxe, int nxd) {
   avrpotd13_(ayz,byz,ffd,&ci,&nx,&nxe,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void gtsmodes1(double pot[], double pott[], int nx, int it, int modesx,
               int nxe, int nt, int modesxd) {
   gtsmodes1_(pot,pott,&nx,&it,&modesx,&nxe,&nt,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void ptsmodes1(double pot[], double pott[], int nx, int it, int modesx,
               int nxe, int nt, int modesxd) {
   ptsmodes1_(pot,pott,&nx,&it,&modesx,&nxe,&nt,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void gtvsmodes1(double vpot[], double vpott[], int nx, int it,
                int modesx, int ndim, int nxv, int nt, int modesxd) {
   gtvsmodes1_(vpot,vpott,&nx,&it,&modesx,&ndim,&nxv,&nt,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void ptvsmodes1(double vpot[], double vpott[], int nx, int it,
                int modesx, int ndim, int nxv, int nt, int modesxd) {
   ptvsmodes1_(vpot,vpott,&nx,&it,&modesx,&ndim,&nxv,&nt,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void lscfguard1(double cus[], double cu[], double q2m0, int nx,
                int nxe) {
   lscfguard1_(cus,cu,&q2m0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void lscfguard1l(double cus[], double cu[], double q2m0, int nx, 
                int nxe) {
   lscfguard1l_(cus,cu,&q2m0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void dcuperpdx13(double complex dcu[], double complex amu[], int nx,
                 int nxv) {
   dcuperpdx13_(dcu,amu,&nx,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void dcuperpd13(double dcu[], double amu[], int nx, int nxe) {
   dcuperpd13_(dcu,amu,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void adcuperpdx13(double complex dcu[], double complex amu[], int nx,
                  int nxv) {
   adcuperpdx13_(dcu,amu,&nx,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void adcuperpd13(double dcu[], double amu[], int nx, int nxe) {
   adcuperpd13_(dcu,amu,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void epoisdx13(double complex dcu[], double complex eyz[], int isign,
               double complex fff[], double ax, double affp, double wp0,
               double ci, double *wf, int nx, int nxv, int nxd) {
   epoisdx13_(dcu,eyz,&isign,fff,&ax,&affp,&wp0,&ci,wf,&nx,&nxv,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void epoisd13(double dcu[], double eyz[], int isign,
              double complex fff[], double ax, double affp, double wp0,
              double ci, double *wf, int nx, int nxe, int nxv) {
   epoisd13_(dcu,eyz,&isign,fff,&ax,&affp,&wp0,&ci,wf,&nx,&nxe,&nxv);
   return;
}
