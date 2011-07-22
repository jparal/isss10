/* 1d PIC library for solving field equations */
/* Wrappers for calling the Fortran routines from a C main program */

#include "field1lib_f.h"

void cguard1_(double *byz, int *nx, int *nxe);

void bguard1_(double *fxyz, int *nx, int *nxe);

void dguard1_(double *fx, int *nx, int *nxe);

void scguard1_(double *cu, double *yj0, double *zj0, int *nx, int *nxe);

void sguard1_(double *q, double *qi0, int *nx, int *nxe);

void acguard1_(double *cu, int *nx, int *nxe);

void aguard1_(double *q, int *nx, int *nxe);

void cguard1l_(double *byz, int *nx, int *nxe);

void bguard1l_(double *fxyz, int *nx, int *nxe);

void dguard1l_(double *fx, int *nx, int *nxe);

void scguard1l_(double *cu, double *yj0, double *zj0, int *nx, int *nxe);

void sguard1l_(double *q, double *qi0, int *nx, int *nxe);

void acguard1l_(double *cu, int *nx, int *nxe);

void aguard1l_(double *q, int *nx, int *nxe);

void poisp1_(double *q, double *fx, int *isign, double *ffc, double *ax,
             double *affp, double *we, int *nx);

void bpois13_(double complex *cu, double complex *byz, int *isign,
              double complex *ffc, double *ax, double *affp, double *ci,
              double *wm, int *nx, int *nxvh, int *nxhd);

void ibpois13_(double complex *cu, double complex *byz,
               double complex *ffc, double *ci, double *wm, int *nx,
               int *nxvh, int *nxhd);

void maxwel1_(double complex *eyz, double complex *byz,
              double complex *cu, double complex *ffc, double *ci,
              double *dt, double *wf, double *wm, int *nx, int *nxvh,
              int *nxhd);

void emfield1_(double complex *fxyz, double complex *fx,
               double complex *eyz, double complex *ffc, int *nx,
               int *nxvh, int *nxhd);

void bmfield1_(double complex *fyz, double complex *eyz,
               double complex *ffc, int *nx, int *nxvh, int *nxhd);

void emfieldr1_(double *fxyz, double *fx, double *eyz,
                double complex *ffc, int *nx, int *nxe, int *nxd);

void bmfieldr1_(double *fyz, double *eyz, double complex *ffc, int *nx,
                int *nxe, int *nxd);

void avpot13_(double complex *byz, double complex *ayz, int *nx,
              int *nxvh);

void avrpot13_(double complex *ayz, double complex *byz,
               double complex *ffc, double *ci, int *nx, int *nxvh,
               int *nxhd);

void gtmodes1_(double *pot, double *pott, int *nx, int *it, int *modesx,
               int *nxe, int *nt2, int *modesxd);

void ptmodes1_(double *pot, double *pott, int *nx, int *it, int *modesx,
               int *nxe, int *nt2, int *modesxd);

void gtvmodes1_(double complex *vpot, double complex *vpott, int *nx,
                int *it, int *modesx, int *ndim, int *nxvh, int *nt,
                int *modesxd);

void ptvmodes1_(double complex *vpot, double complex *vpott, int *nx,
                int *it, int *modesx, int *ndim, int *nxvh, int *nt,
                int *modesxd);

void scfguard1_(double *cus, double *cu, double *q2m0, int *nx,
                int *nxe);

void scfguard1l_(double *cus, double *cu, double *q2m0, int *nx,
                 int *nxe);

void dcuperp13_(double complex *dcu, double complex *amu, int *nx,
                int *nxvh);

void adcuperp13_(double complex *dcu, double complex *amu, int *nx,
                 int *nxvh);

void epois13_(double complex *dcu, double complex *eyz, int *isign,
              double complex *ffe, double *ax, double *affp, double *wp0,
              double *ci, double *wf, int *nx, int *nxvh, int *nxhd);

void wpmxn1_(double *qe, double *qi0, double *qbme, double *wpmax,
             double *wpmin, int *nx, int *nxe);

void baddext1_(double *byz, double *omy, double *omz, int *nx, int *nxe);

void vrcopy1_(double *f, double *g, int *nx, int *ndim, int *nxv);

void vccopy1_(double complex *f, double complex *g, int *nx, int *ndim,
              int *nxv);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cguard1(double byz[], int nx, int nxe) {
   cguard1_(byz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void bguard1(double fxyz[], int nx, int nxe) {
   bguard1_(fxyz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void dguard1(double fx[], int nx, int nxe) {
   dguard1_(fx,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void scguard1(double cu[], double yj0, double zj0, int nx, int nxe) {
   scguard1_(cu,&yj0,&zj0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void sguard1(double q[], double qi0, int nx, int nxe) {
   sguard1_(q,&qi0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void acguard1(double cu[], int nx, int nxe) {
   acguard1_(cu,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void aguard1(double q[], int nx, int nxe) {
   aguard1_(q,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cguard1l(double byz[], int nx, int nxe) {
   cguard1l_(byz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void bguard1l(double fxyz[], int nx, int nxe) {
   bguard1l_(fxyz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void dguard1l(double fx[], int nx, int nxe) {
   dguard1l_(fx,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void scguard1l(double cu[], double yj0, double zj0, int nx, int nxe) {
   scguard1l_(cu,&yj0,&zj0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void sguard1l(double q[], double qi0, int nx, int nxe) {
   sguard1l_(q,&qi0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void acguard1l(double cu[], int nx, int nxe) {
   acguard1l_(cu,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void aguard1l(double q[], int nx, int nxe) {
   aguard1l_(q,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void poisp1(double q[], double fx[], int isign, double ffc[], double ax,
            double affp, double *we, int nx) {
   poisp1_(q,fx,&isign,ffc,&ax,&affp,we,&nx);
   return;
}

/*--------------------------------------------------------------------*/
void bpois13(double complex cu[], double complex byz[], int isign,
             double complex ffc[], double ax, double affp, double ci,
             double *wm, int nx, int nxvh, int nxhd) {
   bpois13_(cu,byz,&isign,ffc,&ax,&affp,&ci,wm,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void ibpois13(double complex cu[], double complex byz[],
              double complex ffc[], double ci, double *wm, int nx,
              int nxvh, int nxhd) {
   ibpois13_(cu,byz,ffc,&ci,wm,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void maxwel1(double complex eyz[], double complex byz[],
             double complex cu[], double complex ffc[], double ci,
             double dt, double *wf, double *wm, int nx, int nxvh,
             int nxhd) {
   maxwel1_(eyz,byz,cu,ffc,&ci,&dt,wf,wm,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void emfield1(double complex fxyz[], double complex fx[],
              double complex eyz[], double complex ffc[], int nx,
              int nxvh, int nxhd) {
   emfield1_(fxyz,fx,eyz,ffc,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void bmfield1(double complex fyz[], double complex eyz[],
              double complex ffc[], int nx, int nxvh, int nxhd) {
   bmfield1_(fyz,eyz,ffc,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void emfieldr1(double fxyz[], double fx[], double eyz[],
               double complex ffc[], int nx, int nxe, int nxd) {
   emfieldr1_(fxyz,fx,eyz,ffc,&nx,&nxe,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void bmfieldr1(double fyz[], double eyz[], double complex ffc[], int nx,
               int nxe, int nxd) {
   bmfieldr1_(fyz,eyz,ffc,&nx,&nxe,&nxd);
   return;
}

/*--------------------------------------------------------------------*/
void avpot13(double complex byz[], double complex ayz[], int nx,
             int nxvh) {
   avpot13_(byz,ayz,&nx,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void avrpot13(double complex ayz[], double complex byz[],
              double complex ffc[], double ci, int nx, int nxvh,
              int nxhd) {
   avrpot13_(ayz,byz,ffc,&ci,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void gtmodes1(double pot[], double pott[], int nx, int it, int modesx,
              int nxe, int nt2, int modesxd) {
   gtmodes1_(pot,pott,&nx,&it,&modesx,&nxe,&nt2,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void ptmodes1(double pot[], double pott[], int nx, int it, int modesx,
              int nxe, int nt2, int modesxd) {
   ptmodes1_(pot,pott,&nx,&it,&modesx,&nxe,&nt2,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void gtvmodes1(double complex vpot[], double complex vpott[], int nx,
               int it, int modesx, int ndim, int nxvh, int nt,
               int modesxd) {
   gtvmodes1_(vpot,vpott,&nx,&it,&modesx,&ndim,&nxvh,&nt,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void ptvmodes1(double complex vpot[], double complex vpott[], int nx,
               int it, int modesx, int ndim, int nxvh, int nt,
               int modesxd) {
   ptvmodes1_(vpot,vpott,&nx,&it,&modesx,&ndim,&nxvh,&nt,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void scfguard1(double cus[], double cu[], double q2m0, int nx,
               int nxe) {
   scfguard1_(cus,cu,&q2m0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void scfguard1l(double cus[], double cu[], double q2m0, int nx,
                int nxe) {
   scfguard1l_(cus,cu,&q2m0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void dcuperp13(double complex dcu[], double complex amu[], int nx,
               int nxvh) {
   dcuperp13_(dcu,amu,&nx,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void adcuperp13(double complex dcu[], double complex amu[], int nx,
                int nxvh) {
   adcuperp13_(dcu,amu,&nx,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void epois13(double complex dcu[], double complex eyz[], int isign,
             double complex ffe[], double ax, double affp, double wp0,
             double ci, double *wf, int nx, int nxvh, int nxhd) {
   epois13_(dcu,eyz,&isign,ffe,&ax,&affp,&wp0,&ci,wf,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void wpmxn1(double qe[], double qi0, double qbme, double *wpmax,
            double *wpmin, int nx, int nxe) {
   wpmxn1_(qe,&qi0,&qbme,wpmax,wpmin,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void baddext1(double byz[], double omy, double omz, int nx, int nxe) {
   baddext1_(byz,&omy,&omz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void vrcopy1(double f[], double g[], int nx, int ndim, int nxv) {
   vrcopy1_(f,g,&nx,&ndim,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void vccopy1(double complex f[], double complex g[], int nx, int ndim,
             int nxv) {
   vccopy1_(f,g,&nx,&ndim,&nxv);
   return;
}
