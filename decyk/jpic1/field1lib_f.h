/* header file for field1lib_f.c */

#include <complex.h>

void cguard1(double byz[], int nx, int nxe);

void bguard1(double fxyz[], int nx, int nxe);

void dguard1(double fx[], int nx, int nxe);

void scguard1(double cu[], double yj0, double zj0, int nx, int nxe);

void sguard1(double q[], double qi0, int nx, int nxe);

void acguard1(double cu[], int nx, int nxe);

void aguard1(double q[], int nx, int nxe);

void cguard1l(double byz[], int nx, int nxe);

void bguard1l(double fxyz[], int nx, int nxe);

void dguard1l(double fx[], int nx, int nxe);

void scguard1l(double cu[], double yj0, double zj0, int nx, int nxe);

void sguard1l(double q[], double qi0, int nx, int nxe);

void acguard1l(double cu[], int nx, int nxe);

void aguard1l(double q[], int nx, int nxe);

void poisp1(double q[], double fx[], int isign, double ffc[], double ax,
            double affp, double *we, int nx);

void bpois13(double complex cu[], double complex byz[], int isign,
             double complex ffc[], double ax, double affp, double ci,
             double *wm, int nx, int nxvh, int nxhd);

void ibpois13(double complex cu[], double complex byz[],
              double complex ffc[], double ci, double *wm, int nx,
              int nxvh, int nxhd);

void maxwel1(double complex eyz[], double complex byz[],
             double complex cu[], double complex ffc[], double ci,
             double dt, double *wf, double *wm, int nx, int nxvh,
             int nxhd);

void emfield1(double complex fxyz[], double complex fx[],
              double complex eyz[], double complex ffc[], int nx,
              int nxvh, int nxhd);

void bmfield1(double complex fyz[], double complex eyz[],
              double complex ffc[], int nx, int nxvh, int nxhd);

void emfieldr1(double fxyz[], double fx[], double eyz[],
               double complex ffc[], int nx, int nxe, int nxd);

void bmfieldr1(double fyz[], double eyz[], double complex ffc[], int nx,
               int nxe, int nxd);

void avpot13(double complex byz[], double complex ayz[], int nx,
             int nxvh);

void avrpot13(double complex ayz[], double complex byz[],
              double complex ffc[], double ci, int nx, int nxvh,
              int nxhd);

void gtmodes1(double pot[], double pott[], int nx, int it, int modesx,
              int nxe, int nt2, int modesxd);

void ptmodes1(double pot[], double pott[], int nx, int it, int modesx,
              int nxe, int nt2, int modesxd);

void gtvmodes1(double complex vpot[], double complex vpott[], int nx,
               int it, int modesx, int ndim, int nxvh, int nt,
               int modesxd);

void ptvmodes1(double complex vpot[], double complex vpott[], int nx,
               int it, int modesx, int ndim, int nxvh, int nt,
               int modesxd);

void scfguard1(double cus[], double cu[], double q2m0, int nx,
               int nxe);

void scfguard1l(double cus[], double cu[], double q2m0, int nx,
                int nxe);

void dcuperp13(double complex dcu[], double complex amu[], int nx,
               int nxvh);

void adcuperp13(double complex dcu[], double complex amu[], int nx,
                int nxvh);

void epois13(double complex dcu[], double complex eyz[], int isign,
             double complex ffe[], double ax, double affp, double wp0,
             double ci, double *wf, int nx, int nxvh, int nxhd);

void wpmxn1(double qe[], double qi0, double qbme, double *wpmax,
            double *wpmin, int nx, int nxe);

void baddext1(double byz[], double omy, double omz, int nx, int nxe);

void vrcopy1(double f[], double g[], int nx, int ndim, int nxv);

void vccopy1(double complex f[], double complex g[], int nx, int ndim,
             int nxv);
