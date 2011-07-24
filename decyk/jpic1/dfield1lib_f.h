/* header file for dfield1lib_f.c */

#include <complex.h>

void lcguard1(double byz[], int nx, int nxe);

void lbguard1(double fxyz[], int nx, int nxe);

void ldguard1(double fx[], int nx, int nxe);

void lscguard1(double cu[], double yj0, double zj0, int nx, int ngx,
               int nxe);

void lsguard1(double q[], double qi0, int nx, int ngx, int nxe);

void lacguard1(double cu[], int nx, int nxe);

void laguard1(double q[], int nx, int nxe);

void lscguard1l(double cu[], double yj0, double zj0, int nx, int ngx,
                int nxe);

void lsguard1l(double q[], double qi0, int nx, int ngx, int nxe);

void dblsin1a(double cu[], double cu2[], int nx, int nxv, int nx2v);

void dblsin1d(double q[], double q2[], int nx, int nxv, int nx2v);

void dblcos1a(double cu[], double cu2[], int nx, int nxv, int nx2v);

void dblcos1d(double q[], double q2[], int nx, int nxv, int nx2v);

void hafdbl1c(double byz[], double byz2[], int nx, int nxe, int nx2v);

void hafdbl1b(double fxyz[], double fxyz2[], int nx, int nxe, int nx2v);

void hafdbl1d(double q[], double q2[], int nx, int nxe, int nx2v);

void poisdx1(double q[], double fx[], int isign, double ffd[],
             double ax, double affp, double *we, int nx);

void poisd1(double q[], double fx[], int isign, double ffd[], double ax,
            double affp, double *we, int nx, int nxe, int nx2v);

void bpoisdx13(double complex cu[], double complex byz[], int isign,
               double complex ffd[], double ax, double affp, double ci,
               double *wm, int nx, int nxv, int nxd);

void bpoisd13(double cu[], double byz[], int isign,
              double complex ffd[], double ax, double affp, double ci,
              double *wm, int nx, int nxe, int nxv);

void ibpoisdx13(double complex cu[], double complex byz[],
                double complex ffd[], double ci, double *wm, int nx,
                int nxv, int nxd);

void ibpoisd13(double cu[], double byz[], double complex ffd[],
                double ci, double *wm, int nx, int nxe, int nxv);

void maxweldx1(double complex eyz[], double complex byz[],
               double complex cu[], double complex ffd[], double ci,
               double dt, double *wf, double *wm, int nx, int nxv,
               int nxd);

void maxweld1(double eyz[], double byz[], double cu[],
              double complex ffd[], double ci, double dt, double *wf,
              double *wm, int nx, int nxe, int nxv);

void dmfieldd1(double complex q2[], double q[], int nx, int nxv,
               int nxe);

void cmfieldd1(double complex cu2[], double cu[], int nx, int nxv,
               int nxe);

void amfieldd1(double complex amu2[], double amu[], int nx, int nxv,
               int nxe);

void emfieldd1(double complex fxyz[], double complex fx[], double eyz[],
               double complex ffd[], int nx, int nxv, int nxe, int nxd);

void bmfieldd1(double complex fxyz[], double eyz[],
               double complex ffd[], int nx, int nxv, int nxe, int nxd);

void avpotdx13(double complex byz[], double complex ayz[], int nx,
               int nxv);

void avpotd13(double byz[], double ayz[], int nx, int nxe);

void avrpotdx13(double complex ayz[], double complex byz[],
                double complex ffd[], double ci, int nx, int nxv,
                int nxd);

void avrpotd13(double ayz[], double byz[], double complex ffd[],
               double ci, int nx, int nxe, int nxd);

void gtsmodes1(double pot[], double pott[], int nx, int it, int modesx,
               int nxe, int nt, int modesxd);

void ptsmodes1(double pot[], double pott[], int nx, int it, int modesx,
               int nxe, int nt, int modesxd);

void gtvsmodes1(double vpot[], double vpott[], int nx, int it,
                int modesx, int ndim, int nxv, int nt, int modesxd);

void ptvsmodes1(double vpot[], double vpott[], int nx, int it,
                int modesx, int ndim, int nxv, int nt, int modesxd);

void lscfguard1(double cus[], double cu[], double q2m0, int nx,
                int nxe);

void lscfguard1l(double cus[], double cu[], double q2m0, int nx, 
                int nxe);

void dcuperpdx13(double complex dcu[], double complex amu[], int nx,
                 int nxv);

void dcuperpd13(double dcu[], double amu[], int nx, int nxe);

void adcuperpdx13(double complex dcu[], double complex amu[], int nx,
                  int nxv);

void adcuperpd13(double dcu[], double amu[], int nx, int nxe);

void epoisdx13(double complex dcu[], double complex eyz[], int isign,
               double complex fff[], double ax, double affp, double wp0,
               double ci, double *wf, int nx, int nxv, int nxd);

void epoisd13(double dcu[], double eyz[], int isign,
              double complex fff[], double ax, double affp, double wp0,
              double ci, double *wf, int nx, int nxe, int nxv);
