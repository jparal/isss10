/* library for simulation post-processors */
/* Wrappers for calling the Fortran routines from a C main program */

#include <string.h>
#include "libanls1_f.h"

void menucr1_(char *code, char *cp, int *ip, double *ap, int *nc,
              int *ncc, int *nci, int *ncr, int *istyle, int *irc);

void menucrv1_(char *code, char *cp, int *ip, double *ap, int *nc,
               int *ncc, int *nci, int *ncr, int *istyle, int *irc);

void helpcr1_(int *nvar, double *ax, double *ay, double *space);

void helpcrv1_(int *nvar, double *ax, double *ay, double *space);

void wpcorr1_(const char *runid, int *indx, int *ntp, int *modes,
              int *psolve, double *t0, double *tend, double *dt,
              int *lts, int *its, int *nts, int *kmin, int *kmax,
              int *ntd, int *ntc, double *wmin, double *wmax,
              double *dw, int *irc, int l);

void wpcorrv1_(const char *runid, int *indx, int *ntp, int *modes,
               int *psolve, double *t0, double *tend, double *dt,
               int *lts, int *its, int *nts, int *kmin, int *kmax,
               int *ntd, int *ntc, int *nvf, double *wmin, double *wmax,
               double *dw, int *irc, int l);

void corre_(double *pot, double complex *f, double complex *g,
            double *c, int *mixup, double complex *sct, int *inft,
            int *nt, int *ntc, int *nt2, int *nft, int *nfth,
            int *ntc2);

void spect_(double *f, double *wm, double *p, double *t0, double *dt,
            int *nt, int *iw, int *nt2, int *iw2);

void xspect_(double *f, double *wm, double *p, double *t0, double *dt,
             int *nt2, int *iw, int *iw2);

void sak1_(double *vpk, int *nx, int *modesx, int *modesxd);

void sq2pot1_(double complex *pott, double *potr, int *nx, int *modesx,
              int *modesxd);

void wel1_(double complex *pott, double *potb, double *affp, double *we,
           int *nx, int *modesx, int *modesxd);

void sq2vpot1_(double complex *vpott, double *vpotr, int *nx,
               int *modesx, int *modesxd);

void wbr1_(double complex *vpott, double *vpotb, double *affp,
           double *ci, double *wm, int *nx, int *modesx, int *modesxd);

void wer1_(double complex *vpotp, double *vpote, double *affp,
           double *wf, int *nx, int *modesx, int *modesxd);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void menucr1(char *code, char *cp, int ip[], double ap[], int nc,
             int ncc, int nci, int ncr, int istyle, int *irc) {
   menucr1_(code,cp,ip,ap,&nc,&ncc,&nci,&ncr,&istyle,irc);
   return;
}

/*--------------------------------------------------------------------*/
void menucrv1(char *code, char *cp, int ip[], double ap[], int nc,
              int ncc, int nci, int ncr, int istyle, int *irc) {
   menucrv1_(code,cp,ip,ap,&nc,&ncc,&nci,&ncr,&istyle,irc);
   return;
}

/*--------------------------------------------------------------------*/
void helpcr1(int nvar, double ax, double *ay, double space) {
   helpcr1_(&nvar,&ax,ay,&space);
   return;
}

/*--------------------------------------------------------------------*/
void helpcrv1(int nvar, double ax, double *ay, double space) {
   helpcrv1_(&nvar,&ax,ay,&space);
   return;
}

/*--------------------------------------------------------------------*/
void wpcorr1(const char *runid, int indx, int ntp, int modes,
             int psolve, double t0, double tend, double dt, int lts,
             int its, int nts, int kmin, int kmax, int ntd, int ntc,
             double wmin, double wmax, double dw, int *irc) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   wpcorr1_(runid,&indx,&ntp,&modes,&psolve,&t0,&tend,&dt,&lts,&its,
            &nts,&kmin,&kmax,&ntd,&ntc,&wmin,&wmax,&dw,irc,
            strlen(runid));
   return;
}

/*--------------------------------------------------------------------*/
void wpcorrv1(const char *runid, int indx, int ntp, int modes,
              int psolve, double t0, double tend, double dt, int lts,
              int its, int nts, int kmin, int kmax, int ntd, int ntc,
              int nvf, double wmin, double wmax, double dw, int *irc) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   wpcorrv1_(runid,&indx,&ntp,&modes,&psolve,&t0,&tend,&dt,&lts,&its,
             &nts,&kmin,&kmax,&ntd,&ntc,&nvf,&wmin,&wmax,&dw,irc,
             strlen(runid));
   return;
}

/*--------------------------------------------------------------------*/
void corre(double pot[], double complex f[], double complex g[],
           double c[], int mixup[], double complex sct[], int inft,
           int nt, int ntc, int nt2, int nft, int nfth,
           int ntc2) {
   corre_(pot,f,g,c,mixup,sct,&inft,&nt,&ntc,&nt2,&nft,&nfth,&ntc2);
   return;
}

/*--------------------------------------------------------------------*/
void spect(double f[], double wm[], double p[], double t0, double dt,
           int nt, int iw, int nt2, int iw2) {
   spect_(f,wm,p,&t0,&dt,&nt,&iw,&nt2,&iw2);
   return;
}

/*--------------------------------------------------------------------*/
void xspect(double f[], double wm[], double p[], double t0, double dt,
            int nt2, int iw, int iw2) {
   xspect_(f,wm,p,&t0,&dt,&nt2,&iw,&iw2);
   return;
}

/*--------------------------------------------------------------------*/
void sak1(double vpk[], int nx, int modesx, int modesxd) {
   sak1_(vpk,&nx,&modesx,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void sq2pot1(double complex pott[], double potr[], int nx, int modesx,
             int modesxd) {
   sq2pot1_(pott,potr,&nx,&modesx,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void wel1(double complex pott[], double potb[], double affp, double *we,
          int nx, int modesx, int modesxd) {
   wel1_(pott,potb,&affp,we,&nx,&modesx,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void sq2vpot1(double complex vpott[], double vpotr[], int nx,
              int modesx, int modesxd) {
   sq2vpot1_(vpott,vpotr,&nx,&modesx,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void wbr1(double complex vpott[], double vpotb[], double affp,
          double ci, double *wm, int nx, int modesx, int modesxd) {
   wbr1_(vpott,vpotb,&affp,&ci,wm,&nx,&modesx,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void wer1(double complex vpotp[], double vpote[], double affp,
          double *wf, int nx, int modesx, int modesxd) {
   wer1_(vpotp,vpote,&affp,wf,&nx,&modesx,&modesxd);
   return;
}
