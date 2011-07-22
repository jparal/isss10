/* header file for libanls1_f.c */

#include <complex.h>

void menucr1(char *code[], char *cp[], int ip[], double ap[], int nc, 
             int ncc, int nci, int ncr, int istyle, int *irc);

void menucrv1(char *code[], char *cp[], int ip[], double ap[], int nc,
              int ncc, int nci, int ncr, int istyle, int *irc);

void helpcr1(int nvar, double ax, double *ay, double space);

void helpcrv1(int nvar, double ax, double *ay, double space);

void wpcorr1(const char *runid, int indx, int ntp, int modes,
             int psolve, double t0, double tend, double dt, int lts,
             int its, int nts, int kmin, int kmax, int ntd, int ntc,
             double wmin, double wmax, double dw, int *irc);

void wpcorrv1(const char *runid, int indx, int ntp, int modes,
              int psolve, double t0, double tend, double dt, int lts,
              int its, int nts, int kmin, int kmax, int ntd, int ntc,
              int nvf, double wmin, double wmax, double dw, int *irc);

void corre(double pot[], double complex f[], double complex g[],
           double c[], int mixup[], double complex sct[], int inft,
           int nt, int ntc, int nt2, int nft, int nfth,
           int ntc2);

void spect(double f[], double wm[], double p[], double t0, double dt,
           int nt, int iw, int nt2, int iw2);

void xspect(double f[], double wm[], double p[], double t0, double dt,
            int nt2, int iw, int iw2);

void sak1(double vpk[], int nx, int modesx, int modesxd) ;

void sq2pot1(double complex pott[], double potr[], int nx, int modesx,
             int modesxd);

void wel1(double complex pott[], double potb[], double affp, double *we,
          int nx, int modesx, int modesxd);

void sq2vpot1(double complex vpott[], double vpotr[], int nx,
              int modesx, int modesxd);

void wbr1(double complex vpott[], double vpotb[], double affp,
          double ci, double *wm, int nx, int modesx, int modesxd);

void wer1(double complex vpotp[], double vpote[], double affp,
          double *wf, int nx, int modesx, int modesxd);
