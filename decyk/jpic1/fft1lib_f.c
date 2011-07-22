/*  1d PIC library for fast fourier transforms */
/* Wrappers for calling the Fortran routines from a C main program */

#include "fft1lib_f.h"

void fft1rx_(double *f, double complex *t, int *isign, int *mixup, 
             double complex *sct, int *indx, int *nxd, int *nxhd);

void fft1r2_(double *f, double complex *t, int *isign, int *mixup,
             double complex *sct, int *indx, int *nxd, int *nxhd);

void fft1r3_(double *f, double complex *t, int *isign, int *mixup,
             double complex *sct, int *indx, int *nxd, int *nxhd);

void fft1c_(double complex *f, int *isign, int *mixup,
            double complex *sct, int *indx, int *nxd, int *nxhd);

void fst1rx_(double *f, double complex *t, int *isign, int *mixup,
             double complex *sctd, int *indx, int *nxd, int *nxhd);

void fct1rx_(double *f, double complex *t, int *isign, int *mixup,
             double complex *sctd, int *indx, int *nxd, int *nxhd);

void fst1r2_(double *f, double complex *t, int *isign, int *mixup,
             double complex *sctd, int *indx, int *nxd, int *nxhd);

void fct1r2_(double *f, double complex *t, int *isign, int *mixup,
             double complex *sctd, int *indx, int *nxd, int *nxhd);
             
void fcst1r3_(double *f, double complex *t, int *isign, int *mixup,
              double complex *sctd, int *indx, int *nxd, int *nxhd);

void fdst1rx_(double *f, double complex *t, int *isign, int *mixup,
              double complex *sctd2, int *indx, int *nxd, int *nxhd,
              int *nx2d);

void fdct1rx_(double *f, double complex *t, int *isign, int *mixup,
              double complex *sctd2, int *indx, int *nxd, int *nxhd,
              int *nx2d);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void fft1rx(double f[], double complex t[], int isign, int mixup[], 
            double complex sct[], int indx, int nxd, int nxhd) {
   fft1rx_(f,t,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fft1r2(double f[], double complex t[], int isign, int mixup[],
            double complex sct[], int indx, int nxd, int nxhd) {
   fft1r2_(f,t,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fft1r3(double f[], double complex t[], int isign, int mixup[],
            double complex sct[], int indx, int nxd, int nxhd) {
   fft1r3_(f,t,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fft1c(double complex f[], int isign, int mixup[],
           double complex sct[], int indx, int nxd, int nxhd) {
   fft1c_(f,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fst1rx(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd) {
   fst1rx_(f,t,&isign,mixup,sctd,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fct1rx(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd) {
   fct1rx_(f,t,&isign,mixup,sctd,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fst1r2(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd) {
   fst1r2_(f,t,&isign,mixup,sctd,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fct1r2(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd) {
   fct1r2_(f,t,&isign,mixup,sctd,&indx,&nxd,&nxhd);
      return;
}

/*--------------------------------------------------------------------*/
void fcst1r3(double f[], double complex t[], int isign, int mixup[],
             double complex sctd[], int indx, int nxd, int nxhd) {
   fcst1r3_(f,t,&isign,mixup,sctd,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void fdst1rx(double f[], double complex t[], int isign, int mixup[],
             double complex sctd2[], int indx, int nxd, int nxhd,
             int nx2d) {
   fdst1rx_(f,t,&isign,mixup,sctd2,&indx,&nxd,&nxhd,&nx2d);
   return;
}

/*--------------------------------------------------------------------*/
void fdct1rx(double f[], double complex t[], int isign, int mixup[],
             double complex sctd2[], int indx, int nxd, int nxhd,
             int nx2d) {
   fdct1rx_(f,t,&isign,mixup,sctd2,&indx,&nxd,&nxhd,&nx2d);
   return;
}
