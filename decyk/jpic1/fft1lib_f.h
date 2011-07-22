/* header file for fft1lib_f.c */

#include <complex.h>

void fft1rx(double f[], double complex t[], int isign, int mixup[], 
            double complex sct[], int indx, int nxd, int nxhd);

void fft1r2(double f[], double complex t[], int isign, int mixup[],
            double complex sct[], int indx, int nxd, int nxhd);

void fft1r3(double f[], double complex t[], int isign, int mixup[],
            double complex sct[], int indx, int nxd, int nxhd);

void fft1c(double complex f[], int isign, int mixup[],
           double complex sct[], int indx, int nxd, int nxhd);

void fst1rx(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd);

void fct1rx(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd);

void fst1r2(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd);

void fct1r2(double f[], double complex t[], int isign, int mixup[],
            double complex sctd[], int indx, int nxd, int nxhd);

void fcst1r3(double f[], double complex t[], int isign, int mixup[],
             double complex sctd[], int indx, int nxd, int nxhd);

void fdst1rx(double f[], double complex t[], int isign, int mixup[],
             double complex sctd2[], int indx, int nxd, int nxhd,
             int nx2d);

void fdct1rx(double f[], double complex t[], int isign, int mixup[],
             double complex sctd2[], int indx, int nxd, int nxhd,
             int nx2d);
