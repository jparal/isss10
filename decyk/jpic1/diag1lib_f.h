/* header file for diag1lib_f.c */

#include <complex.h>

void vdist1(double part[], double fv[], double fvm[], int idimp,
            int np, int nmv, int nmvf);

void vdist13(double part[], double fv[], double fvm[], int idimp,
             int np, int nmv, int nmvf);

void fwrite1(double f[], int nx, int nxv, int *iunit, int *nrec,
             int lrec, const char *name);

void fread1(double f[], int nx, int nxv, int iunit, int *nrec, int lrec,
            const char *name, int *ierr);

void fcwrite1(double f[], int nx, int nxv, int *iunit, int *nrec,
              int lrec, const char *name);

void fcread1(double f[], int nx, int nxv, int iunit, int *nrec,
             int lrec, const char *name, int *ierr);

void fwrite0(double f[], int nxp, int iunit, int *nrec,
             const char *name);

int fgetlrec(double f[], int nx);

int fgetclrec(double complex f[], int nx);

int fgetunit(int iunit);

void fbfopen(double f[], int nx, int iunit, int *nrec,
             const char *fname);

void fbfvopen(double f[], int nx,  int ndim, int iunit, int *nrec,
              const char *fname);

void fbfcopen(double complex f[], int nx, int iunit, int *nrec, 
              const char *fname);

void fbfvcopen(double complex f[], int nx, int ndim, int iunit, 
               int *nrec, const char *fname);

void frclose(int iunit);

void frwnml1(const char *fname, const char *name, char *cn,
             double ddata[], char *sdata ,int iunit, int isign, int nl,
             int ns, int *mcn, int *msd, int *ml, int *ms, int *ierr);
