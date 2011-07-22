/* 1d PIC library for additional diagnostics */
/* Wrappers for calling the Fortran routines from a C main program */

#include <string.h>
#include "diag1lib_f.h"

void vdist1_(double *part, double *fv, double *fvm, int *idimp, int *np,
             int *nmv, int *nmvf);

void vdist13_(double *part, double *fv, double *fvm, int *idimp,
              int *np, int *nmv, int *nmvf);

void fwrite1_(double *f, int *nx, int *nxv, int *iunit, int *nrec,
             int *lrec, const char *name, int l);

void fread1_(double *f, int *nx, int *nxv, int *iunit, int *nrec,
             int *lrec, const char *name, int *ierr, int l);

void fcwrite1_(double *f, int *nx, int *nxv, int *iunit, int *nrec,
               int *lrec, const char *name, int l);

void fcread1_(double *f, int *nx, int *nxv, int *iunit, int *nrec,
              int *lrec, const char *name, int *ierr, int l);

void fwrite0_(double *f, int *nxp, int *iunit, int *nrec,
              const char *name, int l);

int fgetlrec_(double *f, int *nx);

int fgetclrec_(double complex *f, int *nx);

int fgetunit_(int *iunit);

void fbfopen_(double *f, int *nx, int *iunit, int *nrec,
              const char *fname, int l);

void fbfvopen_(double *f, int *nx,  int *ndim, int *iunit, int *nrec,
               const char *fname, int l);

void fbfcopen_(double complex *f, int *nx, int *iunit, int *nrec, 
               const char *fname, int l);

void fbfvcopen_(double complex *f, int *nx, int *ndim, int *iunit, 
                int *nrec, const char *fname, int l);

void frclose_(int *iunit);

void frwnml1_(const char *fname, const char *name, char *cn,
              double *ddata, char *sdata, int *iunit, int *isign,
              int *nl, int *ns, int *mcn, int *msd, int *ml, int *ms,
              int *ierr, int m, int n, int k, int l);
      
/* Interfaces to C */

/*--------------------------------------------------------------------*/
void vdist1(double part[], double fv[], double fvm[], int idimp,
            int np, int nmv, int nmvf) {
   vdist1_(part,fv,fvm,&idimp,&np,&nmv,&nmvf);
   return;
}

/*--------------------------------------------------------------------*/
void vdist13(double part[], double fv[], double fvm[], int idimp,
             int np, int nmv, int nmvf) {
   vdist13_(part,fv,fvm,&idimp,&np,&nmv,&nmvf);
   return;
}

/*--------------------------------------------------------------------*/
void fwrite1(double f[], int nx, int nxv, int *iunit, int *nrec,
             int lrec, const char *name) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fwrite1_(f,&nx,&nxv,iunit,nrec,&lrec,name,strlen(name));
   return;
}

/*--------------------------------------------------------------------*/
void fread1(double f[], int nx, int nxv, int iunit, int *nrec, int lrec,
            const char  *name, int *ierr) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fread1_(f,&nx,&nxv,&iunit,nrec,&lrec,name,ierr,strlen(name));
   return;
}

/*--------------------------------------------------------------------*/
void fcwrite1(double f[], int nx, int nxv, int *iunit, int *nrec,
              int lrec, const char *name) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fcwrite1_(f,&nx,&nxv,iunit,nrec,&lrec,name,strlen(name));
   return;
}

/*--------------------------------------------------------------------*/
void fcread1(double f[], int nx, int nxv, int iunit, int *nrec,
             int lrec, const char *name, int *ierr) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fcread1_(f,&nx,&nxv,&iunit,nrec,&lrec,name,ierr,strlen(name));
   return;
}

/*--------------------------------------------------------------------*/
void fwrite0(double f[], int nxp, int iunit, int *nrec,
             const char *name) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fwrite0_(f,&nxp,&iunit,nrec,name,strlen(name));
   return;
}

/*--------------------------------------------------------------------*/
int fgetlrec(double f[], int nx) {
   return fgetlrec_(f,&nx);
}

/*--------------------------------------------------------------------*/
int fgetclrec(double complex f[], int nx) {
   return fgetclrec_(f,&nx);
}

/*--------------------------------------------------------------------*/
int fgetunit(int iunit) {
   return fgetunit_(&iunit);
}

/*--------------------------------------------------------------------*/
void fbfopen(double f[], int nx, int iunit, int *nrec,
             const char *fname) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fbfopen_(f,&nx,&iunit,nrec,fname,strlen(fname));
   return;
}

/*--------------------------------------------------------------------*/
void fbfvopen(double f[], int nx,  int ndim, int iunit, int *nrec,
              const char *fname) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fbfvopen_(f,&nx,&ndim,&iunit,nrec,fname,strlen(fname));
   return;
}

/*--------------------------------------------------------------------*/
void fbfcopen(double complex f[], int nx, int iunit, int *nrec, 
              const char *fname) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fbfcopen_(f,&nx,&iunit,nrec,fname,strlen(fname));
   return;
}

/*--------------------------------------------------------------------*/
void fbfvcopen(double complex f[], int nx, int ndim, int iunit, 
               int *nrec, const char *fname) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   fbfvcopen_(f,&nx,&ndim,&iunit,nrec,fname,strlen(fname));
   return;
}

/*--------------------------------------------------------------------*/
void frclose(int iunit) {
   frclose_(&iunit);
   return;
}

/*--------------------------------------------------------------------*/
void frwnml1(const char *fname, const char *name, char *cn,
             double ddata[], char *sdata ,int iunit, int isign, int nl,
             int ns, int *mcn, int *msd, int *ml, int *ms, int *ierr) {
/* Passing characters this way may not work on all Fortran compilers */
/* OK with gfortran, g95, ibm, nag, absoft, intel                    */
   frwnml1_(fname,name,cn,ddata,sdata,&iunit,&isign,&nl,&ns,mcn,msd,ml,
            ms,ierr,strlen(fname),strlen(name),0,0);
   return;
}
