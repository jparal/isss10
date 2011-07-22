/* header file for init1lib_f.c */

void distr1(double part[], double vtx, double vdx, int npx, int idimp,
            int nop, int nx, int ipbc);

void distr1h(double part[], double vtx, double vty, double vtz,
             double vdx, double vdy, double vdz, int npx, int idimp,
             int nop, int nx, int ipbc);

void fdistr1(double part[], double (*fnx_)(), double argx1,
             double argx2, double argx3, int npx, int idimp, int nop,
             int nx, int *ierr);

double fldistr1_(double *x, double *anlx, double *anxi, double *shift,
                 int *intg);

double fsdistr1_(double *x, double *ans, double *dkx, double *phase,
                 int *intg);

double fgdistr1_(double *x, double *ang, double *wi, double *x0,
                 int *intg);

double fhdistr1_(double *x, double *anh, double *wi, double *x0,
                 int *intg);

void vdistr1(double part[], double vtx, double vdx, int idimp, int nop);

void vdistr1h(double part[], double vtx, double vty, double vtz,
              double vdx, double vdy, double vdz, int idimp, int nop);

void gbdistr1l(double part[], double byz[], double qbm, int idimp,
               int nop, int nx, int nxv, int ipbc);

double ranorm();

double randum();

double erfn() ;
