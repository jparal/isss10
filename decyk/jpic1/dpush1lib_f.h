/* header file for dpush1lib_f.c */

void gmjpost1(double part[], double amu[], double qm, int nop,
              int idimp, int nxv);

void gsmjpost1(double part[], double amu[], double qm, int nop,
               int idimp, int nxv);

void gdcjpost1(double part[], double fxyz[], double byz[], double cu[],
               double dcu[], double amu[], double omx, double qm,
               double qbm, double dt, int idimp, int nop, int nxv);

void gsdcjpost1(double part[], double fxyz[], double byz[], double cu[],
                double dcu[], double amu[], double omx, double qm,
                double qbm, double dt, int idimp, int nop, int nxv);

void gmjpost1l(double part[], double amu[], double qm, int nop,
               int idimp, int nxv);

void gsmjpost1l(double part[], double amu[], double qm, int nop,
                int idimp, int nxv);

void gdcjpost1l(double part[], double fxyz[], double byz[], double cu[],
                double dcu[], double amu[], double omx, double qm,
                double qbm, double dt, int idimp, int nop, int nxv);

void gsdcjpost1l(double part[], double fxyz[], double byz[],
                 double cu[], double dcu[], double amu[], double omx,
                 double qm, double qbm, double dt, int idimp, int nop,
                 int nxv);
