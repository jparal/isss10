/* header file for bpush1lib_f.c */

void gjpost1(double part[], double cu[], double qm, double dt,
             int nop, int idimp, int nx, int nxv, int ipbc);

void gsjpost1(double part[], double cu[], double qm, double dt,
              int nop, int idimp, int nx, int nxv, int ipbc);

void gsjpost1x(double part[], double cu[], double qm, double dt,
               int nop, int idimp, int nx, int nxv, int ipbc);

void gjpost1l(double part[], double cu[], double qm, double dt,
              int nop, int idimp, int nx, int nxv, int ipbc);

void gsjpost1l(double part[], double cu[], double qm, double dt,
               int nop, int idimp, int nx, int nxv, int ipbc);

void gsjpost1xl(double part[], double cu[], double qm, double dt,
                int nop, int idimp, int nx, int nxv, int ipbc);

void gbpush13(double part[], double fxyz[], double byz[], double omx,
              double qbm, double dt, double dtc, double *ek, int idimp,
              int nop, int nx, int nxv, int ipbc);

void gsbpush13(double part[], double fxyz[], double byz[], double omx,
               double qbm, double dt, double dtc, double *ek, int idimp,
               int nop, int nx, int nxv, int ipbc);

void gbpush13l(double part[], double fxyz[], double byz[], double omx,
               double qbm, double dt, double dtc, double *ek, int idimp,
               int nop, int nx, int nxv, int ipbc);

void gsbpush13l(double part[], double fxyz[], double byz[], double omx,
                double qbm, double dt, double dtc, double *ek,
                int idimp, int nop, int nx, int nxv, int ipbc);

void retard1(double part[], double dtc, int idimp, int nop, int nx,
             int ipbc);
