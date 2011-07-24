//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 library init1lib.f

public class init1d {

   public native static void
      distr(double[] part, int nstart, int nop, double vtx, double vdx,
            int npx, int nx, int ipbc, int idimp);

   public native static void
      distrh(double[] part, int nstart, int nop, double vtx, double vty,
             double vtz, double vdx, double vdy, double vdz, int npx,
             int nx, int ipbc, int idimp);

   public native static void
      fdistr(double[] part, int nstart, int nop, double ampx,
             double scalex, double shiftx, int npx, int nx, int ipbc,
             int ndpro, int idimp);

   public native static void
      vdistr(double[] part, int nstart, int nop, double vtx,
             double vdx, int idimp);

   public native static void
      vdistrh(double[] part, int nstart, int nop, double vtx,
              double vty, double vtz, double vdx, double vdy,
              double vdz, int idimp);

   public native static void
      bdistr(double[] part, double[] byz, int nop, double qbm,
             int nx, int ipbc, int idimp, int ndim, int inorder);

}
