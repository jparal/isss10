//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 library bpush1lib.f

public class bpush1d {

   public native static void
      djpost(double[] part, double[] cu, int nop, double qm, double dt,
             double[] tdjpost, int nx, int ipbc, int idimp, int ndim,
             int inorder, int djopt);

   public native static void
      push3(double[] part, double[] fxyz, double[] byz, double omx,
            int nop, double qbm, double dt, double dtc, double[] ek,
            double[] tpush, int nx, int ipbc, int idimp, int ndim,
            int inorder, int popt);

   public native static void
      retard(double[] part, int nop, double dtc, int nx, int ipbc,
             int idimp);
}
