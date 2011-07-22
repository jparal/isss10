//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 library push1lib.f

public class push1d {

   public native static void
      dpost(double[] part, double[] q, int nop, double qm,
            double[] tdpost, int idimp, int inorder, int dopt);

   public native static void
      push(double[] part, double[] fx, int nop, double qbm, double dt,
           double[] ek, double[] tpush, int nx, int ipbc, int idimp,
           int inorder, int popt);

   public native static void
      sortp(double[] part, double[] pt, int[] ip, int nop, int[] npic,
            double[] tsort, int idimp, int inorder);

   public native static void
      dpostgl(double[] part, double[] q, int nop, double qm, int nx,
              int nxh, double[] tdpost, int idimp);

   public native static void
      pushgl(double[] part, double[] fx, int nop, double qbm, double dt,
             double[] ek, int nx, int nxh, double[] tpush, int idimp);
}
