//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 library dpush1lib.f

public class dpush1d {

   public native static void
      dmjpost(double[] part, double[] amu, int nop, double qm,
              double[] tdcjpost, int idimp, int ndim, int inorder,
              int djopt);

   public native static void
      dcjpost(double[] part, double[] fxyz, double[] byz, double[] cu,
              double[] dcu, double[] amu, double omx, int nop,
              double qm, double qbm, double dt, double[] tdcjpost,
              int idimp,  int ndim, int inorder, int djopt);
}
