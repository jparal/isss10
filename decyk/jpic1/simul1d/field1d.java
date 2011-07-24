//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 library field1lib.f

public class field1d {

   public native static void
      cguard(double[] fxy, int nx, int ndim, int inorder);

   public native static void
      dguard(double[] fx, int nx, int inorder);

   public native static void
      scguard(double[] cu, double yj0, double zj0, int nx, int inorder);

   public native static void
      sguard(double[] q, double qi0, int nx, int inorder);

   public native static void
      acguard(double[] cu, int nx, int inorder);

   public native static void
      aguard(double[] q, int nx, int inorder);

   public native static void
      pois_init(double[] ffc, double ax, double affp, int nx);

   public native static void
      pois(double[] q, double[] fx, int isign, double[] ffc,
           double[] we, int nx, int inorder);

   public native static void
      bpois(double[] cu, double[] byz, double[] ffc, double ci,
            double[] wm, int nx, int ndim, int inorder);

   public native static void
      apois(double[] cu, double[] ayz, double[] ffc, double ci,
            double[] wm, int nx, int ndim, int inorder);

   public native static void
      ibpois(double[] cu, double[] byz, double[] ffc, double ci,
             double[] wm, int nx, int inorder);

   public native static void
      maxwel(double[] eyz, double[] byz, double[] cu, double[] ffc,
             double ci, double dt, double[] wf, double[] wm, int nx,
             int inorder);

   public native static void
      emfield(double[] fxyz, double[] fx, double[] eyz, double[] ffc,
              int nx, int inorder);

   public native static void
      bmfield(double[] fyz, double[] eyz, double[] ffc, int nx,
              int inorder);

   public native static void
      emfieldr(double[] fxyz, double[] fx, double[] eyz, double[] ffc,
               int nx, int inorder);

   public native static void
      bmfieldr(double[] fyz, double[] eyz, double[] ffc, int nx,
               int inorder);

   public native static void
      avpot(double[] byz, double[] ayz, int nx, int inorder);

   public native static void
      avrpot(double[] ayz, double[] byz, double[] ffc, double ci,
             int nx, int inorder);

   public native static void
      gtmodes(double[] pot, double[] pott, int nx, int modesx,
              int order);

   public native static void
      ptmodes(double[] pot, double[] pott, int nx, int modesx,
              int order);

   public native static void
      gtvmodes(double[] vpot, double[] vpott, int nx, int modesx,
               int ndim, int order);

   public native static void
      ptvmodes(double[] vpot, double[] vpott, int nx, int modesx,
               int ndim, int order);

   public native static void
      sfguard(double[] cus, double[] cu, double q2m0, int nx,
              int ndim, int inorder);

   public native static void
      dcuperp(double[] dcu, double[] amu, int nx, int ndim,
              int inorder);

   public native static void
      adcuperp(double[] dcu, double[] amu, int nx, int ndim,
               int inorder);

   public native static void
      epois_init(double[] ffe, double ax, double affp, double wp0,
                 double ci, int nx);

   public native static void
      epois(double[] dcu, double[] eyz, double[] ffe, double ci,
            double[] wf, int nx, int ndim, int inorder);

   public native static void
      iepois(double[] dcu, double[] eyz, double[] ffe, double ci,
             double[] wf, int nx, int ndim, int inorder);

   public native static void
      wpmxn(double[] qe, double qi0, double qbme, double[] wpmax,
            double[] wpmin, int nx, int inorder);

   public native static void
      addqei(double[] qe, double[] qi, int nx, int inorder);

   public native static void
      addqeix(double[] qe, double[] qi, double qbme, double qbmi,
              double[] wpmax, double[] wpmin, int nx, int inorder);

   public native static void
      baddext(double[] byz, double omy, double omz, int nx,
              int ndim, int inorder);

   public native static void
      ivrcopy(double[] f, double[] g, int nx, int ndim, int inorder);

   public native static void
      ivccopy(double[] f, double[] g, int nx, int ndim, int inorder);
}
