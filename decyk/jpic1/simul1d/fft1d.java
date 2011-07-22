//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 libraries fft1lib.f

public class fft1d {

   public native static void
      fft_init(int[] mixup, double[] sct, int indx);

   public native static void
      fft1(double[] f, int isign, int[] mixup, double[] sct,
          double[] tfft, int indx, int order);

   public native static void
      fft2(double[] f, int isign, int[] mixup, double[] sct,
      double[] tfft, int indx, int ndim, int order);

   public native static void
      fst_init(int[] mixup, double[] sctd, int indx);

   public native static void
      fst1(double[] f, int isign, int[] mixup, double[] sctd,
      double[] tfft, int indx, int order);

   public native static void
      fct1(double[] f, int isign, int[] mixup, double[] sctd,
      double[] tfft, int indx, int order);

   public native static void
      fst2(double[] f, int isign, int[] mixup, double[] sctd,
      double[] tfft, int indx, int ndim, int order);

   public native static void
      fct2(double[] f, int isign, int[] mixup, double[] sctd,
      double[] tfft, int indx, int ndim, int order);

   public native static void
      fcst(double[] f, int isign, int[] mixup, double[] sctd,
      double[] tfft, int indx, int ndim, int order); 
}
