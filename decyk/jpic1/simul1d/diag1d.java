//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 library diag1lib.f

import static simul1d.input1d.*;

public class diag1d {

   public native static void
      writebf(double[] f, int nx, int[] iunit, int[] nrec, String name,
              int order);

   public native static void
      readbf(double[] f, int nx, int iunit, int[] nrec, int[] ierr,
             String name, int order);

   public native static void
      writebfv(double[] f, int nx, int[] iunit, int[] nrec, String name,
               int ndim, int order);

   public native static void
      readbfv(double[] f, int nx, int iunit, int[] nrec, int[] ierr,
              String name, int ndim, int order);

   public native static void
      writebfc(double[] f, int nx, int[] iunit, int[] nrec, String name,
               int order);

   public native static void
      readbfc(double[] f, int nx, int iunit, int[] nrec, int[] ierr,
              String name, int order);

   public native static void
      writebfvc(double[] f, int nx, int[] iunit, int[] nrec,
                String name, int ndim, int order);

   public native static void
      readbfvc(double[] f, int nx, int iunit, int[] nrec, int[] ierr,
               String name, int ndim, int order);

   public native static void
      writef(double[] f, int iunit, int[] nrec, String name);

   public native static void
      vdist(double[] part, double[] fv, double[] fvm, int np, int nmv,
            int idimp, int idimv);

   public native static int
      fgetunit(int iunit);

   public native static void
      fbfopen(double[] f, int nx, int iunit, int[] nrec, String fname);

   public native static void
      fbfvopen(double[] f, int nx,  int ndim, int iunit, int[] nrec,
               String fname);

   public native static void
      fbfcopen(double[] f, int nx, int iunit, int[] nrec, String fname);

   public native static void
      fbfvcopen(double[] f, int nx, int ndim, int iunit, int[] nrec,
                String fname);

   public native static void
      frclose(int iunit);

   public native static void
      frwnml1(String fname, String name, String[] cn, double[] ddata,
              String[] sdata, int iunit, int isign, int[] ml, int[] ms,
              int[] ierr);

   public static void jtimer(double[] time, long [] itime, int icntrl) {
// this subroutine performs timing
// input: icntrl, itime
// icntrl = (-1,0,1) = (initialize,ignore,read) clock
// clock should be initialized before it is read!
// time = elapsed time in seconds
      long jclock, nclock;
      double tick = 1.0/1000000000.0;
      if (icntrl==0) 
         return;
// read clock and write time difference from last clock initialization
// calculate time elapsed in nanoseconds                              */
      if (icntrl==1) {
         nclock = System.nanoTime();
         jclock = itime[0];
         itime[0] = nclock;
         time[0] = tick*(double) (nclock - jclock);
      }
      else {
// initialize clock              
         itime[0] = System.nanoTime();
      }
      return;
   }

   public static void jrwnml1(String fname, String name, int iunit,
                              int isign, int[] irc) {
// reads and writes Fortran namelist to and from input1d static data
//   fname = Fortran file name, name = Fortran namelist name
//   iunit = Fortran unit number, isign = (-1,1) = (read,write) namelist
//   ierr = error if not zero

      int nl = 52, ns = 1, sname = 0;
      int i = 0, k = 0, ierr = 0;
      String[] cn = new String[nl+ns];
      double[] ddata = new double[nl];
      String[] sdata = new String[ns];
      int[] ml = new int[1]; int[] ms = new int[1];

      if (name.equals("input1"))
         sname = 1;
      else if (name.equals("pot1d"))
         sname = 2;
      else if (name.equals("vpot1d"))
         sname = 3;
      else if (name.equals("em1d"))
         sname = 4;

      switch (sname) {
      case 1:

// write namelist from packed data
      if (isign > 0) {
// pack data and names
         ddata[i] = idrun; cn[i] = "idrun"; i++;
         ddata[i] = idrun0; cn[i] = "idrun0"; i++;
         ddata[i] = idcode; cn[i] = "idcode"; i++;
         ddata[i] = indx; cn[i] = "indx"; i++;
         ddata[i] = npx; cn[i] = "npx"; i++;
         ddata[i] = npxb; cn[i] = "npxb"; i++;
         ddata[i] = inorder; cn[i] = "inorder"; i++;
         ddata[i] = popt; cn[i] = "popt"; i++;
         ddata[i] = dopt; cn[i] = "dopt"; i++;
         ddata[i] = djopt;  cn[i] = "djopt"; i++;
         ddata[i] = nustrt; cn[i] = "nustrt"; i++;
         ddata[i] = ntr; cn[i] = "ntr"; i++;
         ddata[i] = ntw; cn[i] = "ntw"; i++;
         ddata[i] = ntp; cn[i] = "ntp"; i++;
         ddata[i] = nta; cn[i] = "nta"; i++;
         ddata[i] = ntv; cn[i] = "ntv"; i++;
         ddata[i] = nts; cn[i] = "nts"; i++;
         ddata[i] = nte; cn[i] = "nte"; i++;
//       ddata[i] = ndw; cn[i] = "ndw"; i++;
//       ddata[i] = ndp; cn[i] = "ndp"; i++;
//       ddata[i] = nda; cn[i] = "nda"; i++;
//       ddata[i] = ndv; cn[i] = "ndv"; i++;
//       ddata[i] = nds; cn[i] = "nds"; i++;
//       ddata[i] = nde; cn[i] = "nde"; i++;
         ddata[i] = tend; cn[i] = "tend"; i++;
         ddata[i] = dt; cn[i] = "dt"; i++;
         ddata[i] = qme; cn[i] = "qme"; i++;
         ddata[i] = vtx; cn[i] = "vtx"; i++;
         ddata[i] = vty; cn[i] = "vty"; i++;
         ddata[i] = vtz; cn[i] = "vtz"; i++;
         ddata[i] = vx0; cn[i] = "vx0"; i++;
         ddata[i] = vy0; cn[i] = "vy0"; i++;
         ddata[i] = vz0; cn[i] = "vz0"; i++;
         ddata[i] = vdx; cn[i] = "vdx"; i++;
         ddata[i] = vdy; cn[i] = "vdy"; i++;
         ddata[i] = vdz; cn[i] = "vdz"; i++;
         ddata[i] = vtdx; cn[i] = "vtdx"; i++;
         ddata[i] = vtdy; cn[i] = "vtdy"; i++;
         ddata[i] = vtdz; cn[i] = "vtdz"; i++;
         ddata[i] = psolve; cn[i] = "psolve"; i++;
//       ddata[i] = relativity; cn[i] = "relativity"; i++;
         ddata[i] = omx; cn[i] = "omx"; i++;
         ddata[i] = omy; cn[i] = "omy"; i++;
         ddata[i] = omz; cn[i] = "omz"; i++;
         ddata[i] = ci; cn[i] = "ci"; i++;
         ddata[i] = ax; cn[i] = "ax"; i++;
         ddata[i] = ndim; cn[i] = "ndim"; i++;
         ddata[i] = ndc; cn[i] = "ndc"; i++;
         ddata[i] = movion; cn[i] = "movion"; i++;
         ddata[i] = sortime; cn[i] = "sortime"; i++;
         ddata[i] = nplot; cn[i] = "nplot"; i++;
         ddata[i] = sntasks; cn[i] = "sntasks"; i++;
         ddata[i] = ndprof; cn[i] = "ndprof"; i++;
         ddata[i] = ampdx; cn[i] = "ampdx"; i++;
         ddata[i] = scaledx; cn[i] = "scaledx"; i++;
         ddata[i] = shiftdx; cn[i] = "shiftdx"; i++;
         ddata[i] = modesxp; cn[i] = "modesxp"; i++;
         ddata[i] = modesxa; cn[i] = "modesxa"; i++;
         ddata[i] = modesxe; cn[i] = "modesxe"; i++;
// set number of variables found
         ml[0] = i; ms[0] = k;
// write Fortran namelist from packed data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
      }
// read namelist from packed data
      else if (isign < 0) {
// read Fortran namelist and pack data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
// unpack data and check names
         idrun = (int) ddata[i];
            if (!(cn[i].equals("idrun"))) ierr = i+1; i++;
         idrun0 = (int) ddata[i];
            if (!(cn[i].equals("idrun0"))) ierr = i+1; i++;
         idcode = (int) ddata[i];
            if (!(cn[i].equals("idcode"))) ierr = i+1; i++;
         indx = (int) ddata[i];
            if (!(cn[i].equals("indx"))) ierr = i+1; i++;
         npx = (int) ddata[i];
            if (!(cn[i].equals("npx"))) ierr = i+1; i++;
         npxb = (int) ddata[i];
            if (!(cn[i].equals("npxb"))) ierr = i+1; i++;
         inorder = (int) ddata[i];
            if (!(cn[i].equals("inorder"))) ierr = i+1; i++;
         popt = (int) ddata[i];
            if (!(cn[i].equals("popt"))) ierr = i+1; i++;
         dopt = (int) ddata[i];
            if (!(cn[i].equals("dopt"))) ierr = i+1; i++;
         djopt = (int) ddata[i];
            if (!(cn[i].equals("djopt"))) ierr = i+1; i++;
         nustrt = (int) ddata[i];
            if (!(cn[i].equals("nustrt"))) ierr = i+1; i++;
         ntr = (int) ddata[i];
            if (!(cn[i].equals("ntr"))) ierr = i+1; i++;
         ntw = (int) ddata[i];
            if (!(cn[i].equals("ntw"))) ierr = i+1; i++;
         ntp = (int) ddata[i];
            if (!(cn[i].equals("ntp"))) ierr = i+1; i++;
         nta = (int) ddata[i];
            if (!(cn[i].equals("nta"))) ierr = i+1; i++;
         ntv = (int) ddata[i];
            if (!(cn[i].equals("ntv"))) ierr = i+1; i++;
         nts = (int) ddata[i];
            if (!(cn[i].equals("nts"))) ierr = i+1; i++;
         nte = (int) ddata[i];
            if (!(cn[i].equals("nte"))) ierr = i+1; i++;
//       ndw = (int) ddata[i];
//          if (!(cn[i].equals("ndw"))) ierr = i+1; i++;
//       ndp = (int) ddata[i];
//          if (!(cn[i].equals("ndp"))) ierr = i+1; i++;
//       nda = (int) ddata[i];
//          if (!(cn[i].equals("nda"))) ierr = i+1; i++;
//       ndv = (int) ddata[i];
//          if (!(cn[i].equals("ndv"))) ierr = i+1; i++;
//       nds = (int) ddata[i];
//          if (!(cn[i].equals("nds"))) ierr = i+1; i++;
//       nde = (int) ddata[i];
//          if (!(cn[i].equals("nde"))) ierr = i+1; i++;
         tend = ddata[i];
            if (!(cn[i].equals("tend"))) ierr = i+1; i++;
         dt = ddata[i];
            if (!(cn[i].equals("dt"))) ierr = i+1; i++;
         qme = ddata[i];
            if (!(cn[i].equals("qme"))) ierr = i+1; i++;
         vtx = ddata[i];
            if (!(cn[i].equals("vtx"))) ierr = i+1; i++;
         vty = ddata[i];
            if (!(cn[i].equals("vty"))) ierr = i+1; i++;
         vtz = ddata[i];
            if (!(cn[i].equals("vtz"))) ierr = i+1; i++;
         vx0 = ddata[i];
            if (!(cn[i].equals("vx0"))) ierr = i+1; i++;
         vy0 = ddata[i];
            if (!(cn[i].equals("vy0"))) ierr = i+1; i++;
         vz0 = ddata[i];
            if (!(cn[i].equals("vz0"))) ierr = i+1; i++;
         vdx = ddata[i];
            if (!(cn[i].equals("vdx"))) ierr = i+1; i++;
         vdy = ddata[i];
            if (!(cn[i].equals("vdy"))) ierr = i+1; i++;
         vdz = ddata[i];
            if (!(cn[i].equals("vdz"))) ierr = i+1; i++;
         vtdx = ddata[i];
            if (!(cn[i].equals("vtdx"))) ierr = i+1; i++;
         vtdy = ddata[i];
            if (!(cn[i].equals("vtdy"))) ierr = i+1; i++;
         vtdz = ddata[i];
            if (!(cn[i].equals("vtdz"))) ierr = i+1; i++;
         psolve = (int) ddata[i];
            if (!(cn[i].equals("psolve"))) ierr = i+1; i++;
//       relativity = (int) ddata[i];
//          if (!(cn[i].equals("relativity"))) ierr = i+1; i++;
         omx = ddata[i];
            if (!(cn[i].equals("omx"))) ierr = i+1; i++;
         omy = ddata[i];
            if (!(cn[i].equals("omy"))) ierr = i+1; i++;
         omz = ddata[i];
            if (!(cn[i].equals("omz"))) ierr = i+1; i++;
         ci = ddata[i];
            if (!(cn[i].equals("ci"))) ierr = i+1; i++;
         ax = ddata[i];
            if (!(cn[i].equals("ax"))) ierr = i+1; i++;
         ndim = (int) ddata[i];
            if (!(cn[i].equals("ndim"))) ierr = i+1; i++;
         ndc = (int) ddata[i];
            if (!(cn[i].equals("ndc"))) ierr = i+1; i++;
         movion = (int) ddata[i];
            if (!(cn[i].equals("movion"))) ierr = i+1; i++;
         sortime = (int) ddata[i];
            if (!(cn[i].equals("sortime"))) ierr = i+1; i++;
         nplot = (int) ddata[i];
            if (!(cn[i].equals("nplot"))) ierr = i+1; i++;
         sntasks = (int) ddata[i];
            if (!(cn[i].equals("sntasks"))) ierr = i+1; i++;
         ndprof = (int) ddata[i];
            if (!(cn[i].equals("ndprof"))) ierr = i+1; i++;
         ampdx = (int) ddata[i];
            if (!(cn[i].equals("ampdx"))) ierr = i+1; i++;
         scaledx = ddata[i];
            if (!(cn[i].equals("scaledx"))) ierr = i+1; i++;
         shiftdx = ddata[i];
            if (!(cn[i].equals("shiftdx"))) ierr = i+1; i++;
         modesxp = (int) ddata[i];
            if (!(cn[i].equals("modesxp"))) ierr = i+1; i++;
         modesxa = (int) ddata[i];
            if (!(cn[i].equals("modesxa"))) ierr = i+1; i++;
         modesxe = (int) ddata[i];
            if (!(cn[i].equals("modesxe"))) ierr = i+1; i++;
         if (ierr != 0)
            irc[0] = 1000+ierr;
      }
      break;

      case 2:

// write namelist from packed data
      if (isign > 0) {
// pack data and names
         ddata[i] = idrun; cn[i] = "idrun"; i++;
         ddata[i] = indx; cn[i] = "indx"; i++;
         ddata[i] = ntp; cn[i] = "ntp"; i++;
         ddata[i] = modesxp; cn[i] = "modesxp"; i++;
         ddata[i] = psolve; cn[i] = "psolve"; i++;
         ddata[i] = omx; cn[i] = "omx"; i++;
         ddata[i] = omy; cn[i] = "omy"; i++;
         ddata[i] = omz; cn[i] = "omz"; i++;
         ddata[i] = nprec[0]; cn[i] = "nprec"; i++;
         ddata[i] = t0; cn[i] = "t0"; i++;
         ddata[i] = tend; cn[i] = "tend"; i++;
         ddata[i] = dt; cn[i] = "dt"; i++;
         ddata[i] = ceng; cn[i] = "ceng"; i++;
         sdata[k] = fpname; cn[i+k] = "fpname"; k++;
// set number of variables found
         ml[0] = i; ms[0] = k;
// write Fortran namelist from packed data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
      }
// read namelist from data
      else if (isign < 0) {
// read Fortran namelist and pack data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
// unpack data and check names
         idrun = (int) ddata[i];
            if (!(cn[i].equals("idrun"))) ierr = i+1; i++;
         indx = (int) ddata[i]; 
            if (!(cn[i].equals("indx"))) ierr = i+1; i++;
         ntp = (int) ddata[i];
            if (!(cn[i].equals("ntp"))) ierr = i+1; i++;
         modesxp = (int) ddata[i];
            if (!(cn[i].equals("modesxp"))) ierr = i+1; i++;
         psolve = (int) ddata[i];
            if (!(cn[i].equals("psolve"))) ierr = i+1; i++;
         omx = ddata[i];
            if (!(cn[i].equals("omx"))) ierr = i+1; i++;
         omy = ddata[i];
            if (!(cn[i].equals("omy"))) ierr = i+1; i++;
         omz = ddata[i];
            if (!(cn[i].equals("omz"))) ierr = i+1; i++;
         nprec[0] = (int) ddata[i];
            if (!(cn[i].equals("nprec"))) ierr = i+1; i++;
         t0 = ddata[i];
            if (!(cn[i].equals("t0"))) ierr = i+1; i++;
         tend = ddata[i];
            if (!(cn[i].equals("tend"))) ierr = i+1; i++;
         dt = ddata[i];
            if (!(cn[i].equals("dt"))) ierr = i+1; i++;
         ceng = ddata[i];
            if (!(cn[i].equals("ceng"))) ierr = i+1; i++;
         fpname = sdata[k];
            if (!(cn[i+k].equals("fpname"))) ierr = i+k+1; k++;
         if (ierr != 0)
            irc[0] = 1000+ierr;
      }
      break;

      case 3:

// write namelist from packed data
      if (isign > 0) {
// pack data and names
         ddata[i] = idrun; cn[i] = "idrun"; i++;
         ddata[i] = indx; cn[i] = "indx"; i++;
         ddata[i] = nta; cn[i] = "nta"; i++;
         ddata[i] = modesxa; cn[i] = "modesxa"; i++;
         ddata[i] = psolve; cn[i] = "psolve"; i++;
         ddata[i] = ndim; cn[i] = "ndim"; i++;
         ddata[i] = omx; cn[i] = "omx"; i++;
         ddata[i] = omy; cn[i] = "omy"; i++;
         ddata[i] = omz; cn[i] = "omz"; i++;
         ddata[i] = ci; cn[i] = "ci"; i++;
         ddata[i] = narec[0]; cn[i] = "narec"; i++;
         ddata[i] = t0; cn[i] = "t0"; i++;
         ddata[i] = tend; cn[i] = "tend"; i++;
         ddata[i] = dt; cn[i] = "dt"; i++;
         ddata[i] = ceng; cn[i] = "ceng"; i++;
         sdata[k] = faname; cn[i+k] = "faname"; k++;
// set number of variables found
         ml[0] = i; ms[0] = k;
// write Fortran namelist from packed data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
      }
// read namelist from data
      else if (isign < 0) {
// read Fortran namelist and pack data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
// unpack data and check names
         idrun = (int) ddata[i];
            if (!(cn[i].equals("idrun"))) ierr = i+1; i++;
         indx = (int) ddata[i]; 
            if (!(cn[i].equals("indx"))) ierr = i+1; i++;
         nta = (int) ddata[i];
            if (!(cn[i].equals("nta"))) ierr = i+1; i++;
         modesxa = (int) ddata[i];
            if (!(cn[i].equals("modesxa"))) ierr = i+1; i++;
         psolve = (int) ddata[i];
            if (!(cn[i].equals("psolve"))) ierr = i+1; i++;
         ndim = (int) ddata[i];
            if (!(cn[i].equals("ndim"))) ierr = i+1; i++;
         omx = ddata[i];
            if (!(cn[i].equals("omx"))) ierr = i+1; i++;
         omy = ddata[i];
            if (!(cn[i].equals("omy"))) ierr = i+1; i++;
         omz = ddata[i];
            if (!(cn[i].equals("omz"))) ierr = i+1; i++;
         ci = ddata[i];
            if (!(cn[i].equals("ci"))) ierr = i+1; i++;
         narec[0] = (int) ddata[i];
            if (!(cn[i].equals("narec"))) ierr = i+1; i++;
         t0 = ddata[i];
            if (!(cn[i].equals("t0"))) ierr = i+1; i++;
         tend = ddata[i];
            if (!(cn[i].equals("tend"))) ierr = i+1; i++;
         dt = ddata[i];
            if (!(cn[i].equals("dt"))) ierr = i+1; i++;
         ceng = ddata[i];
            if (!(cn[i].equals("ceng"))) ierr = i+1; i++;
         faname = sdata[k];
            if (!(cn[i+k].equals("faname"))) ierr = i+k+1; k++;
         if (ierr != 0)
            irc[0] = 1000+ierr;
      }
      break;

      case 4:

// write namelist from packed data
      if (isign > 0) {
// pack data and names
         ddata[i] = idrun; cn[i] = "idrun"; i++;
         ddata[i] = indx; cn[i] = "indx"; i++;
         ddata[i] = nte; cn[i] = "nte"; i++;
         ddata[i] = modesxe; cn[i] = "modesxe"; i++;
         ddata[i] = psolve; cn[i] = "psolve"; i++;
         ddata[i] = ndim; cn[i] = "ndim"; i++;
         ddata[i] = omx; cn[i] = "omx"; i++;
         ddata[i] = omy; cn[i] = "omy"; i++;
         ddata[i] = omz; cn[i] = "omz"; i++;
         ddata[i] = ci; cn[i] = "ci"; i++;
         ddata[i] = nerec[0]; cn[i] = "nerec"; i++;
         ddata[i] = t0; cn[i] = "t0"; i++;
         ddata[i] = tend; cn[i] = "tend"; i++;
         ddata[i] = dt; cn[i] = "dt"; i++;
         ddata[i] = ceng; cn[i] = "ceng"; i++;
         sdata[k] = fename; cn[i+k] = "fename"; k++;
// set number of variables found
         ml[0] = i; ms[0] = k;
// write Fortran namelist from packed data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
      }
// read namelist from data
      else if (isign < 0) {
// read Fortran namelist and pack data
         frwnml1(fname,name,cn,ddata,sdata,iunit,isign,ml,ms,irc);
         if (irc[0] != 0)
            return;
// unpack data and check names
         idrun = (int) ddata[i];
            if (!(cn[i].equals("idrun"))) ierr = i+1; i++;
         indx = (int) ddata[i]; 
            if (!(cn[i].equals("indx"))) ierr = i+1; i++;
         nte = (int) ddata[i];
            if (!(cn[i].equals("nte"))) ierr = i+1; i++;
         modesxe = (int) ddata[i];
            if (!(cn[i].equals("modesxe"))) ierr = i+1; i++;
         psolve = (int) ddata[i];
            if (!(cn[i].equals("psolve"))) ierr = i+1; i++;
         ndim = (int) ddata[i];
            if (!(cn[i].equals("ndim"))) ierr = i+1; i++;
         omx = ddata[i];
            if (!(cn[i].equals("omx"))) ierr = i+1; i++;
         omy = ddata[i];
            if (!(cn[i].equals("omy"))) ierr = i+1; i++;
         omz = ddata[i];
            if (!(cn[i].equals("omz"))) ierr = i+1; i++;
         ci = ddata[i];
            if (!(cn[i].equals("ci"))) ierr = i+1; i++;
         nerec[0] = (int) ddata[i];
            if (!(cn[i].equals("nerec"))) ierr = i+1; i++;
         t0 = ddata[i];
            if (!(cn[i].equals("t0"))) ierr = i+1; i++;
         tend = ddata[i];
            if (!(cn[i].equals("tend"))) ierr = i+1; i++;
         dt = ddata[i];
            if (!(cn[i].equals("dt"))) ierr = i+1; i++;
         ceng = ddata[i];
            if (!(cn[i].equals("ceng"))) ierr = i+1; i++;
         fename = sdata[k];
            if (!(cn[i+k].equals("fename"))) ierr = i+k+1; k++;
         if (ierr != 0)
            irc[0] = 1000+ierr;
      }
      break;

      default:
         irc[0] = -18;
      }
      return;
   }
}

