//----------------------------------------------------------------------
package simul1d;

// Java interface to 1d PIC Fortran77 library input1mod.f

import static simul1d.globals.*;

public class input1d {

// Namelist input1

// idrun/idrun0 = run identifier for current/old run
// idcode = code identifier
   public static int idrun = 1, idrun0 = 0, idcode = 0;
// indx = exponent which determines length in x direction, nx=2**indx
// npx = number of background particles distributed in x direction
// npxb = number of beam particles in x direction
// public static int indx =   7, npx =    5120, npxb = 512;;
   public static int indx =  11, npx =   409600, npxb = 40960;
// public static int indx =  14, npx =  4096000, npxb =   409600;
// public static int indx =  18, npx = 40960000, npxb =  4096000;
// inorder = interpolation order
// popt = particle optimization scheme
// dopt = charge deposit optimization scheme
// djopt = current deposit optimization scheme
   public static int inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD;
   public static int djopt = STANDARD;
// nustrt = (0,1,2) = this is an (old start,new start,restart)
// ntr = number of time steps between restart routine
   public static int nustrt = 1, ntr = 0;
// ntw, ndw = number of time steps between energy diagnostic
// ntp, ndp = number of time steps between potential diagnostic
// nta, nda = number of time steps between vector potential diagnostic
// ntv, ndv = number of time steps between velocity-space diagnostic
// nts, nds = number of time steps between phase space diagnostic
   public static int ntw = 1, ntp = 0, nta = 0, ntv = 0, nts = 0;
// public static int ndw = 1, ndp = 0, nda = 0, ndv = 0, nds = 0;
// nte, nde = number of time steps between electromagnetic diagnostic
   public static int nte = 0;
// public static int nde = 0;
// tend = time at end of simulation, in units of plasma frequency
// dt = time interval between successive calculations
   public static double tend =  45.000, dt = 0.2000000e+00;
// qme = charge on electron, in units of e
// vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
   public static double qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0;
// vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
   public static double vx0 = 0.0, vy0 = 0.0, vz0 = 0.0;
// vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
   public static double vdx = 0.0, vdy = 0.0, vdz = 0.0;
// vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
   public static double vtdx = 1.0, vtdy = 1.0, vtdz = 1.0;
// psolve = type of poisson solver = (1,2,3)
   public static int psolve = PERIODIC_2D;
// relativity = (no,yes) = (0,1) = relativity is used
// public static int relativity = 0;
// omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
   public static double omx = 0.0, omy = 0.0, omz = 0.0;
// ci = reciprical of velocity of light
   public static double ci = 0.1;
// ax = half-width of particle in x direction
// public static double ax = .816497;
// public static double ax = .866025;
   public static double ax = .912871;
// ndim = number of velocity dimensions = 1 or 2
   public static int ndim = 2;
// ndc = number of corrections in darwin iteration
   public static int ndc = 2;
// movion = (0,1) = (no,yes) move the ions
   public static int movion = 0;
// sortime = number of time steps between electron sorting
   public static int sortime = 100;
// nplot = maximum number of plots per page
   public static int nplot = 4;
// sntasks = (-1,n) = set maximum number of tasks (-1 = number of cpus-1)
   public static int sntasks = -1;
// ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
//                        hyperbolic secant squared=4)
   public static int ndprof = 0;
// ampdx = amplitude of density compared to uniform in x
// scaledx = scale length for spatial coordinate in x
// shiftdx = shift of spatial coordinate in x
   public static double ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0;
// modesxp = number of modes in x to keep for potential diagnostic
// modesxa = number of modes in x to keep for vector potential diagnostic
// modesxe = number of modes in x to keep for electromagnetic diagnostic
//           (radiative part of vector potential)
   public static int modesxp = 11, modesxa = 11, modesxe = 11;

// t0 = initial time value
// ceng = energy normalization
   public static double t0 = 0.0, ceng = 1.0;

// Namelist output for potential diagnostic
// nprec = current record number for potential writes
   public static int[] nprec = new int[1];
// fpname = file name for potential diagnostic
   public static String fpname = "potk1.0";

// Namelist output for vector potential diagnostic
// narec = current record number for vector potential writes
   public static int[] narec = new int[1];
// faname = file name for vector potential diagnostic
   public static String faname = "vpotk1.0";

// Namelist output for electromagnetic diagnostic
// nerec = current record number for electromagnetic writes
   public static int[] nerec = new int[1];
// fename = file name for electromagnetic diagnostic
   public static String fename = "vpotrk1.0";

}
