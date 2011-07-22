//----------------------------------------------------------------------
// * * * periodic 1d electrostatic particle simulation kernel code * * *
// this is a simple 1d skeleton particle-in-cell code designed for
// exploring new computer architectures.  it contains the critical pieces
// needed for depositing charge, advancing particles, and solving the
// field.  the code moves only electrons, with periodic electrostatic
// forces obtained by solving poisson's equation with fast fourier
// transforms. 
// written by viktor k. decyk, ucla
// Java wrapper to Fortran90 code
// copyright 2011, regents of the university of california
// update: july 7, 2011

import java.io.*;

import static simul1d.globals.*;
import static simul1d.input1d.*;
import static simul1d.init1d.*;
import static simul1d.push1d.*;
import static simul1d.fft1d.*;
import static simul1d.field1d.*;
import static simul1d.diag1d.*;

public class beps1 {

   static {
      System.loadLibrary("CJavaInterface");
   }

   public static void main(String[] args) {
// idimp = dimension of phase space = 2
// nmv = number of segments in v for velocity distribution
      int idimp = 2, nmv = 40, ipbc = 1;
// default unit numbers
      int iuin = 8, iuot = 18, iudm = 19, iup = 11;
      int np, npx1, nx, nxh, nxe, nx1;
      int nloop, itime, itime0, ntime, isign, it;
      int itw = 0, ltw = 0;
      double zero = 0.0, time = 0.0, ts = 0.0;
      int[] jup = new int[1]; int[] irc = new int[1];
      long[] ltime = new long[1];
      double[] tloop = new double[1];
      double[] tpush = new double[1]; double[] tdpost = new double[1];
      double[] tsort = new double[1]; double[] tfft = new double[1];
      double[] we = new double[1]; double[] wke = new double[1];
      double[] ws = new double[1];
      double totpush = 0.0;
      double qbme, affp, qi0;
      double[] part, qe, fxe, pt, sfield, fv, fvm, wt;
      double[] ffc, sct, pott;
      int[] mixup, ip, npic;
      String cdrun, fname, name;
      PrintWriter jot = null;
// read namelist
      iuin = fgetunit(iuin);
      fname = "input1"; name = "input1";
      jrwnml1(fname,name,iuin,-1,irc);
      if (irc[0] != 0)
         System.out.println("namelist input error irc="+irc[0]);
// override input data
      idcode = 1;
      psolve = 1;
// create string from idrun
      cdrun = "" + idrun;
// text output file
      iuot = fgetunit(iuot);
      fname = "output1." + cdrun;
      try {
         jot = new PrintWriter(fname);
      } catch (IOException e) {
         System.err.println(e);
         return;
      }
// np = total number of electrons in simulation
      np = npx + npxb; npx1 = npx + 1;
      nx = 1<<indx; nxh = nx/2;
      nxe = nx + 4;
//    ax = .866025;
      if (inorder==LINEAR) {
//       ax = .912871;
         nxe = nx + 2;
      }
// dimension for index and sorting arrays
      nx1 = nx + 1;
// nloop = number of time steps in simulation
      nloop = (int) (tend/dt + .0001);
// part(1,n) = position x of particle n
// part(2,n) = velocity vx of particle n
      part = new double[idimp*np];
// in real space, qe(j) = charge density at grid point j
// in real space, fxe(j) = longitudinal force/charge at grid point j,
// that is, fxe is the convolution of the longitudinal electric field
// over the particle shape
      qe = new double[nxe]; fxe = new double[nxe];
// ffc = form factor array for poisson solver
      ffc = new double[2*nxh];
// mixup = array of bit reversed addresses for fft
// sct = sine/cosine table for fft
      mixup = new int[nxh]; sct = new double[2*nxh];
//
// open graphics device
//    irc = open_graphs(nplot)
// initialize timer
      jtimer(tloop,ltime,-1);
// initialize constants
      itime0 = 0;
      itime = itime0;
      ntime = itime + itime0;
      qbme = qme;
      affp = (double) nx/(double) np;
// set initial time
      t0 = dt*(double) itime0;
// set default diagnostic file names
      if (ntp > 0)
         fpname = "potk1." + cdrun;
// energy time history
      if (ntw > 0) {
         ltw = (nloop - 1)/ntw - (itime0/ntw) + 1;
         wt = new double[ltw*4];
      }
      else {
         wt = new double[0];
      }
// prepare fft tables
      fft_init(mixup,sct,indx);
// calculate form factors
      pois_init(ffc,ax,affp,nx);
// initialize density profile and velocity distribution
// background electrons
      if (npx > 0)
         distr(part,1,npx,vtx,vx0,npx,nx,ipbc,idimp);
// beam electrons
      if (npxb > 0)
         distr(part,npx1,npxb,vtdx,vdx,npxb,nx,ipbc,idimp);
// initialize charge density to background
      qi0 = -qme/affp;
// sorting arrays
      if (sortime > 0) {
         pt = new double[np]; ip = new int[np];
         npic = new int[nx1];
      }
      else {
         pt = new double[0]; ip = new int[0]; npic = new int[0];
      }
// initialize diagnostics
// diagnostic metafile
      iudm = fgetunit(iudm);
// velocity diagnostic
      if (ntv > 0) {
         fv = new double[2*nmv+2]; fvm = new double[3];
         fv[0] = 8.*vtx;
      }
      else {
         fv = new double[0]; fvm = new double[0];
      }
// potential diagnostic
      if (ntp > 0) {
         sfield = new double[nxe];
         if (modesxp > nxh)
            modesxp = nxh;
         pott = new double[2*modesxp];
// open output file
         if (nprec[0]==0) {
            nprec[0] = -1; iup = fgetunit(iup); jup[0] = iup;
            fbfcopen(pott,modesxp,iup,nprec,fpname);
         }
      }
      else {
         sfield = new double[0]; pott = new double[0];
      }
// record time
      jtimer(tloop,ltime,1);
      time = tloop[0];
      jot.printf("initialization wall clock time = %f sec%n",time);
//
// * * * start main iteration loop * * *
//
      while (nloop > ntime) {
         jot.printf(" T = %d%n",ntime);
// initialize charge density to background
         sguard(qe,qi0,nx,inorder);
// deposit electron charge
         dpost(part,qe,np,qme,tdpost,idimp,inorder,dopt);
// add guard cells
         aguard(qe,nx,inorder);
// velocity diagnostic
         if (ntv > 0) {
            it = ntime/ntv;
            if (ntime==ntv*it) {
// calculate paticle distribution function and moments
               vdist(part,fv,fvm,np,nmv,idimp,1);
// display velocity distributions
//             call displayfv(fv,fvm,' ELECTRON',ntime,nmv,2,irc)
//             if (irc==1) go to 2000
            }
         }
// transform charge to fourier space
         isign = -1;
         fft1(qe,isign,mixup,sct,tfft,indx,inorder);
// potential diagnostic
         if (ntp > 0) {
            it = ntime/ntp;
            if (ntime==ntp*it) {
// calculate potential in fourier space
               isign = 1;
               pois(qe,sfield,isign,ffc,ws,nx,inorder);
// store selected fourier modes
               gtmodes(sfield,pott,nx,modesxp,inorder);
// write diagnostic output
               writebfc(pott,modesxp,jup,nprec,fpname,LINEAR);
// transform potential to real space
               fft1(sfield,isign,mixup,sct,tfft,indx,inorder);
               dguard(sfield,nx,inorder);
// display potential
//             call displays(sfield,' POTENTIAL',ntime,999,0,nx,irc,    &
//   &inorder)
//             if (irc==1) go to 2000
            }
         }
// calculate force/charge in fourier space
         isign = -1;
         pois(qe,fxe,isign,ffc,we,nx,inorder);
// transform force/charge to real space
         isign = 1;
         fft1(fxe,isign,mixup,sct,tfft,indx,inorder);
         cguard(fxe,nx,1,inorder);
// push particles
         wke[0] = 0.0;
         push(part,fxe,np,qbme,dt,wke,tpush,nx,ipbc,idimp,inorder,popt);
// sort electrons
         if (sortime > 0) {
            if (ntime%sortime==0) {
               sortp(part,pt,ip,np,npic,tsort,idimp,inorder);
            }
         }
// energy diagnostic
         if (ntw > 0) {
            it = itime/ntw;
            if (itime==ntw*it) {
               ws[0] = we[0] + wke[0];
               jot.printf(" field, kinetic, total energies = " +
                          "%e %e %e%n",we[0],wke[0],ws[0]);
               wt[itw] = we[0]; wt[itw+ltw] = wke[0];
               wt[itw+2*ltw] = zero; wt[itw+3*ltw] = ws[0]; 
               itw += 1;
            }
         }
         itime += 1;
         ntime = itime + itime0;
         jtimer(tloop,ltime,1);
         time += tloop[0];
      }
//
// * * * end main iteration loop * * *
//
// energy diagnostic
      if (ntw > 0) {
         ts = t0 + dt*(double) (ntw*(itime0-(itime0/ntw)*ntw));
//       call reset_graphs
//       call displayw(wt,ts,dt*real(ntw),itw,irc)
      }
// accumulate timings
      jot.printf("electrostatic Java code beps1%n");
      jot.printf("main wall clock time = %f sec%n",time);
      totpush = tpush[0] + tdpost[0];
      jot.printf("electron push time = %f sec%n",tpush[0]);
      jot.printf("electron charge deposit time = %f sec%n",tdpost[0]);
      jot.printf("total electron push time = %f sec%n",totpush);
      jot.printf("electron sort time = %f sec%n",tsort[0]);
      totpush += tsort[0];
      jot.printf("total electron time = %f sec%n",totpush);
      jot.printf("total fft time = %f sec%n",tfft[0]);
      time -= (totpush + tfft[0]);
      jot.printf("other time = %f sec%n",time);
// write final diagnostic metafile
      fname = "diag1." + cdrun;
// potential diagnostics
      if (ntp > 0) {
         nprec[0] -= 1;
         ceng = affp;
         name = "pot1d";
         jrwnml1(fname,name,iudm,1,irc);
         if (irc[0] != 0)
            System.out.println(name+" output error irc="+irc[0]);
         frclose(iup);
      }
// write out input file
      name = "input1";
      jrwnml1(fname,name,iudm,1,irc);
      if (irc[0] != 0)
         System.out.println(name+" output error irc="+irc[0]);
      jot.printf(" * * * q.e.d. * * *%n");
// debug
      System.out.println(" * * * q.e.d. * * *");
//
      frclose(iudm);
      jot.close();
// close graphics device
//    call close_graphs
   }
}