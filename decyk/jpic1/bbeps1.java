//----------------------------------------------------------------------
// * * * periodic 1d electromagnetic particle simulation kernel code * *
// this is a simple 1d skeleton particle-in-cell code designed for
// exploring new computer architectures.  it contains the critical pieces
// needed for depositing charge and current, advancing particles, and
// solving the fields.  the code moves only electrons, with periodic
// electromagnetic forces obtained by solving maxwell's equation with
// fast fourier transforms.
// written by viktor k. decyk, ucla
// Java wrapper to Fortran90 code
// copyright 2011, regents of the university of california
// update: july 7, 2011

import java.io.*;

import static simul1d.globals.*;
import static simul1d.input1d.*;
import static simul1d.init1d.*;
import static simul1d.bpush1d.*;
import static simul1d.dpush1d.*;
import static simul1d.push1d.*;
import static simul1d.fft1d.*;
import static simul1d.field1d.*;
import static simul1d.diag1d.*;

public class bbeps1 {

   static {
      System.loadLibrary("CJavaInterface");
   }

   public static void main(String[] args) {
// idimp = dimension of phase space = 4
// nmv = number of segments in v for velocity distribution
      int idimp = 4, nmv = 40, ipbc = 1;
// default unit numbers;
      int iuin = 8, iuot = 18, iudm = 19, iup = 11, iua = 15;
      int iue = 26;
      int j;
      int np, npx1, nx, nxh, nxe, nxeh, nx1;
      int nloop, itime, itime0, ntime, isign, it;
      int itw = 0, ltw = 0;
      double zero = 0.0, time = 0.0, ts = 0.0;
      int[] jup = new int[1]; int[] jua = new int[1];
      int[] jue = new int[1]; int[] irc = new int[1];
      long[] ltime = new long[1];
      double[] tloop = new double[1]; double[] tpush = new double[1];
      double[] tdpost = new double[1]; double[] tdjpost = new double[1];
      double[] tsort = new double[1]; double[] tfft = new double[1];
      double[] we = new double[1]; double[] wf = new double[1];
      double[] wm = new double[1]; double[] wef = new double[1];
      double[] wke = new double[1]; double[] ws = new double[1];
      double totpush = 0.0;
      double qbme, affp, dth, qi0, omt, q2m0, wp0;
      double[] part, qe, fxe, cu, amu, fxyze, byze;
      double[] pt, sfield, vfield, fv, fvm, wt;
      double[] eyz, byz, ffc, ffe, sct, pott, vpott, vpotr;
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
      idcode = 2;
      psolve = 1;
      ndim = 2;
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
      nxeh = nxe/2;
// dimension for index and sorting arrays
      nx1 = nx + 1;
// nloop = number of time steps in simulation
      nloop = (int) (tend/dt + .0001);
// part(1,n) = position x of particle n
// part(2,n) = velocity vx of particle n
// part(3,n) = velocity vy of particle n
// part(4,n) = velocity vz of particle n
      part = new double[idimp*np];
// in real space, qe(j) = charge density at grid point j
// in real space, fxe(j) = longitudinal force/charge at grid point j,
// that is, fxe is the convolution of the longitudinal electric field
// over the particle shape
      qe = new double[nxe]; fxe = new double[nxe];
// in real space, fxyze(i,j) = i component of force/charge at grid (j)
// that is, fxyze are the convolutions of the electric field
// over the particle shape
      fxyze = new double[3*nxe];
// cu(i,j) = i component of current at grid (j).
// byze(i,j) = i component of magnetic field at grid (j).  
// byze is the convolution of the magnetic field over the particle shape
      cu = new double[2*nxe]; byze = new double[2*nxe];
// in fourier space, exyz = transverse electric field
// in fourier space, bxyz = magnetic field
      eyz = new double[4*nxeh]; byz = new double[4*nxeh];
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
      omt = Math.sqrt(omy*omy + omz*omz);
      q2m0 = qbme*qme*(double) np/(double) nx;
      dth = 0.0;
// debug
//    dth = .5*dt;
// end debug
      wp0 = q2m0*affp;
// set initial time
      t0 = dt*(double) itime0;
// set default diagnostic file names
      if (ntp > 0)
         fpname = "potk1." + cdrun;
      if (nta > 0)
         faname = "vpotk1." + cdrun;
      if (nte > 0)
         fename = "vpotrk1." + cdrun;
// energy diagnostics
      if (ntw > 0) {
         ltw = (nloop - 1)/ntw - (itime0/ntw) + 1;
         wt = new double[ltw*7];
      }
      else {
         wt = new double[0];
      }
// initialize electromagnetic fields
      if (omt > 0.0) {
         baddext(byze,omy,omz,nx,ndim,inorder);
         cguard(byze,nx,ndim,inorder);
      }
// prepare fft tables
      fft_init(mixup,sct,indx);
// calculate form factors
      pois_init(ffc,ax,affp,nx);
// initialize density profile and velocity distribution
// background electrons
      if (npx > 0)
         distrh(part,1,npx,vtx,vty,vtz,vx0,vy0,vz0,npx,nx,ipbc,idimp);
// beam electrons
      if (npxb > 0)
         distrh(part,npx1,npxb,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,nx,ipbc,
                idimp);
// initialize charge density to background
      qi0 = -qme/affp;
// fix guiding centers for electrons
      if (omt > 0.0)
         bdistr(part,byze,np,qbme,nx,ipbc,idimp,ndim,inorder);
// retard electron positions to deposit current
//     call retard(part,np,dth,nx)
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
         fv = new double[(2*nmv+2)*3]; fvm = new double[3*3];
         fv[0] = 8.*Math.max(vtx,Math.max(vty,vtz));
         fv[2*nmv+2] = fv[0]; fv[(2*nmv+2)*2] = fv[0];
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
// vector potential or electromagnetic diagnostic
      if ((nta > 0) || (nte > 0)) {
         vfield = new double[2*nxe];
      }
      else {
         vfield = new double[0];
      }
// vector potential diagnostic
      if (nta > 0) {
         if (modesxa > nxh)
            modesxa = nxh;
         vpott = new double[4*modesxa];
// open output file
         if (narec[0]==0) {
            narec[0] = -1; iua = fgetunit(iua); jua[0] = iua;
            fbfvcopen(vpott,modesxa,ndim,iua,narec,faname);
         }
      }
      else {
         vpott = new double[0];
      }
// electromagnetic diagnostic
      if (nte > 0) {
         if (modesxe > nxh)
            modesxe = nxh;
         vpotr = new double[4*modesxe];
// open output file
         if (nerec[0]==0) {
            nerec[0] = -1; iue = fgetunit(iue); jue[0] = iue;
            fbfvcopen(vpotr,modesxe,ndim,iue,nerec,fename);
         }
      }
      else {
         vpotr = new double[0];
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
//
// prepare electromagnetic diagnostic
         if (nte > 0) {
            it = ntime/nte;
            if (ntime==nte*it)
               System.arraycopy(cu,0,vfield,0,vfield.length);
         }
// initialize current density to background
         scguard(cu,zero,zero,nx,inorder);
// deposit electron current
         djpost(part,cu,np,qme,dth,tdjpost,nx,ipbc,idimp,ndim,inorder,
                djopt);
// initialize charge density to background
         sguard(qe,qi0,nx,inorder);
// deposit electron charge
         dpost(part,qe,np,qme,tdpost,idimp,inorder,dopt);
// add guard cells for current
         acguard(cu,nx,inorder);
// add guard cells
         aguard(qe,nx,inorder);
// velocity diagnostic
         if (ntv > 0) {
            it = ntime/ntv;
            if (ntime==ntv*it) {
// calculate paticle distribution function and moments
               vdist(part,fv,fvm,np,nmv,idimp,3);
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
               writebfc(pott,modesxp,jup,nprec,fpname,inorder);
// transform potential to real space
               fft1(sfield,isign,mixup,sct,tfft,indx,inorder);
               dguard(sfield,nx,inorder);
// display potential
//             call displays(sfield,' POTENTIAL',ntime,999,0,nx,irc,    &
//   &inorder)
//             if (irc==1) go to 2000
            }
         }
// transform current to fourier space
         isign = -1;
         fft2(cu,isign,mixup,sct,tfft,indx,ndim,inorder);
// electromagnetic diagnostic
         if (nte > 0) {
            it = ntime/nte;
            if (ntime==nte*it) {
// calculate averaged radiative vector potential
               for (j = 0; j < vfield.length; j++) {
                  vfield[j] = 0.5*(vfield[j] + cu[j]);
               }
               avrpot(vfield,byz,ffc,ci,nx,inorder);
// store selected fourier modes
               gtvmodes(vfield,vpotr,nx,modesxe,ndim,inorder);
// write diagnostic output
               writebfvc(vpotr,modesxe,jue,nerec,fename,ndim,LINEAR);
// transform radiative vector potential to real space
               isign = 1;
               fft2(vfield,isign,mixup,sct,tfft,indx,ndim,inorder);
               cguard(vfield,nx,ndim,inorder);
// display radiative vector potential
//             call displayv(vfield,' RADIATIVE VPOTENTIAL',ntime,999,0,&
//   &1,nx,irc,inorder)
//             if (irc==1) go to 2000
            }
         }
// calculate electromagnetic fields in fourier space
         if (ntime==0) {
// calculate initial darwin magnetic field
            ibpois(cu,byz,ffc,ci,wm,nx,inorder);
            wf[0] = 0.0;
// calculate initial darwin electric field
            amu = new double[2*nxe]; ffe = new double[2*nxeh];
// deposit momentum flux
            scguard(amu,zero,zero,nx,inorder);
            dmjpost(part,amu,np,qme,ws,idimp,ndim,inorder,djopt);
            acguard(amu,nx,inorder);
// solve for darwin electric field
            isign = -1;
            fft2(amu,isign,mixup,sct,tfft,indx,ndim,inorder);
            System.arraycopy(cu,0,byze,0,byze.length);
            dcuperp(byze,amu,nx,ndim,inorder);
            epois_init(ffe,ax,affp,wp0,ci,nx);
            iepois(byze,eyz,ffe,ci,wf,nx,ndim,inorder);
            dth = .5*dt;
         }
// calculate electromagnetic fields
         else {
            maxwel(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,inorder);
         }
// vector potential diagnostic
         if (nta > 0) {
            it = ntime/nta;
            if (ntime==nta*it) {
// calculate vector potential in fourier space
               avpot(byz,vfield,nx,inorder);
// store selected fourier modes
               gtvmodes(vfield,vpott,nx,modesxa,ndim,inorder);
// write diagnostic output
               writebfvc(vpott,modesxa,jua,narec,faname,ndim,LINEAR);
// transform vector potential to real space
               isign = 1;
               fft2(vfield,isign,mixup,sct,tfft,indx,ndim,inorder);
               cguard(vfield,nx,ndim,inorder);
// display vector potential
//             call displayv(vfield,' VECTOR POTENTIAL',ntime,999,0,1,  &
//   &nx,irc,inorder)
//             if (irc==1) go to 2000
            }
         }
// calculate longitudinal electric field in fourier space
         isign = -1;
         pois(qe,fxe,isign,ffc,we,nx,inorder);
// add longitudinal and transverse electric fields
         emfield(fxyze,fxe,eyz,ffc,nx,inorder);
// copy magnetic field
         bmfield(byze,byz,ffc,nx,inorder);
// transform force/charge to real space
         isign = 1;
         fft2(fxyze,isign,mixup,sct,tfft,indx,3,inorder);
         cguard(fxyze,nx,3,inorder);
// transform magnetic field to real space
         isign = 1;
         fft2(byze,isign,mixup,sct,tfft,indx,ndim,inorder);
// add external magnetic field
         if (omt > 0.0)
            baddext(byze,omy,omz,nx,ndim,inorder);
         cguard(byze,nx,ndim,inorder);
// push particles
         wke[0] = 0.0;
         push3(part,fxyze,byze,omx,np,qbme,dt,dth,wke,tpush,nx,ipbc,
               idimp,ndim,inorder,popt);
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
               wef[0] = we[0] + wf[0] + wm[0];
               ws[0] = wef[0] + wke[0];
               jot.printf(" field, kinetic, total energies = " +
                          "%e %e %e%n",wef[0],wke[0],ws[0]);
               jot.printf(" electric(l,t), magnetic energies = " +
                          "%e %e %e%n",we[0],wf[0],wm[0]);
               wt[itw] = wef[0]; wt[itw+ltw] = wke[0];
               wt[itw+2*ltw] = zero; wt[itw+3*ltw] = ws[0];
               wt[itw+4*ltw] = we[0]; wt[itw+5*ltw] = wf[0]; 
               wt[itw+6*ltw] = wm[0]; 
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
      jot.printf("electromagnetic Java code bbeps1%n");
      jot.printf("main wall clock time = %f sec%n",time);
      totpush = tpush[0] + tdpost[0] + tdjpost[0];
      jot.printf("electron push time = %f sec%n",tpush[0]);
      jot.printf("electron charge deposit time = %f sec%n",tdpost[0]);
      jot.printf("electron current deposit time = %f sec%n",tdjpost[0]);
      jot.printf("total electron push time = %f sec%n",totpush);
      jot.printf("electron sort time = %f sec%n",tsort[0]);
      totpush += tsort[0];
      jot.printf("total electron time = %f sec%n", totpush);
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
// vector potential diagnostics
      if (nta > 0) {
         narec[0] -= 1;
         ceng = affp;
         name = "vpot1d";
         jrwnml1(fname,name,iudm,1,irc);
         if (irc[0] != 0)
            System.out.println(name+" output error irc="+irc[0]);
         frclose(iua);
      }
// electromagnetic diagnostics
      if (nte > 0) {
         nerec[0] -= 1;
         ceng = affp;
         name = "em1d";
         jrwnml1(fname,name,iudm,1,irc);
         if (irc[0] != 0)
            System.out.println(name+" output error irc="+irc[0]);
         frclose(iue);
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
