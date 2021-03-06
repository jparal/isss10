!-----------------------------------------------------------------------
! * * * periodic 1d electromagnetic particle simulation kernel code * *
! this is a simple 1d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves only electrons, with periodic
! electromagnetic forces obtained by solving maxwell's equation with
! fast fourier transforms.
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G3
! copyright 1999, regents of the university of california
! update: july 14, 2011
      program bbeps1
      use init1d
      use bpush1d
      use dpush1d
      use push1d
      use fft1d
      use field1d
      use diag1d
      implicit none
! idimp = dimension of phase space = 4
! nmv = number of segments in v for velocity distribution
      integer :: idimp = 4, nmv = 40, ipbc = 1
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iup = 11, iua = 15
      integer :: iue = 26
      integer :: np, npx1, nx, nxh, nxe, nxeh, nx1
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: it, itw
      real :: zero = 0.0, time = 0.0, tloop = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tfft = 0.0, totpush = 0.0, ts = 0.0
      real :: qbme, affp, dth, qi0, omt, q2m0, wp0
      real :: we, wf, wm, wef, wke, ws
      real, dimension(:,:), pointer :: part
      real, dimension(:), pointer :: qe, fxe
      real, dimension(:,:), pointer :: cu, amu, fxyze, byze
      complex, dimension(:,:), pointer :: eyz, byz
      complex, dimension(:), pointer :: ffc, ffe
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:), pointer :: sfield
      complex, dimension(:), pointer :: pott
      real, dimension(:,:), pointer :: vfield
      complex, dimension(:,:), pointer :: vpott, vpotr
      real, dimension(:,:), pointer :: fv, fvm
      real, dimension(:,:), pointer :: wt
      character(len=10) :: cdrun
      character(len=32) :: fname
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
  993 format (' electric(l,t), magnetic energies = ',3e14.7)
! read namelist
      iuin = get_funit(iuin)
      open(unit=iuin,file='input1',form='formatted',status='old')
      read (iuin,input1)
! override input data
      idcode = 2
      psolve = 1
      ndim = 2
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file
      iuot = get_funit(iuot)
      fname = 'output1.'//cdrun
      open(unit=iuot,file=trim(fname),form='formatted',status='replace')
! np = total number of electrons in simulation
      np = npx + npxb; npx1 = npx + 1
      nx = 2**indx; nxh = nx/2
      nxe = nx + 4
!     ax = .866025
      if (inorder==LINEAR) then
!        ax = .912871
         nxe = nx + 2
      endif
      nxeh = nxe/2
! dimension for index and sorting arrays
      nx1 = nx + 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! part(1,n) = position x of particle n
! part(2,n) = velocity vx of particle n
! part(3,n) = velocity vy of particle n
! part(4,n) = velocity vz of particle n
      allocate(part(idimp,np))
! in real space, qe(j) = charge density at grid point j
! in real space, fxe(j) = longitudinal force/charge at grid point j,
! that is, fxe is the convolution of the longitudinal electric field
! over the particle shape
      allocate(qe(nxe),fxe(nxe))
! in real space, fxyze(i,j) = i component of force/charge at grid (j)
! that is, fxyze are the convolutions of the electric field
! over the particle shape
      allocate(fxyze(3,nxe))
! cu(i,j) = i component of current at grid (j).
! byze(i,j) = i component of magnetic field at grid (j).  
! byze is the convolution of the magnetic field over the particle shape
      allocate(cu(2,nxe),byze(2,nxe))
! in fourier space, exyz = transverse electric field
! in fourier space, bxyz = magnetic field
      allocate(eyz(2,nxeh),byz(2,nxeh))
! ffc = form factor array for poisson solver
      allocate(ffc(nxh))
! mixup = array of bit reversed addresses for fft
! sct = sine/cosine table for fft
      allocate(mixup(nxh),sct(nxh))
!
! open graphics device
      irc = open_graphs(nplot)
! initialize timer
      call wtimer(time,ltime,-1)
! initialize constants
      itime0 = 0
      itime = itime0
      ntime = itime + itime0
      qbme = qme
      affp = real(nx)/real(np)
      omt = sqrt(omy*omy + omz*omz)
      q2m0 = qbme*qme*real(np)/real(nx)
      dth = 0.0
! debug
!     dth = .5*dt
! end debug
      wp0 = q2m0*affp
! set initial time
      t0 = dt*real(itime0)
! set default diagnostic file names
      if (ntp > 0) fpname = 'potk1.'//cdrun
      if (nta > 0) faname = 'vpotk1.'//cdrun
      if (nte > 0) fename = 'vpotrk1.'//cdrun
! energy diagnostics
      if (ntw > 0) then
         allocate(wt((nloop-1)/ntw-(itime0/ntw)+1,7))
         itw = 0
      endif
! initialize electromagnetic fields
      byze = 0.0
      if (omt > 0.0) then
         call baddext(byze,omy,omz,nx,inorder)
         call cguard(byze,nx,inorder)
      endif
      byz = cmplx(0.0,0.0)
      eyz = cmplx(0.0,0.0)
      cu = 0.0
! prepare fft tables
      call fft_init(mixup,sct,indx)
! calculate form factors
      call pois_init(ffc,ax,affp,nx)
! initialize density profile and velocity distribution
! background electrons
      if (npx > 0) call distr(part,1,npx,vtx,vty,vtz,vx0,vy0,vz0,npx,nx,
     &ipbc)
! beam electrons
      if (npxb > 0) call distr(part,npx1,npxb,vtdx,vtdy,vtdz,vdx,vdy,vdz&
     &,npxb,nx,ipbc)
! initialize charge density to background
      qi0 = -qme/affp
! fix guiding centers for electrons
      if (omt > 0.0) call distr(part,byze,np,qbme,nx,ipbc,inorder)
! retard electron positions to deposit current
!     call retard(part,np,dth,nx)
! sorting arrays
      if (sortime > 0) allocate(pt(np),ip(np),npic(nx1))
! initialize diagnostics
! diagnostic metafile
      iudm = get_funit(iudm)
! velocity diagnostic
      if (ntv > 0) then
         allocate(fv(2*nmv+2,3),fvm(3,3))
         fv(1,:) = 8.*max(vtx,vty,vtz)
      endif
! potential diagnostic
      if (ntp > 0) then
         allocate(sfield(nxe))
         if (modesxp > nxh) modesxp = nxh
         allocate(pott(modesxp))
! open output file
         if (nprec==0) then
            nprec = -1; iup = get_funit(iup)
            call bfopen(pott,modesxp,iup,nprec,trim(fpname))
         endif
      endif
! vector potential or electromagnetic diagnostic
      if ((nta > 0) .or. (nte > 0)) then
         allocate(vfield(2,nxe))
         vfield = 0.0
      endif
! vector potential diagnostic
      if (nta > 0) then
         if (modesxa > nxh) modesxa = nxh
         allocate(vpott(2,modesxa))
! open output file
         if (narec==0) then
            narec = -1; iua = get_funit(iua)
            call bfopen(vpott,modesxa,iua,narec,trim(faname))
         endif
      endif
! electromagnetic diagnostic
      if (nte > 0) then
         if (modesxe > nxh) modesxe = nxh
         allocate(vpotr(2,modesxe))
! open output file
         if (nerec==0) then
            nerec = -1; iue = get_funit(iue)
            call bfopen(vpotr,modesxe,iue,nerec,trim(fename))
         endif
      endif
! record time
      call wtimer(time,ltime)
      write (iuot,*) 'initialization wall clock time = ', time, 'sec'
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
      write (iuot,991) ntime
!
! prepare electromagnetic diagnostic
      if (nte > 0) then
         it = ntime/nte
         if (ntime==nte*it) vfield = cu
      endif
! initialize current density to background
      call sguard(cu,zero,zero,nx,inorder)
! deposit electron current
      call djpost(part,cu,np,qme,dth,tdjpost,nx,ipbc,inorder,djopt)
! initialize charge density to background
      call sguard(qe,qi0,nx,inorder)
! deposit electron charge
      call dpost(part,qe,np,qme,tdpost,inorder,dopt)
! add guard cells for current
      call aguard(cu,nx,inorder)
! add guard cells
      call aguard(qe,nx,inorder)
! velocity diagnostic
      if (ntv > 0) then
         it = ntime/ntv
         if (ntime==ntv*it) then
! calculate paticle distribution function and moments
            call vdist(part,fv,fvm,np,nmv)
! display velocity distributions
            call displayfv(fv,fvm,' ELECTRON',ntime,nmv,2,irc)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,isign,mixup,sct,tfft,indx,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! calculate potential in fourier space
            isign = 1
            call pois(qe,sfield,isign,ffc,ws,nx,inorder)
! store selected fourier modes
            call gtmodes(sfield,pott,nx,modesxp,inorder)
! write diagnostic output
            call writebf(pott,modesxp,iup,nprec,order=LINEAR)
! transform potential to real space
            call fft(sfield,isign,mixup,sct,tfft,indx,inorder)
            call cguard(sfield,nx,inorder)
! display potential
            call displays(sfield,' POTENTIAL',ntime,999,0,nx,irc,inorder&
     &)
            if (irc==1) go to 2000
         endif
      endif
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,inorder)
! electromagnetic diagnostic
      if (nte > 0) then
         it = ntime/nte
         if (ntime==nte*it) then
! calculate averaged radiative vector potential
            vfield = 0.5*(vfield + cu)
            call avrpot(vfield,byz,ffc,ci,nx,inorder)
! store selected fourier modes
            call gtmodes(vfield,vpotr,nx,modesxe,inorder)
! write diagnostic output
            call writebf(vpotr,modesxe,iue,nerec,order=LINEAR)
! transform radiative vector potential to real space
            isign = 1
            call fft(vfield,isign,mixup,sct,tfft,indx,inorder)
            call cguard(vfield,nx,inorder)
! display radiative vector potential
            call displayv(vfield,' RADIATIVE VPOTENTIAL',ntime,999,0,1, &
     &nx,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! calculate electromagnetic fields in fourier space
      if (ntime==0) then
! calculate initial darwin magnetic field
         call ibpois(cu,byz,ffc,ci,wm,nx,inorder)
         wf = 0.
! calculate initial darwin electric field
         allocate(amu(2,nxe),ffe(nxeh))
! deposit momentum flux
         call sguard(amu,zero,zero,nx,inorder)
         call dmjpost(part,amu,np,qme,ts,inorder,djopt)
         call aguard(amu,nx,inorder)
! solve for darwin electric field
         isign = -1
         call fft(amu,isign,mixup,sct,tfft,indx,inorder)
!        byze = cu
         call dcuperp(byze,amu,nx,inorder)
         call epois_init(ffe,ax,affp,wp0,ci,nx)
         call iepois(byze,eyz,ffe,ci,wf,nx,inorder)
         deallocate(amu,ffe)
         dth = .5*dt
! calculate electromagnetic fields
      else
         call maxwel(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,inorder)
      endif
! vector potential diagnostic
      if (nta > 0) then
         it = ntime/nta
         if (ntime==nta*it) then
! calculate vector potential in fourier space
            call avpot(byz,vfield,nx,inorder)
! store selected fourier modes
            call gtmodes(vfield,vpott,nx,modesxa,inorder)
! write diagnostic output
            call writebf(vpott,modesxa,iua,narec,order=LINEAR)
! transform vector potential to real space
            isign = 1
            call fft(vfield,isign,mixup,sct,tfft,indx,inorder)
            call cguard(vfield,nx,inorder)
! display vector potential
            call displayv(vfield,' VECTOR POTENTIAL',ntime,999,0,1,nx,  &
     &irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! calculate longitudinal electric field in fourier space
      isign = -1
      call pois(qe,fxe,isign,ffc,we,nx,inorder)
! add longitudinal and transverse electric fields
      call emfield(fxyze,fxe,eyz,ffc,nx,inorder)
! copy magnetic field
      call emfield(byze,byz,ffc,nx,inorder)
! transform force/charge to real space
      isign = 1
      call fft(fxyze,isign,mixup,sct,tfft,indx,inorder)
      call cguard(fxyze,nx,inorder)
! transform magnetic field to real space
      isign = 1
      call fft(byze,isign,mixup,sct,tfft,indx,inorder)
! add external magnetic field
      if (omt > 0.0) call baddext(byze,omy,omz,nx,inorder)
      call cguard(byze,nx,inorder)
! push particles
      wke = 0.
      call push3(part,fxyze,byze,omx,np,qbme,dt,dth,wke,tpush,nx,ipbc,  &
     &inorder,popt)
! sort electrons
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
            call sortp(part,pt,ip,np,npic,tsort,inorder)
         endif
      endif
! energy diagnostic
      if (ntw > 0) then
         it = itime/ntw
         if (itime==ntw*it) then
            wef = we + wf + wm
            ws = wef + wke
            write (iuot,992) wef, wke, ws
            write (iuot,993) we, wf, wm
            itw = itw + 1
            wt(itw,:) = (/wef,wke,zero,ws,we,wf,wm/)
         endif
      endif
      itime = itime + 1
      ntime = itime + itime0
      call wtimer(tloop,ltime)
      time = time + tloop
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
! energy diagnostic
      if (ntw > 0) then
         ts = t0 + dt*real(ntw)*(itime0-(itime0/ntw)*ntw)
         call reset_graphs
         call displayw(wt,ts,dt*real(ntw),itw,irc)
      endif
! accumulate timings
      write (iuot,*) 'electromagnetic code bbeps1'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost + tdjpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'electron current deposit time = ', tdjpost, 'sec'
      write (iuot,*) 'total electron push time = ', totpush, 'sec'
      write (iuot,*) 'electron sort time = ', tsort
      totpush = totpush + tsort
      write (iuot,*) 'total electron time = ', totpush, 'sec'
      write (iuot,*) 'total fft time=', tfft, 'sec'
      time = time - (totpush + tfft)
      write (iuot,*) 'other time=', time, 'sec'
! write final diagnostic metafile
      fname = 'diag1.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! potential diagnostics
      if (ntp > 0) then
         nprec = nprec - 1
         ceng = affp
         write (iudm,pot1d)
      endif
! vector potential diagnostics
      if (nta > 0) then
         narec = narec - 1
         ceng = affp
         write (iudm,vpot1d)
      endif
! electromagnetic diagnostics
      if (nte > 0) then
         nerec = nerec - 1
         ceng = affp
         write (iudm,em1d)
      endif
! write out input file
      write (iudm,input1)
      write (iuot,*) ' * * * q.e.d. * * *'
      close(unit=iudm)
      close(unit=iuot)
! close graphics device
      call close_graphs
      stop
      end program bbeps1
