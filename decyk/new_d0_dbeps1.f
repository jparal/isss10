!-----------------------------------------------------------------------
! * * * periodic 1d darwin particle simulation kernel code * * *
! this is a simple 1d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves only electronss, with 
! electromagnetic forces obtained by solving non-radiative form of
! maxwell's equation with fast fourier transforms
! for various boundary conditions. 
! using algorithm similar to that described in
! J. Busnardo-Neto, P. L. Pritchett, A. T. Lin, and J. M. Dawson,
! J. Computational Phys. 23, 300 (1977).
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G5
! copyright 1999, regents of the university of california
! update: august 18, 2010
      program d0_dbeps1
      use init1d
      use bpush1d
      use dpush1d
      use push1d
      use fft1d
      use field1d
      use dfield1d
      use diag1d
!     use mp0d, only: mpinit, ncpus
      implicit none
! idimp = dimension of phase space = 4
! nmv = number of segments in v for velocity distribution
      integer :: idimp = 4, nmv = 40
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iup = 11, iua = 15
      integer :: np, npx1, nx, nxh, nxe, nxeh, nx1, ipbc
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: k, it, itw
      real :: zero = 0.0, time = 0.0, tloop = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tfft = 0.0, totpush = 0.0, ts = 0.0
      real :: tdcjpost = 0.0
      real :: qbme, affp, qi0, omt, q2m0, wp0, wpmax, wpmin
      real :: we, wf, wm, wef, wke, ws
      real, dimension(:,:), pointer :: part
      real, dimension(:), pointer :: qe, fxe
      real, dimension(:,:), pointer :: cu, cus, amu, fxyze, eyze, byze
      complex, dimension(:), pointer :: ffc, ffe
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:), pointer :: sfield
      complex, dimension(:), pointer :: pott
      real, dimension(:,:), pointer :: vfield
      complex, dimension(:,:), pointer :: vpott
      real, dimension(:,:), pointer :: fv, fvm
      real, dimension(:,:), pointer :: wt
! dirichlet boundary conditions
!     integer :: indx1, nxv, nx2, nx2v
!     real, dimension(:), pointer :: q2, fx2
!     real, dimension(:,:), pointer :: cu2, amu2, fxyz2, byz2
!     complex, dimension(:), pointer :: ffd
!     integer, dimension(:), pointer :: mixup2
!     complex, dimension(:), pointer :: sct2
!     real, dimension(:), pointer :: sfield2, pots
!     real, dimension(:,:), pointer :: vfield2, vpots, vpotsr
!
      character(len=10) :: cdrun
      character(len=32) :: fname
! debug
      integer :: i, j
      real :: eps, epsmax
! end debug
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
  993 format (' electric(l,t), magnetic energies = ',3e14.7)
! read namelist
      iuin = get_funit(iuin)
      open(unit=iuin,file='input1',form='formatted',status='old')
      read (iuin,input1)
! override input data
      idcode = 6
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
! boundary conditions
      ipbc = psolve
!     if (psolve==VACUUM_1D) ipbc = 2
! initialize for multiprocessing
!     ntasks = mpinit(sntasks)
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
! that is, fxyze are the convolutions of the total electric field
! over the particle shape
      allocate(fxyze(ndim+1,nxe))
! cu(i,j) = i component of current at grid (j).
! cus(i,j) = i component of acceleration density at grid (j).
! amu(i,j) = i component of momentum flux at grid (j). 
      allocate(cu(ndim,nxe),cus(ndim,nxe),amu(ndim,nxe))
! byze(i,j) = i component of magnetic field at grid (j).  
! byze is the convolution of the magnetic field over the particle shape
! eyze(i,j) = i component of transverse electric field at grid (j).  
! eyze is the convolution of the transverse electric field
! over the particle shape
      allocate(byze(ndim,nxe),eyze(ndim,nxe))
! ffc, ffe = form factor arrays for poisson solvers
      allocate(ffc(nxh),ffe(nxh))
! mixup = array of bit reversed addresses for fft
! sct = sine/cosine table for fft
      allocate(mixup(nxh),sct(nxh))
!
! non-periodic boundary conditions
!     indx1 = indx + 1
!     nxv = nx + 2; nx2 = 2*nx; nx2v = 2*nxv
! dirichlet conditions
!     if (psolve==DIRICHLET_2D) then
!        allocate(q2(nx2v),fx2(nx2v),fxyz2(3,nx2v))
!        allocate(cu2(2,nx2v),byz2(2,nx2v),ffd(nxv))
!        allocate(mixup2(nx),sct2(nx))
!     endif
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
      if (ipbc==2) then
         affp = real(nx-2)/real(np)
      endif
      omt = sqrt(omy*omy + omz*omz)
      q2m0 = qbme*qme*real(np)/real(nx)
      wp0 = q2m0*affp
! set initial time
      t0 = dt*real(itime0)
! set default diagnostic file names
      if (ntp > 0) fpname = 'potk1.'//cdrun
      if (nta > 0) faname = 'vpotk1.'//cdrun
! energy diagnostics
      if (ntw > 0) then
         allocate(wt((nloop-1)/ntw-(itime0/ntw)+1,7))
         itw = 0
      endif
! initialize external magnetic field
      byze = 0.0
      if (omt > 0) then
         call baddext(byze,omy,omz,nx,inorder)
         call cguard(byze,nx,inorder)
      endif
! prepare fft tables
      call fft_init(mixup,sct,indx)
! calculate form factors
      call pois_init(ffc,ax,affp,nx)
      call epois_init(ffe,ax,affp,wp0,ci,nx)
! dirichlet boundary conditions
!     if (psolve==DIRICHLET_2D) then
!        call fft_init(mixup2,sct2,indx1)
!        call fst_init(mixup,sct2,indx)
!        call poisd_init(ffd,ax,affp,nx)
!     endif
! debug
      if (ipbc==2) vdx = 0.
! initialize density profile and velocity distribution
! background electrons
      if (npx > 0) call distr(part,1,npx,vtx,vty,vtz,vx0,vy0,vz0,npx,nx,&
     &ipbc)
! beam electrons
      if (npxb > 0) call distr(part,npx1,npxb,vtdx,vtdy,vtdz,vdx,vdy,vdz&
     &,npxb,nx,ipbc)
! initialize charge density to background
      qi0 = -qme/affp
! fix guiding centers for electrons
      if (omt > 0) call distr(part,byze,np,qbme,nx,ipbc,inorder)
! sorting arrays
      if (sortime > 0) allocate(pt(np),ip(np),npic(nx1))
! initialize diagnostics
! diagnostic metafile
      iudm = get_funit(iudm)
! velocity diagnostic
      if (ntv > 0) then
         allocate(fv(2*nmv+2,ndim+1),fvm(3,ndim+1))
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
!        if (psolve==DIRICHLET_2D) then
!           allocate(sfield2(nx2v))
!           if (modesxp > nx) modesxp = nx
!           allocate(pots(modesxp))
! open output file
!           if (nprec==0) then
!              nprec = -1; iup = get_funit(iup)
!              call bfopen(pots,modesxp,iup,nprec,trim(fpname))
!           endif
!        endif
      endif
! vector potential diagnostic
      if (nta > 0) then
         allocate(vfield(ndim,nxe))
         if (modesxa > nxh) modesxa = nxh
         allocate(vpott(ndim,modesxa))
! open output file
         if (narec==0) then
            narec = -1; iua = get_funit(iua)
            call bfopen(vpott,modesxa,iua,narec,trim(faname))
         endif
!        if (psolve==DIRICHLET_2D) then
!           allocate(vfield2(2,nx2v))
!           if (modesxa > nx) modesxa = nx
!           allocate(vpots(2,modesxa))
!           if (narec==0) then
!              narec = -1; iua = get_funit(iua)
!              call bfopen(vpots,modesxa,iua,narec,trim(faname))
!           endif
!        endif
      endif
! record time
      call wtimer(time,ltime)
      write (iuot,*) 'initialization wall clock time = ', time, 'sec'
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
      write (iuot,991) ntime
! initialize current density and momentum flux to background
      if (ipbc==1) then
         call sguard(cu,zero,zero,nx,inorder)
      else if (ipbc==2) then
         call lsguard(cu,zero,zero,nx,inorder)
      endif
      call sguard(amu,zero,zero,nx,inorder)
! deposit current and momentum flux for electrons
      call djpost(part,cu,np,qme,zero,tdjpost,nx,ipbc,inorder,djopt)
      call dmjpost(part,amu,np,qme,tdcjpost,inorder,djopt)
! initialize charge density to background
      if (ipbc==1) then
         call sguard(qe,qi0,nx,inorder)
      else if (ipbc==2) then
         call lsguard(qe,qi0,nx,inorder)
      endif
! deposit electron charge
      call dpost(part,qe,np,qme,tdpost,inorder,dopt)
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
!
! periodic boundary conditions
!
      if (psolve==PERIODIC_2D) then
!
! add guard cells for current and momentum flux
      call aguard(cu,nx,inorder)
      call aguard(amu,nx,inorder)
! add guard cells for charge density
      call aguard(qe,nx,inorder)
! find maximum and minimum plasma frequency
      call wpmxn(qe,qi0,qbme,wpmax,wpmin,nx,inorder)
      wp0 = 0.5*(wpmax + wpmin)
! recalculate form factors
      if ((wp0 > 1.15*q2m0) .or. (wp0 < 0.85*q2m0)) then
         q2m0 = wp0
         wp0 = affp*wp0
         call epois_init(ffe,ax,affp,wp0,ci,nx)
         write (*,*) ntime, 'new shift constants,q2m0,wp0=', q2m0, wp0
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
! calculate longitudinal electric field in fourier space
      isign = -1
      call pois(qe,fxe,isign,ffc,we,nx,inorder)
! transform longitudinal electric field to real space
      isign = 1
      call fft(fxe,isign,mixup,sct,tfft,indx,inorder)
      call cguard(fxe,nx,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,inorder)
! calculate magnetic field in fourier space
      call bpois(cu,byze,ffc,ci,wm,nx,inorder)
! transform magnetic field to real space
      isign = 1
      call fft(byze,isign,mixup,sct,tfft,indx,inorder)
! add external magnetic field
      if (omt > 0) call baddext(byze,omy,omz,nx,inorder)
      call cguard(byze,nx,inorder)
! transform momentum flux to fourier space
      isign = -1
      call fft(amu,isign,mixup,sct,tfft,indx,inorder)
! take time derivative of current
      call dcuperp(cus,amu,nx,inorder)
! calculate convective part of transverse electric field
      call epois(cus,eyze,ffe,ci,wf,nx,inorder)
! transform transverse electric field to real space
      isign = 1
      call fft(eyze,isign,mixup,sct,tfft,indx,inorder)
      call cguard(eyze,nx,inorder)
! merge longitudinal and transverse electric fields
      fxyze(1,:) = fxe
      fxyze(2:3,:) = eyze
!
! inner iteration loop
      do k = 1, ndc
! initialize acceleration density with shift included
      call sguard(cus,eyze,q2m0,nx,inorder)
! initialize current and momentum flux
      call sguard(cu,zero,zero,nx,inorder)
      call sguard(amu,zero,zero,nx,inorder)
! deposit electron current and acceleration density and momentum flux
      call dcjpost(part,fxyze,byze,cu,cus,amu,omx,np,qme,qbme,dt,       &
     &tdcjpost,inorder,popt)
! add guard cells for current, acceleration density, and momentum flux
      call aguard(cu,nx,inorder)
      call aguard(cus,nx,inorder)
      call aguard(amu,nx,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,inorder)
! calculate magnetic field in fourier space
      call bpois(cu,byze,ffc,ci,wm,nx,inorder)
! transform magnetic field to real space
      isign = 1
      call fft(byze,isign,mixup,sct,tfft,indx,inorder)
! add external magnetic field
      if (omt > 0) call baddext(byze,omy,omz,nx,inorder)
      call cguard(byze,nx,inorder)
! transform acceleration density and momentum flux to fourier space
      isign = -1
      call fft(cus,isign,mixup,sct,tfft,indx,inorder)
      call fft(amu,isign,mixup,sct,tfft,indx,inorder)
! take time derivative of current
      call adcuperp(cus,amu,nx,inorder)
! calculate transverse electric field
      call epois(cus,eyze,ffe,ci,wf,nx,inorder)
! transform transverse electric field to real space
      isign = 1
      call fft(eyze,isign,mixup,sct,tfft,indx,inorder)
      call cguard(eyze,nx,inorder)
! copy transverse electric field
      fxyze(2:3,:) = eyze
      enddo
!
! vector potential diagnostic
      if (nta > 0) then
         it = ntime/nta
         if (ntime==nta*it) then
! calculate vector potential in fourier space
            call apois(cu,vfield,ffc,ci,wm,nx,inorder)
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
!
      endif
!
! push particles
      wke = 0.
      call push3(part,fxyze,byze,omx,np,qbme,dt,dt,wke,tpush,nx,ipbc,   &
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
            wef = we + wm
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
      write (iuot,*) 'bounded darwin code d0_dbeps1'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost + tdjpost + tdcjpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'electron current deposit time = ', tdjpost, 'sec'
      write (iuot,*) 'electron acceleration deposit time = ', tdcjpost, &
     &'sec'
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
! write out input file
      write (iudm,input1)
      write (iuot,*) ' * * * q.e.d. * * *'
      close(unit=iudm)
      close(unit=iuot)
! close graphics device
      call close_graphs
      stop
      end program d0_dbeps1
