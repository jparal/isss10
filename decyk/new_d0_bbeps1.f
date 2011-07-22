!-----------------------------------------------------------------------
! * * * periodic 1d electromagnetic particle simulation kernel code * *
! this is a simple 1d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves only electrons, with
! electromagnetic forces obtained by solving maxwell's equation with
! fast fourier transforms for various boundary conditions. 
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G5
! copyright 1999, regents of the university of california
! update: august 18, 2010
      program d0_bbeps1
      use init1d
!     use empush1d
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
      integer :: iue = 26
      integer :: np, npx1, nx, nxh, nxe, nxeh, nx1, ipbc
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
! dirichlet boundary conditions
      integer :: indx1, nxv, nx2, nx2v
      real, dimension(:), pointer :: q2, fx2
      real, dimension(:,:), pointer :: cu2, amu2, fxyz2, byz2
      complex, dimension(:), pointer :: ffd
      integer, dimension(:), pointer :: mixup2
      complex, dimension(:), pointer :: sct2
      real, dimension(:), pointer :: sfield2, pots
      real, dimension(:,:), pointer :: vfield2, vpots, vpotsr
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
      idcode = 5
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
! non-periodic boundary conditions
      indx1 = indx + 1
      nxv = nx + 2; nx2 = 2*nx; nx2v = 2*nxv
! dirichlet conditions
      if (psolve==DIRICHLET_2D) then
         allocate(q2(nx2v),fx2(nx2v),fxyz2(3,nx2v))
         allocate(cu2(2,nx2v),byz2(2,nx2v),ffd(nxv))
         allocate(mixup2(nx),sct2(nx))
      endif
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
      if (omt > 0) then
         call baddext(byze,omy,omz,nx,inorder)
         call cguard(byze,nx,inorder)
      endif
      byz = cmplx(0.0,0.0)
      eyz = cmplx(0.0,0.0)
      cu = 0.0
      if (psolve==DIRICHLET_2D) cu2 = 0.0
! prepare fft tables
      call fft_init(mixup,sct,indx)
! calculate form factors
      call pois_init(ffc,ax,affp,nx)
! dirichlet boundary conditions
      if (psolve==DIRICHLET_2D) then
         call fft_init(mixup2,sct2,indx1)
         call fst_init(mixup,sct2,indx)
         call poisd_init(ffd,ax,affp,nx)
      endif
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
         if (psolve==DIRICHLET_2D) then
            allocate(sfield2(nx2v))
            if (modesxp > nx) modesxp = nx
            allocate(pots(modesxp))
! open output file
            if (nprec==0) then
               nprec = -1; iup = get_funit(iup)
               call bfopen(pots,modesxp,iup,nprec,trim(fpname))
            endif
         endif
      endif
! vector potential or electromagnetic diagnostic
      if ((nta > 0) .or. (nte > 0)) then
         allocate(vfield(2,nxe))
         vfield = 0.0
         if (psolve==DIRICHLET_2D) then
            allocate(vfield2(2,nx2v))
            vfield2 = 0.0
         endif
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
         if (psolve==DIRICHLET_2D) then
            if (modesxa > nx) modesxa = nx
            allocate(vpots(2,modesxa))
            if (narec==0) then
               narec = -1; iua = get_funit(iua)
               call bfopen(vpots,modesxa,iua,narec,trim(faname))
            endif
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
         if (psolve==DIRICHLET_2D) then
            if (modesxe > nx) modesxe = nx
            allocate(vpotsr(2,modesxe))
            if (nerec==0) then
               nerec = -1; iue = get_funit(iue)
               call bfopen(vpotsr,modesxe,iue,nerec,trim(fename))
            endif
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
         if (ntime==nte*it) then
            vfield = cu
            if (psolve==DIRICHLET_2D) vfield2 = cu2
         endif
      endif
! initialize current density to background
      if (ipbc==1) then
         call sguard(cu,zero,zero,nx,inorder)
      else if (ipbc==2) then
         call lsguard(cu,zero,zero,nx,inorder)
      endif
! deposit electron current
      call djpost(part,cu,np,qme,dth,tdjpost,nx,ipbc,inorder,djopt)
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
! dirichlet boundary conditions
!
      if (psolve==DIRICHLET_2D) then
!
! add guard cells
      call laguard(qe,nx,inorder)
! transform charge to fourier space
      isign = -1
      call dblsin(qe,q2,nx,inorder)
      call fft(q2,isign,mixup2,sct2,tfft,indx1,LINEAR)
!     call fst(qe,isign,mixup,sct2,tfft,indx,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! calculate potential in fourier space
            isign = 1
            call poisdx(q2,sfield2,isign,ffd,we,nx)
            call cmfieldd(sfield2,sfield,nx,inorder)
!           call poisd(qe,sfield,isign,ffd,we,nx,inorder)
! store selected fourier modes
!           call gtmodes(sfield2,pott,nx2,modesxp,LINEAR)
            call gtsmodes(sfield,pots,nx,modesxp,inorder)
! write diagnostic output
!           call writebf(pott,modesxp,iup,nprec,order=LINEAR)
            call writebf(pots,modesxp,iup,nprec,order=LINEAR)
! transform potential to real space
            isign = 1
            call fft(sfield2,isign,mixup2,sct2,tfft,indx1,LINEAR)
            call hafdbl(sfield,sfield2,nx,inorder)
!           call fst(sfield,isign,mixup,sct2,tfft,indx,inorder)
            call lcguard(sfield,nx,inorder)
! display potential
            call displays(sfield,' POTENTIAL',ntime,999,0,nx,irc,inorder&
     &)
            if (irc==1) go to 2000
         endif
      endif
! add guard cells for current
      call laguard(cu,nx,inorder)
! transform current to fourier space
      isign = -1
      call dblsin(cu,cu2,nx,inorder)
      call fft(cu2,isign,mixup2,sct2,tfft,indx1,LINEAR)
      call cmfieldd(cu2,cu,nx,inorder)
!     call fst(cu,isign,mixup,sct2,tfft,indx,inorder)
! electromagnetic diagnostic
      if (nte > 0) then
         it = ntime/nte
         if (ntime==nte*it) then
! calculate averaged radiative vector potential
!           vfield2 = 0.5*(vfield2 + cu2)
!           call avrpotdx(vfield2,byz2,ffd,ci,nx)
!           call cmfieldd(vfield2,vfield,nx,inorder)
            vfield = 0.5*(vfield + cu)
            call avrpotd(vfield,byz,ffd,ci,nx,inorder)
! store selected fourier modes
!           call gtmodes(vfield2,vpotr,nx2,modesxe,LINEAR)
            call gtsmodes(vfield,vpotsr,nx,modesxe,inorder)
! write diagnostic output
!           call writebf(vpotr,modesxe,iue,nerec,order=LINEAR)
            call writebf(vpotsr,modesxe,iue,nerec,order=LINEAR)
! transform radiative vector potential to real space
            isign = 1
            call fft(vfield2,isign,mixup2,sct2,tfft,indx1,LINEAR)
            call hafdbl(vfield,vfield2,nx,inorder)
!           call fst(vfield,isign,mixup,sct2,tfft,indx,inorder)
            call lcguard(vfield,nx,inorder)
! display radiative vector potential
            call displayv(vfield,' RADIATIVE VPOTENTIAL',ntime,999,0,1, &
     &nx,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! calculate electromagnetic fields in fourier space
      if (ntime==0) then
! calculate initial darwin magnetic field
         call ibpoisd(cu,byz,ffd,ci,wm,nx,inorder)
         wf = 0.
! calculate initial darwin electric field
!        allocate(amu2(2,nx2v),ffe(nxv))
!        allocate(amu(2,nxe),ffe(nxeh))
! deposit momentum flux
!        call sguard(amu,zero,zero,nx,inorder)
!        call dmjpost(part,amu,np,qme,ts,inorder,djopt)
!        call aguard(amu,nx,inorder)
! solve for darwin electric field
!        isign = -1
!        call fft(amu,isign,mixup,sct,tfft,indx,inorder)
!        byze = cu
!        call dcuperp(byze,amu,nx,inorder)
!        call epois_init(ffe,ax,affp,wp0,ci,nx)
!        call iepois(byze,eyz,ffe,ci,wf,nx,inorder)
!        deallocate(amu2,ffe)
!        deallocate(amu,ffe)
         dth = .5*dt
! calculate electromagnetic fields
      else
         call maxweld(eyz,byz,cu,ffd,ci,dt,wf,wm,nx,inorder)
      endif
! vector potential diagnostic
      if (nta > 0) then
         it = ntime/nta
         if (ntime==nta*it) then
! calculate vector potential in fourier space
!           call avpotdx(byz2,vfield2,nx)
!           call cmfieldd(vfield2,vfield,nx,inorder)
            call avpotd(byz,vfield,nx,inorder)
! store selected fourier modes
!           call gtmodes(vfield2,vpott,nx2,modesxa,LINEAR)
            call gtsmodes(vfield,vpots,nx,modesxa,inorder)
! write diagnostic output
!           call writebf(vpott,modesxa,iua,narec,order=LINEAR)
            call writebf(vpots,modesxa,iua,narec,order=LINEAR)
! transform vector potential to real space
            isign = 1
            call fft(vfield2,isign,mixup2,sct2,tfft,indx1,LINEAR)
            call hafdbl(vfield,vfield2,nx,inorder)
!           call fst(vfield,isign,mixup,sct2,tfft,indx,inorder)
            call lcguard(vfield,nx,inorder)
! display vector potential
            call displayv(vfield,' VECTOR POTENTIAL',ntime,999,0,1,nx,  &
     &irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! calculate longitudinal electric field in fourier space
      isign = -1
      call poisdx(q2,fx2,isign,ffd,we,nx)
!     call poisd(qe,fxe,isign,ffd,we,nx,inorder)
! add longitudinal and transverse electric fields
      call emfieldd(fxyz2,fx2,eyz,ffd,nx)
!     call emfieldr(fxyze,fxe,eyz,ffd,nx,inorder)
! copy magnetic field
      call emfieldd(byz2,byz,ffd,nx)
!     call emfieldr(byze,byz,ffd,nx,inorder)
! transform force/charge to real space
      isign = 1
      call fft(fxyz2,isign,mixup2,sct2,tfft,indx1,LINEAR)
      call hafdbl(fxyze,fxyz2,nx,inorder)
!     call fcst(fxyze,isign,mixup,sct2,tfft,indx,inorder)
      call lcguard(fxyze,nx,inorder)
! transform magnetic field to real space
      isign = 1
      call fft(byz2,isign,mixup2,sct2,tfft,indx1,LINEAR)
      call hafdbl(byze,byz2,nx,inorder)
!     call fct(byze,isign,mixup,sct2,tfft,indx,inorder)
! add external magnetic field
      if (omt > 0) call baddext(byze,omy,omz,nx,inorder)
      call lcguard(byze,nx,inorder)
!
! periodic boundary conditions
!
      else if (psolve==PERIODIC_2D) then
!
! add guard cells
      call aguard(qe,nx,inorder)
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
! add guard cells for current
      call aguard(cu,nx,inorder)
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
         byze = cu
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
      if (omt > 0) call baddext(byze,omy,omz,nx,inorder)
      call cguard(byze,nx,inorder)
!
      endif
!
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
      write (iuot,*) 'bounded electromagnetic code d0_bbeps1'
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
      end program d0_bbeps1
