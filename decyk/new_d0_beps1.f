!-----------------------------------------------------------------------
! * * * periodic 1d electrostatic particle simulation kernel code * * *
! this is a simple 1d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge, advancing particles, and solving the
! field.  the code moves electrons and ions, with electrostatic forces
! obtained by solving poisson's equation with fast fourier transforms
! for various boundary conditions. 
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G5
! copyright 1999, regents of the university of california
! update: july 21, 2010
      program d0_beps1
      use init1d
!     use espush1d
      use push1d
      use fft1d
      use field1d
      use dfield1d
      use diag1d
!     use mp0d, only: mpinit, ncpus
      implicit none
! idimp = dimension of phase space = 2
! nmv = number of segments in v for velocity distribution
      integer :: idimp = 2, nmv = 40
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iup = 11
      integer :: np, npi, npx1, nx, nxh, nxe, nx1, ipbc
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: it, itw
!     integer :: ntasks
      real :: zero = 0.0, wki = 0.0, time = 0.0, tloop = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0, tfft = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tsorti = 0.0
      real :: totpush = 0.0, totpushi = 0.0, ts = 0.0
      real :: qbme, qbmi, affp, qi0, we, wke, ws
      real :: vtxi, vtdxi
      real, dimension(:,:), pointer :: part, parti
      real, dimension(:), pointer :: qe, qi, fxe
      complex, dimension(:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:), pointer :: qt, sfield
      complex, dimension(:), pointer :: dent, pott
      real, dimension(:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
! dirichlet boundary conditions
      integer :: indx1, nxv, nx2, nx2v
!     real, dimension(:), pointer :: q2, fx2, sfield2
      complex, dimension(:), pointer :: ffd
      integer, dimension(:), pointer :: mixup2
      complex, dimension(:), pointer :: sct2
      real, dimension(:), pointer :: dens, pots
!
      character(len=10) :: cdrun
      character(len=32) :: fname
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
! read namelist
      iuin = get_funit(iuin)
      open(unit=iuin,file='input1',form='formatted',status='old')
      read (iuin,input1)
      if (movion==1) then
         rewind iuin
         read (iuin,ions1)
      endif
! override input data
      idcode = 4
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file
      iuot = get_funit(iuot)
      fname = 'output1.'//cdrun
      open(unit=iuot,file=trim(fname),form='formatted',status='replace')
! np = total number of electrons in simulation
      np = npx + npxb
! npi = total number of ions in simulation
      npi = npxi + npxbi
      nx = 2**indx; nxh = nx/2
      nxe = nx + 4
!     ax = .866025
      if (inorder==LINEAR) then
!        ax = .912871
         nxe = nx + 2
      endif
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
      allocate(part(idimp,np))
! in real space, qe(j) = charge density at grid point j
! in real space, qi(j) = ion charge density at grid point j
! in real space, fxe(j) = longitudinal force/charge at grid point j,
! that is, fxe is the convolution of the longitudinal electric field
! over the particle shape
      allocate(qe(nxe),qi(nxe),fxe(nxe))
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
!        allocate(q2(nx2v),fx2(nx2v))
         allocate(ffd(nxv),mixup2(nx),sct2(nx))
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
      if (movion==1) then
         qbmi = qmi/rmass
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
      endif
! set initial time
      t0 = dt*real(itime0)
! set default diagnostic file names
      if (ntd > 0) fdname = 'denk2.'//cdrun
      if (ntp > 0) fpname = 'potk1.'//cdrun
! energy time history
      if (ntw > 0) then
         allocate(wt((nloop-1)/ntw-(itime0/ntw)+1,4))
         itw = 0
      endif
! prepare fft tables
      call fft_init(mixup,sct,indx)
! calculate form factors
      call pois_init(ffc,ax,affp,nx)
! dirichlet boundary conditions
      if (psolve==DIRICHLET_2D) then
!        call fft_init(mixup2,sct2,indx1)
         call fst_init(mixup,sct2,indx)
         call poisd_init(ffd,ax,affp,nx)
      endif
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npi))
      endif
! debug
!     if (ipbc==2) vdx = 0.
! initialize density profile and velocity distribution
      npx1 = npx + 1
! background electrons
!     if (npx > 0) call distr(part,1,npx,vtx,vx0,npx,nx,ipbc)
      if (npx > 0) then
         call fdistr(part,1,npx,ampdx,scaledx,shiftdx,npx,nx,ipbc,ndprof&
     &)
         call vdistr(part,1,npx,vtx,vx0)
      endif
! beam electrons
!     if (npxb > 0) call distr(part,npx1,npxb,vtdx,vdx,npxb,nx,ipbc)
      if (npxb > 0) then
         call fdistr(part,npx1,npxb,ampdx,scaledx,shiftdx,npxb,nx,ipbc, &
     &ndprof)
         call vdistr(part,npx1,npxb,vtdx,vdx)
      endif
! initialize background charge density
      if (movion==0) then
         qi0 = -qme/affp
         if (ipbc==1) then
            call sguard(qi,zero,nx,inorder)
         else if (ipbc==2) then
            call lsguard(qi,zero,nx,inorder)
         endif
         call dpost(part,qi,np,-qme,tdpost,inorder,dopt)
         if (ipbc==1) then
            call aguard(qi,nx,inorder)
         else if (ipbc==2) then
            call laguard(qi,nx,inorder)
         endif
! debug
!        if (ipbc==1) then
!           call sguard(qi,qi0,nx,inorder)
!        else if (ipbc==2) then
!           call lsguard(qi,qi0,nx,inorder)
!        endif
      endif
! initialize ions
      if (movion==1) then
         npx1 = npxi + 1
! background ions
!        if (npxi > 0) call distr(parti,1,npxi,vtxi,vxi0,npxi,nx,ipbc)
         if (npxi > 0) then
            call fdistr(parti,1,npxi,ampdxi,scaledxi,shiftdxi,npxi,nx,  &
     &ipbc,ndprofi)
            call vdistr(parti,1,npxi,vtxi,vxi0)
         endif
! beam ions
!        if (npxbi > 0) call distr(parti,npx1,npxbi,vtdxi,vdxi,npxbi,nx,&
!    &ipbc)
         if (npxbi > 0) then
            call fdistr(parti,npx1,npxbi,ampdxi,scaledxi,shiftdxi,npxbi,&
     &nx,ipbc,ndprofi)
            call vdistr(parti,npx1,npxbi,vtdxi,vdxi)
         endif
      endif
!
! sorting arrays for electrons
      if (sortime > 0) allocate(pt(np),ip(np),npic(nx1))
! initialize diagnostics
! diagnostic metafile
      iudm = get_funit(iudm)
! velocity diagnostic
      if (ntv > 0) then
         allocate(fv(2*nmv+2,1),fvm(3,1))
         fv(1,:) = 8.0*vtx
! ions
         if (movion==1) then
            allocate(fvi(2*nmv+2,1),fvmi(3,1))
            fvi(1,:) = 8.0*vtxi
         endif
      endif
! ion density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0)) then
         allocate(sfield(nxe))
!        if (psolve==DIRICHLET_2D) then
!           allocate(sfield2(nx2v))
!        endif
      endif
! ion density diagnostics
      if (ntd > 0) then
         allocate(qt(nxe))
         if (modesxd > nxh) modesxd = nxh
         allocate(dent(modesxd))
! open output file
         if (ndrec==0) then
            ndrec = -1; iud = get_funit(iud)
            call bfopen(dent,modesxd,iud,ndrec,trim(fdname))
         endif
         if (psolve==DIRICHLET_2D) then
            if (modesxd > nx) modesxd = nx
            allocate(dens(modesxd))
! open output file
            if (ndrec==0) then
               ndrec = -1; iud = get_funit(iud)
               call bfopen(dens,modesxd,iud,ndrec,trim(fdname))
            endif
         endif
      endif
! potential diagnostic
      if (ntp > 0) then
         if (modesxp > nxh) modesxp = nxh
         allocate(pott(modesxp))
! open output file
         if (nprec==0) then
            nprec = -1; iup = get_funit(iup)
            call bfopen(pott,modesxp,iup,nprec,trim(fpname))
         endif
         if (psolve==DIRICHLET_2D) then
            if (modesxp > nx) modesxp = nx
            allocate(pots(modesxp))
! open output file
            if (nprec==0) then
               nprec = -1; iup = get_funit(iup)
               call bfopen(pots,modesxp,iup,nprec,trim(fpname))
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
! deposit electron charge density
! initialize electron charge density to zero
      if (ipbc==1) then
         call sguard(qe,zero,nx,inorder)
      else if (ipbc==2) then
         call lsguard(qe,zero,nx,inorder)
      endif
! deposit electron charge
      call dpost(part,qe,np,qme,tdpost,inorder,dopt)
! add guard cells
      if (ipbc==1) then
          call aguard(qe,nx,inorder)
      else if (ipbc==2) then
          call laguard(qe,nx,inorder)
      endif
! deposit ion charge density
      if (movion==1) then
         if (ipbc==1) then
            call sguard(qi,zero,nx,inorder)
         else if (ipbc==2) then
            call lsguard(qi,zero,nx,inorder)
         endif
         call dpost(parti,qi,npi,qmi,tdposti,inorder,dopt)
         if (ipbc==1) then
             call aguard(qi,nx,inorder)
         else if (ipbc==2) then
             call laguard(qi,nx,inorder)
         endif
      endif
! add electron and ion densities
      call addqei(qe,qi,nx,inorder)
! velocity diagnostic
      if (ntv > 0) then
         it = ntime/ntv
         if (ntime==ntv*it) then
! calculate paticle distribution function and moments
            call vdist(part,fv,fvm,np,nmv)
! display velocity distributions
            call displayfv(fv,fvm,' ELECTRON',ntime,nmv,2,irc)
            if (irc==1) go to 2000
! ion distribution function
            if (movion==1) then
               call vdist(parti,fvi,fvmi,npi,nmv)
               call displayfv(fvi,fvmi,' ION',ntime,nmv,2,irc)
               if (irc==1) go to 2000
            endif
         endif
      endif
!
! dirichlet boundary conditions
!
      if (psolve==DIRICHLET_2D) then
!
! ion density diagnostic
      if (ntd > 0) then
         it = ntime/ntd
         if (ntime==ntd*it) then
            qt = qi
! transform ion density to fourier space
            isign = -1
!           call dblsin(qt,q2,nx,inorder)
!           call fft(q2,isign,mixup2,sct2,tfft,indx1,LINEAR)
            call fst(qt,isign,mixup,sct2,tfft,indx,inorder)
! calculate smoothing in fourier space
            isign = 2
!           call poisdx(q2,sfield2,isign,ffd,ws,nx)
!           call cmfieldd(sfield2,sfield,nx,inorder)
            call poisd(qt,sfield,isign,ffd,ws,nx,inorder)
! store selected fourier modes
!           call gtmodes(sfield2,dent,nx2,modesxd,LINEAR)
            call gtsmodes(sfield,dens,nx,modesxd,inorder)
! write diagnostic output
!           call writebf(dent,modesxd,iud,ndrec,order=LINEAR)
            call writebf(dens,modesxd,iud,ndrec,order=LINEAR)
! transform ion density to real space
            isign = 1
!           call fft(sfield2,isign,mixup2,sct2,tfft,indx1,LINEAR)
!           call hafdbl(sfield,sfield2,nx,inorder)
            call fst(sfield,isign,mixup,sct2,tfft,indx,inorder)
            call lcguard(sfield,nx,inorder)
! display ion density
            call displays(sfield,' ION DENSITY',ntime,999,0,nx,irc,     &
     &inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      isign = -1
!     call dblsin(qe,q2,nx,inorder)
!     call fft(q2,isign,mixup2,sct2,tfft,indx1,LINEAR)
      call fst(qe,isign,mixup,sct2,tfft,indx,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! calculate potential in fourier space
            isign = 1
!           call poisdx(q2,sfield2,isign,ffd,ws,nx)
!           call cmfieldd(sfield2,sfield,nx,inorder)
            call poisd(qe,sfield,isign,ffd,ws,nx,inorder)
! store selected fourier modes
!           call gtmodes(sfield2,pott,nx2,modesxp,LINEAR)
            call gtsmodes(sfield,pots,nx,modesxp,inorder)
! write diagnostic output
!           call writebf(pott,modesxp,iup,nprec,order=LINEAR)
            call writebf(pots,modesxp,iup,nprec,order=LINEAR)
! transform potential to real space
            isign = 1
!           call fft(sfield2,isign,mixup2,sct2,tfft,indx1,LINEAR)
!           call hafdbl(sfield,sfield2,nx,inorder)
            call fst(sfield,isign,mixup,sct2,tfft,indx,inorder)
            call lcguard(sfield,nx,inorder)
! display potential
            call displays(sfield,' POTENTIAL',ntime,999,0,nx,irc,inorder&
     &)
            if (irc==1) go to 2000
         endif
      endif
! calculate force/charge in fourier space
      isign = -1
!     call poisdx(q2,fx2,isign,ffd,we,nx)
      call poisd(qe,fxe,isign,ffd,we,nx,inorder)
! transform force/charge to real space
      isign = 1
!     call fft(fx2,isign,mixup2,sct2,tfft,indx1,LINEAR)
!     call hafdbl(fxe,fx2,nx,inorder)
      call fct(fxe,isign,mixup,sct2,tfft,indx,inorder)
      call lcguard(fxe,nx,inorder)
!
! periodic boundary conditions
!
      else if (psolve==PERIODIC_2D) then
!
! ion density diagnostic
      if (ntd > 0) then
         it = ntime/ntd
         if (ntime==ntd*it) then
            qt = qi
! transform ion density to fourier space
            isign = -1
            call fft(qt,isign,mixup,sct,tfft,indx,inorder)
! calculate smoothing in fourier space
            isign = 2
            call pois(qt,sfield,isign,ffc,ws,nx,inorder)
! store selected fourier modes
            call gtmodes(sfield,dent,nx,modesxd,inorder)
! write diagnostic output
            call writebf(dent,modesxd,iud,ndrec,order=LINEAR)
! transform ion density to real space
            call fft(sfield,isign,mixup,sct,tfft,indx,inorder)
            call cguard(sfield,nx,inorder)
! display ion density
            call displays(sfield,' ION DENSITY',ntime,999,0,nx,irc,     &
     &inorder)
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
! calculate force/charge in fourier space
      isign = -1
      call pois(qe,fxe,isign,ffc,we,nx,inorder)
! transform force/charge to real space
      isign = 1
      call fft(fxe,isign,mixup,sct,tfft,indx,inorder)
      call cguard(fxe,nx,inorder)
!
      endif
!
! push electrons
      wke = 0.0
      call push(part,fxe,np,qbme,dt,wke,tpush,nx,ipbc,inorder,popt)
! push ions
      if (movion==1) then
         wki = 0.0
         call push(parti,fxe,npi,qbmi,dt,wki,tpushi,nx,ipbc,inorder,popt&
     &)
         wki = wki*rmass
      endif
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
            ws = we + wke + wki
            write (iuot,992) we, wke, ws
            itw = itw + 1
            wt(itw,:) = (/we,wke,wki,ws/)
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
      write (iuot,*) 'bounded electrostatic code d0_beps1'
!     write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'total electron push time = ', totpush, 'sec'
      write (iuot,*) 'electron sort time = ', tsort
      totpush = totpush + tsort
      write (iuot,*) 'total electron time = ', totpush, 'sec'
      if (movion==1) then
         totpushi = tpushi + tdposti
         write (iuot,*) 'ion push time = ', tpushi, 'sec'
         write (iuot,*) 'ion charge deposit time = ', tdposti, 'sec'
         write (iuot,*) 'total ion push time = ', totpushi, 'sec'
         write (iuot,*) 'ion sort time = ', tsorti
         totpushi = totpushi + tsorti
         write (iuot,*) 'total ion time = ', totpushi, 'sec'
      endif
      write (iuot,*) 'total fft time=', tfft, 'sec'
      time = time - (totpush + totpushi + tfft)
      write (iuot,*) 'other time=', time, 'sec'
! write final diagnostic metafile
      fname = 'diag1.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! ion density diagnostics
      if (ntd > 0) then
         ndrec = ndrec - 1
         ceng = zero
         write (iudm,den1d)
      endif
! potential diagnostics
      if (ntp > 0) then
         nprec = nprec - 1
         ceng = affp
         write (iudm,pot1d)
      endif
! write out input file
      write (iudm,input1)
      write (iuot,*) ' * * * q.e.d. * * *'
      close(unit=iudm)
      close(unit=iuot)
! close graphics device
      call close_graphs
      stop
      end program d0_beps1
