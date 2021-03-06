!-----------------------------------------------------------------------
! * * * periodic 1d electrostatic particle simulation kernel code * * *
! this is a simple 1d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge, advancing particles, and solving the
! field.  the code moves only electrons, with periodic electrostatic
! forces obtained by solving poisson's equation with fast fourier
! transforms. 
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G3
! copyright 1999, regents of the university of california
! update: june 22, 2011
      program beps1
      use init1d
!     use espush1d
      use push1d
      use fft1d
      use field1d
      use diag1d
!     use mp0d, only: mpinit, ncpus
      implicit none
! idimp = dimension of phase space = 2
! nmv = number of segments in v for velocity distribution
      integer :: idimp = 2, nmv = 40, ipbc = 1
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iup = 11
      integer :: np, npx1, nx, nxh, nxe, nx1
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: it, itw
!     integer :: ntasks
      real :: zero = 0.0, time = 0.0, tloop = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0, tfft = 0.0
      real :: totpush = 0.0, ts = 0.0
      real :: qbme, affp, qi0, we, wke, ws
      real, dimension(:,:), pointer :: part
      real, dimension(:), pointer :: qe, fxe
      complex, dimension(:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:), pointer :: sfield
      complex, dimension(:), pointer :: pott
      real, dimension(:,:), pointer :: fv, fvm
      real, dimension(:,:), pointer :: wt
      character(len=10) :: cdrun
      character(len=32) :: fname
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
! read namelist
      iuin = get_funit(iuin)
      open(unit=iuin,file='input1',form='formatted',status='old')
      read (iuin,input1)
! override input data
      idcode = 1
      psolve = 1
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
! in real space, fxe(j) = longitudinal force/charge at grid point j,
! that is, fxe is the convolution of the longitudinal electric field
! over the particle shape
      allocate(qe(nxe),fxe(nxe))
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
! set initial time
      t0 = dt*real(itime0)
! set default diagnostic file names
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
! initialize density profile and velocity distribution
! background electrons
      if (npx > 0) call distr(part,1,npx,vtx,vx0,npx,nx,ipbc)
! beam electrons
      if (npxb > 0) call distr(part,npx1,npxb,vtdx,vdx,npxb,nx,ipbc)
! initialize charge density to background
      qi0 = -qme/affp
! sorting arrays
      if (sortime > 0) allocate(pt(np),ip(np),npic(nx1))
! initialize diagnostics
! diagnostic metafile
      iudm = get_funit(iudm)
! velocity diagnostic
      if (ntv > 0) then
         allocate(fv(2*nmv+2,1),fvm(3,1))
         fv(1,:) = 8.*vtx
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
! record time
      call wtimer(time,ltime)
      write (iuot,*) 'initialization wall clock time = ', time, 'sec'
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
      write (iuot,991) ntime
! initialize charge density to background
      call sguard(qe,qi0,nx,inorder)
! deposit electron charge
      call dpost(part,qe,np,qme,tdpost,inorder,dopt)
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
! calculate force/charge in fourier space
      isign = -1
      call pois(qe,fxe,isign,ffc,we,nx,inorder)
! transform force/charge to real space
      isign = 1
      call fft(fxe,isign,mixup,sct,tfft,indx,inorder)
      call cguard(fxe,nx,inorder)
! push particles
      wke = 0.
      call push(part,fxe,np,qbme,dt,wke,tpush,nx,ipbc,inorder,popt)
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
            ws = we + wke
            write (iuot,992) we, wke, ws
            itw = itw + 1
            wt(itw,:) = (/we,wke,zero,ws/)
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
      write (iuot,*) 'electrostatic code beps1'
!     write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'total electron push time = ', totpush, 'sec'
      write (iuot,*) 'electron sort time = ', tsort, 'sec'
      totpush = totpush + tsort
      write (iuot,*) 'total electron time = ', totpush, 'sec'
      write (iuot,*) 'total fft time = ', tfft, 'sec'
      time = time - (totpush + tfft)
      write (iuot,*) 'other time = ', time, 'sec'
! write final diagnostic metafile
      fname = 'diag1.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
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
      end program beps1
