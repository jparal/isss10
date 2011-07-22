!-----------------------------------------------------------------------
!
      module input1d
!
! input1mod.f defines namelists containing input and output variables.
! written by viktor k. decyk, ucla
! copyright 2011, regents of the university of california
! update: july 6, 2011
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR,&
     & PERIODIC_2D, DIRICHLET_2D
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PERIODIC_2D, DIRICHLET_2D
      public :: input1
      public :: idrun, idrun0, idcode, indx, npx, npxb, inorder
      public :: popt, dopt, djopt, nustrt, ntr
      public :: ntw, ntp, nta, ntv, nts, nte
!     public :: ndw, ndp, nda, ndv, nds, nde
      public :: tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0
      public :: vdx, vdy, vdz, vtdx, vtdy, vtdz
      public :: psolve
!     public :: relativity
      public :: omx, omy, omz, ci, ax
      public :: ndim
      public :: ndc, movion, sortime, nplot, sntasks
      public :: ndprof, ampdx, scaledx, shiftdx
      public :: modesxp, modesxa, modesxe
      public :: pot1d, vpot1d, em1d
      public :: t0, ceng, nprec, fpname, narec, faname, nerec, fename
!
! Namelist Input
      save
! idrun/idrun0 = run identifier for current/old run
! idcode = code identifier
      integer :: idrun = 1, idrun0 = 0, idcode = 0
! indx = exponent which determines length in x direction, nx=2**indx
! npx = number of background particles distributed in x direction
! npxb = number of beam particles in x direction
!     integer :: indx =   7, npx =    5120, npxb = 512
      integer :: indx =  11, npx =   409600, npxb = 40960
!     integer :: indx =  14, npx =  4096000, npxb =   409600
!     integer :: indx =  18, npx = 40960000, npxb =  4096000
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
      integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
      integer :: djopt = STANDARD
! nustrt = (0,1,2) = this is an (old start,new start,restart)
! ntr = number of time steps between restart routine
      integer :: nustrt = 1, ntr = 0
! ntw, ndw = number of time steps between energy diagnostic
! ntp, ndp = number of time steps between potential diagnostic
! nta, nda = number of time steps between vector potential diagnostic
! ntv, ndv = number of time steps between velocity-space diagnostic
! nts, nds = number of time steps between phase space diagnostic
      integer :: ntw = 1, ntp = 0, nta = 0, ntv = 0, nts = 0
!     integer :: ndw = 1, ndp = 0, nda = 0, ndv = 0, nds = 0
! nte, nde = number of time steps between electromagnetic diagnostic
      integer :: nte = 0
!     integer :: nde = 0
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend =  45.000, dt = 0.2000000e+00
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! psolve = type of poisson solver = (1,2,3)
      integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
!     integer :: relativity = 0
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! ci = reciprical of velocity of light
      real :: ci = 0.1
! ax = half-width of particle in x direction
!     real :: ax = .816497
!     real :: ax = .866025
      real :: ax = .912871
! ndim = number of velocity dimensions = 1 or 2
      integer :: ndim = 2
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
! movion = (0,1) = (no,yes) move the ions
      integer :: movion = 0
! sortime = number of time steps between electron sorting
      integer :: sortime = 100
! nplot = maximum number of plots per page
      integer :: nplot = 4
! sntasks = (-1,n) = set maximum number of tasks (-1 = number of cpus-1)
      integer :: sntasks = -1
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
      integer :: ndprof = 0
! ampdx = amplitude of density compared to uniform in x
! scaledx = scale length for spatial coordinate in x
! shiftdx = shift of spatial coordinate in x
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
! modesxd = number of modes in x to keep for ion density diagnostic
!     integer :: modesxd = 11
! modesxp = number of modes in x to keep for potential diagnostic
! modesxa = number of modes in x to keep for vector potential diagnostic
! modesxe = number of modes in x to keep for electromagnetic diagnostic
!           (radiative part of vector potential)
      integer :: modesxp = 11, modesxa = 11, modesxe = 11
! define namelist
      namelist /input1/ idrun, idrun0, idcode, indx, npx, npxb, inorder,&
     & popt, dopt, djopt, nustrt, ntr, ntw, ntp, nta, ntv, nts, nte,    &
     &tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0, vdx, vdy, vdz, vtdx, &
     &vtdy, vtdz, psolve, omx, omy, omz, ci, ax, ndim, ndc, movion,     &
     &sortime, nplot, sntasks, ndprof, ampdx, scaledx, shiftdx, modesxp,&
     &modesxa, modesxe
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
!
! Namelist output for potential diagnostic
! nprec = current record number for potential writes
      integer :: nprec = 0
! fpname = file name for potential diagnostic
      character(len=32) :: fpname = 'potk1.0'
! define namelist
      namelist /pot1d/ idrun, indx, ntp, modesxp, psolve, omx, omy, omz,&
     & nprec, t0, tend, dt, ceng, fpname
!
! Namelist output for vector potential diagnostic
! narec = current record number for vector potential writes
      integer :: narec = 0
! faname = file name for vector potential diagnostic
      character(len=32) :: faname = 'vpotk1.0'
! define namelist
      namelist /vpot1d/ idrun, indx, nta, modesxa, psolve, ndim, omx,   &
     &omy, omz, ci, narec, t0, tend, dt, ceng, faname
!
! Namelist output for electromagnetic diagnostic
! nerec = current record number for electromagnetic writes
      integer :: nerec = 0
! fename = file name for electromagnetic diagnostic
      character(len=32) :: fename = 'vpotrk1.0'
! define namelist
      namelist /em1d/ idrun, indx, nte, modesxe, psolve, ndim, omx, omy,&
     &omz, ci, nerec, t0, tend, dt, ceng, fename
!
      end module input1d
