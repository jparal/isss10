!-----------------------------------------------------------------------
!
      module init1d
!
! Fortran90 interface to 1d PIC Fortran77 library init1lib.f
! init1mod.f contains interface procedures to initialize particle
!            co-ordinates:
!            defines module init1d
! distr => idistr1 initializes x and vx co-ordinates for 1d code.
!          calls DISTR1
! distr => idistrh1 initializes x and vx,vy,vz co-ordinates for 
!          magnetized 1-2/2d codes.
!          calls DISTR1H
! distr => ibdistr1 calculates guiding centers for magnetized
!          1-2/2d codes.
!          calls GBDISTR1L 
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: july 16, 2011
!
      use input1d
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
      public :: ions1
      public :: ntd, ntj, npxi, npxbi
!     public :: ndd, ndj
      public :: qmi, rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0
      public :: vdxi, vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi
      public :: ndprofi, ampdxi, scaledxi, shiftdxi, modesxd, modesxj
      public :: pot1d, vpot1d, em1d
      public :: t0, ceng, nprec, fpname, narec, faname, nerec, fename
      public :: distr, fdistr, vdistr
      public :: den1d, vcur1d
      public :: ndrec, fdname, njrec, fjname 
!
! define interface to original Fortran77 procedures
      interface
         subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
         implicit none
         integer :: npx, idimp, nop, nx, ipbc
         real :: vtx, vdx
!        real, dimension(*) :: part
         real:: part
         end subroutine
      end interface
      interface
         subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop,n&
     &x,ipbc)
         implicit none
         integer :: npx, idimp, nop, nx, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(*) :: part
         real:: part
         end subroutine
      end interface
      interface
         subroutine FDISTR1(part,fnx,argx1,argx2,argx3,npx,idimp,nop,nx,&
     &ipbc,ierr)
         implicit none
         integer :: npx, idimp, nop, nx, ipbc, ierr
         real :: argx1, argx2, argx3
!        real, dimension(*) :: part
         real:: part
         real, external :: fnx
         end subroutine
      end interface
      interface
         subroutine VDISTR1(part,vtx,vdx,idimp,nop)
         implicit none
         integer :: idimp, nop
         real :: vtx, vdx
!        real, dimension(*) :: part
         real :: part
         end subroutine
      end interface
      interface
         subroutine VDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,idimp,nop)
         implicit none
         integer :: idimp, nop
         real :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(*) :: part
         real :: part
         end subroutine
      end interface
      interface
         subroutine GBDISTR1L(part,byz,qbm,idimp,nop,nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx,  nxv, ipbc
         real :: qbm
         real, dimension(idimp,nop) :: part
!        real, dimension(*) :: byz
         real :: byz
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface distr
         module procedure idistr1
         module procedure idistrh1
         module procedure ibdistr1
      end interface
!
      interface fdistr
         module procedure ifdistr1
      end interface
!
      interface vdistr
         module procedure ivdistr1
         module procedure ivdistrh1
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine idistr1(part,nstart,nop,vtx,vdx,npx,nx,ipbc)
! calculates initial particle co-ordinates and velocities in 1d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: nstart, nop, npx, nx, ipbc
         real :: vtx, vdx
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call DISTR1(part(1,nstart),vtx,vdx,npx,idimp,nop,nx,ipbc)
         end subroutine idistr1
!
         subroutine idistrh1(part,nstart,nop,vtx,vty,vtz,vdx,vdy,vdz,npx&
     &,nx,ipbc)
! calculates initial particle co-ordinates and velocities in 1-2/2d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: nstart, nop, npx, nx, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call DISTR1H(part(1,nstart),vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,n&
     &op,nx,ipbc)
         end subroutine idistrh1
!
         subroutine ifdistr1(part,nstart,nop,ampx,scalex,shiftx,npx,nx,i&
     &pbc,ndpro)
! calculates initial particle co-ordinates in 1d
! with various density profiles
         implicit none
         integer :: nstart, nop, npx, nx, ipbc, ndpro
         real :: ampx, scalex, shiftx
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, ierr
         real sxi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         idimp = size(part,1)
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call FDISTR1(part(1,nstart),FLDISTR1,zero,zero,zero,npx,idim&
     &p,nop,nx,ipbc,ierr)
! linear density
         else if (ndpro==1) then
            call FDISTR1(part(1,nstart),FLDISTR1,ampx,sxi,shiftx,npx,idi&
     &mp,nop,nx,ipbc,ierr)
! sinusoidal density
         else if (ndpro==2) then
            call FDISTR1(part(1,nstart),FSDISTR1,ampx,sxi,shiftx,npx,idi&
     &mp,nop,nx,ipbc,ierr)
! gaussian density
         else if (ndpro==3) then
            call FDISTR1(part(1,nstart),FGDISTR1,ampx,sxi,shiftx,npx,idi&
     &mp,nop,nx,ipbc,ierr)
! hyperbolic secant squared density
         else if (ndpro==4) then
            call FDISTR1(part(1,nstart),FHDISTR1,ampx,sxi,shiftx,npx,idi&
     &mp,nop,nx,ipbc,ierr)
         endif
         end subroutine ifdistr1
!
         subroutine ivdistr1(part,nstart,nop,vtx,vdx)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
         implicit none
         integer :: nstart, nop
         real :: vtx, vdx
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call VDISTR1(part(1,nstart),vtx,vdx,idimp,nop)
         end subroutine ivdistr1
!
         subroutine ivdistrh1(part,nstart,nop,vtx,vty,vtz,vdx,vdy,vdz)
! calculates initial particle velocities in 1-2/2d
! with maxwellian velocity with drift
         implicit none
         integer :: nstart, nop
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call VDISTR1H(part(1,nstart),vtx,vty,vtz,vdx,vdy,vdz,idimp,nop)
         end subroutine ivdistrh1
!
         subroutine ibdistr1(part,byz,nop,qbm,nx,ipbc,inorder)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for 1d
         implicit none
         integer :: nop, nx, ipbc
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: byz
! local data
         integer :: idimp, nxv, order
         idimp = size(part,1)
         nxv = size(byz,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(byz,1))
         case (2)
            if (order==LINEAR) then
               call GBDISTR1L(part,byz(1,1),qbm,idimp,nop,nx,nxv,ipbc)
            else
               call GBDISTR1L(part,byz(1,1),qbm,idimp,nop,nx,nxv,ipbc)
            endif
         end select
         end subroutine ibdistr1
!
      end module init1d
