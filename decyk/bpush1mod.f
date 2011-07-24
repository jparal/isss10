!-----------------------------------------------------------------------
!
      module bpush1d
!
! Fortran90 interface to 1d PIC Fortran77 library bpush1lib.f
! bpush1mod.f contains interface procedures to process particles with
!             magnetic fields:
!             defines module bpush1d
! djpost => igjpost1 deposits current density, with various
!           interpolations and optimizations.
!           calls GJPOST1, GSJPOST1, GSJPOST1X, GJPOST1L, GSJPOST1L, or
!           GSJPOST1XL
! push3 => igbpush13 push particles with magnetic field, with various
!          interpolations and optimizations.
!          calls GBPUSH13, GSBPUSH13, GBPUSH13L, or GSBPUSH13L
! retard => iretard1 retard particle position a half time-step.
!           calls RETARD1
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: july 12, 2011
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use diag1d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer
      public :: djpost, push3, retard
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GJPOST1(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer :: nop, idimp, nx, nxv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST1(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer :: nop, idimp, nx, nxv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST1X(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer :: nop, idimp, nx, nxv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxv) :: cu
         end subroutine
      end interface
      interface
         subroutine GJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer :: nop, idimp, nx, nxv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer :: nop, idimp, nx, nxv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST1XL(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
         implicit none
         integer :: nop, idimp, nx, nxv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxv) :: cu
         end subroutine
      end interface
      interface
         subroutine GBPUSH13(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,n&
     &x,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: omx, qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz
         end subroutine
      end interface
      interface
         subroutine GSBPUSH13(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,&
     &nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: omx, qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz
         end subroutine
      end interface
      interface
         subroutine GBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,&
     &nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: omx, qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz
         end subroutine
      end interface
      interface
         subroutine GSBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop&
     &,nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: omx, qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz
         end subroutine
      end interface
      interface
         subroutine RETARD1(part,dtc,idimp,nop,nx,ipbc)
         implicit none
         integer :: idimp, nop, nx, ipbc
         real :: dtc
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface djpost
         module procedure igjpost1
      end interface
!
      interface push3
         module procedure igbpush13
      end interface
!
      interface retard
         module procedure iretard1
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igjpost1(part,cu,nop,qm,dt,tdjpost,nx,ipbc,inorder,d&
     &jopt)
! deposit current
         implicit none
         integer :: nop, nx, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: cu
! local data
         integer :: idimp, nxv, order, opt, ltime
         real :: tj
         idimp = size(part,1)
         nxv = size(cu,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,ltime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
               else if (opt==VECTOR) then
                  call GSJPOST1XL(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
               else
                  call GJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSJPOST1(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
               else if (opt==VECTOR) then
                  call GSJPOST1X(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
               else
                  call GJPOST1(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine igjpost1
!
         subroutine igbpush13(part,fxyz,byz,omx,nop,qbm,dt,dtc,ek,tpush,&
     &nx,ipbc,inorder,popt)
! push particles with 1-2/2d electromagnetic fields
         implicit none
         integer :: nop, nx, ipbc
         integer, optional :: inorder, popt
         real :: omx, qbm, dt, dtc, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: fxyz, byz
! local data
         integer :: idimp, nxv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(fxyz,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(byz,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,&
     &nop,nx,nxv,ipbc)
               else
                  call GBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,n&
     &op,nx,nxv,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSBPUSH13(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,n&
     &op,nx,nxv,ipbc)
               else
                  call GBPUSH13(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,no&
     &p,nx,nxv,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igbpush13
!
         subroutine iretard1(part,nop,dtc,nx,ipbc)
! retards particle positions half time-step
         implicit none
         integer :: nop, nx
         integer, optional :: ipbc
         real :: dtc
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, lpbc
         idimp = size(part,1)
         lpbc = 1
         if (present(ipbc)) lpbc = ipbc
         call RETARD1(part,dtc,idimp,nop,nx,ipbc)
         end subroutine iretard1
!
      end module bpush1d
