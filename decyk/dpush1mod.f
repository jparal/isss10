!-----------------------------------------------------------------------
!
      module dpush1d
!
! Fortran90 interface to 1d PIC Fortran77 library dpush1lib.f
! dpush1mod.f contains interface procedures to process particles with
!             darwin electric and magnetic fields:
!             defines module dpush1d
! dmjpost => igmjpost1 deposits momentum flux, with various
!            interpolations and optimizations.
!            calls GMJPOST1, GSMJPOST1, GMJPOST1L, or GSMJPOST1L
! dcjpost => igdcjpost1 deposits momentum flux, acceleration density,
!            and current density, with various interpolations and
!            optimizations.
!            calls GDCJPOST1, GSDCJPOST1, GDCJPOST1L, or GSDCJPOST1L
! written by viktor k. decyk, ucla
! copyright 2006, regents of the university of california
! update: july 12, 2011
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use diag1d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, dmjpost, dcjpost
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GMJPOST1(part,amu,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSMJPOST1(part,amu,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST1(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,idi&
     &mp,nop,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: omx, qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz, cu, dcu, amu
         end subroutine
      end interface
      interface
         subroutine GSDCJPOST1(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,id&
     &imp,nop,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: omx, qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz, cu, dcu, amu
         end subroutine
      end interface
      interface
         subroutine GMJPOST1L(part,amu,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSMJPOST1L(part,amu,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,id&
     &imp,nop,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: omx, qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz, cu, dcu, amu
         end subroutine
      end interface
      interface
         subroutine GSDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,i&
     &dimp,nop,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: omx, qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv) :: fxyz
         real, dimension(2,nxv) :: byz, cu, dcu, amu
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dmjpost
         module procedure igmjpost1
      end interface
!
      interface dcjpost
         module procedure igdcjpost1
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igmjpost1(part,amu,nop,qm,tdcjpost,inorder,djopt)
! deposit momentum flux with 1-2/2d electromagnetic fields
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: amu
! local data
         integer :: idimp, nxv, order, opt, ltime
         real :: tdc
         idimp = size(part,1)
         nxv = size(amu,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(amu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSMJPOST1L(part,amu,qm,nop,idimp,nxv)
               else
                  call GMJPOST1L(part,amu,qm,nop,idimp,nxv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSMJPOST1(part,amu,qm,nop,idimp,nxv)
               else
                  call GMJPOST1(part,amu,qm,nop,idimp,nxv)
               endif
            endif
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igmjpost1
!
         subroutine igdcjpost1(part,fxyz,byz,cu,dcu,amu,omx,nop,qm,qbm,d&
     &t,tdcjpost,inorder,djopt)
! deposit momentum flux, acceleration density, and current density
! with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: omx, qm, qbm, dt, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: fxyz, byz, cu, dcu, amu
! local data
         integer :: idimp, nxv, order, opt, ltime
         real :: tdc
         idimp = size(part,1)
         nxv = size(fxyz,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(byz,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,d&
     &t,idimp,nop,nxv)
               else
                  call GDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt&
     &,idimp,nop,nxv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSDCJPOST1(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt&
     &,idimp,nop,nxv)
               else
                  call GDCJPOST1(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,&
     &idimp,nop,nxv)
               endif
            endif
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igdcjpost1
!
      end module dpush1d
