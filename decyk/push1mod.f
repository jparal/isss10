!-----------------------------------------------------------------------
!
      module push1d
!
! Fortran90 interface to 1d PIC Fortran77 library push1lib.f
! push1mod.f contains interface procedures to process particles:
!            defines module push1d
! dpost => igpost1 deposits charge density, with various interpolations
!          and optimizations.
!          calls GPOST1, GSPOST1, GSPOST1X, GPOST1L, GSPOST1L, or
!          GSPOST1XL
! push => igpush1 pushes particles, with various interpolations and
!         optimizations.
!         calls GPUSH1, GSPUSH1, GPUSH1L, or GSPUSH1L
! sortp => isortp1x sort particles by grid with various interpolations.
!          calls SORTP1X, or SORTP1XL
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
      public :: dpost, push, sortp, dpostgl, pushgl
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GPOST1(part,q,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST1(part,q,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST1X(part,q,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         end subroutine
      end interface
      interface
         subroutine GPOST1L(part,q,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST1L(part,q,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST1XL(part,q,qm,nop,idimp,nxv)
         implicit none
         integer :: nop, idimp, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         end subroutine
      end interface
      interface
         subroutine GPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         end subroutine
      end interface
      interface
         subroutine GSPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         end subroutine
      end interface
      interface
         subroutine GPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         end subroutine
      end interface
      interface
         subroutine GSPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
         implicit none
         integer :: idimp, nop, nx, nxv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         end subroutine
      end interface
      interface
         subroutine SORTP1X(part,pt,ip,npic,idimp,nop,nx1)
         implicit none
         integer :: idimp, nop, nx1
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(nx1) :: npic
         end subroutine
      end interface
      interface
         subroutine SORTP1XL(part,pt,ip,npic,idimp,nop,nx1)
         implicit none
         integer :: idimp, nop, nx1
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(nx1) :: npic
         end subroutine
      end interface
      interface
         subroutine DPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxh,nxvh)
         implicit none
         integer :: nop, idimp, nx, nxh, nxvh
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh) :: q
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine PUSH1GL(part,fx,sctx,qbm,dt,ek,idimp,nop,nx,nxh,nxvh&
     &)
         implicit none
         integer :: nop, idimp, nx, nxh, nxvh
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh) :: fx
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dpost
         module procedure igpost1
      end interface
!
      interface push
         module procedure igpush1
      end interface
!
      interface sortp
         module procedure isortp1x
      end interface
!
      interface dpostgl
         module procedure idpost1gl
      end interface
!
      interface pushgl
         module procedure ipush1gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igpost1(part,q,nop,qm,tdpost,inorder,dopt)
! deposit charge
         implicit none
         integer :: nop
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: q
! local data
         integer :: idimp, nxv, order, opt, ltime
         real :: td
         idimp = size(part,1)
         nxv = size(q)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! initialize timer
         call wtimer(td,ltime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call GSPOST1L(part,q,qm,nop,idimp,nxv)
            else if (opt==VECTOR) then
               call GSPOST1XL(part,q,qm,nop,idimp,nxv)
            else
               call GPOST1L(part,q,qm,nop,idimp,nxv)
            endif
         else
            if (opt==LOOKAHEAD) then
               call GSPOST1(part,q,qm,nop,idimp,nxv)
            else if (opt==VECTOR) then
               call GSPOST1X(part,q,qm,nop,idimp,nxv)
            else
               call GPOST1(part,q,qm,nop,idimp,nxv)
            endif
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine igpost1
!
         subroutine igpush1(part,fx,nop,qbm,dt,ek,tpush,nx,ipbc,inorder,&
     &popt)
! push particles with 1d electrostatic fields
         implicit none
         integer :: nop, nx, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: fx
! local data
         integer :: idimp, nxv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(fx)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call GSPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
            else
               call GPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call GSPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
            else
               call GPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
            endif
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igpush1
!
         subroutine isortp1x(part,pt,ip,nop,npic,tsort,inorder)
! sort particles by x grid using memory-conserving bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: pt
         integer, dimension(:), pointer :: ip
         integer, dimension(:), pointer :: npic
! local data
         integer :: idimp, nx1, order, ltime
         real :: ts
         order = QUADRATIC
         if (present(inorder)) order = inorder
         idimp = size(part,1); nx1 = size(npic)
! initialize timer
         call wtimer(ts,ltime,-1)
         if (order==LINEAR) then
            call SORTP1XL(part,pt,ip,npic,idimp,nop,nx1)
         else
            call SORTP1X(part,pt,ip,npic,idimp,nop,nx1)
         endif
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine isortp1x
!
         subroutine idpost1gl(part,q,nop,qm,nx,nxh,tdpost)
! deposit charge using gridless method
         implicit none
         integer :: nop, nx, nxh
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: q
! local data
         integer :: idimp, nxvh, ltime
         real :: td
         complex, dimension(size(q,1)/2) :: sctx
         idimp = size(part,1)
         nxvh = size(q,1)/2
! initialize timer
         call wtimer(td,ltime,-1)
         call DPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxh,nxvh)
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine idpost1gl
!
         subroutine ipush1gl(part,fx,nop,qbm,dt,ek,nx,nxh,tpush)
! push particles with 1d electrostatic fields using gridless method
         implicit none
         integer :: nop, nx, nxh
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: fx
! local data
         integer :: idimp, nxvh, ltime
         real :: tp
         complex, dimension(size(fx,1)/2) :: sctx
         idimp = size(part,1)
         nxvh = size(fx,1)/2
! initialize timer
         call wtimer(tp,ltime,-1)
         call PUSH1GL(part,fx,sctx,qbm,dt,ek,idimp,nop,nx,nxh,nxvh)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ipush1gl
!
      end module push1d
