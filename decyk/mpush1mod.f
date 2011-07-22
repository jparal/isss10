!-----------------------------------------------------------------------
!
      module mpush1d
!
! Fortran90 interface to 1d PIC Fortran77 library mpush1lib.f
! mpush1mod.f contains multi-tasking interface procedures to process
!             particles:
!             defines module mpush1d
! dpost => imgpost1 deposits charge density, with various interpolations
!          and optimizations.
!          calls MGPOST1, MGSPOST1, MGSPOST1X, MGPOST1L, MGSPOST1L, or
!          MGSPOST1XL
! push => imgpush1 push particles, with various interpolations and
!         optimizations.
!         calls MGPUSH1, MGSPUSH1, MGPUSH1L, or MGSPUSH1L
! sortp => imsortp1x sorts particles by x grid using memory-conserving
!          bin sort, with various interpolations.
!          calls MSORTP1X, or MSORTP1XL
! sortp => imdsortp1x sorts particles by x grid using optimized bin sort
!          with various interpolations.
!          calls MDSORTP1X, or MDSORTP1XL
! dpostgl => imdpost1gl deposits charge density using gridless method.
!            calls MDPOST1GL
! pushgl => impush1gl push particles using optimzed gridless method.
!           calls MPUSH1GL
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: january 8, 2010
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use push1d, only: wtimer
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: dpost, push, sortp, dpostgl, pushgl
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MGPOST1(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         real, dimension(nxv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST1(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         real, dimension(nxv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST1X(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr&
     &)
         implicit none
         integer :: nop, idimp, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         real, dimension(nxv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGPOST1L(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         real, dimension(nxv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST1L(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr&
     &)
         implicit none
         integer :: nop, idimp, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         real, dimension(nxv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST1XL(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ier&
     &r)
         implicit none
         integer :: nop, idimp, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: q
         real, dimension(nxv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idta&
     &sk,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idta&
     &sk,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idt&
     &ask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv) :: fx
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSORTP1X(part,pt,ip,npic,idimp,nop,nx1,npicp,idtask,&
     &nmt,ierr)
         implicit none
         integer :: idimp, nop, nx1, nmt, ierr
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(nx1) :: npic
         integer, dimension(nx1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSORTP1XL(part,pt,ip,npic,idimp,nop,nx1,npicp,idtask&
     &,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx1, nmt, ierr
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(nx1) :: npic
         integer, dimension(nx1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDSORTP1X(parta,partb,npic,idimp,nop,nx1,npicp,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx1, nmt, ierr
         real, dimension(idimp,nop) :: parta, partb
         integer, dimension(nx1) :: npic
         integer, dimension(nx1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDSORTP1XL(parta,partb,npic,idimp,nop,nx1,npicp,idta&
     &sk,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx1, nmt, ierr
         real, dimension(idimp,nop) :: parta, partb
         integer, dimension(nx1) :: npic
         integer, dimension(nx1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxh,nxvh,qp,sc&
     &txp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, nxh, nxvh, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh) :: q
         complex, dimension(nxvh) :: sctx
         real, dimension(2*nxvh,nmt) :: qp
         complex, dimension(nxvh,nmt) :: sctxp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPUSH1GL(part,fx,sctx,qbm,dt,ek,idimp,nop,nx,nxh,nxv&
     &h,sctxp,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, nxh, nxvh, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh) :: fx
         complex, dimension(nxvh) :: sctx
         complex, dimension(nxvh,nmt) :: sctxp
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dpost
         module procedure imgpost1
      end interface
!
      interface dpostgl
         module procedure imdpost1gl
      end interface
!
      interface push
         module procedure imgpush1
      end interface
!
      interface pushgl
         module procedure impush1gl
      end interface
!
      interface sortp
         module procedure imsortp1x
         module procedure imdsortp1x
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imgpost1(part,q,nop,qm,tdpost,inorder,dopt)
! multi-tasking charge deposit
         implicit none
         integer :: nop
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: q
! local data
         integer :: idimp, nxv, nmt, ltime, order, opt, ierr
         real :: td
!        real, dimension(size(q,1),ntasks) :: qp
         real, dimension(:,:), allocatable, save :: qp
         integer, save :: szbuf = 0
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(q,1)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! check if size of buffer has changed
         if (szbuf < nxv) then
            if (szbuf /= 0) deallocate(qp)
! allocate buffer
            allocate(qp(nxv,ntasks))
            szbuf = nxv
         endif
! initialize timer
         call wtimer(td,ltime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MGSPOST1L(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr&
     &)
            else if (opt==VECTOR) then
               call MGSPOST1XL(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ier&
     &r)
            else
               call MGPOST1L(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MGSPOST1(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MGSPOST1X(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr&
     &)
            else
               call MGPOST1(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine imgpost1
!
         subroutine imgpush1(part,fx,nop,qbm,dt,ek,tpush,nx,inorder,popt&
     &)
! multi-tasking particle push with 1d electrostatic fields
         implicit none
         integer :: nop, nx
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: fx
! local data
         integer :: idimp, nxv, nmt, ltime, order, opt, ierr
         real :: tp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(fx,1)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MGSPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idt&
     &ask,nmt,ierr)
            else
               call MGPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idta&
     &sk,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MGSPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idta&
     &sk,nmt,ierr)
            else
               call MGPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idtas&
     &k,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgpush1
!
         subroutine imsortp1x(part,pt,ip,nop,npic,tsort,inorder)
! multi-tasking particle sort by x grid using memory-conserving bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: pt
         integer, dimension(:), pointer :: ip, npic
! local data
         integer, dimension(size(npic),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         integer :: idimp, nx1, nmt, ltime, order, ierr
         real :: ts
         idimp = size(part,1)
         nx1 = size(npic)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if (order==LINEAR) then
            call MSORTP1XL(part,pt,ip,npic,idimp,nop,nx1,npicp,idtask,nm&
     &t,ierr)
         else
            call MSORTP1X(part,pt,ip,npic,idimp,nop,nx1,npicp,idtask,nmt&
     &,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine imsortp1x
!
         subroutine imdsortp1x(parta,partb,nop,npic,tsort,inorder)
! multi-tasking particle sort by x grid using optimized bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: parta, partb
         integer, dimension(:), pointer :: npic
         real, dimension(:,:), pointer :: part
! local data
         integer, dimension(size(npic),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         integer :: idimp, nx1, nmt, ltime, order, ierr
         real :: ts
         idimp = size(parta,1)
         nx1 = size(npic)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if (order==LINEAR) then
            call MDSORTP1XL(parta,partb,npic,idimp,nop,nx1,npicp,idtask,&
     &nmt,ierr)
         else
            call MDSORTP1X(parta,partb,npic,idimp,nop,nx1,npicp,idtask,n&
     &mt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
         part => parta
         parta => partb
         partb => part
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine imdsortp1x
!
         subroutine imdpost1gl(part,q,nop,qm,nx,nxh,tdpost)
! multi-tasking deposit charge using gridless method
         implicit none
         integer :: nop, nx, nxh
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: q
! local data
         integer :: idimp, nxvh, nmt, ltime, ierr
         real :: td
         complex, dimension(size(q,1)/2) :: sctx
         real, dimension(size(q,1),ntasks) :: qp
         complex, dimension(size(q,1)/2,ntasks) :: sctxp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxvh = size(q,1)/2
         nmt = ntasks
! initialize timer
         call wtimer(td,ltime,-1)
         call MDPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxh,nxvh,qp,sctxp,id&
     &task,nmt,ierr)
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine imdpost1gl
!
       subroutine impush1gl(part,fx,nop,qbm,dt,ek,nx,nxh,tpush)
! multi-tasking particle push with 1d electrostatic fields
! using gridless method
         implicit none
         integer :: nop, nx, nxh
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: fx
! local data
         integer :: idimp, nxvh, nmt, ltime, ierr
         real :: tp
         complex, dimension(size(fx,1)/2) :: sctx
         complex, dimension(size(fx,1)/2,ntasks) :: sctxp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxvh = size(fx,1)/2
         nmt = ntasks
! initialize timer
         call wtimer(tp,ltime,-1)
         call MPUSH1GL(part,fx,sctx,qbm,dt,ek,idimp,nop,nx,nxh,nxvh,sctx&
     &p,ekp,idtask,nmt,ierr)
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine impush1gl
!
      end module mpush1d
