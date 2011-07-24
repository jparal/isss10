!-----------------------------------------------------------------------
!
      module fft1d
!
! fft1mod.f contains interface procedures to perform ffts:
!           defines module fft1d
! fft_init => ifft1rxinit initializes real to complex fft tables.
!             calls FFT1RX
! fft => ifft1rx performs 1d real to complex fft and its inverse.
!        calls FFT1RX
! fft => ifft1r2 performs multiple 1d real to complex fft and its
!        inverse for scalar or 2 and 3 component vector arrays.
!        calls FFT1RX, FFT1R2, or FFT1R3
! fst_init => ifst1rxinit initializes 1d real sine or cosine transform.
!             tables.
!             calls FST1RX
! fst => ifst1rx performs 1d scalar real sine transform.
!        calls FST1RX
! fst => ifst1r2 performs 1d scalar real sine transform, for scalar
!        or 2 component vector arrays.
!        calls FST1R2
! fct => ifct1rx performs 1d scalar real cosine transform.
!        calls FCT1RX
! fct => ifct1r2 performs 1d scalar real cosine transform, for scalar
!        or 2 component vector arrays.
!        calls FCT1R2
! fcst => ifcst1r3 performs 1d scalar real sine and cosine transform,
!         3 component vector arrays.
!         calls FCST1R3
! Fortran90 interface to 1d PIC Fortran77 library fft1lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: july 12, 2010
!
      use globals, only: LINEAR, QUADRATIC
      use diag1d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: fft_init, fft, fst_init, fst, fct, fcst
!
!
! define interface to original Fortran77 procedures
      interface
         subroutine FFT1RX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(nxhd) :: t, sct
         end subroutine
      end interface
      interface
         subroutine FFT1R2(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(nxhd) :: sct
         complex, dimension(2,nxhd) :: t
         end subroutine
      end interface
      interface
         subroutine FFT1R3(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(nxhd) :: sct
         complex, dimension(3,nxhd) :: t
         end subroutine
      end interface
      interface
         subroutine FST1RX(f,t,isign,mixup,sctd,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(nxhd) :: t
         complex, dimension(nxd) :: sctd
         end subroutine
      end interface
      interface
         subroutine FCT1RX(f,t,isign,mixup,sctd,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(nxhd) :: t
         complex, dimension(nxd) :: sctd
         end subroutine
      end interface
      interface
         subroutine FST1R2(f,t,isign,mixup,sctd,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(2,nxhd) :: t
         complex, dimension(nxd) :: sctd
         end subroutine
      end interface
      interface
         subroutine FCT1R2(f,t,isign,mixup,sctd,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(2,nxhd) :: t
         complex, dimension(nxd) :: sctd
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft_init
         module procedure ifft1rxinit
      end interface
!
      interface fft
         module procedure ifft1rx
         module procedure ifft1r2
      end interface
!
      interface fst_init
         module procedure ifst1rxinit
      end interface
!
      interface fst
         module procedure ifst1rx
         module procedure ifst1r2
      end interface
!
      interface fct
         module procedure ifct1rx
         module procedure ifct1r2
      end interface
!
      interface fcst
         module procedure ifcst1r3
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ifft1rxinit(mixup,sct,indx)
! initialize 1d real to complex fft
         implicit none
         integer :: indx
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         integer :: isign = 0, nxd = 1, nxhd
         real :: f
         complex, dimension(1) :: t
         nxhd = size(mixup)
         call FFT1RX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         end subroutine ifft1rxinit
!
         subroutine ifft1rx(f,isign,mixup,sct,tfft,indx,order)
! perform 1d scalar real to complex fft
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(mixup)) :: t
         nxd = size(f); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call FFT1RX(f(1),t,isign,mixup,sct,indx,nxd,nxhd)
         else
            call FFT1RX(f(2),t,isign,mixup,sct,indx,nxd,nxhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifft1rx
!
         subroutine ifft1r2(f,isign,mixup,sct,tfft,indx,order)
! perform 1d vector real to complex fft for 2 component vectors
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(f,1),size(mixup)) :: t
         nxd = size(f,2); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (1)
            if (inorder==LINEAR) then
               call FFT1RX(f(1,1),t,isign,mixup,sct,indx,nxd,nxhd)
            else
               call FFT1RX(f(1,2),t,isign,mixup,sct,indx,nxd,nxhd)
            endif
         case (2)
            if (inorder==LINEAR) then
               call FFT1R2(f(1,1),t,isign,mixup,sct,indx,nxd,nxhd)
            else
               call FFT1R2(f(1,2),t,isign,mixup,sct,indx,nxd,nxhd)
            endif
         case (3)
            if (inorder==LINEAR) then
               call FFT1R3(f(1,1),t,isign,mixup,sct,indx,nxd,nxhd)
            else
               call FFT1R3(f(1,2),t,isign,mixup,sct,indx,nxd,nxhd)
            endif
         end select
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifft1r2
!
         subroutine ifst1rxinit(mixup,sctd,indx)
! initialize 1d real sine or cosine transforms
         implicit none
         integer :: indx
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
         integer :: isign = 0, nxd, nxhd
         real :: f
         complex, dimension(1) :: t
         nxd = size(sctd); nxhd = size(mixup)
         call FST1RX(f,t,isign,mixup,sctd,indx,nxd,nxhd)
         end subroutine ifst1rxinit
!
         subroutine ifst1rx(f,isign,mixup,sctd,tfft,indx,order)
! perform 1d scalar real sine transform
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(mixup)) :: t
         nxd = size(f); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call FST1RX(f(1),t,isign,mixup,sctd,indx,nxd,nxhd)
         else
            call FST1RX(f(2),t,isign,mixup,sctd,indx,nxd,nxhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifst1rx
!
         subroutine ifct1rx(f,isign,mixup,sctd,tfft,indx,order)
! perform 1d scalar real cosine transform
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(mixup)) :: t
         nxd = size(f); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call FCT1RX(f(1),t,isign,mixup,sctd,indx,nxd,nxhd)
         else
            call FCT1RX(f(2),t,isign,mixup,sctd,indx,nxd,nxhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifct1rx
!
         subroutine ifst1r2(f,isign,mixup,sctd,tfft,indx,order)
! perform 1d scalar real sine transform for 2 component vectors
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(f,1),size(mixup)) :: t
         nxd = size(f); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (1)
            if (inorder==LINEAR) then
               call FST1RX(f(1,1),t,isign,mixup,sctd,indx,nxd,nxhd)
            else
               call FST1RX(f(1,2),t,isign,mixup,sctd,indx,nxd,nxhd)
            endif
         case (2)
            if (inorder==LINEAR) then
               call FST1R2(f(1,1),t,isign,mixup,sctd,indx,nxd,nxhd)
            else
               call FST1R2(f(1,2),t,isign,mixup,sctd,indx,nxd,nxhd)
            endif
         end select
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifst1r2
!
         subroutine ifct1r2(f,isign,mixup,sctd,tfft,indx,order)
! perform 1d scalar real cosine transform for 2 component vectors
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(f,1),size(mixup)) :: t
         nxd = size(f); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (1)
            if (inorder==LINEAR) then
               call FCT1RX(f(1,1),t,isign,mixup,sctd,indx,nxd,nxhd)
            else
               call FCT1RX(f(1,2),t,isign,mixup,sctd,indx,nxd,nxhd)
            endif
         case (2)
            if (inorder==LINEAR) then
               call FCT1R2(f(1,1),t,isign,mixup,sctd,indx,nxd,nxhd)
            else
               call FCT1R2(f(1,2),t,isign,mixup,sctd,indx,nxd,nxhd)
            endif
         end select
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifct1r2
!
         subroutine ifcst1r3(f,isign,mixup,sctd,tfft,indx,order)
! perform 1d scalar real cosine and sine transform for 3 component
! vectors
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(f,1),size(mixup)) :: t
         nxd = size(f); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (3)
            if (inorder==LINEAR) then
               call FCST1R3(f(1,1),t,isign,mixup,sctd,indx,nxd,nxhd)
            else
               call FCST1R3(f(1,2),t,isign,mixup,sctd,indx,nxd,nxhd)
            endif
         end select
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifcst1r3
!
      end module fft1d
