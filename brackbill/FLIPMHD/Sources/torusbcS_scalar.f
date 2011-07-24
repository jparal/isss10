      subroutine torusbcS_scalar(nxp,nyp,nzp,    &
          x,is)
!
!     a routine to impose double periodicity
!     for toroidal geometry for vertex variables
!
!     called by VINIT_GMRES
!
      use vast_kind_param, ONLY :  double
      use corgan_com_M, ONLY : itdim
      use cindex_com_M, ONLY : iwid, jwid, kwid
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z
      implicit real*8 (a-h,o-z)
!
      real(double) ::  x(itdim,*)
      integer :: is
      integer :: ijkr,ijkl,ijkb,ijkf
!
      do 10 k=2,nzp
!
!     periodicity in the toroidal angle
!
      if(periodic_y) then
!
      do 1 i=2,nxp
!
      ijkf=(i-1)*iwid+jwid+(k-1)*kwid+1
      ijkb=(i-1)*iwid+(nyp-1)*jwid+(k-1)*kwid+1
      x(ijkf,is)=x(ijkf,is)+x(ijkb,is)
!
      x(ijkb,is)=x(ijkf,is)
!
    1 continue
!
      endif
!
!
!
!     periodicity in the poloidal angle
!
      if(periodic_x) then
!
      do 2 j=2,nyp
!
      ijkl=iwid+(j-1)*jwid+(k-1)*kwid+1
      ijkr=(nxp-1)*iwid+(j-1)*jwid+(k-1)*kwid+1
      x(ijkl,is)=x(ijkl,is)+x(ijkr,is)
!
      x(ijkr,is)=x(ijkl,is)
!
    2 continue
!
      endif
!
   10 continue
!
!
      return
      end subroutine torusbcS_scalar
