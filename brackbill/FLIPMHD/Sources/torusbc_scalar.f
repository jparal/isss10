      subroutine torusbc_scalar(nxp,nyp,nzp,    &
          x)
!
!     a routine to impose double periodicity
!     for toroidal geometry for vertex variables
!
!     called by VINIT_GMRES
!
      use vast_kind_param, ONLY :  double
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z
      implicit real*8 (a-h,o-z)
!
      real(double) ::  x(nxp,nyp,*)
!
      do 10 k=2,nzp
!
!     periodicity in the toroidal angle
!
      if(periodic_y) then
!
      do 1 i=2,nxp
!
      x(i,2,k)=x(i,nyp,k)+x(i,2,k)
!
      x(i,nyp,k)=x(i,2,k)
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
      x(2,j,k)=x(nxp,j,k)+x(2,j,k)
!
      x(nxp,j,k)=x(2,j,k)
!
    2 continue
!
      endif
!
   10 continue
!
!
      return
      end subroutine torusbc_scalar
