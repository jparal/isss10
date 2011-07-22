      subroutine bc_current(nxp,nyp,nzp,    &
          strait,                           &
          x,y,z)
!
      use vast_kind_param, ONLY : double
      use cophys_com_M, ONLY : cdlt,sdlt,dz
      use cindex_com_M, ONLY : iwid, jwid, kwid
      use blcom_com_M, ONLY : periodic_x, periodic_y, periodic_z
!
      implicit none
      integer, intent(in) :: nxp, nyp, nzp
      real(double) :: x(*),y(*),z(*)
      real(double) :: strait
!
      integer i, j, k, ijkl, ijkr
!
!     a routine to impose double periodicity
!     for toroidal geometry
!     for vertex vector, assembles final result
!
!     called by RESISTIVE_DIFF 
!
!
      do 10 k=2,nzp
!
!     periodicity in the toroidal angle
!
      if(periodic_y) then
!
      do 1 i=2,nxp
!
      ijkl=(i-1)*iwid+jwid+(k-1)*kwid+1
      ijkr=(i-1)*iwid+(nyp-1)*jwid+(k-1)*kwid+1
      x(ijkl)=cdlt*x(ijkr)+sdlt*y(ijkr)+x(ijkl)
      y(ijkl)=-sdlt*x(ijkr)+cdlt*y(ijkr)+y(ijkl)    &
          -strait*dz
      z(ijkl)=z(ijkr)+z(ijkl)
!
      x(ijkr)=cdlt*x(ijkl)-sdlt*y(ijkl)
      y(ijkr)=sdlt*x(ijkl)+cdlt*y(ijkl)    &
          +strait*dz
      z(ijkr)=z(ijkl)
!
    1 continue
!
      else
!
      do i=2,nxp
!
      ijkl=(i-1)*iwid+jwid+(k-1)*kwid+1
      ijkr=(i-1)*iwid+(nyp-1)*jwid+(k-1)*kwid+1
      x(ijkl)=0.
      y(ijkl)=0.
      z(ijkl)=0.
!
      x(ijkr)=0.
      y(ijkr)=0.
      z(ijkr)=0.
!
      enddo
!
      endif
!
!     periodicity in the poloidal angle
!
      if(periodic_x) then
!
      do 2 j=2,nyp
!
      ijkl=iwid+(j-1)*jwid+(k-1)*kwid+1
      ijkr=(nxp-1)*iwid+(j-1)*jwid+(k-1)*kwid+1
      
      x(ijkl)=x(ijkr)+x(ijkl)
      y(ijkl)=y(ijkr)+y(ijkl)
      z(ijkl)=z(ijkr)+z(ijkl)
!
      x(ijkr)=x(ijkl)
      y(ijkr)=y(ijkl)
      z(ijkr)=z(ijkl)
!
    2 continue
!
      endif
!
   10 continue
!
      do i=2,nxp
      do j=2,nyp
!
      ijkl=(i-1)*iwid+(j-1)*jwid+kwid+1
      ijkr=(i-1)*iwid+(j-1)*jwid+(nzp-1)*kwid+1
      x(ijkl)=0.0
      y(ijkl)=0.0
      z(ijkl)=0.0
!
      x(ijkr)=0.0
      y(ijkr)=0.0
      z(ijkr)=0.0
!
      enddo
      enddo
!
      return
      end subroutine bc_current
