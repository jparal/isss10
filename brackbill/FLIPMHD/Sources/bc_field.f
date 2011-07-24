      subroutine bc_field(nxp,nyp,nzp,   &
          x,y,z)
!
      use vast_kind_param, ONLY : double
      implicit none
      integer, intent(in) :: nxp, nyp, nzp
      real(double) :: x(nxp,nyp,*),y(nxp,nyp,*),z(nxp,nyp,*)
      real(double) :: zero
      integer :: i, j, k
!
!     a routine to impose double periodicity
!     for toroidal geometry
!     sets ghost cell values only
!     for cell-centered, scalar variables
!
!     called by resistive_diff
!     sets ghost cell values of B to zero
!
      zero=0.0d0
!
      do 10 k=2,nzp-1
!
      do 1 i=2,nxp-1
!
      x(i,nyp,k)=zero
      y(i,nyp,k)=zero
      z(i,nyp,k)=zero
      x(i,1,k)=zero
      y(i,1,k)=zero
      z(i,1,k)=zero
!
    1 continue
!
!
!     periodicity in the poloidal angle
!
!
      do 2 j=1,nyp
!
      x(nxp,j,k)=zero
      y(nxp,j,k)=zero
      z(nxp,j,k)=zero
      x(1,j,k)=zero
      y(1,j,k)=zero
      z(1,j,k)=zero
!
    2 continue
!
   10 continue
!
      return
      end subroutine bc_field
