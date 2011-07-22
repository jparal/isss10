      subroutine bc_vertex(nxp,nyp,nzp,    &
          x)
!
!     a routine to impose double periodicity
!     for toroidal geometry
!     sets ghost cell values only
!     for vertex-centered, scalar variables
!     either periodic or Dirichlet bc's (x=0)
      use vast_kind_param, ONLY :  double
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z
      implicit none
!----------------------------------------------------
!     D u m m y  A r g u m e n t s
!----------------------------------------------------
      integer , intent(in) :: nxp, nyp, nzp
      real(double) ::  x(nxp,nyp,*)
      integer :: i,j,k
!
!     called by DIVB_PROJECTION, POISSON_CG, VINIT_GMRES
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
      x(i,nyp,k)=x(i,2,k)
      x(i,1,k)=x(i,nyp-1,k)
!
    1 continue
!
      else
!
      do i=2,nxp
!
      x(i,2,k)=0.0d0
      x(i,nyp,k)=0.0d0
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
      x(nxp,j,k)=x(2,j,k)
      x(1,j,k)=x(nxp-1,j,k)
!
    2 continue
!
      else
!
      do j=2,nyp
!
      x(2,j,k)=0.0d0
      x(nxp,j,k)=0.0d0
!
      enddo
      endif
   10 continue
!
      if (periodic_z) then
!
      do i=2,nxp
      do j=2,nyp
!
      x(i,j,nzp)=x(i,j,2)
      x(i,j,1)=x(i,j,nzp-1)
!
      enddo
      enddo
!
      else
!
      do i=2,nxp
      do j=2,nyp
!
      x(i,j,2)=0.0d0
      x(i,j,nzp)=0.0d0
!
      enddo
      enddo
!
      endif
!
      return
      end subroutine bc_vertex
