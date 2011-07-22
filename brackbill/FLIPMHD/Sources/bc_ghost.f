      subroutine bc_ghost(nxp,nyp,nzp,iwid,jwid,kwid,   &
          bx,by,bz)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!DIR$ FREE

      use vast_kind_param, ONLY : double
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nxp, nyp, nzp, iwid, jwid, kwid
      real(double) :: bx(*), by(*), bz(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, ijk
!-----------------------------------------------

!
!     a routine to set ghost cell values of cell-centered
!     vector to zero
!
!     called by VINIT_GMRES
!
!
!     bottom and top ( k=1 and k=nzp)
!
      do i=1,nxp
      do j=1,nyp
!
      ijk=(j-1)*jwid+(i-1)*iwid+1
!
      bx(ijk)=0.0
      by(ijk)=0.0
      bz(ijk)=0.0
!
      ijk=(nzp-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
!
      bx(ijk)=0.0
      by(ijk)=0.0
      bz(ijk)=0.0
!
      enddo
      enddo
!
!
!     right and left ( i=1 and i=nxp)
!
      do j=1,nyp
      do k=1,nzp
!
      ijk=(k-1)*kwid+(j-1)*jwid+1
!
      bx(ijk)=0.0
      by(ijk)=0.0
      bz(ijk)=0.0
!
      ijk=(k-1)*kwid+(j-1)*jwid+(nxp-1)*iwid+1
!
      bx(ijk)=0.0
      by(ijk)=0.0
      bz(ijk)=0.0
!
      enddo
      enddo
!
!     front and back ( j=1 and j=nyp)
!
      do k=1,nzp
      do i=1,nxp
!
      ijk=(k-1)*kwid+(i-1)*iwid+1
!
      bx(ijk)=0.0
      by(ijk)=0.0
      bz(ijk)=0.0
!
      ijk=(k-1)*kwid+(nyp-1)*jwid+(i-1)*iwid+1
!
      bx(ijk)=0.0
      by(ijk)=0.0
      bz(ijk)=0.0
!
      enddo
      enddo
!
      return
      end subroutine bc_ghost
