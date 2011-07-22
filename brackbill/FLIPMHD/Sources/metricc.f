      subroutine metricc(ncells,ijkcell,   &
          vol,   &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz)
!
      use vast_kind_param, ONLY : double
      use geometry_com_M
      implicit none
!
      integer :: ncells,ijkcell(*)
      real(double) ::    &
          vol(*),   &
          tsix(*),tsiy(*),tsiz(*),   &
          etax(*),etay(*),etaz(*),   &
          nux(*),nuy(*),nuz(*)
!
      real(double) :: rvol
      integer :: n,ijk
!
      do 1 n=1,ncells
!
      ijk=ijkcell(n)
!
      rvol=1./vol(ijk)
!
      tsix(ijk)=(c1x(ijk)+c2x(ijk)+c5x(ijk)+c6x(ijk))*rvol
      tsiy(ijk)=(c1y(ijk)+c2y(ijk)+c5y(ijk)+c6y(ijk))*rvol
      tsiz(ijk)=(c1z(ijk)+c2z(ijk)+c5z(ijk)+c6z(ijk))*rvol
!
      etax(ijk)=(c2x(ijk)+c3x(ijk)+c6x(ijk)+c7x(ijk))*rvol
      etay(ijk)=(c2y(ijk)+c3y(ijk)+c6y(ijk)+c7y(ijk))*rvol
      etaz(ijk)=(c2z(ijk)+c3z(ijk)+c6z(ijk)+c7z(ijk))*rvol
!
      nux(ijk)=(c5x(ijk)+c6x(ijk)+c7x(ijk)+c8x(ijk))*rvol
      nuy(ijk)=(c5y(ijk)+c6y(ijk)+c7y(ijk)+c8y(ijk))*rvol
      nuz(ijk)=(c5z(ijk)+c6z(ijk)+c7z(ijk)+c8z(ijk))*rvol
!
    1 continue
!
      return
      end subroutine metricc
