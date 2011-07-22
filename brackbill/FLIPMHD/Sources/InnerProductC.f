      subroutine InnerProductC(Ax,Ay, Az,Bx, By, Bz, WeightC, AdotB)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double
      use blcom_com_M, ONLY: vol, ijkcell
      use cindex_com_M, ONLY: ncells
!
!-----------------------------------------------
!   D u m m y   V a r i a b l e s
!-----------------------------------------------
      real(double), intent(in)  :: Ax(*), Ay(*), Az(*), Bx(*), By(*), Bz(*)
      real(double), intent(in) :: WeightC
      real(double), intent(out) :: AdotB
!
      AdotB=0.0d0
!
      do n=1,ncells
        ijk=ijkcell(n)
        AdotB=AdotB+(Ax(ijk)*Bx(ijk)+Ay(ijk)*By(ijk)+Az(ijk)*Bz(ijk))*vol(ijk)
      enddo
!
      AdotB=AdotB*WeightC
!
      return
!
      end subroutine InnerProductC

