       subroutine eos(p,ro,vol,sie,csq,gm1,ncells,ijkcell)
!
      use vast_kind_param, ONLY : double
      implicit none
!
      integer, intent(in) :: ncells, ijkcell(*)
      real(double) ::  p(*),ro(*),vol(*),sie(*),csq(*)
      real(double), intent(in) :: gm1
!
      integer n, ijk
      real(double) :: zero
      zero=0.0
!
      do 100 n=1,ncells
      ijk=ijkcell(n)
      p(ijk)=gm1*ro(ijk)*sie(ijk)
      p(ijk)=max(zero,p(ijk))
      csq(ijk)=gm1*(gm1+1.0)*sie(ijk)
  100 continue
      return
      end subroutine eos
