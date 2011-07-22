      subroutine stress(ncells,ijkcell,     &
          pixx,pixy,pixz,piyy,piyz,pizz,    &
          exx,exy,exz,eyy,eyz,ezz,mu,lambda)
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      integer :: ijkcell(*)
      real(double) ::     &
          pixx(*),pixy(*),pixz(*),piyy(*),piyz(*),pizz(*),     &
          exx(*),exy(*),exz(*),eyy(*),eyz(*),ezz(*)
!
      real(double) :: mu,lambda
      integer :: n, ijk
!
      do 100 n=1,ncells
      ijk=ijkcell(n)
      pixx(ijk)=2.0*mu*exx(ijk)
      pixy(ijk)=2.0*mu*exy(ijk)
      pixz(ijk)=2.0*mu*exz(ijk)
      piyy(ijk)=2.0*mu*eyy(ijk)
      piyz(ijk)=2.0*mu*eyz(ijk)
      pizz(ijk)=2.0*mu*ezz(ijk)
  100 continue
!
      do 400 n=1,ncells
      ijk=ijkcell(n)
      divu=exx(ijk)+eyy(ijk)+ezz(ijk)
      pixx(ijk)=pixx(ijk)+lambda*divu
      piyy(ijk)=piyy(ijk)+lambda*divu
      pizz(ijk)=pizz(ijk)+lambda*divu
  400 continue
!
      return
      end subroutine stress
