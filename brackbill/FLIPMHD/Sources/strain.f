
      subroutine strain(ncells,ijkcell,iwid,jwid,kwid,     &
          exx,exy,exz,eyy,eyz,ezz,      &
          u,v,w,vol)
!
      use vast_kind_param, ONLY : double
      use corgan_com_M, ONLY : itdim
      use Scratch_com_M, ONLY : gradcx, gradcy, gradcz
      use geometry_com_M
      implicit real*8 (a-h,o-z)
!
      integer :: iwid, jwid, kwid,     &
          ijkcell(*)
      real(double) ::      &
          exx(*),exy(*),exz(*),eyy(*),eyz(*),ezz(*),     &
          u(*),v(*),w(*),vol(*)
!
      allocate (gradcx(itdim), gradcy(itdim), gradcz(itdim))
!
      call gradc(ncells,ijkcell,      &
          u,gradcx,gradcy,gradcz)
!
      do 100 n=1,ncells
      ijk=ijkcell(n)
      exx(ijk)=gradcx(ijk)
      exy(ijk)=0.5*gradcy(ijk)
      exz(ijk)=0.5*gradcz(ijk) 
  100 continue
!
      call gradc(ncells,ijkcell,      &
          v,gradcx,gradcy,gradcz)
!
      do 200 n=1,ncells
      ijk=ijkcell(n)
      exy(ijk)=exy(ijk)+0.5*gradcx(ijk)
      eyy(ijk)=gradcy(ijk)
      eyz(ijk)=0.5*gradcz(ijk) 
  200 continue
!
      call gradc(ncells,ijkcell,      &
          w,gradcx,gradcy,gradcz)
!
      do 300 n=1,ncells
      ijk=ijkcell(n)
      exz(ijk)=exz(ijk)+0.5*gradcx(ijk)
      eyz(ijk)=eyz(ijk)+0.5*gradcy(ijk)
      ezz(ijk)=gradcz(ijk) 
  300 continue
!
      deallocate (gradcx, gradcy, gradcz)
      return
      end subroutine strain
