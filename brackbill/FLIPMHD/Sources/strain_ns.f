      subroutine strain_ns(ncells,ijkcell,iwid,jwid,kwid,     &
          dudx,dudy,dudz,     &
          dvdx,dvdy,dvdz,     &
          dwdx,dwdy,dwdz,     &
          u,v,w,vol)
!
       use vast_kind_param, ONLY : double
       use corgan_com_M, ONLY : itdim
       use Scratch_com_M, ONLY : gradcx, gradcy, gradcz
       use geometry_com_M
!      implicit real*8 (a-h,o-z)
!
      integer ::   &
          ijkcell(*)
      real(double) ::     &
          dudx(*),dudy(*),dudz(*),    &
          dvdx(*),dvdy(*),dvdz(*),    &
          dwdx(*),dwdy(*),dwdz(*),    &
          u(*),v(*),w(*),vol(*)
!
       allocate(gradcx(itdim), gradcy(itdim), gradcz(itdim))
!
!      write(*,*)'ijk = 2664'
!      write(*,*)'gradx,grady,gradz=',gradcx(2664),gradcy(2664),     &
!                gradcz(2664)
      call gradc(ncells,ijkcell,      &
          u,gradcx,gradcy,gradcz)
!      write(*,*)'ijk = 2664'
!      write(*,*)'gradx,grady,gradz=',gradcx(2664),gradcy(2664),     &
!              gradcz(2664)
!
!       write(*,*)'dudz(2664)',dudz(2664)
      do 100 n=1,ncells
      ijk=ijkcell(n)
      dudx(ijk)=gradcx(ijk)
      dudy(ijk)=gradcy(ijk)
      dudz(ijk)=gradcz(ijk) 
  100 continue
!
      call gradc(ncells,ijkcell,      &
          v,gradcx,gradcy,gradcz)
!
      do 200 n=1,ncells
      ijk=ijkcell(n)
      dvdx(ijk)=gradcx(ijk)
      dvdy(ijk)=gradcy(ijk)
      dvdz(ijk)=gradcz(ijk) 
  200 continue
!
      call gradc(ncells,ijkcell,    &
          w,gradcx,gradcy,gradcz)
!
      do 300 n=1,ncells
      ijk=ijkcell(n)
      dwdx(ijk)=gradcx(ijk)
      dwdy(ijk)=gradcy(ijk)
      dwdz(ijk)=gradcz(ijk) 
  300 continue
!
      deallocate (gradcx, gradcy, gradcz)
      return
      end subroutine strain_ns
