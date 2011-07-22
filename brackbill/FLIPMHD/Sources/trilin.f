      subroutine trilin(ncells,ijkcell,itdim,iwid,jwid,kwid,    &
          wate,bxv,byv,bzv,bxpn,bypn,bzpn)
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
      real(double) :: wate(itdim,*),bxv(*),byv(*),bzv(*),    &
          bxpn(*),bypn(*),bzpn(*)
      integer :: ijkcell(*)
!
      do 41 n=1,ncells
!
      ijk=ijkcell(n)
!
!
      bxpn(ijk)=(wate(ijk,1)*bxv(ijk+iwid)   &
               +(wate(ijk,2)*bxv(ijk+iwid+jwid)   &
               +(wate(ijk,3)*bxv(ijk+jwid)   &
               +(wate(ijk,4)*bxv(ijk)   &
               +(wate(ijk,5)*bxv(ijk+iwid+kwid)   &
               +(wate(ijk,6)*bxv(ijk+iwid+jwid+kwid)   &
               +(wate(ijk,7)*bxv(ijk+jwid+kwid)   &
               +(wate(ijk,8)*bxv(ijk+kwid)))))))))
!
      bypn(ijk)=(wate(ijk,1)*byv(ijk+iwid)     &
               +(wate(ijk,2)*byv(ijk+iwid+jwid)     &
               +(wate(ijk,3)*byv(ijk+jwid)     &
               +(wate(ijk,4)*byv(ijk)     &
               +(wate(ijk,5)*byv(ijk+iwid+kwid)     &
               +(wate(ijk,6)*byv(ijk+iwid+jwid+kwid)     &
               +(wate(ijk,7)*byv(ijk+jwid+kwid)     &
               +(wate(ijk,8)*byv(ijk+kwid)))))))))
!
      bzpn(ijk)=(wate(ijk,1)*bzv(ijk+iwid)     &
               +(wate(ijk,2)*bzv(ijk+iwid+jwid)     &
               +(wate(ijk,3)*bzv(ijk+jwid)     &
               +(wate(ijk,4)*bzv(ijk)     &
               +(wate(ijk,5)*bzv(ijk+iwid+kwid)     &
               +(wate(ijk,6)*bzv(ijk+iwid+jwid+kwid)     &
               +(wate(ijk,7)*bzv(ijk+jwid+kwid)     &
               +(wate(ijk,8)*bzv(ijk+kwid)))))))))
!
   41 continue
!
      return
      end subroutine trilin
