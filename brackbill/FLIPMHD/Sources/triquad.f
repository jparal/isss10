      subroutine triquad(ncells,ijkcell,itdim,ijkc,      &
         bx,by,bz,w,bxp,byp,bzp)
!
!     a routine to interpolate cell-centered variables
!     like the magnetic field to the particles
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      integer :: ijkc(27),w(itdim,27),ijkcell(*)
      real(double) ::bx(*),by(*),bz(*),    &
          bxp(*),byp(*),bzp(*)
!
      do n=1,ncells
      ijk=ijkcell(n)
!

      bxp(ijk)=w(ijk,1)*bx(ijk+ijkc(1))     &
        +(w(ijk,2)*bx(ijk+ijkc(2))     &
        +(w(ijk,3)*bx(ijk+ijkc(3))     &
        +(w(ijk,4)*bx(ijk+ijkc(4))     &
        +(w(ijk,5)*bx(ijk+ijkc(5))     &
        +(w(ijk,6)*bx(ijk+ijkc(6))     &
        +(w(ijk,7)*bx(ijk+ijkc(7))     &
        +(w(ijk,8)*bx(ijk+ijkc(8))     &
        +(w(ijk,9)*bx(ijk+ijkc(9))))))))))
      bxp(ijk)=bxp(ijk)      &
        +(w(ijk,10)*bx(ijk+ijkc(10))      &
        +(w(ijk,11)*bx(ijk+ijkc(11))      &
        +(w(ijk,12)*bx(ijk+ijkc(12))      &
        +(w(ijk,13)*bx(ijk+ijkc(13))      &
        +(w(ijk,14)*bx(ijk+ijkc(14))      &
        +(w(ijk,15)*bx(ijk+ijkc(15))      &
        +(w(ijk,16)*bx(ijk+ijkc(16))      &
        +(w(ijk,17)*bx(ijk+ijkc(17))      &
        +(w(ijk,18)*bx(ijk+ijkc(18)))))))))))
      bxp(ijk)=bxp(ijk)     &
        +(w(ijk,19)*bx(ijk+ijkc(19))     &
        +(w(ijk,20)*bx(ijk+ijkc(20))     &
        +(w(ijk,21)*bx(ijk+ijkc(21))     &
        +(w(ijk,22)*bx(ijk+ijkc(22))     &
        +(w(ijk,23)*bx(ijk+ijkc(23))     &
        +(w(ijk,24)*bx(ijk+ijkc(24))     &
        +(w(ijk,25)*bx(ijk+ijkc(25))     &
        +(w(ijk,26)*bx(ijk+ijkc(26))     &
        +(w(ijk,27)*bx(ijk+ijkc(27)))))))))))
!
      byp(ijk)=w(ijk,1)*by(ijk+ijkc(1))      &
        +(w(ijk,2)*by(ijk+ijkc(2))      &
        +(w(ijk,3)*by(ijk+ijkc(3))      &
        +(w(ijk,4)*by(ijk+ijkc(4))      &
        +(w(ijk,5)*by(ijk+ijkc(5))      &
        +(w(ijk,6)*by(ijk+ijkc(6))      &
        +(w(ijk,7)*by(ijk+ijkc(7))      &
        +(w(ijk,8)*by(ijk+ijkc(8))      &
        +(w(ijk,9)*by(ijk+ijkc(9))))))))))
      byp(ijk)=byp(ijk)      &
        +(w(ijk,10)*by(ijk+ijkc(10))      &
        +(w(ijk,11)*by(ijk+ijkc(11))      &
        +(w(ijk,12)*by(ijk+ijkc(12))      &
        +(w(ijk,13)*by(ijk+ijkc(13))      &
        +(w(ijk,14)*by(ijk+ijkc(14))      &
        +(w(ijk,15)*by(ijk+ijkc(15))      &
        +(w(ijk,16)*by(ijk+ijkc(16))      &
        +(w(ijk,17)*by(ijk+ijkc(17))      &
        +(w(ijk,18)*by(ijk+ijkc(18)))))))))))
      byp(ijk)=byp(ijk)      &
        +(w(ijk,19)*by(ijk+ijkc(19))      &
        +(w(ijk,20)*by(ijk+ijkc(20))      &
        +(w(ijk,21)*by(ijk+ijkc(21))      &
        +(w(ijk,22)*by(ijk+ijkc(22))      &
        +(w(ijk,23)*by(ijk+ijkc(23))      &
        +(w(ijk,24)*by(ijk+ijkc(24))      &
        +(w(ijk,25)*by(ijk+ijkc(25))      &
        +(w(ijk,26)*by(ijk+ijkc(26))      &
        +(w(ijk,27)*by(ijk+ijkc(27)))))))))))
!
      bzp(ijk)=w(ijk,1)*bz(ijk+ijkc(1))      &
        +(w(ijk,2)*bz(ijk+ijkc(2))      &
        +(w(ijk,3)*bz(ijk+ijkc(3))      &
        +(w(ijk,4)*bz(ijk+ijkc(4))      &
        +(w(ijk,5)*bz(ijk+ijkc(5))      &
        +(w(ijk,6)*bz(ijk+ijkc(6))      &
        +(w(ijk,7)*bz(ijk+ijkc(7))      &
        +(w(ijk,8)*bz(ijk+ijkc(8))      &
        +(w(ijk,9)*bz(ijk+ijkc(9))))))))))
      bzp(ijk)=bzp(ijk)      &
        +(w(ijk,10)*bz(ijk+ijkc(10))      &
        +(w(ijk,11)*bz(ijk+ijkc(11))      &
        +(w(ijk,12)*bz(ijk+ijkc(12))      &
        +(w(ijk,13)*bz(ijk+ijkc(13))      &
        +(w(ijk,14)*bz(ijk+ijkc(14))      &
        +(w(ijk,15)*bz(ijk+ijkc(15))      &
        +(w(ijk,16)*bz(ijk+ijkc(16))      &
        +(w(ijk,17)*bz(ijk+ijkc(17))      &
        +(w(ijk,18)*bz(ijk+ijkc(18)))))))))))
      bzp(ijk)=bzp(ijk)      &
        +(w(ijk,19)*bz(ijk+ijkc(19))      &
        +(w(ijk,20)*bz(ijk+ijkc(20))      &
        +(w(ijk,21)*bz(ijk+ijkc(21))      &
        +(w(ijk,22)*bz(ijk+ijkc(22))      &
        +(w(ijk,23)*bz(ijk+ijkc(23))      &
        +(w(ijk,24)*bz(ijk+ijkc(24))      &
        +(w(ijk,25)*bz(ijk+ijkc(25))      &
        +(w(ijk,26)*bz(ijk+ijkc(26))      &
        +(w(ijk,27)*bz(ijk+ijkc(27)))))))))))
!
      enddo
!
      return
      end subroutine triquad
