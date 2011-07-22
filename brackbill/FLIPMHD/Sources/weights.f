      subroutine weights(xi,eta,zeta,wght)
!
!     calculate the weights for trilinear interpolation
!
!test      implicit real*8 (a-h,o-z)
!
      use vast_kind_param, ONLY : double
      implicit none
      real(double) :: wght(*),xi,eta,zeta,omx,ome,omz
!
      omx=1.-xi
      ome=1.-eta
      omz=1.-zeta
!
      wght(1)= xi*omz*ome
      wght(2)= xi*omz*eta
      wght(3)= omx*omz*eta
      wght(4)= omx*omz*ome
      wght(5)= xi*zeta*ome
      wght(6)= xi*zeta*eta
      wght(7)= omx*zeta*eta
      wght(8)= omx*zeta*ome
!
      return
      end subroutine weights
