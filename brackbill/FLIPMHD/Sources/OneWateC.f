      subroutine OneWateC(n,the,zeta,nu,wate)
!
      use vast_kind_param, ONLY : double
      use corgan_com_M, ONLY : itdim
!
      real(double) :: wate(itdim,*)
      real(double) :: the(*),zeta(*),nu(*),   &
        wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp
      integer :: n
!
      wi=0.75-the(n)**2
      wim=0.5*(0.5-the(n))**2
      wip=0.5*(0.5+the(n))**2
!
      wj=0.75-zeta(n)**2
      wjm=0.5*(0.5-zeta(n))**2
      wjp=0.5*(0.5+zeta(n))**2
!
      wk=0.75-nu(n)**2
      wkm=0.5*(0.5-nu(n))**2
      wkp=0.5*(0.5+nu(n))**2
!
!     k-plane
!
      wate(n,1)=wi*wj*wk
      wate(n,2)=wip*wj*wk
      wate(n,3)=wip*wjp*wk
      wate(n,4)=wi*wjp*wk
      wate(n,5)=wim*wjp*wk
      wate(n,6)=wim*wj*wk
      wate(n,7)=wim*wjm*wk
      wate(n,8)=wi*wjm*wk
      wate(n,9)=wip*wjm*wk
!
!     k-1 - plane
!
      wate(n,10)=wi*wj*wkm
      wate(n,11)=wip*wj*wkm
      wate(n,12)=wip*wjp*wkm
      wate(n,13)=wi*wjp*wkm
      wate(n,14)=wim*wjp*wkm
      wate(n,15)=wim*wj*wkm
      wate(n,16)=wim*wjm*wkm
      wate(n,17)=wi*wjm*wkm
      wate(n,18)=wip*wjm*wkm
!
!     k+1 - plane
!
      wate(n,19)=wi*wj*wkp
      wate(n,20)=wip*wj*wkp
      wate(n,21)=wip*wjp*wkp
      wate(n,22)=wi*wjp*wkp
      wate(n,23)=wim*wjp*wkp
      wate(n,24)=wim*wj*wkp
      wate(n,25)=wim*wjm*wkp
      wate(n,26)=wi*wjm*wkp
      wate(n,27)=wip*wjm*wkp
!
      return
      end   !subroutine OneWateC

