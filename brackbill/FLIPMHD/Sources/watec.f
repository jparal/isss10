      subroutine watec(ncells,ijkcell,iphead,itdim,     &
                      pxi,peta,pzta,wate)
!
!
!     a routine to calculate the weights for triquadratic
!     interpolation
!
      use vast_kind_param, only:  double
!test      implicit real*8 (a-h,o-z)
      implicit none
!
      integer :: ncells,ijkcell(*),iphead(*),     &
          itdim,np,i,j,k,ijk,n
      real(double) ::  pxi(0:*),peta(0:*),pzta(0:*),   &
          wate(itdim,*),the,zeta,nu,   &
          wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp
!
      do 1 n=1,ncells
!
      ijk=ijkcell(n)
!
      np=iphead(ijk)
      i=int(pxi(np))
      j=int(peta(np))
      k=int(pzta(np))
!
      the=pxi(np)-real(i)-0.5
      zeta=peta(np)-real(j)-0.5
      nu=pzta(np)-real(k)-0.5
!
      wi=0.75-the**2
      wim=0.5*(0.5-the)**2
      wip=0.5*(0.5+the)**2
!
      wj=0.75-zeta**2
      wjm=0.5*(0.5-zeta)**2
      wjp=0.5*(0.5+zeta)**2
!
      wk=0.75-nu**2
      wkm=0.5*(0.5-nu)**2
      wkp=0.5*(0.5+nu)**2
!
!     k-plane
!
      wate(ijk,1)=wi*wj*wk
      wate(ijk,2)=wip*wj*wk
      wate(ijk,3)=wip*wjp*wk
      wate(ijk,4)=wi*wjp*wk
      wate(ijk,5)=wim*wjp*wk
      wate(ijk,6)=wim*wj*wk
      wate(ijk,7)=wim*wjm*wk
      wate(ijk,8)=wi*wjm*wk
      wate(ijk,9)=wip*wjm*wk
!
!     k-1 - plane
!
      wate(ijk,10)=wi*wj*wkm
      wate(ijk,11)=wip*wj*wkm
      wate(ijk,12)=wip*wjp*wkm
      wate(ijk,13)=wi*wjp*wkm
      wate(ijk,14)=wim*wjp*wkm
      wate(ijk,15)=wim*wj*wkm
      wate(ijk,16)=wim*wjm*wkm
      wate(ijk,17)=wi*wjm*wkm
      wate(ijk,18)=wip*wjm*wkm
!
!     k+1 - plane
!
      wate(ijk,19)=wi*wj*wkp
      wate(ijk,20)=wip*wj*wkp
      wate(ijk,21)=wip*wjp*wkp
      wate(ijk,22)=wi*wjp*wkp
      wate(ijk,23)=wim*wjp*wkp
      wate(ijk,24)=wim*wj*wkp
      wate(ijk,25)=wim*wjm*wkp
      wate(ijk,26)=wi*wjm*wkp
      wate(ijk,27)=wip*wjm*wkp
!
    1 continue
!
      return
      end subroutine watec
