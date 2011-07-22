      subroutine watev(ncells,ijkcell,iphead,itdim,     &
         pxi,peta,pzta,wate)
!
!     a routine to calculate the interpolation weights for trilinear
!     interpolation
!
!test      implicit real*8 (a-h,o-z)
      use vast_kind_param, ONLY : double
      implicit none
!
      integer :: ijkcell(*),iphead(*),    &
          inew,jnew,knew,ijk,    &
          ncells,itdim,n
      real(double) :: pxi(0:*),peta(0:*),pzta(0:*),    &
          wate(itdim,*),pxi2,peta2,pzta2,     &
          wi,wip,wj,wjp,wk,wkp
!
      do 372 n=1,ncells
!
      ijk=ijkcell(n)
!
      inew=int(pxi(iphead(ijk)))
      jnew=int(peta(iphead(ijk)))
      knew=int(pzta(iphead(ijk)))
!
!
      pxi2=pxi(iphead(ijk))-real(inew)
      peta2=peta(iphead(ijk))-real(jnew)
      pzta2=pzta(iphead(ijk))-real(knew)
!
!     calculate interpolation weights
!
      wi=1.-pxi2
      wip=pxi2
!
      wj=1.-peta2
      wjp=peta2
!
      wk=1.-pzta2
      wkp=pzta2
!
      wate(ijk,1)=wip*wj*wk
      wate(ijk,2)=wip*wjp*wk
      wate(ijk,3)=wi*wjp*wk
      wate(ijk,4)=wi*wj*wk
!
      wate(ijk,5)=wip*wj*wkp
      wate(ijk,6)=wip*wjp*wkp
      wate(ijk,7)=wi*wjp*wkp
      wate(ijk,8)=wi*wj*wkp
!
  372 continue
!
      return
      end subroutine watev
