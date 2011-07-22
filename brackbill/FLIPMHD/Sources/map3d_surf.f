      subroutine map3d_surf(ijk,ijknew,iwid,jwid,kwid,eps,ifail,   &
                     iphead,itdim,x,y,z,   &
                     xptilde,yptilde,zptilde,   &
                     uptilde,vptilde,wptilde,   &
                     xi,eta,alfa)
!
      use vast_kind_param, ONLY : double
      implicit none
!
      integer :: iphead(*),    &
        ijk,ijknew,iwid,jwid,kwid,ifail,itdim
      real(double) :: x(*),y(*),z(*),     &
     &          a(3,3),r(3),wght(8)
      real(double) :: xptilde,yptilde,zptilde,uptilde,vptilde,wptilde,    &
         xi,eta,alfa
      integer :: l1, l2, iter, np, maxx, l(3)
      real(double) :: xmult,delpxi,delalfa,delpeta,   &
        test,rmax1,rmax2,rmax3,xp,yp,zp,amax,det,eps
!
!
      np=iphead(ijk)
      xi=0.5
      eta=0.5
      alfa=0.0
      iter=0
      ifail=0
!
 3730 continue
!
!
      iter=iter+1
!
      xp=xi*(eta*x(ijknew+kwid+iwid+jwid)      &
            +(1.-eta)*x(ijknew+kwid+iwid))      &
       +(1.-xi)*(eta*x(ijknew+kwid+jwid)      &
                 +(1.-eta)*x(ijknew+kwid))
!
      yp=xi*(eta*y(ijknew+kwid+iwid+jwid)      &
            +(1.-eta)*y(ijknew+kwid+iwid))      &
       +(1.-xi)*(eta*y(ijknew+kwid+jwid)      &
                 +(1.-eta)*y(ijknew+kwid))
!
      zp=xi*(eta*z(ijknew+kwid+iwid+jwid)      &
            +(1.-eta)*z(ijknew+kwid+iwid))      &
       +(1.-xi)*(eta*z(ijknew+kwid+jwid)      &
                 +(1.-eta)*z(ijknew+kwid))
!
      a(1,1)=-(1.-eta)*(x(ijknew+kwid+iwid)-x(ijknew+kwid))      &
            -eta*(x(ijknew+kwid+iwid+jwid)-x(ijknew+kwid+jwid))

      a(1,2)=-(1.-xi)*(x(ijknew+kwid+jwid)-x(ijknew+kwid))      &
            -xi*(x(ijknew+kwid+iwid+jwid)-x(ijknew+kwid+iwid))
!
      a(1,3)=-uptilde
!
      a(2,1)=-(1.-eta)*(y(ijknew+kwid+iwid)-y(ijknew+kwid))      &
            -eta*(y(ijknew+kwid+iwid+jwid)-y(ijknew+kwid+jwid))
!
      a(2,2)=-(1.-xi)*(y(ijknew+kwid+jwid)-y(ijknew+kwid))      &
            -xi*(y(ijknew+kwid+iwid+jwid)-y(ijknew+kwid+iwid))
!
      a(2,3)=-vptilde
!
      a(3,1)=-(1.-eta)*(z(ijknew+kwid+iwid)-z(ijknew+kwid))      &
            -eta*(z(ijknew+kwid+iwid+jwid)-z(ijknew+kwid+jwid))
!
      a(3,2)=-(1.-xi)*(z(ijknew+kwid+jwid)-z(ijknew+kwid))      &
            -xi*(z(ijknew+kwid+iwid+jwid)-z(ijknew+kwid+iwid))
!
      a(3,3)=-wptilde
!
      det=-(a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))      &
          -a(2,3)*(a(1,1)*a(3,2)-a(1,2)*a(3,1))      &
          +a(3,3)*(a(1,1)*a(2,2)-a(1,2)*a(2,1)))
!
      if(det.eq.0.0) then
      ifail=1
      return
      end if
!
      r(1)=xptilde-alfa*uptilde-xp
      r(2)=yptilde-alfa*vptilde-yp
      r(3)=zptilde-alfa*wptilde-zp
!
!     rescale coeeficients
!
      rmax1=1./max(abs(a(1,1)),abs(a(1,2)),abs(a(1,3)))
      rmax2=1./max(abs(a(2,1)),abs(a(2,2)),abs(a(2,3)))
      rmax3=1./max(abs(a(3,1)),abs(a(3,2)),abs(a(3,3)))
!
      a(1,1)=a(1,1)*rmax1
      a(1,2)=a(1,2)*rmax1
      a(1,3)=a(1,3)*rmax1
      r(1)=r(1)*rmax1
!
      a(2,1)=a(2,1)*rmax2
      a(2,2)=a(2,2)*rmax2
      a(2,3)=a(2,3)*rmax2
      r(2)=r(2)*rmax2
!
      a(3,1)=a(3,1)*rmax3
      a(3,2)=a(3,2)*rmax3
      a(3,3)=a(3,3)*rmax3
      r(3)=r(3)*rmax3
!
      l(1)=1
      l(2)=2
      l(3)=3
!
      amax=max(abs(a(1,1)),abs(a(2,1)),abs(a(3,1)))
      maxx=1
      if(abs(a(2,1)).ge.amax) maxx=2
      if(abs(a(3,1)).ge.amax) maxx=3
      l1=l(maxx)
      l(maxx)=l(1)
      l(1)=l1
!
      a(l1,1)=1./a(l1,1)
      xmult=a(l(2),1)*a(l1,1)
      a(l(2),2)=a(l(2),2)-xmult*a(l1,2)
      a(l(2),3)=a(l(2),3)-xmult*a(l1,3)
      r(l(2))=r(l(2))-xmult*r(l1)
!
      xmult=a(l(3),1)*a(l1,1)
      a(l(3),2)=a(l(3),2)-xmult*a(l1,2)
      a(l(3),3)=a(l(3),3)-xmult*a(l1,3)
      r(l(3))=r(l(3))-xmult*r(l1)
!
      maxx=2
      amax=max(abs(a(l(2),2)),abs(a(l(3),2)))
      if(abs(a(l(3),2)).ge.amax) maxx=3
      l2=l(maxx)
      l(maxx)=l(2)
      l(2)=l2
!
      a(l2,2)=1./a(l2,2)
      xmult=a(l(3),2)*a(l2,2)
      a(l(3),3)=a(l(3),3)-xmult*a(l2,3)
      r(l(3))=r(l(3))-xmult*r(l2)
!
      delalfa=-r(l(3))/a(l(3),3)
      delpeta=(-r(l2)-a(l2,3)*delalfa)*a(l2,2)
      delpxi=(-r(l1)-a(l1,3)*delalfa-a(l1,2)*delpeta)*a(l1,1)
!
      xi=xi+delpxi
      eta=eta+delpeta
      alfa=alfa+delalfa
!
!
!
      test=delalfa**2+delpeta**2+delpxi**2
      if(test.le.eps**2) return
      if(iter.le.25) go to 3730
      ifail=2
!
!
  373 continue
!
      return
      end subroutine map3d_surf
