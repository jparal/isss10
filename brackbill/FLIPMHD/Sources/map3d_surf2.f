      subroutine map3d_surf2(ijk,iwid,jwid,kwid,eps,ifail,
     &                iphead,itdim,wate,x,y,z,
     &                xptilde,yptilde,zptilde,
     &                uptilde,vptilde,wptilde,
     &                xi,eta,zta,alfa,dt)
c
      implicit real*8 (a-h,o-z)
c
      dimension iphead(*),wate(itdim,*),x(*),y(*),z(*)
c
c
      iter=0
c
      np=iphead(ijk)
      inew=int(xi)
      jnew=int(eta)
      knew=int(zta)
      ijknew=(knew-1)*kwid+(jnew-1)*jwid+(inew-1)*iwid+1
      alfmin=0.0
      alfmax=dt
      alfold=0.0
c
 3730 continue
c
      iter=iter+1
c
      alfa=0.5*(alfmin+alfmax)      
      delalfa=alfa-alfold
      alfold=alfa
c
      xptilde=xptilde-delalfa*uptilde
      yptilde=yptilde-delalfa*vptilde
      zptilde=zptilde-delalfa*wptilde
c
      call map3d(ijknew,iwid,jwid,kwid,eps,
     &           x,y,z,xptilde,yptilde,zptilde,
     &           delxi,deleta,delnu)
c
      if(abs(delnu).lt.eps) then
        ifail=0
        xi=real(inew)+delxi
        eta=real(jnew)+deleta
        zta=real(knew)+delnu
        return
      end if
c
      if(delnu.gt.0.0) then
        alfmin=alfa
      else
        alfmax=alfa
      endif
c
      if(alfmax-alfmin.lt.eps*dt.or.iter.gt.25) then
      ifail=2
      xi=real(inew)+delxi
      eta=real(jnew)+deleta
      zta=real(knew)+delnu
      return
      endif
c
      go to 3730
c
      end
