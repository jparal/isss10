 
      subroutine area_bc(i1,i2,j1,j2,ks,iwid,jwid,kwid,updown,      &
          x,y,z,ax,ay,az)
!
!     a routine to calculate the directed areas over a surface of constant index ks
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      integer, intent(in) :: i1,i2,j1,j2,ks,iwid,jwid,kwid
      real(double) ::  x(*),y(*),z(*),ax(*),ay(*),az(*)
      real(double), intent(in) :: updown
!
      do 1 j=j1,j2
           do 1 i=i1,i2
!
      ijk=(ks-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
!
      x1=0.5*(x(ijk+iwid)+x(ijk+iwid+jwid)    &
            -x(ijk)-x(ijk+jwid))
!
      x2=0.5*(x(ijk+jwid)+x(ijk+iwid+jwid)    &
            -x(ijk)-x(ijk+iwid))
!
      y1=0.5*(y(ijk+iwid)+y(ijk+iwid+jwid)    &
            -y(ijk)-y(ijk+jwid))
!
      y2=0.5*(y(ijk+jwid)+y(ijk+iwid+jwid)    &
            -y(ijk)-y(ijk+iwid))
!
      z1=0.5*(z(ijk+iwid)+z(ijk+iwid+jwid)    &
            -z(ijk)-z(ijk+jwid))
!
      z2=0.5*(z(ijk+jwid)+z(ijk+iwid+jwid)    &
            -z(ijk)-z(ijk+iwid))
!
      ax(ijk)=(y1*z2-y2*z1)*updown
      ay(ijk)=(z1*x2-z2*x1)*updown
      az(ijk)=(x1*y2-x2*y1)*updown
!
    1 continue
!
      return
      end subroutine area_bc
!     *************************************************************************
      subroutine chek_surf(i1,i2,j1,j2,ksurf,iwid,jwid,kwid,    &
          x,y,z,    &
          px,py,pz,up,vp,wp,ax,ay,az,dt,    &
          nlist,ijktmp3)
!
!     a routine to make a list of those polygons the particle orbit
!     may have intersected
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      integer, intent(in) :: i1,i2,j1,j2,iwid,jwid,kwid
      real(double) :: px,py,pz,up,vp,wp,dt
      real(double) :: ax(*),ay(*),az(*),   &
          x(*),y(*),z(*)
      integer :: nlist,ijktmp3(*)
!
      nlist=0
!
      do 1 i=i1,i2
      do 1 j=j1,j2
!
      ijk=(ksurf-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
!
      testout=(px-x(ijk))*ax(ijk)     &
               +(py-y(ijk))*ay(ijk)     &
               +(pz-z(ijk))*az(ijk)
!
      testin=(up*ax(ijk)+vp*ay(ijk)+wp*az(ijk))*dt
!
!test    ********************************************
!
      if(testout.ge.0..and.testout.lt.testin) then
!
      nlist=nlist+1
      ijktmp3(nlist)=ijk-kwid
      endif
!
    1 continue
!
      return
      end subroutine chek_surf
!     **************************************************************************
      subroutine closest_vertex(i1,i2,j1,j2,k1,k2,    &
          x,y,z,xp,yp,zp,    &
          imin,jmin,kmin,ijkmin)
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid, jwid, kwid
      implicit real*8 (a-h,o-z)
!
      integer, intent(in) :: i1,i2,j1,j2,k1,k2
      integer :: imin,jmin,kmin,ijkmin
      real(double) :: xp,yp,zp
      real(double) :: x(*),y(*),z(*)
      real(double) :: d,dmin
      integer :: i,j,k,ijk
!
      dmin=1.e10
!
      do 1 k=k1,k2
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
!
      d=sqrt((xp-x(ijk))**2+(yp-y(ijk))**2+(zp-z(ijk))**2)
!
      if(d.lt.dmin) then
!
      dmin=d
      imin=i
      jmin=j
      kmin=k
      ijkmin=ijk
!
      endif
!
    1 continue
!
      return
      end subroutine closest_vertex
!     ************************************************************
      subroutine look_evrywher(i1,i2,j1,j2,k1,k2,    &
          eps,succes,    &
          x,y,z,xp,yp,zp,pxi,peta,pzta)
!
!     a routine to examine every cell in the mesh to locate a particle
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      integer, intent(in) :: i1,i2,j1,j2,k1,k2
      real(double) :: x(*),y(*),z(*)
      real(double) :: xp,yp,zp,pxi,peta,pzta
      real(double) :: eps
!
      logical :: succes
      real(double) :: tsi,eta,nu
!
      succes=.false.
!
      do 1 k=k1,k2
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
!
      call map3d(ijk,iwid,jwid,kwid,eps,    &
          x,y,z,    &
          xp,yp,zp,tsi,eta,nu)
!
      if (tsi*(1.-tsi).ge.0.     &
          .and.eta*(1.-eta).ge.0.     &
          .and.nu*(1.-nu).ge.0.) then
!
      succes=.true.
!
      pxi=i+tsi
      peta=j+eta
      pzta=k+nu
!
      go to 601
!
      endif
!
    1 continue
!
  601 continue
!
      return
      end subroutine look_evrywher
!     *********************************************************************************
      subroutine parbctor(ncellsp,ijkctmp,iphead,npart,    &
          cdlt,sdlt,    &
          x,y,z,uptilde,vptilde,wptilde,    &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
!     **********************************************
!
!     a routine to apply boundary conditions to particle motion in a torus
!
!     **********************************************
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid,ibar,jbar,kbar,   &
          ibp1,jbp1,kbp1
      implicit real*8 (a-h,o-z)
!
      integer :: ncellsp,ijkctmp(*),iphead(*),npart
      real(double) ::       &
          x(*),y(*),z(*),uptilde(*),vptilde(*),wptilde(*),     &
          px(0:npart),py(0:npart),pz(0:npart),     &
          up(0:npart),vp(0:npart),wp(0:npart),     &
          pxi(0:npart),peta(0:npart),pzta(0:npart)     
!
      logical :: succes
      real(double) :: tsi,eta,nu
      integer :: inew,jnew,knew
    
!
!
      eps=1.e-4
!
!
      do 374 n=1,ncellsp
!
      ijk=ijkctmp(n)
      np=iphead(ijk)
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
  380 continue
!
      inew=int(pxi(np))
      if((ibp1-inew)*(inew-2).lt.0) then
      if(inew.lt.2) pxi(np)=pxi(np)+real(ibar)
      if(inew.gt.ibp1) pxi(np)=pxi(np)-real(ibar)
      go to 380
      endif
!
!
  390 continue
!
      jnew=int(peta(np))
!
      if((jbp1-jnew)*(jnew-2).lt.0) then
!
      if(jnew.lt.2) then
!
      peta(np)=peta(np)+real(jbar)
!
      wspx=cdlt*px(np)-sdlt*py(np)
      wspy=sdlt*px(np)+cdlt*py(np)
      px(np)=wspx
      py(np)=wspy
!
      wsup=cdlt*up(np)-sdlt*vp(np)
      wsvp=sdlt*up(np)+cdlt*vp(np)
      up(np)=wsup
      vp(np)=wsvp
!
      wsup=cdlt*uptilde(ijk)-sdlt*vptilde(ijk)
      wsvp=sdlt*uptilde(ijk)+cdlt*vptilde(ijk)
      uptilde(ijk)=wsup
      vptilde(ijk)=wsvp
!
      endif
!
      if(jnew.gt.jbp1) then
!
      peta(np)=peta(np)-real(jbar)
!
      wspx=cdlt*px(np)+sdlt*py(np)
      wspy=-sdlt*px(np)+cdlt*py(np)
      px(np)=wspx
      py(np)=wspy
!
      wsup=cdlt*up(np)+sdlt*vp(np)
      wsvp=-sdlt*up(np)+cdlt*vp(np)
      up(np)=wsup
      vp(np)=wsvp
!
      wsup=cdlt*uptilde(ijk)+sdlt*vptilde(ijk)
      wsvp=-sdlt*uptilde(ijk)+cdlt*vptilde(ijk)
      uptilde(ijk)=wsup
      vptilde(ijk)=wsvp
!
      endif
!
      go to 390
!
      endif
!
  374 continue
!
      return
      end subroutine parbctor
!     **************************************************************************
      subroutine parbccyl(ncellsp,ijkctmp,iphead,npart,   &
          dz,   &
          x,y,z,px,py,pz,up,vp,wp,pxi,peta,pzta)
!
!     **********************************************
!
!     a routine to apply boundary conditions to particle motion in a torus
!
!     **********************************************
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY: iwid,jwid,kwid,ibar,jbar,kbar
      implicit real*8 (a-h,o-z)
!
      integer :: ncellsp,ijkctmp(*),npart,iphead(*)
      real(double) ::    &
          x(*),y(*),z(*),     &
          px(0:npart),py(0:npart),pz(0:npart),     &
          up(0:npart),vp(0:npart),wp(0:npart),     &
          pxi(0:npart),peta(0:npart),pzta(0:npart)     
!
      logical :: succes
      real(double) ::  nu
!
      ibp1=ibar+1
      jbp1=jbar+1
      kbp1=kbar+1
!
      eps=1.e-4
!
      do 374 n=1,ncellsp
!
      ijk=ijkctmp(n)
      np=iphead(ijk)
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
  380 continue
!
      inew=int(pxi(np))
      if((ibp1-inew)*(inew-2).lt.0) then
      if(inew.lt.2) pxi(np)=pxi(np)+real(ibar)
      if(inew.gt.ibp1) pxi(np)=pxi(np)-real(ibar)
      go to 380
      endif
!
!
  390 continue
!
      jnew=int(peta(np))
!
      if((jbp1-jnew)*(jnew-2).lt.0) then
!
      if(jnew.lt.2) then
!
      peta(np)=peta(np)+real(jbar)
!
      py(np)=py(np)+dz
!
      endif
!
      if(jnew.gt.jbp1) then
!
      peta(np)=peta(np)-real(jbar)
!
      py(np)=py(np)-dz
!
      endif
!
      go to 390
!
      endif
!
  374 continue
!
      return
      end subroutine parbccyl
!     *******************************************************************
      subroutine parlocat(ncellsp,ijkctmp,iphead,    &
          IJKTMP2,IJKTMP3,IJKTMP4,RMAJ,DZ,ALFA,AX,AY,AZ,DT,    &
          itdim,npart,mgeom,cdlt,sdlt,    &
          wate,x,y,z,BXV,BYV,BZV,    &
          XPTILDE,YPTILDE,ZPTILDE,UPTILDE,VPTILDE,WPTILDE,    &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,    &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
!     ********************************************************
!
!     a routine to locate a particle on a grid
!
!     *******************************************************
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid,ibar,jbar,kbar
      implicit real*8 (a-h,o-z)
!
       integer :: ijkctmp(*),iphead(*)
       real(double) ::          &
          wate(itdim,*),         &
          IJKTMP2(*),IJKTMP3(*),IJKTMP4(*),ALFA(*),         &
          x(*),y(*),z(*),BXV(*),BYV(*),BZV(*),         &
          XPTILDE(*),YPTILDE(*),ZPTILDE(*),         &
          UPTILDE(*),VPTILDE(*),WPTILDE(*),         &
          AX(*),AY(*),AZ(*),         &
          tsix(*),tsiy(*),tsiz(*),         &
          etax(*),etay(*),etaz(*),         &
          nux(*),nuy(*),nuz(*),         &
          px(0:npart),py(0:npart),pz(0:npart),         &
          up(0:npart),vp(0:npart),wp(0:npart),         &
          pxi(0:npart),peta(0:npart),pzta(0:npart)
!
      real(double) :: nu
!
      logical :: succes
!
      one=1.
!
      ibp1=ibar+1
      jbp1=jbar+1
      kbp1=kbar+1
!
      do 380 n=1,ncellsp
!
      ijktmp2(n)=ijkctmp(n)
!
  380 continue
!
      nchek=ncellsp
!
!
!
!     ******************************************************************
!
!     calculate new natural coordinates
!
!     ******************************************************************
      itry=0
!
!
    1 continue
!
      itry=itry+1
!
!
      nchek_nxt=0
!
      do 373 n=1,nchek
!
      ijk=ijktmp2(n)
      np=iphead(ijk)
!
      iold=int(pxi(np))
      iold=min(iold,ibp1)
      iold=max(2,iold)
!
      jold=int(peta(np))
      jold=min(jold,jbp1)
      jold=max(2,jold)
!
      kold=int(pzta(np))
      kold=min(kold,kbp1)
      kold=max(2,kold)
!
!
      ijkold=(kold-1)*kwid+(jold-1)*jwid+(iold-1)*iwid+1
!
      eps=1.e-4
!
      call map3d(ijkold,iwid,jwid,kwid,eps,   &
          x,y,z,   &
          px(np),py(np),pz(np),tsi,eta,nu)
!
!      epsn=1.e-7
      epsn=0.0
      ichange=0.5*(int(sign(one,tsi-1.-epsn)    &
          +sign(one,tsi+epsn)))
!
      jchange=0.5*(int(sign(one,eta-1.-epsn)    &
          +sign(one,eta+epsn)))
!
      kchange=0.5*(int(sign(one,nu-1.-epsn)     &
          +sign(one,nu+epsn)))
!
      inew=iold+ichange
      jnew=jold+jchange
      knew=kold+kchange
!
      ijknew=(knew-1)*kwid+(jnew-1)*jwid+(inew-1)*iwid+1
!
      if(ijknew.ne.ijkold) then
!
      nchek_nxt=nchek_nxt+1
      ijktmp2(nchek_nxt)=ijktmp2(n)
!
!
      pxi(np)=inew+0.5
      peta(np)=jnew+0.5
      pzta(np)=knew+0.5
!
!
      else
!
       pxi(np)=real(iold)+tsi
       peta(np)=real(jold)+eta
       pzta(np)=real(kold)+nu
!
      endif
!
!
  373 continue
!
      nchek=nchek_nxt
!
!
!      check for particles outside the domain
!
       nout=0
!
      do 400 n=1,nchek
      ijk=ijktmp2(n)
      np=iphead(ijk)
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
      if(inew.lt.2.or.inew.gt.ibp1     &
          .or.jnew.lt.2.or.jnew.gt.jbp1     &
          .or.knew.lt.2.or.knew.gt.kbp1) then
!
      nout=nout+1
      ijktmp3(nout)=ijk
      endif
!
  400 continue
!
!     check for axis crossings
!
      do 384 n=1,nout
!
      ijk=ijktmp3(n)
      np=iphead(ijk)
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
         if(jnew.ge.2.and.jnew.le.jbp1    &
          .and.knew.lt.2) then
!
      iaxis=inew
      jaxis=jnew
!
!     the particle has crossed the axis
!
      call map3d_axis(2,ibp1,jnew,jnew,iaxis,jaxis, iwid,jwid,kwid,eps,    &
          x,y,z,    &
          px(np),py(np),pz(np),tsi,eta,nu,succes)
!
      if(.not.succes) then
!
      call map3d_axis(2,ibp1,2,jbp1,iaxis,jaxis,iwid,jwid,kwid,eps,     &
          x,y,z,     &
          px(np),py(np),pz(np),tsi,eta,nu,succes)
!
      endif
      if(succes) then
      pxi(np)=iaxis+tsi
      peta(np)=jaxis+eta
      pzta(np)=2.+nu
      else
      pzta(np)=2.
      endif
!
      endif
!
  384 continue
!
!_parloco      if(nout.gt.0.and.itry.eq.1) then
!_parloco      call parloco(nout,ijktmp3,iphead,     &
!_parloco          RMAJ,DZ,ALFA,DT,ifail,      &
!_parloco          mgeom,cdlt,sdlt,      &
!_parloco          wate,x,y,z,bxv,byv,bzv,      &
!_parloco          xptilde,yptilde,zptilde,uptilde,vptilde,wptilde,      &
!_parloco          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,      &
!_parloco          px,py,pz,up,vp,wp,pxi,peta,pzta)
!_parloco      endif
!
      do 374 n=1,nout
!
      ijk=ijktmp3(n)
      np=iphead(ijk)
!
 3810 continue
!
      inew=int(pxi(np))
!
      if(inew.lt.2.or.inew.gt.ibp1) then
!
      if(inew.lt.2) pxi(np)=pxi(np)+real(ibar)
      if(inew.gt.ibp1) pxi(np)=pxi(np)-real(ibar)
!
      go to 3810
!
      endif
!
  374 continue
!
!     nout is the number of particles for which intersections
!     have still not been found
!
!
!
      if(nout.gt.0) then
!
      nout_nxt=0
!
!     calculate directed areas over k=kbp2 surface
!
      updown=1.
      call area_bc(2,ibp1,2,jbp1,kbp1+1,iwid,jwid,kwid,updown,     &
          x,y,z,ax,ay,az)
!
!     one particle at a time, make list of surface elements
!     which the particle may have passed through
!
      do 510 n=1,nout
      ijk=ijktmp3(n)
      np=iphead(ijk)
!
      succes=.false.
!
      knew=int(pzta(np))
      if(knew.gt.kbp1) then
!
!
      call chek_surf(2,ibp1,2,jbp1,kbp1+1,iwid,jwid,kwid,     &
          x,y,z,     &
          px(np),py(np),pz(np),uptilde(ijk),vptilde(ijk),wptilde(ijk),     &
          ax,ay,az,dt,nlist,ijktmp4)
!
!      if(nlist.gt.0) write(*,*) 'chek_surf, nlist=',nlist
!     calculate intersection with boundary in candidate cells
!
      knew=kbp1
!
      do 520 nl=1,nlist
      ijknew=ijktmp4(nl)
!
      call map3d_surf(ijk,ijknew,iwid,jwid,kwid,eps,ifail,    &
          iphead,itdim,x,y,z,    &
          px(np),py(np),pz(np),    &
          uptilde(ijk),vptilde(ijk),wptilde(ijk),    &
          xi,eta,alfa(ijk))
!
!      write(*,*) 'np=',np,'ijknew=',ijknew,'ifail=',ifail
!
      if(ifail.ne.1) then
      alfa(ijk)=alfa(ijk)/dt
!
      if(xi*(1.-xi).ge.0.     &
          .and.eta*(1.-eta).ge.0.     &
          .and.alfa(ijk)*(1.-alfa(ijk)).ge.0.) then
!
!     particle has been found
!
      succes=.true.
!
      knew=(ijknew-1)/kwid+1
      jnew=(ijknew-(knew-1)*kwid-1)/jwid+1
      inew=(ijknew-(knew-1)*kwid-(jnew-1)*jwid-1)/iwid+1
!
      pxi(np)=xi+inew
      peta(np)=eta+jnew
!
      pzta(np)=kbp1+1-1.e-6
!
      ll=1
!
      call parrefl(ll,ijktmp3(n),iphead,npart,itdim,     &
          iwid,jwid,kwid,x,y,z,wate,bxv,byv,bzv,     &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
!      write(*,*) 'particle found and reflected, np=', np
!
      go to 521
!
      endif
!
      endif
!
  520 continue
!
      endif
!
  521 if(.not.succes) then
!
      nout_nxt=nout_nxt+1
      ijktmp3(nout_nxt)=ijk
!
!
      endif
!
  510 continue
!
      nout=nout_nxt
      endif
!
!
!     check j=jbp2 surface for crossings
!
      if(nout.gt.0) then
      nout_nxt=0
!
      updown=1.
!
      call area_bc(2,kbp1,2,ibp1,jbp1+1,kwid,iwid,jwid,updown,    &
          x,y,z,ax,ay,az)
!
!
      do 610 n=1,nout
!
      ijk=ijktmp3(n)
      np=iphead(ijk)
      jsurf=jbp1+1
!
!      if nlist>0, call periodic bc's
!
      call chek_surf(2,kbp1,2,ibp1,jsurf,kwid,iwid,jwid,    &
          x,y,z,    &
          px(np),py(np),pz(np),uptilde(ijk),vptilde(ijk),wptilde(ijk),    &
          ax,ay,az,dt,nlist,ijktmp4)
!
      succes=.false.
!
      if(nlist.gt.0) then
!
!     ******************************************************************
!     DIAGNOSTIC
      cross=-1.
      ijksurf=5*kwid+jbp1*jwid+5*iwid+1
      phisurf=atan(y(ijksurf)/x(ijksurf))
      parphinp=atan(py(np)/(px(np)+1.e-6))
      parphin=atan((py(np)-vp(np)*dt)/(px(np)-up(np)*dt+1.e-6))
      if(parphinp.gt.phisurf.and.parphin.lt.phisurf) then
!      write(*,*) 'particle np=',np,'has crossed j=jbp2'
      cross=1.
      endif
      do 620 nl=1,nlist
!
      ijknew=ijktmp4(nl)
!
      call map3d_surf(ijk,ijknew,kwid,iwid,jwid,eps,ifail,    &
          iphead,itdim,x,y,z,    &
          px(np),py(np),pz(np),    &
          uptilde(ijk),vptilde(ijk),wptilde(ijk),    &
          zta,xi,alfa(ijk))
!
      if(ifail.eq.0) then
!
      alfa(ijk)=alfa(ijk)/dt
!      write(*,*) 'loop 610, ifail=',ifail,'ijknew,np=',ijknew,np
!
      if(xi*(1.-xi).ge.0.     &
          .and.alfa(ijk)*(1.-alfa(ijk)).ge.0.     &
          .and.zta*(1.-zta).ge.0.) then
!
      succes=.true.
!
!      write(*,*) '610 loop, xi,zta,alfa=', xi,zta,alfa(ijk)
      peta(np)=jsurf
      if(rmaj.gt.0.0) then
!
       call parbctor(1,ijktmp3(n),iphead,npart,     &
          cdlt,sdlt,     &
          x,y,z,uptilde,vptilde,wptilde,     &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
      else
!
       call parbccyl(1,ijktmp3(n),iphead,npart,     &
          dz,     &
          x,y,z,px,py,pz,up,vp,wp,pxi,peta,pzta)
!
      endif
      go to 621
      endif
!
      endif
  620 continue
      endif
!
  621 continue
      if(.not.succes) then
!
      nout_nxt=nout_nxt+1
      ijktmp3(nout_nxt)=ijk
      endif
!
!
  610 continue
!
      nout=nout_nxt
      endif
!
      if(nout.gt.0) then
      nout_nxt=0
!
      updown=-1.
      call area_bc(2,kbp1,2,ibp1,2,kwid,iwid,jwid,updown,    &
          x,y,z,ax,ay,az)
!
!
      do 710 n=1,nout
!
      ijk=ijktmp3(n)
      np=iphead(ijk)
      jsurf=2
!
!      if nlist>0, call periodic bc's
!
      call chek_surf(2,kbp1,2,ibp1,jsurf,kwid,iwid,jwid,     &
          x,y,z,     &
          px(np),py(np),pz(np),uptilde(ijk),vptilde(ijk),wptilde(ijk),     &
          ax,ay,az,dt,nlist,ijktmp4)
!
      succes=.false.
!
      if(nlist.gt.0) then
!
!     ******************************************************************
!     DIAGNOSTIC
      cross=-1.
      ijksurf=5*kwid+jwid+5*iwid+1
      phisurf=atan(y(ijksurf)/x(ijksurf))
      parphinp=atan(py(np)/(px(np)+1.e-6))
      parphin=atan((py(np)-vp(np)*dt)/(px(np)-up(np)*dt+1.e-6))
      if(parphinp.lt.phisurf.and.parphin.gt.phisurf) then
!      write(*,*) 'particle np=',np,'has crossed j=2'
      cross=1.
      endif
      do 720 nl=1,nlist
!
      ijknew=ijktmp4(nl)
!
      call map3d_surf(ijk,ijknew,kwid,iwid,jwid,eps,ifail,   &
          iphead,itdim,x,y,z,   &
          px(np),py(np),pz(np),   &
          uptilde(ijk),vptilde(ijk),wptilde(ijk),   &
          zta,xi,alfa(ijk))
!
      if(ifail.eq.0) then
      alfa(ijk)=alfa(ijk)/dt
!      write(*,*) 'loop 710, ifail=',ifail,'ijknew,np=',ijknew,np
!
      if(xi*(1.-xi).ge.0.     &
          .and.alfa(ijk)*(1.-alfa(ijk)).ge.0.     &
          .and.zta*(1.-zta).ge.0.) then
!
      succes=.true.
!
!      write(*,*) '710 loop, xi,zta,alfa=', xi,zta,alfa(ijk)
      peta(np)=jsurf-1
      if(rmaj.gt.0.0) then
!
       call parbctor(1,ijktmp3(n),iphead,npart,    &
          cdlt,sdlt,    &
          x,y,z,uptilde,vptilde,wptilde,    &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
      else
!
       call parbccyl(1,ijktmp3(n),iphead,npart,     &
          dz,     &
          x,y,z,px,py,pz,up,vp,wp,pxi,peta,pzta)
!
      endif
      go to 721
!
      endif
!
      endif
  720 continue
      endif
!
  721 continue
      if(.not.succes) then
!
      nout_nxt=nout_nxt+1
      ijktmp3(nout_nxt)=ijk
      endif
!
  710 continue
!
      nout=nout_nxt
!
      endif
!
!     **********************************************
      if(nchek.gt.0.and.itry.lt.10) go to 1
!
      nchek_nxt=0
!
!    after an exhaustive search of the boundary,
!
!c      do a global search for the particle
!c
!      if(nchek.gt.0) then
!c
!      do 810 n=1,nchek
!      ijk=ijktmp2(n)
!      np=iphead(ijk)
!c
!      call look_evrywher(2,ibp1,2,jbp1,2,kbp1,  &
!          eps,succes,  &
!          x,y,z,  &
!          px(np),py(np),pz(np),  &
!          pxi(np),peta(np),pzta(np))
!c
!c       after global search
!      if(succes) then
!c
!       write(*,*) 'look_evrywher called, succes=', succes
!       write(*,*)  'pxi,peta,pzta=',pxi(np),peta(np),pzta(np)
!      else
!c
!      nchek_nxt=nchek_nxt+1
!      ijktmp2(nchek_nxt)=ijk
!c
!      endif
!
  810 continue
!
      if(nchek.gt.0) then
!      write(*,*) 'nchek particles',nchek,'not found at ncyc=',ncyc
!      endif
!
      endif
!
 1000 continue
!
!
      return
      end subroutine parlocat
!     ***************************************************************
      subroutine parloco(nout,ijkctmp,iphead,     &
          RMAJ,DZ,ALFA,DT,ifail,     &
          mgeom,cdlt,sdlt,     &
          wate,x,y,z,bxv,byv,bzv,     &
          xptilde,yptilde,zptilde,uptilde,vptilde,wptilde,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
!     ********************************************************
!
!     a routine to locate a particle on a grid
!
!     *******************************************************
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid,ibar,jbar,kbar
      use corgan_com_M, ONLY : npart, itdim
      implicit real*8 (a-h,o-z)
!
      integer :: ijkctmp(*),iphead(*) 
      real(double) ::     &
          wate(itdim,*),          &          
          x(*),y(*),z(*), ALFA(*),          &          
          tsix(*),tsiy(*),tsiz(*),          &
          bxv(*),byv(*),bzv(*),          &
          xptilde(*),yptilde(*),zptilde(*),          &
          uptilde(*),vptilde(*),wptilde(*),          &
          etax(*),etay(*),etaz(*),          &
          nux(*),nuy(*),nuz(*),          &
          px(0:npart),py(0:npart),pz(0:npart),          &
          up(0:npart),vp(0:npart),wp(0:npart),          &
          pxi(0:npart),peta(0:npart),pzta(0:npart)
!

      logical :: succes
!
      ibp1=ibar+1
      jbp1=jbar+1
      kbp1=kbar+1
!
!
      do 1 n=1,nout
!
      ijk=ijkctmp(n)
      np=iphead(ijk)
!
      pzta(np)=real(kbar)+2.-1.e-3
!
      xptilde(ijk)=px(np)
      yptilde(ijk)=py(np)
      zptilde(ijk)=pz(np)
!
!     uptilde(ijk)=up(np)
!     vptilde(ijk)=vp(np)
!     wptilde(ijk)=wp(np)
!
      ALFA(IJK)=0.0
!
    1 continue
!
      eps=1.e-4
!
!
!     **********************************************
!
      nout_nxt=0
      do 2 n=1,nout
!
      ijk=ijkctmp(n)
      np=iphead(ijk)
!
      xi=pxi(np)
      eta=peta(np)
      zta=pzta(np)
!
      iold=int(pxi(np))
      jold=int(peta(np))
      kold=int(pzta(np))
!
      ijkold=(kold-1)*kwid+(jold-1)*jwid+(iold-1)*iwid+1
!
!
      call map3d_surf(ijk,iwid,jwid,kwid,eps,ifail,    &
                     iphead,itdim,wate,x,y,z,    &
                     xptilde(ijk),yptilde(ijk),zptilde(ijk),    &
                     uptilde(ijk),vptilde(ijk),wptilde(ijk),    &
                     xi,eta,kbp1+1,alfa(ijk))
!
!
      delxi=xi-iold
      deleta=eta-jold
      alfa(ijk)=alfa(ijk)/dt
!
      if(delxi*(1.-delxi).ge.0.      &
          .and.deleta*(1.-deleta).ge.0.      &
          .and.alfa(ijk)*(1.-alfa(ijk)).ge.0.      &
          .and.ifail.eq.0) then
!
      pxi(np)=xi
      peta(np)=eta
      pzta(np)=kbp1+.999999
!
      ll=1
      call parrefl(ll,ijkctmp(n),iphead,npart,itdim,     &
          iwid,jwid,kwid,x,y,z,wate,bxv,byv,bzv,     &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
!
      else
      nout_nxt=nout_nxt+1
      ijkctmp(nout_nxt)=ijkctmp(n)
!
      end if
!
!
  2   continue
!
      nout=nout_nxt
!
!
      return
      end subroutine parloco
!
!     *************************************************************
!
      subroutine trilinp(ncells,ijkcell,itdim,iwid,jwid,kwid,     &
          iphead,pxi,peta,pzta,     &
          wate,bxv,byv,bzv,bxpn,bypn,bzpn)
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      real(double) :: wate(itdim,*),bxv(*),byv(*),bzv(*),     &
          bxpn(*),bypn(*),bzpn(*)
      integer :: ijkcell(*),iphead(*)
!
      real(double) :: pxi(0:*),peta(0:*),pzta(0:*)
!
      do 41 n=1,ncells
!
      ijk=ijkcell(n)
      np=iphead(ijk)
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
      ijknew=(knew-1)*kwid+(jnew-1)*jwid+(inew-1)*iwid+1
!
!
!
      bxpn(ijk)=(wate(ijk,1)*bxv(ijknew+iwid)     &
               +(wate(ijk,2)*bxv(ijknew+iwid+jwid)     &
               +(wate(ijk,3)*bxv(ijknew+jwid)     &
               +(wate(ijk,4)*bxv(ijknew)     &
               +(wate(ijk,5)*bxv(ijknew+iwid+kwid)     &
               +(wate(ijk,6)*bxv(ijknew+iwid+jwid+kwid)     &
               +(wate(ijk,7)*bxv(ijknew+jwid+kwid)     &
               +(wate(ijk,8)*bxv(ijknew+kwid)))))))))
!
      bypn(ijk)=(wate(ijk,1)*byv(ijknew+iwid)      &
               +(wate(ijk,2)*byv(ijknew+iwid+jwid)      &
               +(wate(ijk,3)*byv(ijknew+jwid)      &
               +(wate(ijk,4)*byv(ijknew)      &
               +(wate(ijk,5)*byv(ijknew+iwid+kwid)      &
               +(wate(ijk,6)*byv(ijknew+iwid+jwid+kwid)      &
               +(wate(ijk,7)*byv(ijknew+jwid+kwid)      &
               +(wate(ijk,8)*byv(ijknew+kwid)))))))))
!
      bzpn(ijk)=(wate(ijk,1)*bzv(ijknew+iwid)      &
               +(wate(ijk,2)*bzv(ijknew+iwid+jwid)      &
               +(wate(ijk,3)*bzv(ijknew+jwid)      &
               +(wate(ijk,4)*bzv(ijknew)      &
               +(wate(ijk,5)*bzv(ijknew+iwid+kwid)      &
               +(wate(ijk,6)*bzv(ijknew+iwid+jwid+kwid)      &
               +(wate(ijk,7)*bzv(ijknew+jwid+kwid)      &
               +(wate(ijk,8)*bzv(ijknew+kwid)))))))))
!
   41 continue
!
      return
      end subroutine trilinp
