      subroutine plotinit(ibp2,jbp2,kbp2,name)
      character*80 name
      call gplot('u',celplt,6)
      call libdisp
      return
      end
c
      subroutine plotclose
      call donepl
      call gdone
      return
      end
c
*dk output
      subroutine output
c
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M

      include 'cindex.com'
      include 'numpar.com'
      include 'cophys.com'
c
      character*8 zname
c
c     check if it is time to plot data:
            if(numit.ge.100) go to 1
            if(ncyc.eq.0) go to 1
            if(t+1.e-10.lt.tout) return
            ixto=ixto+1
c           make sure plotting interval is nonzero
c           use the last nonzero plotting interval:
            if(dto(ixto).eq.0.) ixto=ixto-1
            tout=tout+dto(ixto)
c
  1   continue
      call plt3dgrid(ijkvtx,nvtx,iwid,jwid,kwid,
     &     2,ibp2,2,jbp2,kbp2,kbp2,x,y,z)
      call endpl(0)
c
      do jplane=2,jbp2
      call plt3dgrid(ijkvtx,nvtx,iwid,jwid,kwid,
     &     2,ibp2,jplane,jplane,2,kbp2,x,y,z)
      enddo
      call endpl(0)
c
c
c     plot vectors
c
      call plt3dvec(ijkvtx,nvtx,iwid,jwid,kwid,
     &     ibar,jbar,kbar,
     &     2,ibp2,2,jbp2,2,kbp2,x,y,z,bxv,byv,bzv)
c
      call endpl(0)
c
      call plt3dvec(ijkvtx,nvtx,iwid,jwid,kwid,
     &     ibar,jbar,kbar,
     &     2,ibp2,2,jbp2,kbp2,kbp2,x,y,z,bxv,byv,bzv)
c
      call endpl(0)
c
c
      do jplane=2,jbp2
      call plt3dgrid(ijkvtx,nvtx,iwid,jwid,kwid,
     &     2,ibp2,jplane,jplane,2,kbp2,x,y,z)
      call endpl(0)
      enddo
c
      do jplane=2,jbp2
      call plt3dvec(ijkvtx,nvtx,iwid,jwid,kwid,
     &     ibar,jbar,kbar,
     &     2,ibp2,jplane,jplane,2,kbp2,x,y,z,ex,ey,ez)
      call endpl(0)
      enddo
c
      call pltpar3d(ncells,ijkcell,nsp,itdim,iwid,jwid,kwid,
     &     ijkvtx,nvtx,ibar,jbar,kbar,
     &     iphead,iphd2,link,x,y,z,
     &     xmin,xmax,ymin,ymax,zmin,zmax,
     &     ico,pxi,peta,pzta,up,vp,wp,a11,a12,a13,
     &     wate,mask)
c
      call endpl(0)
c
      nh=mod(ncyc,nhst)
c
      if(nh.gt.1) then
      zname='efnrg'
c
      call plthist(1,nh,thistry,efnrg,zname)
c
      zname='ebnrg'
c
      call plthist(1,nh,thistry,ebnrg,zname)
c
      zname='eknrg'
c
      call plthist(1,nh,thistry,eknrg,zname)
c
      zname='charge'
c
      call plthist(1,nh,thistry,charge,zname)
c
      endif
c
      return
      end
c
      subroutine plt3dgrid(ijkvtx,nvtx,iwid,jwid,kwid,
     &     i1,i2,j1,j2,k1,k2,x,y,z)
c
c     a routine to plot a grid in 3 dimensions
c     calls disspla routines
c     note: plots surfaces of constant i, j, or k
c     so that one pair of range indices must be equal, i.e. i1=i2
c
      dimension ijkvtx(*),x(*),y(*),z(*)
c
c     define 2-dimensional plot boundaries
c
      call page(11.,8.5)
      call area2d(10.,7.5)
c
c
c     calculate the range of values in each coordinate direction
c
      lend=ijkvtx(nvtx)
c
      ijkmax=ismax(lend,x,1)
      xmax=x(ijkmax)
c
      ijkmin=ismin(lend,x,1)
      xmin=x(ijkmin)
c
      ijkmax=ismax(lend,y,1)
      ymax=y(ijkmax)
c
      ijkmin=ismin(lend,y,1)
      ymin=y(ijkmin)
c
      ijkmax=ismax(lend,z,1)
      zmax=z(ijkmax)
c
      ijkmin=ismin(lend,z,1)
      zmin=z(ijkmin)
c
c     define workbox
c     with arguments 0.,0.,0. ratio of relative to absolute units
c     same for each axis
c
      call volm3d(0.,0.,0.)
c
c     define origin, step size, maximum axis values
c
      call graf3d(xmin,'scale',xmax,
     &            ymin,'scale',ymax,
     &            zmin,'scale',zmax)
c
      call x3name('x-axis',6)
      call y3name('y-axis',6)
      call z3name('z-axis',6)
c
c     use default viewing pyramid, no calls needed
c
c      default (-1.5*xaxis,-1.5*yaxis,-1.5*zaxis)
c
c     plot grid
c
      if(i1.eq.i2) then
c
c     plot a surface of constant i
c
      do 11 k=k1,k2
      do 11 j=j1,j2-1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i1-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+jwid),y(ijk+jwid),z(ijk+jwid),0)
c
   11 continue
c
      do 12 j=j1,j2
      do 12 k=k1,k2-1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i1-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &     x(ijk+kwid),y(ijk+kwid),z(ijk+kwid),0)
c
   12 continue
c
      endif
c
      if(j1.eq.j2) then
c
c     plot a surface of constant j
c
      do 21 k=k1,k2
      do 21 i=i1,i2-1
c
      ijk=(k-1)*kwid+(j1-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+iwid),y(ijk+iwid),z(ijk+iwid),0)
c
   21 continue
c
      do 22 i=i1,i2
      do 22 k=k1,k2-1
c
      ijk=(k-1)*kwid+(j1-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &     x(ijk+kwid),y(ijk+kwid),z(ijk+kwid),0)
c
   22 continue
c
      endif
c
      if(k1.eq.k2) then
c
c     plot a surface of constant k
c
      do 31 j=j1,j2
      do 31 i=i1,i2-1
c
      ijk=(k1-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+iwid),y(ijk+iwid),z(ijk+iwid),0)
c
   31 continue
c
      do 32 i=i1,i2
      do 32 j=j1,j2-1
c
      ijk=(k1-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &     x(ijk+jwid),y(ijk+jwid),z(ijk+jwid),0)
c
   32 continue
c
      endif
c
      call endgr(0)
c
       return
       end
      subroutine plt3dvec(ijkvtx,nvtx,iwid,jwid,kwid,
     &     ibar,jbar,kbar,
     &     i1,i2,j1,j2,k1,k2,x,y,z,u,v,w)
c
c     a routine to plot a grid in 3 dimensions
c     calls disspla routines
c     note: plots surfaces of constant i, j, or k
c     so that one pair of range indices must be equal, i.e. i1=i2
c
      dimension ijkvtx(*),x(*),y(*),z(*),u(*),v(*),w(*)
c
c     define 2-dimensional plot boundaries
c
      call page(11.,8.5)
      call area2d(10.,7.5)
c
c
c     calculate the range of values in each coordinate direction
c
      lend=ijkvtx(nvtx)
c
      ijkmax=ismax(lend,x,1)
      xmax=x(ijkmax)
c
      ijkmin=ismin(lend,x,1)
      xmin=x(ijkmin)
c
      ijkmax=ismax(lend,y,1)
      ymax=y(ijkmax)
c
      ijkmin=ismin(lend,y,1)
      ymin=y(ijkmin)
c
      ijkmax=ismax(lend,z,1)
      zmax=z(ijkmax)
c
      ijkmin=ismin(lend,z,1)
      zmin=z(ijkmin)
c
c
c     calculate the maximumm vector length
c
      vmax=0.0
      do n=1,nvtx
      ijk=ijkvtx(n)
      vlensq=u(ijk)**2+v(ijk)**2+w(ijk)**2
      vmax=amax1(vmax,vlensq)
      enddo
      vmax=sqrt(vmax)
c
c     calculate average mesh spacing in x,y,z
c
      dxbar=(xmax-xmin)/float(ibar)
      dybar=(ymax-ymin)/float(jbar)
      dzbar=(zmax-zmin)/float(kbar)
c
      rvmax=1./(vmax+1.e-20)
c
      drou=0.75*dxbar*rvmax
      drov=0.75*dybar*rvmax
      drow=0.75*dzbar*rvmax
c
c     define workbox
c     with arguments 0.,0.,0. ratio of relative to absolute units
c     same for each axis
c
      call volm3d(0.,0.,0.)
c
c     define origin, step size, maximum axis values
c
      call graf3d(xmin,'scale',xmax,
     &            ymin,'scale',ymax,
     &            zmin,'scale',zmax)
c
      call x3name('x-axis',6)
      call y3name('y-axis',6)
      call z3name('z-axis',6)
c
c     use default viewing pyramid, no calls needed
c
c      default (-1.5*xaxis,-1.5*yaxis,-1.5*zaxis)
c
c     plot grid
c
c     plot a surface of constant i
c
      do 11 k=k1,k2
      do 11 j=j1,j2
      do 11 i=i1,i2
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &   x(ijk)+u(ijk)*drou,y(ijk)+v(ijk)*drov,z(ijk)+w(ijk)*drow,1101)
c
   11 continue
c
c
c
      call endgr(0)
c
       return
       end
      subroutine pltfram3d(ijkvtx,nvtx,iwid,jwid,kwid,
     &     xmin,xmax,ymin,ymax,zmin,zmax,
     &     i1,i2,j1,j2,k1,k2,x,y,z)
c
c     a routine to plot a grid in 3 dimensions
c     calls disspla routines
c     note: plots surfaces of constant i, j, or k
c     so that one pair of range indices must be equal, i.e. i1=i2
c
      dimension ijkvtx(*),x(*),y(*),z(*)
c
c     define 2-dimensional plot boundaries
c
      call page(11.,8.5)
      call area2d(10.,7.5)
c
c
c     calculate the range of values in each coordinate direction
c
      lend=ijkvtx(nvtx)
c
      ijkmax=ismax(lend,x,1)
      xmax=x(ijkmax)
c
      ijkmin=ismin(lend,x,1)
      xmin=x(ijkmin)
c
      ijkmax=ismax(lend,y,1)
      ymax=y(ijkmax)
c
      ijkmin=ismin(lend,y,1)
      ymin=y(ijkmin)
c
      ijkmax=ismax(lend,z,1)
      zmax=z(ijkmax)
c
      ijkmin=ismin(lend,z,1)
      zmin=z(ijkmin)
c
c     define workbox
c     with arguments 0.,0.,0. ratio of relative to absolute units
c     same for each axis
c
      call volm3d(0.,0.,0.)
c
c     define origin, step size, maximum axis values
c
      call graf3d(xmin,'scale',xmax,
     &            ymin,'scale',ymax,
     &            zmin,'scale',zmax)
c
      call x3name('x-axis',6)
      call y3name('y-axis',6)
      call z3name('z-axis',6)
c
c     use default viewing pyramid, no calls needed
c
c      default (-1.5*xaxis,-1.5*yaxis,-1.5*zaxis)
c
c     plot grid
c
      do 11 i=i1,i2-1
c
      j=j1
      k=k1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+iwid),y(ijk+iwid),z(ijk+iwid),0)
c
      j=j2
      k=k1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+iwid),y(ijk+iwid),z(ijk+iwid),0)
c
      j=j2
      k=k2
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+iwid),y(ijk+iwid),z(ijk+iwid),0)
c
      j=j1
      k=k2
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+iwid),y(ijk+iwid),z(ijk+iwid),0)
c
   11 continue
c
      do 21 j=j1,j2-1
c
      i=i1
      k=k1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+jwid),y(ijk+jwid),z(ijk+jwid),0)
c
      i=i2
      k=k1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+jwid),y(ijk+jwid),z(ijk+jwid),0)
c
      i=i2
      k=k2
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+jwid),y(ijk+jwid),z(ijk+jwid),0)
c
      i=i2
      k=k1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+jwid),y(ijk+jwid),z(ijk+jwid),0)
c
   21 continue
c
      do 31 k=k1,k2-1
c
      i=i1
      j=j1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+kwid),y(ijk+kwid),z(ijk+kwid),0)
c
      i=i2
      j=j1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+kwid),y(ijk+kwid),z(ijk+kwid),0)
c
      i=i2
      j=j2
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+kwid),y(ijk+kwid),z(ijk+kwid),0)
c
      i=i2
      j=j1
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
      call rlvec3(x(ijk),y(ijk),z(ijk),
     &           x(ijk+kwid),y(ijk+kwid),z(ijk+kwid),0)
c
   31 continue
c
c
c
      call endgr(0)
c
       return
       end
      subroutine plthist(n1,n2,t,y,zname)
c
c     a routine to plot history variables
c
      character*8 zname
c
      dimension t(*),y(*)
c
      call page(11.,8.5)
      call nobrdr
      call grace(3.0)
      call area2d(10.,7.5)
c
      nmax=ismax(n2,t,1)
      tmax=t(nmax)
      nmin=ismin(n2,t,1)
      tmin=t(nmin)
c
      nmax=ismax(n2,y,1)
      ymax=y(nmax)
      nmin=ismin(n2,y,1)
      ymin=y(nmin)
c
      call xname("t,",1)
      call yname(zname,8)
c
      call graf(tmin,'scale',tmax,
     &          ymin,'scale',ymax)
c
      call curve(t,y,n2,0)
c
      call endpl(0)
c
      return
      end
      subroutine pltpar3d(ncells,ijkcell,nsp,itdim,iwid,jwid,kwid,
     &     ijkvtx,nvtx,ibar,jbar,kbar,
     &     iphead,iphd2,link,x,y,z,
     &     xmin,xmax,ymin,ymax,zmin,zmax,
     &     ico,pxi,peta,pzta,up,vp,wp,xp,yp,zp,
     &     wate,mask)
c
      dimension ico(0:*),link(0:*),
     &          pxi(0:*),peta(0:*),pzta(0:*),
     &          up(0:*),vp(0:*),wp(0:*),
     &       x(*),y(*),z(*),
     &     xp(*),yp(*),zp(*),
     &          wate(itdim,27),mask(*),
     &     ijkcell(*),iphead(*),iphd2(*)
c
      dimension ijkvtx(*)
c
c
      logical nomore
c
        real nu
c
      call pltfram3d(ijkvtx,nvtx,iwid,jwid,kwid,
     &     xmin,xmax,ymin,ymax,zmin,zmax,
     &     2,ibar+2,2,jbar+2,2,kbar+2,x,y,z)
c
c     a routine to plot the particles
c
c
c
c     define 2-dimensional plot boundaries
c
      call page(11.,8.5)
      call area2d(10.,7.5)
c
c
c     define workbox
c     with arguments 0.,0.,0. ratio of relative to absolute units
c     same for each axis
c
      call volm3d(0.,0.,0.)
c
c     define origin, step size, maximum axis values
c
      call graf3d(xmin,'scale',xmax,
     &            ymin,'scale',ymax,
     &            zmin,'scale',zmax)
c
c
c     scale particle markeres
c
      call sclpic(0.15)
c
c
      do 20 n=1,ncells
c
        iphd2(ijkcell(n))=0
   20   continue
c
c
c
    1   continue
        nomore=.true.
c
cdir$ ivdep
c
        newcell=0
        do 100 n=1,ncells
        if(iphead(ijkcell(n)).ne.0) then
        newcell=newcell+1
        endif
  100   continue
c
        if(newcell.ne.0) then
        nomore=.false.
c
cdir$ ivdep
c
        do 200 n=1,ncells
c
        ijk=ijkcell(n)
        np=iphead(ijk)
c
        k=ifix(pzta(np))
        j=ifix(peta(np))
        i=ifix(pxi(np))
c
        the=pxi(np)-i
        zeta=peta(np)-j
        nu=pzta(np)-k
c
      wi=1.-the
      wip=the
c
      wj=1.-zeta
      wjp=zeta
c
      wk=1.-nu
      wkp=nu
c
c     k-plane
c
      wate(ijk,1)=wip*wj*wk
      wate(ijk,2)=wip*wjp*wk
      wate(ijk,3)=wi*wjp*wk
      wate(ijk,4)=wi*wj*wk
c
      wate(ijk,5)=wip*wj*wkp
      wate(ijk,6)=wip*wjp*wkp
      wate(ijk,7)=wi*wjp*wkp
      wate(ijk,8)=wi*wj*wkp
c
  200 continue
c
cdir$ ivdep
c
      do 205 n=1,ncells
c
      ijk=ijkcell(n)
c
      xp(ijk)=(wate(ijk,1)*x(ijk+iwid)
     &          +(wate(ijk,2)*x(ijk+iwid+jwid)
     &          +(wate(ijk,3)*x(ijk+jwid)
     &          +(wate(ijk,4)*x(ijk)
     &          +(wate(ijk,5)*x(ijk+iwid+kwid)
     &          +(wate(ijk,6)*x(ijk+iwid+jwid+kwid)
     &          +(wate(ijk,7)*x(ijk+jwid+kwid)
     &          +(wate(ijk,8)*x(ijk+kwid)))))))))
c
      yp(ijk)=(wate(ijk,1)*y(ijk+iwid)
     &          +(wate(ijk,2)*y(ijk+iwid+jwid)
     &          +(wate(ijk,3)*y(ijk+jwid)
     &          +(wate(ijk,4)*y(ijk)
     &          +(wate(ijk,5)*y(ijk+iwid+kwid)
     &          +(wate(ijk,6)*y(ijk+iwid+jwid+kwid)
     &          +(wate(ijk,7)*y(ijk+jwid+kwid)
     &          +(wate(ijk,8)*y(ijk+kwid)))))))))
c
      zp(ijk)=(wate(ijk,1)*z(ijk+iwid)
     &          +(wate(ijk,2)*z(ijk+iwid+jwid)
     &          +(wate(ijk,3)*z(ijk+jwid)
     &          +(wate(ijk,4)*z(ijk)
     &          +(wate(ijk,5)*z(ijk+iwid+kwid)
     &          +(wate(ijk,6)*z(ijk+iwid+jwid+kwid)
     &          +(wate(ijk,7)*z(ijk+jwid+kwid)
     &          +(wate(ijk,8)*z(ijk+kwid)))))))))
c
  205 continue
c
c     plot particles
c
      do 210 n=1,ncells
c
      ijk=ijkcell(n)
      if(iphead(ijk).ne.0) then
c
      is=ico(iphead(ijk))
c
      call marker(is)
      call curv3d(xp(ijk),yp(ijk),zp(ijk),1,-1)
c
      endif
c
  210 continue
c
c
c
c     ************************************************
c
c
c
      do 500 n=1,ncells
      ijk=ijkcell(n)
      np=iphead(ijk)
      if(np.gt.0) then
      iphead(ijk)=link(np)
      link(np)=iphd2(ijk)
      iphd2(ijk)=np
c
      endif
  500 continue
      endif
c
      if(.not.nomore) go to 1
c
      do 900 n=1,ncells
c
      ijk=ijkcell(n)
c
      iphead(ijk)=iphd2(ijk)
c
      iphd2(ijk)=0
c
  900 continue
c
      return
      end
