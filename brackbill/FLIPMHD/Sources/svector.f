      subroutine bcperv(nxp,nyp,nzp,iwid,jwid,kwid,    &
          sdlt,cdlt,    &
          ax,ay,az)
!
!      a routine to impose periodic boundary conditions
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      real(double) :: ax(*),ay(*),az(*)
!
      one=1.0
!
!     set reference vectors for i=1 and i=nxp
!
      do 1 k=2,nzp
!
      do 11 j=1,nyp
!
      ijk1=1+(j-1)*jwid+(k-1)*kwid
      ijknxp=1+(nxp-1)*iwid+(j-1)*jwid+(k-1)*kwid
!
      ax(ijk1)=ax(ijknxp-iwid)
      ax(ijknxp)=ax(ijk1+iwid)
!
      ay(ijk1)=ay(ijknxp-iwid)
      ay(ijknxp)=ay(ijk1+iwid)
!
      az(ijk1)=az(ijknxp-iwid)
      az(ijknxp)=az(ijk1+iwid)
!
   11 continue
!
!     set reference vectors for j=1 and j=nyp
!
      do 12 i=1,nxp
!
      ijk1=1+(i-1)*iwid+(k-1)*kwid
      ijknyp=1+(i-1)*iwid+(nyp-1)*jwid+(k-1)*kwid
!
      ax(ijk1)=ax(ijknyp-jwid)*cdlt+ay(ijknyp-jwid)*sdlt
      ay(ijk1)=-ax(ijknyp-jwid)*sdlt+ay(ijknyp-jwid)*cdlt
      az(ijk1)=az(ijknyp-jwid)
!
      ax(ijknyp)=ax(ijk1+jwid)*cdlt-ay(ijk1+jwid)*sdlt
      ay(ijknyp)=ax(ijk1+jwid)*sdlt+ay(ijk1+jwid)*cdlt
      az(ijknyp)=az(ijk1+jwid)
!
   12 continue
!
    1 continue
!
      return
      end subroutine bcperv
!     ***********************************************************************
      subroutine bvec(nxp,nyp,nzp,ijkcell,ncells,   &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,   &
          rmaj,RWALL,strait,toroid,   &
          q0,bzi,PSI,   &
          x,y,z,   &
          bx,by,bz)
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      integer :: ijkcell(*)
      real(double) :: x(*),y(*),z(*),bx(*),by(*),bz(*),     &
          tsix(*),tsiy(*),tsiz(*),     &
          etax(*),etay(*),etaz(*),     &
          nux(*),nuy(*),nuz(*)
!
      ibar=nxp-2
      jbar=nyp-2
      kbar=nzp-2
      pi=acos(-1.)
!
      do 1 n=1,ncells
!
      ijk=ijkcell(n)
!
         yv=y(ijk)*toroid
         ws= toroid*sqrt(x(ijk)**2+yv**2)-rmaj+strait*x(ijk)
!
!     wsr is the minor radius
!
         wsr=sqrt(ws**2+z(ijk)**2)
!
!     wsmr is the major radius
!
      wsmr=sqrt(x(ijk)**2+y(ijk)**2+z(ijk)**2)
!     Ref: Bauer, Betancourt, and Garabedian
!     "A Computational Method in Plasma Physics"
!
!
!
!     calculate toroidal angle
!
      one=1.
      cosphi=x(ijk)/wsmr
      cosphi=min(one,cosphi)
      cosphi=max(-one,cosphi)
      phi=acos(cosphi)
!
      wsq0=q0*real(ibar)*cos(2.*phi)*wsr/real(jbar)
      denom=2.*pi/(real(ibar)*real(kbar))
!
      bx(ijk)=bzi*(nuy(ijk)*(wsq0*etaz(ijk)+tsiz(ijk))    &
                  -nuz(ijk)*(wsq0*etay(ijk)+tsiy(ijk)))*denom    &
          *(rmaj/wsmr+strait)*wsr
!
      by(ijk)=bzi*(nuz(ijk)*(wsq0*etax(ijk)+tsix(ijk))     &
                  -nux(ijk)*(wsq0*etaz(ijk)+tsiz(ijk)))*denom     &
          *(rmaj/wsmr+strait)*wsr
!
      bz(ijk)=bzi*(nux(ijk)*(wsq0*etay(ijk)+tsiy(ijk))     &
                  -nuy(ijk)*(wsq0*etax(ijk)+tsix(ijk)))*denom     &
          *(rmaj/wsmr+strait)*wsr
!
      rbmag=1./(sqrt(bx(ijk)**2+by(ijk)**2+bz(ijk)**2)+1.e-20)
!
      bx(ijk)=bx(ijk)*rbmag
      by(ijk)=by(ijk)*rbmag
      bz(ijk)=bz(ijk)*rbmag
!
    1 continue
!
      return
      end subroutine bvec
!     *****************************************************************
      subroutine fluxsurf(i1,i2,j1,j2,k1,k2,iwid,jwid,kwid,    &
          x,y,z,rwall,rmaj,toroid,strait,delta,    &
          psi)
!
!     a routine to define a flux function to drive the grid generator
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      real(double) :: x(*),y(*),z(*),psi(*)
!
      do 1 k=k1,k2
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      ijk=1+(i-1)*iwid+(j-1)*jwid+(k-1)*kwid
!
!     calculate toroidal coordinates from cartesian coordinates
!
      ws=toroid*sqrt(x(ijk)**2+y(ijk)**2)-rmaj+strait*x(ijk)
!
      rminor=sqrt(ws**2+z(ijk)**2)
!
      xv=x(ijk)*toroid+strait
      yv=y(ijk)*toroid
!
      rmr=1./sqrt(xv**2+yv**2)
!
      sphi=yv*rmr
      cphi=xv*rmr
!
      stheta=ws/(rminor+1.e-35)
      ctheta=z(ijk)/(rminor+1.e-35)
!
!      deltar=delta*rminor*(rwall-rminor)/rwall**2
      deltar=delta*rminor**3*(rwall-rminor)/rwall**4*256./27.
!
!     multiple angle formulae
!
!     s2theta=2.*stheta*ctheta
      c2theta=2.*ctheta**2-1.
!
!     s3theta=3.*stheta-4.*stheta**2
      c3theta=4.*ctheta**3-3.*ctheta
!
      psi(ijk)=rminor*(1.-deltar*c2theta)
!
    1 continue
!
      return
      end subroutine fluxsurf
!     *******************************************************************
      subroutine refvec(i1,i2,j1,j2,k1,k2,iwid,jwid,kwid,    &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,    &
          x,y,z,    &
          TOROID,STRAIT,RWALL,RMAJ,    &
          psi,cgx,cgy,cgz,AGX,AGY,AGZ,BGX,BGY,BGZ)
!
!     a routine to construct the unit normal to surfaces of constant flux
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      real(double) :: tsix(*),tsiy(*),tsiz(*),     &
          etax(*),etay(*),etaz(*),     &
          nux(*),nuy(*),nuz(*),     &
          x(*),y(*),z(*),     &
          psi(*),cgx(*),cgy(*),cgz(*),     &
          agx(*),agy(*),agz(*),     &
          bgx(*),bgy(*),bgz(*)
!
      do 1 k=k1,k2-1
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      ijk=1+(i-1)*iwid+(j-1)*jwid+(k-1)*kwid
!
      psi1=0.5*(psi(ijk+iwid)-psi(ijk-iwid))
      psi2=0.5*(psi(ijk+jwid)-psi(ijk-jwid))
      psi3=0.5*(psi(ijk+kwid)-psi(ijk-kwid))
!
      cgx(ijk)=psi1*tsix(ijk)+psi2*etax(ijk)+psi3*nux(ijk)
      cgy(ijk)=psi1*tsiy(ijk)+psi2*etay(ijk)+psi3*nuy(ijk)
      cgz(ijk)=psi1*tsiz(ijk)+psi2*etaz(ijk)+psi3*nuz(ijk)
!
      rcg=1./(sqrt(cgx(ijk)**2+cgy(ijk)**2+cgz(ijk)**2)+1.e-30)
!c
      cgx(ijk)=cgx(ijk)*rcg
      cgy(ijk)=cgy(ijk)*rcg
      cgz(ijk)=cgz(ijk)*rcg
!
!     construct a vector that is orthogonal to b (calculated in bvec)
!     and c, which is normal to flux surfaces
!
      agx(ijk)=bgy(ijk)*cgz(ijk)-bgz(ijk)*cgy(ijk)
      agy(ijk)=bgz(ijk)*cgx(ijk)-bgx(ijk)*cgz(ijk)
      agz(ijk)=bgx(ijk)*cgy(ijk)-bgy(ijk)*cgx(ijk)
!
      rag=1./(sqrt(agx(ijk)**2+agy(ijk)**2+agz(ijk)**2)+1.e-20)
!
      agx(ijk)=agx(ijk)*rag
      agy(ijk)=agy(ijk)*rag
      agz(ijk)=agz(ijk)*rag
!
    1 continue
!
!     calculate reference vectors on the axis and the outer wall
!
      do 2 j=j1,j2
      do 2 i=i1,i2
!
      ijk=1+(i-1)*iwid+(j-1)*jwid+(k2-1)*kwid
!
      ws=toroid*sqrt(x(ijk)**2+y(ijk)**2)-rmaj+strait*x(ijk)
!
      xv=x(ijk)*toroid+strait
      yv=y(ijk)*toroid
!
      rmr=1./sqrt(xv**2+yv**2)
!
      sphi=yv*rmr
      cphi=xv*rmr
      stheta=ws/rwall
      ctheta=z(ijk)/rwall
!
      cgx(ijk)=sphi*ctheta
      cgy(ijk)=cphi*ctheta
      cgz(ijk)=stheta
      agx(ijk)=agx(ijk-kwid)
      agy(ijk)=agy(ijk-kwid)
      agz(ijk)=agz(ijk-kwid)
!
      bgx(ijk)=bgx(ijk-kwid)
      bgy(ijk)=bgy(ijk-kwid)
      bgz(ijk)=bgz(ijk-kwid)
!
!      cgx(ijk)=cgx(ijk-kwid)
!      cgy(ijk)=cgy(ijk-kwid)
!      cgz(ijk)=cgz(ijk-kwid)
!
!
!     set unit vectors on the axis equal to unit vectors at the wall
!
      ijka=1+(i-1)*iwid+(j-1)*jwid
!
      cgx(ijka)=cgx(ijk)
      cgy(ijka)=cgy(ijk)
      cgz(ijka)=cgz(ijk)
!
!
    2 continue
!
!
      return
      end subroutine refvec
