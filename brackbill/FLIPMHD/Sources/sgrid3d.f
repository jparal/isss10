      subroutine adaptw(nxp,nyp,nzp,i1,i2,j1,j2,k1,k2,     &
          toroid,strait,rmaj,     &
          x,y,z,wgrid)
!
!     a routine to calculate the weight function for a mesh in toroidal
!     geometry
!
!     the weight function is the Jacobian for toroidal coordinates
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      real(double) ::  x(nxp,nyp,*),y(nxp,nyp,*),z(nxp,nyp,*),     &
              wgrid(nxp,nyp,*)
!
      one=1.0
!
      do 1 k=k1,k2
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      ws=toroid*sqrt(x(i,j,k)**2+y(i,j,k)**2)-rmaj+strait*x(i,j,k)
      rminor=sqrt(ws**2+z(i,j,k)**2)
      rmajor=toroid*sqrt(x(i,j,k)**2+y(i,j,k)**2)+strait
!
      sinphi=y(i,j,k)/rmajor
      cosphi=x(i,j,k)/rmajor
!
      stheta=ws/(rminor+1.e-30)
      ctheta=z(i,j,k)/(rminor+1.e-30)
      ctheta=min(ctheta,one)
      ctheta=max(ctheta,-one)
      cosphi=min(cosphi,one)
      cosphi=max(cosphi,-one)
!
!
      phi=acos(cosphi)
      theta=acos(ctheta)
!
!      wgrid(i,j,k)=wgrid(i,j,k)*(1.0+1.e3*cosphi**4)
!      wgrid(i,j,k)=wgrid(i,j,k)
!     &     *(1.+025.*exp(-100.*(rminor-0.5)**2))
!      wgrid(i,j,k)=wgrid(i,j,k)
!      wgrid(i,j,k)=wgrid(i,j,k)*(1.+0.75*cos(2.*theta-3.*phi))
!
    1 continue
!
      return
      end subroutine adaptw
!     ********************************************************************************
      subroutine contras(iwid,jwid,kwid,ncells,ijkcell,    &
          x,y,z,    &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz)
!
!     a routine to calculate the contravariant mesh
!     vectors at the vertices of the mesh
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
     real(double) ::  x(*),y(*),z(*), & 
           tsix(*),tsiy(*),tsiz(*),    &
          etax(*),etay(*),etaz(*),    &
          nux(*),nuy(*),nuz(*)
!
      real(double) ::  jacob
!
!
      do 10 n=1,ncells
!
      ijk=ijkcell(n)
      ipjk=ijk+iwid
      ijpk=ijk+jwid
      imjk=ijk-iwid
      ijmk=ijk-jwid
      ijkp=ijk+kwid
      ijkm=ijk-kwid
!
!     compute covariant base vectors
!
      x1=0.5*(x(ipjk)-x(imjk))
      x2=0.5*(x(ijpk)-x(ijmk))
      x3=0.5*(x(ijkp)-x(ijkm))
!
      y1=0.5*(y(ipjk)-y(imjk))
      y2=0.5*(y(ijpk)-y(ijmk))
      y3=0.5*(y(ijkp)-y(ijkm))
!
      z1=0.5*(z(ipjk)-z(imjk))
      z2=0.5*(z(ijpk)-z(ijmk))
      z3=0.5*(z(ijkp)-z(ijkm))
!
!     compute determinant of the metric tensor
!
      jacob=x1*(y2*z3-y3*z2)   &
         +y1*(z2*x3-z3*x2)     &
         +z1*(x2*y3-x3*y2)
!
      rjacob=1./jacob
!
!     compute contravariant base vectors
!
      tsix(ijk)=(y2*z3-y3*z2)*rjacob
      tsiy(ijk)=(z2*x3-z3*x2)*rjacob
      tsiz(ijk)=(x2*y3-x3*y2)*rjacob
!
      etax(ijk)=(y3*z1-y1*z3)*rjacob
      etay(ijk)=(z3*x1-z1*x3)*rjacob
      etaz(ijk)=(x3*y1-x1*y3)*rjacob
!
      nux(ijk)=(y1*z2-y2*z1)*rjacob
      nuy(ijk)=(z1*x2-z2*x1)*rjacob
      nuz(ijk)=(x1*y2-x2*y1)*rjacob
!
!
   10 continue
!
!
      return
      end subroutine contras
!     ******************************************************************
      subroutine diagnos(ijkcell,ncells,      &
          ax,ay,az,bx,by,bz,cx,cy,cz,      &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,      &
          divaa1,divbb2,divcc3,diva,divb,divc,      &
          adottsi,bdoteta,cdotnu)
!
!    a routine to check minimization
!
      implicit real*8 (a-h,o-z)
!
      dimension ax(*),ay(*),az(*),    &
               bx(*),by(*),bz(*),    &
               cx(*),cy(*),cz(*),    &
               tsix(*),tsiy(*),tsiz(*),    &
               etax(*),etay(*),etaz(*),    &
               nux(*),nuy(*),nuz(*),    &
          divaa1(*),divbb2(*),divcc3(*),    &
          ijkcell(*)
!
      adottsi=0.0
      bdoteta=0.0
      cdotnu=0.0
!
      diva=0.0
      divb=0.0
      divc=0.0
!
      do 1 n=1,ncells
      ijk=ijkcell(n)
!
      adottsi=adottsi     &
          +(ax(ijk)*tsix(ijk)+ay(ijk)*tsiy(ijk)+az(ijk)*tsiz(ijk))**2     &
         /(tsix(ijk)**2+tsiy(ijk)**2+tsiz(ijk)**2)
!
      bdoteta=bdoteta     &
          +(bx(ijk)*etax(ijk)+by(ijk)*etay(ijk)+bz(ijk)*etaz(ijk))**2     &
          /(etax(ijk)**2+etay(ijk)**2+etaz(ijk)**2)
!
      cdotnu=cdotnu       &
          +(cx(ijk)*nux(ijk)+cy(ijk)*nuy(ijk)+cz(ijk)*nuz(ijk))**2      &
          /(nux(ijk)**2+nuy(ijk)**2+nuz(ijk)**2)
!
      diva=diva+divaa1(ijk)**2
      divb=divb+divbb2(ijk)**2
      divc=divc+divcc3(ijk)**2
!
    1 continue
!
      rncells=1./real(ncells)
      adottsi=adottsi*rncells
      bdoteta=bdoteta*rncells
      cdotnu=cdotnu*rncells
!
      return
      end subroutine diagnos
!     ***********************************************************************
      subroutine divaa(iwid,jwid,kwid,ncells,ijkcell,      &
          ax,ay,az,      &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,      &
          projx,projy,projz,divaa1)
!
!     calculate the contravariant components of divaa
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      integer :: & 
          ijkcell(*)
      real(double) :: divaa1(*),      &
          ax(*),ay(*),az(*),      &
          projx(*),projy(*),projz(*),      &
          tsix(*),tsiy(*),tsiz(*),      &
          etax(*),etay(*),etaz(*),      &
          nux(*),nuy(*),nuz(*)      
!
!
!
!
!dir$ ivdep
      do 10 n=1,ncells
!
      ijk=ijkcell(n)
!
      ax1=0.5*(ax(ijk+iwid)-ax(ijk-iwid))
      ax2=0.5*(ax(ijk+jwid)-ax(ijk-jwid))
      ax3=0.5*(ax(ijk+kwid)-ax(ijk-kwid))
!
      ay1=0.5*(ay(ijk+iwid)-ay(ijk-iwid))
      ay2=0.5*(ay(ijk+jwid)-ay(ijk-jwid))
      ay3=0.5*(ay(ijk+kwid)-ay(ijk-kwid))
!
      az1=0.5*(az(ijk+iwid)-az(ijk-iwid))
      az2=0.5*(az(ijk+jwid)-az(ijk-jwid))
      az3=0.5*(az(ijk+kwid)-az(ijk-kwid))
!
      adottsi=ax(ijk)*tsix(ijk)+(ay(ijk)*tsiy(ijk)+(az(ijk)*tsiz(ijk)))
!
      adoteta=ax(ijk)*etax(ijk)+(ay(ijk)*etay(ijk)+(az(ijk)*etaz(ijk)))
!
      adotnu=ax(ijk)*nux(ijk)+(ay(ijk)*nuy(ijk)+(az(ijk)*nuz(ijk)))
!
      agradx=adottsi*ax1+adoteta*ax2+adotnu*ax3
      agrady=adottsi*ay1+adoteta*ay2+adotnu*ay3
      agradz=adottsi*az1+adoteta*az2+adotnu*az3
!
      diva=tsix(ijk)*ax1+tsiy(ijk)*ay1+tsiz(ijk)*az1      &
         +etax(ijk)*ax2+etay(ijk)*ay2+etaz(ijk)*az2      &
         +nux(ijk)*ax3+nuy(ijk)*ay3+nuz(ijk)*az3
!
      adotproj=projx(ijk)*ax(ijk)+projy(ijk)*ay(ijk)+projz(ijk)*az(ijk)
!
      a1dota=ax1*ax(ijk)+(ay1*ay(ijk)+(az1*az(ijk)))
      a2dota=ax2*ax(ijk)+(ay2*ay(ijk)+(az2*az(ijk)))
      a3dota=ax3*ax(ijk)+(ay3*ay(ijk)+(az3*az(ijk)))
!
      gradax=tsix(ijk)*a1dota+etax(ijk)*a2dota+nux(ijk)*a3dota
      graday=tsiy(ijk)*a1dota+etay(ijk)*a2dota+nuy(ijk)*a3dota
      gradaz=tsiz(ijk)*a1dota+etaz(ijk)*a2dota+nuz(ijk)*a3dota
!
!
      divaa1(ijk)=projx(ijk)*(2.*gradax-agradx)      &
                +projy(ijk)*(2.*graday-agrady)      &
                +projz(ijk)*(2.*gradaz-agradz)      &
                -diva*adotproj
!
   10 continue
!
      return
      end subroutine divaa
!     *****************************************************************************
      subroutine gridsolv(x,y,z,w,     &
          AX,AY,AZ,LAMA,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          PROJX,PROJY,PROJZ,     &
          DIVAA1,     &
          iwid,jwid,kwid,ijkcell,ncells,     &
          DX1)
!
!     a mesh generator based on
!     Winslow's variable diffusion
!     functional
!
      use vast_kind_param, ONLY: double
      implicit real*8 (a-h,o-z)
!
      real(double) ::     &
             y(*),     &
             z(*),     &
             w(*),     &
          ax(*),ay(*),az(*),     &
          divaa1(*),     &
          tsix(*),tsiy(*),tsiz(*),     &
          etax(*),etay(*),etaz(*),     &
          nux(*),nuy(*),nuz(*),     &
          projx(*),projy(*),projz(*),     &
          ijkcell(*),DX1(*)
!
!
!
!dir$ ivdep
      do 10 n=1,ncells
!
      ijk=ijkcell(n)
      ipjk=ijk+iwid
      ipjpk=ijk+iwid+jwid
      ijpk=ijk+jwid
      imjpk=ijk-iwid+jwid
      imjk=ijk-iwid
      imjmk=ijk-iwid-jwid
      ijmk=ijk-jwid
      ipjmk=ijk+iwid-jwid
!
      ijkp=ijk+kwid
      ipjkp=ijk+iwid+kwid
      ipjpkp=ijk+iwid+jwid+kwid
      ijpkp=ijk+jwid+kwid
      imjpkp=ijk-iwid+jwid+kwid
      imjkp=ijk-iwid+kwid
      imjmkp=ijk-iwid-jwid+kwid
      ijmkp=ijk-jwid+kwid
      ipjmkp=ijk+iwid-jwid+kwid
!
      ijkm=ijk-kwid
      ipjkm=ijk+iwid-kwid
      ipjpkm=ijk+iwid+jwid-kwid
      ijpkm=ijk+jwid-kwid
      imjpkm=ijk-iwid+jwid-kwid
      imjkm=ijk-iwid-kwid
      imjmkm=ijk-iwid-jwid-kwid
      ijmkm=ijk-jwid-kwid
      ipjmkm=ijk+iwid-jwid-kwid
!
!     compute covariant base vectors
!
!
!
!
!     compute covariant base vectors
!
      x1=0.5*(x(ipjk)-x(imjk))
      x2=0.5*(x(ijpk)-x(ijmk))
      x3=0.5*(x(ijkp)-x(ijkm))
!
      y1=0.5*(y(ipjk)-y(imjk))
      y2=0.5*(y(ijpk)-y(ijmk))
      y3=0.5*(y(ijkp)-y(ijkm))
!
      z1=0.5*(z(ipjk)-z(imjk))
      z2=0.5*(z(ijpk)-z(ijmk))
      z3=0.5*(z(ijkp)-z(ijkm))
!
      w1=(w(ipjk)-w(imjk))/(w(ipjk)+w(imjk)+1.e-10)
      w2=(w(ijpk)-w(ijmk))/(w(ijpk)+w(ijmk)+1.e-10)
      w3=(w(ijkp)-w(ijkm))/(w(ijkp)+w(ijkm)+1.e-10)
!
!     calculate elements of the metric tensor
!
      g11=tsix(ijk)**2+tsiy(ijk)**2+tsiz(ijk)**2
      g12=tsix(ijk)*etax(ijk)+tsiy(ijk)*etay(ijk)+tsiz(ijk)*etaz(ijk)
      g13=tsix(ijk)*nux(ijk)+tsiy(ijk)*nuy(ijk)+tsiz(ijk)*nuz(ijk)
      g22=etax(ijk)**2+etay(ijk)**2+etaz(ijk)**2
      g23=etax(ijk)*nux(ijk)+etay(ijk)*nuy(ijk)+etaz(ijk)*nuz(ijk)
      g33=nux(ijk)**2+nuy(ijk)**2+nuz(ijk)**2
!
!     calculate contravariant components of reference vectors
!
      adottsi=ax(ijk)*tsix(ijk)+(ay(ijk)*tsiy(ijk)+(az(ijk)*tsiz(ijk)))
!
      adoteta=ax(ijk)*etax(ijk)+(ay(ijk)*etay(ijk)+(az(ijk)*etaz(ijk)))
!
      adotnu=ax(ijk)*nux(ijk)+(ay(ijk)*nuy(ijk)+(az(ijk)*nuz(ijk)))
      asq=ax(ijk)**2+ay(ijk)**2+az(ijk)**2
!
!
!     calculate elements of directional metric
!
      ga11=(1.-lama)*g11+lama*(asq*g11-adottsi**2)
      ga12=(1.-lama)*g12+lama*(asq*g12-adottsi*adoteta)
      ga13=(1.-lama)*g13+lama*(asq*g13-adottsi*adotnu)
      ga22=(1.-lama)*g22+lama*(asq*g22-adoteta**2)
      ga23=(1.-lama)*g23+lama*(asq*g23-adoteta*adotnu)
      ga33=(1.-lama)*g33+lama*(asq*g33-adotnu**2)
!
!     compute derivatives
!
      x11=x(ipjk)-2.*x(ijk)+x(imjk)
      x12=0.25*(x(ipjpk)-x(imjpk)+x(imjmk)-x(ipjmk))
      x13=0.25*(x(ipjkp)-x(imjkp)+x(imjkm)-x(ipjkm))
      x22=x(ijpk)-2.*x(ijk)+x(ijmk)
      x23=0.25*(x(ijpkp)-x(ijmkp)+x(ijmkm)-x(ijpkm))
      x33=x(ijkp)-2.*x(ijk)+x(ijkm)
!
      y11=y(ipjk)-2.*y(ijk)+y(imjk)
      y12=0.25*(y(ipjpk)-y(imjpk)+y(imjmk)-y(ipjmk))
      y13=0.25*(y(ipjkp)-y(imjkp)+y(imjkm)-y(ipjkm))
      y22=y(ijpk)-2.*y(ijk)+y(ijmk)
      y23=0.25*(y(ijpkp)-y(ijmkp)+y(ijmkm)-y(ijpkm))
      y33=y(ijkp)-2.*y(ijk)+y(ijkm)
!
      z11=z(ipjk)-2.*z(ijk)+z(imjk)
      z12=0.25*(z(ipjpk)-z(imjpk)+z(imjmk)-z(ipjmk))
      z13=0.25*(z(ipjkp)-z(imjkp)+z(imjkm)-z(ipjkm))
      z22=z(ijpk)-2.*z(ijk)+z(ijmk)
      z23=0.25*(z(ijpkp)-z(ijmkp)+z(ijmkm)-z(ijpkm))
      z33=z(ijkp)-2.*z(ijk)+z(ijkm)
!
      rax=ga11*(x11+w1*x1)     &
          +ga12*(2.*x12+(w1*x2+w2*x1))     &
          +ga13*(2.*x13+(w1*x3+w3*x1))     &
          +ga22*(x22+w2*x2)     &
          +ga23*(2.*x23+(w2*x3+w3*x2))     &
          +ga33*(x33+w3*x3)
!
      ray=ga11*(y11+w1*y1)     &
          +ga12*(2.*y12+(w1*y2+w2*y1))     &
          +ga13*(2.*y13+(w1*y3+w3*y1))     &
          +ga22*(y22+w2*y2)     &
          +ga23*(2.*y23+(w2*y3+w3*y2))     &
          +ga33*(y33+w3*y3)
!
      raz=ga11*(z11+w1*z1)     &
          +ga12*(2.*z12+(w1*z2+w2*z1))     &
          +ga13*(2.*z13+(w1*z3+w3*z1))     &
          +ga22*(z22+w2*z2)     &
          +ga23*(2.*z23+(w2*z3+w3*z2))     &
          +ga33*(z33+w3*z3)
!
!     assemble the contravariant components of r
!
      r1=rax*projx(ijk)+(ray*projy(ijk)+(raz*projz(ijk)))    &
         -lama*divaa1(ijk)
!
      ws1=-1./(2.*(ga11+ga22+ga33))
!
      dx1(ijk)=-r1*ws1
!
   10 continue
!
!
      return
      end subroutine gridsolv
!     ************************************************************************
      subroutine meshgen(nxp,nyp,nzp,x,y,z,w,vol,dt,taug,numit,     &   
          DX1,DX2,DX3,     &
          DELTA,PSI,     &
          DIVAA1,DIVBB2,DIVCC3,     &
          ijkcell,ISTEP,IOTA,     &
          periodic,     &
          toroid,strait,rmaj,RWALL,     &
          BZI,Q0)
!
!     a mesh generator based on
!     Winslow's variable diffusion
!     functional
!
      use vast_kind_param, ONLY : double
      use blcom_com_M, ONLY:      &
        tsix,tsiy,tsiz,     &
        etax,etay,etaz,      &
        nux,nuy,nuz
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use cophys_com_M, ONLY : cdlt,sdlt,dz
      use corgan_com_M, ONLY : itdim
      implicit real*8 (a-h,o-z)
!
      real(double) ::  x(*),     &
             y(*),     &
             z(*),     &
             w(*),     &
          ax(itdim),ay(itdim),az(itdim),     &
          bx(itdim),by(itdim),bz(itdim),     &
          cx(itdim),cy(itdim),cz(itdim),     &
          dx1(*),dx2(*),dx3(*),     &
          divaa1(*),divbb2(*),divcc3(*),     &
          vol(*),     &
          ijkcell(*),ISTEP(*),IOTA(*)
!
      real(double) :: jacob,lama,lamb,lamc, zero
!
      logical :: periodic
!
      istart=3
      jstart=3
      if(periodic) then
      istart=2
      jstart=2
      endif
      zero=0.0d0
!
!     construct an array of cell indices
!
      ncells=0
      do 1 i=istart,nxp-1
      do 1 j=jstart,nyp-1
      do 1 k=3,nzp-1
!
      ncells=ncells+1
!
      ijkcell(ncells)=(i-1)*iwid+(j-1)*jwid+(k-1)*kwid+1
!
      ijk=ijkcell(ncells)
      divaa1(ijk)=0.0
      divbb2(ijk)=0.0
      divcc3(ijk)=0.0
!
    1 continue
!
!     calculate the divergence of the reference vectors,
!     store in divaa1, divbb2, divcc3
!
      if(lama.ne.0.) then
!
      call divaa(iwid,jwid,kwid,ncells,ijkcell,     &
          ax,ay,az,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          tsix,tsiy,tsiz,divaa1)
!
      end if
!
      if(lamb.ne.0.0) then
!
      call divaa(iwid,jwid,kwid,ncells,ijkcell,    &
          bx,by,bz,    &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,    &
          etax,etay,etaz,divbb2)
!
      endif
!
!
      do 1000 iter=1,numit
!
!     compute weight function for toroidal mesh
!
      call torusw(nxp,nyp,nzp,1,nxp,1,nyp,2,nzp,     &
          x,y,z,w)
!
      call adaptw(nxp,nyp,nzp,1,nxp,1,nyp,2,nzp,    &
          toroid,strait,rmaj,    &
          x,y,z,w)
!
      call contras(iwid,jwid,kwid,ncells,ijkcell,     &
          x,y,z,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz)
!
!    monitor minimization
!
      call diagnos(ijkcell,ncells,     &
          ax,ay,az,bx,by,bz,cx,cy,cz,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          divaa1,divbb2,divcc3,diva,divb,divc,     &
          adot,bdot,cdot)
!
!     call routines to cause grid to follow flux surfaces
!
      if(lama+lamb+lamc.gt.0.0) then
!
!
      call fluxsurf(1,nxp,1,nyp,1,nzp,iwid,jwid,kwid,     &
          x,y,z,rwall,rmaj,toroid,strait,delta,     &
          psi)
!
      call bvec(nxp,nyp,nzp,ijkcell,ncells,IWID,JWID,KWID,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          rmaj,RWALL,strait,toroid,     &
          Q0,BZI,PSI,     &
          x,y,z,     &
          bx,by,bz)
!
      call refvec(2,nxp,2,nyp,2,nzp,iwid,jwid,kwid,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          x,y,z,     &
          TOROID,STRAIT,RWALL,RMAJ,     &
          psi,cx,cy,cz,AX,AY,AZ,BX,BY,BZ)
!
      call torusbc(nxp,nyp,nzp,     &
          zero,                     &
          ax,ay,az)
!
      call torusbc(nxp,nyp,nzp,     &
          zero,                     &
          bx,by,bz)
!
      call torusbc(nxp,nyp,nzp,      &
          zero,                      &
          cx,cy,cz)
!
      if(lama.ne.0.0) then
!
      call divaa(iwid,jwid,kwid,ncells,ijkcell,      &
          cx,cy,cz,      &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,      &
          nux,nuy,nuz,divcc3)

      endif

      if(lamb.ne.0.0) then
!
      call divaa(iwid,jwid,kwid,ncells,ijkcell,      &
          ax,ay,az,      &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,      &
          tsix,tsiy,tsiz,divaa1)
!
      endif
!
      if(lamc.ne.0.0) then
!
      call divaa(iwid,jwid,kwid,ncells,ijkcell,     &
          bx,by,bz,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          etax,etay,etaz,divbb2)
!
      endif
!
      endif
!
!     calculate change in contravariant components of mesh vector
!
      call gridsolv(x,y,z,w,     &
          AX,AY,AZ,LAMA,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          tsix,tsiy,tsiz,     &
          DIVAA1,     &
          iwid,jwid,kwid,ijkcell,ncells,     &
          DX1)
!
      call gridsolv(x,y,z,w,     &
          BX,BY,BZ,LAMB,     &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,     &
          etax,etay,etaz,     &
          DIVBB2,     &
          iwid,jwid,kwid,ijkcell,ncells,     &
          DX2)
!
      call gridsolv(x,y,z,w,      &
          CX,CY,CZ,LAMC,      &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,      &
          nux,nuy,nuz,      &
          DIVCC3,      &
          iwid,jwid,kwid,ijkcell,ncells,      &
          DX3)
!
      do 10 n=1,ncells
!
      ijk=ijkcell(n)
!
      x1=0.5*(x(ijk+iwid)-x(ijk-iwid))
      x2=0.5*(x(ijk+jwid)-x(ijk-jwid))
      x3=0.5*(x(ijk+kwid)-x(ijk-kwid))
!
      y1=0.5*(y(ijk+iwid)-y(ijk-iwid))
      y2=0.5*(y(ijk+jwid)-y(ijk-jwid))
      y3=0.5*(y(ijk+kwid)-y(ijk-kwid))
!
      z1=0.5*(z(ijk+iwid)-z(ijk-iwid))
      z2=0.5*(z(ijk+jwid)-z(ijk-jwid))
      z3=0.5*(z(ijk+kwid)-z(ijk-kwid))
!
      dx=dx1(ijk)*x1+(dx2(ijk)*x2+(dx3(ijk)*x3))
      dy=dx1(ijk)*y1+(dx2(ijk)*y2+(dx3(ijk)*y3))
      dz=dx1(ijk)*z1+(dx2(ijk)*z2+(dx3(ijk)*z3))
!
      x(ijk)=x(ijk)+dx
      y(ijk)=y(ijk)+dy
      z(ijk)=z(ijk)+dz
!
   10 continue
!
!     impose periodic boundary conditions for torus
!
      call torusbg(nxp,nyp,nzp,istep,IOTA,     &
          x,y,z)
!
 1000 continue
!
      return
      end subroutine meshgen
!     ************************************************************************
      subroutine torusbc(nxp,nyp,nzp,     &
          strait,                         &
          x,y,z)
!
!     a routine to impose double periodicity
!     for toroidal geometry
!     sets ghost cell values for cell-centered vector
!     no assembly
!
!     called by SGRID3D, GEOM_JUB, DIVB_PROJECTION,
!     ACCEL_3DMHD, DEBUG, BUDGET VINIT
!
      use vast_kind_param, ONLY : double
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z
      use cophys_com_M, ONLY : cdlt,sdlt,dz
      implicit real*8 (a-h,o-z)
!
      real(double) ::  x(nxp,nyp,*),y(nxp,nyp,*),z(nxp,nyp,*)
      real(double), intent(in) :: strait
!
      do 10 k=2,nzp
!
!     periodicity in the toroidal angle
!
      if(periodic_y) then
!
      do 1 i=1,nxp
!
      x(i,1,k)=cdlt*x(i,nyp-1,k)+sdlt*y(i,nyp-1,k)
      y(i,1,k)=-sdlt*x(i,nyp-1,k)+cdlt*y(i,nyp-1,k)     &
          -strait*dz
      z(i,1,k)=z(i,nyp-1,k)
!
      x(i,nyp,k)=cdlt*x(i,2,k)-sdlt*y(i,2,k)
      y(i,nyp,k)=sdlt*x(i,2,k)+cdlt*y(i,2,k)     &
          +strait*dz
      z(i,nyp,k)=z(i,2,k)
!
    1 continue
!
      else
!
      do i=1,nxp
!
      x(i,1,k)=x(i,2,k)
      y(i,1,k)=y(i,2,k)
      z(i,1,k)=z(i,2,k)
!
      x(i,nyp,k)=x(i,nyp-1,k)
      y(i,nyp,k)=y(i,nyp-1,k)
      z(i,nyp,k)=z(i,nyp-1,k)
!
      enddo
      endif

!
!     periodicity in the poloidal angle
!
      if(periodic_x) then
!
      do 2 j=1,nyp
!
      x(1,j,k)=x(nxp-1,j,k)
      y(1,j,k)=y(nxp-1,j,k)
      z(1,j,k)=z(nxp-1,j,k)
!
      x(nxp,j,k)=x(2,j,k)
      y(nxp,j,k)=y(2,j,k)
      z(nxp,j,k)=z(2,j,k)
!
    2 continue
!
      else
!
      do j=1,nyp
!
      x(1,j,k)=x(2,j,k)
      y(1,j,k)=y(2,j,k)
      z(1,j,k)=z(2,j,k)
!
      x(nxp,j,k)=x(nxp-1,j,k)
      y(nxp,j,k)=y(nxp-1,j,k)
      z(nxp,j,k)=z(nxp-1,j,k)
!
      enddo
!
      endif
!
   10 continue
!
      if(periodic_z) then
!
      do i=1,nxp
      do j=1,nyp
!
      x(i,j,1)=x(i,j,nzp-1)
      y(i,j,1)=y(i,j,nzp-1)
      z(i,j,1)=z(i,j,nzp-1)
!
      x(i,j,nzp)=x(i,j,2)
      y(i,j,nzp)=y(i,j,2)
      z(i,j,nzp)=z(i,j,2)
!
      enddo
      enddo
!
      else
!
      do i=1,nxp
      do j=1,nyp
!
      x(i,j,1)=x(i,j,2)
      y(i,j,1)=y(i,j,2)
      z(i,j,1)=z(i,j,2)
!
      x(i,j,nzp)=x(i,j,nzp-1)
      y(i,j,nzp)=y(i,j,nzp-1)
      z(i,j,nzp)=z(i,j,nzp-1)
!
      enddo
      enddo
!
      endif
      return
      end subroutine torusbc
!     ********************************************************************
      subroutine torusbg(nxp,nyp,nzp,istep,IOTA,     &
          x,y,z)
!
!     a routine to impose double periodicity
!     for toroidal geometry
!
      use vast_kind_param, ONLY : double
      use cophys_com_M, ONLY : cdlt,sdlt,dz,rwall,rmaj
      use corgan_com_M, ONLY : strait,toroid
      implicit real*8 (a-h,o-z)
!
      real(double) x(nxp,nyp,*),y(nxp,nyp,*),z(nxp,nyp,*),    &
          IOTA(*)
      integer :: istep(*)
!
!     istep is recalculated as function of minor radius
!
      delt=acos(cdlt)
      pi=acos(-1.)
      ibar=nxp-2
      kbar=nzp-2
      dr=rwall/kbar
!
      do 11 k=2,nzp
      ws=toroid*sqrt(x(2,2,k)**2+y(2,2,k)**2)-rmaj+strait*x(2,2,k)
      wsr=sqrt(ws**2+z(2,2,k)**2)+1.e-20
      kt=int(wsr/dr)+2
      istep(k)=int(iota(kt)*ibar+0.5)*delt/(2.*pi)
  11  continue
!
!
      do 10 k=2,nzp
!
!     periodicity in the toroidal angle
!
      do 1 i=1,nxp
!
      i1=i+istep(k)
      if(i1.lt.1) i1=i1+nxp-2
      if(i1.gt.nxp) i1=i1-nxp+2
!
      x(i1,1,k)=cdlt*x(i,nyp-1,k)+sdlt*y(i,nyp-1,k)
      y(i1,1,k)=-sdlt*x(i,nyp-1,k)+cdlt*y(i,nyp-1,k)      &
          -strait*dz
      z(i1,1,k)=z(i,nyp-1,k)
!
      x(i,nyp,k)=cdlt*x(i1,2,k)-sdlt*y(i1,2,k)
      y(i,nyp,k)=sdlt*x(i1,2,k)+cdlt*y(i1,2,k)       &
          +strait*dz
      z(i,nyp,k)=z(i1,2,k)
!
    1 continue
!
!
!     periodicity in the poloidal angle
!
      do 2 j=1,nyp
!
      x(1,j,k)=x(nxp-1,j,k)
      y(1,j,k)=y(nxp-1,j,k)
      z(1,j,k)=z(nxp-1,j,k)
!
      x(nxp,j,k)=x(2,j,k)
      y(nxp,j,k)=y(2,j,k)
      z(nxp,j,k)=z(2,j,k)
!
    2 continue
!
   10 continue
!
      return
      end subroutine torusbg
!     *****************************************************************
      subroutine torusw(nxp,nyp,nzp,i1,i2,j1,j2,k1,k2,     &
          x,y,z,wgrid)
!
!     a routine to calculate the weight function for a mesh in toroidal
!     geometry
!
!     the weight function is the Jacobian for toroidal coordinates
!
      use vast_kind_param, ONLY : double
      use cophys_com_M, ONLY : rmaj
      use corgan_com_M, ONLY : toroid,strait
      implicit real*8 (a-h,o-z)
!
      real(double) :: x(nxp,nyp,*),y(nxp,nyp,*),z(nxp,nyp,*),     &
               wgrid(nxp,nyp,*)
!
      do 1 k=k1,k2
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      ws=toroid*sqrt(x(i,j,k)**2+y(i,j,k)**2)-rmaj+strait*x(i,j,k)
      wsr=sqrt(ws**2+z(i,j,k)**2)+1.e-20
      wsmr=toroid*sqrt(x(i,j,k)**2+y(i,j,k)**2)+strait
!
!     to prevent collapse near the origin
!
      wsrp=1./(1./wsr+0.25/wsr**2)
!
      wgrid(i,j,k)=wsrp*wsmr
!
    1 continue
!
      return
      end subroutine torusw
!     ****************************************************************************
