*dk gridinit
      subroutine gridinit(ibar,jbar,kbar,iwid,jwid,kwid,
     &     delt,dphi,dtheta,rwall,rmaj,dzstr,istep,
     &     del1,del2,del3,del4,del5,del6,del7,
     &     x,y,z)
c
      implicit real*8 (a-h,o-z)
c
      dimension x(*),y(*),z(*),istep(*)
c
c     impose a twist on the mesh
c     one rotation over the length of the mesh in j
c
      rfibar=1./real(ibar)
      rfjbar=1./real(jbar)
      rfkbar=1./real(kbar)
c
      pi=acos(-1.)
c
      do 90 k=1,kbar+2
c
      dtdj=real(istep(k))*dtheta*rfjbar
      dt0=-0.5*real(istep(k))*dtheta
c
      do 88 j=1,jbar+2
c
      phiang=(j-2)*dphi
c
      do 86 i=1,ibar+2
c
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
c
      theta=(i-2)*dtheta+(j-2)*dtdj+dt0
      theta0=theta-0.5*pi
      delr=rwall*rfkbar*(1.+del1*cos(theta0)
     &                 +del2*cos(2.*(theta0+0.5*pi))
     &                 +del3*cos(3.*theta0)
     &                 +del4*cos(4.*theta0)
     &                 +del5*cos(5.*theta0)
     &                 +del7*cos(7.*theta0))
c
         rr = (k-2)*delr
         ws = rmaj+rr*sin(theta)
         xx =ws*cos(phiang)
         x(ijk) = xx
         yy =ws*sin(phiang)+real(j-2)*dzstr
         y(ijk) = yy
         z(ijk) = rr*cos(theta)
c
   86 continue
c
   88 continue
c
   90 continue
c
c
      return
      end
