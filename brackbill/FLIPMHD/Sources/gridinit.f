      subroutine gridinit(ibar,jbar,kbar,iwid,jwid,kwid,   &
          delt,dphi,dtheta,dx,dy,dz,rwall,rmaj,dzstr,istep,   &
          del1,del2,del3,del4,del5,del6,del7,   &
          x,y,z,cartesian,xl,xr,yb,yt,ze,zf)
!
      use vast_kind_param, ONLY : double
      implicit none
      real(double) :: x(*),y(*),z(*)
      real(double) :: delt,dphi,dtheta,dx,dy,dz,rwall,rmaj,dzstr,   &
       del1,del2,del3,del4,del5,del6,del7,   &
       xl,xr,yb,yt,ze,zf
      integer :: ibar,jbar,kbar,iwid,jwid,kwid,istep(*)
      integer :: i,j,k,ijk
      real(double) :: delr,theta0,theta,rr,xx,yy,ws,dtdj,pi,phiang,dt0
      logical cartesian
!
      real(double) :: rfibar,rfjbar,rfkbar
!
!     impose a twist on the mesh
!     one rotation over the length of the mesh in j
!
      if (cartesian) then
	call gridinit_cart(ibar,jbar,kbar,iwid,jwid,kwid,   &
        dx,dy,dz,x,y,z,xl,xr,yb,yt,ze,zf)
	return
      endif
      rfibar=1./float(ibar)
      rfjbar=1./float(jbar)
      rfkbar=1./float(kbar)
!
      pi=acos(-1.)
!
      do 90 k=1,kbar+2
!
      dtdj=float(istep(k))*dtheta*rfjbar
      dt0=-0.5*float(istep(k))*dtheta
!
      do 88 j=1,jbar+2
!
      phiang=(j-2)*dphi
!
      do 86 i=1,ibar+2
!
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
!
      theta=(i-2)*dtheta+(j-2)*dtdj+dt0
      theta0=theta-0.5*pi
      delr=rwall*rfkbar*(1.+del1*cos(theta0)     &
                      +del2*cos(2.*(theta0+0.5*pi))     &
                      +del3*cos(3.*theta0)     &
                      +del4*cos(4.*theta0)     &
                      +del5*cos(5.*theta0)     &
                      +del7*cos(7.*theta0))
!
         rr = (k-2)*delr
         ws = rmaj+rr*sin(theta)
         xx =ws*cos(phiang)
         x(ijk) = xx
         yy =ws*sin(phiang)+float(j-2)*dzstr
         y(ijk) = yy
         z(ijk) = rr*cos(theta)
!
   86 continue
!
   88 continue
!
   90 continue
!
!
      return
      end subroutine gridinit

      subroutine gridinit_cart(ibar,jbar,kbar,iwid,jwid,kwid,     &
        dx,dy,dz,x,y,z,xl,xr,yb,yt,ze,zf)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),z(*)

      do k=1,kbar+2
      do j=1,jbar+2
      do i=1,ibar+2
      
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
    
      x(ijk)=(i-2)*dx
      y(ijk)=(j-2)*dy
      z(ijk)=(k-2)*dz

      enddo
      enddo
      enddo

      xl=0.
      xr=float(ibar)*dx
      yb=0.
      yt=float(jbar)*dy
      ze=0.
      zf=float(kbar)*dz
      return
      end subroutine gridinit_cart
