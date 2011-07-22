      subroutine map3d(ijk,iwid,jwid,kwid,eps,   &
          x,y,z,   &
          xp,yp,zp,tsi,eta,nu)
!
!     a routine to calculate the natural coordinates
!     of a particle
!
      use vast_kind_param, ONLY : double
      implicit none
      real(double) :: x(*),y(*),z(*)
!
      real(double) :: jacob
      integer :: iwid, jwid, kwid
      real(double) :: tsi, eta, nu
      real(double) :: tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz
      real(double) :: dxdxi,dxdeta,dxdnu,dydxi,dydeta,dydnu,dzdxi,dzdeta,dzdnu
      real(double) :: dtsi,deta,dnu
      real(double) :: xp,yp,zp
      real(double) :: eps,rjacob, test, xl,yl,zl
      integer :: itmax, iter, ijk
!
!     initial guess:  particle at the center of the cell
!
      tsi=0.5
      eta=0.5
      nu=0.5
      itmax=10
!
!
      iter=0
    1 continue
!
      iter=iter+1
!
!     compute forward transformation
!
      xl=(1.-nu)*(tsi*(1.-eta)*x(ijk+iwid)    &
                +tsi*eta*x(ijk+iwid+jwid)    &
                +(1.-tsi)*eta*x(ijk+jwid)    &
                +(1.-tsi)*(1.-eta)*x(ijk))    &
          +nu*(tsi*(1.-eta)*x(ijk+iwid+kwid)    &
              +tsi*eta*x(ijk+iwid+jwid+kwid)    &
              +(1.-tsi)*eta*x(ijk+jwid+kwid)    &
              +(1.-tsi)*(1.-eta)*x(ijk+kwid))
!
      yl=(1.-nu)*(tsi*(1.-eta)*y(ijk+iwid)    &
                +tsi*eta*y(ijk+iwid+jwid)    &
                +(1.-tsi)*eta*y(ijk+jwid)    &
                +(1.-tsi)*(1.-eta)*y(ijk))    &
          +nu*(tsi*(1.-eta)*y(ijk+iwid+kwid)    &
              +tsi*eta*y(ijk+iwid+jwid+kwid)    &
              +(1.-tsi)*eta*y(ijk+jwid+kwid)    &
              +(1.-tsi)*(1.-eta)*y(ijk+kwid))
!
      zl=(1.-nu)*(tsi*(1.-eta)*z(ijk+iwid)    &
                +tsi*eta*z(ijk+iwid+jwid)    &
                +(1.-tsi)*eta*z(ijk+jwid)    &
                +(1.-tsi)*(1.-eta)*z(ijk))    &
          +nu*(tsi*(1.-eta)*z(ijk+iwid+kwid)    &
              +tsi*eta*z(ijk+iwid+jwid+kwid)    &
              +(1.-tsi)*eta*z(ijk+jwid+kwid)    &
              +(1.-tsi)*(1.-eta)*z(ijk+kwid))
!
!     compute covariant base vectors
!
      dxdxi=(1.-nu)*((1.-eta)*(x(ijk+iwid)-x(ijk))    &
                  +eta*(x(ijk+iwid+jwid)-x(ijk+jwid)))    &
             +nu*((1.-eta)*(x(ijk+iwid+kwid)-x(ijk+kwid))    &
                     +eta*(x(ijk+iwid+jwid+kwid)-x(ijk+jwid+kwid)))
!
      dxdeta=(1.-nu)*(tsi*(x(ijk+iwid+jwid)-x(ijk+iwid))    &
                  +(1.-tsi)*(x(ijk+jwid)-x(ijk)))    &
               +nu*(tsi*(x(ijk+iwid+jwid+kwid)-x(ijk+iwid+kwid))    &
                +(1.-tsi)*(x(ijk+jwid+kwid)-x(ijk+kwid)))
!
      dxdnu=(1.-tsi)*(eta*(x(ijk+jwid+kwid)-x(ijk+jwid))    &
                 +(1.-eta)*(x(ijk+kwid)-x(ijk)))    &
              +tsi*(eta*(x(ijk+iwid+jwid+kwid)-x(ijk+iwid+jwid))    &
                 +(1.-eta)*(x(ijk+iwid+kwid)-x(ijk+iwid)))
!
      dydxi=(1.-nu)*((1.-eta)*(y(ijk+iwid)-y(ijk))    &
                  +eta*(y(ijk+iwid+jwid)-y(ijk+jwid)))    &
             +nu*((1.-eta)*(y(ijk+iwid+kwid)-y(ijk+kwid))    &
                     +eta*(y(ijk+iwid+jwid+kwid)-y(ijk+jwid+kwid)))
!
      dydeta=(1.-nu)*(tsi*(y(ijk+iwid+jwid)-y(ijk+iwid))    &
                  +(1.-tsi)*(y(ijk+jwid)-y(ijk)))    &
               +nu*(tsi*(y(ijk+iwid+jwid+kwid)-y(ijk+iwid+kwid))    &
                +(1.-tsi)*(y(ijk+jwid+kwid)-y(ijk+kwid)))
!
      dydnu=(1.-tsi)*(eta*(y(ijk+jwid+kwid)-y(ijk+jwid))    &
                 +(1.-eta)*(y(ijk+kwid)-y(ijk)))    &
              +tsi*(eta*(y(ijk+iwid+jwid+kwid)-y(ijk+iwid+jwid))    &
                 +(1.-eta)*(y(ijk+iwid+kwid)-y(ijk+iwid)))
!
      dzdxi=(1.-nu)*((1.-eta)*(z(ijk+iwid)-z(ijk))    &
                  +eta*(z(ijk+iwid+jwid)-z(ijk+jwid)))    &
             +nu*((1.-eta)*(z(ijk+iwid+kwid)-z(ijk+kwid))    &
                     +eta*(z(ijk+iwid+jwid+kwid)-z(ijk+jwid+kwid)))
!
      dzdeta=(1.-nu)*(tsi*(z(ijk+iwid+jwid)-z(ijk+iwid))    &
                  +(1.-tsi)*(z(ijk+jwid)-z(ijk)))    &
               +nu*(tsi*(z(ijk+iwid+jwid+kwid)-z(ijk+iwid+kwid))    &
                +(1.-tsi)*(z(ijk+jwid+kwid)-z(ijk+kwid)))
!
      dzdnu=(1.-tsi)*(eta*(z(ijk+jwid+kwid)-z(ijk+jwid))    &
                 +(1.-eta)*(z(ijk+kwid)-z(ijk)))    &
              +tsi*(eta*(z(ijk+iwid+jwid+kwid)-z(ijk+iwid+jwid))    &
                 +(1.-eta)*(z(ijk+iwid+kwid)-z(ijk+iwid)))
!
!
!
!     calculate the determinant of the metric tensor
!
      jacob=dxdxi*(dydeta*dzdnu-dydnu*dzdeta)    &
          +dydxi*(dzdeta*dxdnu-dzdnu*dxdeta)    &
          +dzdxi*(dxdeta*dydnu-dxdnu*dydeta)
!
      rjacob=1./jacob
!
!     calculate the contravariant base vectors
!
      tsix=(dydeta*dzdnu-dydnu*dzdeta)*rjacob
      tsiy=(dzdeta*dxdnu-dzdnu*dxdeta)*rjacob
      tsiz=(dxdeta*dydnu-dxdnu*dydeta)*rjacob
!
      etax=(dydnu*dzdxi-dydxi*dzdnu)*rjacob
      etay=(dzdnu*dxdxi-dzdxi*dxdnu)*rjacob
      etaz=(dxdnu*dydxi-dxdxi*dydnu)*rjacob
!
      nux=(dydxi*dzdeta-dydeta*dzdxi)*rjacob
      nuy=(dzdxi*dxdeta-dzdeta*dxdxi)*rjacob
      nuz=(dxdxi*dydeta-dxdeta*dydxi)*rjacob
!
!     calculate new value for the natural coordinates
!
      dtsi=(xp-xl)*tsix    &
         +(yp-yl)*tsiy    &
         +(zp-zl)*tsiz
!
      deta=(xp-xl)*etax    &
         +(yp-yl)*etay    &
         +(zp-zl)*etaz
!
      dnu=(xp-xl)*nux    &
         +(yp-yl)*nuy    &
         +(zp-zl)*nuz
!
!    test for convergence
!
      test=dtsi**2+deta**2+dnu**2
!
      if(test.gt.eps**2.and.iter.le.itmax)then
      tsi=tsi+dtsi
      eta=eta+deta
      nu=nu+dnu
!
      go to 1
!
      endif
!
      if(iter.gt.itmax) then
      write(*,*) 'map3d didnt converge..test=',test
      stop 'map3d'
      endif
!
!
!
      return
      end subroutine map3d 
