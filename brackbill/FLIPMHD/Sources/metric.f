      subroutine metric(ncells,ijkcell,iwid,jwid,kwid,    &
          x,y,z,    &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,vol)
!
!     a routine to calculate the contravariant
!     base vectors for a 3d Cartesian mesh
!
      use vast_kind_param, ONLY : double
      implicit none
!
      real(double) :: x(1),y(1),z(1),    &
          tsix(1),tsiy(1),tsiz(1),    &
          etax(1),etay(1),etaz(1)
      real(double) ::  vol(1)
      real(double) :: nux(1),nuy(1),nuz(1)
      real(double) :: dxdnu,dxdxi,dxdeta,   &
                      dydnu,dydxi,dydeta,  &
                      dzdnu,dzdxi,dzdeta
!
      real(double) :: jacob,rjacob
      integer ::  ncells,n,ijkcell(*),iwid,jwid,kwid,  &
         ijk,ipjk,ipjpk,ijpk,ipjkp,ipjpkp,ijpkp,ijkp
!
!
      do 1 n=1,ncells
!
!     compute covariant base vectors
!
      ijk=ijkcell(n)
      ipjk=ijk+iwid
      ipjpk=ijk+iwid+jwid
      ijpk=ijk+jwid
!
      ipjkp=ijk+iwid+kwid
      ipjpkp=ijk+iwid+jwid+kwid
      ijpkp=ijk+jwid+kwid
      ijkp=ijk+kwid
!
      dxdxi=0.25*(x(ipjk)+x(ipjpk)+x(ipjkp)+x(ipjpkp)   &
     &        -x(ijk)-x(ijpk)-x(ijkp)-x(ijpkp))
      dxdeta=0.25*(x(ijpk)+x(ipjpk)+x(ijpkp)+x(ipjpkp)   &
     &        -x(ijk)-x(ipjk)-x(ijkp)-x(ipjkp))
      dxdnu=0.25*(x(ijkp)+x(ipjkp)+x(ipjpkp)+x(ijpkp)   &
     &        -x(ijk)-x(ipjk)-x(ipjpk)-x(ijpk))
!
      dydxi=0.25*(y(ipjk)+y(ipjpk)+y(ipjkp)+y(ipjpkp)   &
     &        -y(ijk)-y(ijpk)-y(ijkp)-y(ijpkp))
      dydeta=0.25*(y(ijpk)+y(ipjpk)+y(ijpkp)+y(ipjpkp)   &
     &        -y(ijk)-y(ipjk)-y(ijkp)-y(ipjkp))
      dydnu=0.25*(y(ijkp)+y(ipjkp)+y(ipjpkp)+y(ijpkp)   &
     &        -y(ijk)-y(ipjk)-y(ipjpk)-y(ijpk))
!
      dzdxi=0.25*(z(ipjk)+z(ipjpk)+z(ipjkp)+z(ipjpkp)   &
     &        -z(ijk)-z(ijpk)-z(ijkp)-z(ijpkp))
      dzdeta=0.25*(z(ijpk)+z(ipjpk)+z(ijpkp)+z(ipjpkp)   &
     &        -z(ijk)-z(ipjk)-z(ijkp)-z(ipjkp))
      dzdnu=0.25*(z(ijkp)+z(ipjkp)+z(ipjpkp)+z(ijpkp)   &
     &        -z(ijk)-z(ipjk)-z(ipjpk)-z(ijpk))
!
!
!     calculate the determinant of the metric tensor
!
      jacob=dxdxi*(dydeta*dzdnu-dydnu*dzdeta)   &
     &     +dydxi*(dzdeta*dxdnu-dzdnu*dxdeta)   &
     &     +dzdxi*(dxdeta*dydnu-dxdnu*dydeta)
!
      vol(ijk)=jacob
      rjacob=1./jacob
!
!     calculate the contravariant base vectors
!
      tsix(ijk)=(dydeta*dzdnu-dydnu*dzdeta)*rjacob
      tsiy(ijk)=(dzdeta*dxdnu-dzdnu*dxdeta)*rjacob
      tsiz(ijk)=(dxdeta*dydnu-dxdnu*dydeta)*rjacob
!
      etax(ijk)=(dydnu*dzdxi-dydxi*dzdnu)*rjacob
      etay(ijk)=(dzdnu*dxdxi-dzdxi*dxdnu)*rjacob
      etaz(ijk)=(dxdnu*dydxi-dxdxi*dydnu)*rjacob
!
      nux(ijk)=(dydxi*dzdeta-dydeta*dzdxi)*rjacob
      nuy(ijk)=(dzdxi*dxdeta-dzdeta*dxdxi)*rjacob
      nuz(ijk)=(dxdxi*dydeta-dxdeta*dydxi)*rjacob
!
!     calculate elements of the metric tensor
!
!
    1 continue
!
      return
      end subroutine metric
