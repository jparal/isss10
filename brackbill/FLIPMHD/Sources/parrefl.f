      subroutine parrefl(ncellsp,ijkctmp,iphead,npart,itdim,   &
          iwid,jwid,kwid,x,y,z,wate,bxv,byv,bzv,   &
          px,py,pz,up,vp,wp,pxi,peta,pzta)
!
!     ***************************************************************
!
!     a routine to reflect a particle at a surface
!     px,py,pz are the coordinates of a point outside the surface
!     pxi,peta,pzta are the natural coordinates of the intersection of the particle
!                orbit with the surface
!     the new position is calculated by reflecting in the normal
!     Naitou's reflection conditions are applied to the velocity
!
!     **************************************************************
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      integer :: ijkctmp(*),iphead(*)
      real(double) ::    &
          x(*),y(*),z(*),bxv(*),byv(*),bzv(*),    &
          px(0:npart),py(0:npart),pz(0:npart),    &
          up(0:npart),vp(0:npart),wp(0:npart),    &
          pxi(0:npart),peta(0:npart),pzta(0:npart),    &
          wate(itdim,*)
!
      real(double) :: normx,normy,normz
!
!
!dir$ ivdep
!
      do 1 n=1,ncellsp
!
      ijk_host=ijkctmp(n)
      np=iphead(ijk_host)
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
      ijk_chek=(knew-1)*kwid+(jnew-1)*jwid+(inew-1)*iwid+1
!
      tsi=pxi(np)-inew
      eta=peta(np)-jnew
!
      dx1=(1.-eta)*(x(ijk_chek+iwid+kwid)-x(ijk_chek+kwid))   &
          +eta*(x(ijk_chek+iwid+jwid+kwid)-x(ijk_chek+jwid+kwid))
!
      dx2=(1.-tsi)*(x(ijk_chek+jwid+kwid)-x(ijk_chek+kwid))    &
          +tsi*(x(ijk_chek+iwid+jwid+kwid)-x(ijk_chek+iwid+kwid))
!
      dy1=(1.-eta)*(y(ijk_chek+iwid+kwid)-y(ijk_chek+kwid))     &
          +eta*(y(ijk_chek+iwid+jwid+kwid)-y(ijk_chek+jwid+kwid))
!
      dy2=(1.-tsi)*(y(ijk_chek+jwid+kwid)-y(ijk_chek+kwid))     &
          +tsi*(y(ijk_chek+iwid+jwid+kwid)-y(ijk_chek+iwid+kwid))
!
      dz1=(1.-eta)*(z(ijk_chek+iwid+kwid)-z(ijk_chek+kwid))     &
          +eta*(z(ijk_chek+iwid+jwid+kwid)-z(ijk_chek+jwid+kwid))
!
      dz2=(1.-tsi)*(z(ijk_chek+jwid+kwid)-z(ijk_chek+kwid))     &
          +tsi*(z(ijk_chek+iwid+jwid+kwid)-z(ijk_chek+iwid+kwid))
!
      normx=dy1*dz2-dy2*dz1
      normy=dz1*dx2-dz2*dx1
      normz=dx1*dy2-dx2*dy1
!
      x_intersect=tsi*(x(ijk_chek+iwid+kwid)*(1.-eta)    &
             +x(ijk_chek+iwid+jwid+kwid)*eta)    &
       +(1.-tsi)*(x(ijk_chek+kwid)*(1.-eta)    &
                  +x(ijk_chek+jwid+kwid)*eta)
!
      y_intersect=tsi*(y(ijk_chek+iwid+kwid)*(1.-eta)   &
             +y(ijk_chek+iwid+jwid+kwid)*eta)   &
       +(1.-tsi)*(y(ijk_chek+kwid)*(1.-eta)   &
                  +y(ijk_chek+jwid+kwid)*eta)
!
      z_intersect=tsi*(z(ijk_chek+iwid+kwid)*(1.-eta)    &
             +z(ijk_chek+iwid+jwid+kwid)*eta)    &
       +(1.-tsi)*(z(ijk_chek+kwid)*(1.-eta)    &
                  +z(ijk_chek+jwid+kwid)*eta)
!
!
!
      rnorm=1./(normx**2+normy**2+normz**2)
!
!
      dxdotn=((px(np)-x_intersect)*normx    &
            +(py(np)-y_intersect)*normy    &
            +(pz(np)-z_intersect)*normz)    &
          *rnorm
!
      px(np)=px(np)-2.*normx*dxdotn
      py(np)=py(np)-2.*normy*dxdotn
      pz(np)=pz(np)-2.*normz*dxdotn
!
      udotn=(up(np)*normx+vp(np)*normy+wp(np)*normz)  &
         *rnorm
!
!      up(np)=+up(np)-2.*normx*udotn
!      vp(np)=+vp(np)-2.*normy*udotn
!      wp(np)=+wp(np)-2.*normz*udotn
!
    1 continue
!
      return
      end subroutine parrefl
