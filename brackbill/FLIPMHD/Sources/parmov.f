      subroutine parmov(epsuv)
!
      use vast_kind_param, ONLY : double
      use corgan_com_M, ONLY :  itdim, npart, mgeom
      use cindex_com_M, ONLY : ncells,iwid,jwid,kwid,  &
           ibp1,jbp1,kbp1
      use cophys_com_M, ONLY : dx, dy, dz,   &
          rmaj,cdlt,sdlt
      use blcom_com_M, ONLY : dudx,dudy,dudz,   &
           dvdx,dvdy,dvdz,    &
           dwdx,dwdy,dwdz,    &
           area_x,area_y,area_z,    &
           cartesian,divu,         &
           CurlEx,CurlEy,CurlEz,number,vol,    &
           iphead,link,iphd2,pxi,peta,pzta,  &
           ijkParticle,               &
           wate,ul,vl,wl,ax,ay,az,    &
           px,py,pz,mass,     &
           x,y,z,             &
           tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,   &
           mupx,mupy,mupz,    &
           dbvx,dbvy,dbvz,    &
           ep,work,pdv,       &
           ijkcell,           &
           xr,xl,yb,yt,ze,zf, &
           up,vp,wp,          &
           periodic_x,periodic_y,periodic_z
      use numpar_com_M, ONLY : dt
      use Scratch_com_M, ONLY : upv1, vpv1, wpv1,   &
         dudxp, dudyp, dudzp,     &
         dvdxp, dvdyp, dvdzp,     &
         dwdxp, dwdyp, dwdzp,     &
         axv1, ayv1, azv1,        &
         the, zeta, nu ,          &
         workp
      use Timing
      implicit none
!
      integer :: OMP_GET_THREAD_NUM
      integer ::    &
          ijknew,index,ll,np,ijk,n,    &
          ijkstep,   &
          inew,jnew,knew,l,    &
          i, j, k, parcount
!
!
      real(double) :: epsuv,    &
         xcount, ycount, zcount,    &
         b22,b23,b32,b33,rhs2,rhs3,  &
         wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp
      real(double) :: parmass, parnrg, momx, momy, momz,   &
         mux, muy,muz
      real(double) :: Tstart, Tfinish
!
      dimension ijkstep(27)
!
      allocate (upv1(itdim), vpv1(itdim), wpv1(itdim))
      allocate (axv1(itdim), ayv1(itdim), azv1(itdim))
      allocate (workp(itdim))
      allocate (the(itdim), zeta(itdim), nu(itdim))
      allocate (dudxp(itdim), dudyp(itdim), dudzp(itdim))
      allocate (dvdxp(itdim), dvdyp(itdim), dvdzp(itdim))
      allocate (dwdxp(itdim), dwdyp(itdim), dwdzp(itdim))
!
!     *******************************************
!
!     an explicit particle mover 
!
!     *********************************************
!
      call cpu_time(Tstart)
!
!
      call celstep(iwid,jwid,kwid,ijkstep)
!
      do 1111 n=1,ncells
         ijk=ijkcell(n)
         if(iphead(ijk).gt.0) then
!
!     normalize the contribution from resistive diffusion
!
            CurlEx(ijk)=CurlEx(ijk)*vol(ijk)*dt/number(ijk)
            CurlEy(ijk)=CurlEy(ijk)*vol(ijk)*dt/number(ijk)
            CurlEz(ijk)=CurlEz(ijk)*vol(ijk)*dt/number(ijk)
         else
            CurlEx(ijk)=0.0d0
            CurlEy(ijk)=0.0d0
            CurlEz(ijk)=0.0d0
         endif
 1111 continue
!
        call OMP_SET_NUM_THREADS(8)
!$omp parallel do private(n,ijk,np,ll,i,j,k,   &
!$omp inew,jnew,knew,                 &
!$omp wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp)       &
!$omp shared(dt,xl,xr,yb,yt,ze,zf,             &
!$omp periodic_x,periodic_y,periodic_z),       &
!$omp reduction(+:xcount,ycount,zcount)
!
     do n=1,ncells
        ijk=ijkcell(n)
        np=iphead(ijk)
        do while (np.gt.0)
!
!
!       calculate weights for trilinear interpolation
!
           i=int(pxi(np))
           j=int(peta(np))
           k=int(pzta(np))
!
           the(n)=pxi(np)-real(i)
           zeta(n)=peta(np)-real(j)
           nu(n)=pzta(np)-real(k)

!
!       compute interpolation weights
!
           wate(n,1) = the(n)*(1.-zeta(n))*(1.-nu(n))
           wate(n,2) = the(n)*zeta(n)*(1.-nu(n))
           wate(n,3) = (1.-the(n))*zeta(n)*(1.-nu(n))
           wate(n,4) = (1.-the(n))*(1.-zeta(n))*(1.-nu(n))
!
           wate(n,5) = the(n)*(1.-zeta(n))*nu(n)
           wate(n,6) = the(n)*zeta(n)*nu(n)
           wate(n,7) = (1.-the(n))*zeta(n)*nu(n)
           wate(n,8) = (1.-the(n))*(1.-zeta(n))*nu(n)

!
!    interpolate grid velocity to particle
!
             upv1(n)=   wate(n,1)*ul(ijk+iwid)                &
                      + wate(n,2)*ul(ijk+iwid+jwid)          &
                      + wate(n,3)*ul(ijk+jwid)               &
                      + wate(n,4)*ul(ijk)                    &
                      + wate(n,5)*ul(ijk+iwid+kwid)          &
                      + wate(n,6)*ul(ijk+iwid+jwid+kwid)     &
                      + wate(n,7)*ul(ijk+jwid+kwid)          &
                      + wate(n,8)*ul(ijk+kwid)

             vpv1(n)=   wate(n,1)*vl(ijk+iwid)                &
                      + wate(n,2)*vl(ijk+iwid+jwid)          &
                      + wate(n,3)*vl(ijk+jwid)               &
                      + wate(n,4)*vl(ijk)                    &
                      + wate(n,5)*vl(ijk+iwid+kwid)          &
                      + wate(n,6)*vl(ijk+iwid+jwid+kwid)     &
                      + wate(n,7)*vl(ijk+jwid+kwid)          &
                      + wate(n,8)*vl(ijk+kwid)

              wpv1(n)=   wate(n,1)*wl(ijk+iwid)                &
                      + wate(n,2)*wl(ijk+iwid+jwid)          &
                      + wate(n,3)*wl(ijk+jwid)               &
                      + wate(n,4)*wl(ijk)                    &
                      + wate(n,5)*wl(ijk+iwid+kwid)          &
                      + wate(n,6)*wl(ijk+iwid+jwid+kwid)     &
                      + wate(n,7)*wl(ijk+jwid+kwid)          &
                      + wate(n,8)*wl(ijk+kwid)

!     interpolate grid acceleration to particles
!
               axv1(n)=   wate(n,1)*ax(ijk+iwid)                &
                        + wate(n,2)*ax(ijk+iwid+jwid)          &
                        + wate(n,3)*ax(ijk+jwid)               &
                        + wate(n,4)*ax(ijk)                    &
                        + wate(n,5)*ax(ijk+iwid+kwid)          &
                        + wate(n,6)*ax(ijk+iwid+jwid+kwid)     &
                        + wate(n,7)*ax(ijk+jwid+kwid)          &
                        + wate(n,8)*ax(ijk+kwid)

               ayv1(n)=   wate(n,1)*ay(ijk+iwid)                &
                        + wate(n,2)*ay(ijk+iwid+jwid)          &
                        + wate(n,3)*ay(ijk+jwid)               &
                        + wate(n,4)*ay(ijk)                    &
                        + wate(n,5)*ay(ijk+iwid+kwid)          &
                        + wate(n,6)*ay(ijk+iwid+jwid+kwid)     &
                        + wate(n,7)*ay(ijk+jwid+kwid)          &
                        + wate(n,8)*ay(ijk+kwid)

               azv1(n)=   wate(n,1)*az(ijk+iwid)                &
                        + wate(n,2)*az(ijk+iwid+jwid)          &
                        + wate(n,3)*az(ijk+jwid)               &
                        + wate(n,4)*az(ijk)                    &
                        + wate(n,5)*az(ijk+iwid+kwid)          &
                        + wate(n,6)*az(ijk+iwid+jwid+kwid)     &
                        + wate(n,7)*az(ijk+jwid+kwid)          &
                        + wate(n,8)*az(ijk+kwid)

               workp(n)=   wate(n,1)*work(ijk+iwid)                &
                         + wate(n,2)*work(ijk+iwid+jwid)          &
                         + wate(n,3)*work(ijk+jwid)               &
                         + wate(n,4)*work(ijk)                    &
                         + wate(n,5)*work(ijk+iwid+kwid)          &
                         + wate(n,6)*work(ijk+iwid+jwid+kwid)     &
                         + wate(n,7)*work(ijk+jwid+kwid)          &
                         + wate(n,8)*work(ijk+kwid)
!
!     advance particle positions and velocity
!
!
               px(np)=px(np)+upv1(n)*dt
               py(np)=py(np)+vpv1(n)*dt
               pz(np)=pz(np)+wpv1(n)*dt
!
               up(np)=up(np)+axv1(n)
               vp(np)=vp(np)+ayv1(n)
               wp(np)=wp(np)+azv1(n)
!
!
!     adjust particle internal energy
!
              ep(np)=ep(np)+mass(np)*(workp(n)    &
                  -0.5*(axv1(n)**2+ayv1(n)**2+azv1(n)**2))
!
!     update magnetic moment
!     include contribution of resistive diffusion
!
              mupx(np)=mupx(np)+mass(np)*dbvx(ijk)
              mupy(np)=mupy(np)+mass(np)*dbvy(ijk)
              mupz(np)=mupz(np)+mass(np)*dbvz(ijk)
!
!
     the(n)=the(n)-0.5
     zeta(n)=zeta(n)-0.5
     nu(n)=nu(n)-0.5
!             call OneWateC(n,the,zeta,nu,wate)
!
!
      wi=0.75-the(n)**2
      wim=0.5*(0.5-the(n))**2
      wip=0.5*(0.5+the(n))**2
!
      wj=0.75-zeta(n)**2
      wjm=0.5*(0.5-zeta(n))**2
      wjp=0.5*(0.5+zeta(n))**2
!
      wk=0.75-nu(n)**2
      wkm=0.5*(0.5-nu(n))**2
      wkp=0.5*(0.5+nu(n))**2
!
!     k-plane
!
      wate(n,1)=wi*wj*wk
      wate(n,2)=wip*wj*wk
      wate(n,3)=wip*wjp*wk
      wate(n,4)=wi*wjp*wk
      wate(n,5)=wim*wjp*wk
      wate(n,6)=wim*wj*wk
      wate(n,7)=wim*wjm*wk
      wate(n,8)=wi*wjm*wk
      wate(n,9)=wip*wjm*wk
!
!     k-1 - plane
!
      wate(n,10)=wi*wj*wkm
      wate(n,11)=wip*wj*wkm
      wate(n,12)=wip*wjp*wkm
      wate(n,13)=wi*wjp*wkm
      wate(n,14)=wim*wjp*wkm
      wate(n,15)=wim*wj*wkm
      wate(n,16)=wim*wjm*wkm
      wate(n,17)=wi*wjm*wkm
      wate(n,18)=wip*wjm*wkm
!
!     k+1 - plane
!
      wate(n,19)=wi*wj*wkp
      wate(n,20)=wip*wj*wkp
      wate(n,21)=wip*wjp*wkp
      wate(n,22)=wi*wjp*wkp
      wate(n,23)=wim*wjp*wkp
      wate(n,24)=wim*wj*wkp
      wate(n,25)=wim*wjm*wkp
      wate(n,26)=wi*wjm*wkp
      wate(n,27)=wip*wjm*wkp
!

!
             do 350 ll=1,27
                ep(np)=ep(np)+mass(np)*wate(n,ll)*pdv(ijk+ijkstep(ll))
  350        continue
!
!     calculate new natural coordinates for the particle
!
   21 if(px(np).lt.xl.or.px(np).gt.xr) then
       xcount=xcount+1
!
        if(px(np).lt.xl) then
         if(periodic_x) then
           px(np)=px(np)+xr-xl
         else
           px(np)=xl+2.*(xl-px(np))
         endif
        endif
!
        if(px(np).gt.xr) then
         if(periodic_x) then
           px(np)=px(np)-(xr-xl)
         else
           px(np)=xr-2.*(px(np)-xr)
         endif
        endif
!
        go to 21
        endif
   22 if(py(np).lt.yb.or.py(np).gt.yt) then
       ycount=ycount+1
!
        if(py(np).lt.yb) then
         if(periodic_y) then
           py(np)=py(np)+(yt-yb)
         else
           py(np)=py(np)+2.*(yb-py(np))
         endif
        endif
!
        if(py(np).gt.yt) then
         if(periodic_y) then
           py(np)=py(np)-(yt-yb)
         else
           py(np)=yt-2.*(py(np)-yt)
         endif
        endif
        go to 22
        endif
!

   23  if(pz(np).lt.ze.or.pz(np).gt.zf) then
        zcount=zcount+1
        if(pz(np).lt.ze) then
         if(periodic_z) then
           pz(np)=pz(np)+(zf-ze)
         else
           pz(np)=ze+2.*(ze-pz(np))
         endif
        endif
!
        if(pz(np).gt.zf) then
         if(periodic_z) then
           pz(np)=pz(np)-(zf-ze)
         else
           pz(np)=zf-2.*(pz(np)-zf)
         endif
        endif
!
!
        go to 23
        endif
!
!       calculate the new logical coordinates of the particle
!
         pxi(np)=2.+(px(np)-xl)/dx
         peta(np)=2.+(py(np)-yb)/dy
         pzta(np)=2.+(pz(np)-ze)/dz
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
      if((inew-2)*(ibp1-inew).lt.0)  write(6,*) 'np=',np,'inew=',inew
      if((jnew-2)*(jbp1-jnew).lt.0)  write(6,*) 'np=',np,'jnew=',jnew
      if((knew-2)*(kbp1-knew).lt.0)  write(6,*) 'np=',np,'knew=',knew
!
      inew=min(ibp1,max(2,inew))
      jnew=min(jbp1,max(2,jnew))
      knew=min(kbp1,max(2,knew))
!
      ijkParticle(np)=(knew-1)*kwid+(jnew-1)*jwid+(inew-1)*iwid+1
!
      np=link(np)
!
      enddo
!
   enddo
!
        call OMP_SET_NUM_THREADS(1)

      parcount=0

      iphd2=0
      do n=1,ncells
         ijk=ijkcell(n)
         np=iphead(ijk) 
         do while (np.gt.0)
            parcount=parcount+1
            iphead(ijk)=link(np)
            link(np)=iphd2(ijkParticle(np))
            iphd2(ijkParticle(np))=np
            np=iphead(ijk)
         enddo
       enddo
!
      do n=1,ncells
         ijk=ijkcell(n)
         iphead(ijk)=iphd2(ijk)
         iphd2(ijk)=0
      enddo
!
      deallocate (upv1, vpv1, wpv1)
      deallocate (axv1, ayv1, azv1)
      deallocate (workp)
      deallocate(the,zeta,nu)
      deallocate (dudxp,dudyp,dudzp,dvdxp,dvdyp,dvdzp,dwdxp,dwdyp,dwdzp)
!
      call cpu_time(Tfinish)
      do l=1,20
         if(RoutineName(l).eq.'parmov') CpuTime(l)=CpuTime(l)+Tfinish-Tstart
      enddo
!
      return
      end subroutine parmov
