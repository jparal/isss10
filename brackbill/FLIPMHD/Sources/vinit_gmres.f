      subroutine vinit
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use logical_com_M
      use blcom_com_M
      use geometry_com_M
      use cplot_com_M, ONLY : iplot

      use cindex_com_M
      use numpar_com_M
      use cophys_com_M
      use gmres_com_M, ONLY : srce, phi
      use DXFiles
!
      real(double) ::  rnorm,bnorm,refnorm,vvolaxis,    &
         fmo,zero,factor,dummy,totalvolume,vtxmass,cellmass
      real(double) :: Half, One, Bsq, GradPhisq, BGradPhi
      integer :: n,itest,istart,istop,jstart,jstop,kstart,kstop,    &
          itsub,ijkr,ijkl,nvtxtmp
!     a routine to prepare varibles for the next computation cycle
!     ****************************************************
!
      dummy=0.0
      Half=0.5d0
      One=1.d0
!
       if(.not.cartesian) then
!
!      this stuff is a carryover from the toroidal version
!
          if(rmaj.gt.0.0) then
             call torusbcv(ibp1+1,jbp1+1,kbp1+1,     &
                cdlt,sdlt,DUMMY,dz,     &
                periodic_x,periodic_y,periodic_z,     &
             umom,vmom,wmom)
          else
             call torusbc_scalar(ibp2,jbp2,kbp2,     &
                umom)
             call torusbc_scalar(ibp2,jbp2,kbp2,     &
                vmom)
             call torusbc_scalar(ibp2,jbp2,kbp2,     &
                wmom)
          end if
!
   vtxmass=0.0
      do n=1,nvtx
         ijk=ijkvtx(n)
         mv(ijk)=0.0
         do is=1,nsp
            mv(ijk)=mv(ijk)+mv_s(ijk,is)
         enddo
       vtxmass=vtxmass+mv(ijk)
      enddo

      do is=1,nsp
        call torusbcS_scalar(ibp1+1,jbp1+1,kbp1+1,    &
           mv_s,is)
      enddo
!
      call torusbc_scalar(ibp1+1,jbp1+1,kbp1+1,    &
          mv)
!
      call torusbc_scalar(ibp1+1,jbp1+1,kbp1+1,    &
          numberv)
!
      call torusbc_scalar(ibp1+1,jbp1+1,kbp1+1,    &
          color)
!
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,     &
          umom)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,     &
          vmom)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,     &
          wmom)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,     &
          mv)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,     &
          numberv)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,     &
          color)
!
!     *********************************************
!
      else
!
!      this routine imposes periodicity in x and y
!
      do n=1,nvtx
         ijk=ijkvtx(n)
         mv(ijk)=0.0
         do is=1,nsp
            mv(ijk)=mv(ijk)+mv_s(ijk,is)
         enddo
      enddo
!
       call torusbc_scalar(ibp2,jbp2,kbp2,     &
          umom)
!
       call torusbc_scalar(ibp2,jbp2,kbp2,     &
          vmom)
!
       call torusbc_scalar(ibp2,jbp2,kbp2,     &
          wmom)
!
       call torusbc_scalar(ibp2,jbp2,kbp2,     &
          mv)
!
       do is=1,nsp
          call torusbcS_scalar(ibp2,jbp2,kbp2,     &
             mv_s,is)
       enddo
       call torusbc_scalar(ibp2,jbp2,kbp2,     &
          numberv)
!
       call torusbc_scalar(ibp2,jbp2,kbp2,     &
          color)
!
      endif
!
      do 6 n=1,nvtx
      ijk=ijkvtx(n)
      factor=mv(ijk)/(mv(ijk)+1.e-10)**2
!
      u(ijk)=umom(ijk)*factor
      v(ijk)=vmom(ijk)*factor
      w(ijk)=wmom(ijk)*factor
!
      ul(ijk)=u(ijk)
      vl(ijk)=v(ijk)
      wl(ijk)=w(ijk)
!
      factor=numberv(ijk)/(numberv(ijk)+1.e-10)**2
!
      color(ijk)=color(ijk)*factor
!
    6 continue
!
!  apply rigid, free-slip wall conditions
!
      call bc_wall(ibp1,jbp1,kbp1,iwid,jwid,kwid,    &
                  c5x,c6x,c7x,c8x,    &
                  c5y,c6y,c7y,c8y,    &
                  c5z,c6z,c7z,c8z,    &
                  ul,vl,wl)
!
      if(.not.periodic_y) then
!
!     impose no slip conditions on j=2,j=jbp2 boundaries
!
      call list(2,ibp2,2,2,2,kbp2,iwid,jwid,kwid,      &
          nvtxtmp,ijktmp2)
!
      call bc_noslip(nvtxtmp,ijktmp2,     &
          ul,vl,wl)
!
      call list(2,ibp2,jbp2,jbp2,2,kbp2,iwid,jwid,kwid,    &
          nvtxtmp,ijktmp2)
!
      call bc_noslip(nvtxtmp,ijktmp2,     &
          ul,vl,wl)
!
      endif
!
!
      refnorm=0.0
      totalvolume=0.0
      cellmass=0.0
!
      do n=1,ncells
         ijk=ijkcell(n)
!
         rho(ijk)=mc(ijk)/(vol(ijk)+1.e-10)
         cellmass=cellmass+mc(ijk)
!
         bmagx(ijk)=bmagx(ijk)/(vol(ijk)+1.e-10)
         bmagy(ijk)=bmagy(ijk)/(vol(ijk)+1.e-10)
         bmagz(ijk)=bmagz(ijk)/(vol(ijk)+1.e-10)
!
!        initialize the magnetic field
!
         bxn(ijk)=bmagx(ijk)
         byn(ijk)=bmagy(ijk)
         bzn(ijk)=bmagz(ijk)
!
         refnorm=refnorm+(bxn(ijk)**2+byn(ijk)**2+bzn(ijk)**2)
         totalvolume=totalvolume+vol(ijk)
!
!        compute the Alfven speed
!
         vasq(ijk)=(bxn(ijk)**2+byn(ijk)**2+bzn(ijk)**2)     &
          /(rho(ijk)+1.e-20)
!
         if(mc(ijk).ne.0.0) then
            sie(ijk)=sie1p(ijk)/(mc(ijk)+1.e-10)
         else
            sie(ijk)=0.0
         end if
!
      enddo
!
      refnorm=dsqrt(refnorm)/(totalvolume)**(.33)
!
!     set ghost cell values of b to zero
!
      call bc_ghost(ibp1+1,jbp1+1,kbp1+1,iwid,jwid,kwid,    &
          bxn,byn,bzn)
!
      zero=0.0
!
!      call torusbc(ibp2,jbp2,kbp2,     &
!          zero,                        &
!          bxn,byn,bzn)
!
      if(.not.cartesian) then
!    adjust rho on the k=2 and k=kbp1 boundaries
      fmo=1./(real(kbar)+0.5)
      do 60 i=1,ibp2
      do 60 j=1,jbp2
      ijkl=1+(i-1)*iwid+(j-1)*jwid
      ijkr=1+(i-1)*iwid+(j-1)*jwid+kbp1*kwid
      rho(ijkl+kwid)=rho(ijkl+kwid)     &
             *(mc(ijkl+kwid)-2.*mc(ijkl))/mc(ijkl+kwid)
      rho(ijkr-kwid)=rho(ijkr-kwid)     &
             *(mc(ijkr-kwid)+mc(ijkr)*fmo)/mc(ijkr-kwid)
   60 continue
!
      endif
!
!
!     zero working arrays
!
      call list(1,ibp2,1,jbp2,1,kbp2,iwid,jwid,kwid,     &
          nvtxtmp,ijktmp2)
!
      call SetZeroVector(nvtxtmp,ijktmp2,gradx,grady,gradz)

!     calculate the magnetic field intensity from the magnetization
!     by projection
!
      if(MAGNETIZED) then
!
!     1.  Calculate the divergence of M
!
      call bc_ghost(ibp1+1,jbp1+1,kbp1+1,iwid,jwid,kwid,   &
          bxn,byn,bzn)

      call divv(    &
         bxn,byn,bzn,srce)
!x
       call torusbc_scalar(ibp2,jbp2,kbp2,     &
          srce)
!
!
      call setzero(nvtx,ijkvtx,wate(1,7))
!
!     impose Dirichlet boundary conditions as default
!     by setting range of do loops
!
      istart=3
      istop=ibp1
      jstart=3
      jstop=jbp1
      kstart=3
      kstop=kbp1
!
      if(periodic_x) istart=2
      if(periodic_y) jstart=2
      if(periodic_z) kstart=2
!
      call list(istart,istop,jstart,jstop,kstart,kstop,     &
          iwid,jwid,kwid,      &
          nvtxtmp,ijktmp2)
!
      rnorm=0.0
      do n=1,nvtxtmp
         ijk=ijktmp2(n)
         rnorm=rnorm+dabs(srce(ijk))
      enddo

      phi=0.0

      do n=1,nvtxtmp
      ijk=ijktmp2(n)
!
!     set initial guess equal to source
!
!!!!      if(ncyc.eq.1) then
      phi(ijk)=srce(ijk)
!!!!      endif
      enddo
!
      call bc_scalar(ibp1+1,jbp1+1,kbp1+1,     &
          periodic_x,periodic_y,periodic_z,    &
          phi)
!
      itest = 1
      itmax=10
      itsub=5
      error=1.e-6
      if(itest .eq. 1)then
         if(rnorm.gt.error*refnorm) then
            call poisson_vtx(nvtxtmp,ijktmp2,refnorm,    &
               PRECON,    &
               itsub,error,rnorm,srce,    &
               bnorm,iter_pois)
         endif
!
      else
!
     itmax=100
      call poisson_cg(ncells,ijkcell,nvtx,ijkvtx,    &
          periodic_x,periodic_y,periodic_z,    &
          iwid,jwid,kwid,PRECON,    &
          ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,    &
          vol,vvol,vvolaxis,    &
          wate(1,10),wate(1,11),wate(1,12),itmax,error,srce,    &
          wate(1,3),wate(1,4),gradx,grady,gradz,wate(1,7),    &
          wate(1,9),wate(1,2))
!
      endif

!
!     calculate divergence free magnetic field
!
      call gradc(ncells,ijkcell,    &
          phi,gradx,grady,gradz)
!
      call InnerProductC(bxn,byn,bzn,bxn,byn,bzn,Half,Bsq)

      call InnerProductC(gradx,grady,gradz,gradx,grady,gradz,Half,GradPhisq)

      call InnerProductC(bxn,byn,bzn, gradx,grady,gradz,One,BGradPhi)
!
      do n=1,ncells
      ijk=ijkcell(n)
      bxn(ijk)=bxn(ijk)-gradx(ijk)
      bxl(ijk)=bxn(ijk)
      byn(ijk)=byn(ijk)-grady(ijk)
      byl(ijk)=byn(ijk)
      bzn(ijk)=bzn(ijk)-gradz(ijk)
      bzl(ijk)=bzn(ijk)
      enddo
!
      call divv(    &
         bxn,byn,bzn,srce)
!
             call torusbc_scalar(ibp2,jbp2,kbp2,     &
                srce)

      rnorm=0.0
      do n=1,nvtxtmp
         ijk=ijktmp2(n)
         rnorm=rnorm+dabs(srce(ijk))
      enddo


      call bc_ghost(ibp1+1,jbp1+1,kbp1+1,iwid,jwid,kwid,    &
          bxn,byn,bzn)
!
      call bc_ghost(ibp1+1,jbp1+1,kbp1+1,iwid,jwid,kwid,    &
          bxl,byl,bzl)
!
!
      endif
!
      return
      end subroutine vinit
