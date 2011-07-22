      program flip3d 
!
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use logical_com_M
      use blcom_com_M
      use nk_com_M
      use nkmhd_com_M
      use numpar_com_M
      use geometry_com_M
      use cophys_com_M
      use cindex_com_M
!      use gmres_com_M
      use cplot_com_M, ONLY : iplot
      use DXFiles
      use Timing

!
      character*80 :: name
      integer :: nvtxtmp,itsub,n,itlwd
      real(double) :: epsuv,NewVolume,    &
          dtcntr,totnrg,t3,dtsq,vvolaxis,xx,    &
          told,tmin,trl,psi
      real(double) :: Tstart_flip, Tfinish_flip, Tused_flip
!------------------------------------------------------------------
!               main program
!------------------------------------------------------------------
!
!ll   1. open files:
!        --------------
      call begin(name)
!
!ll   2. initialize run:
!        --------------
      call cpu_time(Tstart_flip)
      call initiv
      call OpenDXFiles
      call DXStart
!
      iplot=0
!      call plotinit(ibp2,jbp2,kbp2,name)
!
!
!     do j=1,3
!     do i=1,itdim
!     wate(i,j)=0.0
!     enddo
!     enddo
!     call debug(ncells,ijkcell,iwid,jwid,kwid,    &
!         cdlt,sdlt,strait,dz,    &
!         ibp1,jbp1,kbp1,nvtx,ijkvtx,wate(1,1),wate(1,2),wate(1,3),    &
!         x,y,z,vol,vvol,gradx,grady,gradz,divu)
!
      call parset
!
!     interpolate particle quantities to the grid
!
!
      call parcelc(bmagx,bmagy,bmagx)
!
!
       call parcelv
!      call parcelv(ncells,ijkcell,nvtx,ijkvtx,nrg,iwid,jwid,kwid,dt,    &
!          ijkctmp,    &
!          itdim,iphead,iphd2,link,    &
!          ico,mass,pxi,peta,pzta,up,vp,wp,    &
!          wate,    &
!          mv,umom,vmom,wmom,numberv,color)
!
!     ******************************************************************
!
!     generate an adaptive grid
!
!     *****************************************************************
!
      if(ADAPTIVEGRID) then
!
      call meshgen(ibp2,jbp2,kbp2,x,y,z,wgrid,vol,dt,taug,numgrid,    &
          WATE(1,7),WATE(1,8),WATE(1,9),    &
          WATE(1,10),WATE(1,11),WATE(1,12),    &
          WATE(1,13),WATE(1,14),WATE(1,15),    &
          LAMA,LAMB,LAMC,    &
          WATE(1,1),WATE(1,2),WATE(1,3),    &
          DELTA,PSI,    &
          WATE(1,4),WATE(1,5),WATE(1,6),    &
          ijkcell,ISTEP,IOTA,    &
          periodic,    &
          toroid,strait,rmaj,RWALL,    &
          BTORUS,Q0,    &
          cdlt,sdlt,dz)
!
      call celindex(2,ibp1,2,jbp1,2,kbp1,iwid,jwid,kwid,    &
          ijkcell,ijkctmp,ncells,ijkvtx,nvtxkm,nvtx)
!
      call metric(ncells,ijkcell,iwid,jwid,kwid,    &
          x,y,z,    &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,divpix)
!
      endif
!
!
!dak      call second(stim)
!dak      call getjtl (itlm)
!
!
      call budget
!
      itlwd=0
      t1=stim
      t2=t1
      trl=t1
!      tmin = (.02 + .005) * fibar*fjbar*fkbar
      tmin=0.0
!dak      tlm = itlm - tmin
      tlm = 300000.0
      write(*,*)'time limit, tlm=',tlm
      if (tlm .gt. 0) go to 1
!        (time limit too small for problem)
         write (kpt,9020) tmin
         write (kfm,9020) tmin
         stop 1
!
!ll   3. time-step: integrate equations in time:
!        --------------------------------------
    1 continue
!
!     >>>>>>>>>>>>>>> ncyc-loop <<<<<<<<<<<<<<<
!
      do ncyc= 1,ncyc2
!
      nh=mod(ncyc,nhst)
!       (calculate timestep)
        call timstp
!
!
   30    continue
!
         t=t+dt
         told=t2
!dak         call second(t2)
!
!      compute compute time/cycle/cell in millisec
!
         xx=(t2-told)*rbjbkb*1000.
         grind=xx
!     ****************************************************
!
!     prepare for the next computation cycle
!
!     ****************************************************
!
      call vinit
!
!         if(lpr.eq.0) go to 50
!
            if((ncyc).ge.1) then
               if(t+1.e-10.ge.tout) then
                  ixto=ixto+1
 !
!         call output
          call DXOutput
       iplot=iplot+1
                  if(dto(ixto).eq.0.0) ixto=ixto-1
                  tout=tout+dto(ixto)
               endif
             endif
!
!   50    continue
!
!     <<<<<<<<<<<<<Solve Equations of Motion>>>>>>>>>>>>>>>>>
!
!
      call strain(ncells,ijkcell,iwid,jwid,kwid,    &
          exx,exy,exz,eyy,eyz,ezz,    &
          ul,vl,wl,vol)
!
      call stress(ncells,ijkcell,   &
          pixx,pixy,pixz,piyy,piyz,pizz,    &
          exx,exy,exz,eyy,eyz,ezz,mu,lam)
!
!
      if(.not.MAGNETIZED) then
!
         if (IMPLICITPRESSURE) then
!
!----------------------------------
! implicit pressure solution
!----------------------------------
!
         do n=1,ncells
           p(ijkcell(n))=0.0
         end do
!
         call accel_3dmhd
!
         call strain(ncells,ijkcell,iwid,jwid,kwid,    &
           exx,exy,exz,eyy,eyz,ezz,    &
           ul,vl,wl,vol)
!
         do n=1,ncells
           ijk=ijkcell(n)
           divu(ijk)=exx(ijk)+eyy(ijk)+ezz(ijk)
         enddo
!
         call eos(p,rho,vol,sie,csq,gm1,ncells,ijkcell)
!
!       this version of the source is contained in implctp_cr.f
!
        numit=itmax
        call poisson(ncells,ijkcell,nvtxkm,nvtx,ijkvtx,    &
          iwid,jwid,kwid,PRECON,mv,    &
          ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,    &
          vol,vvol,vvolaxis,    &
          wate(1,1),wate(1,2),wate(1,3),numit,error,DT,CNTR,wate(1,4),    &
          rho,divu,csq,wate(1,5),wate(1,6),wate(1,7),wate(1,8),    &
          wate(1,9),wate(1,10),p,wate(1,11),wate(1,12))
!
         itsub=20
!
!
         dtsq=cntr*dt**2
         dtcntr=dt*cntr
!
!        do n=1,ncells
!c      ijk=ijkcell(n)
!      wate(ijk,26)=1./(csq(ijk)*rho(ijk)*dtsq+1.e-10)
!      wate(ijk,27)=p(ijk)*wate(ijk,26)-divu(ijk)/dtcntr
!c      enddo
!c
!      call poisson(ncells,ijkcell,nvtxkm,nvtx,ijkvtx,    &
!          iwid,jwid,kwid,PRECON,gm1,mv,    &
!          ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,    &
!          vol,vvol,vvolaxis,    &
!          itsub,numit,error,rnorm,DT,CNTR,wate(1,27),rho,divu,csq,    &
!          wate(1,26),wate(1,25),wate(1,24),wate(1,23),wate(1,1),bnorm,    &
!          gradx,grady,gradz,p,wate(1,22))
!
!
      else
!
!-------------------------------------------
!     explicit pressure
!-------------------------------------------
!
      call eos(p,rho,vol,sie,csq,gm1,ncells,ijkcell)
!
      endif
!
      call accel_3dmhd
      endif
!
      if(MAGNETIZED) then
!
!     solve implicit MHD equations
!
 1111 continue
!
      call eos(p,rho,vol,sie,csq,gm1,ncells,ijkcell)
!
      do n=1,ncells
!
      ijk=ijkcell(n)
!
      pl(ijk)=p(ijk)
      bxl(ijk)=bxn(ijk)
      byl(ijk)=byn(ijk)
      bzl(ijk)=bzn(ijk)
!
      enddo
!
      call list(1,ibp2,1,jbp2,1,kbp2,iwid,jwid,kwid,   &
          nvtxtmp,ijktmp2)
!
      call SetZeroVector(nvtxtmp,ijktmp2,jx,jy,jz)
!
      if(RESISTIVE) call resistive_diff
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      error=3.e-2
      eps=1.e-5
!
      call bc_ghost(ibp2,jbp2,kbp2,iwid,jwid,kwid,   &
          exx,exy,exz)
!
      call bc_ghost(ibp2,jbp2,kbp2,iwid,jwid,kwid,   &
          eyy,eyz,ezz)
!
!     check density array
!
!     calculate L1 norm of RHS
!
      ynorm=0.0

      do n=1,ncells
         ijk=ijkcell(n)
         ijktmp2(ijk)=n
         ynorm=ynorm    &
          +(abs(rho(ijk))+abs(bxn(ijk))+abs(byn(ijk))+abs(bzn(ijk)))    &
          *vol(ijk)
      enddo
!
!
      call mfnk_3dmhd(    &
          error,eps,    &
          numit,itmax)
!

!
      if(numit.gt.itmax) then
      t=t-0.5*dt
      dt=0.5*dt
      write(*,*) 'flip3d:  time step cut by itercg at ncyc=', ncyc
      go to 1111
      endif
!
      endif ! implicit MHD solver
!
!     compute strain on the grid
!
      call strain_ns(ncells,ijkcell,iwid,jwid,kwid,    &
          dudx,dudy,dudz,    &
          dvdx,dvdy,dvdz,    &
          dwdx,dwdy,dwdz,    &
          ul,vl,wl,vol)
!
      do n=1,ncells
         ijk=ijkcell(n)
         NewVolume=vol(ijk)    &
           *(1.d0+(dudx(ijk)+dvdy(ijk)+dwdz(ijk))*dt)
         dudx(ijk)=dudx(ijk)*vol(ijk)/NewVolume
         dudy(ijk)=dudy(ijk)*vol(ijk)/NewVolume
         dudz(ijk)=dudz(ijk)*vol(ijk)/NewVolume
         dvdx(ijk)=dvdx(ijk)*vol(ijk)/NewVolume
         dvdy(ijk)=dvdy(ijk)*vol(ijk)/NewVolume
         dvdz(ijk)=dvdz(ijk)*vol(ijk)/NewVolume
         dwdx(ijk)=dwdx(ijk)*vol(ijk)/NewVolume
         dwdy(ijk)=dwdy(ijk)*vol(ijk)/NewVolume
         dwdz(ijk)=dwdz(ijk)*vol(ijk)/NewVolume
      enddo
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dudx)
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dudy)
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dudz)
!
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dvdx)
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dvdy)
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dvdz)
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dwdx)
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dwdy)
!
      call bc_periodic(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          dwdz)
!
      call divphi_3dmhd
!
      call budget
!
!

      epsuv=1.d-10
!
      call parmov(epsuv)
!
      call parcelv
!      call parcelv(ncells,ijkcell,nvtx,ijkvtx,nrg,iwid,jwid,kwid,dt,    &
!          ijkctmp,    &
!          itdim,iphead,iphd2,link,    &
!          ico,mass,pxi,peta,pzta,up,vp,wp,    &
!          wate,    &
!          mv,umom,vmom,wmom,numberv,color)
!
!
      call parcelc(bmagx,bmagy,bmagx)
!
!dak         call second (t3)
    if(mod(ncyc,10).eq.1) then
   write(*,*) 'ncyc      t          dt   iterB iterP  IntlEnrg    KEnrg     BEnrg     TotalEnrg'
   write(*,*) ' '
    endif
      write(*,1131) ncyc,t,dt,numit,iter_pois, efnrg(nh), eknrg(nh), ebnrg(nh), efnrg(nh)+eknrg(nh)+ebnrg(nh)
      write(20,1131) ncyc,t,dt,numit,iter_pois, efnrg(nh), eknrg(nh), ebnrg(nh), efnrg(nh)+eknrg(nh)+ebnrg(nh)
 1131 format(i4,2f11.4,2i5,4f11.4)
    
         if (t.ge.twfin) then
            write(*,*) 'flip3d:  run ending'
            write(*,*) 'flip3d: t, twfin=',t,twfin
            goto 1010
         endif
!
         if (t2-t1.ge.dcttd)  go to 60
!
         if(t+1.e-10.lt.pttd) go to 70
!
         pttd=pttd+dpttd
!
   60    continue
         t1= t2
!
   70    continue
!
    enddo  ! main loop
!
!     >>>>>>>>>>>>>>> end of ncyc-loop <<<<<<<<<<<<<<<
!
      ncyc= ncyc2
!
!ll   4. terminate run:
!        -------------
 1010 continue
!
      tout = t
!
!  NEED GRAPHICS HERE
!
      if (lpr.ne.0) call DXOutput
      write(*,*) 'flip3d_nol:  iplot=',iplot
      call DXFinish
!
      write(*,*) 'flip3d_nol:  iplot=',iplot
!
      call cpu_time(Tfinish_flip)
      Tused_flip=Tfinish_flip-Tstart_flip
      call TimingOutput(Tused_flip)
      call CloseFiles
      stop
!
 9000 format(5h ncyc,i4,3h t=,1pe9.2,4h dt=,e9.2,8h courant,1pe9.2,   &
     9h lagrange,e9.2,   &
     7h numit=,i2,   &
      7h grind=,0pf7.3,1x,a1)
 9010 format (1x,i4,1pe9.2,1e11.3,1x,3e8.1,   &
        2e10.3,3e8.1,e10.3,i4,i3,e8.1,4x,e8.1)
 9020 format (34h>>> error: time limit must exceed ,1pe9.2,   &
             11hseconds <<<)
!
      end program flip3d
