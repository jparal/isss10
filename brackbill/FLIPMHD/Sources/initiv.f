!*************************************************************
!*                                                           *
      subroutine initiv
!*                                                           *
!*************************************************************
!
      USE vast_kind_param, ONLY:  double
      use logical_com_M
      use corgan_com_M
      use cindex_com_M
      use blcom_com_M
      use cindex_com_M
      use numpar_com_M
      use cophys_com_M
      use cplot_com_M
      implicit none
!

!
      real(double), dimension(8*nreg) :: xvitmp, yvitmp, zvitmp
      character ::  name*80
      integer :: itrlmp,kxi,nomnit,itlwd,itlwtt,ncr3,nwrtp,   &
          instab,nflm,itouch,nr,l
      real(double) :: sxth,thrd,wtp,wbtm,vbck,vfrnt,urght,ulft,a0fac,   &
          drbar,deigen(100),eigen(100),rmax,rota,amp,zk,zm,dtdum,   &
          del1,del2,del3,del4,del5,del6,del7,h
! 
!     NAMELIST Input
!
      namelist / celin / PRECON,ITMAX,itmag,cartesian,   &
      ADAPTIVEGRID, IMPLICITPRESSURE, MAGNETIZED, RESISTIVE,     &
      ibar,jbar,kbar,kplas,instab,nflm,EXPLICIT,error,precon,   &
      dx,dy,dz,   &
        eps,   &
         resist,resistivity,   &
       NOMAG,periodic_x,periodic_y,periodic_z,   &
          del1,del2,del3,del4,del5,del7,   &
      btorus ,bvertcl, efld0, q0,nphist,   &
       numgrid,taug,delta,   &
      tmax,tmax0,tmax1,tmax2,tmax3,   &
      lpr,itouch,nwrtp,ncyc2,   &
      a0,b0,   &
      dtdum,tlimd,dcttd,dpttd,   &
      dt,tlimd,dcttd,dpttd,twfin,dtmax,   &
      dto,dtoc,iout,   &
      lama,lamb,lamc,adaptg,newstuf,fields,   &
      rhoi,siei,beta,resist,   &
      eps,mgeom,mu,lam,   &
      rpl,rvac,rwall,rmaj,dz,iota,gx,gy,gz,gm1,imp   
!
      namelist / species / nsp, SpeciesName
!
      namelist / reg / nrg,npcelx,npcely,npcelz,  &
          rhr,siep,icoi,utordrft,qom,&
          bxr,byr,bzr,&
          xcenter,ycenter,zcenter,rex,rey,rez,rix,riy,riz,    &
          xvitmp,yvitmp,zvitmp,uvi,vvi,wvi,    &
          modex, modey, perturbx, perturby, perturbz
!
!l    *** additional declarations from initiv ***
!
      logical  cntrun
      data amp/0./,rota/0./,rmax/.95/
!
!---------------------------------------------------------------
!
!  default values:
      IMPLICITPRESSURE=.false.
      ADAPTIVEGRID=.false.
      RESISTIVE=.false.
      MAGNETIZED=.true.
      mgeom=1
      cntrun=.false.
      modex=0
      modey=0
      perturbx=0
      perturby=0
      perturbz=0
!
      imp=0
!
!     Cartesian geometry is the default
!
      cartesian=.true.
      resist=.false.
      resistivity=0.0d0
      eps=1.e-5
      itmag=25
      lama=0.
      lamb=0.
      lamc=0.
      adaptg=.false.
      fields=.false.
!dak      call dateh  (ddateh)
!dak      call timeh (hour)
      gm1= .6666667
      r4pi=1.
      dtpdv=1.e20
      periodic_x=.true.
      periodic_y=.true.
      periodic_z=.false.
!
!ll   1.  read parameters from cards:
!         --------------------------
!
      read(5,celin)
      write(59,celin)
!
!     default plot interval
      ixto=1
      if(dto(ixto).eq.0.) dto(ixto)=twfin/4.0
 
!
!
!ll   2.  check if continuation of previous run:
!         -------------------------------------
      if (ibar.lt.0 .and. jbar.eq.0) stop
      if (ibar.ne.0) goto 30
!
!l    2.1 continuation -- read parameters and all data from tape:
!         ------------------------------------------------------
      cntrun = .true.
!
!
!
!
!
   20 continue
!
      ncyc1= ncyc+1
!
      if (touch)  itouch= 1
!
!
      if (wrtp)   nwrtp = 1
!
      go to 40
!
!l    2.2 read further parameters from data cards:
!         ---------------------------------------
   30 continue
!
      read(5,species)
      write(59,species)
!
      read(5,reg)
      do l=1,8
         do nr=1,nrg
         xvi(l,nr)=xvitmp((nr-1)*8+l)
         yvi(l,nr)=yvitmp((nr-1)*8+l)
         zvi(l,nr)=zvitmp((nr-1)*8+l)
         enddo
      enddo

!
      write(*,*) 'initiv: dtmax=',dtmax
      strait=0.0
      toroid=1.0
      if(rmaj.ge.1.0e+05)  rmaj= 0.
      if(rmaj.eq.0.0)      strait=1.0
      if(rmaj.eq.0.0)      toroid=0.0
!
      touch= .false.
      if (itouch.gt.0)             touch= .true.
!
      wrtp= .false.
      if (nwrtp.gt.0)  wrtp= .true.
!
!
      kvac= kplas+2
      if(.not.cartesian) then
      dx  = rwall
      dy  = rmaj
      a   = rvac
      endif
!
      ncyc1= 1
      if (ncyc2.lt.0)  ncyc2= 77777
!
!
!
   40 continue
!
      if (tmax0.le.0.)  tmax0= 1.
      if (tmax1.le.0.)  tmax1= 1.
      if (tmax2.le.0.)  tmax2= 1.
      if (tmax3.le.0.)  tmax3= 1.
      if (tmax .le.0.)  tmax= (tmax0+tmax1+tmax2)/3.
!
!
!ll   3.  print input data
!         ----------------
      kpr= kpt
!
      write(*,*) 'initiv: writing input data'
      write (59,1111) name
 1111 format(a80)
      write (59,celin)
      write(59,reg)
         kpr= kfm
!
      write(*,*) 'initiv: cntrun=',cntrun
      if (cntrun)  go to 120
!
!
!ll   4.  calculate data for run-setup:
!         -----------------------------
      ibp1= ibar+1
      jbp1= jbar+1
      kbp1= kbar+1
      ibp2= ibar+2
      jbp2= jbar+2
      kbp2= kbar+2
      kvac= kplas+2
      kvp = kvac+1
      kvm = kvac-1
      iwid = 1
      jwid = ibp2*iwid
      kwid = jbp2*jwid
      iper = ibar*iwid
      jper = jbar*jwid
      fibar=real(ibar)
      rfibar= 1./fibar
      fjbar=real(jbar)
      rfjbar= 1./fjbar
      fkbar=real(kbar)
      rfkbar= 1./fkbar
      rbjbkb=1./(ibar*jbar*kbar)
      pi=acos(-1.)
      ilim= kbp2*kwid
!dak      call second(time)
      if (ilim .le. itdim) goto 105
      write(*,*) 'initiv:  ilim, itdim=',ilim, itdim
         stop 2
  105 continue
!
      if (pars.gt.0.)  itlwd= itlwd+ncr3
      itlwtt= 1.1*itlwd
!
      t= 0.
      write(*,*) 'initiv: t,dt=',t,dt
      mu= 0.
      lam= 0.
      tout=t+dto(ixto)
      pttd=t+dpttd
      ncyc=0
      numtd=0
      drbar=0.9*rwall
!      dtpos=dt*0.1
!      dt=dtpos
      dtc= dt
      dtv= dt
      nomnit= 10
      numit = 10
      a0fac=a0/(2.*(1.+a0**2))
      rdt=1./(dt+1.e-20)
      kxi=-1
      colamu=(1.+1.67)/(lam+mu+mu+1.e-10)
      urght=0.
      ulft=0.
      vfrnt=1.
      vbck=1.
      vfrnt=0.
      vbck=0.
      wbtm=0.
      wtp=0.
      thrd=1.0/3.0
      sxth=1.0/6.0
      itrlmp=25
!
!
      h = 2.0*pi/dz
!
!
  120 continue
!
!
!ll   5.  define mesh and initialize diagnostics:
!         ---------------------------------------
!
      call celindex(2,ibp1,2,jbp1,2,kbp1,iwid,jwid,kwid,    &
          ijkcell,ijkctmp,ncells,ijkvtx,nvtxkm,nvtx)
!
!
!l    5.1 calculate mesh:
!         --------------
      call mshset
      write(*,*) 'initiv: return from mshset,dt,t=',dt,t
!
!
      return
      end subroutine initiv
