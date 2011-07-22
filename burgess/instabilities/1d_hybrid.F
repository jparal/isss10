      program hybrid
c#######################################################################
c
c-----------------------------------------------------------------------
c23456789012345678921234567893123456789412345678951234567896123456789712
c-----------------------------------------------------------------------
c
c 1-D periodic version of a modern predictor-corrector hybrid code
c
c   *******************************************************************
c   *                                                                 *
c   * This code has been modified for ISSS-9 and is not a             *
c   * production version.  Last update: May 2009.                     *
c   *                                                                 *
c   * Code updates will be available online.                          *
c   *                                                                 *
c   * For questions, contact the author at varban@ssl.berkeley.edu    *
c   *                                                                 *
c   * A reference for this code is:                                   *
c   *                                                                 *
c   * Krauss-Varban, D.,                                              *
c   * From theoretical foundation to invaluable research tool:        *
c   * Modern hybrid simulations, Proceedings of the 7th International * 
c   * Symposium for Space Simulations (ISSS-7), pp. 15-18,            *
c   * Kyoto University, arXiv:physics/0610133, 2005.                  *
c   *                                                                 *
c   *******************************************************************
c
c                  --------------------------------
c                  Note for ISSS-9 Software Library
c                  --------------------------------
c
c The copyrights of the simulation codes in the ISSS-9 Software Library
c are retained by the author(s) of each simulation code.
c You are encouraged to contact the author(s) if you have any trouble
c in installing and/or applying the simulation code.  
c
c You are requested to acknowledge the author(s) of the code(s) if you
c publish the results of computer simulations using these codes or
c their modified versions, in scientific journals or otherwise.  We 
c recommend that you collaborate with the author(s) of the code(s) and 
c publish the simulation work jointly with them as much as possible, 
c at least, for the first work under the usage of these codes.
c
c-----------------------------------------------------------------------
c
c   - fast implementation of predictor-corrector algorithm
c   - substepping of e-m fields
c   - periodic boundary conditions, no injection 
c   - arbitrary number of species, constant number of particles
c
c-----------------------------------------------------------------------
c
c
      INCLUDE 'block'
c
      character text*20
c
c************
c    defaults for namelist variables
      data it /0/, t /0./
      data nsp /1/
      data nx /400/, xmax /400./, dtwci /.05/, resis /2.e-6/
      data text/' '/
      data npx /800,800/
      data bete /.5/
      data wpiwci /2000./, iterf /6/, icycle /0/
      data gammae/1.666666/
      data nwrtf /200/
      data ntimes /1/
      data theta /90./
      data iseed/0/
      namelist /datum/ ntimes,xmax,dtwci,wspec,nwrtf
     1 ,resis,bete,npx,nx,anspec
     2 ,frac,vxspec,iterf,wpiwci,qspec,vyspec,vzspec
     3 ,btspec,gammae,nsp,theta,ibxy,nspecm
     4 ,intwrt,text
      ns=nss
      nc=ncx
c
c****************
c initialize species-dependent defaults
c****************
      do  5 k = 1, ns
        anspec(k) = 1.0
        btspec(k) = 1.0
        frac(k) = 1.0
        qspec(k) = 1.0
        vxspec(k) = 0.0
        vyspec(k) = 0.0
        vzspec(k) = 0.0
        wspec(k) = 1.0
  5   continue      
c****************
c    open files and read namelist
c****************
      open (7,file='output.txt',form='formatted',access='sequential')
      open (5,file='input.txt',form='formatted',status='old',
     &      access='sequential')
c make sure the current platform supports as large an integer as needed
c else declare long integer
c the following assumes IEEE format for these files (4 byte)
      lrec =  4 *ncx
      lrec1 = 4 *(ncx-1)
      lrec2 = 4 *(ncx-2)
c Normally, write unformatted direct access files.
c Here, text files for GnuPlot access.
c
c     open (4,file='bx.gda',form='unformatted',
c    &      access='direct', recl=lrec2)
c
      open (8,file='by.dat',form='formatted',
     &      access='sequential')
      open (9,file='bz.dat',form='formatted',
     &      access='sequential')
      open (12,file='dn.dat',form='formatted',
     &      access='sequential')
      open (13,file='tpar2.dat',form='formatted',
     &      access='sequential')
      open (14,file='tper2.dat',form='formatted',
     &      access='sequential')
      open (15,file='tpar1.dat',form='formatted',
     &      access='sequential')
      open (16,file='tper1.dat',form='formatted',
     &      access='sequential')
      open (24,file='vx.dat',form='formatted',
     &      access='sequential')
      open (25,file='vy.dat',form='formatted',
     &      access='sequential')
      open (26,file='vz.dat',form='formatted',
     &      access='sequential')
c
      open (31,file='epara1.dat',form='formatted',
     &      access='sequential')
      open (32,file='eperp1.dat',form='formatted',
     &      access='sequential')
      open (33,file='epara2.dat',form='formatted',
     &      access='sequential')
      open (34,file='eperp2.dat',form='formatted',
     &      access='sequential')
c
      read (5,datum)
c***************************
c    check dimensions
c***************************
      npartm=0
      npart=0
      do k=1,ns
        if(nspecm(k).lt.npx(k)) call endrun ('npartm')
        npartm=npartm+nspecm(k)
        npart=npart+npx(k)
      enddo
      if (npartm.gt.nbb) call endrun ('nbb too small')
      if(nsp.gt.ns) call endrun('species')
      if((nx+2).gt.nc) call endrun ('ncc')
c
c******************
c    could attach date and time stamp to label 'title', here
c******************
c
      title =  ' ' // text // ' '
      write(6,*) title
c
c**********************
c    initialize variables etc.
c**********************
      idarec=0
      call init
      it0=1
c
c**********************
c     main time loop
c**********************
c
      do it = it0,ntimes
c        write(6,*) it,it0
c
        t=t+dt
c
c*******************************
c        start subcycle loop
c
        do icycle=1, iterf
c
c - - - Predictor-Corrector, 2 particle pushes per large time step - -
c     N +1/2 step
          ipc = 1
          if(icycle.eq.1)  then
            call trans
            if (  mod(it,nwrtf).eq.0) call pptmp
          endif
          call field
c      N +1   predictor
          ipc = 2
          call field
c      N +3/2 predictor
          ipc = 3
          if(icycle .eq. iterf) then
            call trans
          endif
          call field
c      N +1   corrector
          ipc = 4
          call field
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        enddo
c
c        end subcycle loop
c*******************************
c
        call diag2
c
        if(mod(it,intwrt).eq.0) write(6,*) it
c
      enddo
c
c**********************
c  end main time loop
c  close any files etc.
c*********************
c
      write(*,*) "regular end"
      call exit
c
      end
c#######################################################################
      subroutine init
c#######################################################################
c
c*************
c
      INCLUDE 'block'
c
c***************
c    define some useful variables from input
c***************
c
      pifac=180./3.1415926
      nx1=nx+1
      nx2=nx+2
      hx=xmax/float(nx)
      hxi=1./hx
c
c      iterf must be even or one:
c
      if( (mod(iterf,2) .ne. 0) .and. (iterf .ne. 1) )  then
        write(*,*) " iterf must be even or one.  iterf = ", iterf
        stop
      endif
c
      ntot=0
      do 10 k=1,nsp
c
      nspec0(k)=npx(k)
c
      tx0(k)=btspec(k)/(2.*wpiwci**2)
      vbx(k)=vxspec(k)/wpiwci
      vby(k)=vyspec(k)/wpiwci
      vbz(k)=vzspec(k)/wpiwci
c
      x0=0.0
      x1=xmax
c
   10 ntot=ntot+nspec0(k)
      te0=bete/(2.*wpiwci**2)
      bxc=cos(theta/pifac)/wpiwci
      bzc=sin(theta/pifac)/wpiwci
c
      dt=wpiwci*dtwci
c
c********************
c   set up B (byc,bzc)  based on
c   whether B is in the xz (ibxy=0) or the xy plane
c*******************
      if(ibxy.eq.0) then
        write(6,*) 'can t have main field out of plane here'
        stop
      else
        byc=bzc
      endif
c*******************
c  write out input variables
      write (7,110)
      write (7,120) nx
      write (7,130) xmax
      write (7,160) dtwci
      write (7,170) wpiwci
      write (7,180) resis
      write (7,190) bete
      write (7,200)
      write (7,210) (is,frac(is),vxspec(is),btspec(is),anspec(is),
     &               wspec(is),nspec0(is),is=1,nsp)
      write (7,220) theta
      write (7,290) nwrtf
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      denmin = 0.0
c
      nparh(1)=1
      nspec2(1)=npx(1)
      do is=2,nss
        nparh(is)=nparh(is-1)+nspecm(is-1)
        nspec2(is)=nparh(is)+npx(is)-1
      enddo
c
      kspc=0
      isp =0
c
20    call getb (isp,l3,l4,nspec2,nss,nparh)
      if (isp.ne.-1) then
      if (kspc.ne.isp) then
      nppx=npx(isp)
      pinvx=1./nppx
c**********************
c here, T is defined as the total temperature
c
      cc = sqrt(3.0/(2.0*anspec(isp) +1.0))
      vpar(isp)=sqrt(btspec(isp) *cc/(wspec(isp)*wpiwci**2))
      vper(isp)=vpar(isp)*sqrt(anspec(isp))
      dfac(isp)=(x1-x0)/(hx*float(nspec0(isp)))
c
      endif
40    do 60 l=l3,l4
      if (nppx.le.0) go to 60
c 
      cs=1.-pinvx*(nppx-0.5)
      x(l)=x0+(x1-x0)*cs
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
      vmag=sqrt(-alog(1.-.999999*random(iseed)))
      th=6.2831853*random(iseed)
c*****************************
c Initialize particle velocity with possible temperature
c anisotropy. Note that for the periodic case, here vbx
c is added to vx and vy since the drift is along B
c*****************************
c NOTE: Not valid when field rotates ...
c       here: input vbx is really v_parallel
c*************
c parallel velocity :
c
      vxa=vpar(isp)*vmag*cos(th) +vbx(isp)
      vmag=sqrt(-alog(1.-.999999*random(iseed)))
      th=6.2831853*random(iseed)
c
c the two perpendicular velocities :
c
      vya=vper(isp)*vmag*cos(th)
      vza=vper(isp)*vmag*sin(th)
c
c vx is the same independent of whether B lies in X-Z
c or X-Y plane. Initialize vx :
c
      vx(l)=vxa*cos(theta/pifac)-vza*sin(theta/pifac)
c
c vy and vz are initialized based on the plane of B.
c
      if(ibxy.eq.0) then
c
c B is in the X-Z plane :
c
         vy(l)=vya+vby(isp)
         vz(l)=vza*cos(theta/pifac)+vxa*sin(theta/pifac)+vbz(isp)
      else
c 
c B is in the X-Y plane :
c
         vz(l)=vya+vbz(isp)
         vy(l)=vza*cos(theta/pifac)+vxa*sin(theta/pifac)+vby(isp)

      endif
c
c***********************************
c
   60 nppx=nppx-1
c
c --- end per species initialization
c
      kspc=isp
      go to 20
c
c --- end all species initialization
c
      endif
c
c     set the minimum allowed density (e.g., 1/2 particle per cell)
c
      denmin = 0.5 *dfac(1)
      write(7,*) ' minimum density:      ', denmin
c
c***********************
c    initialize field arrays
c***********************
c
c    fields are initialized throughout domain
c
c***********************
c   set the fields in the entire domain
c   e-fields are not defined here ...
c***********************
c
      do i=1,nx2
          bz(i) = 0.0
          by(i) = byc 
c         bx(i) = bxc  constant in 1-D from div B = 0
          ex(i) = 0.0
          ey(i) = 0.0
          ez(i) = 0.0
      enddo
c
c***********************
c   set the friction force and resitivity
c   friction force can be zero at t=0
c**********************
c
      do 100 i=1,nx2
c       fox(i)=0.     no resistive force in x in 1-D
        foy(i)=0.
        foz(i)=0.
c
        eta(i) = resis
c
  100 continue
c                  set den and curl to zero so it is defined
      do i=1,nx2
          den(i) = 0.0
c         curlex(i)=0.0    there is no B_x advance in 1-D
          curley(i)=0.0
          curlez(i)=0.0
      enddo
c******************
c    run one time step with dt=0 to initialize fields and diagnostics
c******************
      write (7,300)
      ipc = 0
      dtsav=dt
      dt=0.
      call trans
      call pptmp
      write (7,310)
      dt=dtsav
      write (7,315)
      call field
      write (7,320)
      call diag2
      return
c
  110 format (' namelist input parameters/values:      ')
  120 format ('       nx=   no. of cells in x          ',i15)
  130 format ('     xmax=   system length in x         ',f15.1)
  160 format ('    dtwci=   time step                  ',f15.3)
  170 format ('   wpiwci=   wpi/wci                    ',f15.1)
  180 format ('    resis=   resistivity                ',e15.5)
  190 format ('     bete=   electron beta              ',f15.3)
  200 format ('   species','   density','     vb/va','      beta','   an
     1isotropy','   mi/mp',' #particles')
  210 format (i10,f10.3,f10.2,f10.3,f10.2,f10.1,i10)
  220 format ('    theta=                              ',f15.2)
  290 format ('    nwrtf=   increment of field writes  ',i15)
  300 format ('  init call to trans  ')
  310 format ('  init call to field  ')
  315 format ('  init call to inject ')
  320 format ('  init call to diagnos')
      end
c
c#######################################################################
      subroutine trans
c#######################################################################
c
c*************
      INCLUDE 'block'
c************
c      zero source arrays
c      remember density and velocities of previous substep
c      (N or N+1, respectively)
c      write(*,*) "icycle = ", icycle, ",  ipc = ", ipc
      do 20 i=1,nx2
        if(ipc .eq. 3) then
c          deno(i)=den(i)      not in this sub-step
          vixo(i)=vix(i)
          viyo(i)=viy(i)
          vizo(i)=viz(i)
        else
          deno(i)=den(i)
          den(i)=0.
        endif
        vix(i)=0.
        viy(i)=0.
        viz(i)=0.
   20 continue
      do 30 k=1,nsp
      do 30 i=1,nx2
        if(ipc .ne. 3)   dns(i,k)=1.e-20
        vxs(i,k)=0.
        vys(i,k)=0.
        vzs(i,k)=0.
   30 continue
c******************
c in non-periodic systems, 
c allow for injection of particles in each time step.
c
c      only inject and recalculate species numbers in first part
c      of predictor-corrector:
c
      if((ipc .eq. 0) .or. (ipc .eq. 1)) then
c
c prepare properties of newly injected particles:
c injection rate could be calculated here
c
c******************
c  
c push particles and plot phase space
c
c               read and write in parmov1; need moments at N+1/2
c               density at N+1/2 via average old, new
        isp =0
40      call getb (isp,l3,l4,nspec2,nss,nparh)
        if (isp.ne.-1) then
          call parmov1
c         call diag1  --- particle diagnostics could go here
          go to 40
        endif
c
c - - - - - - - - - - 
c
      else
c               read only in parmov2; need v moments at N+1
c               moments at larger than N+1 via extrapolation
        isp =0
48      call getb (isp,l3,l4,nspec2,nss,nparh)
        if (isp.ne.-1) then
          call parmov2
          go to 48
        endif
      endif
c
c - - - - - - - - - - 
c
c******************
c    put ghost cell contributions back into real domain
c    
c    This version is periodic in x.
c
      if(ipc .ne. 3)  then
c                             do regular case, including density:
        do 80 k=1,nsp
            dns(2,k)=dns(2,k)+dns(nx2,k)
            vxs(2,k)=vxs(2,k)+vxs(nx2,k)
            vys(2,k)=vys(2,k)+vys(nx2,k)
            vzs(2,k)=vzs(2,k)+vzs(nx2,k)
            dns(nx1,k)=dns(nx1,k)+dns(1,k)
            vxs(nx1,k)=vxs(nx1,k)+vxs(1,k)
            vys(nx1,k)=vys(nx1,k)+vys(1,k)
            vzs(nx1,k)=vzs(nx1,k)+vzs(1,k)
            dns(1,k)=dns(nx1,k)
            vxs(1,k)=vxs(nx1,k)
            vys(1,k)=vys(nx1,k)
            vzs(1,k)=vzs(nx1,k)
            dns(nx2,k)=dns(2,k)
            vxs(nx2,k)=vxs(2,k)
            vys(nx2,k)=vys(2,k)
            vzs(nx2,k)=vzs(2,k)
   80   continue
c
c    normalize density and velocities of each species
c    note that velocities are 1/2 step off, so there is not
c    an exact match  with the densities
c
        do 90 k=1,nsp
            do 90 i=1,nx2      
              denavg=dns(i,k)
              if(denavg .lt. 0.2)   denavg = 1.0e10
              vxs(i,k)=vxs(i,k)/denavg
              vys(i,k)=vys(i,k)/denavg
              vzs(i,k)=vzs(i,k)/denavg
              dns(i,k)=dns(i,k)*dfac(k) *frac(k)
   90   continue
c    calculate total density
        do 100 k=1,nsp
            do 100 i=1,nx2
              den(i)=den(i)+dns(i,k)*qspec(k)
  100   continue
c
          do 110 i=1,nx2
            if (den(i).lt.denmin) den(i)=denmin
  110   continue
c
      else
c            do not re-calculate density moment
        do k=1,nsp
            vxs(2,k)=vxs(2,k)+vxs(nx2,k)
            vys(2,k)=vys(2,k)+vys(nx2,k)
            vzs(2,k)=vzs(2,k)+vzs(nx2,k)
            vxs(nx1,k)=vxs(nx1,k)+vxs(1,k)
            vys(nx1,k)=vys(nx1,k)+vys(1,k)
            vzs(nx1,k)=vzs(nx1,k)+vzs(1,k)
            vxs(1,k)=vxs(nx1,k)
            vys(1,k)=vys(nx1,k)
            vzs(1,k)=vzs(nx1,k)
            vxs(nx2,k)=vxs(2,k)
            vys(nx2,k)=vys(2,k)
            vzs(nx2,k)=vzs(2,k)
        enddo
c
c    normalize velocities of each species
c    note that velocities are 1/2 step off, so there is not
c    an exact match with the densities, needs protection here
c    against collected v-moment without (or with tiny) n-moment
c
        do k=1,nsp
            do i=1,nx2      
              denavg=dns(i,k)/(dfac(k) *frac(k))
              if(denavg .lt. 0.2)   denavg = 1.0e10
              vxs(i,k)=vxs(i,k)/denavg
              vys(i,k)=vys(i,k)/denavg
              vzs(i,k)=vzs(i,k)/denavg
            enddo
        enddo
c
      endif
c
c    calculate total ion velocity
      do 140 k=1,nsp
      do 130 i=1,nx2
c
c at N+1/2 (IPC = 1), use averaged density as approximation
c for true moment matching time and place of v-eval:
c
      if(ipc .eq. 1)   then
        denavg = 0.5 *(deno(i) +den(i))
      else
        denavg = den(i)
      endif
      vix(i)=vix(i)+qspec(k)*dns(i,k)*vxs(i,k)/denavg
      viy(i)=viy(i)+qspec(k)*dns(i,k)*vys(i,k)/denavg
      viz(i)=viz(i)+qspec(k)*dns(i,k)*vzs(i,k)/denavg
  130 continue
  140 continue
c
      if (it.eq.0) then
          do i=1,nx2
            deno(i)=den(i)
            vixo(i)=vix(i)
            viyo(i)=viy(i)
            vizo(i)=viz(i)
          enddo
      endif
c
c************
c
      return
      end
c
c#######################################################################
      subroutine parmov1
c#######################################################################
c
c******
      INCLUDE 'block'
c******
      h=dt*qspec(isp)/wspec(isp)
      hh=.5*h
      dt2 = 0.5 *dt
c
      do 100 l=l3,l4
c
c     could inject recycled particles here (non-periodic systems)
c
c*******************
c   calculate particle weights for magnetic field
c*******************
c   changed to 1.99999 from 1.9999999999 (PC limited accuracy)
      rxb=hxi*x(l)+1.99999
      ixb=rxb
      ixb1=ixb+1
      fxb=rxb-ixb
      w1b=(1.-fxb)
      w2b=fxb
c*******************
c   calculate B force on the particle
c*******************
      bxa=bxc
      bya=by(ixb)*w1b+by(ixb1)*w2b
      bza=bz(ixb)*w1b+bz(ixb1)*w2b
c*******************
c   calculate particle weights for electric field
c*******************
      rxe=hxi*x(l)+1.50000001
      ixe=rxe
      ixe1=ixe+1
      fxe=rxe-ixe
      w1e=(1.-fxe)
      w2e=fxe
c******************
c   calculate E and drag force on the particle
c******************
      exa=ex(ixe)*w1e+ex(ixe1)*w2e
      eya=ey(ixe)*w1e+ey(ixe1)*w2e+foy(ixe)*w1e+foy(ixe1)*w2e
      eza=ez(ixe)*w1e+ez(ixe1)*w2e+foz(ixe)*w1e+foz(ixe1)*w2e
c*****************
c   calculate velocity
c*****************
      ff=2./(1.+hh*hh*(bxa**2+bya**2+bza**2))
      vex=vx(l)+exa*hh
      vey=vy(l)+eya*hh
      vez=vz(l)+eza*hh
      p2xs=vex+(vey*bza-vez*bya)*hh
      p2ys=vey+(vez*bxa-vex*bza)*hh
      p2zs=vez+(vex*bya-vey*bxa)*hh
      vx(l)=vex+ff*(p2ys*bza-p2zs*bya)*hh+exa*hh
      vy(l)=vey+ff*(p2zs*bxa-p2xs*bza)*hh+eya*hh
      vz(l)=vez+ff*(p2xs*bya-p2ys*bxa)*hh+eza*hh
c*****************
c    push the ions
c*****************
      x(l)=x(l)+dt*vx(l)
c*****************
c    put the particles outside of the box back in.
c*****************
c
        if(x(l).ge.xmax)  then
          x(l) = x(l) -xmax
        endif
        if(x(l).le.0.)    then
          x(l) = xmax +x(l)
        endif
c
  100 continue
c*****************
c
c   inject new particles here that have not been taken care
c   of via recycling.
c
c***********************
c    gather particle density at new position
c    and velocity at 1/2 step
c***********************
      do 210 l=l3,l4
      if((x(l) .le. 0.0).or.(x(l) .ge. xmax)) go to 210
c
      xl=x(l) -dt2*vx(l)
c
      if((xl .le. 0.0).or.(xl .ge. xmax))  then
c                 make sure a moment is still evaluated
        xl = x(l)
      endif
      rx=hxi*xl+1.50000001
      ix=rx
      ix1=ix+1
      fx=rx-ix
      w1=(1.-fx)
      w2=fx
      vxs(ix, isp)=w1*vx(l)+vxs(ix, isp)
      vxs(ix1,isp)=w2*vx(l)+vxs(ix1,isp)
      vys(ix, isp)=w1*vy(l)+vys(ix, isp)
      vys(ix1,isp)=w2*vy(l)+vys(ix1,isp)
      vzs(ix, isp)=w1*vz(l)+vzs(ix, isp)
      vzs(ix1,isp)=w2*vz(l)+vzs(ix1,isp)
c
      rx=hxi*x(l)+1.50000001
      ix=rx
      ix1=ix+1
      fx=rx-ix
      w1=(1.-fx)
      w2=fx
      dns(ix, isp)=w1+dns(ix, isp)
      dns(ix1,isp)=w2+dns(ix1,isp)
c
 210  continue
 110  continue
      return
      end
c
c#######################################################################
      subroutine parmov2
c#######################################################################
c
c******
      INCLUDE 'block'
c******
c        for second part of predictor-corrector
c        read only, do not change values nor inject
c        push one-half step
C       here: push full step, but then cut in half
C       to have good interpolation values
c******
      h=dt*qspec(isp)/wspec(isp)
c******
      hh=.5*h
      do 100 l=l3,l4
      if((x(l).le. 0.0).or.(x(l).ge. xmax)) then
        goto 100
      endif
c*******************
c   calculate particle weights for magnetic field
c*******************
c  changed to 1.99999 from 1.9999999999 (PC limited accuracy)
      rxb=hxi*x(l)+1.99999
      ixb=rxb
      ixb1=ixb+1
      fxb=rxb-ixb
      w1b=(1.-fxb)
      w2b=fxb
c*******************
c   calculate B force on the particle
c*******************
      bxa=bxc
      bya=by(ixb)*w1b+by(ixb1)*w2b
      bza=bz(ixb)*w1b+bz(ixb1)*w2b
c*******************
c   calculate particle weights for electric field
c*******************
      rxe=hxi*x(l)+1.50000001
      ixe=rxe
      ixe1=ixe+1
      fxe=rxe-ixe
      w1e=(1.-fxe)
      w2e=fxe
c******************
c   calculate E and drag force on the particle
c******************
      exa=ex(ixe)*w1e+ex(ixe1)*w2e
      eya=ey(ixe)*w1e+ey(ixe1)*w2e+foy(ixe)*w1e+foy(ixe1)*w2e
      eza=ez(ixe)*w1e+ez(ixe1)*w2e+foz(ixe)*w1e+foz(ixe1)*w2e
c*****************
c   calculate velocity (note differences to parmov1)
c*****************
c     ff=2./(1.+hh*hh*(bxa**2+bya**2+bza**2))
      ff=1./(1.+hh*hh*(bxa**2+bya**2+bza**2))
      vex=vx(l)+exa*hh
      vey=vy(l)+eya*hh
      vez=vz(l)+eza*hh
      p2xs=vex+(vey*bza-vez*bya)*hh
      p2ys=vey+(vez*bxa-vex*bza)*hh
      p2zs=vez+(vex*bya-vey*bxa)*hh
c     vx(l)=vex+ff*(p2ys*bza-p2zs*bya)*hh+exa*hh
c     vy(l)=vey+ff*(p2zs*bxa-p2xs*bza)*hh+eya*hh
c     vz(l)=vez+ff*(p2xs*bya-p2ys*bxa)*hh+eza*hh
      vex=vex+ff*(p2ys*bza-p2zs*bya)*hh
      vey=vey+ff*(p2zs*bxa-p2xs*bza)*hh
      vez=vez+ff*(p2xs*bya-p2ys*bxa)*hh
c*****************
c    push the ions not even one half step and
c    gather particle velocity only, at current position
c*****************
      ix=ixe
      ix1=ix+1
      vxs(ix, isp)=w1e*vex+vxs(ix, isp)
      vxs(ix1,isp)=w2e*vex+vxs(ix1,isp)
      vys(ix, isp)=w1e*vey+vys(ix, isp)
      vys(ix1,isp)=w2e*vey+vys(ix1,isp)
      vzs(ix, isp)=w1e*vez+vzs(ix, isp)
      vzs(ix1,isp)=w2e*vez+vzs(ix1,isp)
c
 100  continue
      return
      end
c
c#######################################################################
      subroutine diagnos
c#######################################################################
c
c*************
      INCLUDE 'block'
c************
      entry diag1
c************
c    place desired particle diagnostics here
c************
      return
c
      entry diag2
c***********
c    go through this part after sources and fields are known
c***********
c
c    write b fields and density to disk files
c
      if (nwrtf.le.0) go to 130
      if (mod(it,nwrtf).ne.0) go to 130
c
      idarec = idarec +1
c
c usually, direct access files are written instead of text files...
c     write (4, rec=idarec)  
c    &      (bxc*wpiwci,         i = 2, nx1)
c
c normally, last data point is not written 
c for periodic system (same as first)...
c
c here: simply scale subsequent time slices
c
      scale = 2.0 *float(it/nwrtf)
c
      do i = 2, nx2
        xi = float(i-2) *hx
        write (8,*)  xi, '  ', it,  '  ', by(i)*wpiwci
        write (9,*)  xi, '  ', it,  '  ', bz(i)*wpiwci
        write (12,*) xi, '  ', it,  '  ', den(i)
        write (13,*) xi, '  ', it,  '  ', tpar(i,2)
        write (14,*) xi, '  ', it,  '  ', tperp(i,2)
        write (15,*) xi, '  ', it,  '  ', tpar(i,1)
        write (16,*) xi, '  ', it,  '  ', tperp(i,1)
        write (24,*) xi, '  ', it,  '  ', vix(i)*wpiwci
        write (25,*) xi, '  ', it,  '  ', viy(i)*wpiwci
        write (26,*) xi, '  ', it,  '  ', viz(i)*wpiwci
c
      enddo
c
c     write (8,*)  '  ' 
c     write (9,*)  '  '
c     write (12,*) '  '
c     write (13,*) '  '
c     write (14,*) '  '
c     write (24,*) '  '
c     write (25,*) '  '
c     write (26,*) '  '
c
c     write (99, rec=idarec) 
c    &      (ez(i)*wpiwci**2,  i = 2, nx1)

c     write (15, rec=idarec) 
c    &      (tpar(i,2),   i = 2, nx1)
c     write (16, rec=idarec) 
c    &      (tperp(i,2),   i = 2, nx1)
c
  130 continue
      return
c
      end
c***********************************************************************
c
c#######################################################################
      subroutine field
c#######################################################################
c
c=======================================================================
c
c************
      include "block"
c************
c
c=======================================================================
c
c       define substep timestep
c
      odt = dt
      dt = dt /float(iterf)
c
c       - order and tasks depend on predictor-corrector step ipc - 
c
      if(ipc .eq. 0) then
c       call bcalc
        call ecalc
        call focalc
      elseif(ipc .eq. 1) then
        call bcalc
        call ecalc  
      elseif(ipc .eq. 2) then
        call ecalc
        call bcalc
        call focalc
      elseif(ipc .eq. 3) then
        call bcalc
        call ecalc 
      elseif(ipc .eq. 4) then
        call ecalc
        call bcalc
        call focalc
      else
        write(*,*) ' *** incorrect ipc in field ***' 
        stop
      endif 
c
      dt = odt
c
c=======================================================================
c
      return
      end
c
c#######################################################################
      subroutine ecalc
c#######################################################################
c
c************
      include "block"
c************
c
c calculation depends on step within  predictor-corrector sequence
c
c - - - - - - - - - -
c
      if(ipc .eq. 1) then
        w2 = float(iterf +2*icycle -1) /float(2 *iterf)
        w1 = 1.0 -w2
        v2 = float(2*icycle -1) /float(2 *iterf)
        v1 = 1.0 -v2
        do i=1,nx2
             denave = v1*deno(i) +v2*den(i)
             pe(i) = te0 *(denave **gammae)
        enddo
        do i=2,nx1
          by1=byt(i+1)
          by2=byt(i  )
          bz1=bzt(i+1)
          bz2=bzt(i  )
          vixa = w1*vixo(i) +w2*vix(i)
          viya = w1*viyo(i) +w2*viy(i)
          viza = w1*vizo(i) +w2*viz(i)
          dena = v1*deno(i) +v2*den(i)
          a=1./dena
          dxa=a/hx
          dbydx=by1-by2
          dbzdx=bz1-bz2
          bxav=bxc
          byav=.5*(by1+by2)
          bzav=.5*(bz1+bz2)
          dpedx=0.5*(pe(i+1)-pe(i-1))
          ext(i)=(viza*byav -viya*bzav)
     1               -dxa*(dbydx*byav +dbzdx*bzav+dpedx)
          eyt(i)=(vixa*bzav -viza*bxav)
     1               +dxa*dbydx*bxav -eta(i)*dbzdx*hxi
          ezt(i)=(viya*bxav -vixa*byav)
     1               +dxa*dbzdx*bxav +eta(i)*dbydx*hxi
        enddo
c
c******************************************
c
c     boundary conditions
c
c**********************
c     set E at the left hand boundary
c**********************
        ext(1)=ext(nx1)
        eyt(1)=eyt(nx1)
        ezt(1)=ezt(nx1)
c**********************
c     set E at the right hand boundary
c**********************
        ext(nx2)=ext(2)
        eyt(nx2)=eyt(2)
        ezt(nx2)=ezt(2)
c****************************
c
      elseif(ipc .eq. 2) then
c                              use extrapolation for N+1 predictor
          do i=1,nx2
            ex(i) = -ex(i) +2.0 *ext(i) 
            ey(i) = -ey(i) +2.0 *eyt(i)
            ez(i) = -ez(i) +2.0 *ezt(i)
          enddo
c
c      - boundary conditions are automatically satisfied -
c      - need curl E for this substep: - 
c
          do i=2,nx2
            deydx=ey(i) -ey(i-1)
            dezdx=ez(i) -ez(i-1)
            curley(i)=-dezdx*hxi
            curlez(i)= deydx*hxi
          enddo
c
      elseif((ipc .eq. 0) .or. (ipc .eq. 3)) then
        w2 = float(iterf +2*icycle +1) /float(2 *iterf)
          if(icycle .eq. iterf) w2 = float(2*icycle +1) /float(2 *iterf)
        w1 = 1.0 -w2
        v2 = float(2*icycle +1) /float(2 *iterf)
        v1 = 1.0 -v2
          do i=1,nx2
             denave = v1*deno(i) +v2*den(i)
             pe(i) = te0 *(denave **gammae)
          enddo
        do i=2,nx1
          by1=by(i+1)
          by2=by(i  )
          bz1=bz(i+1)
          bz2=bz(i  )
          vixa = w1*vixo(i) +w2*vix(i)
          viya = w1*viyo(i) +w2*viy(i)
          viza = w1*vizo(i) +w2*viz(i)
          dena = v1*deno(i) +v2*den(i)
          a=1./dena
          dxa=a/hx
          dbydx=by1-by2
          dbzdx=bz1-bz2
          bxav=bxc
          byav=.5*(by1+by2)
          bzav=.5*(bz1+bz2)
          dpedx=0.5*(pe(i+1)-pe(i-1))
          ex(i)=(viza*byav -viya*bzav)
     1               -dxa*(dbydx*byav +dbzdx*bzav+dpedx)
          ey(i)=(vixa*bzav -viza*bxav)
     1               +dxa*dbydx*bxav -eta(i)*dbzdx*hxi
          ez(i)=(viya*bxav -vixa*byav)
     1               +dxa*dbzdx*bxav +eta(i)*dbydx*hxi
        enddo
c
c******************************************
c
c     boundary conditions
c
c**********************
c     set E at the left hand boundary
c**********************
        ex(1)=ex(nx1)
        ey(1)=ey(nx1)
        ez(1)=ez(nx1)
c**********************
c     set E at the right hand boundary
c**********************
        ex(nx2)=ex(2)
        ey(nx2)=ey(2)
        ez(nx2)=ez(2)
c****************************
c
c      - need curl E  for it=0, only- 
c
      if(ipc .eq. 0) then
          do i=2,nx2
            deydx=ey(i) -ey(i-1)
            dezdx=ez(i) -ez(i-1)
            curley(i)=-dezdx*hxi
            curlez(i)= deydx*hxi
          enddo
      endif
c
c - - - - - - - - - -
c
      elseif(ipc .eq. 4) then
c                              use interpolation for N+1 corrector
          do i=1,nx2
            ex(i) = 0.5*( ext(i) +ex(i) ) 
            ey(i) = 0.5*( eyt(i) +ey(i) )
            ez(i) = 0.5*( ezt(i) +ez(i) )
          enddo
c
c      - boundary conditions are automatically satisfied -
c      - need curl E for this substep: - 
c
          do i=2,nx2
            deydx=ey(i) -ey(i-1)
            dezdx=ez(i) -ez(i-1)
            curley(i)=-dezdx*hxi
            curlez(i)= deydx*hxi
          enddo
c
      else
      write(*,*) ' *** incorrect ipc in ecalc *** '
      stop
      endif
c
c - - - - - - - - - -
c
      return
      end
c#######################################################################
      subroutine bcalc
c#######################################################################
c
c*********************
      include "block"
c*********************
c
       dt2=0.5 *dt
c
c calculation depends on substep within predictor-corrector sequence
c
c - - - - - - - - - -
c
      if(ipc .eq. 1) then
          do i=2,nx2
            byt(i)=by(i)-dt2*curley(i)
            bzt(i)=bz(i)-dt2*curlez(i)
          enddo
c        ************************
c            set ghost cell values
c        ************************
          byt(1)=byt(nx1)
          bzt(1)=bzt(nx1)
c
c - - - - - - - - - -
c
      elseif(ipc .eq. 2) then
          do i=2,nx2
            by(i)=byt(i)-dt2*curley(i)
            bz(i)=bzt(i)-dt2*curlez(i)
          enddo
c        ************************
c            set ghost cell values
c        ************************
          by(1)=by(nx1)
          bz(1)=bz(nx1)
c
c - - - - - - - - - -
c
      elseif(ipc .eq. 3) then
          do i=2,nx2
            by(i)=by(i)-dt2*curley(i)
            bz(i)=bz(i)-dt2*curlez(i)
          enddo
c        ************************
c            set ghost cell values
c        ************************
          by(1)=by(nx1)
          bz(1)=bz(nx1)
c
c - - - - - - - - - -
c
      elseif(ipc .eq. 4) then
          do i=2,nx2
            by(i)=byt(i)-dt2*curley(i)
            bz(i)=bzt(i)-dt2*curlez(i)
          enddo
c        ************************
c            set ghost cell values
c        ************************
          by(1)=by(nx1)
          bz(1)=bz(nx1)
c
c - - - - - - - - - -
c
      endif
      return
      end
c
c#######################################################################
      subroutine focalc
c#######################################################################
c
c*************
      include "block"
c*************
      do i=2,nx1
          dbydx=by(i+1) -by(i)
          dbzdx=bz(i+1) -bz(i)
          foy(i)= eta(i)*dbzdx*hxi
          foz(i)=-eta(i)*dbydx*hxi
      enddo
c**************************
c       first and last values (periodic) 
c*************************
        foy(1)  =foy(nx1)
        foy(nx2)=foy(2)
        foz(1)  =foz(nx1)
        foz(nx2)=foz(2)
      return
      end
c
c*************
c
c#######################################################################
      subroutine pptmp
c#######################################################################
c
c Note: this version assumes that all 'species' are in reality
c       the same and have to be treated in conjunction with each other
c*************
       INCLUDE 'block'
c*************
      do 10 k=1,nsp
      do 10 i=1,nx2
       tpar(i,k)=0.
       tperp(i,k)=0.
       trat(i,k)=0.
10    densit(i,k)=1.e-20
      isp =0
20    call getb (isp,l3,l4,nspec2,nss,nparh)
      if (isp.ne.-1) then
      do 25 l=l3,l4
      if((x(l).le.0.0).or.(x(l).ge.xmax)) go to 25
c
      rx=x(l)*hxi+1.50000001
      ix=min(nx1,int(rx))
      fx=rx-ix
      w1=(1.-fx)
      w2=fx
c replace species average with total average if the
c 'species' are really the same (approximate, only):
c (works with varying frac if incorporated in dns;
c  here, for two species, only)
      if(isp .lt. 3)   then
        dns1 = w1*dns(ix,1) +w2*dns(ix+1,1)
        dns2 = w1*dns(ix,2) +w2*dns(ix+1,2)
        dnst = dns1 +dns2
        vxavg=(  dns1 *(w1*vxs(ix,1) +w2*vxs(ix+1,1)) 
     &          +dns2 *(w1*vxs(ix,2) +w2*vxs(ix+1,2))  )/ dnst
        vyavg=(  dns1 *(w1*vys(ix,1) +w2*vys(ix+1,1)) 
     &          +dns2 *(w1*vys(ix,2) +w2*vys(ix+1,2))  )/ dnst
        vzavg=(  dns1 *(w1*vzs(ix,1) +w2*vzs(ix+1,1)) 
     &          +dns2 *(w1*vzs(ix,2) +w2*vzs(ix+1,2))  )/ dnst
c
      else
        vxavg = w1*vxs(ix,isp) +w2*vxs(ix+1,isp)  
        vyavg = w1*vys(ix,isp) +w2*vys(ix+1,isp)
        vzavg = w1*vzs(ix,isp) +w2*vzs(ix+1,isp)
      endif
      vxa=vx(l)-vxavg
      vya=vy(l)-vyavg
      vza=vz(l)-vzavg
c      for diagnostics, B is 1/2 cell off, here, for simplicity
      bxa = bxc
      bya = w1*by(ix) +w2*by(ix+1)
      bza = w1*bz(ix) +w2*bz(ix+1)
      argum=(bxa**2+bya**2+bza**2)
      if (argum .lt. 1.0e-29) argum = 1.0e-29
      btota = sqrt(argum)
      if (btota.lt.1.e-20) btota=1.e-20
      wpar=(vxa*bxa+vya*bya+vza*bza)/btota
      wperp2=vxa**2+vya**2+vza**2-wpar**2
      tpar(ix,isp)=tpar(ix,isp)+wpar**2*w1
      tpar(ix+1,isp)=tpar(ix+1,isp)+wpar**2*w2
      tperp(ix,isp)=tperp(ix,isp)+wperp2*w1
      tperp(ix+1,isp)=tperp(ix+1,isp)+wperp2*w2
      densit(ix,isp)=densit(ix,isp)+w1
      densit(ix+1,isp)=densit(ix+1,isp)+w2
c
c HERE: Add species 1 and 2 in (1)
c       - assumes same wspec for these species
c       - take care of weighting due to frac
c       That is, here the results in collected (1) are for
c       entire plasma, if there are only two species.
c
      if(isp.eq.2)  then
        rfrac = frac(2)/ frac(1)
        w1 = w1 *rfrac
        w2 = w2 *rfrac
        tpar(ix,1)=tpar(ix,1)+wpar**2*w1
        tpar(ix+1,1)=tpar(ix+1,1)+wpar**2*w2
        tperp(ix,1)=tperp(ix,1)+wperp2*w1
        tperp(ix+1,1)=tperp(ix+1,1)+wperp2*w2
        densit(ix,1)=densit(ix,1)+w1
        densit(ix+1,1)=densit(ix+1,1)+w2
      endif
25    continue
      go to 20
      endif
      do 30 k=1,nsp
c  periodic boundary condition
c***************************
      tpar(2,k)=tpar(2,k)+tpar(nx2,k)
      tperp(2,k)=tperp(2,k)+tperp(nx2,k)
      densit(2,k)=densit(2,k)+densit(nx2,k)
      tpar(nx1,k)=tpar(nx1,k)+tpar(1,k)
      tperp(nx1,k)=tperp(nx1,k)+tperp(1,k)
      densit(nx1,k)=densit(nx1,k)+densit(1,k)
      tpar(1,k)=tpar(nx1,k)
      tperp(1,k)=tperp(nx1,k)
      densit(1,k)=densit(nx1,k)
      tpar(nx2,k)=tpar(2,k)
      tperp(nx2,k)=tperp(2,k)
      densit(nx2,k)=densit(2,k)
c
30    continue
      do 60 k=1,nsp
      do 60 i=1,nx2
      tpar(i,k)=tpar(i,k)/
     &            (densit(i,k) +denmin)*wspec(k)/tx0(k)
      tperp(i,k)=tperp(i,k)/
     &            (densit(i,k) +denmin)*wspec(k)/tx0(k)/2.
      if (tpar(i,k).lt.0.00001) tpar(i,k)=0.00001
      trat(i,k)=tperp(i,k)/tpar(i,k)
60    continue
c
c     average values for energy conservation
c
c     k = 1,2: (1) includes T of all species here
      epara1 = 0.0
      eperp1 = 0.0
      epara2 = 0.0
      eperp2 = 0.0
        do i=2,nx1
          epara1 = epara1  +tpar(i,1)  *den(i)
          eperp1 = eperp1  +tperp(i,1) *den(i)
          epara2 = epara2  +tpar(i,2)  *dns(i,2)
          eperp2 = eperp2  +tperp(i,2) *dns(i,2)
        enddo
      epara1 = epara1/ float(nx)
      eperp1 = eperp1/ float(nx)
      epara2 = epara2/ float(nx)
      eperp2 = eperp2/ float(nx)
c
      write(31,*) t/wpiwci, epara1
      write(32,*) t/wpiwci, eperp1
      write(33,*) t/wpiwci, epara2
      write(34,*) t/wpiwci, eperp2
c
      write(7,*) t/wpiwci, epara1
      write(7,*) t/wpiwci, eperp1
      write(7,*) t/wpiwci, epara2
      write(7,*) t/wpiwci, eperp2
c
      return
      end
c
c#######################################################################
      subroutine getb(is,l3,l4,npart,nss,nparh)
c#######################################################################
c
      integer npart(*),nparh(*)
c****************
c writing of particle data to disk could go here.
c simply advance species in this version.
c****************
      is=is+1
      if(is.gt.nss) then
          is=-1
          return
      else
          l3=nparh(is)
          l4=npart(is)
          return
      endif
c
      end
c
c#######################################################################
      subroutine endrun (h)
c#######################################################################
c
c****************
c  writes message to output if run terminates
c****************
      character *(*) h
      write (6,10) h
      write (7,10) h
      nchar=len(h)
c
      call exit (1)
c
   10 format (///,'run terminated - - -',a10)
      end
c#######################################################################
      function random(is)
c#######################################################################
c
c convert to appropriate random generator:
c
c     random = rand()
c
      call random_number(rr)
      random = rr
c
      return
      end
