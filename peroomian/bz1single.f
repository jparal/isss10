                                                                                
c                                                                               
c    Calculate ion trajectories and velocities using                            
c         ZWING magnetic field model                                            
c                                                                               
c                                                                               
                                                                                
       implicit double precision (a-h,o-z)                                      
       parameter(npn=100)                                        
       dimension dvz1(npn),dvz2(npn),dvz3(npn)                                  
       dimension dx(npn),dy(npn),dz(npn),dvx(npn),dvy(npn),dvz(npn)             
       dimension dx1(npn),dy1(npn),dz1(npn),dvx1(npn),dvy1(npn)                 
       dimension dx2(npn),dy2(npn),dz2(npn),dvx2(npn),dvy2(npn)                 
       dimension dx3(npn),dy3(npn),dz3(npn),dvx3(npn),dvy3(npn)                 
       dimension rri(npn),inpa(npn),earthdist(npn),ener(npn)                              
       dimension x(npn),y(npn),z(npn),vx(npn),vy(npn),vz(npn)                   
       dimension xp(npn),yp(npn),zp(npn),dt(npn),xnewdt(npn)                    
       dimension x1(npn),y1(npn),z1(npn),vx1(npn),vy1(npn),vz1(npn)             
       dimension x2(npn),y2(npn),z2(npn),vx2(npn),vy2(npn),vz2(npn)             
       logical swapout                                 
       double precision itime(npn)                                              
       character*80 fname1,rec
                                                                                
       common /extra1/dvz1,dvz2,dvz3,dx,dy,dz,dvx,dvy,dvz                       
       common /extra2/dx1,dy1,dz1,dvx1,dvy1,dx2,dy2,dz2,dvx2,dvy2               
       common /extra3/dx3,dy3,dz3,dvx3,dvy3,itime,rri,inpa                      
       common /extra4/x,y,z,vx,vy,vz,xp,yp,zp,dt,xnewdt                         
       common /extra5/x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2                 
                                                                                
       common/zwingpar/ th,q,alen,eps,xl,ac
                                                                                
       pi = 4.*atan(1.)

       nabc = 1
       itsteps = 1000000
       ipart=0


       open(7,file='part_input.in',status='old',err=616)
       read(7,*) ipart, xo, zo, apmass, vstream, wkevo, nskiip
       close(7)
       goto 620
616    write(6,*) 'particle input file missing'
       stop

620    continue

      open(7,file='zw_parameters.in',status='old',err=716)
      read(7,*) th,q,alen1,eps,xl,ac
      alen=alen1/q
      close(7)
      goto 720
716   write(6,*) 'field input file missing'
       stop

720    continue




c      call getarg(1,rec)
c      ipart=int(gets(rec)+0.2)

c      call getarg(2,rec)
c      xo=gets(rec)

c      call getarg(3,rec)
c      zo=gets(rec)

c      call getarg(4,rec)
c      apmass=gets(rec)

c      call getarg(5,rec)
c      vstream=gets(rec)

c      call getarg(6,rec)
c      wkevo=gets(rec)

c      call getarg(7,rec)
c      nskiip=int(gets(rec)+0.2)


       write(6,*) 'got it! '
       write(6,*) 'ipart, apmass, nskiip = ',ipart,apmass,nskiip
       write(6,*) 'xo, zo = ',xo,zo
       write(6,*) 'vstream, wkevo = ',vstream,wkevo

c
                                                                                
       dt0 = 0.005
       do 10 i=1,npn                                                            
          dt(i)=dt0                                                             
10     continue                                                                 
                                                                                
                                                                                
c------- set up calculation constants                                           
      ey = 0.1                                                                  
      re = 6371.2                                                               
      vsw0=vstream                                                                 
      xdt0=1.                                                                   
                                                                                
c----- calculate magnetic field components at launching position                
      xxx=xo                                                                    
      zzz=zo                                                                    
      call zwing(xxx,yjunk,zzz,bx,byjunk,bz)                                    
      bx0=abs(bx)                                                               
      bn0=abs(bz)                                                               
      btot=sqrt(bx0*bx0+bn0*bn0)                                                
      factor=1.0/bn0                                                        
      timec=(1.6e-19*1.e-9/1.67e-27)*bn0                                        
                                                                                
c--------------------------------------------------                             
      vswx=vsw0*bx0/btot                                                        
      vswz=vsw0*bn0/btot                                                        
      if (zo .lt. 0.0) vswz = -vswz                                             
                                                                                
      rmiqi = apmass*1.e-8                                                       
      wkev = wkevo*1.e3                                                         
      bequat = bn0*1.e-9                                                        
      blobe  = bx0*1.e-9                                                        
      vth = sqrt(2.*wkev/rmiqi)                                                 
      ebarcoefx = 1.e-3/(vth*bequat)                                            
      ebarcoefz = 1.e-3/(vth*blobe )                                            
      rhoperp = vth/(bequat/rmiqi)                                              
      vth = 1.e-3*vth                                                           
      vswx = vswx/vth                                                           
      vswz = vswz/vth                                                           
      rhoperp = 1.e-3*rhoperp/re                                                
      xbar = rhoperp                                                            
      rey = -ebarcoefx*ey                                                       
      rez = -ebarcoefz*ey                                                       
      oneoth = 1.0/th                                                           
                                                                                
c      vconvx = -rey                                                            
c      vconvz = rez                                                             
                                                                                
c     vconvx = ((-ey*bn0)/(btot*btot))*1000.                                    
c     vconvz = (-1.*(-ey)*bx0/(btot*btot))*1000.                                
                                                                                
       vconvx = 0.0                                                             
       vconvz = 0.0                                                             
                                                                                
c     vconvx = -1.0*vconvx/vth                                                  
c     vconvz = -1.0*vconvz/vth                                                  
                                                                                
c     print*, 'vconvx=',vconvx,' vconvz=',vconvz                                
c     print*, 'vswx=',vswx,' vswz=',vswz                                        
c     print*, 'bequat=',bequat,' blobe =',blobe ,' vth=',vth                    
                                                                                
                                                                                
c----------------------------------------------------                           
c   Open files for input / output
c-----------------------------------------------------                          

       write(fname1,97) ipart
97     format('particles/particle.',i5.5)
       open(31,file=fname1,status='unknown')                        
                                                                                
                                                                                
       open(47,file='./phil2.v2',status='old')                        
       do 17 jj=1,nskiip
	  read(47,*) rj,rj,rj
17     continue
       do 18 jj=1,nabc
123	  read(47,*) vx(jj),vy(jj),vz(jj)
18     continue
       close(47)
 

                                                                                
                                                                                
c-----------------------------------------------------------                    
                                                                                
       do 650 i=1,nabc                                                          
          inpa(i)=ipart                                                       
          vx(i)=vx(i)+vconvx+vswx                                               
          vy(i) = 0.
          vz(i)=vz(i)+vconvz+vswz                                               
          itime(i)=0.                                                
                                                                                
                                                                                
          x(i)=xo/xbar                                                          
          xp(i)=xo/xbar                                                         
          y(i)=0.                                                              
          yp(i) = y(i)                                                          
          zp(i)=zo/xbar                                                         
          z(i)=zo/xbar                                                          
                                                                                
c-----------------------------------------------------------                    
c      Calculating initial phases and pitchangles.                              
c-----------------------------------------------------------                    
c                                                                               
c         xxx=x(i)*xbar                                                         
c         zzz=z(i)*xbar                                                         
c         call zwing(xxx,yjunk,zzz,bx,byjunk,bz)                                
c         bx=bx*factor                                                          
c         bz=bz*factor                                                          
c         bmod=sqrt(bx*bx+bz*bz)                                                
c         ttt = bx/bmod                                             
c                                                                               
c         if (ttt .le. -1.0) ttt = -0.99999                                     
c         if (ttt .ge. 1.0)  ttt = 0.9999999                                    
c         thetab = acos(ttt)                                                    
c                                                                               
          if (bz .lt. 0.0) thetab=2.0*pi-thetab                                 
c                                                                               
c         tx=x(i)*xbar                                                          
c         ty=y(i)*xbar                                                          
c         tz=z(i)*xbar                                                          
c                                                                               
c         energy=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)                            
c         if (abs(energy) .lt. 0.00001) then                                    
c            vx(i)=-1.0                                                         
c            pphi=0.0                                                           
c            rpitch=0.0                                                         
c            ibadpart=ibadpart+1                                                
c            goto 3344                                                          
c         endif                                                                 
                                                                                
c         denom = energy*(bx*bx+bz*bz)                                          
cc          rpitch=acos((vx(i)*bx+vz(i)*bz)/sqrt(denom))                         
c         ttt=((vx(i)*bx+vz(i)*bz)/sqrt(denom))                                 
c         if (ttt .le. -1.0) ttt = -0.9999999                                   
c         if (ttt .ge. 1.0)  ttt = 0.9999999                                    
c         rpitch = acos(ttt)                                                    
                                                                                
                                                                                
c----- phase calculation ---------------------------------                      
                                                                                
c         vnorm=sin(rpitch)*sqrt(energy)                                        
cc          pphi=acos((vy(i)/abs(vnorm)))                                        
c          ttt=((vy(i)/abs(vnorm)))                                              
c          if (ttt .le. -1.0) ttt = -0.99999990                                  
c          if (ttt .ge. 1.0)  ttt = 0.9999999                                    
c          pphi = acos(ttt)                                                      
c                                                                                
c          vprody=bz*vx(i)-bx*vz(i)                                              
c          if (vprody.gt.0.0) pphi=2*pi-pphi                                     
c                                                                                
c3344      write(40,662) inpa(i)+ipadd,itime(i),tx,ty,tz,                 
c     +          vx(i),vy(i),vz(i),pphi,rpitch,energy                            
                                                                                
c-------------------------------------------------------------------            
                                                                                
                                                                                
                                                                                
          np=i                                                                  
                                                                                
          xnewdt(i)=xdt0                                                        
                                                                                
650    continue                                                                 
662   format(1x,i10,1x,f15.3,1x,10f11.5)                                        
      close(40)                                                                 
                                                                                
                                                                                
       zchange=3./xbar                                                          
c      print*, z(1),'zchange ',zchange                                          
                                                                                
       itmax=0                                                                  
       do 20 it=itmax+1,itmax+itsteps                                           
                                                                                
c        itimer=int(itime/120.)                                                 
c        np=5.+5.*float(itimer)                                                 
c        if (np .gt. nabc) np=nabc                                              
                                                                                
                                                                                
c------------------------------------------------------------------------       
c--- change dt  -------------------------------------------                     
c----------------------------------------------------------                     
       do 39 inp=1,np                                                           
         tx=x(inp)*xbar                                                         
         tz=z(inp)*xbar                                                         
         call zwing(tx,0.,tz,bx,by,bz)                                          
         btnt=sqrt(bx*bx+by*by+bz*bz)                                           
         if (btnt.eq.0.) btnt=100.
         trr=sqrt(tx*tx+tz*tz)                                                  
         if (trr .gt. 5.) then                                                  
           dt(inp)=dt0*btot/btnt                                                
         endif                                                                  
         if (trr .le. 5.  .or.  tx .lt.  2.) then                               
           dt(inp)=dt0*0.01*btot/btnt                                           
         endif                                                                  
                                                                                
         if (dt(inp) .lt.  0.000001) dt(inp)=0.000001                           
         if (dt(inp) .gt.  0.01) dt(inp)=0.01                                   
                                                                                
39     continue                                                                 
                                                                                
       do 40 inp=1,np                                                           
                                                                                
      xxx=x(inp)*xbar                                                           
      yyy=y(inp)*xbar                                                           
      zzz=z(inp)*xbar                                                           
      call zwing(xxx,yjunk,zzz,bx,byjunk,bz)                                    
      bx=bx*factor                                                              
      bxtemp=bx                                                                 
      bn=bz*factor                                                              
c      bntemp=bz                                                                
                                                                                
                                                                                
c-----------------------------------------------------------------              
                                                                                
c-------------------------------------------------------------                  
       dx(inp)=vx(inp)                                                          
       dy(inp)=vy(inp)                                                          
       dz(inp)=vz(inp)                                                          
       dvx(inp)=vy(inp)*bn                                                      
       dvy(inp)=vz(inp)*bx-vx(inp)*bn+rey                                       
       dvz(inp)=-vy(inp)*bx                                                     
c----------------------------------------------------------                     
c---- start Runge Katta --------------------------------                        
c----------------------------------------------------------                     
c---update position and vel ----------------------------------                  
c----------------------------------------------------------                     
c                                                                               
       x1(inp)=x(inp)+dt(inp)*dx(inp)*0.5                                       
       y1(inp)=y(inp)+dt(inp)*dy(inp)*0.5                                       
       z1(inp)=z(inp)+dt(inp)*dz(inp)*0.5                                       
       vx1(inp)=vx(inp)+dt(inp)*dvx(inp)*0.5                                    
       vy1(inp)=vy(inp)+dt(inp)*dvy(inp)*0.5                                    
       vz1(inp)=vz(inp)+dt(inp)*dvz(inp)*0.5                                    
                                                                                
40     continue                                                                 
                                                                                
                                                                                
                                                                                
       do 41 inp=1,np                                                           
c----------------------------------------------------------                     
c---- call derivs (x1) ----- 1 ---------------------------                      
c----------------------------------------------------------                     
cc     f1val = (x1(inp)*xbar*oneoxl)**(-q)*oneoth                               
cc     f2val = -q*f1val/(x1(inp)*xbar)                                          
cc     bx1= -ac*f1val*tanh(f1val*z1(inp)*xbar)                                  
cc     bn1=  ac*f2val*(z1(inp)*xbar*tanh(f1val*z1(inp)*xbar)-1./f1val)          
                                                                                
                                                                                
c-----------------------------------------------------------                    
c     Zwingmann model                                                          
c     Giving x & z this section of code will return bx & bz.                    
c     Assignment:   xxx=x    &    zzz=z                                         
c     Return:       bx=bx    &    bn=bz                                         
                                                                                
                                                                                
      xxx=x1(inp)*xbar                                                          
      zzz=z1(inp)*xbar                                                          
      call zwing(xxx,yjunk,zzz,bx,byjunk,bz)                                    
      bx1=bx*factor                                                             
      bn1=bz*factor                                                             
                                                                                
cc     modify on 08-08-90                                                       
c-------------------------------------------------------------                  
                                                                                
       dx1(inp)=vx1(inp)                                                        
       dy1(inp)=vy1(inp)                                                        
       dz1(inp)=vz1(inp)                                                        
       dvx1(inp)=vy1(inp)*bn1                                                   
       dvy1(inp)=vz1(inp)*bx1-vx1(inp)*bn1+rey                                  
       dvz1(inp)=-vy1(inp)*bx1                                                  
c-------------------------------------------------------------                  
c---update position and vel -- 1 -----------------------------                  
c-------------------------------------------------------------                  
                                                                                
       x1(inp)=x(inp)+dt(inp)*dx1(inp)*0.5                                      
       y1(inp)=y(inp)+dt(inp)*dy1(inp)*0.5                                      
       z1(inp)=z(inp)+dt(inp)*dz1(inp)*0.5                                      
       vx1(inp)=vx(inp)+dt(inp)*dvx1(inp)*0.5                                   
       vy1(inp)=vy(inp)+dt(inp)*dvy1(inp)*0.5                                   
       vz1(inp)=vz(inp)+dt(inp)*dvz1(inp)*0.5                                   
                                                                                
41     continue                                                                 
                                                                                
                                                                                
       do 42 inp=1,np                                                           
c---------------------------------------------------------                      
c---- call derivs (x1) --------- 2 ------------------------                     
c----------------------------------------------------------                     
cc     f1val = (x1(inp)*xbar*oneoxl)**(-q)*oneoth                               
cc     f2val = -q*f1val/(x1(inp)*xbar)                                          
cc     bx1= -ac*f1val*tanh(f1val*z1(inp)*xbar)                                  
cc     bn1=  ac*f2val*(z1(inp)*xbar*tanh(f1val*z1(inp)*xbar)-1./f1val)          
                                                                                
c-----------------------------------------------------------                    
c     Zwingmann model                                                          
c     Giving x & z this section of code will return bx & bz.                    
c     Assignment:   xxx=x    &    zzz=z                                         
c     Return:       bx=bx    &    bn=bz                                         
                                                                                
      xxx=x1(inp)*xbar                                                          
      zzz=z1(inp)*xbar                                                          
      call zwing(xxx,yjunk,zzz,bx,byjunk,bz)                                    
      bx1=bx*factor                                                             
      bn1=bz*factor                                                             
                                                                                
c      End of Tsyganenko model                                                  
c-----------------------------------------------------------------              
                                                                                
cc     modify on 08-08-90                                                       
c-------------------------------------------------------------                  
                                                                                
       dx2(inp)=vx1(inp)                                                        
       dy2(inp)=vy1(inp)                                                        
       dz2(inp)=vz1(inp)                                                        
       dvx2(inp)=vy1(inp)*bn1                                                   
       dvy2(inp)=vz1(inp)*bx1-vx1(inp)*bn1+rey                                  
       dvz2(inp)=-vy1(inp)*bx1                                                  
c-------------------------------------------------------------                  
c---update position and vel ------ 2 -------------------------                  
c-------------------------------------------------------------                  
                                                                                
       x1(inp)=x(inp)+dt(inp)*dx2(inp)                                          
       y1(inp)=y(inp)+dt(inp)*dy2(inp)                                          
       z1(inp)=z(inp)+dt(inp)*dz2(inp)                                          
       vx1(inp)=vx(inp)+dt(inp)*dvx2(inp)                                       
       vy1(inp)=vy(inp)+dt(inp)*dvy2(inp)                                       
       vz1(inp)=vz(inp)+dt(inp)*dvz2(inp)                                       
                                                                                
                                                                                
       dx2(inp)=dx1(inp)+dx2(inp)                                               
       dy2(inp)=dy1(inp)+dy2(inp)                                               
       dz2(inp)=dz1(inp)+dz2(inp)                                               
       dvx2(inp)=dvx1(inp)+dvx2(inp)                                            
       dvy2(inp)=dvy1(inp)+dvy2(inp)                                            
       dvz2(inp)=dvz1(inp)+dvz2(inp)                                            
                                                                                
42     continue                                                                 
                                                                                
                                                                                
       do 43 inp=1,np                                                           
c----------------------------------------------------------                     
c---- call derivs (x1) ------------------------------------                     
c----------------------------------------------------------                     
cc     f1val = (x1(inp)*xbar*oneoxl)**(-q)*oneoth                               
cc     f2val = -q*f1val/(x1(inp)*xbar)                                          
cc     bx1= -ac*f1val*tanh(f1val*z1(inp)*xbar)                                  
cc     bn1=  ac*f2val*(z1(inp)*xbar*tanh(f1val*z1(inp)*xbar)-1./f1val)          
                                                                                
c-----------------------------------------------------------                    
c     Zwingmann model                                                          
c     Giving x & z this section of code will return bx & bz.                    
c     Assignment:   xxx=x    &    zzz=z                                         
c     Return:       bx=bx    &    bn=bz                                         
                                                                                
      xxx=x1(inp)*xbar                                                          
      zzz=z1(inp)*xbar                                                          
      call zwing(xxx,yjunk,zzz,bx,byjunk,bz)                                    
      bx1=bx*factor                                                             
      bn1=bz*factor                                                             
                                                                                
c      End of Zwingmann model                                                  
c-----------------------------------------------------------------              
                                                                                
c-------------------------------------------------------------                  
                                                                                
       dx3(inp)=vx1(inp)                                                        
       dy3(inp)=vy1(inp)                                                        
       dz3(inp)=vz1(inp)                                                        
       dvx3(inp)=vy1(inp)*bn1                                                   
       dvy3(inp)=vz1(inp)*bx1-vx1(inp)*bn1+rey                                  
       dvz3(inp)=-vy1(inp)*bx1                                                  
c-----------------------------------------------------------------              
c-----------------------------------------------------------------              
       x(inp)=x(inp)+(dx(inp)+dx3(inp)+2.*dx2(inp))*dt(inp)/6.                  
       y(inp)=y(inp)+(dy(inp)+dy3(inp)+2.*dy2(inp))*dt(inp)/6.                  
       z(inp)=z(inp)+(dz(inp)+dz3(inp)+2.*dz2(inp))*dt(inp)/6.                  
       vx(inp)=vx(inp)+(dvx(inp)+dvx3(inp)+2.*dvx2(inp))*dt(inp)/6.             
       vy(inp)=vy(inp)+(dvy(inp)+dvy3(inp)+2.*dvy2(inp))*dt(inp)/6.             
       vz(inp)=vz(inp)+(dvz(inp)+dvz3(inp)+2.*dvz2(inp))*dt(inp)/6.             
                                                                                
c      itime(inp)=it                                                            
       itime(inp)=itime(inp)+dt(inp)/timec                                      
43     continue                                                                 
c------------------------------------------------------------                   
                                                                                
       inp=0                                                                    
101    continue                                                                 
       inp = inp+1                                                              
                                                                                
109    continue                                                                 
                                                                                
       txp= xp(inp)*xbar                                                        
       typ= yp(inp)*xbar                                                        
       tzp= zp(inp)*xbar                                                        
       tx = x(inp)*xbar                                                         
       ty = y(inp)*xbar                                                         
       tz = z(inp)*xbar                                                         
  
       ener(inp)=(vx(inp)*vx(inp)+vy(inp)*vy(inp)+
     >     vz(inp)*vz(inp))*wkevo

                                                                                
          write(31,333) inpa(inp),itime(inp),                         
     >              tx,ty,tz,vx(inp),vy(inp),vz(inp),
     >              ener(inp)                            
                                                                                
                                                                                
                                                                                
      earthdist(inp)=sqrt(tx*tx+tz*tz)                                          
      xp(inp) = x(inp)                                                          
      yp(inp) = y(inp)                                                          
      zp(inp) = z(inp)                                                          
                                                                                
                                                                                
      if (tx .le. 0.0 .or. tx .ge. 110.0 .or.                                   
     +    tz .ge. 15.0 .or. tz .le. -15.0 .or.                                  
     +    earthdist(inp) .lt. 1.5 .or.                                          
     +    ty .ge. 25.0 .or. ty .le. -25.0) then
                                                                                
                                                                                
          call swap(x(inp),x(np))                                               
          call swap(xp(inp),xp(np))                                             
          call swap(yp(inp),yp(np))                                             
          call swap(zp(inp),zp(np))                                             
          call swap(y(inp),y(np))                                               
          call swap(z(inp),z(np))                                               
          call swap(vx(inp),vx(np))                                             
          call swap(vy(inp),vy(np))                                             
          call swap(vz(inp),vz(np))                                             
          call swap(rri(inp),rri(np))                                           
          call swap(earthdist(inp),earthdist(np))                       
          call swap(itime(inp),itime(np))                                       
          call iswap(inpa(inp),inpa(np))                                        
                                                                                
          np=np - 1                                                             
          if (np .le. 0) then                                                   
             goto 22                                                            
          endif                                                                 
                                                                                
          goto 109                                                              
                                                                                
      endif                                                                     
                                                                                
      if (inp .lt. np) goto 101                                                 
                                                                                
20    continue                                                                  
22    continue                                                                  
333   format(i4,1x,f8.3,7(1x,g15.4))                                           
                                                                                
c       if (nabc .eq. 1) close(22)                                              
                                                                                
      close(32)                                                                 
      close(33)                                                                 
      close(41)                                                                 
      close(42)                                                                 
      close(36)                                                                 
                                                                                
      print*,'bz1.f terminate normally'                                         
                                                                                
      stop                                                                      
      end                                                                       
                                                                                
                                                                                
      subroutine zwing(xa,ya,za,bx,by,bz)                                       

       implicit double precision (a-h,o-z)                                      
       common/zwingpar/ th,q,alen,eps,xl,ac
                                                                                
      x=xa                                                                      
      z=za                                                                      
                                                                                
c     th=4.0                                                                    
c     q=0.5                                                                     
c     alen=100./q                                                               
c     eps=0.02                                                                  
c     xl=45.                                                                    
c     ac=75.                                                                    
                                                                                
      if (abs(z).gt.10.0) return
                                                                                
      f1=(1./th)*((x/xl)**(-q))*exp(x/alen)                                     
c     f2=-q*f1/x                                                                
      f2=(1./th)*exp(x/alen)*((-q/xl)*(x/xl)**(-q-1.)+                          
     >      (1./alen)*(x/xl)**(-q))                                             
      bz=f2*(-z*tanh(f1*z)+1./f1)                                               
      bz=-ac*bz                                                                 
      bx=f1*tanh(f1*z)                                                          
      bx=-ac*bx                                                                 
      by=0.                                                                     
                                                                                
      return                                                                    
      end                                                                       
                                                                                
                                                                                
       subroutine iswap(ia,ib)                                                  
       implicit double precision (a-h,o-z)                                      
                                                                                
       it=ia                                                                    
       ia=ib                                                                    
       ib=it                                                                    
                                                                                
       return                                                                   
       end                                                                      
                                                                                
                                                                                
       subroutine swap(a,b)                                                     
       implicit double precision (a-h,o-z)                                      
                                                                                
       t=a                                                                      
       a=b                                                                      
       b=t                                                                      
                                                                                
       return                                                                   
       end                                                                      
                                                                                
                                                                                
       subroutine lswap(la,lb)                                                  
       implicit double precision (a-h,o-z)                                      
                                                                                
       logical la,lb,lt                                                         
       lt=la                                                                    
       la=lb                                                                    
       lb=lt                                                                    
                                                                                
       return                                                                   
       end                                                                      
                                                                                
                                                                                
                                                                                
                                                                                
c------------------------------------------
      function gets(cc)
c------------------------------------------
c
c  decodes integer or floating f format
c  from character string
c
      character*(*) cc
      nn=len(cc)
c
      fak=1.
      ief=0
      l1=0
      l2=0
      do 200 i=1,nn
      if(cc(i:i).eq.'e'.or.cc(i:i).eq.'E')ief=i
      if(cc(i:i).eq.'d'.or.cc(i:i).eq.'D')ief=i
      if(l1.eq.0.and.cc(i:i).ne.' ')l1=i
      if(cc(i:i).ne.' ')l2=i
      if(cc(i:i).eq.' '.and.l2.gt.0) goto 201
200   continue
201   continue
      nn=l2
      if(ief.gt.0) then
      lex=l2-ief
      iex=-9999999
      if(lex.eq.1)read(cc(ief+1:l2),'(i1)',err=900) iex
      if(lex.eq.2)read(cc(ief+1:l2),'(i2)',err=900) iex
      if(lex.eq.3)read(cc(ief+1:l2),'(i3)',err=900) iex
      if(lex.eq.4)read(cc(ief+1:l2),'(i4)',err=900) iex
      if(lex.eq.5)read(cc(ief+1:l2),'(i5)',err=900) iex
      if(iex.gt.-999999) then
      if(iex.lt.0)fak=1./( 10.**(-iex) )
      if(iex.gt.0)fak=10.**iex
      else
      write(0,*)'gets: cannot read ',cc
      endif
      nn=ief-1
      endif
c
      sig=1.
      ss=0.
      tt=1.
      ip=0
      do 100 l=1,nn
      if(cc(l:l).ne.' ')  then
      if(cc(l:l).eq.'.') then
      ip=1
      else if(cc(l:l).eq.'-')  then
      sig=-1.
      else
c      read(cc(l:l),'(i1)',err=900) ii
      ii=ichar(cc(l:l))-48
      if(ir.lt.0.or.ir.gt.9) goto 109
      if(ip.eq.0)  then
      ss=10.*ss+float(ii)
      else
      tt=0.1*tt
      ss=ss+tt*float(ii)
      endif
      endif
      endif
100   continue
109   continue
      gets=ss*sig*fak
      return
900   continue
      write(0,*)' gets: error reading formatted integer:'
      write(0,*)nn
      write(0,'(a,a,a)')'$',cc,'$'
      stop
      end

       function lenc(line,lenmax)
c   finds length of line of text
       character *(*) line
c
       lenc=0
       do ii=1,lenmax
        iout=ichar(line(ii:ii))
        if((iout.le.32).or.(iout.gt.126)) then
         lenc=ii-1
         return
         endif
        enddo
       lenc=lenmax
       return
       end


