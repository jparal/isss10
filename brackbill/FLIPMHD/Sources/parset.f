      subroutine parset
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M
      use cindex_com_M
      use numpar_com_M
      use cophys_com_M
      use Timing

!
      integer :: ijkreg,kr,kx,ky,kz,      &
          istp,jstep,kstep,iref,l,n,np
      real(double) :: bxdif,bxbar,bydif,bybar,bzdif,bzbar,      &
          epsreg,xi,eta,zta,xi1,eta1,zta1,      &
          hz,dhzdz,fkx,fky,fkz,      &
          sinkx,sinky,coskx,cosky,      &
          profile,pvolume,rhobar,rhodif,rse,rsi,siebar,siedif,      &
          rnpcx,rnpcy,rnpcz,ws,wsnc,volmin,      &
          wspx,wspy,wspz,wsxi,wseta,wsnu,perturbx,      &
          c1,c2,c3,c4,c5,c6,c7,c8
      real(double) ::      &
          x1,x2,x3,x4,x5,x6,x7,x8,      &
          y1,y2,y3,y4,y5,y6,y7,y8,      &
          z1,z2,z3,z4,z5,z6,z7,z8
      real(double) ::       &
         x1278,      &
         x1476,      &
         x1573,      &
         x2385,      &
         x2684,      &
         x3456,      &
         x13 ,      &
         x16,      &
         x18,      &
         x24,      &
         x25,      &
         x27,      &
         x36,      &
         x38,      &
         x45,      &
         x47,      &
         x57,      &
         x68

!
      real(double) :: Tstart, Tfinish
      real(double) :: regvol(8)
!
      real(double) :: ijkv(8)
      integer :: ijkc(27)
      real(double) ::wghtc(27)
      real(double) :: wght1(8),wght2(8),wght3(8),wght4(8),    &
              wght5(8),wght6(8),wght7(8),wght8(8)
!
       call cpu_time(Tstart)
       pi=acos(-1.) 
!
!     set up linked list
!
      do 25 n=1,ncells
      iphead(ijkcell(n))=0
   25 continue
!
      iphead(1)=1
      do 26 np=1,npart
      pxi(np)=1.
      peta(np)=1.
      pzta(np)=1.
!
      link(np)=np+1
   26 continue
      link(npart)=0
!
!
!     set all particle variables for np=0
!
      link(0)=0
      px(0)=0.0
      py(0)=0.0
      pz(0)=0.0
      pxi(0)=0.0
      peta(0)=0.0
      pzta(0)=0.0
      mupx(0)=0.0
      mupy(0)=0.0
      mupz(0)=0.0
      up(0)=0.0
      vp(0)=0.0
      wp(0)=0.0
      ico(0)=0
      ep(0)=0.0
!
      write(8,*)'npcelx(1)', npcelx(1)
      write(8,*)'npcelx(2)', npcelx(2)
      write(8,*)'npcely(1)', npcely(1)
      write(8,*)'npcely(2)', npcely(2)
      write(8,*)'npcelz(1)', npcelz(1)
      write(8,*)'npcelz(2)', npcelz(2)
      do 1000 kr=1,nrg
!
      numtot(kr)=0
!
!     check convexity of region
!
      do l=1,8
         write(*,*) 'parset: l,kr,xvi,yvi,zvi=',l,kr,xvi(l,kr),yvi(l,kr),zvi(l,kr)
         enddo
      iref=8*(kr-1)
      istp=1
      jstep=2
      kstep=4
!
      call triple(xvi,yvi,zvi,iref+1,+1,+3,+4,regvol(1))
      call triple(xvi,yvi,zvi,iref+2,+1,-1,+4,regvol(2))
      call triple(xvi,yvi,zvi,iref+3,+1,-1,+4,regvol(3))
      call triple(xvi,yvi,zvi,iref+4,-3,-1,+4,regvol(4))
      call triple(xvi,yvi,zvi,iref+5,+1,+3,-4,regvol(5))
      call triple(xvi,yvi,zvi,iref+6,+1,-1,-4,regvol(6))
      call triple(xvi,yvi,zvi,iref+7,+1,-1,-4,regvol(7))
      call triple(xvi,yvi,zvi,iref+8,-3,-1,-4,regvol(8))
!
      regvol(5:8) = -regvol(5:8)
      volmin=0.0
      do l=1,8
      volmin=min(regvol(l),volmin)
      write(*,*) 'parset: l,regvol=',l,regvol(l)
      end do
!
      if (volmin.lt.0.0) then
      write(*,*) 'region', nrg, ' is not convex'
      go to 1000
      endif
!
      wsnc=1./(npcelx(kr)*npcely(kr)*npcelz(kr))
      rnpcx=1./real(npcelx(kr))
      rnpcy=1./real(npcely(kr))
      rnpcz=1./real(npcelz(kr))
!
!     permute list of region vertices to lexicographic order
!
      ws=xvi(3,kr)
      xvi(3,kr)=xvi(4,kr)
      xvi(4,kr)=ws
!
      ws=yvi(3,kr)
      yvi(3,kr)=yvi(4,kr)
      yvi(4,kr)=ws
!
      ws=zvi(3,kr)
      zvi(3,kr)=zvi(4,kr)
      zvi(4,kr)=ws
!
      ws=xvi(7,kr)
      xvi(7,kr)=xvi(8,kr)
      xvi(8,kr)=ws
!
      ws=yvi(7,kr)
      yvi(7,kr)=yvi(8,kr)
      yvi(8,kr)=ws
!
      ws=zvi(7,kr)
      zvi(7,kr)=zvi(8,kr)
      zvi(8,kr)=ws
!
!
      do 100 kx=1,npcelx(kr)
      do 100 ky=1,npcely(kr)
      do 100 kz=1,npcelz(kr)
!
      xi=(0.5+(kx-1))*rnpcx
      eta=(0.5+(ky-1))*rnpcy
      zta=(0.5+(kz-1))*rnpcz
!
      call weights(xi,eta,zta,wght)
!
      xi1=(kx-1)*rnpcx
      eta1=(ky-1)*rnpcy
      zta1=(kz-1)*rnpcz
      call weights(xi1,eta1,zta1,wght4)
      xi1=xi1+rnpcx
      call weights(xi1,eta1,zta1,wght1)
      eta1=eta1+rnpcy
      call weights(xi1,eta1,zta1,wght2)
      xi1=xi1-rnpcx
      call weights(xi1,eta1,zta1,wght3)
      eta1=eta1-rnpcy
      zta1=zta1+rnpcz
      call weights(xi1,eta1,zta1,wght8)
      xi1=xi1+rnpcx
      call weights(xi1,eta1,zta1,wght5)
      eta1=eta1+rnpcy
      call weights(xi1,eta1,zta1,wght6)
      xi1=xi1-rnpcx
      call weights(xi1,eta1,zta1,wght7)
!
      ijkreg=8*(kr-1)+1
      istp=1
      jstep=2
      kstep=4
      epsreg=1.e-4
!
      do 1 n=1,ncells
      np=iphead(1)
!      write(*,*)'kr=',kr
!      write(*,*)'kx=',kx
!      write(*,*)'ky=',ky
!      write(*,*)'kz=',kz
!      write(*,*)'n=',n
      if(np.eq.0) go to 5000
!
      ijk=ijkcell(n)
      ipjk=ijk+iwid
      ipjpk=ijk+iwid+jwid
      ijpk=ijk+jwid
!
      ijkp=ijk+kwid
      ipjkp=ijk+iwid+kwid
      ijpkp=ijk+jwid+kwid
      ipjpkp=ijk+iwid+jwid+kwid
!
      wspx=wght(1)*x(ipjk)     &
           +(wght(2)*x(ipjpk)     &
           +(wght(3)*x(ijpk)     &
           +(wght(4)*x(ijk)     &
           +(wght(5)*x(ipjkp)     &
           +(wght(6)*x(ipjpkp)     &
           +(wght(7)*x(ijpkp)     &
           +(wght(8)*x(ijkp))))))))
!
      wspy=wght(1)*y(ipjk)    &
           +(wght(2)*y(ipjpk)    &
           +(wght(3)*y(ijpk)    &
           +(wght(4)*y(ijk)    &
           +(wght(5)*y(ipjkp)    &
           +(wght(6)*y(ipjpkp)    &
           +(wght(7)*y(ijpkp)    &
           +(wght(8)*y(ijkp))))))))
!
      wspz=wght(1)*z(ipjk)    &
           +(wght(2)*z(ipjpk)    &
           +(wght(3)*z(ijpk)    &
           +(wght(4)*z(ijk)    &
           +(wght(5)*z(ipjkp)    &
           +(wght(6)*z(ipjpkp)    &
           +(wght(7)*z(ijpkp)    &
           +(wght(8)*z(ijkp))))))))
!
!     check whether particle is in region or not
!
      call map3d(ijkreg,istp,jstep,kstep,epsreg,     &
         xvi,yvi,zvi,     &
         wspx,wspy,wspz,wsxi,wseta,wsnu)
!
      if(wsxi*(1.-wsxi).lt.0.0     &
        .or.wseta*(1.-wseta).lt.0.0     &
        .or.wsnu*(1.-wsnu).lt.0.0)  go to 1
!
      rse=((wspx-xcenter(kr))/(rex(kr)+1.e-10))**2     &
        +((wspy-ycenter(kr))/(rey(kr)+1.e-10))**2     &
        +((wspz-zcenter(kr))/(rez(kr)+1.e-10))**2     
      rsi=((wspx-xcenter(kr))/(rix(kr)+1.e-10))**2     &
        +((wspy-ycenter(kr))/(riy(kr)+1.e-10))**2     &
        +((wspz-zcenter(kr))/(riz(kr)+1.e-10))**2
      if(rsi.lt.1..or.rse.gt.1.) go to 1
!
!     particle will be generated
!
      numtot(kr)=numtot(kr)+1
!
      iphead(1)=link(np)
      link(np)=iphead(ijk)
      iphead(ijk)=np
!
      px(np)=wspx
      py(np)=wspy
      pz(np)=wspz
!
      k=1+(ijk-1)/kwid
      j=1+(ijk-1-(k-1)*kwid)/jwid
      i=1+(ijk-1-(j-1)*jwid-(k-1)*kwid)/iwid
!
      pxi(np)=i+xi
      peta(np)=j+eta
      pzta(np)=k+zta
!
      x1=wght1(1)*x(ipjk)      &
           +(wght1(2)*x(ipjpk)      &
           +(wght1(3)*x(ijpk)      &
           +(wght1(4)*x(ijk)      &
           +(wght1(5)*x(ipjkp)      &
           +(wght1(6)*x(ipjpkp)      &
           +(wght1(7)*x(ijpkp)      &
           +(wght1(8)*x(ijkp))))))))
      x2=wght2(1)*x(ipjk)     &
           +(wght2(2)*x(ipjpk)     &
           +(wght2(3)*x(ijpk)     &
           +(wght2(4)*x(ijk)     &
           +(wght2(5)*x(ipjkp)     &
           +(wght2(6)*x(ipjpkp)     &
           +(wght2(7)*x(ijpkp)     &
           +(wght2(8)*x(ijkp))))))))
      x3=wght3(1)*x(ipjk)      &
           +(wght3(2)*x(ipjpk)      &
           +(wght3(3)*x(ijpk)      &
           +(wght3(4)*x(ijk)      &
           +(wght3(5)*x(ipjkp)      &
           +(wght3(6)*x(ipjpkp)      &
           +(wght3(7)*x(ijpkp)      &
           +(wght3(8)*x(ijkp))))))))
      x4=wght4(1)*x(ipjk)      &
           +(wght4(2)*x(ipjpk)      &
           +(wght4(3)*x(ijpk)      &
           +(wght4(4)*x(ijk)      &
           +(wght4(5)*x(ipjkp)      &
           +(wght4(6)*x(ipjpkp)      &
           +(wght4(7)*x(ijpkp)      &
           +(wght4(8)*x(ijkp))))))))
      x5=wght5(1)*x(ipjk)     &
           +(wght5(2)*x(ipjpk)     &
           +(wght5(3)*x(ijpk)     &
           +(wght5(4)*x(ijk)     &
           +(wght5(5)*x(ipjkp)     &
           +(wght5(6)*x(ipjpkp)     &
           +(wght5(7)*x(ijpkp)     &
           +(wght5(8)*x(ijkp))))))))
      x6=wght6(1)*x(ipjk)     &
           +(wght6(2)*x(ipjpk)     &
           +(wght6(3)*x(ijpk)     &
           +(wght6(4)*x(ijk)     &
           +(wght6(5)*x(ipjkp)     &
           +(wght6(6)*x(ipjpkp)     &
           +(wght6(7)*x(ijpkp)     &
           +(wght6(8)*x(ijkp))))))))
      x7=wght7(1)*x(ipjk)    &
           +(wght7(2)*x(ipjpk)    &
           +(wght7(3)*x(ijpk)    &
           +(wght7(4)*x(ijk)    &
           +(wght7(5)*x(ipjkp)    &
           +(wght7(6)*x(ipjpkp)    &
           +(wght7(7)*x(ijpkp)    &
           +(wght7(8)*x(ijkp))))))))
      x8=wght8(1)*x(ipjk)    &
           +(wght8(2)*x(ipjpk)    &
           +(wght8(3)*x(ijpk)    &
           +(wght8(4)*x(ijk)    &
           +(wght8(5)*x(ipjkp)    &
           +(wght8(6)*x(ipjpkp)    &
           +(wght8(7)*x(ijpkp)    &
           +(wght8(8)*x(ijkp))))))))
!
      y1=wght1(1)*y(ipjk)    &
           +(wght1(2)*y(ipjpk)    &
           +(wght1(3)*y(ijpk)    &
           +(wght1(4)*y(ijk)    &
           +(wght1(5)*y(ipjkp)    &
           +(wght1(6)*y(ipjpkp)    &
           +(wght1(7)*y(ijpkp)    &
           +(wght1(8)*y(ijkp))))))))
      y2=wght2(1)*y(ipjk)    &
           +(wght2(2)*y(ipjpk)    &
           +(wght2(3)*y(ijpk)    &
           +(wght2(4)*y(ijk)    &
           +(wght2(5)*y(ipjkp)    &
           +(wght2(6)*y(ipjpkp)    &
           +(wght2(7)*y(ijpkp)    &
           +(wght2(8)*y(ijkp))))))))
      y3=wght3(1)*y(ipjk)    &
           +(wght3(2)*y(ipjpk)    &
           +(wght3(3)*y(ijpk)    &
           +(wght3(4)*y(ijk)    &
           +(wght3(5)*y(ipjkp)    &
           +(wght3(6)*y(ipjpkp)    &
           +(wght3(7)*y(ijpkp)    &
           +(wght3(8)*y(ijkp))))))))
      y4=wght4(1)*y(ipjk)    &
           +(wght4(2)*y(ipjpk)    &
           +(wght4(3)*y(ijpk)    &
           +(wght4(4)*y(ijk)    &
           +(wght4(5)*y(ipjkp)    &
           +(wght4(6)*y(ipjpkp)    &
           +(wght4(7)*y(ijpkp)    &
           +(wght4(8)*y(ijkp))))))))
      y5=wght5(1)*y(ipjk)    &
           +(wght5(2)*y(ipjpk)    &
           +(wght5(3)*y(ijpk)    &
           +(wght5(4)*y(ijk)    &
           +(wght5(5)*y(ipjkp)    &
           +(wght5(6)*y(ipjpkp)    &
           +(wght5(7)*y(ijpkp)    &
           +(wght5(8)*y(ijkp))))))))
      y6=wght6(1)*y(ipjk)    &
           +(wght6(2)*y(ipjpk)    &
           +(wght6(3)*y(ijpk)    &
           +(wght6(4)*y(ijk)    &
           +(wght6(5)*y(ipjkp)    &
           +(wght6(6)*y(ipjpkp)    &
           +(wght6(7)*y(ijpkp)    &
           +(wght6(8)*y(ijkp))))))))
      y7=wght7(1)*y(ipjk)    &
           +(wght7(2)*y(ipjpk)    &
           +(wght7(3)*y(ijpk)    &
           +(wght7(4)*y(ijk)    &
           +(wght7(5)*y(ipjkp)    &
           +(wght7(6)*y(ipjpkp)    &
           +(wght7(7)*y(ijpkp)    &
           +(wght7(8)*y(ijkp))))))))
      y8=wght8(1)*y(ipjk)    &
           +(wght8(2)*y(ipjpk)    &
           +(wght8(3)*y(ijpk)    &
           +(wght8(4)*y(ijk)    &
           +(wght8(5)*y(ipjkp)    &
           +(wght8(6)*y(ipjpkp)    &
           +(wght8(7)*y(ijpkp)    &
           +(wght8(8)*y(ijkp))))))))
!
      z1=wght1(1)*z(ipjk)    &
           +(wght1(2)*z(ipjpk)    &
           +(wght1(3)*z(ijpk)    &
           +(wght1(4)*z(ijk)    &
           +(wght1(5)*z(ipjkp)    &
           +(wght1(6)*z(ipjpkp)    &
           +(wght1(7)*z(ijpkp)    &
           +(wght1(8)*z(ijkp))))))))
      z2=wght2(1)*z(ipjk)    &
           +(wght2(2)*z(ipjpk)    &
           +(wght2(3)*z(ijpk)    &
           +(wght2(4)*z(ijk)    &
           +(wght2(5)*z(ipjkp)    &
           +(wght2(6)*z(ipjpkp)    &
           +(wght2(7)*z(ijpkp)    &
           +(wght2(8)*z(ijkp))))))))
      z3=wght3(1)*z(ipjk)    &
           +(wght3(2)*z(ipjpk)    &
           +(wght3(3)*z(ijpk)    &
          +(wght3(4)*z(ijk)    &
           +(wght3(5)*z(ipjkp)    &
           +(wght3(6)*z(ipjpkp)    &
           +(wght3(7)*z(ijpkp)    &
           +(wght3(8)*z(ijkp))))))))
      z4=wght4(1)*z(ipjk)    &
           +(wght4(2)*z(ipjpk)    &
           +(wght4(3)*z(ijpk)    &
           +(wght4(4)*z(ijk)    &
           +(wght4(5)*z(ipjkp)    &
           +(wght4(6)*z(ipjpkp)    &
           +(wght4(7)*z(ijpkp)    &
           +(wght4(8)*z(ijkp))))))))
      z5=wght5(1)*z(ipjk)    &
           +(wght5(2)*z(ipjpk)    &
           +(wght5(3)*z(ijpk)    &
           +(wght5(4)*z(ijk)    &
           +(wght5(5)*z(ipjkp)    &
           +(wght5(6)*z(ipjpkp)    &
           +(wght5(7)*z(ijpkp)    &
           +(wght5(8)*z(ijkp))))))))
      z6=wght6(1)*z(ipjk)    &
           +(wght6(2)*z(ipjpk)    &
           +(wght6(3)*z(ijpk)    &
           +(wght6(4)*z(ijk)    &
           +(wght6(5)*z(ipjkp)    &
           +(wght6(6)*z(ipjpkp)    &
           +(wght6(7)*z(ijpkp)    &
           +(wght6(8)*z(ijkp))))))))
      z7=wght7(1)*z(ipjk)    &
           +(wght7(2)*z(ipjpk)    &
           +(wght7(3)*z(ijpk)    &
           +(wght7(4)*z(ijk)    &
           +(wght7(5)*z(ipjkp)    &
           +(wght7(6)*z(ipjpkp)    &
           +(wght7(7)*z(ijpkp)    &
           +(wght7(8)*z(ijkp))))))))
      z8=wght8(1)*z(ipjk)    &
           +(wght8(2)*z(ipjpk)    &
           +(wght8(3)*z(ijpk)    &
           +(wght8(4)*z(ijk)    &
           +(wght8(5)*z(ipjkp)    &
           +(wght8(6)*z(ipjpkp)    &
           +(wght8(7)*z(ijpkp)    &
           +(wght8(8)*z(ijkp))))))))
!
         x1278 = y1*z2-y2*z1+y7*z8-y8*z7
         x1476 = y1*z4-y4*z1+y7*z6-y6*z7
         x1573 = y1*z5-y5*z1+y7*z3-y3*z7
         x2385 = y2*z3-y3*z2+y8*z5-y5*z8
         x2684 = y2*z6-y6*z2+y8*z4-y4*z8
         x3456 = y3*z4-y4*z3+y5*z6-y6*z5
         x13   = y1*z3-y3*z1
         x16   = y1*z6-y6*z1
         x18   = y1*z8-y8*z1
         x24   = y2*z4-y4*z2
         x25   = y2*z5-y5*z2
         x27   = y2*z7-y7*z2
         x36   = y3*z6-y6*z3
         x38   = y3*z8-y8*z3
         x45   = y4*z5-y5*z4
         x47   = y4*z7-y7*z4
         x57   = y5*z7-y7*z5
         x68   = y6*z8-y8*z6
         c1 =(x25-x45-x24-x2385+x2684-x3456)/12.
         c2 =(x13+x36-x16+x1476-x1573-x3456)/12.
         c3 =(x24+x47-x27-x1278+x1476-x2684)/12.
         c4 =(x18-x38-x13-x1278+x1573-x2385)/12.
         c5 =(x16+x68-x18+x1278-x1476+x2684)/12.
         c6 =(x27-x57-x25+x1278-x1573+x2385)/12.
         c7 =(x38-x68-x36+x2385-x2684+x3456)/12.
         c8 =(x45+x57-x47-x1476+x1573+x3456)/12.
!
      bxdif=0.5*(bxr(2)-bxr(1))
      bxbar=0.5*(bxr(2)+bxr(1))

      bydif=0.5*(byr(2)-byr(1))
      bybar=0.5*(byr(2)+byr(1))

      bzdif=0.5*(bzr(2)-bzr(1))
      bzbar=0.5*(bzr(2)+bzr(1))

      rhodif=0.5*(rhr(2)-rhr(1))
      rhobar=0.5*(rhr(2)+rhr(1))
!
      siedif=0.5*(siep(2)-siep(1))
      siebar=0.5*(siep(2)+siep(1))
!
      profile=tanh((wspz-zcenter(kr))/0.1)
!
      pvolume=c1*x1+c2*x2+c3*x3+c4*x4     &
            +c5*x5+c6*x6+c7*x7+c8*x8
      mass(np)=pvolume*(rhodif*profile+rhobar)
      mupx(np)=pvolume*(bxdif*profile+bxbar)
      mupy(np)=pvolume*(bydif*profile+bybar)
      mupz(np)=pvolume*(bzdif*profile+bzbar)
!
      ep(np)=mass(np)*(siedif*profile+siebar)
!
      ico(np)=icoi(kr)
!
      fkx=modex(kr)*2.*pi/(float(ibar)*dx)
      sinkx=sin(fkx*wspx)
      coskx=cos(fkx*wspx)
      hz=exp(-abs(wspz-zcenter(kr))*fkx)
      if(wspz.gt.zcenter(kr)) then
      dhzdz=-fkx*hz
      else
      dhzdz=fkx*hz
      endif
      fky=modey(kr)*2.*pi/(float(jbar)*dy)
      sinky=sin(fky*wspy)
      cosky=cos(fky*wspy)
!
    perturbx=0.05
    perturby=0.05
!
      up(np)=(uvi(kr)-perturbx*coskx*dhzdz)*cosky
      wp(np)=(wvi(kr)+perturbx*sinkx*hz)*cosky
      vp(np)=vvi(kr)      !+perturbx*coskx*dhzdz
!
!
    1 continue
  100 continue
      write(*,*) 'Parset: in region', kr, numtot(kr), 'particles created'
 1000 continue
!
      call cpu_time(Tfinish)
      do l=1,20
         if(RoutineName(l).eq.'parset')   &
             CPUTime(l)=CPUTime(l)+Tfinish-Tstart
      enddo
      return
!
 5000 continue
      write(6,*) "parset: not enough particles"
      stop 'parset'
      write(*,*)'parset: particle initialization completed'
!
      return
      end subroutine parset
