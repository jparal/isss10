      subroutine parvol(kx,ky,kz,rnpcx,rnpcy,rnpcz,
     c                  x,y,z,ijk,iwid,jwid,kwid,pvolume)
c
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),z(*)
      dimension wght1(8),wght2(8),wght3(8),wght4(8),
     &          wght5(8),wght6(8),wght7(8),wght8(8)
c
c  compute particle volume
c
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
c    
      ipjk=ijk+iwid
      ijpk=ijk+jwid
      ipjpk=ijpk+iwid
      ijkp=ijk+kwid
      ipjpkp=ipjpk+kwid
      ipjkp=ipjk+kwid
      ijpkp=ijpk+kwid
c
      x1=wght1(1)*x(ipjk)
     &      +(wght1(2)*x(ipjpk)
     &      +(wght1(3)*x(ijpk)
     &      +(wght1(4)*x(ijk)
     &      +(wght1(5)*x(ipjkp)
     &      +(wght1(6)*x(ipjpkp)
     &      +(wght1(7)*x(ijpkp)
     &      +(wght1(8)*x(ijkp))))))))
      x2=wght2(1)*x(ipjk)
     &      +(wght2(2)*x(ipjpk)
     &      +(wght2(3)*x(ijpk)
     &      +(wght2(4)*x(ijk)
     &      +(wght2(5)*x(ipjkp)
     &      +(wght2(6)*x(ipjpkp)
     &      +(wght2(7)*x(ijpkp)
     &      +(wght2(8)*x(ijkp))))))))
      x3=wght3(1)*x(ipjk)
     &      +(wght3(2)*x(ipjpk)
     &      +(wght3(3)*x(ijpk)
     &      +(wght3(4)*x(ijk)
     &      +(wght3(5)*x(ipjkp)
     &      +(wght3(6)*x(ipjpkp)
     &      +(wght3(7)*x(ijpkp)
     &      +(wght3(8)*x(ijkp))))))))
      x4=wght4(1)*x(ipjk)
     &      +(wght4(2)*x(ipjpk)
     &      +(wght4(3)*x(ijpk)
     &      +(wght4(4)*x(ijk)
     &      +(wght4(5)*x(ipjkp)
     &      +(wght4(6)*x(ipjpkp)
     &      +(wght4(7)*x(ijpkp)
     &      +(wght4(8)*x(ijkp))))))))
      x5=wght5(1)*x(ipjk)
     &      +(wght5(2)*x(ipjpk)
     &      +(wght5(3)*x(ijpk)
     &      +(wght5(4)*x(ijk)
     &      +(wght5(5)*x(ipjkp)
     &      +(wght5(6)*x(ipjpkp)
     &      +(wght5(7)*x(ijpkp)
     &      +(wght5(8)*x(ijkp))))))))
      x6=wght6(1)*x(ipjk)
     &      +(wght6(2)*x(ipjpk)
     &      +(wght6(3)*x(ijpk)
     &      +(wght6(4)*x(ijk)
     &      +(wght6(5)*x(ipjkp)
     &      +(wght6(6)*x(ipjpkp)
     &      +(wght6(7)*x(ijpkp)
     &      +(wght6(8)*x(ijkp))))))))
      x7=wght7(1)*x(ipjk)
     &      +(wght7(2)*x(ipjpk)
     &      +(wght7(3)*x(ijpk)
     &      +(wght7(4)*x(ijk)
     &      +(wght7(5)*x(ipjkp)
     &      +(wght7(6)*x(ipjpkp)
     &      +(wght7(7)*x(ijpkp)
     &      +(wght7(8)*x(ijkp))))))))
      x8=wght8(1)*x(ipjk)
     &      +(wght8(2)*x(ipjpk)
     &      +(wght8(3)*x(ijpk)
     &      +(wght8(4)*x(ijk)
     &      +(wght8(5)*x(ipjkp)
     &      +(wght8(6)*x(ipjpkp)
     &      +(wght8(7)*x(ijpkp)
     &      +(wght8(8)*x(ijkp))))))))
c
      y1=wght1(1)*y(ipjk)
     &      +(wght1(2)*y(ipjpk)
     &      +(wght1(3)*y(ijpk)
     &      +(wght1(4)*y(ijk)
     &      +(wght1(5)*y(ipjkp)
     &      +(wght1(6)*y(ipjpkp)
     &      +(wght1(7)*y(ijpkp)
     &      +(wght1(8)*y(ijkp))))))))
      y2=wght2(1)*y(ipjk)
     &      +(wght2(2)*y(ipjpk)
     &      +(wght2(3)*y(ijpk)
     &      +(wght2(4)*y(ijk)
     &      +(wght2(5)*y(ipjkp)
     &      +(wght2(6)*y(ipjpkp)
     &      +(wght2(7)*y(ijpkp)
     &      +(wght2(8)*y(ijkp))))))))
      y3=wght3(1)*y(ipjk)
     &      +(wght3(2)*y(ipjpk)
     &      +(wght3(3)*y(ijpk)
     &      +(wght3(4)*y(ijk)
     &      +(wght3(5)*y(ipjkp)
     &      +(wght3(6)*y(ipjpkp)
     &      +(wght3(7)*y(ijpkp)
     &      +(wght3(8)*y(ijkp))))))))
      y4=wght4(1)*y(ipjk)
     &      +(wght4(2)*y(ipjpk)
     &      +(wght4(3)*y(ijpk)
     &      +(wght4(4)*y(ijk)
     &      +(wght4(5)*y(ipjkp)
     &      +(wght4(6)*y(ipjpkp)
     &      +(wght4(7)*y(ijpkp)
     &      +(wght4(8)*y(ijkp))))))))
      y5=wght5(1)*y(ipjk)
     &      +(wght5(2)*y(ipjpk)
     &      +(wght5(3)*y(ijpk)
     &      +(wght5(4)*y(ijk)
     &      +(wght5(5)*y(ipjkp)
     &      +(wght5(6)*y(ipjpkp)
     &      +(wght5(7)*y(ijpkp)
     &      +(wght5(8)*y(ijkp))))))))
      y6=wght6(1)*y(ipjk)
     &      +(wght6(2)*y(ipjpk)
     &      +(wght6(3)*y(ijpk)
     &      +(wght6(4)*y(ijk)
     &      +(wght6(5)*y(ipjkp)
     &      +(wght6(6)*y(ipjpkp)
     &      +(wght6(7)*y(ijpkp)
     &      +(wght6(8)*y(ijkp))))))))
      y7=wght7(1)*y(ipjk)
     &      +(wght7(2)*y(ipjpk)
     &      +(wght7(3)*y(ijpk)
     &      +(wght7(4)*y(ijk)
     &      +(wght7(5)*y(ipjkp)
     &      +(wght7(6)*y(ipjpkp)
     &      +(wght7(7)*y(ijpkp)
     &      +(wght7(8)*y(ijkp))))))))
      y8=wght8(1)*y(ipjk)
     &      +(wght8(2)*y(ipjpk)
     &      +(wght8(3)*y(ijpk)
     &      +(wght8(4)*y(ijk)
     &      +(wght8(5)*y(ipjkp)
     &      +(wght8(6)*y(ipjpkp)
     &      +(wght8(7)*y(ijpkp)
     &      +(wght8(8)*y(ijkp))))))))
c
c
      z1=wght1(1)*z(ipjk)
     &      +(wght1(2)*z(ipjpk)
     &      +(wght1(3)*z(ijpk)
     &      +(wght1(4)*z(ijk)
     &      +(wght1(5)*z(ipjkp)
     &      +(wght1(6)*z(ipjpkp)
     &      +(wght1(7)*z(ijpkp)
     &      +(wght1(8)*z(ijkp))))))))
      z2=wght2(1)*z(ipjk)
     &      +(wght2(2)*z(ipjpk)
     &      +(wght2(3)*z(ijpk)
     &      +(wght2(4)*z(ijk)
     &      +(wght2(5)*z(ipjkp)
     &      +(wght2(6)*z(ipjpkp)
     &      +(wght2(7)*z(ijpkp)
     &      +(wght2(8)*z(ijkp))))))))
      z3=wght3(1)*z(ipjk)
     &      +(wght3(2)*z(ipjpk)
     &      +(wght3(3)*z(ijpk)
     &      +(wght3(4)*z(ijk)
     &      +(wght3(5)*z(ipjkp)
     &      +(wght3(6)*z(ipjpkp)
     &      +(wght3(7)*z(ijpkp)
     &      +(wght3(8)*z(ijkp))))))))
      z4=wght4(1)*z(ipjk)
     &      +(wght4(2)*z(ipjpk)
     &      +(wght4(3)*z(ijpk)
     &      +(wght4(4)*z(ijk)
     &      +(wght4(5)*z(ipjkp)
     &      +(wght4(6)*z(ipjpkp)
     &      +(wght4(7)*z(ijpkp)
     &      +(wght4(8)*z(ijkp))))))))
      z5=wght5(1)*z(ipjk)
     &      +(wght5(2)*z(ipjpk)
     &      +(wght5(3)*z(ijpk)
     &      +(wght5(4)*z(ijk)
     &      +(wght5(5)*z(ipjkp)
     &      +(wght5(6)*z(ipjpkp)
     &      +(wght5(7)*z(ijpkp)
     &      +(wght5(8)*z(ijkp))))))))
      z6=wght6(1)*z(ipjk)
     &      +(wght6(2)*z(ipjpk)
     &      +(wght6(3)*z(ijpk)
     &      +(wght6(4)*z(ijk)
     &      +(wght6(5)*z(ipjkp)
     &      +(wght6(6)*z(ipjpkp)
     &      +(wght6(7)*z(ijpkp)
     &      +(wght6(8)*z(ijkp))))))))
      z7=wght7(1)*z(ipjk)
     &      +(wght7(2)*z(ipjpk)
     &      +(wght7(3)*z(ijpk)
     &      +(wght7(4)*z(ijk)
     &      +(wght7(5)*z(ipjkp)
     &      +(wght7(6)*z(ipjpkp)
     &      +(wght7(7)*z(ijpkp)
     &      +(wght7(8)*z(ijkp))))))))
      z8=wght8(1)*z(ipjk)
     &      +(wght8(2)*z(ipjpk)
     &      +(wght8(3)*z(ijpk)
     &      +(wght8(4)*z(ijk)
     &      +(wght8(5)*z(ipjkp)
     &      +(wght8(6)*z(ipjpkp)
     &      +(wght8(7)*z(ijpkp)
     &      +(wght8(8)*z(ijkp))))))))
c
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
c
      pvolume=c1*x1+c2*x2+c3*x3+c4*x4
     &       +c5*x5+c6*x6+c7*x7+c8*x8
c
      return
      end
