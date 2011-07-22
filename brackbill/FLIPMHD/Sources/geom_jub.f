      subroutine geom(ncells,ijkcell,nvtx,ijkvtx,  &
          nvtxkm,  &
          x,y,z,  &
          vol,vvol)
!
      use vast_kind_param, ONLY : double
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z
      use geometry_com_M
      use cindex_com_M, ONLY : iwid,jwid,kwid,ibar,jbar,kbar
      use cophys_com_M, ONLY : cdlt,sdlt,dz
!
      implicit none
!
      integer, intent(in) :: ncells, ijkcell(*),nvtx, ijkvtx(*),  &
         nvtxkm
      real(double) :: dummy,   &
          x(*),y(*),z(*),  &
          vol(*),vvol(*)
      integer n,ijk,ipjk,ipjpk,ijpk, ijkp,ipjkp,ipjpkp,ijpkp,   &
          i,j,k,ijkl,ijkr,ijkb,ijkt
      real(double) ::  x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,   &
        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,    &
        x1278,y1278,z1278 ,x1476,y1476,z1476,   &
        x1573,y1573,z1573,x2385,y2385,z2385,    &
        x2684,y2684,z2684,x3456,y3456,z3456,    &
        x13,y13,z13,x16,y16,z16,x18,y18,z18,x24,y24,z24,   &
        x25,y25,z25,x27,y27,z27,x36,y36,z36,x38,y38,z38,   &
        x45,y45,z45,x47,y47,z47,x57,y57,z57,x68,y68,z68

!
!ll   calculate geometric coefficients:
      do 381 n=1,ncells
!
      ijk=ijkcell(n)
!
      ipjk=ijk+iwid
      ipjpk=ijk+iwid+jwid
      ijpk=ijk+jwid
!
      ijkp=ijk+kwid
      ipjkp=ijk+iwid+kwid
      ipjpkp=ijk+iwid+jwid+kwid
      ijpkp=ijk+jwid+kwid
!
         x1=x(ipjk)
         y1=y(ipjk)
         z1=z(ipjk)
         x2=x(ipjpk)
         y2=y(ipjpk)
         z2=z(ipjpk)
         x3=x(ijpk)
         y3=y(ijpk)
         z3=z(ijpk)
         x4=x(ijk)
         y4=y(ijk)
         z4=z(ijk)
         x5=x(ipjkp)
         y5=y(ipjkp)
         z5=z(ipjkp)
         x6=x(ipjpkp)
         y6=y(ipjpkp)
         z6=z(ipjpkp)
         x7=x(ijpkp)
         y7=y(ijpkp)
         z7=z(ijpkp)
         x8=x(ijkp)
         y8=y(ijkp)
         z8=z(ijkp)
         x1278 = y1*z2-y2*z1+y7*z8-y8*z7
         y1278 = z1*x2-z2*x1+z7*x8-z8*x7
         z1278 = x1*y2-x2*y1+x7*y8-x8*y7
         x1476 = y1*z4-y4*z1+y7*z6-y6*z7
         y1476 = z1*x4-z4*x1+z7*x6-z6*x7
         z1476 = x1*y4-x4*y1+x7*y6-x6*y7
         x1573 = y1*z5-y5*z1+y7*z3-y3*z7
         y1573 = z1*x5-z5*x1+z7*x3-z3*x7
         z1573 = x1*y5-x5*y1+x7*y3-x3*y7
         x2385 = y2*z3-y3*z2+y8*z5-y5*z8
         y2385 = z2*x3-z3*x2+z8*x5-z5*x8
         z2385 = x2*y3-x3*y2+x8*y5-x5*y8
         x2684 = y2*z6-y6*z2+y8*z4-y4*z8
         y2684 = z2*x6-z6*x2+z8*x4-z4*x8
         z2684 = x2*y6-x6*y2+x8*y4-x4*y8
         x3456 = y3*z4-y4*z3+y5*z6-y6*z5
         y3456 = z3*x4-z4*x3+z5*x6-z6*x5
         z3456 = x3*y4-x4*y3+x5*y6-x6*y5
         x13   = y1*z3-y3*z1
         y13   = z1*x3-z3*x1
         z13   = x1*y3-x3*y1
         x16   = y1*z6-y6*z1
         y16   = z1*x6-z6*x1
         z16   = x1*y6-x6*y1
         x18   = y1*z8-y8*z1
         y18   = z1*x8-z8*x1
         z18   = x1*y8-x8*y1
         x24   = y2*z4-y4*z2
         y24   = z2*x4-z4*x2
         z24   = x2*y4-x4*y2
         x25   = y2*z5-y5*z2
         y25   = z2*x5-z5*x2
         z25   = x2*y5-x5*y2
         x27   = y2*z7-y7*z2
         y27   = z2*x7-z7*x2
         z27   = x2*y7-x7*y2
         x36   = y3*z6-y6*z3
         y36   = z3*x6-z6*x3
         z36   = x3*y6-x6*y3
         x38   = y3*z8-y8*z3
         y38   = z3*x8-z8*x3
         z38   = x3*y8-x8*y3
         x45   = y4*z5-y5*z4
         y45   = z4*x5-z5*x4
         z45   = x4*y5-x5*y4
         x47   = y4*z7-y7*z4
         y47   = z4*x7-z7*x4
         z47   = x4*y7-x7*y4
         x57   = y5*z7-y7*z5
         y57   = z5*x7-z7*x5
         z57   = x5*y7-x7*y5
         x68   = y6*z8-y8*z6
         y68   = z6*x8-z8*x6
         z68   = x6*y8-x8*y6
         c1x(ijk) =(x25-x45-x24-x2385+x2684-x3456)/12.
         c1y(ijk) =(y25-y45-y24-y2385+y2684-y3456)/12.
         c1z(ijk) =(z25-z45-z24-z2385+z2684-z3456)/12.
         c2x(ijk) =(x13+x36-x16+x1476-x1573-x3456)/12.
         c2y(ijk) =(y13+y36-y16+y1476-y1573-y3456)/12.
         c2z(ijk) =(z13+z36-z16+z1476-z1573-z3456)/12.
         c3x(ijk) =(x24+x47-x27-x1278+x1476-x2684)/12.
         c3y(ijk) =(y24+y47-y27-y1278+y1476-y2684)/12.
         c3z(ijk) =(z24+z47-z27-z1278+z1476-z2684)/12.
         c4x(ijk) =(x18-x38-x13-x1278+x1573-x2385)/12.
         c4y(ijk) =(y18-y38-y13-y1278+y1573-y2385)/12.
         c4z(ijk) =(z18-z38-z13-z1278+z1573-z2385)/12.
         c5x(ijk) =(x16+x68-x18+x1278-x1476+x2684)/12.
         c5y(ijk) =(y16+y68-y18+y1278-y1476+y2684)/12.
         c5z(ijk) =(z16+z68-z18+z1278-z1476+z2684)/12.
         c6x(ijk) =(x27-x57-x25+x1278-x1573+x2385)/12.
         c6y(ijk) =(y27-y57-y25+y1278-y1573+y2385)/12.
         c6z(ijk) =(z27-z57-z25+z1278-z1573+z2385)/12.
         c7x(ijk) =(x38-x68-x36+x2385-x2684+x3456)/12.
         c7y(ijk) =(y38-y68-y36+y2385-y2684+y3456)/12.
         c7z(ijk) =(z38-z68-z36+z2385-z2684+z3456)/12.
         c8x(ijk) =(x45+x57-x47-x1476+x1573+x3456)/12.
         c8y(ijk) =(y45+y57-y47-y1476+y1573+y3456)/12.
         c8z(ijk) =(z45+z57-z47-z1476+z1573+z3456)/12.
!
!     calculate cell volume:
         vol(ijk)= c1x(ijk)*x(ipjk)+c2x(ijk)*x(ipjpk)  &
                 +c3x(ijk)*x(ijpk)+c4x(ijk)*x(ijk)  &
                 +c5x(ijk)*x(ipjkp)+c6x(ijk)*x(ipjpkp)  &
                 +c7x(ijk)*x(ijpkp)+c8x(ijk)*x(ijkp)
!
  381 continue
!
      do 4 n=1,nvtx
!
      ijk=ijkvtx(n)
!
!
!
!      vvol(ijk)=0.125*(vol(ijk)  &
!                    +vol(ijk-iwid)  &
!                    +vol(ijk-iwid-jwid)  &
!                    +vol(ijk-jwid)  &
!                    +vol(ijk-kwid)  &
!                    +vol(ijk-iwid-kwid)  &
!                    +vol(ijk-iwid-jwid-kwid)  &
!                    +vol(ijk-jwid-kwid))
!c
      vvol(ijk)=0.0
    4 continue
!
      write(*,*) 'geom: after 4'
      call volume_vtx(ncells,ijkcell,iwid,jwid,kwid,&
         x,y,z,    &
          vvol)
!
!     totlvolc=0.0
!     do 41 n=1,ncells
!     ijk=ijkcell(n)
!     totlvolc=totlvolc+vol(ijk)
!  41 continue
!
!     write(*,*) 'mshset: totlvolc=',totlvolc
!
!     totlvolv=0.0
!     do 42 n=1,nvtx
!     ijk=ijkvtx(n)
!     totlvolv=totlvolv+vvol(ijk)
!  42 continue
!
!     write(*,*) 'mshset: totlvolv=',totlvolv
!
      do 5 k=2,kbar+2
!
      if(periodic_y) then
!
      do 51 i=2,ibar+2
!
      ijkb=1+(i-1)*iwid+jwid+(k-1)*kwid
      ijkt=1+(i-1)*iwid+(jbar+1)*jwid+(k-1)*kwid
!
      vvol(ijkb)=vvol(ijkb)+vvol(ijkt)
      vvol(ijkt)=vvol(ijkb)
!
   51 continue
!
      endif
!
      if(periodic_x) then
!
      do 52 j=2,jbar+2
!
      ijkl=1+iwid+(j-1)*jwid+(k-1)*kwid
      ijkr=1+(ibar+1)*iwid+(j-1)*jwid+(k-1)*kwid
!
      vvol(ijkl)=vvol(ijkl)+vvol(ijkr)
      vvol(ijkr)=vvol(ijkl)
!
   52 continue
!
      endif
    5 continue
!
!
!     TORUSBC IS SET UP TO IMPOSE
!     PERIODIC BOUNDARY CONDTIONS ON COORDINATE
!     LIKE VECTORS
!     FOR TRULY PERIODIC VARIABLES, STRAIT SHOULD BE ZERO
!
      DUMMY=0.0
!
     write(*,*) 'geom: calling torusbc'
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c1x,c1y,c1z)
!
     write(*,*) 'geom: returning from torusbc'
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c2x,c2y,c2z)
!
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c3x,c3y,c3z)
!
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c4x,c4y,c4z)
!
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c5x,c5y,c5z)
!
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c6x,c6y,c6z)
!
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c7x,c7y,c7z)
!
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c8x,c8y,c8z)
!
      call torusbc(ibar+2,jbar+2,kbar+2,    &
          DUMMY,                            &
          c8x,c8y,vol)
!
!
      return
      end subroutine geom
