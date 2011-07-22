      subroutine axisavg(ibp1,jbp1,iwid,jwid,kwid,    &
          gradxf)
!
      implicit real*8 (a-h,o-z)
!
      dimension gradxf(*)
!
!     calculate gradient at the axis by averaging gradients in the poloidal angle
!
      do 1 j=2,jbp1+1
!
      gradx=0.0
!
      do 11 i=2,ibp1
      ijk=1+(i-1)*iwid+(j-1)*jwid+kwid
      gradx=gradx+gradxf(ijk)
   11 continue
!
      do 12 i=2,ibp1+1
      ijk=1+(i-1)*iwid+(j-1)*jwid+kwid
      gradxf(ijk)=gradx
   12 continue
!
    1 continue
!
      return
      end
