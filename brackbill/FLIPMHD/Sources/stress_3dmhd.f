      subroutine stress_3dmhd
!
!     a routine to calculate the Maxwell stress
!     and combine it with the viscous stress
!
      use vast_kind_param, only:  double
      use blcom_com_M, ONLY : pixx, pixy, pixz,  &
         piyy, piyz, pizz,                       &
         bxl, byl, bzl, p,                       &
         exx, exy, exz, eyy, eyz, ezz,           &
         ijkcell
      use numpar_com_M, ONLY : r4pi
      use cindex_com_M, ONLY : ncells
  
      integer :: n,     &
          ijk
!
!     Maxwell stress evaluated in active cells
!     ghost cell stress set too zero before calling mfnk_3dmhd
!
      do n=1,ncells
!
      ijk=ijkcell(n)
!
      exx(ijk)=-p(ijk)     &
          -0.5*(bxl(ijk)**2+byl(ijk)**2+bzl(ijk)**2)*r4pi     &
          +bxl(ijk)**2*r4pi+pixx(ijk)
      exy(ijk)=bxl(ijk)*byl(ijk)*r4pi+pixy(ijk)
      exz(ijk)=bxl(ijk)*bzl(ijk)*r4pi+pixz(ijk)
!
      eyy(ijk)=-p(ijk)      &
          -0.5*(bxl(ijk)**2+byl(ijk)**2+bzl(ijk)**2)*r4pi    &
          +byl(ijk)**2*r4pi+piyy(ijk)
      eyz(ijk)=byl(ijk)*bzl(ijk)*r4pi+piyz(ijk)      
!
      ezz(ijk)=-p(ijk)    &
          -0.5*(bxl(ijk)**2+byl(ijk)**2+bzl(ijk)**2)*r4pi    &
         +bzl(ijk)**2*r4pi+pizz(ijk)
!
      enddo
!
      return
      end subroutine stress_3dmhd
