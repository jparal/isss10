      subroutine resistive_diff
!
!     a routine to calculate diffusion of the magnetic field due to
!     resistance
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M
      use geometry_com_M
      use cindex_com_M
      use numpar_com_M
      use cophys_com_M
      use Scratch_com_M, ONLY : upv1, vpv1, wpv1
!
      implicit real*8 (a-h,o-z)
!
      real(double) ::  rvvol,jdotb,bsq,zero
      integer :: n
!
!     set b=0 in ghost cells
!
      call bc_field(ibp1+1,jbp1+1, kbp1+1,   &
          bxl,byl,bzl)
!
!     calculate the electric field due to resistance
!     to current flow
!     note:  curlv returns curl times vertex volume
!
      call curlv(nvtx,ijkvtx,      &
          bxl,byl,bzl,jx,jy,jz)
!
!     impose periodic boundary conditions
!
      zero=0.0
!
      call bc_current(ibp1+1,jbp1+1,kbp1+1,    &
          zero,                                &
          jx,jy,jz)
!
!     calculate the electric field
!
      do n=1,nvtx
!
         ijk=ijkvtx(n)
!
         rvvol=1./vvol(ijk)
         ex(ijk)=resistivity*jx(ijk)*rvvol
         ey(ijk)=resistivity*jy(ijk)*rvvol
         ez(ijk)=resistivity*jz(ijk)*rvvol
!
      enddo
   
!
!     calculate the curl of the electric field
!     note:  curlc divides by cell volume
!
      call curlc(ncells,ijkcell,     &
          vol,     &
          ex,ey,ez,CurlEx,CurlEy,CurlEz)
!
!
      do n=1,nvtx
      ijk=ijkvtx(n)
!
!     calculate the energy dissipated by resistive diffusion
!
      EdotJ(ijk)=(ex(ijk)*jx(ijk)      &
          +ey(ijk)*jy(ijk)     &
          +ez(ijk)*jz(ijk))*vvol(ijk)
!
      enddo
!
      if(ncyc.eq.1) JouleHeating=0.0d0
!
      do n=1,ncells
      ijk=ijkcell(n)
!
      Ohmic_heating(ijk)=0.125*(EdotJ(ijk)     &
          +EdotJ(ijk+iwid)     &
          +EdotJ(ijk+iwid+jwid)     &
          +EdotJ(ijk+jwid)     &
          +EdotJ(ijk+kwid)     &
          +EdotJ(ijk+iwid+kwid)     &
          +EdotJ(ijk+iwid+jwid+kwid)     &
          +EdotJ(ijk+jwid+kwid))     &
          *dt
!
!      Ohmic_heating(ijk)=0.5*((bxl(ijk)+bxn(ijk))*CurlEx(ijk)  &
!                            +(byl(ijk)+byn(ijk))*CurlEy(ijk)  &
!                            +(bzl(ijk)+bzn(ijk))*CurlEz(ijk))*vol(ijk)*dt
      JouleHeating=JouleHeating+Ohmic_heating(ijk)
!
      enddo
!
      return
      end subroutine resistive_diff
