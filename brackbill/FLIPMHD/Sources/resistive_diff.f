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
      allocate (upv1(itdim), vpv1(itdim), wpv1(itdim))
!     calculate the electric field due to resistance
!     to current flow
!
      call curlv(nvtx,ijkvtx,      &
          bxn,byn,bzn,jx,jy,jz)
!jub     &     bmagx,bmagy,bmagz,jx,jy,jz)
!
!     impose periodic boundary conditions
!
      zero=0.0
!
      call bc_current(ibp1+1,jbp1+1,kbp1+1,    &
          zero,                                &
          jx,jy,jz)
!
!     calculate the values of the magnetic field
!     at the vertices of the mesh, and store in CurlEx,CurlEy,CurlEz
!
      call b_vtx(ncells,ijkcell,iwid,jwid,kwid,    &
         nvtx,ijkvtx,    &
         bxl,byl,bzl,    &
         vol,upv1,    &
         bxv,byv,bzv)
!
!     calculate the electric field
!
      do n=1,nvtx
!
      ijk=ijkvtx(n)
!
      rvvol=1./vvol(ijk)
      jx(ijk)=jx(ijk)*rvvol
      jy(ijk)=jy(ijk)*rvvol
      jz(ijk)=jz(ijk)*rvvol
!
      jdotb=bxv(ijk)*jx(ijk)+byv(ijk)*jy(ijk)+bzv(ijk)*jz(ijk)
      bsq=bxv(ijk)**2+byv(ijk)**2+bzv(ijk)**2+1.d-20
!
!     the parallel current will not contribute to resistive diffusion
!
!jub      ex(ijk)=resistivity*(jx(ijk)-bxv(ijk)*jdotb/bsq)
!jub      ey(ijk)=resistivity*(jy(ijk)-byv(ijk)*jdotb/bsq)
!jub      ez(ijk)=resistivity*(jz(ijk)-bzv(ijk)*jdotb/bsq)
!
      ex(ijk)=resistivity*jx(ijk)
      ey(ijk)=resistivity*jy(ijk)
      ez(ijk)=resistivity*jz(ijk)
!
!     only the parallel current contributes to resistive diffusion
!
!jub      ex(ijk)=resistivity*bxv(ijk)*jdotb/bsq
!jub      ey(ijk)=resistivity*byv(ijk)*jdotb/bsq
!jub      ez(ijk)=resistivity*bzv(ijk)*jdotb/bsq
!
      enddo
   
!
!     calculate the curl of the electric field
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
      EdotJ(ijk)=ex(ijk)*jx(ijk)      &
          +ey(ijk)*jy(ijk)     &
          +ez(ijk)*jz(ijk)
!
      enddo
!
      JouleHeating=0.0d0
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
          *dt*vol(ijk)
!
      JouleHeating=JouleHeating+Ohmic_heating(ijk)
!
      enddo
!
      deallocate (upv1, vpv1, wpv1)
      return
      end subroutine resistive_diff
