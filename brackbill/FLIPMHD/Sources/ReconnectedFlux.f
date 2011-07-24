      subroutine ReconnectedFlux
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
      use Scratch_com_M, ONLY : rhov_s, gradcx, gradcy, gradcz,   &
         uc, vc, wc
!
      implicit real*8 (a-h,o-z)
!
      real(double) ::  rvvol,zero, rbsq, RecFlux, OhmicHeating, Growth
      integer :: is,n
      integer, dimension(8) :: ijkv
!
      allocate(rhov_s(itdim))
      allocate(gradcx(itdim), gradcy(itdim), gradcz(itdim))
      allocate(uc(itdim), vc(itdim), wc(itdim))
!
      call vtxindx(iwid,jwid,kwid,ijkv)
!
      is=1
      do n=1,nvtx
         ijk=ijkvtx(n)
         rhov_s(ijk)=mv_s(ijk,is)/vvol(ijk)
      enddo
!
      call gradc(ncells,ijkcell,    &
          rhov_s,gradcx,gradcy,gradcz)
!
      RecFlux=0.0
      Growth=0.0
      uc=0.0
      vc=0.0
      wc=0.0
      do n=1,ncells
!
         ijk=ijkcell(n)
!
         RecFlux=RecFlux    &
            +dabs(bxn(ijk)*gradcx(ijk)   &
                 +byn(ijk)*gradcy(ijk)   &
                 +bzn(ijk)*gradcz(ijk))*vol(ijk)
       
!
         Growth=Growth+0.5*mv(ijk)*w(ijk)**2
      enddo
!
      OhmicHeating=0.0
      do n=1,nvtx
         ijk=ijkvtx(n)
         OhmicHeating=OhmicHeating    &
           +resistivity*(jx(ijk)**2+jy(ijk)**2+jz(ijk)**2)*vvol(ijk)
      enddo
!
      if(ncyc.eq.1) write(21,*) '    t    RecFlux OhmicHeating Growth'
      write(21,*)  t, RecFlux, OhmicHeating, Growth
!
       deallocate (gradcx, gradcy, gradcz)
       deallocate (uc,vc,wc)
       deallocate(rhov_s)
   
      return
      end subroutine ReconnectedFlux
