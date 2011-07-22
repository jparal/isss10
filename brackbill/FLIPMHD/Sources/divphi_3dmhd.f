     subroutine divphi_3dmhd
!
     use vast_kind_param, ONLY : double
     use cindex_com_M, ONLY : ncells
     use blcom_com_M, ONLY : ijkcell,    &
       dudx, dudy, dudz,                 &
       dvdx, dvdy, dvdz,                 &
       dwdx, dwdy, dwdz,                 &
       bxn, byn, bzn,                    &
       bmagx, bmagy, bmagz,              &
       vol,                              &
       divuphix, divuphiy, divuphiz
!
       implicit none
!
       real(double) :: divu, gradphix, gradphiy, gradphiz
       integer n, ijk
!
       do n=1,ncells
         ijk=ijkcell(n)
         gradphix=bxn(ijk)-bmagx(ijk)
         gradphiy=byn(ijk)-bmagy(ijk)
         gradphiz=bzn(ijk)-bmagz(ijk)
!
         divu=dudx(ijk)+dvdy(ijk)+dwdz(ijk)
!
         divuphix(ijk)=(dudx(ijk)*gradphix   &
                       +dvdx(ijk)*gradphiy   &
                       +dwdx(ijk)*gradphiz   &
                       -divu*gradphix)*vol(ijk)
!
         divuphiy(ijk)=(dudy(ijk)*gradphix   &
                       +dvdy(ijk)*gradphiy   &
                       +dwdy(ijk)*gradphiz   &
                       -divu*gradphiy)*vol(ijk)
!
         divuphiz(ijk)=(dudz(ijk)*gradphix   &
                       +dvdz(ijk)*gradphiy   &
                       +dwdz(ijk)*gradphiz   &
                       -divu*gradphiz)*vol(ijk)
!
     enddo
!
     return
     end
