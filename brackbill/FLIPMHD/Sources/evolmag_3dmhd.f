     subroutine evolmag_3dmhd
!
     use vast_kind_param, ONLY : double
     use cindex_com_M, ONLY : ncells
     use blcom_com_M, ONLY : ijkcell,    &
       dudx, dudy, dudz,                 &
       dvdx, dvdy, dvdz,                 &
       dwdx, dwdy, dwdz,                 &
       bxn, byn, bzn,                    &
       bxl, byl, bzl,                    &
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
!
!
         divuphix(ijk)=(bxl(ijk)*bmagx(ijk))*vol(ijk)/(bxn(ijk)+1.d-20)
!
         divuphiy(ijk)=(byl(ijk)*bmagy(ijk))*vol(ijk)/(byn(ijk)+1.d-20)
!
         divuphiz(ijk)=(bzl(ijk)*bmagz(ijk))*vol(ijk)/(bzn(ijk)+1.d-20)
!
     enddo
!
     return
     end
