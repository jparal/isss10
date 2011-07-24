     subroutine divphi_3dmhd
!
     use vast_kind_param, ONLY : double
     use cindex_com_M, ONLY : ncells
     use numpar_com_M, ONLY : dt
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
       real(double) :: gradphix1,gradphiy1,gradphiz1
       real(double) :: a11, a12, a13,   &
                       a21, a22, a23,   &
                       a31, a32, a33,   &
                       b22, b23,        &
                       b32, b33,        &
                       c33,             &
                       B2, B3, C3
       integer n, ijk
!
!
       do n=1,ncells
         ijk=ijkcell(n)
         gradphix=bxn(ijk)-bmagx(ijk)
         gradphiy=byn(ijk)-bmagy(ijk)
         gradphiz=bzn(ijk)-bmagz(ijk)
!
         a11=1.+dudx(ijk)*dt
         a12=dvdx(ijk)*dt
         a13=dwdx(ijk)*dt
!
         a21=dudy(ijk)*dt
         a22=1.+dvdy(ijk)*dt
         a23=dwdy(ijk)*dt
!
         a31=dudz(ijk)*dt
         a32=dvdz(ijk)*dt
         a33=1.+dwdz(ijk)*dt
!
         b22=a22-a21*a12/a11
         b23=a23-a21*a13/a11
!
         b32=a32-a31*a12/a11
         b33=a33-a31*a13/a11
!
         c33= b33-b32*b23/b22
!
         B2=gradphiy-a21*gradphix/a11
         B3=gradphiz-a31*gradphix/a11
!
         C3=B3-b32*B2/b22
!
       gradphiz1=C3/c33
       gradphiy1=(B2-b23*gradphiz1)/b22
       gradphix1=(gradphix-a12*gradphiy1-a13*gradphiz1)/a11
!
         divuphix(ijk)=(dudx(ijk)*gradphix1   &
                       +dvdx(ijk)*gradphiy1  &
                       +dwdx(ijk)*gradphiz1   &
                        )
!
         divuphiy(ijk)=(dudy(ijk)*gradphix1   &
                       +dvdy(ijk)*gradphiy1   &
                       +dwdy(ijk)*gradphiz1   &
                       )
!
         divuphiz(ijk)=(dudz(ijk)*gradphix1   &
                       +dvdz(ijk)*gradphiy1   &
                       +dwdz(ijk)*gradphiz1   &
                       )
!
     enddo
!
     return
     end
