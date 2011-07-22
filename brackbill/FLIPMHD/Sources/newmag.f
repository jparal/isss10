     subroutine newmag
!
     use vast_kind_param, ONLY : double
     use cindex_com_M, ONLY : ncells
     use blcom_com_M, ONLY : ijkcell,    &
       dudx, dudy, dudz,                 &
       dvdx, dvdy, dvdz,                 &
       dwdx, dwdy, dwdz,                 &
       bxn, byn, bzn,                    &
       bmagx, bmagy, bmagz,              &
       vol
     use numpar_com_M, ONLY : dt
!
       implicit none
!
       real(double) :: divu,      &
         a11,a12,a13,   &
         a21,a22,a23,   &
         a31,a32,a33,   &
         r1,r2,r3,      &
         b22,b23,       &
         b32,b33,       &
         s2,s3,         &
         c33,t3
!
       integer n, ijk
!
       do n=1,ncells
         ijk=ijkcell(n)
!
         divu=dudx(ijk)+dvdy(ijk)+dwdz(ijk)
!
         a11=1.+(divu-dudx(ijk))*dt
         a12=-dudy(ijk)*dt
         a13=-dudz(ijk)*dt
!
         a21=-dvdx(ijk)*dt
         a22=1.+(divu-dvdy(ijk))*dt
         a23=-dvdz(ijk)*dt
!
         a31=-dwdx(ijk)*dt
         a32=-dwdy(ijk)*dt
         a33=1.+(divu-dwdz(ijk))*dt
!
         r1=bmagx(ijk)
         r2=bmagy(ijk)
         r3=bmagz(ijk)
!
         b22=a22-(a21*a12)/a11
         b23=a23-(a21*a13)/a11
!
         b32=a32-(a31*a12)/a11
         b33=a33-(a31*a13)/a11
!
         s2=r2-a21*r1/a11
         s3=r3-a31*r1/a11
!
         c33=b33-b32*b23/b22
         t3=s3-b32*s2/b22
!
         bzn(ijk)=t3/c33
         byn(ijk)=(s2-b23*bzn(ijk))/b22
         bxn(ijk)=(r1-a12*byn(ijk)-a13*bzn(ijk))/a11
!
     enddo
!
     return
     end
