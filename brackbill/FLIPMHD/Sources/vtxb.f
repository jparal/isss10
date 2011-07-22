      subroutine vtxb(i1,i2,j1,j2,k1,k2,     &
          vol,bxn,byn,bzn,bxv,byv,bzv)
!
!     a routine to calculate the magnetic field at a vertex
!     by averaging cell-centered values
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      real(double) :: bxn(*),byn(*),bzn(*),     &
               bxv(*),byv(*),bzv(*),vol(*)
!
      do 1 k=k1,k2
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
!
      rvvol=1./(vol(ijk)	&
              +vol(ijk-iwid)	&
              +vol(ijk-iwid-jwid)	&
              +vol(ijk-jwid)	&
              +vol(ijk-kwid)	&
              +vol(ijk-iwid-kwid)	&
              +vol(ijk-iwid-jwid-kwid)	&
              +vol(ijk-jwid-kwid))
!
      bxv(ijk)=(bxn(ijk)*vol(ijk)	&
                    +bxn(ijk-iwid)*vol(ijk-iwid)	&
                    +bxn(ijk-iwid-jwid)*vol(ijk-iwid-jwid)	&
                    +bxn(ijk-jwid)*vol(ijk-jwid)	&
                    +bxn(ijk-kwid)*vol(ijk-kwid)	&
                    +bxn(ijk-iwid-kwid)*vol(ijk-iwid-kwid)	&
                    +bxn(ijk-iwid-jwid-kwid)*vol(ijk-iwid-jwid-kwid)	&
                    +bxn(ijk-jwid-kwid)*vol(ijk-jwid-kwid))*rvvol
!
      byv(ijk)=(byn(ijk)*vol(ijk)	&
                    +byn(ijk-iwid)*vol(ijk-iwid)	&
                    +byn(ijk-iwid-jwid)*vol(ijk-iwid-jwid)	&
                    +byn(ijk-jwid)*vol(ijk-jwid)	&
                    +byn(ijk-kwid)*vol(ijk-kwid)	&
                    +byn(ijk-iwid-kwid)*vol(ijk-iwid-kwid)	&
                    +byn(ijk-iwid-jwid-kwid)*vol(ijk-iwid-jwid-kwid)	&
                    +byn(ijk-jwid-kwid)*vol(ijk-jwid-kwid))*rvvol
!
      bzv(ijk)=(bzn(ijk)*vol(ijk)	&
                    +bzn(ijk-iwid)*vol(ijk-iwid)	&
                    +bzn(ijk-iwid-jwid)*vol(ijk-iwid-jwid)	&
                    +bzn(ijk-jwid)*vol(ijk-jwid)	&
                    +bzn(ijk-kwid)*vol(ijk-kwid)	&
                    +bzn(ijk-iwid-kwid)*vol(ijk-iwid-kwid)	&
                    +bzn(ijk-iwid-jwid-kwid)*vol(ijk-iwid-jwid-kwid)	&
                    +bzn(ijk-jwid-kwid)*vol(ijk-jwid-kwid))*rvvol
!
    1 continue
!
      return
      end subroutine vtxb
