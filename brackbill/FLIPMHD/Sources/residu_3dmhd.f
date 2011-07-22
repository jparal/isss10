      subroutine residu_3dmhd(    &
          solution,s)
!
!     a routine to calculate the residual error
!     for 3d mhd
!
      use vast_kind_param, only:  double
      use corgan_com_M, ONLY : itdim
      use blcom_com_M, ONLY:    &
          p,bxn,byn,bzn,    &
          pl,bxl,byl,bzl,    &
          dudx,dudy,dudz,    &
          dvdx,dvdy,dvdz,    &
          dwdx,dwdy,dwdz,    &
          ijkcell,ijkvtx,    &
          vol,    &
          ul,vl,wl,    &
          CurlEx,CurlEy,CurlEz,    &
          resist,resistivity,jx,jy,jz,vvol
      use geometry_com_M
      use numpar_com_M
      use cophys_com_M
      use cindex_com_M
!
      integer :: n 
!
      real(double) ::    &
          divu, solution(itdim,*), s(itdim,*)
!
      do n=1,ncells
         ijk=ijkcell(n)
         pl(ijk)=solution(n,1)
         bxl(ijk)=solution(n,2)
         byl(ijk)=solution(n,3)
         bzl(ijk)=solution(n,4)
      enddo
!
      call stress_3dmhd
!
      call accel_3dmhd
!
      call strain_ns(ncells,ijkcell,iwid,jwid,kwid,    &
          dudx,dudy,dudz,    &
          dvdx,dvdy,dvdz,    &
          dwdx,dwdy,dwdz,    &
          ul,vl,wl,vol)
!
      do n=1,ncells
!
      ijk=ijkcell(n)
!
      divu=(dudx(ijk)+dvdy(ijk)+dwdz(ijk))
!
      s(n,1)=(pl(ijk)-p(ijk))*vol(ijk)    &
          +(gm1+1.0)*pl(ijk)*divu*vol(ijk)*dt
!
      s(n,2)=(bxl(ijk)-bxn(ijk))*vol(ijk)    &
          +(bxl(ijk)*divu     &
          -(bxl(ijk)*dudx(ijk)+byl(ijk)*dudy(ijk)+bzl(ijk)*dudz(ijk))     &
          +CurlEx(ijk))     &
!           )      &
          *vol(ijk)*dt      
!
      s(n,3)=(byl(ijk)-byn(ijk))*vol(ijk)    &
          +(byl(ijk)*divu    &
          -(bxl(ijk)*dvdx(ijk)+byl(ijk)*dvdy(ijk)+bzl(ijk)*dvdz(ijk))    &
         +CurlEy(ijk))    &
!           )    &
          *vol(ijk)*dt
!
      s(n,4)=(bzl(ijk)-bzn(ijk))*vol(ijk)    &
          +(bzl(ijk)*divu    &
          -(bxl(ijk)*dwdx(ijk)+byl(ijk)*dwdy(ijk)+bzl(ijk)*dwdz(ijk))    &
         +CurlEz(ijk))    &
!          )    &
          *vol(ijk)*dt
!
      enddo
!
      return 
      end subroutine residu_3dmhd
