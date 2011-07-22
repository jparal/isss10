      subroutine axisgrad(ibp1,jbp1,    &
          gradxf,gradyf,gradzf)
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      real(double) :: gradxf(*),gradyf(*),gradzf(*)
!
!     calculate gradient at the axis by averaging gradients in the poloidal angle
!
      do 1 j=2,jbp1+1
!
      gradx=0.0
      grady=0.0
      gradz=0.0
!
      do 11 i=2,ibp1
      ijk=1+(i-1)*iwid+(j-1)*jwid+kwid
      gradx=gradx+gradxf(ijk)
      grady=grady+gradyf(ijk)
      gradz=gradz+gradzf(ijk)
   11 continue
!
      do 12 i=2,ibp1+1
      ijk=1+(i-1)*iwid+(j-1)*jwid+kwid
      gradxf(ijk)=gradx
      gradyf(ijk)=grady
      gradzf(ijk)=gradz
   12 continue
!
    1 continue
!
      return
      end subroutine axisgrad
!     *******************************************************************
      subroutine axisvec(ibp1,jbp1,     &
          vvol,vvolaxis,vecx,vecy,vecz)
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      real(double) :: vecx(*),vecy(*),vecz(*),     &
          vvol(*),vvolaxis(*)
!
!     calculate average vector at the axis by volume averaging
!
      do 1 j=2,jbp1+1
!
      avgx=0.0
      avgy=0.0
      avgz=0.0
!
      do 11 i=2,ibp1
      ijk=1+(i-1)*iwid+(j-1)*jwid+kwid
      avgx=avgx+vecx(ijk)*vvolaxis(ijk)
      avgy=avgy+vecy(ijk)*vvolaxis(ijk)
      avgz=avgz+vecz(ijk)*vvolaxis(ijk)
   11 continue
!
      do 12 i=2,ibp1+1
      ijk=1+(i-1)*iwid+(j-1)*jwid+kwid
      rvvol=1./vvol(ijk)
      vecx(ijk)=avgx*rvvol
      vecy(ijk)=avgy*rvvol
      vecz(ijk)=avgz*rvvol
   12 continue
!
    1 continue
!
      return
      end subroutine axisvec
!     ****************************************************************
      subroutine axisvol(ibp1,jbp1,    & 
          VVOL,NVTX,IJKVTX,     &
          vol,vvolaxis)
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      real(double) :: vol(*),vvolaxis(*),     &
          vvol(*)
      integer :: ijkvtx(*)
!
      do 11 n=1,nvtx
      vvolaxis(ijkvtx(n))=vvol(ijkvtx(n))
   11 continue
!
!
      do 1 j=2,jbp1+1
      do 1 i=2,ibp1+1
!
      ijk=1+(i-1)*iwid+(j-1)*jwid+kwid
!
      vvolaxis(ijk)=0.125*(vol(ijk)      &
                         +vol(ijk-iwid)      &
                         +vol(ijk-iwid-jwid)      &
                         +vol(ijk-jwid))
!
    1 continue
!
      return
      end subroutine axisvol
!     *********************************************************************
      subroutine bcphi(ibp1,jbp1,kbp1,phi)
!
!     a routine to impose doubly periodic boundary conditions on phi
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      real(double) :: phi(*)
!
      do 5 k=2,kbp1+1
!
      do 51 j=1,jbp1+1
      ijkl=1+(j-1)*jwid+(k-1)*kwid
      ijkr=1+ibp1*iwid+(j-1)*jwid+(k-1)*kwid
      phi(ijkr)=phi(ijkl+iwid)
      phi(ijkl)=phi(ijkr-iwid)
   51 continue
!
      do 52 i=1,ibp1+1
      ijkb=1+(i-1)*iwid+(k-1)*kwid
      ijkt=1+(i-1)*iwid+jbp1*jwid+(k-1)*kwid
      phi(ijkb)=phi(ijkt-jwid)
      phi(ijkt)=phi(ijkb+jwid)
   52 continue
!
    5 continue
!
!
!     k=kbp2 plane
!
      do 1 j=2,jbp1
      do 1 i=2,ibp1
      ijk=1+(i-1)*iwid+(j-1)*jwid+kbp1*kwid
      phi(ijk)=phi(ijk-kwid)
    1 continue
!
!     minor axis of torus
!
      do 4 j=2,jbp1
      do 42 i=2,ibp1
      ijk=1+(i-1)*iwid+(j-1)*jwid
      phi(ijk)=0.0
   42 continue
!
    4 continue
!
      return
      end subroutine bcphi
!     ******************************************************************

      subroutine curlc(ncells,ijkcell,     &
          vol,     &
          vx,vy,vz,curlcx,curlcy,curlcz)
!
!     a routine to calculate a cell-centered curl from a vertex-centered
!     vector with components vx,vy,vz
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use geometry_com_M
      implicit none
!
      integer :: ncells,ijkcell(*),n,ijk
!
      real(double) ::    &
          vol(*),    &
          vx(*),vy(*),vz(*),    &
          curlcx(*),curlcy(*),curlcz(*)
      real(double) :: rvol
!
      do 100 n=1,ncells
!
      ijk=ijkcell(n)
!
      rvol=1./vol(ijk)
!
      curlcx(ijk)=(c1y(ijk)*vz(ijk+iwid)      &
                +c2y(ijk)*vz(ijk+iwid+jwid)      &
                +c3y(ijk)*vz(ijk+jwid)      &
                +c4y(ijk)*vz(ijk)      &
                +c5y(ijk)*vz(ijk+iwid+kwid)      &
                +c6y(ijk)*vz(ijk+iwid+jwid+kwid)      &
                +c7y(ijk)*vz(ijk+jwid+kwid)      &
                +c8y(ijk)*vz(ijk+kwid)      &
                -c1z(ijk)*vy(ijk+iwid)      &
                -c2z(ijk)*vy(ijk+iwid+jwid)      &
                -c3z(ijk)*vy(ijk+jwid)      &
                -c4z(ijk)*vy(ijk)      &
                -c5z(ijk)*vy(ijk+iwid+kwid)      &
                -c6z(ijk)*vy(ijk+iwid+jwid+kwid)      &
                -c7z(ijk)*vy(ijk+jwid+kwid)      &
                -c8z(ijk)*vy(ijk+kwid))      &
                  *rvol
!
      curlcy(ijk)=(c1z(ijk)*vx(ijk+iwid)     &
                +c2z(ijk)*vx(ijk+iwid+jwid)     &
                +c3z(ijk)*vx(ijk+jwid)     &
                +c4z(ijk)*vx(ijk)     &
                +c5z(ijk)*vx(ijk+iwid+kwid)     &
                +c6z(ijk)*vx(ijk+iwid+jwid+kwid)     &
                +c7z(ijk)*vx(ijk+jwid+kwid)     &
                +c8z(ijk)*vx(ijk+kwid)     &
                -c1x(ijk)*vz(ijk+iwid)     &
                -c2x(ijk)*vz(ijk+iwid+jwid)     &
                -c3x(ijk)*vz(ijk+jwid)     &
                -c4x(ijk)*vz(ijk)     &
                -c5x(ijk)*vz(ijk+iwid+kwid)     &
                -c6x(ijk)*vz(ijk+iwid+jwid+kwid)     &
                -c7x(ijk)*vz(ijk+jwid+kwid)     &
                -c8x(ijk)*vz(ijk+kwid))     &
                  *rvol
!
      curlcz(ijk)=(c1x(ijk)*vy(ijk+iwid)       &
                +c2x(ijk)*vy(ijk+iwid+jwid)       &
                +c3x(ijk)*vy(ijk+jwid)       &
                +c4x(ijk)*vy(ijk)       &
                +c5x(ijk)*vy(ijk+iwid+kwid)       &
                +c6x(ijk)*vy(ijk+iwid+jwid+kwid)       &
                +c7x(ijk)*vy(ijk+jwid+kwid)       &
                +c8x(ijk)*vy(ijk+kwid)       &
                -c1y(ijk)*vx(ijk+iwid)       &
                -c2y(ijk)*vx(ijk+iwid+jwid)       &
                -c3y(ijk)*vx(ijk+jwid)       &
                -c4y(ijk)*vx(ijk)       &
                -c5y(ijk)*vx(ijk+iwid+kwid)       &
                -c6y(ijk)*vx(ijk+iwid+jwid+kwid)       &
                -c7y(ijk)*vx(ijk+jwid+kwid)       &
                -c8y(ijk)*vx(ijk+kwid))       &
                  *rvol
!
  100 continue
!
      return
      end subroutine curlc
!     ************************************************************************
      subroutine curlv(nvtx,ijkvtx,   &
          vx,vy,vz,curlvx,curlvy,curlvz)
!
!     a routine to calculate a vertex-centered curl
!     from a cell-centered vector
!     vector with components vx,vy,vz
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use geometry_com_M
      implicit real*8 (a-h,o-z)
!
      integer ::      &
     &     ijkvtx(*)
      real(double) ::    &
          vx(*),vy(*),vz(*),     &
          curlvx(*),curlvy(*),curlvz(*)
!
      do 100 n=1,nvtx
!
      ijk=ijkvtx(n)
!
!
      curlvx(ijk)=-(c1y(ijk-iwid)*vz(ijk-iwid)     &
                +c2y(ijk-iwid-jwid)*vz(ijk-iwid-jwid)     &
                +c3y(ijk-jwid)*vz(ijk-jwid)     &
                +c4y(ijk)*vz(ijk)      &
                +c5y(ijk-iwid-kwid)*vz(ijk-iwid-kwid)      &
                +c6y(ijk-iwid-jwid-kwid)*vz(ijk-iwid-jwid-kwid)      &
                +c7y(ijk-jwid-kwid)*vz(ijk-jwid-kwid)      &
                +c8y(ijk-kwid)*vz(ijk-kwid)      &
                -c1z(ijk-iwid)*vy(ijk-iwid)      &
                -c2z(ijk-iwid-jwid)*vy(ijk-iwid-jwid)      &
                -c3z(ijk-jwid)*vy(ijk-jwid)      &
                -c4z(ijk)*vy(ijk)      &
                -c5z(ijk-iwid-kwid)*vy(ijk-iwid-kwid)      &
                -c6z(ijk-iwid-jwid-kwid)*vy(ijk-iwid-jwid-kwid)      &
                -c7z(ijk-jwid-kwid)*vy(ijk-jwid-kwid)      &
                -c8z(ijk-kwid)*vy(ijk-kwid))      
!
      curlvy(ijk)=-(c1z(ijk-iwid)*vx(ijk-iwid)    &
                +c2z(ijk-iwid-jwid)*vx(ijk-iwid-jwid)    &
                +c3z(ijk-jwid)*vx(ijk-jwid)    &
                +c4z(ijk)*vx(ijk)    &
                +c5z(ijk-iwid-kwid)*vx(ijk-iwid-kwid)    &
                +c6z(ijk-iwid-jwid-kwid)*vx(ijk-iwid-jwid-kwid)    &
                +c7z(ijk-jwid-kwid)*vx(ijk-jwid-kwid)    &
                +c8z(ijk-kwid)*vx(ijk-kwid)    &
                -c1x(ijk-iwid)*vz(ijk-iwid)    &
                -c2x(ijk-iwid-jwid)*vz(ijk-iwid-jwid)    &
                -c3x(ijk-jwid)*vz(ijk-jwid)    &
                -c4x(ijk)*vz(ijk)    &
                -c5x(ijk-iwid-kwid)*vz(ijk-iwid-kwid)    &
                -c6x(ijk-iwid-jwid-kwid)*vz(ijk-iwid-jwid-kwid)    &
                -c7x(ijk-jwid-kwid)*vz(ijk-jwid-kwid)    &
                -c8x(ijk-kwid)*vz(ijk-kwid))
!
      curlvz(ijk)=-(c1x(ijk-iwid)*vy(ijk-iwid)      &
                +c2x(ijk-iwid-jwid)*vy(ijk-iwid-jwid)      &
                +c3x(ijk-jwid)*vy(ijk-jwid)      &
                +c4x(ijk)*vy(ijk)      &
                +c5x(ijk-iwid-kwid)*vy(ijk-iwid-kwid)      &
                +c6x(ijk-iwid-jwid-kwid)*vy(ijk-iwid-jwid-kwid)      &
                +c7x(ijk-jwid-kwid)*vy(ijk-jwid-kwid)      &
                +c8x(ijk-kwid)*vy(ijk-kwid)      &
                -c1y(ijk-iwid)*vx(ijk-iwid)      &
                -c2y(ijk-iwid-jwid)*vx(ijk-iwid-jwid)      &
                -c3y(ijk-jwid)*vx(ijk-jwid)      &
                -c4y(ijk)*vx(ijk)      &
                -c5y(ijk-iwid-kwid)*vx(ijk-iwid-kwid)      &
                -c6y(ijk-iwid-jwid-kwid)*vx(ijk-iwid-jwid-kwid)      &
                -c7y(ijk-jwid-kwid)*vx(ijk-jwid-kwid)      &
                -c8y(ijk-kwid)*vx(ijk-kwid))
!
  100 continue
!
      return
      end subroutine curlv
!     ***********************************************************************

      subroutine divc(ncells,ijkcell,    &
                        vol,ex,ey,ez,dive)
!
!     a routine to calculate the divergence of a vector e
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use geometry_com_M
      implicit real*8 (a-h,o-z)
!
      integer :: ijkcell(*)
      real(double) ::     &
          vol(*),ex(*),ey(*),ez(*),dive(*)
!
      do 1 n=1,ncells
!
      ijk=ijkcell(n)
!
      divex=(c1x(ijk)*ex(ijk+iwid)     &
          +(c2x(ijk)*ex(ijk+iwid+jwid)     &
          +(c3x(ijk)*ex(ijk+jwid)     &
          +(c4x(ijk)*ex(ijk)     &
          +(c5x(ijk)*ex(ijk+iwid+kwid)     &
          +(c6x(ijk)*ex(ijk+iwid+jwid+kwid)     &
          +(c7x(ijk)*ex(ijk+jwid+kwid)     &
          +(c8x(ijk)*ex(ijk+kwid)))))))))
!
      divey=(c1y(ijk)*ey(ijk+iwid)     &
          +(c2y(ijk)*ey(ijk+iwid+jwid)     &
          +(c3y(ijk)*ey(ijk+jwid)     &
          +(c4y(ijk)*ey(ijk)     &
          +(c5y(ijk)*ey(ijk+iwid+kwid)     &
          +(c6y(ijk)*ey(ijk+iwid+jwid+kwid)     &
          +(c7y(ijk)*ey(ijk+jwid+kwid)     &
          +(c8y(ijk)*ey(ijk+kwid)))))))))
!
      divez=(c1z(ijk)*ez(ijk+iwid)     &
          +(c2z(ijk)*ez(ijk+iwid+jwid)     &
          +(c3z(ijk)*ez(ijk+jwid)     &
          +(c4z(ijk)*ez(ijk)     &
          +(c5z(ijk)*ez(ijk+iwid+kwid)     &
          +(c6z(ijk)*ez(ijk+iwid+jwid+kwid)     &
          +(c7z(ijk)*ez(ijk+jwid+kwid)     &
          +(c8z(ijk)*ez(ijk+kwid)))))))))
!
      dive(ijk)=(divex+divey+divez)/vol(ijk)
!
    1 continue
!
      return
      end subroutine divc
!     *******************************************************************
      subroutine divv(     &
          vx,vy,vz,diverge)
!
!     a routine to calculate a vertex-centered divergence
!     from a cell-centered vector
!     vector with components vx,vy,vz
!
      use vast_kind_param, only:  double
      use corgan_com_M, ONLY : itdim
      use geometry_com_M
      use cindex_com_M, ONLY : nvtx, iwid, jwid, kwid
      use blcom_com_M, ONLY : ijkvtx
      implicit none
!
      real(double) ::                        &
          divergex, divergey, divergez,      &
          vx(*),vy(*),vz(*),                 &
          diverge(*)
      integer :: n, ijk
!
      do 100 n=1,nvtx
!
      ijk=ijkvtx(n)
!
!
!
      divergex=-(c1x(ijk-iwid)*vx(ijk-iwid)     &
                +c2x(ijk-iwid-jwid)*vx(ijk-iwid-jwid)     &
                +c3x(ijk-jwid)*vx(ijk-jwid)     &
                +c4x(ijk)*vx(ijk)     &
                +c5x(ijk-iwid-kwid)*vx(ijk-iwid-kwid)     &
                +c6x(ijk-iwid-jwid-kwid)*vx(ijk-iwid-jwid-kwid)     &
                +c7x(ijk-jwid-kwid)*vx(ijk-jwid-kwid)     &
                +c8x(ijk-kwid)*vx(ijk-kwid))
!
!
      divergey=-(c1y(ijk-iwid)*vy(ijk-iwid)     &
                +c2y(ijk-iwid-jwid)*vy(ijk-iwid-jwid)     &
                +c3y(ijk-jwid)*vy(ijk-jwid)     &
                +c4y(ijk)*vy(ijk)     &
                +c5y(ijk-iwid-kwid)*vy(ijk-iwid-kwid)     &
                +c6y(ijk-iwid-jwid-kwid)*vy(ijk-iwid-jwid-kwid)     &
                +c7y(ijk-jwid-kwid)*vy(ijk-jwid-kwid)     &
                +c8y(ijk-kwid)*vy(ijk-kwid))

      divergez=-(c1z(ijk-iwid)*vz(ijk-iwid)     &
                +c2z(ijk-iwid-jwid)*vz(ijk-iwid-jwid)     &
                +c3z(ijk-jwid)*vz(ijk-jwid)     &
                +c4z(ijk)*vz(ijk)     &
               +c5z(ijk-iwid-kwid)*vz(ijk-iwid-kwid)     &
                +c6z(ijk-iwid-jwid-kwid)*vz(ijk-iwid-jwid-kwid)     &
                +c7z(ijk-jwid-kwid)*vz(ijk-jwid-kwid)     &
                +c8z(ijk-kwid)*vz(ijk-kwid))
!
      diverge(ijk)=divergex+divergey+divergez

!
  100 continue
!
      return
      end subroutine divv
!     **********************************************************************
      subroutine divpi(nvtx,ijkvtx,     &
          pixx,pixy,pixz,piyy,piyz,pizz,     &
          divpix,divpiy,divpiz)
!
!     a rouTINE to calculate the divergence of a tensor
!     the tensor is cell-centered
!     the divergence, a vector, is vertex-centered
!
      use vast_kind_param, ONLY : double
      use geometry_com_M
      use cindex_com_M, ONLY : iwid,jwid,kwid
      implicit real*8 (a-h,o-z)
!
      integer ::   &
          ijkvtx(*)
      real(double) ::      &
          pixx(*),pixy(*),pixz(*),piyy(*),piyz(*),pizz(*),     &
          divpix(*),divpiy(*),divpiz(*)
!
      do 100 n=1,nvtx
!
      ijk=ijkvtx(n)
!
!
      divxx=-(c1x(ijk-iwid)*pixx(ijk-iwid)     &
                +(c2x(ijk-iwid-jwid)*pixx(ijk-iwid-jwid)     &
                +(c3x(ijk-jwid)*pixx(ijk-jwid)     &
                +(c4x(ijk)*pixx(ijk)     &
                +(c5x(ijk-iwid-kwid)*pixx(ijk-iwid-kwid)     &
                +(c6x(ijk-iwid-jwid-kwid)*pixx(ijk-iwid-jwid-kwid)     &
                +(c7x(ijk-jwid-kwid)*pixx(ijk-jwid-kwid)     &
                +(c8x(ijk-kwid)*pixx(ijk-kwid)))))))))
!
      divyx=-(c1y(ijk-iwid)*pixy(ijk-iwid)    &
                +(c2y(ijk-iwid-jwid)*pixy(ijk-iwid-jwid)    &
                +(c3y(ijk-jwid)*pixy(ijk-jwid)    &
                +(c4y(ijk)*pixy(ijk)    &
                +(c5y(ijk-iwid-kwid)*pixy(ijk-iwid-kwid)    &
                +(c6y(ijk-iwid-jwid-kwid)*pixy(ijk-iwid-jwid-kwid)    &
                +(c7y(ijk-jwid-kwid)*pixy(ijk-jwid-kwid)    &
                +(c8y(ijk-kwid)*pixy(ijk-kwid)))))))))
!
      divzx=-(c1z(ijk-iwid)*pixz(ijk-iwid)     &
                +(c2z(ijk-iwid-jwid)*pixz(ijk-iwid-jwid)     &
                +(c3z(ijk-jwid)*pixz(ijk-jwid)     &
                +(c4z(ijk)*pixz(ijk)     &
                +(c5z(ijk-iwid-kwid)*pixz(ijk-iwid-kwid)     &
                +(c6z(ijk-iwid-jwid-kwid)*pixz(ijk-iwid-jwid-kwid)     &
                +(c7z(ijk-jwid-kwid)*pixz(ijk-jwid-kwid)     &
                +(c8z(ijk-kwid)*pixz(ijk-kwid)))))))))
!
      divpix(ijk)=divxx+divyx+divzx
!
  100 continue
!
      do 200 n=1,nvtx
!
      ijk=ijkvtx(n)
!
      divxy=-(c1x(ijk-iwid)*pixy(ijk-iwid)     &
                +(c2x(ijk-iwid-jwid)*pixy(ijk-iwid-jwid)     &
                +(c3x(ijk-jwid)*pixy(ijk-jwid)     &
                +(c4x(ijk)*pixy(ijk)     &
                +(c5x(ijk-iwid-kwid)*pixy(ijk-iwid-kwid)     &
                +(c6x(ijk-iwid-jwid-kwid)*pixy(ijk-iwid-jwid-kwid)     &
                +(c7x(ijk-jwid-kwid)*pixy(ijk-jwid-kwid)     &
                +(c8x(ijk-kwid)*pixy(ijk-kwid)))))))))
!
      divyy=-(c1y(ijk-iwid)*piyy(ijk-iwid)    &
                +(c2y(ijk-iwid-jwid)*piyy(ijk-iwid-jwid)    &
                +(c3y(ijk-jwid)*piyy(ijk-jwid)    &
                +(c4y(ijk)*piyy(ijk)    &
                +(c5y(ijk-iwid-kwid)*piyy(ijk-iwid-kwid)    &
                +(c6y(ijk-iwid-jwid-kwid)*piyy(ijk-iwid-jwid-kwid)    &
                +(c7y(ijk-jwid-kwid)*piyy(ijk-jwid-kwid)    &
                +(c8y(ijk-kwid)*piyy(ijk-kwid)))))))))
!
      divzy=-(c1z(ijk-iwid)*piyz(ijk-iwid)     &
                +(c2z(ijk-iwid-jwid)*piyz(ijk-iwid-jwid)     &
                +(c3z(ijk-jwid)*piyz(ijk-jwid)     &
                +(c4z(ijk)*piyz(ijk)     &
                +(c5z(ijk-iwid-kwid)*piyz(ijk-iwid-kwid)     &
                +(c6z(ijk-iwid-jwid-kwid)*piyz(ijk-iwid-jwid-kwid)     &
                +(c7z(ijk-jwid-kwid)*piyz(ijk-jwid-kwid)     &
                +(c8z(ijk-kwid)*piyz(ijk-kwid)))))))))
!
      divpiy(ijk)=divxy+divyy+divzy
!
  200 continue
!
      do 300 n=1,nvtx
!
      ijk=ijkvtx(n)
!
      divxz=-(c1x(ijk-iwid)*pixz(ijk-iwid)      &
                +(c2x(ijk-iwid-jwid)*pixz(ijk-iwid-jwid)      &
                +(c3x(ijk-jwid)*pixz(ijk-jwid)      &
                +(c4x(ijk)*pixz(ijk)      &
                +(c5x(ijk-iwid-kwid)*pixz(ijk-iwid-kwid)      &
                +(c6x(ijk-iwid-jwid-kwid)*pixz(ijk-iwid-jwid-kwid)      &
                +(c7x(ijk-jwid-kwid)*pixz(ijk-jwid-kwid)      &
                +(c8x(ijk-kwid)*pixz(ijk-kwid)))))))))
!
      divyz=-(c1y(ijk-iwid)*piyz(ijk-iwid)     &
                +(c2y(ijk-iwid-jwid)*piyz(ijk-iwid-jwid)     &
                +(c3y(ijk-jwid)*piyz(ijk-jwid)     &
                +(c4y(ijk)*piyz(ijk)     &
                +(c5y(ijk-iwid-kwid)*piyz(ijk-iwid-kwid)     &
                +(c6y(ijk-iwid-jwid-kwid)*piyz(ijk-iwid-jwid-kwid)     &
                +(c7y(ijk-jwid-kwid)*piyz(ijk-jwid-kwid)     &
                +(c8y(ijk-kwid)*piyz(ijk-kwid)))))))))
!
      divzz=-(c1z(ijk-iwid)*pizz(ijk-iwid)     &
                +(c2z(ijk-iwid-jwid)*pizz(ijk-iwid-jwid)     &
                +(c3z(ijk-jwid)*pizz(ijk-jwid)     &
                +(c4z(ijk)*pizz(ijk)     &
                +(c5z(ijk-iwid-kwid)*pizz(ijk-iwid-kwid)     &
                +(c6z(ijk-iwid-jwid-kwid)*pizz(ijk-iwid-jwid-kwid)     &
                +(c7z(ijk-jwid-kwid)*pizz(ijk-jwid-kwid)     &
                +(c8z(ijk-kwid)*pizz(ijk-kwid)))))))))
!
      divpiz(ijk)=divxz+divyz+divzz
!
  300 continue
!
      return
      end subroutine divpi
!     **********************************************************************
      subroutine gradc(ncells,ijkcell,    &
          psi,gradcx,gradcy,gradcz)
!
!     a routine to calculate a cell-centered grad from a vertex-centered
!     scalar
!
      use vast_kind_param, ONLY : double
      use blcom_com_M, ONLY : vol
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use geometry_com_M
      implicit none
!
      integer :: ncells,ijkcell(*)
      real(double) ::    &
          psi(*),    &
          gradcx(*),gradcy(*),gradcz(*)
      integer :: n,ijk
      real(double) :: rvol
!
      do 100 n=1,ncells
!
      ijk=ijkcell(n)
!
      rvol=1./(vol(ijk)+1.e-20)
!
      gradcx(ijk)=(c1x(ijk)*psi(ijk+iwid)     &
                +(c2x(ijk)*psi(ijk+iwid+jwid)     &
                +(c3x(ijk)*psi(ijk+jwid)     &
                +(c4x(ijk)*psi(ijk)     &
                +(c5x(ijk)*psi(ijk+iwid+kwid)     &
                +(c6x(ijk)*psi(ijk+iwid+jwid+kwid)     &
                +(c7x(ijk)*psi(ijk+jwid+kwid)     &
                +(c8x(ijk)*psi(ijk+kwid)))))))))    &
                  *rvol
!
      gradcy(ijk)=(c1y(ijk)*psi(ijk+iwid)     &
                +(c2y(ijk)*psi(ijk+iwid+jwid)     &
                +(c3y(ijk)*psi(ijk+jwid)     &
                +(c4y(ijk)*psi(ijk)     &
                +(c5y(ijk)*psi(ijk+iwid+kwid)     &
                +(c6y(ijk)*psi(ijk+iwid+jwid+kwid)     &
                +(c7y(ijk)*psi(ijk+jwid+kwid)     &
                +(c8y(ijk)*psi(ijk+kwid)))))))))    &
                  *rvol
!
      gradcz(ijk)=(c1z(ijk)*psi(ijk+iwid)      &
                +(c2z(ijk)*psi(ijk+iwid+jwid)      &
                +(c3z(ijk)*psi(ijk+jwid)      &
                +(c4z(ijk)*psi(ijk)      &
                +(c5z(ijk)*psi(ijk+iwid+kwid)      &
                +(c6z(ijk)*psi(ijk+iwid+jwid+kwid)      &
                +(c7z(ijk)*psi(ijk+jwid+kwid)      &
                +(c8z(ijk)*psi(ijk+kwid)))))))))    &
                  *rvol
!
  100 continue
!
      return
      end subroutine gradc
!     ********************************************************************
      subroutine gradf(nvtx,ijkvtx,     &
          f,gradxf,gradyf,gradzf)
!
!     a routine to calculate the gradient of f
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use geometry_com_M
      implicit real*8 (a-h,o-z)
!
      integer ::  &
          ijkvtx(*)
      real(double)  ::    &
          f(*),gradxf(*),gradyf(*),gradzf(*)
!
!      tiny=0.125*1.e-3
      tiny=0.0
!
      do 100 n=1,nvtx
!
      ijk=ijkvtx(n)
!
!
      gradxf(ijk)=-(c1x(ijk-iwid)*(f(ijk-iwid))     &
                +(c2x(ijk-iwid-jwid)*(f(ijk-iwid-jwid))     &
                +(c3x(ijk-jwid)*(f(ijk-jwid))     &
                +(c4x(ijk)*(f(ijk))     &
                +(c5x(ijk-iwid-kwid)*(f(ijk-iwid-kwid))     &
                +(c6x(ijk-iwid-jwid-kwid)*(f(ijk-iwid-jwid-kwid))     &
                +(c7x(ijk-jwid-kwid)*(f(ijk-jwid-kwid))     &
                +(c8x(ijk-kwid)*(f(ijk-kwid))))))))))
!
      gradyf(ijk)=-(c1y(ijk-iwid)*(f(ijk-iwid))     &
                +(c2y(ijk-iwid-jwid)*(f(ijk-iwid-jwid))     &
                +(c3y(ijk-jwid)*(f(ijk-jwid))     &
                +(c4y(ijk)*(f(ijk))     &
                +(c5y(ijk-iwid-kwid)*(f(ijk-iwid-kwid))     &
                +(c6y(ijk-iwid-jwid-kwid)*(f(ijk-iwid-jwid-kwid))     &
                +(c7y(ijk-jwid-kwid)*(f(ijk-jwid-kwid))     &
                +(c8y(ijk-kwid)*(f(ijk-kwid))))))))))
!
      gradzf(ijk)=-(c1z(ijk-iwid)*(f(ijk-iwid))     &
                +(c2z(ijk-iwid-jwid)*(f(ijk-iwid-jwid))     &
                +(c3z(ijk-jwid)*(f(ijk-jwid))     &
                +(c4z(ijk)*(f(ijk))     &
                +(c5z(ijk-iwid-kwid)*(f(ijk-iwid-kwid))     &
                +(c6z(ijk-iwid-jwid-kwid)*(f(ijk-iwid-jwid-kwid))     &
                +(c7z(ijk-jwid-kwid)*(f(ijk-jwid-kwid))     &
                +(c8z(ijk-kwid)*(f(ijk-kwid))))))))))
!
  100 continue
!
      return
      end subroutine gradf
!     *******************************************************************
      subroutine graddiv(ncells,ijkcell,   &
          ibp1,jbp1,kbp1,cdlt,sdlt,    &
          vol,    &
          ex,ey,ez,divetil,gdivx,gdivy,gdivz)
!
!     a routine to calculate the right-hand-side for Faraday's law
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use geometry_com_M
      implicit real*8 (a-h,o-z)
!
      integer :: ijkcell(*)
      real(double) ::    &
          vol(*),ex(*),ey(*),ez(*),divetil(*),    &
          gdivx(*),gdivy(*),gdivz(*)
!
!
!       ************************************************************************
!
!     CALCULATE GRAD(DIV(E))
!
      call list(1,ibp1,1,jbp1,2,kbp1,iwid,jwid,kwid,     &
          ncells,ijkcell)
!
      call divc(ncells,ijkcell,             &
                        vol,             &
                        ex,ey,ez,DIVETIL)
!
      call list(2,ibp1,2,jbp1,2,kbp1,iwid,jwid,kwid,     &
          ncells,ijkcell)
!
      call gradf(ncells,ijkcell,     &
          DIVETIL,GDIVX,GDIVY,GDIVZ)
!
      call bcperv(ibp1+1,jbp1+1,kbp1+1,iwid,jwid,kwid,	&
          sdlt,cdlt,	&
          GDIVX,GDIVY,GDIVZ)
!
!     it is assumed cartesian=.t.
!      call axisgrad(ibp1,jbp1,	&
!          GDIVX,GDIVY,GDIVZ)
!
!
!     END CALCULATION OF GRAD(DIV(E))
!
!     ********************************************************************
!
      return
      end subroutine graddiv
!     *******************************************************************
      subroutine list(i1,i2,j1,j2,k1,k2,iwid,jwid,kwid,     &
          nlist,ijklist)
!
!     a routine to form a singly indexed array of indices of active cells
!
!     nlist ... total number of active cells
!     ijklist ... cell indices
!
      implicit real*8 (a-h,o-z)
!
      integer :: ijklist(*)
!
      nlist=0
!
      do 1 k=k1,k2
      do 1 j=j1,j2
      do 1 i=i1,i2
!
      nlist=nlist+1
!
      ijklist(nlist)=1+(i-1)*iwid+(j-1)*jwid+(k-1)*kwid
!
    1 continue
!
      return
      end subroutine list
!     **************************************************************************
      subroutine torusbcv(nxp,nyp,nzp,    &
          cdlt,sdlt,strait,dz,    &
          periodic_x,periodic_y,periodic_z,    &
          x,y,z)
!
!     a routine to impose double periodicity
!     for toroidal geometry for vertex vector variables
!     assembles data
!
!     called by ACCEL_3DMHD, VINIT_GMRES
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
      logical periodic_x,periodic_y,periodic_z
!
      real(double) :: x(nxp,nyp,*),y(nxp,nyp,*),z(nxp,nyp,*)
!
      do 10 k=2,nzp
!
!     periodicity in the toroidal angle
!
      if(periodic_y) then
!
      do 1 i=2,nxp
!
      x(i,2,k)=cdlt*x(i,nyp,k)+sdlt*y(i,nyp,k)+x(i,2,k)
      y(i,2,k)=-sdlt*x(i,nyp,k)+cdlt*y(i,nyp,k)+y(i,2,k)    &
          -strait*dz
      z(i,2,k)=z(i,nyp,k)+z(i,2,k)
!
      x(i,nyp,k)=cdlt*x(i,2,k)-sdlt*y(i,2,k)
      y(i,nyp,k)=sdlt*x(i,2,k)+cdlt*y(i,2,k)    &
          +strait*dz
      z(i,nyp,k)=z(i,2,k)
!
    1 continue
!
      endif
!
!     periodicity in the poloidal angle
!
      if(periodic_x) then
!
      do 2 j=2,nyp
!
      x(2,j,k)=x(nxp,j,k)+x(2,j,k)
      y(2,j,k)=y(nxp,j,k)+y(2,j,k)
      z(2,j,k)=z(nxp,j,k)+z(2,j,k)
!
      x(nxp,j,k)=x(2,j,k)
      y(nxp,j,k)=y(2,j,k)
      z(nxp,j,k)=z(2,j,k)
!
    2 continue
!
      endif
!
   10 continue
!
      return
      end subroutine torusbcv
!     ********************************************************************
