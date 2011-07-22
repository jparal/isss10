      subroutine jacobi(nxp,nyp,nzp,   &
          residu,diag,   &
          phi)
!
      use vast_kind_param, ONLY : double
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z
      use cindex_com_M, ONLY : iwid,jwid,kwid
      real(double) :: residu(*),diag(*),phi(*)
!
      integer :: nxp,nyp,nzp,    &
          istart,iend,jstart,jend,kstart,kend
!
      if(periodic_x) then
         istart=2
      else
         istart=3
      endif
         iend=nxp-1
!
      if(periodic_y) then
         jstart=2
      else
         jstart=3
      endif
         jend=nyp-1
!
      if(periodic_z) then
         kstart=2
      else
         kstart=3
      endif
         kend=nzp-1
!
      do i=istart,iend
         do j=jstart,jend
            do k=kstart,kend
!
            ijk=(k-1)*kwid+(j-1)*jwid+(i-1)*iwid+1
            phi(ijk)=phi(ijk)-0.20*residu(ijk)/diag(ijk)
!
            enddo
         enddo
      enddo
!
      return
      end subroutine jacobi
