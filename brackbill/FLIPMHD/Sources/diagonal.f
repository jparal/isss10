      subroutine diagonal_vtx(nvtx,ijkvtx)
      use vast_kind_param, ONLY : double
      use blcom_com_M, ONLY : vol
      use cindex_com_M, ONLY : iwid,jwid,kwid
      use geometry_com_M
      use gmres_com_M, ONLY : diag
      implicit none
!
!     a routine to calculate the diagonal elements of the
!     matrix representing the Laplacian
!
      integer :: nvtx, ijkvtx(*)
      integer ijk, n
!
      write(*,*) 'diagonal_vtx: nvtx=',nvtx
      do 1 n=1,nvtx
      ijk=ijkvtx(n)
!
      diag(ijk)=    &
         +((c1x(ijk-iwid)**2+c1y(ijk-iwid)**2+c1z(ijk-iwid)**2)    &
         /(vol(ijk-iwid)+1.d-50)    &
         +(c2x(ijk-iwid-jwid)**2+c2y(ijk-iwid-jwid)**2    &
         +c2z(ijk-iwid-jwid)**2)/(vol(ijk-iwid-jwid)+1.d-50)    &
         +(c3x(ijk-jwid)**2+c3y(ijk-jwid)**2+c3z(ijk-jwid)**2)    &
         /(vol(ijk-jwid)+1.d-50)    &
         +(c4x(ijk)**2+c4y(ijk)**2+c4z(ijk)**2)/(vol(ijk)+1.d-50)    &
         +(c5x(ijk-iwid-kwid)**2+c5y(ijk-iwid-kwid)**2    &
         +c5z(ijk-iwid-kwid)**2)/(vol(ijk-iwid-kwid)+1.d-50)    &
         +(c6x(ijk-iwid-jwid-kwid)**2    &
         +c6y(ijk-iwid-jwid-kwid)**2+c6z(ijk-iwid-jwid-kwid)**2)    &
         /(vol(ijk-iwid-jwid-kwid)+1.d-50)    &
         +(c7x(ijk-jwid-kwid)**2+c7y(ijk-jwid-kwid)**2    &
         +c7z(ijk-jwid-kwid)**2)/(vol(ijk-jwid-kwid)+1.d-50)    &
         +(c8x(ijk-kwid)**2+c8y(ijk-kwid)**2    &
         +c8z(ijk-kwid)**2)/(vol(ijk-kwid)+1.d-50))
!
    1 continue
!
      return
       end subroutine diagonal_vtx
