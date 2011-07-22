      subroutine bc_noslip(nvtx,ijkvtx,u,v,w)
      use vast_kind_param, ONLY : double
      implicit none
      integer, intent(in) :: nvtx, ijkvtx(*)
      real(double) :: u(*),v(*),w(*)
!
      integer ijk, n
!     a routine to impose no slip conditions on a surface
!     of the mesh
!
!
      do n=1,nvtx
!
      ijk=ijkvtx(n)
!
      u(ijk)=0.0
      v(ijk)=0.0
      w(ijk)=0.0
!
      enddo
!
      return
      end subroutine bc_noslip
