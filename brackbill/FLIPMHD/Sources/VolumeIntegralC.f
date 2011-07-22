      subroutine VolumeIntegralC(q, IntegralC)
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : ncells
      use blcom_com_M, ONLY : vol, ijkcell
!-----------------------------------------------------------------
!L o c a l  V a r i a b l e s
!-----------------------------------------------------------------
      integer :: ijk, n
      real(double) :: q(*), IntegralC
!     a routine for computing the volume integral of vertex-centered variables
!
      IntegralC=0.0d0
           do n=1, ncells
              ijk=ijkcell(n)
              IntegralC=IntegralC+q(ijk)*vol(ijk)
           enddo
     return
     end subroutine VolumeIntegralC

