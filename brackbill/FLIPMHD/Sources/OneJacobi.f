      subroutine OneJacobi(GMit,eps,phi,du, F,Jdu)
!
      use vast_kind_param, ONLY: double
      use blcom_com_M, ONLY : vol, ijkvtx
      use corgan_com_M, ONLY : itdim
      use cindex_com_M, ONLY : nvtx
      use gmres_com_M, ONLY : srce
     
      implicit none
      integer :: ij,n, GMit
      real(double) eps,reps
      real(double) :: Jdu(*)
      real(double) :: phi(*),    &
          du(itdim,*), F(*)
      logical :: complement

!     where we evaluate

!            F(u + eps.du) - F(u)
!     J.du = -------------------
!                    eps

!     and the new residual r = -F(u)-J.du

!     F(u) from the current Newton iteration, find F(u+eps.du) = Feps
!
      complement=.false.
!
!
      do n = 1,nvtx
         ij = ijkvtx(n)
         phi(ij) = phi(ij)+eps*du(ij,GMit)
      end do
!
!
      call residue_vtx(nvtx,ijkvtx,srce, phi, Jdu)
!
      reps=1.d0/eps
!
      do n=1,nvtx
         ij=ijkvtx(n)
         phi(ij) = phi(ij)-eps*du(ij, GMit)
         Jdu(ij) = (Jdu(ij)-F(ij))*reps
      enddo

!
      return
      end

