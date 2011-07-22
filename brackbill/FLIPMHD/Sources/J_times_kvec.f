      subroutine J_times_kvec(pq)
!
!     a routine to approximate the Jacobian  times
!     a krylov vector (pq,bxq,byq,bzq) for the 3D MHD equations.
!
!     itdim is the maximum array length set in corgan_com
!     if reversed order of indices in solution, for example, could eliminate
!     need for itdim an argument list
!     pq is krylov vector
!
      use vast_kind_param, only:  double
      use cindex_com_M, ONLY : ncells
      use corgan_com_M, ONLY : itdim
      use nkmhd_com_M, ONLY : num_eq
      use nk_com_M, ONLY : solution, apq, ptil, s
      use blcom_com_M, ONLY : p,csq,ijktmp2
      implicit none
!
      integer ::   &
         n,m,jj
!
      real(double) ::     &
          pq(itdim,num_eq)
!
      real(double) :: eps_NK, vec_scale
!
!
!****************************************************************************
!
!     calculate eps_NK
!
      eps_NK=0.0e0
      vec_scale=0.0e0
!
      do  n=1,ncells
         do jj=1,num_eq
            eps_NK= eps_NK +     &
               (abs(solution(n,jj)))
            vec_scale= vec_scale +    &
               (abs(pq(n,jj)))
         enddo
      enddo
!
      eps_NK = 1.0e-6*eps_NK/(vec_scale+1.d-20) 
!
!     calculate physics solution vector plus
!     perturbed Krylov vector. 
!
      do  n=1,ncells
         do jj=1,num_eq
            ptil(n,jj)  = solution(n,jj) + eps_NK*pq(n,jj)
         enddo
      enddo
!
!     calculate perturbed values for residual errors
!     and store them temporarily in apq
!

      call residu_3dmhd(    &
          ptil,apq) 

!     Now approximate Jacobian matrix times the Krylov vector, q,
!     and store in apq
!
      do m=1,num_eq
         do n=1,ncells
         apq(n,m) = (apq(n,m) - s(n,m))/eps_NK
         enddo
      enddo
!
!   Return to GMRES with J x q stored in
!   apq
!
      return
      end subroutine J_times_kvec




