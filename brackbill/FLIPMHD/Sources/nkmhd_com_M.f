      module nkmhd_com_M
      use vast_kind_param, only:  double
      use corgan_com_M, only: itdim
!     TYPE DECLARATIONS
!
!
      integer, parameter :: num_eq=4
      integer, parameter :: itsub_mhd=20
!
      real(double), dimension(num_eq*itdim) ::   &
               Jq
      real(double), dimension(num_eq*itdim,itsub_mhd) ::   &
               q
      real(double) :: g(itsub_mhd+1,2),s_local(itsub_mhd+1),   &
               h(itsub_mhd+1,itsub_mhd+1)     
      real(double) :: ynorm,dxnorm,rnorm

      end module nkmhd_com_M

! 
