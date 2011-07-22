      module nk_com_M
      use vast_kind_param, only:  double
      use corgan_com_M
      use nkmhd_com_M
      real(double) ::    &
          pq(itdim,num_eq),    &
          ptil(itdim,num_eq),  &
          apq(itdim,num_eq),   &
          solution(itdim,num_eq),s(itdim,num_eq),   &
          dp(itdim,num_eq),    &
          celerr(itdim)
      end module nk_com_M
