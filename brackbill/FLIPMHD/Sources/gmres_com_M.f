       module gmres_com_M
!
       use vast_kind_param, ONLY : double
       use corgan_com_M 
!
       real(double), dimension(itdim) :: phi,diag,residu,wKrylov,Ax,  &
         rhs
       real(double) :: srce(itdim), Jdu(itdim)
       real(double) :: q(itdim,GMitmax)
!
       end module gmres_com_M
