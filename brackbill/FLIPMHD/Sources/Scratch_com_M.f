     module Scratch_com_M
!
     use vast_kind_param, ONLY : double
     use corgan_com_M, ONLY : itdim, npart
     real(double), dimension(:), allocatable :: upv1, vpv1, wpv1
     real(double), dimension(:), allocatable :: axv1, ayv1, azv1
     real(double), dimension(:), allocatable :: workp
     real(double), dimension(:), allocatable :: gradcx, gradcy, gradcz
     real(double), dimension(:), allocatable :: plotx, ploty, plotz
     real(double), dimension(:), allocatable :: dudxp, dudyp, dudzp
     real(double), dimension(:), allocatable :: dvdxp, dvdyp, dvdzp
     real(double), dimension(:), allocatable :: dwdxp, dwdyp, dwdzp
     real(double), dimension(:), allocatable :: the, zeta, nu
     real(double), dimension(:,:), allocatable :: mc_tmp, sie1p_tmp,  &
        bxn_tmp, byn_tmp, bzn_tmp, number_tmp
     real(double), dimension(:,:), allocatable :: mv_tmp, umom_tmp,  &
        vmom_tmp, wmom_tmp, numberv_tmp, color_tmp
!
     end module Scratch_com_M
