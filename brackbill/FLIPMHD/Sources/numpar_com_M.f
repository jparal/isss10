      module numpar_com_M
      use vast_kind_param, only : double
      character*1  :: iddt
      integer :: numgrididdt
      real(double)  ::    &
             mu    , lam   ,   &
             a0    , b0    , colamu,   &
             t     , dt    , rdt   , dtv   , dtc   , dtvsav,  &
             pi    ,     &
             taug,r4pi,   &
             dtcon,dtvis,dtgrow,dtmin,dtpdv,dtmax,   &
             dtiter,dtresist
      integer ::    &
             numgrid
!
      end module numpar_com_M
