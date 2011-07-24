      module corgan_com_M
      use vast_kind_param, only:  double
!
      integer, parameter :: nhst=5000,ijmx=2000
      integer, parameter :: nsampl=1024
      integer, parameter :: idx=32,idy=32,idz=64
      integer, parameter :: idxg=idx+2,idyg=idy+2,idzg=idz+2
      integer, parameter :: itdim=idxg*idyg*idzg
      integer, parameter :: idxyzg=idxg*idyg*idzg
      integer, parameter :: npart=1769472,npart_object=5
      integer, parameter :: nspecies_object=1
      integer, parameter :: nreg_object=3
      integer, parameter :: nlist=150
      integer, parameter :: nreg=3
      integer, parameter :: nspecies=2
      integer, parameter :: GMitmax=10
      integer, dimension(4) :: idisp
      integer ::             &
             krd   , kpr   , kpt   , kfm   , kpl   , ktpin , ktpout,             &
             jbnm  , itlm  , numtd , mgeom ,              &
             lpr   ,              &
             ncyc  , ncyc1 , ncyc2 ,              &
             ixto
      logical ::    touch,frctn,outflg,wrtp,periodic,   &
          newstuf,PRECON,     &
          nomag,expnd
!
      real(double), parameter ::  clite=1.,fourpi=1.,cntr=1.0
      real(double) ::                       &
             dto(10),dtoc(10),             &
             dcttd , dpttd , pttd  ,             &
             tout  , tmovi,             &
             pars  , strait, toroid,             &
             tlimd , twfin , tt,             &
             ddate , hour  ,             &
             iota(idzg) ,             &
             tlm   , stim  , t1    , t2    ,             &
             pflag
!
      end module corgan_com_M
