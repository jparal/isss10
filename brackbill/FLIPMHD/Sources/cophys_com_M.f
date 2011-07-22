          module cophys_com_M
!
          use vast_kind_param, only:  double
          real(double), dimension (100) :: rwz
          real(double) ::      &
             gm1   , rho0  , beta  , asq   , a     , wsa   ,     &
             tke   , tie   , etotl , tmass , tmomx , tmomy , tmomz ,    &
             time  , flux,   surface_area,   &
             toti  , totk,   totkz,   JouleHeating,    &
             totmb,  totmm,   totbb, diffbb,  tote,    &
             tparke, tgheat, gcur_max,    &
             tmom  , tgheat_par, tgheat_per,     &
             tbe   , tle   , grind,    &
             dx    , dy    , dz   ,    &
             rpl   , rvac  , rwall , rmaj  , cdlhf , sdlhf ,    &
             cdlt  , sdlt  , cdph  , sdph  , cdphhf, sdphhf, dphi  ,    &
             tmax  , tmax0 , tmax1 , tmax2 , tmax3 ,      &
             btorus, bvertcl,rhoi  , siei  ,      &
             racc  , zacc  , tflux , pflux ,      &
             gx    , gy    , gz    , tiesav, tkesav, tbesav,    &
             rq0   , q0,  coll,  efld0,     &
             c5last, dvbmax

        end module cophys_com_M
