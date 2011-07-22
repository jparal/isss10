      subroutine rotate(ibar,jbar,kbar,    &
          rmaj,dz,strait,iota,istep,    &
          delt,dtheta,dzstr,dphi,    &
          cdlt,sdlt,cdlhf,sdlhf,    &
          cdph,sdph,cdphhf,sdphhf)
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      real(double) :: iota(*)
      integer :: istep(*)
!     calculate elements of rotation matrix for toroidal direction
!
      pi=acos(-1.)
!
      rfibar=1./real(ibar)
      rfjbar=1./real(jbar)
!
         dtheta=2.*pi*rfibar
         delt=0.0
      if(rmaj.gt.0.0)   delt=dz/rmaj
         dphi = delt*rfjbar
         dzstr=strait*dz*rfjbar
!
         cdlt = cos(delt)
         sdlt = sin(delt)
         cdlhf= sqrt(.5*(1.+cdlt))
         sdlhf= sqrt(.5*(1.-cdlt))
         cdph = cos(dphi)
         sdph = sin(dphi)
         cdphhf = cos(dphi/2.)
         sdphhf = sin(dphi/2.)
!
!
!     define for quasi-ballooning coordinates
!
      do 1 k=1,kbar+2
!
      istep(k)=int(iota(k)*ibar+0.5)*delt/(2.*pi)
!
    1 continue
!
      return
      end subroutine rotate
