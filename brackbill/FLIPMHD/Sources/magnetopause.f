      subroutine magnetopause(ncells,ijkcell,iwid,jwid,kwid,    &
          vol,    &
          color,gradx,grady,gradz,    &
          bxn,byn,bzn,cjump,    &
          reconnected_flux,surface_area)
!
      use vast_kind_param, only:  double
      integer ncells,ijkcell(*),iwid,jwid,kwid
!
      real(double) ::     &
          color(*),gradx(*),grady(*),gradz(*),    &
          bxn(*),byn(*),bzn(*),vol(*),cjump,    &
          reconnected_flux,surface_area

!
!     a routine to calculate the reconnected flux at the magnetopause
!     the magnetic field is projected on to the gradient of color,
!     a vertex-centered characteristic function that has different constant
!     values in the magnetosphere and solar wind
!
      if(cjump.eq.0.d0) cjump=1.0d0
      gradcsq=0.0d0
      bzsq=0.0d0
!
!     compute the gradient of color
!
      call gradc(ncells,ijkcell,    &
          vol,    &
          color,gradx,grady,gradz)
!
!     project the magnetic field on to the gradient of color,
!     and sum the reconnected flux
!
      reconnected_flux=0.0d0
      surface_area=0.0d0
!
      do n=1,ncells
      ijk=ijkcell(n)
!
      reconnected_flux=reconnected_flux    &
          +(dabs(bxn(ijk))*gradx(ijk)    &
           +dabs(byn(ijk))*grady(ijk)    &
           +dabs(bzn(ijk))*gradz(ijk))*vol(ijk)

      surface_area=surface_area    &
         +dsqrt(gradx(ijk)**2+grady(ijk)**2+gradz(ijk)**2)    &
         *vol(ijk)

!
      enddo
!
      return
      end subroutine magnetopause
