      subroutine volume_vtx(ncells,ijkcell,iwid,jwid,kwid,    &
         x,y,z,    &
          vvol)
!
!     a routine to calculate a semi-correct vertex volume
!
      use vast_kind_param, ONLY : double
      use geometry_com_M
      implicit real*8 (a-h,o-z)
!
      integer :: ijkcell(*)
      real(double) ::    &
          x(*),y(*),z(*),vvol(*)
!
      do 1 n=1,ncells
      ijk=ijkcell(n)
!
      xc=0.125*(x(ijk+iwid)    &
              +x(ijk+iwid+jwid)    &
              +x(ijk+jwid)    &
              +x(ijk)    &
              +x(ijk+iwid+kwid)    &
              +x(ijk+iwid+jwid+kwid)    &
              +x(ijk+jwid+kwid)    &
              +x(ijk+kwid))
!
      vvol(ijk+iwid)=vvol(ijk+iwid)     &
           +c1x(ijk)*(x(ijk+iwid)-xc)
!
      vvol(ijk+iwid+jwid)=vvol(ijk+iwid+jwid)    &
           +c2x(ijk)*(x(ijk+iwid+jwid)-xc)
!
      vvol(ijk+jwid)=vvol(ijk+jwid)     &
           +c3x(ijk)*(x(ijk+jwid)-xc)
!
      vvol(ijk)=vvol(ijk)    &
           +c4x(ijk)*(x(ijk)-xc)
!
      vvol(ijk+iwid+kwid)=vvol(ijk+iwid+kwid)    &
          +c5x(ijk)*(x(ijk+iwid+kwid)-xc)
!
      vvol(ijk+iwid+jwid+kwid)=vvol(ijk+iwid+jwid+kwid)    &
           +c6x(ijk)*(x(ijk+iwid+jwid+kwid)-xc)
!
      vvol(ijk+jwid+kwid)=vvol(ijk+jwid+kwid)     &
           +c7x(ijk)*(x(ijk+jwid+kwid)-xc)
!
      vvol(ijk+kwid)=vvol(ijk+kwid)     &
           +c8x(ijk)*(x(ijk+kwid)-xc)
!
    1 continue
!
      return
      end subroutine volume_vtx
