      subroutine triple(x,y,z,ijk,iwid,jwid,kwid,vol)
!
!     a routine to calculate the sub-volume of a hexagon
!
      use vast_kind_param, ONLY : double
      implicit real*8 (a-h,o-z)
!
      real(double) :: x(*),y(*),z(*)
!
      vol=((x(ijk+iwid)-x(ijk))*(y(ijk+jwid)-y(ijk))     &
         -(y(ijk+iwid)-y(ijk))*(x(ijk+jwid)-x(ijk)))     &
         *(z(ijk+kwid)-z(ijk))     &
!
        +((y(ijk+iwid)-y(ijk))*(z(ijk+jwid)-z(ijk))     &
         -(z(ijk+iwid)-z(ijk))*(y(ijk+jwid)-y(ijk)))     &
         *(x(ijk+kwid)-x(ijk))     &
!
        +((z(ijk+iwid)-z(ijk))*(x(ijk+jwid)-x(ijk))     &
         -(x(ijk+iwid)-x(ijk))*(z(ijk+jwid)-z(ijk)))     &
         *(y(ijk+kwid)-y(ijk))
!
      return
      end subroutine triple
