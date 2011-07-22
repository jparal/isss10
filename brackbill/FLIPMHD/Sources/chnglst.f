      subroutine chnglst(np,iphed0,iphed1,link)
!
      implicit none
      integer, intent(in) :: np
      integer iphed0, iphed1, link
!
!
!     a routine for moving a particle from the beginning of one
!     list to the beginning of another
!
!     the arguments of chnglst are:
!     np...the index of a particle changing cells
!     iphed0...header for the donating list
!     iphed1...header for the receiving list
!     link...the pointer for particle np
!
      iphed0=link
      link=iphed1
      iphed1=np
!
      return
      end subroutine chnglst
