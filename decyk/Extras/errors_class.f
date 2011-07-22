!-----------------------------------------------------------------------
!
      module errors_class
!
! this class provides error handling and debugging utilities
! written by viktor k. decyk, ucla
! copyright 2001, regents of the university of california
! update: december 28, 2001
!
      implicit none
      private
      public :: NOERR, LOGIT, QUIT, EXCEPTION, EDEFAULT
      public :: eunit, erstr, monitor
      public :: INTS, REALS, DOUBLES, CMPLXS, PTRS
      public :: set_eunit, set_edefault, ehandler, werrfl, set_monitor
      public :: clr_alloc, count_alloc, pr_alloc
!
      integer, parameter :: NOERR = 0, LOGIT = 0, QUIT = 1
! EXCEPTION = flag if error has occurred
! EDEFAULT = default behavior for error handler ehandler
      integer, save :: EXCEPTION = NOERR, EDEFAULT = QUIT
! default Fortran unit number for error file
      integer, save :: eunit = 2
      character(len=80), save :: erstr = ' '
! monitor = (0,1,2) = (disable,normal,extended) error processing
      integer, save :: monitor = 1
! constants and counters for memory allocation of various types
      integer, parameter :: INTS = 1, REALS = 2, DOUBLES = 3
      integer, parameter :: CMPLXS = 4, PTRS = 5
      integer, save :: int_alloc = 0, real_alloc = 0, double_alloc = 0
      integer, save :: cmplx_alloc = 0, ptr_alloc = 0, unkwn_alloc = 0
!
      contains
!
         subroutine set_eunit(iunit)
! this subroutine sets the unit number for writing error file
! iunit = new unit number
         integer, intent(in) :: iunit
         if (iunit.ge.0) eunit = iunit
         end subroutine set_eunit
!
         subroutine set_edefault(newdefault)
! this subroutine sets default behavior for error handler
! newdefault= new default behavior
         integer, intent(in) :: newdefault
         if ((newdefault >= 0) .and. (newdefault <= 1)) then
            EDEFAULT = newdefault
         endif
         end subroutine set_edefault
!
         subroutine ehandler(how,estr)
! this subroutine handles errors and optionally logs error message
! how = keyword on how to handle error
! estr = optional error message string to be logged in error file
         integer, intent(in) :: how
         character(len=*), intent(in), optional :: estr
         if (present(estr)) write (eunit,*) trim(estr)
         if (how==QUIT) stop
         EXCEPTION = NOERR
         end subroutine ehandler
!
         subroutine werrfl(estr)
! this subroutine prints errors 
! estr = optional error message string to be logged in error file
         character(len=*), intent(in) :: estr
         write (eunit,*) trim(estr)
         end subroutine werrfl
!
         subroutine set_monitor(monval)
! this subroutine sets new monitor value
! monval = new monitor value
         implicit none
         integer, intent(in) :: monval
! reset monitor value
         if (monval > 1) then
            monitor = 2
         else if (monval < 1) then
            monitor = 0
         else
            monitor = 1
         endif
         end subroutine set_monitor
!
         subroutine clr_alloc()
! this subroutine clears memory allocation counters
         int_alloc = 0
         real_alloc = 0
         double_alloc = 0
         cmplx_alloc = 0
         ptr_alloc = 0
         unkwn_alloc = 0
         end subroutine clr_alloc
!
         subroutine count_alloc(num,atype)
! this subroutine keeps track of memory allocations for debugging
! num = number of entities allocated, atype = type of entity allocated
! num should be negative for deallocations
         integer, intent(in) :: num, atype
         if (atype==INTS) then
            int_alloc = int_alloc + num
         else if (atype==REALS) then
            real_alloc = real_alloc + num
         else if (atype==DOUBLES) then
            double_alloc = double_alloc + num
         else if (atype==CMPLXS) then
            cmplx_alloc = cmplx_alloc + num
         else if (atype==PTRS) then
            ptr_alloc = ptr_alloc + num
         else
            unkwn_alloc = unkwn_alloc + num
         endif
         end subroutine count_alloc
!
         subroutine pr_alloc()
! this subroutine prints out current allocation status
         integer :: atotal
         atotal = int_alloc + real_alloc + double_alloc + cmplx_alloc + &
     &ptr_alloc + unkwn_alloc
         if (atotal==0) then
            write (6,*) 'Current Memory Allocation:', atotal
         else
            write (6,*) 'Current Memory Allocation:'
            if (int_alloc /= 0) write (6,*) 'integers = ', int_alloc
            if (real_alloc /= 0) write (6,*) 'reals = ', real_alloc
            if (double_alloc /= 0) write (6,*) 'doubles = ',double_alloc
            if (cmplx_alloc /= 0) write (6,*) 'cmplxs = ', cmplx_alloc
            if (ptr_alloc /= 0) write (6,*) 'ptrs = ', ptr_alloc
            if (unkwn_alloc /= 0) write (6,*) 'unknown = ', unkwn_alloc
         endif
         end subroutine pr_alloc
!
      end module errors_class
