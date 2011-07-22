!-----------------------------------------------------------------------
!
      module graf1d_class
!
! class library for describing 1d gks graphics
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: august 16, 2002
!
      use errors_class
      implicit none
      private
      public :: EXCEPTION, GRAF1D_ERR
      public :: GROPEN, GRCLOSE, SETNPLT
      public :: graf1d, new_graf1d, del_graf1d, display  
!
! nx = number of points plotted in each subarray
! real(scale)/aimag(scale) = power of 2 scale of x/y coordinate for plot
! real(clip)/aimag(clip) = flag for positive and/or negative x/y values
! label = long character string label for plot
! nsubs = number of subarrays to be plotted
! mksmarkerstyle = (1,2) descriptor is for (real,complex) display
! label2 = additional character string comment for plot
! lsubs = short labels for subarrays being plotted
      type graf1d
          integer :: nx
          complex :: scale, clip
          character(len=80) :: label
          integer :: style, nsubs, marker
          character(len=80) :: label2
          character(len=10), dimension(:), pointer :: lsubs
      end type graf1d
!
      character(len=6), private, save :: self = 'graf1d'
! Error codes, 1 = invalid number of points, 2 = too few arrays
! 3 = incorrect descriptor, 4 = array too small
! 5 = non-conforming arrays
      integer, save :: GRAF1D_ERR = 0
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GROPEN
         implicit none
         end subroutine
      end interface
      interface
         subroutine GRCLOSE
         implicit none
         end subroutine
      end interface
      interface
         subroutine SETNPLT(nplt,irc)
         implicit none
         integer :: nplt, irc
         end subroutine
      end interface
      interface
         subroutine DISPR(f,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs,chr,c&
     &hrs,irc)
         implicit none
         integer :: isc, ist, mks, nx, nxv, ngs, irc
         real :: xmin, xmax
         character(len=*) :: label, chr
         character(len=*), dimension(ngs) :: chrs
         real, dimension(nxv,ngs) :: f
         end subroutine
      end interface
      interface
         subroutine DISPC(f,g,label,zsc,zst,mks,nx,nxv,ngs,chr,chrs,irc)
         implicit none
         integer :: mks, nx, nxv, ngs, irc
         complex :: zsc, zst
         character(len=*) :: label, chr
         character(len=*), dimension(ngs) :: chrs
         real, dimension(nxv,ngs) :: f, g
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface new_graf1d
         module procedure new_graf1d_r
         module procedure new_graf1d_c
      end interface
!
      interface display
         module procedure idispr
         module procedure idispc
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine new_graf1d_r(this,nx,xmin,xmax,label,nsubs,scale,cli&
     &p,marker,label2,lsubs)
! this subroutine creates a real descriptor for 1d graphics
! this = graf1d descriptor of 1d graphics data
! nx = number of points plotted in each subarray
! xmin/xmax = range of x values in plot
! label = long character string label for plot
! nsubs = number of subarrays to be plotted 
! scale = power of 2 scale of y coordinate for plot,
! if abs(scale) > 116, then the program finds the minimum value of scale
! clip = (1,0,-1) = flag for choosing (positive,all,negative) values
! clip = 2, uses largest power of 2 between minimum and maxmium values
! marker = flag to determine whether lines or markers are used
! label2 = additional character string comment for plot
! lsubs = array of ngs short character labels
! default values: nsubs=1,scale=999,clip=0,marker=0,label2=' ',lsubs=' '
         implicit none
         type (graf1d) :: this
         integer :: nx
         real :: xmin, xmax
         character(len=*) :: label
         integer, optional :: nsubs, scale, clip, marker
         character(len=*), optional :: label2
         character(len=*), dimension(:), optional :: lsubs
         if (monitor==2) call werrfl('new_graf1d_r started')
! check for errors
         if (monitor > 0) then
         if (nx < 1) then
            this%style = 0
            write (erstr,*) 'new_graf1d_r: invalid number of points, nx &
     &= ', nx
            GRAF1D_ERR = 1; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         endif
! save required descriptor values
         this%nx = nx
         this%scale = cmplx(xmin,999.)
         this%clip = cmplx(xmax,0.)
         this%label = label
         this%style = 1
! set optional descriptor values, reset to defaults if out of range
         this%nsubs = 1
         if (present(nsubs)) then
            if (nsubs >= 1) then
               this%nsubs = nsubs
            else if (monitor > 0) then
               this%style = 0
               write (erstr,*) 'new_graf1d_r: too few arrays, nsubs = ',
     & nsubs
               GRAF1D_ERR = 2; EXCEPTION = EXCEPTION + 1
               call ehandler(EDEFAULT,erstr); return
            endif
         endif
         if (present(scale)) this%scale = cmplx(xmin,real(scale))
         if (present(clip)) then
            if ((clip >= (-1)) .and. (clip <= 2)) then
               this%clip = cmplx(xmax,real(clip))
            endif
         endif
         this%marker = 0
         if (present(marker)) this%marker= marker
         this%label2 = ' '
         if (present(label2)) this%label2 = label2
! allocate character array
         allocate(this%lsubs(this%nsubs))
         call count_alloc((size(this%lsubs)),PTRS)
         if (present(lsubs)) then
            this%lsubs = lsubs
         else
            this%lsubs = ' '
         endif
         if (monitor==2) call werrfl('new_graf1d_r complete')
         end subroutine new_graf1d_r
!
         subroutine new_graf1d_c(this,nx,label,nsubs,scale,clip,marker,l&
     &abel2,lsubs)
! this subroutine creates a complex descriptor for 1d graphics
! this = graf1d descriptor of 1d graphics data
! nx = number of points plotted in each subarray
! label = long character string label for plot
! nsubs = number of subarrays to be plotted
! real(scale)/aimag(scale) = power of 2 scale of x/y coordinate for plot
! if abs(scale) > 116, then the program finds the minimum value of scale
! real(clip)/aimag(clip) = flag for positive and/or negative x/y values
! clip = (1,0,-1) = flag for choosing (positive,all,negative) values
! clip = 2, uses largest power of 2 between minimum and maxmium values
! marker = flag to determine whether lines or markers are used
! label2 = additional character string comment for plot
! lsubs = array of ngs short character labels
! default values: nsubs=1,scale=(999,999),clip=(0,0),marker=0,label2=' '
! lsubs=' '
         implicit none
         type (graf1d) :: this
         integer :: nx
         character(len=*) :: label
         integer, optional :: nsubs, marker
         complex, optional :: scale, clip
         character(len=*), optional :: label2
         character(len=*), dimension(:), optional :: lsubs
         if (monitor==2) call werrfl('new_graf1d_c started')
! check for errors
         if (monitor > 0) then
         if (nx < 1) then
            this%style = 0
            write (erstr,*) 'new_graf1d_c: invalid number of points, nx &
     &= ', nx
            GRAF1D_ERR = 1; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         endif
! save required descriptor values
         this%nx = nx
         this%label = label
         this%style = 2
! set optional descriptor values, reset to defaults if out of range
         this%nsubs = 1
         if (present(nsubs)) then
            if (nsubs >= 1) then
               this%nsubs = nsubs
            else if (monitor > 0) then
               this%style = 0
               write (erstr,*) 'new_graf1d_r: too few arrays, nsubs = ',&
     & nsubs
               GRAF1D_ERR = 2; EXCEPTION = EXCEPTION + 1
               call ehandler(EDEFAULT,erstr); return
            endif
         endif
         this%scale = cmplx(999.,999.)
         if (present(scale)) this%scale = scale
         this%clip = cmplx(0.,0.)
         if (present(clip)) then
            if ((real(clip) >= (-1.)) .and. (real(clip) <= 2.)) then
               this%clip = cmplx(real(clip),aimag(this%clip))
            endif
            if ((aimag(clip) >= (-1.)) .and. (aimag(clip) <= 2.)) then
               this%clip = cmplx(real(this%clip),aimag(clip))
            endif
         endif
         this%marker = 0
         if (present(marker)) this%marker = marker
         this%label2 = ' '
         if (present(label2)) this%label2 = label2
! allocate character array
         allocate(this%lsubs(this%nsubs))
         call count_alloc((size(this%lsubs)),PTRS)
         if (present(lsubs)) then
            this%lsubs = lsubs
         else
            this%lsubs = ' '
         endif
         if (monitor==2) call werrfl('new_graf1d_c complete')
         end subroutine new_graf1d_c
!
         subroutine del_graf1d(this)
! delete descriptor for 1d graphics
! this = graf1d descriptor of 1d graphics data
         implicit none
         type (graf1d) :: this
         if (monitor==2) call werrfl('del_graf1d started')
! deallocate descriptor
         if (associated(this%lsubs)) then
            call count_alloc(-(size(this%lsubs)),PTRS)
            deallocate(this%lsubs)
         endif
         if (monitor==2) call werrfl('del_graf1d complete')
         end subroutine del_graf1d
!
         subroutine idispr(this,f,irc,scale,clip,marker,label2)
! this subroutine plots ngs subarrays of the array f, on a common graph,
! each plot with nx points, versus a linear function in x,
! where xmin < x < xmax.
! this = graf1d descriptor of 1d graphics data
! f = array containing subarrays to be plotted
! irc = return code (0 = normal return)
! scale = power of 2 scale of y coordinate for plot,
! if abs(scale) > 116, then the program finds the minimum value of scale
! clip = (1,0,-1) = flag for choosing (positive,all,negative) values
! clip = 2, uses largest power of 2 between minimum and maxmium values
! marker = flag to determine whether lines or markers are used
! label2 = additional character string comment for plot
         type (graf1d), intent(in) :: this
         integer :: irc
         real, dimension(:,:) :: f
         integer, optional :: scale, clip, marker
         character(len=*), optional :: label2
! local data
         integer :: nxv, isc, ist, mks
         character(len=80) :: chr
         real :: xmin, xmax
         if (monitor==2) call werrfl('idispr started')
         nxv = size(f,1)
! check for errors
         if (monitor > 0) then
         if (this%style /= 1) then
            write (erstr,*) 'idispr: incorrect descriptor for real displ&
     &ay'
            GRAF1D_ERR = 3; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         if (this%nx > nxv) then
            write (erstr,*) 'idispr: array too small, nx, nxv=', this%nx&
     &, nxv
            GRAF1D_ERR = 4; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         if (this%nsubs > size(f,2)) then
            write (erstr,*) 'idispr: too few arrays, nsubs, ndim = ', th&
     &is%nsubs, size(f,2)
            GRAF1D_ERR = 2; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         endif
! unpack arguments
         xmin = real(this%scale); xmax = real(this%clip)
         isc = aimag(this%scale); if (present(scale)) isc = scale
         ist = aimag(this%clip);  if (present(clip)) ist = clip
         mks = this%marker; if (present(marker)) mks = marker
         chr = this%label2; if (present(label2)) chr = label2
! call function
         call DISPR(f,trim(this%label),xmin,xmax,isc,ist,mks,this%nx,nxv&
     &,this%nsubs,trim(chr),this%lsubs,irc)
         if (monitor==2) call werrfl('idispr complete')
         end subroutine idispr
!
         subroutine idispc(this,f,g,irc,scale,clip,marker,label2)
! this subroutine plots ngs subarrays of the array f, on a common graph,
! each plot with nx points, versus the corresponding subarray of the
! array g
! this = graf1d descriptor of 1d graphics data
! f, g = arrays containing subarrays to be plotted
! irc = return code (0 = normal return)
! real(scale)/aimag(scale) = power of 2 scale of x/y coordinate for plot
! if abs(scale) > 116, then the program finds the minimum value of scale
! real(clip)/aimag(clip) = flag for positive and/or negative x/y values
! clip = (1,0,-1) = flag for choosing (positive,all,negative) values
! clip = 2, uses largest power of 2 between minimum and maxmium values
! marker = flag to determine whether lines or markers are used
! label2 = additional character string comment for plot
         implicit none
         type (graf1d), intent(in) :: this
         integer :: irc
         real, dimension(:,:) :: f, g
         complex, optional :: scale, clip
         integer, optional :: marker
         character(len=*), optional :: label2
! local data
         integer :: nxv, mks
         complex :: zsc, zst
         character(len=80) :: chr
         if (monitor==2) call werrfl('idispc started')
         nxv = size(f,1)
! check for errors
         if (monitor > 0) then
         if (this%style /= 2) then
            write (erstr,*) 'idispc: incorrect descriptor for complex di&
     &splay'
            GRAF1D_ERR = 3; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         if (this%nx > nxv) then
            write (erstr,*) 'idispc: array too small, nx, nxv=', this%nx&
     &, nxv
            GRAF1D_ERR = 4; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         if (this%nsubs > size(f,2)) then
            write (erstr,*) 'idispc: too few f arrays, nsubs, ndim = ', &
     &this%nsubs, size(f,2)
            GRAF1D_ERR = 2; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         if (nxv /= size(g,1)) then
            write (erstr,*) 'idispc: non-conforming arrays'
            GRAF1D_ERR = 5; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         if (this%nsubs > size(g,2)) then
            write (erstr,*) 'idispc: too few g arrays, nsubs, ndim = ', &
     &this%nsubs, size(g,2)
            GRAF1D_ERR = 4; EXCEPTION = EXCEPTION + 1
            call ehandler(EDEFAULT,erstr); return
         endif
         endif
! unpack arguments
         zsc = aimag(this%scale); if (present(scale)) zsc = scale
         zst = aimag(this%clip);  if (present(clip)) zst = clip
         mks = this%marker; if (present(marker)) mks = marker
         chr = this%label2; if (present(label2)) chr = label2
! call function
         call DISPC(f,g,trim(this%label),zsc,zst,mks,this%nx,nxv,this%ns&
     &ubs,trim(chr),this%lsubs,irc)
         if (monitor==2) call werrfl('idispc complete')
         end subroutine idispc
!
         subroutine idisplayw(wt,dtw,nt,irc)
! this subroutine displays time history of electric field, kinetic, and
! total energies
! wt = time history array for energies
! dtw = time between energy values
! nt = number of energy values to be displayed
! irc = return code (0 = normal return)
         implicit none
         integer :: nt, irc
         real :: dtw
         real, dimension(:,:) :: wt
! local data
         real :: tmn, tmx
         character*36 lbl
         character(len=10), dimension(3) :: chrs = (/'ELEC FIELD','ELECT&
     & ENRG','TOTAL ENRG'/)
         type (graf1d) :: line, mline
         if (monitor==2) call werrfl('idisplayw started')
! tmn/tmx = range of time values in plot
         tmn = 0.
         tmx = dtw*(nt - 1)
         call new_graf1d(line,nt,tmn,tmx,' ',clip=1)
! electric field energy
         lbl = ' ELECTRIC FIELD ENERGY VERSUS TIME '
         call display(line,wt(:,1:1),irc,label2=lbl)
         if (irc==1) return
! electron kinetic energy
         lbl = ' ELECTRON KINETIC ENERGY VERSUS TIME'
         call display(line,wt(:,2:2),irc,label2=lbl)
         if (irc==1) return
! total energy
         lbl = ' TOTAL ENERGY VERSUS TIME '
         call display(line,wt(:,3:3),irc,label2=lbl)
         if (irc==1) return
! all energies on common plot
         call new_graf1d(mline,nt,tmn,tmx,' ',nsubs=3,clip=1,lsubs=chrs)
         lbl = ' ENERGIES VERSUS TIME '
         call display(mline,wt,irc,label2=lbl)
         if (monitor==2) call werrfl('idisplayw complete')
         return
         end subroutine idisplayw


!
      end module graf1d_class

