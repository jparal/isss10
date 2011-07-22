!-----------------------------------------------------------------------
! 1d PIC library for additional diagnostics
! Wrappers to allow some procedures or parts from diag1mod.f to be
! called from C and Java
! written by viktor k. decyk, ucla
! copyright 2011, regents of the university of california
! update: july 6, 2011
!-----------------------------------------------------------------------
!
      integer function fgetlrec(f,nx)
! find record length for real data for direct access files
      implicit none
      integer :: nx
      real, dimension(*) :: f
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1)
      fgetlrec = nx*lrec
      end function
!
      integer function fgetclrec(f,nx)
! find record length for complex data for direct access files
      implicit none
      integer :: nx
      complex, dimension(*) :: f
! local data
      integer :: lrec
      inquire(iolength=lrec) f(1)
      fgetclrec = nx*lrec
      end function
!
      integer function fgetunit(iunit)
! find an unconnected fortran unit number
      implicit none
      integer :: iunit
! local data
      integer :: i
      logical :: connected
      fgetunit = -1
! check connection status
      do i = iunit, 99
         inquire(unit=i,opened=connected)
         if (.not.connected) then
            fgetunit = i
            exit
         endif
      enddo
      end function
!
      subroutine fbfopen(f,nx,iunit,nrec,fname)
! open direct access binary file for real 1d scalar data
      implicit none
      integer :: nx, iunit, nrec
      real, dimension(*) :: f
      character(len=*) :: fname
! local data
      integer :: lrec, ierr
      if (nrec > 0) return
      inquire(iolength=lrec) f(1)
      lrec = nx*lrec
      if (nrec==0) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old',iostat=ierr)
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace',iostat=ierr)
      endif
      if (ierr==0) nrec = 1
      end subroutine
!
      subroutine fbfvopen(f,nx,ndim,iunit,nrec,fname)
! open direct access binary file for complex 1d vector data
      implicit none
      integer :: nx, ndim, iunit, nrec
      real, dimension(*) :: f
      character(len=*) :: fname
! local data
      integer :: lrec, nnx, ierr
      if (nrec > 0) return
      nnx = ndim*nx
      inquire(iolength=lrec) f(1)
      lrec = nnx*lrec
      if (nrec==0) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old',iostat=ierr)
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace',iostat=ierr)
      endif
      if (ierr==0) nrec = 1
      end subroutine
!
      subroutine fbfcopen(f,nx,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 1d scalar data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
      implicit none
      integer :: nx, iunit, nrec
      complex, dimension(:), pointer :: f
      character(len=*) :: fname
! local data
      integer :: lrec, ierr
      if (nrec > 0) return
      inquire(iolength=lrec) f(1)
      lrec = nx*lrec
      if (nrec==0) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old',iostat=ierr)
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace',iostat=ierr)
      endif
      if (ierr==0) nrec = 1
      end subroutine
!
      subroutine fbfvcopen(f,nx,ndim,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 1d vector data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
      implicit none
      integer :: nx, ndim, iunit, nrec
      complex, dimension(*) :: f
      character(len=*) :: fname
! local data
      integer :: lrec, nnx, ierr
      if (nrec > 0) return
      nnx = ndim*nx
      inquire(iolength=lrec) f(1)
      lrec = nnx*lrec
      if (nrec==0) then
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='old',iostat=ierr)
      else
         open(unit=iunit,file=fname,form='unformatted',access='direct', &
     &recl=lrec,status='replace',iostat=ierr)
      endif
      if (ierr==0) nrec = 1
      end subroutine
!
      subroutine frclose(iunit)
      integer iunit
      close(unit=iunit)
      end subroutine
!
      subroutine frwnml1(fname,name,cn,ddata,sdata,iunit,isign,nl,ns,mcn&
     &,msd,ml,ms,ierr)
! reads and writes namelist and packs and unpacks namelist variables
! into double precision array to passing to C or Java
! input: fname = Fortran file name, name = Fortran namelist name
!    iunit = Fortran unit number, isign = (-1,1) = (read,write) namelist
!    nl = size of variable array, ns = size of string array
!    if isign > 0, the following are also input:
!       cn = names of variable, ddata = values of variables
!       sdata = values of strings
! output: ml = number of variables found, ms = number of strings found
!    mcn = size of cn string, msd = size of sdata string
!    ierr = error if not zero
!    if isign < 0, the following are also output:
!    cn = names of variable, ddata = values of variables
!    sdata = values of strings,
! isign = 0 = return mcn, msd only
      use input1d
      implicit none
      integer :: iunit, isign, nl, ns, mcn, msd, ml, ms, ierr
      character(len=*) :: fname, name
      character(len=16), dimension(nl+ns) :: cn
      double precision, dimension(nl) :: ddata
      character(len=32), dimension(ns) :: sdata
! local data
      integer :: i, j, k, l, ios
      logical :: connected, namef
      character(len=32) :: lname
      character(len=1), parameter :: z = char(0)
      interface rmnullc
         function rmnullc(str) result(sstr)
         implicit none
         character(len=*) :: str
         character(len=len(str)) :: sstr
         end function
      end interface
! set string sizes
      mcn = len(cn); msd = len(sdata)
      ierr = 0
! return variable string size only
      if (isign==0) return
      i = 0; j = 0; k = 0; l = 0
!
      select case(name)
      case ('input1')
!
! read namelist and write data
      if (isign.lt.0) then
! open namelist for input
         open(unit=iunit,file=fname,form='formatted',status='old',      &
     &iostat=ios)
         if (ios /= 0) then
            ierr = -2
            return
         endif
         read (iunit,input1,iostat=ios)
         if (ios /= 0) then
            ierr = -3
            return
         endif
         close(unit=iunit)
! pack data and names
         i = i + 1; j = min(nl,i); ddata(j) = idrun
                                     cn(j) = 'idrun'//z
         i = i + 1; j = min(nl,i); ddata(j) = idrun0
                                     cn(j) = 'idrun0'//z
         i = i + 1; j = min(nl,i); ddata(j) = idcode
                                     cn(j) = 'idcode'//z
         i = i + 1; j = min(nl,i); ddata(j) = indx
                                     cn(j) = 'indx'//z
         i = i + 1; j = min(nl,i); ddata(j) = npx
                                     cn(j) = 'npx'//z
         i = i + 1; j = min(nl,i); ddata(j) = npxb
                                     cn(j) = 'npxb'//z
         i = i + 1; j = min(nl,i); ddata(j) = inorder
                                     cn(j) = 'inorder'//z
         i = i + 1; j = min(nl,i); ddata(j) = popt
                                     cn(j) = 'popt'//z
         i = i + 1; j = min(nl,i); ddata(j) = dopt
                                     cn(j) = 'dopt'//z
         i = i + 1; j = min(nl,i); ddata(j) = djopt
                                     cn(j) = 'djopt'//z
         i = i + 1; j = min(nl,i); ddata(j) = nustrt
                                     cn(j) = 'nustrt'//z
         i = i + 1; j = min(nl,i); ddata(j) = ntr
                                     cn(j) = 'ntr'//z
         i = i + 1; j = min(nl,i); ddata(j) = ntw
                                     cn(j) = 'ntw'//z
         i = i + 1; j = min(nl,i); ddata(j) = ntp
                                     cn(j) = 'ntp'//z
         i = i + 1; j = min(nl,i); ddata(j) = nta
                                     cn(j) = 'nta'//z
         i = i + 1; j = min(nl,i); ddata(j) = ntv
                                     cn(j) = 'ntv'//z
         i = i + 1; j = min(nl,i); ddata(j) = nts
                                     cn(j) = 'nts'//z
         i = i + 1; j = min(nl,i); ddata(j) = nte
                                     cn(j) = 'nte'//z
!        i = i + 1; j = min(nl,i); ddata(j) = ndw
!                                    cn(j) = 'ndw'//z
!        i = i + 1; j = min(nl,i); ddata(j) = ndp
!                                    cn(j) = 'ndp'//z
!        i = i + 1; j = min(nl,i); ddata(j) = nda
!                                    cn(j) = 'nda'//z
!        i = i + 1; j = min(nl,i); ddata(j) = ndv
!                                    cn(j) = 'ndv'//z
!        i = i + 1; j = min(nl,i); ddata(j) = nds
!                                    cn(j) = 'nds'//z
!        i = i + 1; j = min(nl,i); ddata(j) = nde
!                                    cn(j) = 'nde'//z
         i = i + 1; j = min(nl,i); ddata(j) = tend
                                     cn(j) = 'tend'//z
         i = i + 1; j = min(nl,i); ddata(j) = dt
                                     cn(j) = 'dt'//z
         i = i + 1; j = min(nl,i); ddata(j) = qme
                                     cn(j) = 'qme'//z
         i = i + 1; j = min(nl,i); ddata(j) = vtx
                                     cn(j) = 'vtx'//z
         i = i + 1; j = min(nl,i); ddata(j) = vty
                                     cn(j) = 'vty'//z
         i = i + 1; j = min(nl,i); ddata(j) = vtz
                                     cn(j) = 'vtz'//z
         i = i + 1; j = min(nl,i); ddata(j) = vx0
                                     cn(j) = 'vx0'//z
         i = i + 1; j = min(nl,i); ddata(j) = vy0
                                     cn(j) = 'vy0'//z
         i = i + 1; j = min(nl,i); ddata(j) = vz0
                                     cn(j) = 'vz0'//z
         i = i + 1; j = min(nl,i); ddata(j) = vdx
                                     cn(j) = 'vdx'//z
         i = i + 1; j = min(nl,i); ddata(j) = vdy
                                     cn(j) = 'vdy'//z
         i = i + 1; j = min(nl,i); ddata(j) = vdz
                                     cn(j) = 'vdz'//z
         i = i + 1; j = min(nl,i); ddata(j) = vtdx
                                     cn(j) = 'vtdx'//z
         i = i + 1; j = min(nl,i); ddata(j) = vtdy
                                     cn(j) = 'vtdy'//z
         i = i + 1; j = min(nl,i); ddata(j) = vtdz
                                     cn(j) = 'vtdz'//z
         i = i + 1; j = min(nl,i); ddata(j) = psolve
                                     cn(j) = 'psolve'//z
!        i = i + 1; j = min(nl,i); ddata(j) = relativity
!                                    cn(j) = 'relativity'//z
         i = i + 1; j = min(nl,i); ddata(j) = omx
                                     cn(j) = 'omx'//z
         i = i + 1; j = min(nl,i); ddata(j) = omy
                                     cn(j) = 'omy'//z
         i = i + 1; j = min(nl,i); ddata(j) = omz
                                     cn(j) = 'omz'//z
         i = i + 1; j = min(nl,i); ddata(j) = ci
                                     cn(j) = 'ci'//z
         i = i + 1; j = min(nl,i); ddata(j) = ax
                                     cn(j) = 'ax'//z
         i = i + 1; j = min(nl,i); ddata(j) = ndim
                                     cn(j) = 'ndim'//z
         i = i + 1; j = min(nl,i); ddata(j) = ndc
                                     cn(j) = 'ndc'//z
         i = i + 1; j = min(nl,i); ddata(j) = movion
                                     cn(j) = 'movion'//z
         i = i + 1; j = min(nl,i); ddata(j) = sortime
                                     cn(j) = 'sortime'//z
         i = i + 1; j = min(nl,i); ddata(j) = nplot
                                     cn(j) = 'nplot'//z
         i = i + 1; j = min(nl,i); ddata(j) = sntasks
                                     cn(j) = 'sntasks'//z
         i = i + 1; j = min(nl,i); ddata(j) = ndprof
                                     cn(j) = 'ndprof'//z
         i = i + 1; j = min(nl,i); ddata(j) = ampdx
                                     cn(j) = 'ampdx'//z
         i = i + 1; j = min(nl,i); ddata(j) = scaledx
                                     cn(j) = 'scaledx'//z
         i = i + 1; j = min(nl,i); ddata(j) = shiftdx
                                     cn(j) = 'shiftdx'//z
         i = i + 1; j = min(nl,i); ddata(j) = modesxp
                                     cn(j) = 'modesxp'//z
         i = i + 1; j = min(nl,i); ddata(j) = modesxa
                                     cn(j) = 'modesxa'//z
         i = i + 1; j = min(nl,i); ddata(j) = modesxe
                                     cn(j) = 'modesxe'//z
! check for overflow
         if (i.gt.nl) ierr = i
         if (k.gt.ns) ierr = nl + k
! set number of variables found
         ml = j; ms = l;
! write namelist from data
      else if (isign.gt.0) then
! unpack data and check names
         i = min(nl,i+1); idrun = ddata(i)
            if (cn(i) /= 'idrun'//z) ierr = i
         i = min(nl,i+1); idrun0 = ddata(i)
            if (cn(i) /= 'idrun0'//z) ierr = i
         i = min(nl,i+1); idcode = ddata(i)
            if (cn(i) /= 'idcode'//z) ierr = i
         i = min(nl,i+1); indx = ddata(i)
            if (cn(i) /= 'indx'//z) ierr = i
         i = min(nl,i+1); npx = ddata(i)
            if (cn(i) /= 'npx'//z) ierr = i
         i = min(nl,i+1); npxb = ddata(i)
            if (cn(i) /= 'npxb'//z) ierr = i
         i = min(nl,i+1); inorder = ddata(i)
            if (cn(i) /= 'inorder'//z) ierr = i
         i = min(nl,i+1); popt = ddata(i)
            if (cn(i) /= 'popt'//z) ierr = i
         i = min(nl,i+1); dopt = ddata(i)
            if (cn(i) /= 'dopt'//z) ierr = i
         i = min(nl,i+1); djopt = ddata(i)
            if (cn(i) /= 'djopt'//z) ierr = i
         i = min(nl,i+1); nustrt = ddata(i)
            if (cn(i) /= 'nustrt'//z) ierr = i
         i = min(nl,i+1); ntr = ddata(i)
            if (cn(i) /= 'ntr'//z) ierr = i
         i = min(nl,i+1); ntw = ddata(i)
            if (cn(i) /= 'ntw'//z) ierr = i
         i = min(nl,i+1); ntp = ddata(i)
            if (cn(i) /= 'ntp'//z) ierr = i
         i = min(nl,i+1); nta = ddata(i)
            if (cn(i) /= 'nta'//z) ierr = i
         i = min(nl,i+1); ntv = ddata(i)
            if (cn(i) /= 'ntv'//z) ierr = i
         i = min(nl,i+1); nts = ddata(i)
            if (cn(i) /= 'nts'//z) ierr = i
         i = min(nl,i+1); nte = ddata(i)
            if (cn(i) /= 'nte'//z) ierr = i
!        i = min(nl,i+1); ndw = ddata(i)
!           if (cn(i) /= 'ndw'//z) ierr = i
!        i = min(nl,i+1); ndp = ddata(i)
!           if (cn(i) /= 'ndp'//z) ierr = i
!        i = min(nl,i+1); nda = ddata(i)
!           if (cn(i) /= 'nda'//z) ierr = i
!        i = min(nl,i+1); ndv = ddata(i)
!           if (cn(i) /= 'ndv'//z) ierr = i
!        i = min(nl,i+1); nds = ddata(i)
!           if (cn(i) /= 'nds'//z) ierr = i
!        i = min(nl,i+1); nde = ddata(i)
!           if (cn(i) /= 'nde'//z) ierr = i
         i = min(nl,i+1); tend = ddata(i)
            if (cn(i) /= 'tend'//z) ierr = i
         i = min(nl,i+1); dt = ddata(i)
            if (cn(i) /= 'dt'//z) ierr = i
         i = min(nl,i+1); qme = ddata(i)
            if (cn(i) /= 'qme'//z) ierr = i
         i = min(nl,i+1); vtx = ddata(i)
            if (cn(i) /= 'vtx'//z) ierr = i
         i = min(nl,i+1); vty = ddata(i)
            if (cn(i) /= 'vty'//z) ierr = i
         i = min(nl,i+1); vtz = ddata(i)
            if (cn(i) /= 'vtz'//z) ierr = i
         i = min(nl,i+1); vx0 = ddata(i)
            if (cn(i) /= 'vx0'//z) ierr = i
         i = min(nl,i+1); vy0 = ddata(i)
            if (cn(i) /= 'vy0'//z) ierr = i
         i = min(nl,i+1); vz0 = ddata(i)
            if (cn(i) /= 'vz0'//z) ierr = i
         i = min(nl,i+1); vdx = ddata(i)
            if (cn(i) /= 'vdx'//z) ierr = i
         i = min(nl,i+1); vdy = ddata(i)
            if (cn(i) /= 'vdy'//z) ierr = i
         i = min(nl,i+1); vdz = ddata(i)
            if (cn(i) /= 'vdz'//z) ierr = i
         i = min(nl,i+1); vtdx = ddata(i)
            if (cn(i) /= 'vtdx'//z) ierr = i
         i = min(nl,i+1); vtdy = ddata(i)
            if (cn(i) /= 'vtdy'//z) ierr = i
         i = min(nl,i+1); vtdz = ddata(i)
            if (cn(i) /= 'vtdz'//z) ierr = i
         i = min(nl,i+1); psolve = ddata(i)
            if (cn(i) /= 'psolve'//z) ierr = i
!        i = min(nl,i+1); relativity = ddata(i)
!           if (cn(i) /= 'relativity'//z) ierr = i
         i = min(nl,i+1); omx = ddata(i)
            if (cn(i) /= 'omx'//z) ierr = i
         i = min(nl,i+1); omy = ddata(i)
            if (cn(i) /= 'omy'//z) ierr = i
         i = min(nl,i+1); omz = ddata(i)
            if (cn(i) /= 'omz'//z) ierr = i
         i = min(nl,i+1); ci = ddata(i)
            if (cn(i) /= 'ci'//z) ierr = i
         i = min(nl,i+1); ax = ddata(i)
            if (cn(i) /= 'ax'//z) ierr = i
         i = min(nl,i+1); ndim = ddata(i)
            if (cn(i) /= 'ndim'//z) ierr = i
         i = min(nl,i+1); ndc = ddata(i)
            if (cn(i) /= 'ndc'//z) ierr = i
         i = min(nl,i+1); movion = ddata(i)
            if (cn(i) /= 'movion'//z) ierr = i
         i = min(nl,i+1); sortime = ddata(i)
            if (cn(i) /= 'sortime'//z) ierr = i
         i = min(nl,i+1); nplot = ddata(i)
            if (cn(i) /= 'nplot'//z) ierr = i
         i = min(nl,i+1); sntasks = ddata(i)
            if (cn(i) /= 'sntasks'//z) ierr = i
         i = min(nl,i+1); ndprof = ddata(i)
            if (cn(i) /= 'ndprof'//z) ierr = i
         i = min(nl,i+1); ampdx = ddata(i)
            if (cn(i) /= 'ampdx'//z) ierr = i
         i = min(nl,i+1); scaledx = ddata(i)
            if (cn(i) /= 'scaledx'//z) ierr = i
         i = min(nl,i+1); shiftdx = ddata(i)
            if (cn(i) /= 'shiftdx'//z) ierr = i
         i = min(nl,i+1); modesxp = ddata(i)
            if (cn(i) /= 'modesxp'//z) ierr = i
         i = min(nl,i+1); modesxa = ddata(i)
            if (cn(i) /= 'modesxa'//z) ierr = i
         i = min(nl,i+1); modesxe = ddata(i)
            if (cn(i) /= 'modesxe'//z) ierr = i
! open namelist for output
         inquire(unit=iunit,opened=connected,named=namef,name=lname)
         if (connected) then
            if (.not.namef) lname = 'undefined'
         endif
         if ((.not.connected).or.(fname.ne.lname)) then
            open(unit=iunit,file=fname,form='formatted',status='replace'&
     &,iostat=ios)
            if (ios /= 0) then
               ierr = -4
               return
            endif
         endif
         write (iunit,input1,iostat=ios)
         if (ios /= 0) then
            ierr = -5
            return
         endif
! set number of variables found
         ml = i; ms = k
      endif
!
      case ('pot1d')
!
! read namelist and write data
      if (isign.lt.0) then
! open namelist for input
         open(unit=iunit,file=fname,form='formatted',status='old',      &
     &iostat=ios)
         if (ios /= 0) then
            ierr = -6
            return
         endif
         read (iunit,pot1d,iostat=ios)
         if (ios /= 0) then
            ierr = -7
            return
         endif
         close(unit=iunit)
! pack data and names
         i = i + 1; j = min(nl,i); ddata(j) = idrun
                                     cn(j) = 'idrun'//z
         i = i + 1; j = min(nl,i); ddata(j) = indx
                                     cn(j) = 'indx'//z
         i = i + 1; j = min(nl,i); ddata(j) = ntp
                                     cn(j) = 'ntp'//z
         i = i + 1; j = min(nl,i); ddata(j) = modesxp
                                     cn(j) = 'modesxp'//z
         i = i + 1; j = min(nl,i); ddata(j) = psolve
                                     cn(j) = 'psolve'//z
         i = i + 1; j = min(nl,i); ddata(j) = omx
                                     cn(j) = 'omx'//z
         i = i + 1; j = min(nl,i); ddata(j) = omy
                                     cn(j) = 'omy'//z
         i = i + 1; j = min(nl,i); ddata(j) = omz
                                     cn(j) = 'omz'//z
         i = i + 1; j = min(nl,i); ddata(j) = nprec
                                     cn(j) = 'nprec'//z
         i = i + 1; j = min(nl,i); ddata(j) = t0
                                     cn(j) = 't0'//z
         i = i + 1; j = min(nl,i); ddata(j) = tend
                                     cn(j) = 'tend'//z
         i = i + 1; j = min(nl,i); ddata(j) = dt
                                     cn(j) = 'dt'//z
         i = i + 1; j = min(nl,i); ddata(j) = ceng
                                     cn(j) = 'ceng'//z
         k = k + 1; l = min(ns,k); sdata(l) = trim(fpname)//z
                                   cn(j+l) = 'fpname'//z
! check for overflow
         if (i.gt.nl) ierr = i
         if (k.gt.ns) ierr = nl + k
! set number of variables found
         ml = j; ms = l;
! write namelist from data
      else if (isign.gt.0) then
! unpack data and check names
         i = min(nl,i+1); idrun = ddata(i)
            if (cn(i) /= 'idrun'//z) ierr = i
         i = min(nl,i+1); indx = ddata(i)
            if (cn(i) /= 'indx'//z) ierr = i
         i = min(nl,i+1); ntp = ddata(i)
            if (cn(i) /= 'ntp'//z) ierr = i
         i = min(nl,i+1); modesxp = ddata(i)
            if (cn(i) /= 'modesxp'//z) ierr = i
         i = min(nl,i+1); psolve = ddata(i)
            if (cn(i) /= 'psolve'//z) ierr = i
         i = min(nl,i+1); omx = ddata(i)
            if (cn(i) /= 'omx'//z) ierr = i
         i = min(nl,i+1); omy = ddata(i)
            if (cn(i) /= 'omy'//z) ierr = i
         i = min(nl,i+1); omz = ddata(i)
            if (cn(i) /= 'omz'//z) ierr = i
         i = min(nl,i+1); nprec = ddata(i)
            if (cn(i) /= 'nprec'//z) ierr = i
         i = min(nl,i+1); t0 = ddata(i)
            if (cn(i) /= 't0'//z) ierr = i
         i = min(nl,i+1); tend = ddata(i)
            if (cn(i) /= 'tend'//z) ierr = i
         i = min(nl,i+1); dt = ddata(i)
            if (cn(i) /= 'dt'//z) ierr = i
         i = min(nl,i+1); ceng = ddata(i)
            if (cn(i) /= 'ceng'//z) ierr = i
         k = min(ns,k+1); fpname = rmnullc(sdata(k))
            if (cn(i+k) /= 'fpname'//z) ierr = i + k
! open namelist for output
         inquire(unit=iunit,opened=connected,named=namef,name=lname)
         if (connected) then
            if (.not.namef) lname = 'undefined'
         endif
         if ((.not.connected).or.(fname.ne.lname)) then
            open(unit=iunit,file=fname,form='formatted',status='replace'&
     &,iostat=ios)
            if (ios /= 0) then
               ierr = -8
               return
            endif
         endif
         write (iunit,pot1d,iostat=ios)
         if (ios /= 0) then
            ierr = -9
            return
         endif
! set number of variables found
         ml = i; ms = k
      endif
!
      case ('vpot1d')
!
! read namelist and write data
      if (isign.lt.0) then
! open namelist for input
         open(unit=iunit,file=fname,form='formatted',status='old',      &
     &iostat=ios)
         if (ios /= 0) then
            ierr = -10
            return
         endif
         read (iunit,vpot1d,iostat=ios)
         if (ios /= 0) then
            ierr = -11
            return
         endif
         close(unit=iunit)
! pack data and names
         i = i + 1; j = min(nl,i); ddata(j) = idrun
                                     cn(j) = 'idrun'//z
         i = i + 1; j = min(nl,i); ddata(j) = indx
                                     cn(j) = 'indx'//z
         i = i + 1; j = min(nl,i); ddata(j) = nta
                                     cn(j) = 'nta'//z
         i = i + 1; j = min(nl,i); ddata(j) = modesxa
                                     cn(j) = 'modesxa'//z
         i = i + 1; j = min(nl,i); ddata(j) = psolve
                                     cn(j) = 'psolve'//z
         i = i + 1; j = min(nl,i); ddata(j) = ndim
                                     cn(j) = 'ndim'//z
         i = i + 1; j = min(nl,i); ddata(j) = omx
                                     cn(j) = 'omx'//z
         i = i + 1; j = min(nl,i); ddata(j) = omy
                                     cn(j) = 'omy'//z
         i = i + 1; j = min(nl,i); ddata(j) = omz
                                     cn(j) = 'omz'//z
         i = i + 1; j = min(nl,i); ddata(j) = ci
                                     cn(j) = 'ci'//z
         i = i + 1; j = min(nl,i); ddata(j) = narec
                                     cn(j) = 'narec'//z
         i = i + 1; j = min(nl,i); ddata(j) = t0
                                     cn(j) = 't0'//z
         i = i + 1; j = min(nl,i); ddata(j) = tend
                                     cn(j) = 'tend'//z
         i = i + 1; j = min(nl,i); ddata(j) = dt
                                     cn(j) = 'dt'//z
         i = i + 1; j = min(nl,i); ddata(j) = ceng
                                     cn(j) = 'ceng'//z
         k = k + 1; l = min(ns,k); sdata(l) = trim(faname)//z
                                   cn(j+l) = 'faname'//z
! check for overflow
         if (i.gt.nl) ierr = i
         if (k.gt.ns) ierr = nl + k
! set number of variables found
         ml = j; ms = l;
! write namelist from data
      else if (isign.gt.0) then
! unpack data and check names
         i = min(nl,i+1); idrun = ddata(i)
            if (cn(i) /= 'idrun'//z) ierr = i
         i = min(nl,i+1); indx = ddata(i)
            if (cn(i) /= 'indx'//z) ierr = i
         i = min(nl,i+1); nta = ddata(i)
            if (cn(i) /= 'nta'//z) ierr = i
         i = min(nl,i+1); modesxa = ddata(i)
            if (cn(i) /= 'modesxa'//z) ierr = i
         i = min(nl,i+1); psolve = ddata(i)
            if (cn(i) /= 'psolve'//z) ierr = i
         i = min(nl,i+1); ndim = ddata(i)
            if (cn(i) /= 'ndim'//z) ierr = i
         i = min(nl,i+1); omx = ddata(i)
            if (cn(i) /= 'omx'//z) ierr = i
         i = min(nl,i+1); omy = ddata(i)
            if (cn(i) /= 'omy'//z) ierr = i
         i = min(nl,i+1); omz = ddata(i)
            if (cn(i) /= 'omz'//z) ierr = i
         i = min(nl,i+1); ci = ddata(i)
            if (cn(i) /= 'ci'//z) ierr = i
         i = min(nl,i+1); narec = ddata(i)
            if (cn(i) /= 'narec'//z) ierr = i
         i = min(nl,i+1); t0 = ddata(i)
            if (cn(i) /= 't0'//z) ierr = i
         i = min(nl,i+1); tend = ddata(i)
            if (cn(i) /= 'tend'//z) ierr = i
         i = min(nl,i+1); dt = ddata(i)
            if (cn(i) /= 'dt'//z) ierr = i
         i = min(nl,i+1); ceng = ddata(i)
            if (cn(i) /= 'ceng'//z) ierr = i
         k = min(ns,k+1); faname = rmnullc(sdata(k))
            if (cn(i+k) /= 'faname'//z) ierr = i + k
! open namelist for output
         inquire(unit=iunit,opened=connected,named=namef,name=lname)
         if (connected) then
            if (.not.namef) lname = 'undefined'
         endif
         if ((.not.connected).or.(fname.ne.lname)) then
            open(unit=iunit,file=fname,form='formatted',status='replace'&
     &,iostat=ios)
            if (ios /= 0) then
               ierr = -12
               return
            endif
         endif
         write (iunit,vpot1d,iostat=ios)
         if (ios /= 0) then
            ierr = -13
            return
         endif
! set number of variables found
         ml = i; ms = k
      endif
!
      case ('em1d')
!
! read namelist and write data
      if (isign.lt.0) then
! open namelist for input
         open(unit=iunit,file=fname,form='formatted',status='old',      &
     &iostat=ios)
         if (ios /= 0) then
            ierr = -14
            return
         endif
         read (iunit,em1d,iostat=ios)
         if (ios /= 0) then
            ierr = -15
            return
         endif
         close(unit=iunit)
! pack data and names
         i = i + 1; j = min(nl,i); ddata(j) = idrun
                                     cn(j) = 'idrun'//z
         i = i + 1; j = min(nl,i); ddata(j) = indx
                                     cn(j) = 'indx'//z
         i = i + 1; j = min(nl,i); ddata(j) = nte
                                     cn(j) = 'nte'//z
         i = i + 1; j = min(nl,i); ddata(j) = modesxe
                                     cn(j) = 'modesxe'//z
         i = i + 1; j = min(nl,i); ddata(j) = psolve
                                     cn(j) = 'psolve'//z
         i = i + 1; j = min(nl,i); ddata(j) = ndim
                                     cn(j) = 'ndim'//z
         i = i + 1; j = min(nl,i); ddata(j) = omx
                                     cn(j) = 'omx'//z
         i = i + 1; j = min(nl,i); ddata(j) = omy
                                     cn(j) = 'omy'//z
         i = i + 1; j = min(nl,i); ddata(j) = omz
                                     cn(j) = 'omz'//z
         i = i + 1; j = min(nl,i); ddata(j) = ci
                                     cn(j) = 'ci'//z
         i = i + 1; j = min(nl,i); ddata(j) = nerec
                                     cn(j) = 'nerec'//z
         i = i + 1; j = min(nl,i); ddata(j) = t0
                                     cn(j) = 't0'//z
         i = i + 1; j = min(nl,i); ddata(j) = tend
                                     cn(j) = 'tend'//z
         i = i + 1; j = min(nl,i); ddata(j) = dt
                                     cn(j) = 'dt'//z
         i = i + 1; j = min(nl,i); ddata(j) = ceng
                                     cn(j) = 'ceng'//z
         k = k + 1; l = min(ns,k); sdata(l) = trim(fename)//z
                                   cn(j+l) = 'fename'//z
! check for overflow
         if (i.gt.nl) ierr = i
         if (k.gt.ns) ierr = nl + k
! set number of variables found
         ml = j; ms = l;
! write namelist from data
      else if (isign.gt.0) then
! unpack data and check names
         i = min(nl,i+1); idrun = ddata(i)
            if (cn(i) /= 'idrun'//z) ierr = i
         i = min(nl,i+1); indx = ddata(i)
            if (cn(i) /= 'indx'//z) ierr = i
         i = min(nl,i+1); nte = ddata(i)
            if (cn(i) /= 'nte'//z) ierr = i
         i = min(nl,i+1); modesxe = ddata(i)
            if (cn(i) /= 'modesxe'//z) ierr = i
         i = min(nl,i+1); psolve = ddata(i)
            if (cn(i) /= 'psolve'//z) ierr = i
         i = min(nl,i+1); ndim = ddata(i)
            if (cn(i) /= 'ndim'//z) ierr = i
         i = min(nl,i+1); omx = ddata(i)
            if (cn(i) /= 'omx'//z) ierr = i
         i = min(nl,i+1); omy = ddata(i)
            if (cn(i) /= 'omy'//z) ierr = i
         i = min(nl,i+1); omz = ddata(i)
            if (cn(i) /= 'omz'//z) ierr = i
         i = min(nl,i+1); ci = ddata(i)
            if (cn(i) /= 'ci'//z) ierr = i
         i = min(nl,i+1); nerec = ddata(i)
            if (cn(i) /= 'nerec'//z) ierr = i
         i = min(nl,i+1); t0 = ddata(i)
            if (cn(i) /= 't0'//z) ierr = i
         i = min(nl,i+1); tend = ddata(i)
            if (cn(i) /= 'tend'//z) ierr = i
         i = min(nl,i+1); dt = ddata(i)
            if (cn(i) /= 'dt'//z) ierr = i
         i = min(nl,i+1); ceng = ddata(i)
            if (cn(i) /= 'ceng'//z) ierr = i
         k = min(ns,k+1); fename = rmnullc(sdata(k))
            if (cn(i+k) /= 'fename'//z) ierr = i + k
! open namelist for output
         inquire(unit=iunit,opened=connected,named=namef,name=lname)
         if (connected) then
            if (.not.namef) lname = 'undefined'
         endif
         if ((.not.connected).or.(fname.ne.lname)) then
            open(unit=iunit,file=fname,form='formatted',status='replace'&
     &,iostat=ios)
            if (ios /= 0) then
               ierr = -16
               return
            endif
         endif
         write (iunit,em1d,iostat=ios)
         if (ios /= 0) then
            ierr = -17
            return
         endif
! set number of variables found
         ml = i; ms = k
      endif
!
      case default
         ierr = -18
      end select
      end subroutine
!
      function rmnullc(str) result(sstr)
! this function replaces null characters with blanks
      implicit none
      character(len=*) :: str
! local data
      integer :: i
      character(len=len(str)) :: sstr
      sstr = str
      do i = 1, len(str)
      if (ichar(str(i:i))==0) sstr(i:i) = ' '
      enddo
      end function
