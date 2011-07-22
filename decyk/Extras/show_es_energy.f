!-----------------------------------------------------------------------
      program show_es_energy
      use graf1d_class
      implicit none
      integer :: nloop, irc, itime, l, nplot, ntw, it
      real :: we, wke, wt, time, dt
      character(len=60) :: prompt
      character(len=64) :: label
      character(len=128), dimension(3) :: labels
      character(len=32) :: fname, ftname
      character(len=5) :: line1
      character(len=34) :: line2
      character(len=10), dimension(3) :: chrs
      real, dimension(:,:), allocatable :: w
      type (graf1d) :: single, total
  991 format (5h t = ,i7)
  992 format (19h * * * q.e.d. * * *)
  993 format (34h field, kinetic, total energies = ,3e14.7)
! open graphics device
      call GROPEN0(1)
! get file name
      label = ' '
    5 prompt = 'enter name of em energy file, q to quit:'
      call GTINPUT(label,prompt,fname,irc)
      if ((fname=='q').or.(fname=='Q')) go to 210
! clear screen
      call CLRSCRN
! open es energy file
      open(unit=18,file=trim(fname),form='formatted',status='old',iostat&
     &=irc)
      if (irc /= 0) then
         label = 'open error for '//trim(fname)
         go to 5
      endif
      read (18,*) line1
      nloop = 0
      do
         read (18,'(a5,i7)',iostat=irc) line1, itime
         if (irc /= 0) exit
         read (18,'(a34,3e14.7)',iostat=irc) line2, we, wke, wt
         if (irc /= 0) exit
         nloop = nloop + 1
      end do
      if (nloop==0) go to 5
      prompt = ' '
      write (label,*) 'finished reading file, nloop = ', nloop
! get file name
    6 prompt = 'enter ntw and dt to write energies (0 to omit):'
      call GTINPUT(label,prompt,ftname,irc)
      if ((ftname=='q').or.(ftname=='Q')) go to 210
      if (irc==1) then
         label = ' '
         go to 6
      endif
      read (ftname,*,iostat=irc) ntw, dt
      if (irc /= 0) ntw = 0
      allocate(w(nloop,3))
      time = real(nloop-1)
      rewind 18
      read (18,*) line1
      if (ntw > 0) then
         fname = 'eng'
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
         write (19,*) 'time wel ke wtot dke'
         fname = 'leng'
         open(unit=20,file=trim(fname),form='formatted',status='replace'&
     &)
         write (20,*) 'time logwel logke logwtot'
      endif
      do l = 1, nloop
         read (18,'(a5,i7)') line1, itime
         read (18,'(a34,3e14.7)') line2, w(l,1), w(l,2), w(l,3)
! energy diagnostic
         if (ntw > 0) then
            it = itime/ntw
            if (itime==ntw*it) then
               time = dt*real(itime)
               wt = w(l,3) - w(1,3)
               we = w(l,2) - w(1,2)
               write (19,*) time, w(l,1), w(l,2), wt, we
               write (20,*) time, log(w(l,1)), log(w(l,2))
            endif
         endif
      end do
      we = minval(w(:,3)); wt = maxval(w(:,3)); wke = wt - we
      prompt = ' Hit carrriage return or enter key to continue '
      write (labels(1),*) 'delta/init energy=', wke, w(1,3)
      we = sum(w(:,1))/real(nloop)
      write (labels(2),*) 'average el energy=', we
      wke = sum(w(:,2))/real(nloop)
      write (labels(3),*) 'average eke energy=', wke
      call GTMINPUT(labels,prompt,fname,3,irc)
      if (irc==1) go to 20
      nplot = 4
! clear screen
      call CLRSCRN
! set number of plots per page
      call SETNPLT(nplot,irc)
      call new_graf1d(single,nloop,0.,time,'ENERGY',clip=2)
      chrs(1) = 'T  FIELD  '
      chrs(2) = 'KINETIC   '
      chrs(3) = 'TOTAL     '
      call new_graf1d(total,nloop,0.,time,'ENERGY',nsubs=3,clip=2,lsubs=&
     &chrs(1:3))
   10 call display(single,w(:,1:1),irc,label2='TOTAL FIELD')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,2:2),irc,label2='KINETIC ENERGY')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,3:3),irc,label2='TOTAL ENERGY')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(total,w(:,1:3),irc,label2='TOTAL')
      if (irc==6) go to 10
   20 call del_graf1d(single)
      call del_graf1d(total)
! close graphics device
  210 call GRCLOSE
      end program
