!-----------------------------------------------------------------------
      program showmomentum
      use graf1d_class
      implicit none
      integer :: nloop, irc, itime, l, nplot
      real :: px, py, pz, sx, sy, sz, pxi, pyi, pzi, tx, ty, tz, time
      character(len=60) :: prompt
      character(len=64) :: label
      character(len=128), dimension(3) :: labels
      character(len=32) :: fname
      character(len=5) :: line1
      character(len=21) :: line2
      character(len=18) :: line3
      character(len=16) :: line4
      character(len=10), dimension(12) :: chrs
      real, dimension(:,:), allocatable :: w
      type (graf1d) :: single, particle, field, total
  991 format (5h t = ,i7)
  994 format (21h electron momentum = ,3e14.7)
  995 format (16h ion momentum = ,3e14.7)
  996 format (18h field momentum = ,3e14.7)
  997 format (18h total momentum = ,3e14.7)
! open graphics device
      call GROPEN0(1)
! get file name
      label = ' '
    5 prompt = 'enter name of momentum file, q to quit:'
      call GTINPUT(label,prompt,fname,irc)
      if ((fname=='q').or.(fname=='Q')) go to 210
! clear screen
      call CLRSCRN
! open momentum file
      open(unit=19,file=trim(fname),form='formatted',status='old',iostat&
     &=irc)
      if (irc /= 0) then
         label = 'open error for '//trim(fname)
         go to 5
      endif
      nloop = 0
      do
         read (19,'(a5,i7)',iostat=irc) line1, itime
         if (irc /= 0) exit
         read (19,'(a21,3e14.7)',iostat=irc) line2, px, py, pz
         if (irc /= 0) exit
         read (19,'(a18,3e14.7)',iostat=irc) line3, sx, sy, sz
         if (irc /= 0) exit
         read (19,'(a16,3e14.7)',iostat=irc) line4, pxi, pyi, pzi
         if (irc /= 0) exit
         read (19,'(a18,3e14.7)',iostat=irc) line3, tx, ty, tz
         if (irc /= 0) exit
         nloop = nloop + 1
      end do
      if (nloop==0) go to 5
      prompt = ' '
      write (label,*) 'finished reading file, nloop = ', nloop
      call PTOTPUT(label,prompt)
      allocate(w(nloop,12))
      time = real(nloop-1)
      rewind 19
      do l = 1, nloop
         read (19,'(a5,i7)') line1, itime
         read (19,'(a21,3e14.7)') line2, w(l,1), w(l,2), w(l,3)
         read (19,'(a18,3e14.7)') line3, w(l,4), w(l,5), w(l,6)
         read (19,'(a16,3e14.7)') line4, w(l,7), w(l,8), w(l,9)
         read (19,'(a18,3e14.7)') line3, w(l,10), w(l,11), w(l,12)
      end do
      prompt = ' Hit carrriage return or enter key to continue '
      tx = minval(w(:,10)); ty = maxval(w(:,10)); tz = ty - tx
      write (labels(1),*) 'delta/init x momentum=', tz, w(1,10)
      tx = minval(w(:,11)); ty = maxval(w(:,11)); tz = ty - tx
      write (labels(2),*) 'delta/init y momentum=', tz, w(1,11)
      tx = minval(w(:,12)); ty = maxval(w(:,12)); tz = ty - tx
      write (labels(3),*) 'delta/init z momentum=', tz, w(1,12)
      call GTMINPUT(labels,prompt,fname,3,irc)
      if (irc==1) go to 20
      nplot = 4
! clear screen
      call CLRSCRN
! set number of plots per page
      call SETNPLT(nplot,irc)
      call new_graf1d(single,nloop,0.,time,'MOMENTUM',clip=2)
      chrs(1) = 'ELECTRON X'
      chrs(2) = 'ELECTRON Y'
      chrs(3) = 'ELECTRON Z'
      call new_graf1d(particle,nloop,0.,time,'MOMENTUM',nsubs=3,clip=2,l&
     &subs=chrs(1:3))
      chrs(4) = 'FIELD X   '
      chrs(5) = 'FIELD Y   '
      chrs(6) = 'FIELD Z   '
      call new_graf1d(single,nloop,0.,time,'MOMENTUM',clip=2)
      chrs(7) = 'ION X     '
      chrs(8) = 'ION Y     '
      chrs(9) = 'ION Z     '
      call new_graf1d(field,nloop,0.,time,'MOMENTUM',nsubs=3,clip=2,lsub&
     &s=chrs(4:6))
      chrs(10) = 'TOTAL X   '
      chrs(11) = 'TOTAL Y   '
      chrs(12) = 'TOTAL Z   '
      call new_graf1d(total,nloop,0.,time,'MOMENTUM',nsubs=3,clip=2,lsub&
     &s=chrs(7:9))
   10 call display(single,w(:,1:1),irc,label2='ELECTRON X MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,2:2),irc,label2='ELECTRON Y MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,3:3),irc,label2='ELECTRON Z MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(particle,w(:,1:3),irc,label2='ELECTRON MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,4:4),irc,label2='FIELD X MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,5:5),irc,label2='FIELD Y MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,6:6),irc,label2='FIELD Z MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(field,w(:,4:6),irc,label2='FIELD MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,7:7),irc,label2='ION X MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,8:8),irc,label2='ION Y MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,9:9),irc,label2='ION Z MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(particle,w(:,7:9),irc,label2='ION MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,10:10),irc,label2='TOTAL X MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,11:11),irc,label2='TOTAL Y MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(single,w(:,12:12),irc,label2='TOTAL Z MOMENTUM')
      if (irc==6) go to 10
      if (irc==1) go to 20
      call display(total,w(:,10:12),irc,label2='TOTAL MOMENTUM')
      if (irc==6) go to 10
   20 call del_graf1d(single)
      call del_graf1d(particle)
      call del_graf1d(field)
      call del_graf1d(total)
! close graphics device
  210 call GRCLOSE
      end program
