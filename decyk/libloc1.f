      subroutine starch
c this subroutine sets machine architecture description
c these settings are for apple macintosh
      common /march/ iebc, irvb, longi, mtype
      save /march/
c iebc = (0,1) = (no,yes) characters are in ebcdic
      iebc = 0
c irvb = (0,1) = (no,yes) integers are stored in reverse order
      irvb = 0
c longi = (0,1) = (no,yes) integers are 64 bits
      longi = 0
c mtype = machine type
c 1=rs/6000, 2=sun, 3=cray c90, 4=ibm es/9000, 5=vax, 6=dec, 7=hp, 8=sgi
c 9=ibm pc, 10=mac, 11=paragon, 12=cray t3d, 13=fujitsu vpp500
      mtype = 1
      return
      end
c     subroutine TIMERA(icntrl,chr,time)
c this subroutine performs timing
c input: icntrl, chr 
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c chr = character variable for labeling timings
c time = elapsed time in seconds
c written for macintosh
c     character*8 chr
c     integer TickCount
c     external TickCount
c     save jclock
c  91 format (1x,a8,1x,'wall clock time = ',e14.7,1x,'sec')
c     data jclock /0/
c     if (icntrl.eq.0) return
c     if (icntrl.eq.1) go to 10
c initialize clock
c get current system tick count
c     jclock = TickCount()
c     return
c read clock and write time difference from last clock initialization
c get current system tick count
c  10 nclock = TickCount()
c     time = real(nclock - jclock)/60.
c     write (6,91) chr, time
c     return
c     end
