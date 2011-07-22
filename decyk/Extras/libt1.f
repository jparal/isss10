c interactive plot10 driver package for tektronix 4014
c viktor k. decyk, ucla
c copyright 1990, regents of the university of california
c update: february 15, 1996
      subroutine gitwks(idcon,iwtype)
c this is a site-dependent subroutine which returns
c connection identifier and workstation type
c version for plot10 tektronix graphics
c iwtype = workstation type
c 4010 = Tektronix 4010, 4014 = Tektronix 4014, 4016 = Tektronix 4014EGM
c 4105 = Tektronix 4105, 4207 = Tektronix 4207
      character*1 c
      dimension itw(5)
   91 format (' enter tektronix type: (1=4010,2=4014,3=4014egm,4=4105,5=
     14207,q=quit,?=help)')
   92 format (a1)
c itw = workstation type conversion table
      data itw /4010,4014,4016,4105,4207/
c special case for vax
c     open(unit=6,file='sys$output',status='old',carriagecontrol='none')
c write prompt
   10 write (6,91)
      read (5,92,end=20) c
c help requested
      if (c.ne.'?') go to 20
      write (6,*) ' Tektronix terminals supported'
      write (6,*) ' 1 = 4010 = 1024x780 screen, monochrome'
      write (6,*) ' 2 = 4014 = 1024x780 screen, monochrome'
      write (6,*) ' 3 = 4014egm = 4096x3120 screen, monochrome'
      write (6,*) ' 4 = 4105 = 1024x780 screen, 8 colors'
      write (6,*) ' 5 = 4207 = 1024x780 screen, 8 colors'
      go to 10
c convert to integer
   20 iterm = ichar(c) - ichar('0')
c iterm = 0 means abort
      if ((c.eq.'q').or.(c.eq.'Q').or.(iterm.eq.0)) stop 1
c request again if terminal type is invalid
      if ((iterm.lt.0).or.(iterm.gt.5)) go to 10
c convert iterm to workstation type
      iwtype = itw(iterm)
c idcon = connection identifier, 1 seems to work
      idcon = 1
      return
      end
      subroutine tputc(lout,len)
c this program sends characters stored as integers to terminal
c using fortran io
c written for the ibm rs/6000
c input: all
c lout = buffer of particles packed as integers
c len = number of characters to be sent
c lw = number of bytes per word
      parameter(lw=4)
c iebc = (0,1) = (no,yes) output characters are in ebcdic
c mtype = machine type
c 1=rs/6000, 2=sun, 3=cray c90, 4=ibm es/9000, 5=vax, 6=dec, 7=hp, 8=sgi
c 9=ibm pc, 10=mac, 11=paragon, 12=cray t3d, 13=fujitsu vpp500
      common /march/ iebc, irvb, longi, mtype
      dimension lout(*)
c cs = format string
      character*9 cs
c c = output character array
      character*1 c(240)
c lb = integer array for extracting packed characters from integers
      dimension lb(lw)
c iate = ascii to ebcdic conversion table
      dimension iate(128)
      save iate
c ebcdic code for ascii 124 is non-standard
      data iate /0,1,2,3,55,45,46,47,22,5,37,11,12,13,14,15,16,17,18,19,
     160,61,50,38,24,25,63,39,34,29,53,31,64,90,127,123,91,108,80,125,77
     2,93,92,78,107,96,75,97,240,241,242,243,244,245,246,247,248,249,122
     3,94,76,126,110,111,124,193,194,195,196,197,198,199,200,201,209,210
     4,211,212,213,214,215,216,217,226,227,228,229,230,231,232,233,173,2
     524,189,95,109,121,129,130,131,132,133,134,135,136,137,145,146,147,
     6148,149,150,151,152,153,162,163,164,165,166,167,168,169,192,79,208
     7,161,7/
c     if (len.gt.0) call itputg(lout,len,irc)
c return if no characters
      if (len.lt.1) return
c use default format string
      if (len.eq.240) then
         cs = '(240a1,$)'
c find format for exact number of characters to write
      else
         cs = '(   a1,$)'
         iz = ichar('0')
c hundreds place
         ih = len/100
         is = len - 100*ih
c tens place
         id = is/10
c ones place
         is = is - 10*id
c fix format string
         if (ih.gt.0) then
            cs(2:2) = char(ih+iz)
            cs(3:3) = char(id+iz)
         elseif (id.gt.0) then
            cs(3:3) = char(id+iz)
         endif
         cs(4:4) = char(is+iz)
      endif
c l = number of integer words to be processed
      l = (len - 1)/lw + 1
      lw1 = lw - 1
      lwp = lw + 1
c copy characters from integer array to character array
      do 40 i = 1, l
      i1 = lw*(i - 1)
      lb(1) = lout(i)
c loop for extracting packed characters from integer
      do 10 j = 1, lw1
      lb(j+1) = lb(j)/256
      lb(j) = lb(j) - lb(j+1)*256
   10 continue
c copy to character array
      if (iebc.eq.0) then
         do 20 j = 1, lw
         c(i1+j) = char(lb(lwp-j))
   20    continue
c translate to ebcdic
      else
         do 30 j = 1, lw
         c(i1+j) = char(iate(lb(lwp-j)+1))
   30    continue
      endif
   40 continue
c write out character array
      write (6,cs) (c(j),j=1,len)
c ibm rs/6000 and sgi need to rewind 6
      if ((mtype.eq.1).or.(mtype.eq.8)) rewind 6
      return
      end
