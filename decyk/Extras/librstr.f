      subroutine gitwks(idcon,iwtype)
c this is a site-dependent subroutine which returns
c connection identifier and workstation type
c version for nersc raster library
c iwtype = workstation type
c 1 = cga, 2 = ega, 3 = vga, 4 = mac1, 5 = mac2, 6 = mac3, 7 = mac4
c 8 = ascii, 9 = tif1, 10 = tif2, 11 = tif3, 12 = tif4
      character*1 c
   91 format (' enter format: (1=cga,2=ega,3=vga,4=mac1,5=mac2,6=mac3,7=
     1mac4,q=quit,?=help)')
   92 format (a1)
c write prompt
   10 write (6,91)
      read (5,92,end=20) c
c help requested
      if (c.ne.'?') go to 20
      write (6,*) ' Compressed raster image file with NERSC compression'
      write (6,*) ' 1 = cga = 320x200 pixels, 4 colors'
      write (6,*) ' 2 = ega = 640x350 pixels, monochrome'
      write (6,*) ' 3 = vga = 320x200 pixels, 256 colors'
      write (6,*) ' 4 = mac1 = 320x200 pixels, 256 colors'
      write (6,*) ' 5 = mac2 = 320x320 pixels, 256 colors'
      write (6,*) ' 6 = mac3 = 400x400 pixels, 256 colors'
      write (6,*) ' 7 = mac4 = 640x400 pixels, 256 colors'
      write (6,*) ' Plus five additional formats'
      write (6,*) ' 8 = ascii = 78x22 pixels, 64 printable characters'
      write (6,*) ' 9 = tif1 = 320x200 pixels, monochrome TIFF'
      write (6,*) ' a = tif2 = 320x200 pixels, 256 color TIFF'
      write (6,*) ' b = tif3 = 640x400 pixels, monochrome TIFF'
      write (6,*) ' c = tif4 = 640x400 pixels, 256 color TIFF'
      go to 10
c convert to integer
   20 id = ichar(c) - ichar('0')
c id = 0 means abort
      if ((c.eq.'q').or.(c.eq.'Q').or.(id.eq.0)) stop 1
      if ((c.eq.'a').or.(c.eq.'A')) id = 10
      if ((c.eq.'b').or.(c.eq.'B')) id = 11
      if ((c.eq.'c').or.(c.eq.'C')) id = 12
c request again if format type is invalid
      if ((id.lt.0).or.(id.gt.12)) go to 10
c open file for raster images
      if (id.ne.8) then
         open(unit=20,file='mgraph',form='formatted',status='unknown')
c special case for macintosh
c        open(unit=20,file='mgraph',form='unformatted',status='unknown')
c special case for vax
c        open(unit=20,file='mgraph',form='formatted',status='unknown',ca
c    1rriagecontrol='none')
      endif
c return workstation type
      iwtype = id
c idcon = connection identifier, 20 seems to work
      idcon = 20
      return
      end
c gks device driver for nersc raster library
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c version for ibm rs/6000
c update: august 14, 1999
      block data
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      save /gksrstr/
      data kosv /0/
      end
      subroutine gqops(istat)
c inquire operating state value
c input arguments: none
c istat = operating state (0=gks closed,1=gks open,2=workstation open,
c 3=workstation active,4=segment open)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      istat = kosv
      return
      end
      subroutine gopks(nerrfl,meml)
c open gks
c input arguments: all
c nerrfl = error file unit number, 6 for terminal
c meml = storage limit for one segment, in bytes
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 1
      return
      end
      subroutine gopwk(idwk,idcon,iwtype)
c open workstation
c input arguments: all
c idwk = workstation identifier
c idcon = connection identifier
c iwtype = workstation type
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c iwt = current workstation type
      kosv = 2
      iwk = idwk
      idc = idcon
      iwt = iwtype
      return
      end
      subroutine gqopwk(n,ierr,now,idwk)
c inquire set number of open workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c now = number of open workstations
c idwk = workstation identifier of the nth member
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      now = 0
      ierr = 0
      if (kosv.ge.2) now = 1
      if ((n.lt.0).or.(n.gt.1)) ierr = -12
      if (n.eq.1) idwk = iwk
      return
      end
      subroutine gqwkca(iwtype,ierr,iwkca)
c inquire workstation category
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator
c iwkca = workstation category
c (0 = output, 1 = input, 2 = outin, 3 = wiss, 4 = mo, 5 = mi)
      ierr = 0
      iwkca = 0
      if (iwtype.eq.8) iwkca = 2
      return
      end
      subroutine gscnid(idcon,connam)
c set connection identifier
c input arguments: all
c idcon = connection identifier
c connam = connection name
      character*8 connam
c     open(unit=idcon,file=connam,form='formatted',status='unknown')
      return
      end
      subroutine gqwkc(idwk,ierr,idcon,iwtype)
c inquire workstation connection and type
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c idcon = connection identifier
c iwtype = workstation type
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
c idc = current connection identifier
c iwt = current workstation type
      ierr = 0
      idcon = idc
      if (idwk.eq.iwk) then
         iwtype = iwt
      else
         ierr = 20
      endif
      return
      end
      subroutine gacwk(idwk)
c activate workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwt = current workstation type
c initialize nersc raster device
      call gminit(iwt)
c create raster defaults
      call rsdflts
      kosv = 3
      return
      end
      subroutine gqacwk(n,ierr,naw,idwk)
c inquire set of active workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c naw = number of active workstations
c idwk = workstation identifier of the nth member
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      naw = 0
      ierr = 0
      if (kosv.ge.3) naw = 1
      if ((n.lt.0).or.(n.gt.1)) ierr = -12
      if (n.eq.1) idwk = iwk
      return
      end
      subroutine gqwks(idwk,ierr,istat)
c inquire workstation state
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c istat = workstation state (0 = inactive, 1 = active)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      ierr = 0
      istat = 0
      if (kosv.lt.2) ierr = 25
      if (kosv.ge.3) istat = 1
      if (idwk.ne.iwk) then
         istat = 0
         ierr = 20
      endif
      return
      end
      subroutine gqcf(iwtype,ierr,ncoli,iscol,npci)
c inquire color facilities
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator (0=inquiry successful)
c ncoli = number of colors available
c iscol = color availability indicator, 0 = monochrome, 1 = color
c npci = number of predefined color indices
      ierr = 0
c cga format
      if (iwtype.eq.1) then
         ncoli = 4
         iscol = 1
         npci = 4
c ega format
      elseif (iwtype.eq.2) then
         ncoli = 2
         iscol = 0
         npci = 2
c printer plots
      elseif (iwtype.eq.8) then
         ncoli = 64
         iscol = 1
         npci = 64
c monochrome tiff format
      elseif ((iwtype.eq.9).or.(iwtype.eq.11)) then
         ncoli = 2
         iscol = 0
         npci = 2
c 256 color formats
      else
         ncoli = 256
         iscol = 1
         npci = 256
      endif
      return
      end
      subroutine gscr(idwk,ic,cr,cg,cb)
c set color representation
c input arguments: all
c idwk = workstation identifier
c ic = color index
c cr/cg/cb = red/green/blue component (0 < cr,cg,cb < 1)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c icc = color index conversion table
      dimension reds(4), greens(4), blues(4)
      data reds /0.,0.,1.,1./
      data greens /0.,0.,0.,1./
      data blues /0.,1.,0.,1./
c return if index is out of range
      if ((ic.lt.0).or.(ic.gt.255)) return
c cga display
      if (iwt.eq.1) then
c find which cga color index most closely matches gks index
c and set corresponding element of color index conversion table
         imin = 1
         cdm = cr*cr + cg*cg + cb*cb
         do 10 i = 1, 4
         cd = (cr - reds(i))**2 + (cg - greens(i))**2 + (cb - blues(i))*
     1*2
         if (cd.lt.cdm) then
            imin = i
            cdm = cd
         endif
   10    continue
         if ((ic.ge.0).and.(ic.le.3)) icc(ic+1) = imin - 1
c other displays
      endif
c convert and store palette entry
      call swpal(ic,cr,cg,cb)
      return
      end
      subroutine gqeci(idwk,n,ierr,ncoli,icol)
c inquire list element of color indices
c input arguments: idwk, n
c idwk = workstation identifier
c n = requested list element
c ierr = error indicator
c ncoli = number of indices currently defined
c icol = color index of requested element
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c iwk = current workstation identifier
c cga format
      if (iwt.eq.1) then
         ncoli = 4
c ega format
      elseif (iwt.eq.2) then
         ncoli = 2
c printer plots
      elseif (iwt.eq.8) then
         ncoli = 64
c monochrome tiff format
      elseif ((iwt.eq.9).or.(iwt.eq.11)) then
         ncoli = 2
c 256 color formats
      else
         ncoli = 256
      endif
      ierr = 0
      if (idwk.ne.iwk) ierr = -12
      if ((n.lt.0).or.(n.gt.ncoli)) ierr = -12
      if ((n.gt.0).and.(n.le.ncoli)) icol = n - 1
      return
      end
      subroutine gqdsp(iwtype,ierr,idcun,dcx,dcy,lx,ly)
c inquire maximum display surface size
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator (0=inquiry successful)
c idcun = device coordinate units, (0=meters,1=other)
c dcx/dcy = width/height in device coordinate units
c lx/ly = width/height in device (raster) units
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwt = current workstation type
      dimension lxs(12), lys(12)
c lxs, lys = screen sizes for nersc format raster
      data lxs /320,640,320,320,320,400,640,78,320,320,640,640/
      data lys /200,350,200,200,320,400,400,22,200,200,400,400/
      idcun = 1
c if workstation requested is open, query actual size
      if ((kosv.eq.3).and.(iwtype.eq.iwt)) then
c get screen size
         call gscreen(lx,ly)
c return maximum possible size, if not actually opened
      else
         lx = lxs(iwt)
         ly = lys(iwt)
      endif
      dcx = float(lx) - 1.
      dcy = float(ly) - 1.
      ierr = 0
      return
      end
      subroutine gswkwn(idwk,xmin,xmax,ymin,ymax)
c set workstation window
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = window x coordinates in ndc
c ymin/ymax = window y coordinates in ndc
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c trans = normalization transformations
c clip to maximum window size
      wwnvp(1) = amax1(xmin,0.)
      wwnvp(2) = amin1(xmax,1.)
      wwnvp(3) = amax1(ymin,0.)
      wwnvp(4) = amin1(ymax,1.)
c redefine viewport of transformation 0
      trans(1,1) = wwnvp(1)
      trans(2,1) = wwnvp(2)
      trans(3,1) = wwnvp(3)
      trans(4,1) = wwnvp(4)
      return
      end
      subroutine gswkvp(idwk,xmin,xmax,ymin,ymax)
c set workstation viewport
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = viewport x coordinates in device coordinates
c ymin/ymax = viewport y coordinates in device coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c get screen size
      call gscreen(lx,ly)
      dcx = float(lx) - 1.
      dcy = float(ly) - 1.
c clip to maximum screen size
      wwnvp(5) = amax1(xmin,0.)
      wwnvp(7) = amax1(ymin,0.)
      wwnvp(6) = amin1(xmax,dcx)
      wwnvp(8) = amin1(ymax,dcy)
      return
      end
      subroutine gqsts(idwk,idstr,mldr,ierr,mode,iesw,lstr,str,ipet,eare
     1a,lenb,ipos,ldr,datar)
c inquire string device state
c input arguments: idwk, idstr, mldr
c idwk = workstation identifier
c idstr = string device number
c mldr = maximum dimension of data record array
c ierr = error indicator (0=inquiry successful)
c mode = operating mode (0=request,1=sample,2=event)
c iesw = echo switch (0=no echo,1=echo)
c lstr = number of characters in initial string
c str = initial string
c ipet = prompt/echo type (1=normal)
c earea(1)/earea(2) = echo area x coordinates in device coordinates
c earea(3)/earea(4) = echo area y coordinates in device coordinates
c lenb = input buffer size
c ipos = inital edit position
c ldr = length of data record array
c datar = data record array
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c cea = character echo area
      character*(*) str
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 1
      lstr = 1
      str(1:1) = ' '
      ipet = 1
      do 10 j = 1, 4
      earea(j) = cea(j)
   10 continue
      lenb = 80
      ipos = 1
      ldr = 1
c fortran input supported
      if (iwt.eq.8) then
         ierr = 0
c no input supported
      else
         ierr = 140
      endif
      return
      end
      subroutine gqlcs(idwk,idloc,itype,mldr,ierr,mode,iesw,nrt,pxi,pyi,
     1ipet,earea,ldr,datar)
c inquire locator device state
c input arguments: idwk, idloc, itype, mldr
c idwk = workstation identifier
c idloc = locator device number
c itype = type of returned value (0=set,1=realized)
c mldr = maximum dimension of data record array
c ierr = error indicator (0=inquiry successful)
c mode = operating mode (0=request,1=sample,2=event)
c iesw = echo switch (0=no echo,1=echo)
c nrt = normalization transformation number for initial position
c pxi/pyi = initial position, in world coordinates
c ipet = prompt/echo type (1=normal)
c earea(1)/earea(2) = echo area x coordinates in device coordinates
c earea(3)/earea(4) = echo area y coordinates in device coordinates
c ldr = length of data record array
c datar = data record array
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 0
      nrt = 0
      pxi = .5
      pyi = .5
      ipet = 1
      earea(1) = 0.
      earea(2) = 1.
      earea(3) = 0.
      earea(4) = 1.
      ldr = 1
      ierr = 140
      return
      end
      subroutine ginst(idwk,idstr,lstr,str,ipet,xmin,xmax,ymin,ymax,lenb
     1,ipos,ldr,datar)
c initialize string device
c input arguments: all
c idwk = workstation identifier
c idstr = string device number
c lstr = number of characters in initial string
c str = initial string
c ipet = prompt/echo type (1=normal)
c xmin/xmax = echo area x coordinates in device coordinates
c ymin/ymax = echo area y coordinates in device coordinates
c lenb = input buffer size
c ipos = inital edit position
c ldr = length of data record array
c datar = data record array
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c cea = character echo area
      character*80 datar(ldr)
      character*(*) str
c save characer echo area
      cea(1) = xmin
      cea(2) = xmax
      cea(3) = ymin
      cea(4) = ymax
      return
      end
      subroutine ginlc(idwk,idloc,nrt,px,py,ipet,xmin,xmax,ymin,ymax,ldr
     1,datar)
c initialize locator device
c input arguments: all
c idwk = workstation identifier
c idloc = locator device number
c nrt = initial normalization transformation number
c px/py = initial locator position, in world coordinates
c ipet = prompt/echo type (1=normal)
c xmin/xmax = echo area x coordinates in device coordinates
c ymin/ymax = echo area y coordinates in device coordinates
c ldr = length of data record array
c datar = data record array
      character*80 datar(ldr)
      return
      end
      subroutine gdawk(idwk)
c deactivate workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwt = current workstation type
c terminate nersc raster device
      if (iwt.lt.8) call gmclose
c terminate tiff raster device
      if (iwt.gt.8) call gmclose
      kosv = 2
      return
      end
      subroutine gclwk(idwk)
c close workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 1
      return
      end
      subroutine gclks
c close gks
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 0
      return
      end
      subroutine gclrwk(idwk,icofl)
c clear workstation
c input arguments: all
c idwk = workstation identifier
c icofl = control flag (0=conditionally,1=always)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c printer plots
      if (iwt.eq.8) then
c        call system('clear')
         call zimage(' ')
c raster images
      else
         call zimage(char(0))
      endif
      return
      end
      subroutine gqcntn(ierr,nrt)
c inquire current normalization transformation number
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c nrt = transformation number (0 <= nrt <= 25)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
      ierr = 0
      nrt = nrtn
      return
      end
      subroutine gqnt(nrt,ierr,window,viewpt)
c inquire normalization transformation
c input arguments: nrt
c nrt = transformation number (0 <= nrt <= 25)
c ierr = error indicator (0=inquiry successful)
c window(1)/window(2) = window x coordinates in world coordinates
c window(3)/window(4) = window y coordinates in world coordinates
c viewpt(1)/viewpt(2) = viewport x coordinates in ndc
c viewpt(3)/viewpt(4) = viewport y coordinates in ndc
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      dimension window(4), viewpt(4)
      ierr = 0
      if ((nrt.lt.0).or.(nrt.ge.maxt)) then
         nrt = 0
         ierr = 1
      endif
      n = nrt + 1
      do 10 j = 1, 4
      window(j) = trans(j+4,n)
      viewpt(j) = trans(j,n)
   10 continue
      return
      end
      subroutine gswn(nrt,xmin,xmax,ymin,ymax)
c set window
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = window x coordinates in world coordinates
c ymin/ymax = window y coordinates in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      n = nrt + 1
c store transformation
      trans(5,n) = xmin
      trans(6,n) = xmax
      trans(7,n) = ymin
      trans(8,n) = ymax
      return
      end
      subroutine gsvp(nrt,xmin,xmax,ymin,ymax)
c set viewport
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = viewport x coordinates in ndc
c ymin/ymax = viewport y coordinates in ndc
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      n = nrt + 1
c store transformation
      trans(1,n) = xmin
      trans(2,n) = xmax
      trans(3,n) = ymin
      trans(4,n) = ymax
      return
      end
      subroutine gselnt(nrt)
c select normalization transformation
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      nrtn = nrt
      n = nrt + 1
c clip to workstation window
      xmin = amax1(wwnvp(1),trans(1,n))
      xmax = amin1(wwnvp(2),trans(2,n))
      ymin = amax1(wwnvp(3),trans(3,n))
      ymax = amin1(wwnvp(4),trans(4,n))
c convert from viewport to screen coordinates
      scx = (wwnvp(6) - wwnvp(5))/(wwnvp(2) - wwnvp(1))
      dvt(1) = (xmin - wwnvp(1))*scx + wwnvp(5)
      dvt(2) = (xmax - wwnvp(1))*scx + wwnvp(5)
      scy = (wwnvp(8) - wwnvp(7))/(wwnvp(4) - wwnvp(3))
      dvt(3) = (ymin - wwnvp(3))*scy + wwnvp(7)
      dvt(4) = (ymax - wwnvp(3))*scy + wwnvp(7)
      return
      end
      subroutine gqln(ierr,ltype)
c inquire linetype
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c ltc = current line type
      ierr = 0
      ltype = ltc
      if ((ltype.lt.1).or.(ltype.gt.4)) then
         ierr = 1
         ltype = 1
      endif
      return
      end
      subroutine gslwsc(alwsc)
c set linewidth scale factor
c input arguments: all
c alwsc = linewidth scale factor, (alwsc > 1.0)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c lws = current linewidth scale factor
      lws = alwsc
      if (lws.lt.1) lws = 1
      if (lws.gt.3) lws = 3
      return
      end
      subroutine gsplci(icol)
c set polyline color index
c input arguments: all
c icol = color index
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c lcl = current line color
c icc = color index conversion table
c cga display
      if (iwt.eq.1) then
c convert color index to cga index
         if ((icol.ge.0).and.(icol.le.3)) then
            lcl = icc(icol+1)
         else
            lcl = icc(2)
         endif
c other displays
      else
         lcl = icol
c convert to printable character
         if (iwt.eq.8) call ppc(lcl)
      endif
      return
      end
      subroutine gsln(ltype)
c set linetype
c input arguments: all
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c ltc = current line type
      if ((ltype.ge.1).and.(ltype.le.4)) then
         ltc = ltype
      else
         ltc = 1
      endif
      return
      end
      subroutine gpl(n,px,py)
c draw polyline
c input arguments: all
c n = number of points to be connected by a polyline
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c kcf = clipping flag
c lcl = current line color
c ltc = current line type
c lws = current linewidth scale factor
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
c clipping is on
      if (kcf.eq.1) then
         xmin = dvt(1)
         xmax = dvt(2)
         ymin = dvt(3)
         ymax = dvt(4)
c clipping is off
      else
         xmin = wwnvp(5)
         xmax = wwnvp(6)
         ymin = wwnvp(7)
         ymax = wwnvp(8)
      endif
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy
c convert to screen coordinates
      ax = px(1)*scx + aminx
      ay = py(1)*scy + aminy
c clip to viewport
      ax0 = amin1(amax1(ax,xmin),xmax)
      ay0 = amin1(amax1(ay,ymin),ymax)
      ix = ax0 + .5
      iy = ay0 + .5
c move point
      call dmove(ix,iy)
      do 10 j = 2, n
c convert to screen coordinates
      ax = px(j)*scx + aminx
      ay = py(j)*scy + aminy
      ierr = 0
c clip x coordinate
      if (ax.lt.xmin) then
         ax = xmin
         if (ax0.eq.xmin) ierr = 1
      elseif (ax.gt.xmax) then
         ax = xmax
         if (ax0.eq.xmax) ierr = 1
      endif
      ax0 = ax
c clip y coordinate
      if (ay.lt.ymin) then
         ay = ymin
         if (ay0.eq.ymin) ierr = 1
      elseif (ay.gt.ymax) then
         ay = ymax
         if (ay0.eq.ymax) ierr = 1
      endif
      ay0 = ay
      ix = ax + .5
      iy = ay + .5
c draw line
      if (ierr.eq.0) then
         call dashln(ix,iy,lcl,ltc,lws)
c move cursor
      else
         call dmove(ix,iy)
      endif
   10 continue
      return
      end
      subroutine gqmk(ierr,mtype)
c inquire marker type
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c mtype = marker type
c 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross, 6 = square,
c 7 = square with cross, 8 = diamond, 9 = diamond with cross,
c 10 = filled circle
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      ierr = 0
      mtype = mtc
      if ((mtype.lt.1).or.(mtype.gt.10)) then
         ierr = 1
         mtype = 3
      endif
      return
      end
      subroutine gsmksc(amksc)
c set marker size scale factor
c input arguments: all
c amksc = linewidth scale factor, (amksc > 1.0)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c ams = current marker size scale factor
      ams = amksc
      return
      end
      subroutine gspmci(imcol)
c set polymarker color index
c input arguments: all
c imcol = polymarker color index
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c lcm = current marker color
c icc = color index conversion table
c cga display
      if (iwt.eq.1) then
c convert color index to cga index
         if ((imcol.ge.0).and.(imcol.le.3)) then
            lcm = icc(imcol+1)
         else
            lcm = icc(2)
         endif
c other displays
      else
         lcm = imcol
      endif
      return
      end
      subroutine gsmk(mtype)
c set marker type
c input arguments: all
c mtype = marker type
c 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross, 6 = square,
c 7 = square with cross, 8 = diamond, 9 = diamond with cross,
c 10 = filled circle
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      if ((mtype.ge.1).and.(mtype.le.10)) then
         mtc = mtype
      else
c set marker symbol to star
         mtc = 3
      endif
      return
      end
      subroutine gpm(n,px,py)
c draw polymarker
c input arguments: all
c n = number of points to be connected by a polymarker
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c nrtn = current normalization transformation number
c kcf = clipping flag
c lcm = current marker color
c mtc = current marker type
c ams = current marker size scale factor
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
      character*1 amks(10)
      save amks
      data amks /'.','+','*','o','x','#','H','V','W','@'/
c clipping is on
      if (kcf.eq.1) then
         xmin = dvt(1)
         xmax = dvt(2)
         ymin = dvt(3)
         ymax = dvt(4)
c clipping is off
      else
         xmin = wwnvp(5)
         xmax = wwnvp(6)
         ymin = wwnvp(7)
         ymax = wwnvp(8)
      endif
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy
c calculate scaling factors
      sy = .4*ams
      sx = sy
c set direction cosines of character up vector
      upx = 0.0
      upy = 1.0
c find center of character
      dx = 5.5*sy
      dy = 6.5*sy
c set integer clipping range
      minx = xmin + .5
      maxx = xmax + .5
      miny = ymin + .5
      maxy = ymax + .5
      do 10 j = 1, n
c convert to screen coordinates
      ax = px(j)*scx + aminx
      ay = py(j)*scy + aminy
c plot if within clipping window
      if ((ax.ge.xmin).and.(ax.le.xmax).and.(ay.ge.ymin).and.(ay.le.ymax
     1)) then
c printer plots
         if (iwt.eq.8) then
c convert to integer 
            ix = ax + .5
            iy = ay + .5
c write marker
            call cwrite(amks(mtc),ix,iy)
c raster images
         else
c special case of dot markers
            if (mtc.eq.1) then
c convert to integer 
               ix = ax + .5
               iy = ay + .5
c draw dot
               call dpnt(ix,iy,lcm)
c other markers
            else
c convert to integer 
               ix = (ax - dx) + .5
               iy = (ay - dy) + .5
c draw character string
               call csdraw(amks(mtc),ix,iy,sx,sy,upx,upy,minx,maxx,miny,
     1maxy,lcm)
            endif
         endif
      endif
   10 continue
      return
      end
      subroutine gqtxal(ierr,itxalh,itxalv)
c inquire text alignment
c input arguments: none
c ierr = error indicator
c itxalh = horizontal component
c 0 = normal, 1 = left, 2 = center, 3 = right
c itxalv = vertical component:
c = 0 normal, 1 = top, 2 = cap, 3 = half, 4 = base, 5 = bottom
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c itx = current horizontal text alignment
c ity = current vertical text alignment
      ierr = 0
      itxalh = itx
      itxalv = ity
      return
      end
      subroutine gqchh(ierr,chh)
c inquire character height
c input arguments: none
c ierr = error indicator
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c cch = current character height
      chh = cch
      ierr = 0
      return
      end
      subroutine gqtxp(ierr,itxp)
c inquire text path
c input arguments: none
c ierr = error indicator
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c itxtp = current text path
      itxp = itxtp
      ierr = 0
      return
      end
      subroutine gqchup(ierr,chux,chuy)
c inquire character up vector
c input arguments: none
c ierr = error indicator
c chux/chuy = up vector, in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c chuv = character up vector
      chux = chuv(1)
      chuy = chuv(2)
      ierr = 0
      return
      end
      subroutine gstxal(itxalh,itxalv)
c set text alignment
c input arguments: all
c itxalh = horizontal component:
c 0 = normal, 1 = left, 2 = center, 3 = right
c itxalv = vertical component:
c = 0 normal, 1 = top, 2 = cap, 3 = half, 4 = base, 5 = bottom
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c itx = current horizontal text alignment
c ity = current vertical text alignment
      if ((itxalh.ge.0).and.(itxalh.le.3)) itx = itxalh
      if ((itxalv.ge.0).and.(itxalh.le.3)) ity = itxalv
      return
      end
      subroutine gstxfp(nfont,iprec)
c set text font
c input arguments: all
c nfont = character font number
c iprec = text precision (0=string,1=character,2=stroke)
      return
      end
      subroutine gstxp(itxp)
c set text path
c input arguments: all
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c itxtp = current text path
      itxtp = itxp
      return
      end
      subroutine gstxci(itcol)
c set text color index
c input arguments: all
c itcol = text color index
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c lct = current text color
c icc = color index conversion table
c cga display
      if (iwt.eq.1) then
c convert color index to cga index
         if ((itcol.ge.0).and.(itcol.le.3)) then
            lct = icc(itcol+1)
         else
            lct = icc(2)
         endif
c other displays
      else
         lct = itcol
      endif
      return
      end
      subroutine gschh(chh)
c set character height
c input arguments: all
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c cch = current character height
      cch = chh
      return
      end
      subroutine gschup(chux,chuy)
c set character up vector
c input arguments: all
c chux/chuy = up vector, in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c chuv = character up vector
      chuv(1) = chux
      chuv(2) = chuy
      return
      end
      subroutine gschxp(chxp)
c set character expansion factor
c input arguments: all
c chxp = character expansion factor (>0.)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c chx = character width expansion factor
      if (chxp.gt.0.) then
         chx = chxp
      else
         chx = 1.0
      endif
      return
      end
      subroutine gtx(px,py,chars)
c display text
c input arguments: all
c px/py = starting x/y position of text, in world coordinates
c chars = test string to be displayed
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c kcf = clipping flag
c lct = current text color
c cch = current character height
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
c chuv = character up vector
      character*(*) chars
      n = len(chars)
c clipping is on
      if (kcf.eq.1) then
         xmin = dvt(1)
         xmax = dvt(2)
         ymin = dvt(3)
         ymax = dvt(4)
c clipping is off
      else
         xmin = wwnvp(5)
         xmax = wwnvp(6)
         ymin = wwnvp(7)
         ymax = wwnvp(8)
      endif
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy
c calculate character box size, in user coordinates
      cbx = cch*chx*(scy/scx)
      cby = cch
c calculate scaling factor
      sy = cch*scy/15.
      sx = sy*chx
c set direction cosines of character up angle
      achu = 1./sqrt(chuv(1)*chuv(1) + chuv(2)*chuv(2))
      upx = chuv(1)*achu
      upy = chuv(2)*achu
c determine horizontal offset and termination
      cx = px
      dx = float(n)*cbx
      if (itx.eq.2) then
         cx = cx - .5*dx
      elseif (itx.eq.3) then
         cx = cx - dx
      endif
c determine vertical offset
      cy = py
      if (ity.eq.3) then
         cy = cy - .5*cby
      elseif ((ity.eq.1).or.(ity.eq.2)) then
         cy = cy - cby
      endif
c convert to screen coordinates
      ax = cx*scx + aminx
      ay = cy*scy + aminy
c set integer clipping range
      minx = xmin + .5
      maxx = xmax + .5
      miny = ymin + .5
      maxy = ymax + .5
c plot if within clipping window
      if ((ax.le.xmax).and.(ay.ge.ymin).and.(ay.le.ymax)) then
c convert to integer
         ix = ax + .5
         iy = ay + .5
c write character
         if (iwt.eq.8) then
            call cwrite(chars,ix,iy)
c draw character string
         else
            call csdraw(chars,ix,iy,sx,sy,upx,upy,minx,maxx,miny,maxy,lc
     1t)
         endif
      endif
      return
      end
      subroutine gsfaci(ifcol)
c set fill area color index
c input arguments: all
c ifcol = color index
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c lcf = current fill color
c icc = color index conversion table
c cga display
      if (iwt.eq.1) then
c convert color index to cga index
         if ((ifcol.ge.0).and.(ifcol.le.3)) then
            lcf = icc(ifcol+1)
         else
            lcf = icc(2)
         endif
c other displays
      else
         lcf = ifcol
c convert to printable character
         if (iwt.eq.8) call ppc(lcf)
      endif
      return
      end
      subroutine gsfais(ints)
c set fill area interior style
c input arguments: all
c ints = desired interior style:
c 0 = hollow (default), 1 = solid, 2 = pattern, 3 = hatch
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c infs = interior fill style
      infs = ints
      return
      end
      subroutine gfa(n,px,py)
c fill area
c input arguments: all
c n = number of points in fill area
c px,py = arrays of points, in world coordinates
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c lcf = current fill color
c kcf = clipping flag
c infs = interior fill style
c lws = current linewidth scale factor
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
c clipping is on
      if (kcf.eq.1) then
         xmin = dvt(1)
         xmax = dvt(2)
         ymin = dvt(3)
         ymax = dvt(4)
c clipping is off
      else
         xmin = wwnvp(5)
         xmax = wwnvp(6)
         ymin = wwnvp(7)
         ymax = wwnvp(8)
      endif
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy
c convert to screen coordinates
      ax = px(1)*scx + aminx
      ay = py(1)*scy + aminy
c clip to viewport
      ax0 = amin1(amax1(ax,xmin),xmax)
      ay0 = amin1(amax1(ay,ymin),ymax)
      ix = ax0 + .5
      iy = ay0 + .5
c move point
      call dmove(ix,iy)
      do 10 j = 2, n
c convert to screen coordinates
      ax = px(j)*scx + aminx
      ay = py(j)*scy + aminy
      ierr = 0
c clip x coordinate
      if (ax.lt.xmin) then
         ax = xmin
         if (ax0.eq.xmin) ierr = 1
      elseif (ax.gt.xmax) then
         ax = xmax
         if (ax0.eq.xmax) ierr = 1
      endif
      ax0 = ax
c clip y coordinate
      if (ay.lt.ymin) then
         ay = ymin
         if (ay0.eq.ymin) ierr = 1
      elseif (ay.gt.ymax) then
         ay = ymax
         if (ay0.eq.ymax) ierr = 1
      endif
      ay0 = ay
      ix = ax + .5
      iy = ay + .5
c draw line
      if (ierr.eq.0) then
         call dline(ix,iy,lcf,lws)
c move cursor
      else
         call dmove(ix,iy)
      endif
   10 continue
c convert to screen coordinates
      ax = px(1)*scx + aminx
      ay = py(1)*scy + aminy
      ierr = 0
c clip x coordinate
      if (ax.lt.xmin) then
         ax = xmin
         if (ax0.eq.xmin) ierr = 1
      elseif (ax.gt.xmax) then
         ax = xmax
         if (ax0.eq.xmax) ierr = 1
      endif
c clip y coordinate
      if (ay.lt.ymin) then
         ay = ymin
         if (ay0.eq.ymin) ierr = 1
      elseif (ay.gt.ymax) then
         ay = ymax
         if (ay0.eq.ymax) ierr = 1
      endif
      ix = ax + .5
      iy = ay + .5
c draw line
      if (ierr.eq.0) then
         call dline(ix,iy,lcf,lws)
c move cursor
      else
         call dmove(ix,iy)
      endif
      return
      end
      subroutine gca(px,py,qx,qy,icxd,icyd,ncs,nrs,idx,idy,icola)
c cell array
c input arguments: all
c px,py = lower-left cell corner, in world coordinates
c qx,qy = upper-right cell corner, in world coordinates
c icxd,icyd = color index array dimensions
c ncs,nrs = starting column and row in the color index array
c idx,idy = number of columns and rows in the cell array
c icola = color index array
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c nrtn = current normalization transformation number
c icc = color index conversion table
c wwnvp = workstation window and viewport
c dvt = device transformation
c trans = normalization transformations
      dimension icola(icxd,icyd)
c lxm, lym = maximum address of pixels in x, y
      lxm = wwnvp(6)
      lym = wwnvp(8)
c find location of upper left and lower right hand corner of image
      xu = px
      xl = qx
      yu = amax1(py,qy)
      yl = amin1(py,qy)
c calculate transformation factors
      m = nrtn + 1
      scx = (dvt(2) - dvt(1))/(trans(6,m) - trans(5,m))
      aminx = dvt(1) - trans(5,m)*scx
      scy = (dvt(4) - dvt(3))/(trans(8,m) - trans(7,m))
      aminy = dvt(3) - trans(7,m)*scy
c convert to screen coordinates
      ax = xu*scx + aminx
      ay = yu*scy + aminy
      bx = xl*scx + aminx
      by = yl*scy + aminy
c clip to workstation viewport
      ix0 = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
      iy0 = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
      ix1 = amin1(amax1(bx,wwnvp(5)),wwnvp(6)) + .5
      iy1 = amin1(amax1(by,wwnvp(7)),wwnvp(8)) + .5
c calculate initial offsets
      ipx = ix0 - 1
      ipy = iy0 + 1
      apx = float(ix0) + .5
      apy = 1.5
c calculate maximum index
      lxs = idx
c     if ((ipx+idx).gt.lxm) lxs = lxm - ipx
      lys = idy
c     if ((ipy-idy).lt.0) lys = ipy
c calculate scalings for pixels (dx = dy = 1., for no rescaling)
      dx = float(ix1 - ix0)/float(lxs - 1)
      dy = float(lys - 1)/float(iy0 - iy1)
c invf = (0,1) = (no,yes) image should be inverted vertically
      if (py.ge.qy) then
         invf = 0
         koff = nrs - 1
      else
         invf = 1
         koff = lys - nrs + 2
      endif
      km = iy0 - iy1 + 1
      joff = ncs - 1
c outer loop over rows
      do 20 kk = 1, km
      k = dy*float(kk - 1) + apy
c normal image
      if (invf.eq.0) then
         k1 = k + koff
c inverted image
      else
         k1 = koff - k
      endif
      iy = ipy - kk
c reset previous color code to background
      idr = 0
c cga display
      if (iwt.eq.1) then
         idr = icc(idr+1)
c convert to printable character
      elseif (iwt.eq.8) then
         call ppc(idr)
      endif
c move cursor
      call dmove(ix0,iy)
c next loop over bytes in row of color plane
      do 10 j = 1, lxs
      ix = dx*float(j - 1) + apx
      itc = icola(j+joff,k1)
c cga display
      if (iwt.eq.1) then
         itc = icc(itc+1)
c convert to printable character
      elseif (iwt.eq.8) then
         call ppc(itc)
      endif
c no change in color
      if (itc.eq.idr) then
c draw line in current color if at end of picture
         if (ix.eq.ix1) call dline(ix,iy,idr,1)
c color change
      else
c draw to previous point
         if (ix.gt.ix0) call dline(ix-1,iy,idr,1)
c reset previous color code to current color
         idr = itc
c move to current point
         call dmove(ix,iy)
c draw line in current color if at end of picture
         if (ix.eq.ix1) call dline(ix,iy,idr,1)
      endif
   10 continue
   20 continue
      return
      end
      subroutine gqclip(ierr,indcl,clrect)
c inquire clipping indicator
c input arguments: none
c ierr = error indicator
c indcl = clipping indicator (0=no clip, 1=clip)
c clrect = clipping rectangle, in ndc
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c kcf = clipping flag
c trans = normalization transformations
      dimension clrect(4)
      indcl = kcf
c find clipping rectangle
      n = nrtn + 1
      do 10 j = 1, 4
      clrect(j) = trans(j,n)
   10 continue
      ierr = 0
      return
      end
      subroutine gsclip(iclsw)
c set clipping indicator
c input arguments: all
c iclsw = clipping switch (0=no clip,1=clip)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kcf = clipping flag
      kcf = iclsw
      return
      end
      subroutine guwk(idwk,iregfl)
c update workstation
c input arguments: all
c idwk = workstation identifier
c iregfl = regeneration flag (0=postponed,1=perform)
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
c printer plots
      if (iwt.eq.8) then
         call wimage
c raster images
      else
         call sgraph(iwt)
      endif
      return
      end
      subroutine grqst(idwk,idstr,istat,lostr,str)
c request string
c input arguments: idwk, idstr
c idwk = workstation identifier
c idstr = string device number
c istat = return status (0=none,1=ok)
c lostr = number of characters in string
c str = returned string
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
      character*(*) str
c input = scratch input array
      character*80 input
   91 format (a80)
c fortran input supported
      if (iwt.eq.8) then
         n = len(str)
         lostr = n
         i = 0
c pad input with blanks
         do 10 j = 1, 80
         input(j:j) = ' '
   10    continue
         read (5,91,end=20) input
         str = input
c determine how many characters in input
   20    i = i + 1
         if (i.gt.n) go to 30
         its = ichar(str(i:i))
c try again if printable character found
         if ((its.ge.32).and.(its.ne.ichar(' '))) go to 20
c control character found
         if (its.lt.32) then
            lostr = i
         elseif (i.lt.n) then
c two blank characters in a row will be considered end of record
            if (str(i+1:i+1).eq.' ') then
               lostr = i - 1
c if only one blank is found, try again
            else
               go to 20
            endif
         endif
   30    istat = 1
c fortran input not supported
      else
         istat = 0
         lostr = 0
      endif
      return
      end
      subroutine grqlc(idwk,idloc,istat,nrt,px,py)
c request locator
c input arguments: idwk, idloc
c idwk = workstation identifier
c idloc = locator device number
c istat = return status (0=none,1=ok)
c nrt = normalization transformation number
c px/py = position, in world coordinates
      istat = 0
      nrt = 0
      return
      end
      subroutine gsstm(idwk,idstr,mode,iesw)
c set string mode
c input arguments: all
c idwk = workstation identifier
c idstr = string device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      return
      end
      subroutine gslcm(idwk,idloc,mode,iesw)
c set locator mode
c input arguments: all
c idwk = workstation identifier
c idloc = locator device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      return
      end
      subroutine gwait(tout,idwk,icl,idnr)
c await event
c input arguments: tout
c tout = timeout period, seconds
c idwk = workstation identifier
c icl = input class:
c (0=no class,1=locator,2=stroke,3=valuator,4=choice,5=pick,6=string)
c idnr = device number
      idwk = 0
      icl = 0
      idnr = 0
      return
      end
      subroutine ggtst(lostr,str)
c get string
c input arguments: none
c lostr = number of characters in string
c str = input string
      character*(*) str
      lostr = 0
      return
      end
      subroutine ggtlc(nrt,px,py)
c get locator
c input arguments: none
c nrt = normalization transformation number
c px/py = position, in world coordinates
      nrt = 0
      return
      end
      subroutine gesc(idfct,ldi,datai,mldr,ldr,datar)
c escape function
c input arguments: idfct, ldi, datai, mldr
c idfct = escape function identifier
c ldi = length of input data record
c datai = input data record
c mldr = maximum length of output data record
c ldr = length of output data record
c datar = output data record
      character*80 datai(ldi), datar(mldr)
      ldr = 0
c escape functions not supported
      ierr = 180
      return
      end
      subroutine rsdflts
c this subroutine creates default tables for raster driver for gks
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c iwt = current workstation type
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c kcf = clipping flag
c lcl = current line color
c lcm = current marker color
c lct = current text color
c lcf = current fill color
c ltc = current line type
c mtc = current marker type
c infs = interior fill style
c itxtp = current text path
c lws = current linewidth scale factor
c ams = current marker size scale factor
c chx = character width expansion factor
c cch = current character height
c icc = color index conversion table
c wwnvp = workstation window and viewport
c dvt = device transformation
c cea = character echo area
c trans = normalization transformations
c chuv = character up vector
c create color index conversion table
c maps gks indices to cga indices: cga_index = icc(gks_index+1)
      icc(1) = 0
      icc(2) = 3
      icc(3) = 1
      icc(4) = 2
c set default workstation window to square
      wwnvp(1) = 0.
      wwnvp(2) = 1.
      wwnvp(3) = 0.
      wwnvp(4) = 1.
c set default workstation viewport to square
c get screen size
      call gscreen(lx,ly)
      if (lx.gt.ly) then
         wwnvp(5) = .5*float(lx - ly)
         wwnvp(8) = float(ly - 1)
         wwnvp(6) = wwnvp(5) + wwnvp(8)
         wwnvp(7) = 0.
      else
         wwnvp(6) = float(lx - 1)
         wwnvp(7) = .5*float(ly - lx)
         wwnvp(5) = 0.
         wwnvp(8) = wwnvp(6) + wwnvp(7)
      endif
c make default character echo area equal to workstation viewport
      do 10 j = 1, 4
      cea(j) = wwnvp(j+4)
   10 continue
c default transformation
      nrtn = 0
      trans(1,1) = 0.
      trans(2,1) = 1.
      trans(3,1) = 0.
      trans(4,1) = 1.
      trans(5,1) = 0.
      trans(6,1) = 1.
      trans(7,1) = 0.
      trans(8,1) = 1.
      do 30 k = 2, maxt
      do 20 j = 1, 8
      trans(j,k) = trans(j,1)
   20 continue
   30 continue
c convert from viewport to screen coordinates
      dvt(1) = wwnvp(5)
      dvt(2) = wwnvp(6)
      dvt(3) = wwnvp(7)
      dvt(4) = wwnvp(8)
c set default clipping state to on
      kcf = 1
c default colors
      lcl = 1
      lcm = 1
      lct = 1
      lcf = 1
c set default line type to solid
      ltc = 1
c set default marker symbol to star
      mtc = 3
c default text alignment
      itx = 0
      ity = 0
c default linewidth scale factor
      lws = 1
c default marker size scale factor
      ams = 1.0
c set current text path to right
      itxtp = 0
c default character width expansion factor
      chx = 1.0
c default character height
      cch = 0.01
c set default character up vector
      chuv(1) = 0.
      chuv(2) = 1.
c set default fill area interior style to hollow
      infs = 0
      return
      end
c internal library for nersc raster images
      subroutine gminit(igtype)
c this subroutine initializes nersc raster device
c input arguments: all
c igtype = (1,2,3,4,5,6,7) = (cga,ega,vga,mac1,mac2,mac3,mac4) format
c igtype = 8 = character plots, igtype = (9,11) = (tif1,tif3) format
c igtype = (10,12) = (tif2,tif4) format
c default is vga format
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c g = current image array
      character*1 chr
      dimension lxs(12), lys(12)
c lxs, lys = screen sizes for nersc format raster
      data lxs /320,640,320,320,320,400,640,78,320,320,640,640/
      data lys /200,350,200,200,320,400,400,22,200,200,400,400/
      save /grsdev/
      id = igtype
c default is vga if format type is invalid
      if ((id.lt.1).or.(id.gt.12)) id = 3
c lx, ly = the size of the image
      lx = lxs(id)
      ly = lys(id)
c check if requested image overflows image array
      if ((lx*ly).gt.(lxm*lym)) then
         write (6,*) 'image too large: lx,ly,lxm,lym=',lx,ly,lxm,lym
         stop 1
      endif
c nbit = number of bits per pixel, nbit < 9
      nbit = 8
c cga format has 4 colors
      if (id.eq.1) nbit = 2
c ega format is monochrome 
      if (id.eq.2) nbit = 1
c monochrome tiff
      if ((id.eq.9).or.(id.eq.11)) nbit = 1
c intrl = (0,1) = (no,yes) interlace image
      intrl = 0
c cga format is interlaced
      if (id.eq.1) intrl = 1
c if ixor = 1, current image is xored with previous image
      ixor = 0
c ipmx =  maximum color value in the palette
      ipmx = 0
c vga format has maximum value of only 64
      if (id.eq.3) ipmx = 64
c mac formats have maximum values of 256
      if ((id.ge.4).and.(id.le.7)) ipmx = 256
c color tiff format has maximum value of 256
      if ((id.eq.10).or.(id.eq.12)) ipmx = 256
c determine header code
      chr = char(240+id)
c write out header for mfe format movies
      if (igtype.lt.8) call buffwr(chr,1)
c write out tiff header
      if (igtype.ge.9) call tfhead
c read machine architecture description
      call starch
      return
      end
      subroutine swpal(ic,cr,cg,cb)
c if ic >= 0, convert palette entry to integer ip in the range,
c 0<=ichar(ip)<=ipmx-1, and save palette entry in character array pal
c if ic < 0, write palette to nersc output device
c if ic >= npald, write palette to tiff output device
c input arguments: all
c ic = color index
c cr/cg/cb = red/green/blue component (0 < cr,cg,cb < 1)
c npald = maximum number of palette entries
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      parameter(npald=256,lpald=3*npald)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c ipmx =  maximum color value in the palette
c pal(i+3*j) = ith component of jth palette entry
      character*1 pal(lpald)
c npal = number of palette entries (0 = use default)
      data npal /-1/
      save pal,npal
c on first entry, set default palette
      if (npal.lt.0) then
c clear all entries
         do 10 j = 1, lpald
         pal(j) = char(0)
   10    continue
c set index 1 to foreground
         pal(4) = char(ipmx-1)
         pal(5) = char(ipmx-1)
         pal(6) = char(ipmx-1)
         npal = 0
      endif
c write palette to nersc output device, if vga or mac raster format
      if ((ic.lt.0).and.(ipmx.gt.0)) then
         if (npal.gt.0) npal = npald
c write palette
         call wrpal(pal,npal,lpald,ipmx)
      endif
c write palette to tiff output device
      if ((ic.ge.npald).and.(ipmx.gt.0)) then
         if (npal.gt.0) npal = npald
c write palette
         call wrtpal(pal,npal,lpald,ipmx)
      endif
c return if index is out of range
      if ((ic.lt.0).or.(ic.ge.npald)) return
      jc = 3*ic
      apmx = float(ipmx-1)
c store palette entry
      ip = apmx*cr
      pal(jc+1) = char(ip)
      ip = apmx*cg
      pal(jc+2) = char(ip)
      ip = apmx*cb
      pal(jc+3) = char(ip)
c count entries
      npal = npal + 1
      return
      end
      subroutine wrpal(pal,npal,lpald,ipmx)
c this subroutine writes palette for mfe vga or mac format files
c pal is a character array, with rgb values in successive bytes
c pal(i+3*j) = ith component of jth palette entry
c npal = (0,n) = (default,n) palette entries
c lpald = size of palette array
c ipmx =  maximum color value in the palette
      character*1 pal(lpald)
c color = scratch array for 8 entries
      character*1 color(24)
c icolor = integer array to define 8 primary colors
      dimension icolor(24)
      save icolor
c data for 8 primary colors
      data icolor /0,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,0,1,1,1/
c write palette code
      color(1) = char(240)
      call buffwr(color,1)
      ipm = ipmx - 1
c write default palette
      if (npal.le.0) then
c first 8 entries define primary default colors
         do 10 i = 1, 24
         color(i) = char(ipm*icolor(i))
   10    continue
         call buffwr(color,24)
c pad remaining 248 entries with foreground (white)
         do 20 i = 1, 24
         color(i) = char(ipm)
   20    continue
         do 30 i = 1, 31
         call buffwr(color,24)
   30    continue
c write specified palette
      else
         len = 3*npal
         if (len.gt.768) len = 768
         call buffwr(pal,len)
c pad remaining entries in palette with foreground (white)
         if (len.lt.768) then
            nl = (768 - len)
c n = number of full groups of 8 entries
            n = nl/24
c nl = remaining number of characters to write
            nl = nl - 24*n
c pad character array
            do 40 i = 1, 24
            color(i) = char(ipm)
   40       continue
c write full groups of 8 entries
            if (n.gt.0) then
               do 50 i = 1, n
               call buffwr(color,24)
   50          continue
            endif
c write remainder
            if (nl.gt.0) call buffwr(color,nl)
         endif
      endif
      return
      end
      subroutine gscreen(lxd,lyd)
c this subroutine returns information about screen size
c input arguments: none
c lxd, lyd = the size of the image
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
      lxd = lx
      lyd = ly
      return
      end
      subroutine gmclose
c this subroutine terminates nersc raster device
c input arguments: none
      character*1 chr(4)
c write four null characters
      do 10 i = 1, 4
      chr(i) = char(0)
   10 continue
      call buffwr(chr,4)
c flush buffers
      call buffwr(chr,0)
      return
      end
      subroutine zimage(blank)
c this subroutine clears current image
c input arguments: all
c blank = character used to clear image
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
c g = current image array
      character*1 blank
c fill current image array with null characters
      lxy = lx*ly
      do 10 j = 1, lxy
      g(j) = blank
   10 continue
      return
      end
      subroutine dashln(i,j,ict,lnt,lwtype)
c this subroutine draws a dashed line from (icx,icy) to (i,j)
c input arguments: all
c i, j = location of point
c ict = color index
c lnt = line type
c 1 = solid, 2 = short-dash, 3 = dot, 4 = dash-dot, 5 = long-dash
c lwtype = linewidth
      dimension lt(4,4)
      save ns,nl,ltype
c ltype = current line type
      data ltype /0/
c lt = software dash type patterns
      data lt /9,6,9,6,5,5,5,5,14,6,4,6,23,7,23,7/
      l = lnt - 1
c use current style if line style is known
      if ((l.lt.0).or.(l.gt.7)) l = ltype
c reset line style, if necessary
      if (l.ne.ltype) then
         ltype = l
c reset software line style
         if ((ltype.ge.1).and.(ltype.le.4)) then
         ns = 0
         nl = lt(ns+1,ltype)
         endif
      endif
c use solid style if requested or unknown
      if ((ltype.eq.0).or.(ltype.gt.4)) then
c draw line
         call dline(i,j,ict,lwtype)
         return
      endif
c icx, icy = current location of cursor
      call gtloc(icx,icy)
c software dashed line
      cost = float(i - icx)
      sint = float(j - icy)
      alen = sqrt(cost*cost + sint*sint)
      len = alen + .5
c current pattern segment is longer than line length
      if (nl.ge.len) go to 20
c find starting location and direction cosines
      ax0 = float(icx) + .5
      ay0 = float(icy) + .5
      cost = cost/alen
      sint = sint/alen
c iterate pattern segments
   10 anl = float(nl)
c find end coordinate of next segment
      ix = ax0 + anl*cost
      iy = ay0 + anl*sint
c dark or bright vector flag
      it = ns - (ns/2)*2
c draw line
      if (it.eq.0) then
         call dline(ix,iy,ict,lwtype)
c move cursor
      else
         call dmove(ix,iy)
      endif
c increment pattern segment index
      ns = ns + 1
      if (ns.eq.4) ns = 0
c add length of next pattern segment
      nl = nl + lt(ns+1,ltype)
c do next segment
      if (nl.lt.len) go to 10
c finish up last segment, which may be incomplete
   20 it = ns - (ns/2)*2
c draw line
      if (it.eq.0) then
         call dline(i,j,ict,lwtype)
c move cursor
      else
         call dmove(i,j)
      endif
c adjust length of next pattern segment
      nl = nl - len
c if segment complete, reset to next pattern segment
      if (nl.eq.0) then
         ns = ns + 1
         if (ns.eq.4) ns = 0
         nl = lt(ns+1,ltype)
      endif
      return
      end
      subroutine dmove(i,j)
c this subroutine moves cursor to location (i,j)
c input arguments: all
c i, j = cursor location
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c icx, icy = current location of cursor
      icx = i
      icy = j
      return
      end
      subroutine gtloc(i,j)
c this subroutine returns the current cursor location
c input arguments: none
c i, j = cursor location
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c icx, icy = current location of cursor
      i = icx
      j = icy
      return
      end
      subroutine dline(i,j,ict,lwtype)
c this subroutine draws line from (icx,icy) to (i,j)
c input arguments: all
c i, j = location of end of line
c ict = color index
c lwtype = linewidth
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
c icx, icy = current location of cursor
c g = current image array
c do not draw zero width lines, but save cursor
      if (lwtype.lt.1) then
c save last point
         icx = i
         icy = j
         return
      endif
c find extent of line segment
      ii = i - icx
      jj = j - icy
c default direction and slope of line
      kk = 1
      aji = 1.
c half-width of line
      alwtype = .5*float(lwtype)
c line closer to horizontal
      if (iabs(jj).le.iabs(ii)) then
c line decreases in x
         if (ii.lt.0) kk = -1
c find slope of line
         if (ii.ne.0) aji = float(jj)/float(ii)
c outer loop over number of pixels in width of line
         do 20 n = 1, lwtype
         aj = float(icy + n) - alwtype
         do 10 l = icx, i, kk
c clip in x direction
         if ((l.ge.0).and.(l.lt.lx)) then
            m = aji*float(l - icx) + aj
c clip in y direction
            if ((m.ge.0).or.(m.lt.ly)) then
               g(lx*m+l+1) = char(ict)
            endif
         endif
   10    continue
   20    continue
c line closer to vertical
      else
c line decreases in y
         if (jj.lt.0) kk = -1
c find slope of line
         aji = float(ii)/float(jj)
c outer loop over number of pixels in width of line
         do 40 n = 1, lwtype
         ai = float(icx + n) - alwtype
         do 30 m = icy, j, kk
c clip in y direction
         if ((m.ge.0).or.(m.lt.ly)) then
            l = aji*float(m - icy) + ai
c clip in x direction
            if ((l.ge.0).and.(l.lt.lx)) then
               g(lx*m+l+1) = char(ict)
            endif
         endif
   30    continue
   40    continue
      endif
c save last point
      icx = i
      icy = j
      return
      end
      subroutine dpnt(i,j,ict)
c this subroutine draws point at location (i,j)
c input arguments: all
c i, j = location of point
c ict = color index
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
c icx, icy = current location of cursor
c g = current image array
c clip in y direction
      if ((j.ge.0).or.(j.lt.ly)) then
c clip in x direction
         if ((i.ge.0).and.(i.lt.lx)) then
            g(lx*j+i+1) = char(ict)
         endif
      endif
c save last point
      icx = i
      icy = j
      return
      end
      subroutine csdraw(chr,ix,iy,sx,sy,upx,upy,minx,maxx,miny,maxy,ict)
c this subroutine draws a string at location ix,iy using a scalable,
c rotatable tektronix "stick" font
c input: all
c chr = input string to be drawn
c ix, iy = x, y coordinate of lower left hand corner of character
c sx, sy = scaling factor in x, y relative to a 11 x 13 raster
c upx, upy = x, y direction cosines of character up angle
c minx = the minimum horizontal screen coordinate
c maxx = the maximum horizontal screen coordinate
c miny = the minimum vertical screen coordinate
c maxy = the maximum vertical screen coordinate
c ict = color index used in drawing characters
c the font is defined in an 11 x 13 raster by draws and moves with 4 bit
c addresses.  the initial location ix,iy is the lower left hand corner
c of the character.  addresses are stored as (x,y) pairs (0<x<10,0<y<12)
c and are preceded by a move flag (=15) if the following coordinate
c represents a move rather than a draw.  the 11 x 13 raster is enlarged
c or reduced by the scaling factors (sx,sy).  the character rotation is
c determined by the up vector direction cosines upx, upy, so that the
c actual address used is given by:
c xa = sy*y*upx + sx*x*upy, ya = sy*y*upy - sx*x*upx.
c if any part of the character address is outside the range
c minx < ix+xa < maxx and miny < iy+ya < maxy, it is clipped.
c the 4 bit addresses are stored in packed integers in the array icfont.
c the location of the addresses in icfont for a particular character are
c stored in the array icfloc, and the number of addresses for that
c character are stored in the array icflen.
c lw = number of 4 bit font coordinate addresses per integer word
      parameter(lw=8)
c iebc = (0,1) = (no,yes) input characters are in ebcdic
      common /march/ iebc, irvb, longi, mtype
      character*(*) chr
c lb = scratch array of 4 bit font coordinate addresses
      dimension lb(8)
      dimension icflen(94), icfloc(94)
c icfont = character font coordinate address array
      dimension icfon1(114), icfon2(114), icfont(228)
      dimension ieta(256)
      equivalence (icfon1(1), icfont(1)), (icfon2(1), icfont(115))
      save nw,icflen,icfloc,icfont,ieta
c nw = -2**(4*lw-1) + 1, used for adding and subtracting high bit
      data nw /-2147483647/
c icflen = array of numbers of coordinate addresses in font array
      data icflen /10,10,20,35,43,25,5,9,9,15,10,7,5,5,5,32,12,21,28,9,1
     19,23,9,35,23,10,12,7,10,7,18,47,13,23,25,16,14,9,32,12,15,16,9,7,8
     2,6,27,12,32,14,25,10,13,11,24,7,16,9,9,5,9,7,2,5,24,19,17,22,21,18
     3,28,13,10,11,9,11,24,13,19,22,22,11,21,18,16,11,24,7,22,9,23,10,23
     4,13/
c icfloc = array of locations of characters in font coordinate array
      data icfloc /1,3,5,8,13,19,23,24,26,28,30,32,33,34,35,36,40,42,45,
     149,51,54,57,59,64,67,69,71,72,74,75,78,84,86,89,93,95,97,99,103,10
     25,107,109,111,112,113,114,118,120,124,126,130,132,134,136,139,140,
     3142,144,146,147,149,150,151,152,155,158,161,164,167,170,174,176,17
     48,180,182,184,187,189,192,195,198,200,203,206,208,210,213,214,217,
     5219,222,224,227/
c first half of icfont array
      data icfon1 /-184213676,1543503872,-225259652,2030043136,-21785408
     14,1895448655,145227776,-266266598,977822304,-2137868103,-135079911
     27,62914560,-266682549,1008470793,406341963,-164148918,974718726,37
     32244480,-100622159,-1010518880,1612842506,1073741824,-174735360,-1
     467562346,-1073741824,-200853612,-1073741824,-266682532,1358565552,
     5-261460132,1342177280,-217766608,-244752384,-184217600,-266686464,
     6-264207692,-959923574,1209402370,273613227,-224017137,545259520,-2
     757767222,-1433902496,1074397184,-257767222,-1433902225,-2036030848
     8,537001984,-133644182,1610612736,-88031128,1783244801,1048576,-896
     903392,-1608382454,709386336,-255146715,0,-121454432,-2107086262,67
     a1219744,1114139274,-1463812096,-94216064,-1563899222,671219744,-17
     b4747820,1392508928,-174747820,1395720192,-134059840,-259354716,671
     c08864,-234331456,-257767222,-1433055407,1342177280,-145258860,-182
     d0109754,1196968266,1519957700,-1028616142,335939600,123512736,-261
     e464064,210545320,-2046363542,1244135424,-90655036,-1028620222,3359
     f39610,805306368,209492648,-1533906944,212660327,1872756736,2126621
     g12,1610612736,-90655036,-1028620222,335939610,974103395,217082479,
     h-1398800384,-217641136,1559480256,-266205688,684682412,217759850,0
     i,217710592,207137952,211856384,-264207692/
c second half of icfont array
      data icfon2 /-959923574,1209402370,272629760,210545320,-2046427136
     1,-264207692,-959923574,1209402370,273638560,210545320,-2046386176,
     2-267319286,709386848,-2136815158,-1342177280,-184168692,-140928614
     34,-255839736,171622400,-255830774,1522532352,-255839741,86335314,1
     4887478700,-1393505792,-255814299,257337772,-1408435712,0,-18869452
     54,0,-256241664,-154482170,0,-265979360,-1610612736,-171483136,-244
     6152182,1779409956,35684513,217084552,-1972754430,1048576,-92700032
     7,1612843274,268435456,-99954777,-2010642942,545300736,-263566744,-
     82105515998,134811648,1089602227,-976830715,1157627904,-89603392,-1
     9604171702,1605018240,268500992,217084552,-1972764672,-184184997,15
     a43503872,-137943807,1048576,217743434,0,-238894556,83886080,149975
     b683,-2056974506,2022221472,149975688,-1972764672,-266313080,-19727
     c54430,2097152,-268382453,747416230,-2078014208,-99954773,-19432709
     d06,612672768,149963895,-1990197248,-250476534,675430498,-198741606
     e4,-205318906,135868168,1744830464,-91615327,-2145385976,-260025078
     f,1518338048,-260034045,86327122,1887478696,-1460631040,-87414783,2
     g034694,612672768,-259358710,0,-137968204,-1804377004,890636032,-17
     h1602092,1342177280,-205208138,-1770559914,890503936,-261979257,121
     i2833792/
c ebcdic/ascii translation with conventions at ucla oac's ibm 3090vf.
c ascii codes for ebcdic 74,79,95,113,139,155 are non-standard
c ebcdic codes 34,53,106,161,192,208,224 are added for ibm compatibility
      data ieta /0,1,2,3,-1,9,-1,127,-1,-1,-1,11,12,13,14,15,16,17,18,19
     1,-1,-1,8,-1,24,25,-1,-1,28,29,30,31,-1,-1,28,-1,-1,10,23,27,-1,-1,
     2-1,-1,-1,5,6,7,-1,-1,22,-1,-1,30,-1,4,-1,-1,-1,-1,20,21,-1,26,32,-
     31,-1,-1,-1,-1,-1,-1,-1,-1,92,46,60,40,43,124,38,-1,-1,-1,-1,-1,-1,
     4-1,-1,-1,33,36,42,41,59,126,45,47,-1,-1,-1,-1,-1,-1,-1,-1,124,44,3
     57,95,62,63,-1,94,-1,-1,-1,-1,-1,-1,-1,96,58,35,64,39,61,34,-1,97,9
     68,99,100,101,102,103,104,105,-1,123,-1,-1,-1,-1,-1,106,107,108,109
     7,110,111,112,113,114,-1,125,-1,-1,-1,-1,-1,126,115,116,117,118,119
     8,120,121,122,-1,-1,-1,91,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
     9,-1,93,-1,-1,123,65,66,67,68,69,70,71,72,73,-1,-1,-1,-1,-1,-1,125,
     a74,75,76,77,78,79,80,81,82,-1,-1,-1,-1,-1,-1,92,-1,83,84,85,86,87,
     b88,89,90,-1,-1,-1,-1,-1,-1,48,49,50,51,52,53,54,55,56,57,-1,-1,-1,
     c-1,-1,-1/
      n = len(chr)
c clip to viewport
      ix0 = min0(max0(ix,minx),maxx)
      iy0 = min0(max0(iy,miny),maxy)
c move cursor
      call dmove(ix0,iy0)
c main loop over characters in string
      do 50 k = 1, n
c ic = ascii code of character to be drawn (33 <= ic <= 126)
c input is in ascii
      if (iebc.eq.0) then
         ic = ichar(chr(k:k))
c input is in ebcdic
      else
         ic = ieta(ichar(chr(k:k))+1)
      endif
c is = offset index into character font coordinate address array
      is = ic - 32
c blanks
      if (is.eq.0) go to 40
c skip if characters are not printable
      if ((is.lt.0).or.(is.gt.94)) go to 50
      id = 0
c set vertical offset for characters with descenders (g,j,p,q,y)
      if ((is.eq.71).or.(is.eq.74).or.(is.eq.80).or.(is.eq.81).or.(is.eq
     1.89)) id = 4
c mw = high-bit operator
      mw = nw - 1
      lwm = lw - 1
      lwp = lw + 1
c lena = number of 4 bit font coordinates addresses for character ic 
      lena = icflen(is)
c lenw = number of integer words for addresses needed by character ic
      lenw = (lena - 1)/lw + 1
c ioff = offset into character font coordinate array for character ic
      ioff = icfloc(is) - 1
c il = number of addresses in packed integer word
      il = lw
c my = (0,1) = processing (x,y) coordinate address
      my = 0
c im = (0,1) = (move,draw) to new coordinate
      im = 0
c clear clipping flag
      ierr = 0
c main loop over packed integers
      do 30 i = 1, lenw
      it = icfont(i+ioff)
c unpack 4 bit addresses from packed integer
      lb(1) = it
c remove high-bit if present 
      if (it.lt.0) lb(1) = lb(1) - mw
      do 10 j = 1, lwm
      lb(j+1) = lb(j)/16
      lb(j) = lb(j) - lb(j+1)*16
   10 continue
c restore high-bit if necessary
      if (it.lt.0) lb(lw) = lb(lw) + 8
c adjust number of addreses in last packed integer
      if (i.eq.lenw) il = lena - (i - 1)*lw
c draw character
      do 20 j = 1, il
      it = lb(lwp-j)
c x coordinate or move flag
      if (my.eq.0) then
c set move flag
         if (it.eq.15) then
            im = 1
c x coordinate
         else
            at = sx*float(it)
            ax = upy*at + float(ix) + .5
            ay = -upx*at + float(iy) + .5
c set y coordinate flag
            my = 1
         endif
c y coordinate
      else
         at = sy*float(it - id)
         jx = ax + upx*at
         jy = ay + upy*at
c clip x coordinate to viewport
         if (jx.lt.minx) then
            jx = minx
c set clipping flag
            if (ix0.eq.minx) ierr = 1
         elseif (jx.gt.maxx) then
            jx = maxx
c set clipping flag
            if (ix0.eq.maxx) ierr = 1
         endif
         ix0 = jx
c clip y coordinate to viewport
         if (jy.lt.miny) then
            jy = miny
c set clipping flag
            if (iy0.eq.miny) ierr = 1
         elseif (jy.gt.maxy) then
            jy = maxy
c set clipping flag
            if (iy0.eq.maxy) ierr = 1
         endif
         iy0 = jy
c move or draw if visible
         if (ierr.eq.0) then
c draw line
            if (im.eq.0) then
               call dline(jx,jy,ict,1)
c move cursor
            else
               call dmove(jx,jy)
            endif
c move if invisible
         else
            call dmove(jx,jy)
         endif
c reset (x,y) coordinate flag
         my = 0
c reset (move,draw) flag
         im = 0
c clear clipping flag
         ierr = 0
      endif
   20 continue
   30 continue
c find location of next character
   40 at = sx*float(14)
      ix = upy*at + float(ix) + .5
      iy = -upx*at + float(iy) + .5
c clip to viewport
      ix0 = min0(max0(ix,minx),maxx)
      iy0 = min0(max0(iy,miny),maxy)
c move cursor
      call dmove(ix0,iy0)
   50 continue
      return
      end
      subroutine dimagx(image,lx,ly,lz,lenb,npix,nf)
c this subroutine displays raster image stored in character array
c an identity transformation has been assumed
c image = uncompressed single image
c lx, ly = the size of the image, in pixels
c lz = width of picture, in bytes
c lenb = size of picture, in bytes
c npix = number of pixels per byte
c nf = current frame number being processed
c optimized version for nersc raster library
c npald = number of palette entries
      parameter(npald=256)
      parameter(maxt=3)
c ifrg = index of foreground color
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c iwt = current workstation type
c icc = color index conversion table
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      character*1 image(lenb)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = scratch integer image array
      dimension img8(8)
c lbl = character variable for numerical label
      character*8 lbl
      save alx
c alx = x coordinate for numerical labels
      data alx /.86/
c nbit = the number of colors, pixel depth
      nbit = 8/npix
c calculate maximum size of image
      lxs = lx
      if (lxs.gt.isx) lxs = isx
      lys = ly
      if (lys.gt.isy) lys = isy
      ix0 = 0
      iy0 = 0
      ix1 = ix0 + lxs - 1
      iy1 = iy0 + lys - 1
      ipx = ix0 - 1
      ipy = iy1 + 1
c normalize label location
      lafx = alx*float(isx)
      lafy = 0
c clear image
c printer plots
      if (iwt.eq.8) then
         call zimage(' ')
c raster images
      else
         call zimage(char(0))
      endif
c eight bit color
      if (nbit.eq.8) then
c outer loop over rows
         do 20 k = 1, lys
         ioff = lz*(k - 1)
         iy = ipy - k
c reset previous color code to background
         idr = 0
c cga display
         if (iwt.eq.1) then
            idr = icc(idr+1)
c convert to printable character
         elseif (iwt.eq.8) then
            call ppc(idr)
         endif
c move cursor
         call dmove(ix0,iy)
c loop over bytes in row of color plane
         do 10 j = 1, lxs
         ix = j + ipx
c do not use lookup table
         if (lupt.eq.0) then
            itc = ichar(image(j+ioff))
c use lookup table
         else
            itc = ipal(ichar(image(j+ioff))+1)
         endif
c cga display
         if (iwt.eq.1) then
            itc = icc(itc+1)
c convert to printable character
         elseif (iwt.eq.8) then
            call ppc(itc)
         endif
c no change in color
         if (itc.eq.idr) then
c draw line in current color if at end of picture
            if (ix.eq.ix1) call dline(ix,iy,idr,1)
c color change
         else
c draw to previous point
            if (ix.gt.ix0) call dline(ix-1,iy,idr,1)
c reset previous color code to current color
            idr = itc
c move to current point
            call dmove(ix,iy)
c draw line in current color if at end of picture
            if (ix.eq.ix1) call dline(ix,iy,idr,1)
         endif
   10    continue
   20    continue
c nbits per pixel
      else
c maximum width
         lzs = lxs/npix
c convert from nbits per pixel to 8 bits per pixel
         ntc = 2**nbit
         npixm = npix - 1
         npixp = npix + 1
c outer loop over rows
         do 60 k = 1, lys
         ioff = lz*(k - 1)
         iy = ipy - k
c reset previous color code to background
         idr = 0
c cga display
         if (iwt.eq.1) then
            idr = icc(idr+1)
c convert to printable character
         elseif (iwt.eq.8) then
            call ppc(idr)
         endif
c move cursor
         call dmove(ix0,iy)
c loop over bytes in row of color plane
         do 50 j = 1, lzs
         joff = (j - 1)*npix + ipx
         itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
         do 30 i = 1, npixm
         it1 = itc/ntc
         img8(npixp-i) = itc-it1*ntc
         itc = it1
   30    continue
         img8(1) = itc
         do 40 i = 1, npix
         ix = i + joff
c do not use lookup table
         if (lupt.eq.0) then
            itc = img8(i)
c use lookup table
         else
            itc = ipal(img8(i)+1)
         endif
c cga display
         if (iwt.eq.1) then
            itc = icc(itc+1)
c convert to printable character
         elseif (iwt.eq.8) then
            call ppc(itc)
         endif
c no change in color
         if (itc.eq.idr) then
c draw line in current color if at end of picture
            if (ix.eq.ix1) call dline(ix,iy,idr,1)
c color change
         else
c draw to previous point
            if (ix.gt.ix0) call dline(ix-1,iy,idr,1)
c reset previous color code to current color
            idr = itc
c move to current point
            call dmove(ix,iy)
c draw line in current color if at end of picture
            if (ix.eq.ix1) call dline(ix,iy,idr,1)
         endif
   40    continue
   50    continue
   60    continue
      endif
c add label
c first find how many digits in nf
      id = 0
      n = nf
   70 id = id + 1
      n = n/10
      if (n.gt.0) go to 70
c create label template
      lbl = '#       '
c create left justified label
      is = ichar('0')
      if (id.gt.7) id = 7
      ls = 10**(id - 1)
      nt = nf
      do 80 i = 1, id
      i1 = i + 1
      n = nt/ls
      lbl(i1:i1) = char(n+is)
      nt = nt - n*ls
      ls = ls/10
   80 continue
c set color to foreground
      lct = ifrg
c cga display
      if (iwt.eq.1) then
         lct = icc(lct+1)
c convert to printable character
      elseif (iwt.eq.8) then
         call ppc(lct)
      endif
c set length of string
      n = id + 1
c write character
      if (iwt.eq.8) then
         call cwrite(lbl(1:n),lafx,lafy)
c draw character string
      else
         call csdraw(lbl(1:n),lafx,lafy,1.,1.,0.,1.,0,isx-1,0,isy-1,lct)
      endif
c write out plot to device
c printer plots
      if (iwt.eq.8) then
         call wimage
c raster images
      else
         call sgraph(iwt)
      endif
      return
      end
      subroutine sgraph(igtype)
c this subroutine writes out plot to device
c image stored in g array is reduced first to multiple pixels
c per character, then compressed with a run length encoding scheme.
c igtype = (1,2,3,4,5,6,7) = (cga,ega,vga,mac1,mac2,mac3,mac4) format
c igtype = 8 = character plots, igtype = (9,11) = (tif1,tif3) format
c igtype = (10,12) = (tif2,tif4) format
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
c lenbm/lengm = maximum size of reduced/compressed image
c lzm/lzgm = maximum size of reduced/compressed image line
      parameter(lzm=lxm,lenbm=lzm*lym,lengm=lenbm+lzm,lzgm=lzm+1)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
c nbit = number of bits per pixel, nbit < 9
c intrl = (0,1) = (no,yes) interlace image
c if ixor = 1, current image is xored with previous image
c g = current image array
c image, jmage = reduced image arrays
c     character*1 image(lenbm), jmage(lenbm)
      character*1 image(lenbm), jmage(1)
c line/limg = reduced/compressed image line
      character*1 line(lzm), limg(lzgm)
c istrt = (0,1) = (first,subsequent) write to device
c irev = (0,1) = (no,yes) flip image vertically
c ifrmt = format type (1,2,3=ucla format, 4=mfe format)
      data istrt,irev,ifrmt /0,1,4/
      save image, jmage, line, limg, istrt, irev
c npix = number of pixels per character
      npix = 8/nbit
c lz = (lx - 1)/npix + 1, number of characters per image line
      lz = (lx - 1)/npix + 1
c lenb = number of characters in image
      lenb = lz*ly
c leng = maximum length of compressed image
      leng = ly*(lz + (lz - 1)/128 + 1)
c lzg = maximum length of compressed line
      lzg = lz + 1
c reduce image data
      call gimage (image,lz,irev,lenbm)
c write tiff file
      if ((igtype.ge.9).and.(igtype.le.12)) then
         call timage (image,g,lz,ly,nbit,lenb,leng)
c compress image data
      else
c write palette on first write to device
         if (istrt.eq.0) then
            call swpal(-1,0.,0.,0.)
            istrt = 1
         endif
         call pctsav(image,jmage,g,line,limg,ixor,lz,ly,lenb,leng,lzg,if
     1rmt,intrl)
      endif
      return
      end
      subroutine gimage (image,lz,irev,lenb)
c this subroutine reduces image data, which is stored in a character
c array g with 8 bit pixels per character, of which only nbit bits are
c valid. and packs the data into an array image with multiple pixels
c per each character.
c image = output reduced image array
c lz = (lx - 1)/npix + 1, is number of characters per image line,
c where npix = number of pixels per character
c irev = (0,1) = (no,yes) flip image vertically
c lenb = number of characters in image
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
c nbit = number of bits per pixel, nbit < 9
c intrl = (0,1) = (no,yes) interlace image
c g = current image array
      character*1 image(lenb)
c npix = number of pixels per character
      npix = 8/nbit
c ntc = number of valid colors
      ntc = 2**nbit
      ly1 = ly + 1
      lyh = (ly - 1)/2 + 1
c loop over rows
      do 30 k = 1, ly
      k1 = k
      if (irev.eq.1) k1 = ly1 - k
      k2 = k
      k3 = (k - 1)/2
      if (intrl.eq.1) k2 = lyh*((k2 - 1) - 2*k3) + k3 + 1
      j2off = lx*(k1 - 1)
      joff = lz*(k2 - 1)
c loop over columns
      do 20 j = 1, lz
      j1 = (j - 1)*npix
      itc = 0
      jj = j + joff
c extract low order nbits from g array
      do 10 i = 1, npix
      it = 0
      j2 = j1 + i
      if (j2.le.lx) then
c        if (g(j2+j2off).ne.char(0)) it = 1
         is = ichar(g(j2+j2off))
         it = is - (is/ntc)*ntc
      endif
      itc = ntc*itc + it
   10 continue
      image(jj) = char(itc)
   20 continue
   30 continue
      return
      end
      subroutine pctsav(image,jmage,img,line,limg,ixor,lz,ly,lenb,leng,l
     1zg,ifrmt,intrl)
c this subroutine performs compression of successive images.
c input is in array image, and compressed output is in array img.
c if ixor = 1, then the current image is xored with the previous image
c before compression, and the current image is saved in the array jmage.
c line and limg are scratch arrays needed by subroutine comprs
c image = uncompressed input image
c jmage = uncompressed previous image
c img = compressed output image
c line, limg = uncompressed, compressed image line
c ixor = (0,1) = (no,yes) successive frames are differenced
c lz, ly = width, height of picture, in bytes
c lenb = size of picture, in bytes
c leng = maximum size of compressed image
c lzg = maximum size of compressed line
c ifrmt = format type (1,2,3=ucla format, 4=mfe format)
c intrl = (0,1) = (no,yes) image is interlaced
      character*1 image(lenb), jmage(lenb)
      character*1 img(leng), line(lz), limg(lzg)
      character*1 chr(4)
   91 format(36h compressed image overflow, ibmax = ,i20,8h leng = ,i20)
c compress interlaced image
      if (intrl.eq.1) then
         lyh = ly/2
         lyh1 = (ly - 1)/2 + 1
         len = lz*lyh1
c compress odd lines in image
         call comprs(image,img,line,limg,lz,lyh1,lz,len,leng,lzg,ibmax)
         img(ibmax) = char(240)
         jbmax = ibmax
         it1 = len + 1
         it2 = ibmax + 1
         len = lz*lyh
c compress even lines in image
         call comprs(image(it1),img(it2),line,limg,lz,lyh,lz,len,leng,lz
     1g,ibmax)
         ibmax = ibmax + jbmax
         if (ibmax.gt.leng) write (6,91) ibmax, leng
         go to 50
      endif
c compress successive image
      if (ixor.eq.0) then
c copy uncompressed image
         if (ifrmt.eq.3) then
            ibmax = lenb
            do 10 i = 1, ibmax
            img(i) = image(i)
   10       continue
c compress image
         else
            call comprs(image,img,line,limg,lz,ly,lz,lenb,leng,lzg,ibmax
     1)
         endif
         if (ibmax.gt.leng) write (6,91) ibmax, leng
      else
c xor and compress successive images
         do 20 i = 1, lenb
         jmage(i) = char(ieor(ichar(image(i)),ichar(jmage(i))))
c        jmage(i) = char(ichar(image(i)).xor.ichar(jmage(i)))
   20    continue
c copy uncompressed image
         if (ifrmt.eq.3) then
            ibmax = lenb
            do 30 i = 1, ibmax
            img(i) = jmage(i)
   30       continue
c compress image
         else
            call comprs(jmage,img,line,limg,lz,ly,lz,lenb,leng,lzg,ibmax
     1)
         endif
         if (ibmax.gt.leng) write (6,91) ibmax, leng
c save old image
         do 40 i = 1, lenb
         jmage(i) = image(i)
   40    continue
      endif
c write length of compressed image to disk
      if (ifrmt.ne.4) then
         call convic(ibmax,chr,4,1)
         call buffwr(chr,4)
      endif
c write compressed image to disk
   50 call buffwr(img,ibmax)
      return
      end
      subroutine convic(lin,chr,len,n)
c this subroutine converts packed integer to character data
c lin = input array of integers
c chr = output array of characters
c len = number of characters output, should have len = n*lw
c n = number of integers to be converted
c lw = number of bytes per word
      parameter (lw=4)
      character*1 chr(len)
      dimension lin(n)
      dimension lb(lw)
      save nw
c nw = -2**(8*lw-1) + 1, used for adding and subtracting high bit
      data nw /-2147483647/
c exit if len < 1
      if (len.lt.1) return
      mw = nw - 1
c l = number of integers
      l = (len - 1)/lw + 1
c m = number of 4 byte integers
      m = len/lw
c mr = number of bytes left over
      mr = len - m*lw
      lwm = lw - 1
      lwp = lw + 1
c treat 4 byte integers
      if (m.gt.0) then
         do 30 i = 1, m
         i1 = (i - 1)*lw
         lb(1) = lin(i)
c subtract out high bit of first byte
         if (lb(1).lt.0) lb(1) = lb(1) - mw
         do 10 j = 1, lwm
         lb(j+1) = lb(j)/256
         lb(j) = lb(j) - lb(j+1)*256
   10    continue
c add high bit back in
         if (lin(i).lt.0) lb(lw) = lb(lw) + 128
         do 20 j = 1, lw
         chr(i1 + j) = char(lb(lwp - j))
   20    continue
   30    continue
      endif
c treat left over bytes
      if (mr.gt.0) then 
         i1 = m*lw
         lb(1) = lin(l)
c subtract out high bit of first byte
         if (lb(1).lt.0) lb(1) = lb(1) - mw
         do 40 j = 1, lwm
         lb(j+1) = lb(j)/256
         lb(j) = lb(j) - lb(j+1)*256
   40    continue
c add high bit back in
         if (lin(l).lt.0) lb(lw) = lb(lw) + 128
         do 50 j = 1, mr
         chr(i1 + j) = char(lb(lwp - j))
   50    continue
      endif
      return
      end
      subroutine comprs (image,img,line,limg,lenl,numl,lenv,lmax,lmaxg,l
     1enlg,ib)
c this subroutine compresses binary data, using the algorithm
c by r. h. frobose, jr., lawrence livermore lab report ucrl-51858
c written by viktor k. decyk, ucla
c input is in array image, of size lmax = lenv*numl, where
c lenv = spacing between rows in bytes, numl = number of lines
c lenl = length of line in bytes
c output is in array img, of maximum size lmaxg = lmax+(lmax-1)/3+1,
c and actual size given in variable ib
c line and limg are scratch arrays, of size lenl and lenlg respectively,
c where lenlg = lenl+(lenl-1)/3+1
c image = uncompressed input image
c img = compressed output image
c line, limg = uncompressed, compressed image line
c lenl = length of line in bytes
c numl = number of lines
c lenv = spacing between rows in bytes
c lmax = length of output image array = lenv*numl
c lmaxg = maximum size of img array = lmax+(lmax-1)/3+1,
c lenlg = maximum size of compressed line
c ib = actual size of img array
      character*1 image(lmax), img(lmaxg), line(lenl), limg(lenlg)
      integer ic0,ic1,ic2,ic3
      save icn,ic0,ic1,ic2,ic3,istyle
      data icn,ic0,ic1,ic2,ic3 /64,0,64,128,192/
c istyle = (0,1,2) = (no xor,xor,both) successive lines
      data istyle /2/
      icnr = icn/2 - 2
c prevent xor on first line
      ib = lmax
      ib0 = 0
      il = 0
      go to 120
   10 if (ixor.eq.1) go to 100
      if ((istyle.gt.0).and.(lib.gt.(ib-ib0))) go to 30
c xored line is longer
      do 20 j = 1, lib
      img(j+ib0) = limg(j)
   20 continue
      ib = ib0 + lib
c start new line
   30 il = il + lenv
      if (il.ge.lmax) go to 450
c looking for identical lines
      irl = 0
      il1 = il - lenv
      i = 0
   40 if (i.eq.lenl) go to 50
      i = i + 1
      if (image(i+il).eq.image(i+il1)) go to 40
c line not identical
      if (irl.eq.0) go to 60
c repeat previous line irl times.  different line follows.
      ib = ib + 1
      itc = ic3 + irl
      img(ib) = char(itc)
      go to 60
c line identical
   50 irl = irl + 1
      il = il + lenv
      i = 0
      if ((il.lt.lmax).and.(irl.lt.icnr)) go to 40
c repeat previous line irl times.  buffer full.
      ib = ib + 1
      itc = ic3 + irl
      img(ib) = char(itc)
      if (il.ge.lmax) go to 450
      irl = 0
      i = 0
      go to 40
c test whether to skip xor
   60 ib0 = ib
      if (istyle.eq.0) go to 120
c  xor current line
      il1 = il - lenv
      do 70 i = 1, lenl
      line(i) = char(ieor(ichar(image(i+il)),ichar(image(i+il1))))
c     line(i) = char(ichar(image(i+il)).xor.ichar(image(i+il1)))
   70 continue
      lib = 1
      itc = ic3 + icnr + 1
      limg(lib) = char(itc)
      ixor = 1
      go to 140
c second pass
c save previous compressed line
  100 do 110 j = 1, lib
      img(j+ib) = limg(j)
  110 continue
      ib = ib + lib
c test whether to perform xor
      if (istyle.eq.1) go to 30
  120 do 130 i = 1, lenl
      line(i) = image(i+il)
  130 continue
      lib = 0
      ixor = 0
  140 i = 1
      iz = 0
      id = 0
  150 if (line(i).eq.char(0)) go to 400
c non-zero bytes
c looking for different bytes
  200 id = id + 1
      if (i.eq.lenl) go to 260
      if (id.eq.icn) go to 240
      i = i + 1
      if (line(i).ne.char(0)) go to 205
      if ((i.eq.lenl).or.(line(i+1).ne.char(0))) go to 200
      go to 220
  205 if (line(i).ne.line(i-1)) go to 200
      if (id.eq.1) go to 300
      if ((i.eq.lenl).or.(line(i+1).ne.line(i))) go to 200
      id = id - 1
c plot next id bytes as they appear. pair follows.
      lib = lib + 1
      itc = ic2 + id - 1
      limg(lib) = char(itc)
      i1 = (i - id - 2)
      do 210 j = 1, id
      limg(j+lib) = line(j+i1)
  210 continue
      lib = lib + id
      go to 300
c plot next id bytes as they appear. zero follows.
  220 lib = lib + 1
      itc = ic2 + id - 1
      limg(lib) = char(itc)
      i1 = (i - id - 1)
      do 230 j = 1, id
      limg(j+lib) = line(j+i1)
  230 continue
      lib = lib + id
      id = 0
      go to 400
c plot next id bytes as they appear.  buffer full.
  240 lib = lib + 1
      itc = ic2 + id - 1
      limg(lib) = char(itc)
      i1 = (i - id)
      do 250 j = 1, id
      limg(j+lib) = line(j+i1)
  250 continue
      lib = lib + id
      i = i + 1
      id = 0
      go to 150
c plot next id bytes as they appear.  line full.
  260 lib = lib + 1
      itc = ic2 + id - 1
      limg(lib) = char(itc)
      i1 = (i - id)
      do 270 j = 1, id
      limg(j+lib) = line(j+i1)
  270 continue
      lib = lib + id
      go to 10
c non-zero bytes
c looking for identical bytes
  300 ir = 1
  310 ir = ir + 1
      if (i.eq.lenl) go to 330
      if (ir.eq.icn) go to 320
      i = i + 1
      if (line(i).eq.line(i-1)) go to 310
c repeat next byte ir times.  different byte follows.
      lib = lib + 1
      itc = ic0 + ir - 1
      limg(lib) = char(itc)
      lib = lib + 1
      limg(lib) = line(i-1)
      id = 0
      go to 150
c repeat next byte ir times.  buffer full.
  320 lib = lib + 1
      itc = ic0 + ir - 1
      limg(lib) = char(itc)
      lib = lib + 1
      limg(lib) = line(i)
      i = i + 1
      id = 0
      go to 150
c repeat next byte ir times.  line full.
  330 lib = lib + 1
      itc = ic0 + ir - 1
      limg(lib) = char(itc)
      lib = lib + 1
      limg(lib) = line(i)
      go to 10
c zero bytes
  400 iz = iz + 1
      if (i.eq.lenl) go to 430
      i = i + 1
      if (line(i).eq.char(0)) go to 400
  410 if (iz.le.icn) go to 420
c skip next iz bytes, and plot the following byte.  buffer full.
      lib = lib + 1
      itc = ic1 + icn - 1
      limg(lib) = char(itc)
      lib = lib + 1
      limg(lib) = char(0)
      iz = iz - (icn + 1)
      if (iz.eq.0) go to 200
      go to 410
c skip next iz bytes, and plot the following byte.  next byte non-zero.
  420 lib = lib + 1
      itc = ic1 + iz - 1
      limg(lib) = char(itc)
      lib = lib + 1
      limg(lib) = line(i)
      if (i.eq.lenl) go to 10
      i = i + 1
      iz = 0
      go to 150
c plot end-of-line sentinel
  430 lib = lib + 1
      limg(lib) = char(ic3)
      go to 10
c end-of-frame byte
  450 ib = ib + 1
      img(ib) = char(0)
      return
      end
      subroutine buffwr(line,n)
c this subroutine packs image data into buffer and writes when full
c input arguments: all
c line = input character array to be written
c n = number of characters to be packed (0 = flush buffer)
c lmax = maximum number of characters in lout buffer
      parameter(lmax=80)
      character*1 c0
      character*1 line(*)
      character*1 lout(lmax)
      save len,lout
c len = number of characters currently stored in lout buffer
      data len /0/
c nc = number of characters to be written to lout
      nc = n
c ncr = number of characters remaining to be written
      ncr = n
      i = 0
c flush buffer requested
      if (n.eq.0) then
         if (len.ge.0) go to 40
c buffer already flushed
         return
      endif
c calculate position pointers in character arrays
      ncr = (nc + len) - lmax
      if (ncr.ge.0) nc = lmax - len
c pack output buffer with special null character
   10 if (len.eq.0) then
         c0 = char(232)
         do 20 j = 1, lmax
         lout(j) = c0
   20    continue
      endif
c pack input array into output buffer
      do 30 j = 1, nc
      lout(j+len) = line(i+j)
   30 continue
c update ofset in input array
      i = i + nc
c update last position written in lout buffer
      len = len + nc
c no more characters left to write
      if (ncr.lt.0) return
c flush buffer if full
   40 if (len.gt.0) call wrpct(lout,lmax)
c reset last position writteen in lout buffer
      len = 0
c force flush buffer
      if (n.eq.0) call wrpct(lout,len)
c more characters remaining to be written
      if (ncr.gt.0) then
c update position pointers in character arrays
         nc = ncr
         ncr = ncr - lmax
         if (ncr.ge.0) nc = lmax
         go to 10
      endif
      return
      end
      subroutine wrpct(lout,len)
c this subroutine writes compressed raster file to file
c input arguments: all
c lout = character data to write
c len = number of characters to write
c mtype = machine type
c 1=rs/6000, 2=sun, 3=cray c90, 4=ibm es/9000, 5=vax, 6=dec, 7=hp, 8=sgi
c 9=ibm pc, 10=mac, 11=paragon, 12=cray t3d, 13=fujitsu vpp500
      common /march/ iebc, irvb, longi, mtype
      character*1 lout(*)
      character*8 cs
c set format string
      cs = '(80a1,$)'
c special case for ibm es/9000
      if (mtype.eq.4) cs = '(80a1)  '
c special case for macintosh
      if (mtype.eq.10) then
         if (len.gt.0) write (20) (lout(j),j=1,len)
      else
         if (len.gt.0) write (20,cs) (lout(j),j=1,len)
      endif
      return
      end
c internal library for printer plot images
      subroutine ppc(ict)
c this subroutine converts color indices to printable characters
c input and output argument: ict
c ict = color index to be converted
      character*64 pp
      save pp
c pp = printable character conversion table
      data pp /' *abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQR
     1STUVWXYZ'/
      it = ict + 1
      if (it.lt.1) it = 1
      if (it.gt.64) it = 64
      ict = ichar(pp(it:it))
      return
      end
      subroutine cwrite(chr,ix,iy)
c this subroutine writes character string chr at location (ix,iy)
c input: all
c chr = input string to be written
c ix, iy = x, y location to write character
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
c icx, icy = current location of cursor
c g = current image array
      character*(*) chr
      n = len(chr)
c clip length of string, if necessary
      nn = min0(ix+n,lx) - ix
      joff = ix + lx*iy
c write string
      do 10 j = 1, nn
      g(j+joff) = chr(j:j)
   10 continue
c save last location
      icx = ix + nn - 1
      icy = iy + nn - 1
      return
      end
      subroutine wimage
c this subroutine writes image g to terminal
c lxm, lym = maximum size of g image array
      parameter(lxm=720,lym=540)
      common /grsdev/ lx,ly,nbit,intrl,ixor,ipmx,icx,icy,g
      character*1 g(lym*(lxm+(lxm-1)/128+1))
c lx, ly = the size of the current image
c g = current image array
c nfc = (0,1) = (no,yes) use fortran style control
      data nfc /0/
      save nfc
   91 format ('1',128a1)
   92 format (1x,128a1)
      do 10 k = 1, ly
      joff = lx*(ly - k)
      if ((nfc.eq.1).and.(k.eq.1)) then
         write (6,91) (g(j+joff),j=1,lx)
      else
         write (6,92) (g(j+joff),j=1,lx)
      endif
   10 continue
      return
      end
c internal library for tiff raster images
      subroutine tfhead
c this subroutines write Image File Header (IFH) for tiff files
      character*1 chr(4)
c clear all entries
      do 10 j = 1, 4
      chr(j) = char(0)
   10 continue
c code for Big-Endian, location first IDF not written
      chr(1) = char(77)
      chr(2) = char(77)
      chr(4) = char(42)
      call buffwr(chr,4)
      return
      end
      subroutine timage (image,img,lz,ly,nbit,lenb,leng)
c this subroutine writes image to tiff file with header
c image = reduced image array
c img = compressed output image
c lz, ly = the size of the reduced image
c nbit = number of bits per pixel, nbit < 9
c lenb = number of characters in image
c leng = maximum length of compressed image
      character*1 image(lenb), img(leng)
c nifd = length of image file directory
c maxs = maximum number of strips
      parameter(nifd=170,maxs=24)
c chr = Image File Directory (IFD) record
      character*1 chr(nifd)
c ioff, noff = strip location and size
      dimension ioff(maxs), noff(maxs)
      save last
c last = location of last write
c kmprs = (0,1) = (no,yes) compress file
      data last,kmprs /4,1/
c largest possible number of strips in image
      nstrips = min(3*nbit,ly)
c nrps = maximum number of rows/strip
      nrps = (ly - 1)/nstrips + 1
c actual maximum number of strips in image
      nstrips = (ly - 1)/nrps + 1
c imgloc = location of image data
      imgloc = last + 4
c calculate size of compressed image
      nrow = nrps
      lrow = nrps*lz
      lost = 0
      lrows = 1
      do 10 j = 1, nstrips
      if (j.eq.nstrips) then
         nrow = ly -  nrps*(nstrips - 1)
         lrow = nrow*lz
      endif
      ioff(j) = imgloc + lost
c no compression
      if (kmprs.eq.0) then
         ib = lrow
c PackBits compression
      else
         call pkbits (image(lrows),img(1+lost),lz,nrow,lz,lrow,leng,ib)
      endif
      noff(j) = ib
      lost = lost + ib
      lrows = lrows + lrow
   10 continue
c ifdloc = location of IFD data
      ifdloc = imgloc + lost
c update location of last write to end of IFD data
      last = ifdloc + (16 + 6*nstrips)
      if (nbit.eq.8) last = last + 1536
c convert offset to character
      call convic(last,chr,4,1)
c write offset to IFD record
      call buffwr(chr,4)
c write image
c no compression
      if (kmprs.eq.0) then
         call buffwr(image,lost)
c PackBits compression
      else
         call buffwr(img,lost)
      endif
c write IDF data
c clear 16 entries in character buffer
      do 20 j = 1, 16
      chr(j) = char(0)
   20 continue
c XResolution Data: 72, 1
      chr(4) = char(72)
      chr(8) = char(1)
c XResolution Data: 72, 1
      chr(12) = char(72)
      chr(16) = char(1)
c write Resolution data
      call buffwr(chr,16)
c write StripOffsets
      do 30 j = 1, nstrips
      call convic(ioff(j),chr,4,1)
      call buffwr(chr,4)
   30 continue
c write StripByteCounts
      do 40 j = 1, nstrips
      lh = noff(j)/256
      chr(1) = char(lh)
      chr(2) = char(noff(j) - lh*256)
      call buffwr(chr,2)
   40 continue
c write Palette
      if (nbit.eq.8) call swpal(256,0.,0.,0.)
c write IFD record
c clear all entries
      do 50 j = 1, nifd
      chr(j) = char(0)
   50 continue
c number of tags in IFD
c 13 if monochrome image
      if (nbit.eq.1) then
         ntags = 13
c 14 if color palette image
      else
         ntags = 14
      endif
c lidf = length of IDF record without last pointer
      lidf = 2 + 12*ntags
c update location of last write to end of pointerless idf record
      last = last + lidf
      chr(2) = char(ntags)
c set default values for TifTags
c TadId = 256, DataType = SHORT, and DataCount = 1
      do 60 j = 1, ntags
      joff = 12*(j - 1)
      chr(3+joff) = char(1)
      chr(6+joff) = char(3)
      chr(10+joff) = char(1)
   60 continue
c NewSubfiletype Tag: multi-page image
      chr(3) = char(0)
      chr(4) = char(254)
      chr(6) = char(4)
      chr(14) = char(2)
c ImageWidth Tag: lx
      lx = (8*lz)/nbit
      lh = lx/256
      chr(23) = char(lh)
      chr(24) = char(lx - lh*256)
c ImageLength Tag: ly
      chr(28) = char(1)
      lh = ly/256
      chr(35) = char(lh)
      chr(36) = char(ly - lh*256)
c BitsPerSample: nbit
      chr(40) = char(2)
      chr(48) = char(nbit)
c Compression
      chr(52) = char(3)
c no compression
      if (kmprs.eq.0) then
         chr(60) = char(1)
c PackBits compression
      else
         chr(59) = char(128)
         chr(60) = char(5)
      endif
c PhotometricInterpretation
      chr(64) = char(6)
c BlackisZero
      if (nbit.eq.1) then
         chr(72) = char(1)
C Palette color
      else
         chr(72) = char(3)
      endif
c StripOffsets
      chr(76) = char(17)
      chr(78) = char(4)
      chr(82) = char(nstrips)
c store location of last write to end of Resolution data
      lost = ifdloc + 16
c convert offset to character
      call convic(lost,chr(83),4,1)
c SamplesPerPixel: one
      chr(88) = char(21)
      chr(96) = char(1)
c RowsPerStrip
      chr(100) = char(22)
      chr(102) = char(4)
      lh = nrps/256
      chr(109) = char(lh)
      chr(110) = char(nrps - lh*256)
c StripByteCounts
      chr(112) = char(23)
      chr(118) = char(nstrips)
c store location of last write to end of StripOffsets record
      lost = lost + 4*nstrips
c convert offset to character
      call convic(lost,chr(119),4,1)
c XResolution
      chr(124) = char(26)
      chr(126) = char(5)
c convert offset to character
      call convic(ifdloc,chr(131),4,1)
c YResolution
      chr(136) = char(27)
      chr(138) = char(5)
c update location of last write to end of XResolution Data
      lost = ifdloc + 8
c convert offset to character
      call convic(lost,chr(143),4,1)
c ResolutionUnit: inch
      chr(148) = char(40)
      chr(156) = char(2)
c ColorMap
      if (nbit.eq.8) then
         chr(160) = char(64)
         chr(165) = char(3)
         chr(166) = char(0)
c update location of last write to end of StripByteCounts Data
         lost = ifdloc + (16 + 6*nstrips)
c convert offset to character
         call convic(lost,chr(167),4,1)
      endif
c write IDF record without next pointer
      call buffwr(chr,lidf)
      return
      end
      subroutine pkbits (image,img,lenl,numl,lenv,lmax,lmaxg,ib)
c this subroutine compresses binary data, using the algorithm
c known as Apple Macintosh PackBits
c written by viktor k. decyk, ucla
c input is in array image, of size lmax = lenv*numl, where
c lenv = spacing between rows in bytes, numl = number of lines
c lenl = length of line in bytes
c output is in array img, of maximum size lmaxg
c and actual size given in variable ib
c image = uncompressed input image
c img = compressed output image
c lenl = length of line in bytes
c numl = number of lines
c lenv = spacing between rows in bytes
c lmax = length of input image array = lenv*numl
c lmaxg = maximum size of img array = lmax + numl*((lenv-1)/128+1)
c ib = actual size of img array
      character*1 image(lmax), img(lmaxg)
      integer ic0,ic2
      save icn,ic0,ic2
      data icn,ic0,ic2 /128,256,0/
      ib = 0
      il = -lenv
c start new line
   10 il = il + lenv
      if (il.ge.lmax) go to 120
      linl = il + lenl
      i = il + 1
      id = 0
c looking for different bytes
   20 id = id + 1
      if (i.eq.linl) go to 60
      if (id.eq.icn) go to 40
      i = i + 1
      if (i.eq.linl) go to 20
      if (image(i).ne.image(i-1)) go to 20
      if (id.eq.1) go to 80
      if ((i.eq.linl).or.(image(i+1).ne.image(i))) go to 20
      id = id - 1
c plot next id bytes as they appear. pair follows.
      ib = ib + 1
      itc = ic2 + (id - 1)
      img(ib) = char(itc)
      i1 = (i - id - 2)
      do 30 j = 1, id
      img(j+ib) = image(j+i1)
   30 continue
      ib = ib + id
      go to 80
c plot next id bytes as they appear.  buffer full.
   40 ib = ib + 1
      itc = ic2 + (id - 1)
      img(ib) = char(itc)
      i1 = (i - id)
      do 50 j = 1, id
      img(j+ib) = image(j+i1)
   50 continue
      ib = ib + id
      i = i + 1
      id = 0
      go to 20
c plot next id bytes as they appear.  line full.
   60 ib = ib + 1
      itc = ic2 + (id - 1)
      img(ib) = char(itc)
      i1 = (i - id)
      do 70 j = 1, id
      img(j+ib) = image(j+i1)
   70 continue
      ib = ib + id
      go to 10
c looking for identical bytes
   80 ir = 1
   90 ir = ir + 1
      if (i.eq.linl) go to 110
      if (ir.eq.icn) go to 100
      i = i + 1
      if (image(i).eq.image(i-1)) go to 90
c repeat next byte ir times.  different byte follows.
      ib = ib + 1
      itc = ic0 - (ir - 1)
      img(ib) = char(itc)
      ib = ib + 1
      img(ib) = image(i-1)
      id = 0
      go to 20
c repeat next byte ir times.  buffer full.
  100 ib = ib + 1
      itc = ic0 - (ir - 1)
      img(ib) = char(itc)
      ib = ib + 1
      img(ib) = image(i)
      i = i + 1
      id = 0
      go to 20
c repeat next byte ir times.  line full.
  110 ib = ib + 1
      itc = ic0 - (ir - 1)
      img(ib) = char(itc)
      ib = ib + 1
      img(ib) = image(i)
      go to 10
c check if compression overflowed img array
  120 if (ib.gt.lmaxg) write (6,*) 'overflow: ib,lmaxg=',ib,lmaxg
      return
      end
      subroutine wrtpal(pal,npal,lpald,ipmx)
c this subroutine writes palette for tiff format files
c pal is a character array, with rgb values in successive bytes
c pal(i+3*j) = ith component of jth palette entry
c npal = (0,n) = (default,n) palette entries
c lpald = size of palette array
c ipmx =  maximum color value in the palette
      character*1 pal(lpald)
c color = scratch array for 8 entries
      character*1 color(48)
c icolor = integer array to define 8 primary colors
      dimension icolor(24)
      save icolor
c data for 8 primary colors
      data icolor /0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,0,1,0,1,0,1/
      ipm = 65535
c write default palette
      if (npal.le.0) then
c first 8 entries define primary default colors
         do 10 i = 1, 24
         ip = ipm*icolor(i)
         lh = ip/256
         color(2*i-1) = char(lh)
         color(2*i) = char(ip - lh*256)
   10    continue
         call buffwr(color,48)
c pad remaining 248 entries with foreground (white)
         do 20 i = 1, 48
         color(i) = char(255)
   20    continue
         do 30 i = 1, 31
         call buffwr(color,48)
   30    continue
c write specified palette
      else
         len = npal
         if (len.gt.256) len = 256
         lb = (len - 1)/24 + 1
         nbf = len - 24*(lb - 1)
         apmx = 65535./float(ipmx - 1)
c cycle through red, green, and blue
         do 60 l = 1, 3
         nb = 24
c cycle through palette in blocks of 24
         do 50 i = 1, lb
         ioff = 24*(i - 1) - 1
         if (i.eq.lb) nb = nbf
         do 40 ii = 1, nb
         ip = apmx*float(ichar(pal(l+3*(ii+ioff)))) + .5
         lh = ip/256
         color(2*ii-1) = char(lh)
         color(2*ii) = char(ip - lh*256)
   40    continue
         call buffwr(color,2*nb)
   50    continue
   60    continue
c pad remaining entries in palette with foreground (white)
         if (len.lt.256) then
            nl = (768 - 3*len)
c n = number of full groups of 8 entries
            n = nl/24
c nl = remaining number of characters to write
            nl = nl - 24*n
c pad character array
            do 70 i = 1, 48
            color(i) = char(ipm)
   70       continue
c write full groups of 8 entries
            if (n.gt.0) then
               do 80 i = 1, n
               call buffwr(color,48)
   80          continue
            endif
c write remainder
            if (nl.gt.0) call buffwr(color,2*nl)
         endif
      endif
      return
      end
      integer function kywait()
c special function to request keystroke
c returns: ascii code for keystroke
      parameter(maxt=3)
      common /gksrstr/ kosv,iwk,idc,iwt,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf
     1,ltc,mtc,infs,itxtp,lws,icc,ams,chx,cch,wwnvp,dvt,cea,trans,chuv
      dimension icc(5), wwnvp(8), dvt(4), cea(4), trans(8,maxt), chuv(2)
c iwt = current workstation type
      character*1 str
c fortran input supported
      if (iwt.eq.8) then
         read (5,*) str
         kywait = ichar(str)
      else
         kywait = 0
      endif
      return
      end
