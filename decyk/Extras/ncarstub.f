      subroutine gitwks(idcon,iwtype)
c this is a site-dependent subroutine which returns
c connection identifier and workstation type
c version for ncargraphics, version 3.2
c iwtype = workstation type
c 1 = CGM, 3 = WISS, 7 = pre-existing X11 window,
c 8 = non-existing X11 window, 10 = text dump of graphics output
      iwtype = 1
c idcon = connection identifier, 
c except for iwtype = 7, this is a fortran unit number or is not needed
      idcon = 1
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
c dummy version for gks libraries
      character*1 image(lenb)
      call fimage(image,lx,ly,lz,lenb,npix,nf)
      return
      end
c gks null device driver for missing subroutines in ncar graphics
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: february 25, 1995
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
      character*(*) str
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 0
      lstr = 1
      ipet = 1
      earea(1) = 0.
      earea(2) = 1.
      earea(3) = 0.
      earea(4) = 1.
      lenb = 1
      ipos = 1
      ldr = 1
      ierr = 140
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
      subroutine gscnid(idcon,connam)
c set connection identifier
c input arguments: all
c idcon = connection identifier
c connam = connection name
      character*8 connam
      return
      end
