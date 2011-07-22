c-----------------------------------------------------------------------
c 1d PIC library for fast fourier transforms
c fft1lib.f contains procedures to perform ffts:
c FFT1RX performs real to complex fft and its inverse for scalar array.
c FFT1R2 performs real to complex fft and its inverse for 2 component
c        vector array.
c FFT1R3 performs real to complex fft and its inverse for 3 component
c        vector array.
c FFT1C performs complex to complex fft and its inverse.
c FST1RX performs fast sine transform and its inverse for scalar array.
c FCT1RX performs fast cosine transform and its inverse for scalar array
c FST1R2 performs fast sine transform and its inverse for 2 component
c        vector array.
c FCT1R2 performs fast cosine transform and its inverse for 2 component
c        vector array.
c FCST1R3 performs fast cosine and sine transforms and its inverse for
c         3 component vector array.
c FDST1RX performs fast sine DST-III transform and its inverse.
c FDCT1RX performs fast cosine DCT-III transform and its inverse.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: july 8, 2010
c-----------------------------------------------------------------------
      subroutine FFT1RX(f,t,isign,mixup,sct,indx,nxd,nxhd)
c this subroutine performs a one dimensional real to complex fast
c fourier transform and its inverse, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(n) = (1/nx)*sum(f(j)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(j) = sum(f(n)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c fourier coefficients are stored as follows:
c f(1) = real part of mode 0, f(2) = real part of mode nx/2
c f(2*j-1),f(2*j) = real,imaginary part of mode j-1, 0 < j < nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(nxd), t(nxhd), mixup(nxhd), sct(nxhd)
c local data
      integer nx, nxh, nxhh, nxh2, indx1, j, k, l, lb, ll, jb, it, j1
      integer nxs, nxs2, km, km2, k1, k2
      real dnx, arg, ani
      complex t1, t2
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign) 50,10,120
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxh
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nx
      dnx = 6.28318530717959/real(nx)
      do 40 j = 1, nxh
      arg = dnx*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
c bit-reverse array elements to complex temporary
   50 do 60 j = 1, nxh
      j1 = mixup(j)
      t(j) = cmplx(f(2*j1-1),f(2*j1))
   60 continue
c reduction
      do 90 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 80 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 70 j = 1, nxs
      t1 = sct(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
   70 continue
   80 continue
   90 continue
c unscramble coefficients and normalize result
      ani = 1./real(2*nx)
      do 100 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),-real(sct(j)))
      t(j) = ani*(t1 + t2)
      t(nxh2-j) = ani*conjg(t1 - t2)
  100 continue
      ani = 2.*ani
      t(nxhh+1) = ani*conjg(t(nxhh+1))
      t(1) = ani*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1))
     1)
c move to real destination
      do 110 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
  110 continue
      return
c forward fourier transform
c move to complex temporary
  120 do 130 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
  130 continue
c scramble coefficients
      do 140 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sct(j)),real(sct(j)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
  140 continue
      t(nxhh+1) = 2.*conjg(t(nxhh+1))
      t(1) = cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1)))
c bit-reverse array elements to real destination
      do 150 j = 1, nxh
      j1 = mixup(j)
      f(2*j-1) = real(t(j1))
      f(2*j) = aimag(t(j1))
  150 continue
c move back to complex temporary
      do 160 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
  160 continue
c reduction
      do 190 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 180 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 170 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  170 continue
  180 continue
  190 continue
c move to real destination
      do 200 j = 1, nxh
      f(2*j-1) = real(t(j))
      f(2*j) = aimag(t(j))
  200 continue
      return
      end
      subroutine FFT1R2(f,t,isign,mixup,sct,indx,nxd,nxhd)
c this subroutine performs two one dimensional real to complex fast
c fourier transforms and their inverses, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(1:2,n) = (1/nx)*sum(f(1:2,j)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(1:2,j) = sum(f(1:2,n)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c fourier coefficients are stored as follows:
c f(1) = real part of mode 0, f(2) = real part of mode nx/2
c f(2*j-1),f(2*j) = real,imaginary part of mode j-1, 0 < j < nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(2,nxd), t(2,nxhd), mixup(nxhd), sct(nxhd)
c local data
      integer nx, nxh, nxhh, nxh2, indx1, j, k, l, lb, ll, jb, it, j1
      integer nxs, nxs2, km, km2, k1, k2, jj
      real dnx, arg, ani
      complex t1, t2, t3
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign) 50,10,140
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxh
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nx
      dnx = 6.28318530717959/real(nx)
      do 40 j = 1, nxh
      arg = dnx*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
c bit-reverse array elements to complex temporary
   50 do 60 j = 1, nxh
      j1 = mixup(j)
      t(1,j) = cmplx(f(1,2*j1-1),f(1,2*j1))
      t(2,j) = cmplx(f(2,2*j1-1),f(2,2*j1))
   60 continue
c reduction
      do 90 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 80 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 70 j = 1, nxs
      t1 = sct(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
   70 continue
   80 continue
   90 continue
c unscramble coefficients and normalize result
      ani = 1./real(2*nx)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),-real(sct(j)))
      do 100 jj = 1, 2
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = ani*(t1 + t2)
      t(jj,nxh2-j) = ani*conjg(t1 - t2)
  100 continue
  110 continue
      ani = 2.*ani
      do 120 jj = 1, 2
      t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
      t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) -
     1 aimag(t(jj,1)))
  120 continue
c move to real destination
c     do 130 j = 1, nxh
c     f(1,2*j-1) = real(t(1,j))
c     f(1,2*j) = aimag(t(1,j))
c     f(2,2*j-1) = real(t(2,j))
c     f(2,2*j) = aimag(t(2,j))
c 130 continue
c move to complex destination
      do 130 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(2,2*j-1) = aimag(t(1,j))
      f(1,2*j) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
  130 continue
      return
c forward fourier transform
c move real source to complex temporary
c 140 do 150 j = 1, nxh
c     t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
c     t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
c 150 continue
c move complex source to complex temporary
  140 do 150 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(2,2*j-1))
      t(2,j) = cmplx(f(1,2*j),f(2,2*j))
  150 continue
c scramble coefficients
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),real(sct(j)))
      do 160 jj = 1, 2
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = t1 + t2
      t(jj,nxh2-j) = conjg(t1 - t2)
  160 continue
  170 continue
      do 180 jj = 1, 2
      t(jj,nxhh+1) = 2.*conjg(t(jj,nxhh+1))
      t(jj,1) = cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) - aim
     1ag(t(jj,1)))
  180 continue
c bit-reverse array elements to real destination
      do 190 j = 1, nxh
      j1 = mixup(j)
      f(1,2*j-1) = real(t(1,j1))
      f(1,2*j) = aimag(t(1,j1))
      f(2,2*j-1) = real(t(2,j1))
      f(2,2*j) = aimag(t(2,j1))
  190 continue
c move back to complex temporary
      do 200 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
      t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
  200 continue
c reduction
      do 230 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 220 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 210 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
  210 continue
  220 continue
  230 continue
c move to real destination
      do 240 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = aimag(t(1,j))
      f(2,2*j-1) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
  240 continue
      return
      end
      subroutine FFT1R3(f,t,isign,mixup,sct,indx,nxd,nxhd)
c this subroutine performs three one dimensional real to complex fast
c fourier transforms and their inverses, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(1:3,n) = (1/nx)*sum(f(1:3,j)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(1:3,j) = sum(f(1:3,n)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c fourier coefficients are stored as follows:
c f(1) = real part of mode 0, f(2) = real part of mode nx/2
c f(2*j-1),f(2*j) = real,imaginary part of mode j-1, 0 < j < nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sct
      dimension f(3,nxd), t(3,nxhd), mixup(nxhd), sct(nxhd)
c local data
      integer nx, nxh, nxhh, nxh2, indx1, j, k, l, lb, ll, jb, it, j1
      integer nxs, nxs2, km, km2, k1, k2, jj
      real dnx, arg, ani
      complex t1, t2, t3, t4
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      indx1 = indx - 1
      if (isign) 50,10,140
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxh
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nx
      dnx = 6.28318530717959/real(nx)
      do 40 j = 1, nxh
      arg = dnx*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
c bit-reverse array elements to complex temporary
   50 do 60 j = 1, nxh
      j1 = mixup(j)
      t(1,j) = cmplx(f(1,2*j1-1),f(1,2*j1))
      t(2,j) = cmplx(f(2,2*j1-1),f(2,2*j1))
      t(3,j) = cmplx(f(3,2*j1-1),f(3,2*j1))
   60 continue
c reduction
      do 90 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 80 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 70 j = 1, nxs
      t1 = sct(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t4 = t1*t(3,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(3,j+k2) = t(3,j+k1) - t4
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
      t(3,j+k1) = t(3,j+k1) + t4
   70 continue
   80 continue
   90 continue
c unscramble coefficients and normalize result
      ani = 1./real(2*nx)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),-real(sct(j)))
      do 100 jj = 1, 3
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = ani*(t1 + t2)
      t(jj,nxh2-j) = ani*conjg(t1 - t2)
  100 continue
  110 continue
      ani = 2.*ani
      do 120 jj = 1, 3
      t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
      t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) -
     1 aimag(t(jj,1)))
  120 continue
c move to real destination
c     do 130 j = 1, nxh
c     f(1,2*j-1) = real(t(1,j))
c     f(1,2*j) = aimag(t(1,j))
c     f(2,2*j-1) = real(t(2,j))
c     f(2,2*j) = aimag(t(2,j))
c     f(3,2*j-1) = real(t(3,j))
c     f(3,2*j) = aimag(t(3,j))
c 130 continue
c move to complex destination
      do 130 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(2,2*j-1) = aimag(t(1,j))
      f(3,2*j-1) = real(t(2,j))
      f(1,2*j) = aimag(t(2,j))
      f(2,2*j) = real(t(3,j))
      f(3,2*j) = aimag(t(3,j))
  130 continue
      return
c forward fourier transform
c move real source to complex temporary
c 140 do 150 j = 1, nxh
c     t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
c     t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
c     t(3,j) = cmplx(f(3,2*j-1),f(3,2*j))
c 150 continue
c move complex source to complex temporary
  140 do 150 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(2,2*j-1))
      t(2,j) = cmplx(f(3,2*j-1),f(1,2*j))
      t(3,j) = cmplx(f(2,2*j),f(3,2*j))
  150 continue
c scramble coefficients
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(j)),real(sct(j)))
      do 160 jj = 1, 3
      t2 = conjg(t(jj,nxh2-j))
      t1 = t(jj,j) + t2
      t2 = (t(jj,j) - t2)*t3
      t(jj,j) = t1 + t2
      t(jj,nxh2-j) = conjg(t1 - t2)
  160 continue
  170 continue
      do 180 jj = 1, 3
      t(jj,nxhh+1) = 2.*conjg(t(jj,nxhh+1))
      t(jj,1) = cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)) - aim
     1ag(t(jj,1)))
  180 continue
c bit-reverse array elements to real destination
      do 190 j = 1, nxh
      j1 = mixup(j)
      f(1,2*j-1) = real(t(1,j1))
      f(1,2*j) = aimag(t(1,j1))
      f(2,2*j-1) = real(t(2,j1))
      f(2,2*j) = aimag(t(2,j1))
      f(3,2*j-1) = real(t(3,j1))
      f(3,2*j) = aimag(t(3,j1))
  190 continue
c move back to complex temporary
      do 200 j = 1, nxh
      t(1,j) = cmplx(f(1,2*j-1),f(1,2*j))
      t(2,j) = cmplx(f(2,2*j-1),f(2,2*j))
      t(3,j) = cmplx(f(3,2*j-1),f(3,2*j))
  200 continue
c reduction
      do 230 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = km + km
      do 220 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 210 j = 1, nxs
      t1 = conjg(sct(1+km2*(j-1)))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t4 = t1*t(3,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(3,j+k2) = t(3,j+k1) - t4
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
      t(3,j+k1) = t(3,j+k1) + t4
  210 continue
  220 continue
  230 continue
c move to real destination
      do 240 j = 1, nxh
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = aimag(t(1,j))
      f(2,2*j-1) = real(t(2,j))
      f(2,2*j) = aimag(t(2,j))
      f(3,2*j-1) = real(t(3,j))
      f(3,2*j) = aimag(t(3,j))
  240 continue
      return
      end
      subroutine FFT1C(f,isign,mixup,sct,indx,nxd,nxhd)
c this subroutine performs a one dimensional complex to complex fast
c fourier transform and its inverse, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = (-1,1), approximate flop count: 5*nx*log2(nx)
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(n) = (1/nx)*sum(f(j)*exp(-sqrt(-1)*2pi*n*j/nx))
c if isign = 1, a forward fourier transform is performed
c f(j) = sum(f(n)*exp(sqrt(-1)*2pi*n*j/nx))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxd = dimension of f
c nxhd = nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      complex f, sct
      dimension f(nxd), mixup(nxd), sct(nxhd)
c local data
      integer nx, nxh, j, k, l, lb, ll, jb, it, j1
      integer nxs, nxs2, km, k1, k2
      real dnx, arg, ani
      complex t
      nx = 2**indx
      nxh = nx/2
      if (isign.ne.0) go to 40
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nx
      lb = j - 1
      ll = 0
      do 10 k = 1, indx
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nx
      dnx = 6.28318530717959/real(nx)
      do 30 j = 1, nxh
      arg = dnx*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
c bit-reverse array elements
   40 do 50 j = 1, nx
      j1 = mixup(j)
      if (j.ge.j1) go to 50
      t = f(j1)
      f(j1) = f(j)
      f(j) = t
   50 continue
      if (isign.gt.0) go to 100
c inverse fourier transform
      do 80 l = 1, indx
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxh/nxs
      do 70 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 60 j = 1, nxs
      t = sct(1+km*(j-1))*f(j+k2)
      f(j+k2) = f(j+k1) - t
      f(j+k1) = f(j+k1) + t
   60 continue
   70 continue
   80 continue
c normalize result
      ani = 1./real(nx)
      do 90 j = 1, nx
      f(j) = f(j)*ani
   90 continue
      return
c forward fourier transform
  100 do 130 l = 1, indx
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxh/nxs
      do 120 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 110 j = 1, nxs
      t = conjg(sct(1+km*(j-1)))*f(j+k2)
      f(j+k2) = f(j+k1) - t
      f(j+k1) = f(j+k1) + t
  110 continue
  120 continue
  130 continue
      return
      end
      subroutine FST1RX(f,t,isign,mixup,sctd,indx,nxd,nxhd)
c this subroutine performs a one dimensional fast real sine transform
c and its inverse, using complex arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = 0, input: all except f, output: mixup, sctd
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 23), where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse sine transform is performed
c f(n) = (1/nx)*sum(f(j)*sin(pi*n*j/nx))
c if isign = 1, a forward sine transform is performed
c f(j) = 2*sum(f(n)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxd >= nx + 1
c nxhd >= nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sctd
      dimension f(nxd), t(nxhd), mixup(nxhd), sctd(nxd)
c local data
      integer nx, nxh, nxhh, nxh2, nx2, indx1, j, k, l, jb, it, lb, ll
      integer j1, nxs, nxs2, km, km2, k1, k2
      real dnx, arg, ani, at1, at2, at3, at4, sum1
      complex t1, t2
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      nx2 = nx + 2
      indx1 = indx - 1
      if (isign.ne.0) go to 40
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/real(nx+nx)
      do 30 j = 1, nx
      arg = dnx*real(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
c create auxiliary array and bit-reverse to complex temporary
   40 do 50 j = 2, nxh
      j1 = mixup(j)
      at2 = -aimag(sctd(2*j1-1))
      at1 = (at2 + .5)*f(2*j1-1)
      at2 = (at2 - .5)*f(nx2-2*j1+1)
      at4 = -aimag(sctd(2*j1))
      at3 = (at4 + .5)*f(2*j1)
      at4 = (at4 - .5)*f(nx2-2*j1)
      t(j) = cmplx(at1+at2,at3+at4)
   50 continue
      at4 = -aimag(sctd(2))
      at3 = (at4 + .5)*f(2)
      at4 = (at4 - .5)*f(nx2-2)
      t(1) = cmplx(0.,at3+at4)
c reduction
      do 80 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 4*km
      do 70 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 60 j = 1, nxs
      t1 = sctd(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize result
c inverse fourier transform
      if (isign.lt.0) then
         ani = 1./real(2*nx)
         do 90 j = 2, nxhh
         t2 = conjg(t(nxh2-j))
         t1 = t(j) + t2
         t2 = (t(j) - t2)*cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         t(j) = ani*(t1 + t2)
         t(nxh2-j) = ani*conjg(t1 - t2)
   90    continue
         ani = 2.*ani
         t(nxhh+1) = ani*conjg(t(nxhh+1))
         t(1) = ani*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(
     11)))
c forward fourier transform
      else if (isign.gt.0) then
         do 100 j = 2, nxhh
         t2 = conjg(t(nxh2-j))
         t1 = t(j) + t2
         t2 = (t(j) - t2)*cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         t(j) = t1 + t2
         t(nxh2-j) = conjg(t1 - t2)
  100    continue
         t(nxhh+1) = 2.0*conjg(t(nxhh+1))
         t(1) = 2.0*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(
     11)))
      endif
c perform recursion and move to real destination
      sum1 = .5*real(t(1))
      f(1) = 0.0
      f(2) = sum1
      do 110 j = 2, nxh
      sum1 = sum1 + real(t(j))
      f(2*j-1) = -aimag(t(j))
      f(2*j) = sum1
  110 continue
      f(nx+1) = 0.0
      return
      end
      subroutine FCT1RX(f,t,isign,mixup,sctd,indx,nxd,nxhd)
c this subroutine performs a one dimensional fast real cosine transform
c and its inverse, using complex arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = 0, input: all except f, output: mixup, sctd
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 23), where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse cosine transform is performed
c f(n) = (1/nx)*(.5*f(1) + ((-1)**n)*f(nx+1) + sum(f(j)*cos(pi*n*j/nx)))
c if isign = 1, a forward cosine transform is performed
c f(j) = 2*(.5*f(1) + ((-1)**j)*f(n+1) + sum(f(n)*cos(pi*n*j/nx)))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxd >= nx + 1
c nxhd >= nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sctd
      dimension f(nxd), t(nxhd), mixup(nxhd), sctd(nxd)
c local data
      integer nx, nxh, nxhh, nxh2, nx2, indx1, j, k, l, jb, it, lb, ll
      integer j1, nxs, nxs2, km, km2, k1, k2
      real dnx, arg, ani, at1, at2, at3, at4, sum1
      complex t1, t2
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      nx2 = nx + 2
      indx1 = indx - 1
      if (isign.ne.0) go to 40
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/real(nx+nx)
      do 30 j = 1, nx
      arg = dnx*real(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
c create auxiliary array and bit-reverse to complex temporary
   40 sum1 = .5*(f(1) - f(nx+1))
      do 50 j = 2, nxh
      j1 = mixup(j)
      at2 = aimag(sctd(2*j1-1))
      at1 = (.5 + at2)*f(2*j1-1)
      at2 = (.5 - at2)*f(nx2-2*j1+1)
      sum1 = sum1 + real(sctd(2*j1-1))*f(2*j1-1)
      at4 = aimag(sctd(2*j1))
      at3 = (.5 + at4)*f(2*j1)
      at4 = (.5 - at4)*f(nx2-2*j1)
      sum1 = sum1 + real(sctd(2*j1))*f(2*j1)
      t(j) = cmplx(at1+at2,at3+at4)
   50 continue
      at4 = aimag(sctd(2))
      at3 = (.5 + at4)*f(2)
      at4 = (.5 - at4)*f(nx)
      sum1 = sum1 + real(sctd(2))*f(2)
      t(1) = cmplx(.5*(f(1)+f(nx+1)),at3+at4)
c reduction
      do 80 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 4*km
      do 70 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 60 j = 1, nxs
      t1 = sctd(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize result
c inverse fourier transform
      if (isign.lt.0) then
         ani = 1./real(2*nx)
         do 90 j = 2, nxhh
         t2 = conjg(t(nxh2-j))
         t1 = t(j) + t2
         t2 = (t(j) - t2)*cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         t(j) = ani*(t1 + t2)
         t(nxh2-j) = ani*conjg(t1 - t2)
   90    continue
         ani = 2.*ani
         t(nxhh+1) = ani*conjg(t(nxhh+1))
         t(1) = ani*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(
     11)))
         sum1 = ani*sum1
c forward fourier transform
      else if (isign.gt.0) then
         do 100 j = 2, nxhh
         t2 = conjg(t(nxh2-j))
         t1 = t(j) + t2
         t2 = (t(j) - t2)*cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         t(j) = t1 + t2
         t(nxh2-j) = conjg(t1 - t2)
  100    continue
         t(nxhh+1) = 2.0*conjg(t(nxhh+1))
         t(1) = 2.0*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(
     11)))
         sum1 = sum1 + sum1
      endif
c perform recursion and move to real destination
      f(1) = real(t(1))
      f(2) = sum1
      do 110 j = 2, nxh
      sum1 = sum1 - aimag(t(j))
      f(2*j-1) = real(t(j))
      f(2*j) = sum1
  110 continue
      f(nx+1) = aimag(t(1))
      return
      end
      subroutine FST1R2(f,t,isign,mixup,sctd,indx,nxd,nxhd)
c this subroutine performs two one dimensional fast real sine transforms
c and their inverses, using complex arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = 0, input: all except f, output: mixup, sctd
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 23), where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, inverse sine transforms are performed
c f(1:2,n) = (1/nx)*sum(f(1:2,j)*sin(pi*n*j/nx))
c if isign = 1, forward sine transforms are performed
c f(1:2,j) = 2*sum(f(1:2,n)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxd >= nx + 1
c nxhd >= nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sctd
      dimension f(2,nxd), t(2,nxhd), mixup(nxhd), sctd(nxd)
c local data
      integer nx, nxh, nxhh, nxh2, nx2, indx1, j, k, l, jb, it, lb, ll
      integer j1, nxs, nxs2, km, km2, k1, k2, jj
      real dnx, arg, ani, at1, at2, at3, at4, at5, at6, sum1, sum2
      complex t1, t2, t3
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      nx2 = nx + 2
      indx1 = indx - 1
      if (isign.ne.0) go to 40
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/real(nx+nx)
      do 30 j = 1, nx
      arg = dnx*real(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
c create auxiliary array and bit-reverse to complex temporary
   40 do 50 j = 2, nxh
      j1 = mixup(j)
      at5 = -aimag(sctd(2*j1-1))
      at1 = (at5 + .5)*f(1,2*j1-1)
      at2 = (at5 - .5)*f(1,nx2-2*j1+1)
      at6 = -aimag(sctd(2*j1))
      at3 = (at6 + .5)*f(1,2*j1)
      at4 = (at6 - .5)*f(1,nx2-2*j1)
      t(1,j) = cmplx(at1+at2,at3+at4)
      at1 = (at5 + .5)*f(2,2*j1-1)
      at2 = (at5 - .5)*f(2,nx2-2*j1+1)
      at3 = (at6 + .5)*f(2,2*j1)
      at4 = (at6 - .5)*f(2,nx2-2*j1)
      t(2,j) = cmplx(at1+at2,at3+at4)
   50 continue
      at6 = -aimag(sctd(2))
      at3 = (at6 + .5)*f(1,2)
      at4 = (at6 - .5)*f(1,nx2-2)
      t(1,1) = cmplx(0.,at3+at4)
      at3 = (at6 + .5)*f(2,2)
      at4 = (at6 - .5)*f(2,nx2-2)
      t(2,1) = cmplx(0.,at3+at4)
c reduction
      do 80 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 4*km
      do 70 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 60 j = 1, nxs
      t1 = sctd(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize result
c inverse fourier transform
      if (isign.lt.0) then
         ani = 1./real(2*nx)
         do 100 j = 2, nxhh
         t3 = cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         do 90 jj = 1, 2
         t2 = conjg(t(jj,nxh2-j))
         t1 = t(jj,j) + t2
         t2 = (t(jj,j) - t2)*t3
         t(jj,j) = ani*(t1 + t2)
         t(jj,nxh2-j) = ani*conjg(t1 - t2)
   90    continue
  100    continue
         ani = 2.*ani
         do 110 jj = 1, 2
         t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
         t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)
     1) - aimag(t(jj,1)))
  110    continue
c forward fourier transform
      else if (isign.gt.0) then
         do 130 j = 2, nxhh
         t3 = cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         do 120 jj = 1, 2
         t2 = conjg(t(jj,nxh2-j))
         t1 = t(jj,j) + t2
         t2 = (t(jj,j) - t2)*t3
         t(jj,j) = t1 + t2
         t(jj,nxh2-j) = conjg(t1 - t2)
  120    continue
  130    continue
         do 140 jj = 1, 2
         t(jj,nxhh+1) = 2.0*conjg(t(jj,nxhh+1))
         t(jj,1) = 2.0*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)
     1) - aimag(t(jj,1)))
  140    continue
      endif
c perform recursion and move to real destination
      sum1 = .5*real(t(1,1))
      sum2 = .5*real(t(2,1))
      f(1,1) = 0.0
      f(2,1) = 0.0
      f(1,2) = sum1
      f(2,2) = sum2
      do 150 j = 2, nxh
      sum1 = sum1 + real(t(1,j))
      sum2 = sum2 + real(t(2,j))
      f(1,2*j-1) = -aimag(t(1,j))
      f(1,2*j) = sum1
      f(2,2*j-1) = -aimag(t(2,j))
      f(2,2*j) = sum2
  150 continue
      f(1,nx+1) = 0.0
      f(2,nx+1) = 0.0
      return
      end
      subroutine FCT1R2(f,t,isign,mixup,sctd,indx,nxd,nxhd)
c this subroutine performs two one dimensional fast real cosine
c transforms and their inverses, using complex arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = 0, input: all except f, output: mixup, sctd
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 23), where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, inverse cosine transforms are performed
c f(1:2,n) = (1/nx)*(.5*f(1:2,1) + ((-1)**n)*f(1:2,nx+1)
c            + sum(f(1:2,j)*cos(pi*n*j/nx)))
c if isign = 1, forward cosine transforms are performed
c f(1:2,j) = 2*(.5*f(1:2,1) + ((-1)**j)*f(1:2,n+1)
c            + sum(f(1:2,n)*cos(pi*n*j/nx)))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxd >= nx + 1
c nxhd >= nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sctd
      dimension f(2,nxd), t(2,nxhd), mixup(nxhd), sctd(nxd)
c local data
      integer nx, nxh, nxhh, nxh2, nx2, indx1, j, k, l, jb, it, lb, ll
      integer j1, nxs, nxs2, km, km2, k1, k2, jj
      real dnx, arg, ani, at1, at2, at3, at4, at5, at6, at7, at8
      real sum1, sum2
      complex t1, t2, t3
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      nx2 = nx + 2
      indx1 = indx - 1
      if (isign.ne.0) go to 40
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/real(nx+nx)
      do 30 j = 1, nx
      arg = dnx*real(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
c create auxiliary array and bit-reverse to complex temporary
   40 sum1 = .5*(f(1,1) - f(1,nx+1))
      sum2 = .5*(f(2,1) - f(2,nx+1))
      do 50 j = 2, nxh
      j1 = mixup(j)
      at7 = real(sctd(2*j1-1))
      at5 = aimag(sctd(2*j1-1))
      at1 = (.5 + at5)*f(1,2*j1-1)
      at2 = (.5 - at5)*f(1,nx2-2*j1+1)
      sum1 = sum1 + at7*f(1,2*j1-1)
      at8 = real(sctd(2*j1))
      at6 = aimag(sctd(2*j1))
      at3 = (.5 + at6)*f(1,2*j1)
      at4 = (.5 - at6)*f(1,nx2-2*j1)
      sum1 = sum1 + at8*f(1,2*j1)
      t(1,j) = cmplx(at1+at2,at3+at4)
      at1 = (.5 + at5)*f(2,2*j1-1)
      at2 = (.5 - at5)*f(2,nx2-2*j1+1)
      sum2 = sum2 + at7*f(2,2*j1-1)
      at3 = (.5 + at6)*f(2,2*j1)
      at4 = (.5 - at6)*f(2,nx2-2*j1)
      sum2 = sum2 + at8*f(2,2*j1)
      t(2,j) = cmplx(at1+at2,at3+at4)
   50 continue
      at8 = real(sctd(2))
      at6 = aimag(sctd(2))
      at3 = (.5 + at6)*f(1,2)
      at4 = (.5 - at6)*f(1,nx)
      sum1 = sum1 + at8*f(1,2)
      t(1,1) = cmplx(.5*(f(1,1)+f(1,nx+1)),at3+at4)
      at3 = (.5 + at6)*f(2,2)
      at4 = (.5 - at6)*f(2,nx)
      sum2 = sum2 + at8*f(2,2)
      t(2,1) = cmplx(.5*(f(2,1)+f(2,nx+1)),at3+at4)
c reduction
      do 80 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 4*km
      do 70 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 60 j = 1, nxs
      t1 = sctd(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize result
c inverse fourier transform
      if (isign.lt.0) then
         ani = 1./real(2*nx)
         do 100 j = 2, nxhh
         t3 = cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         do 90 jj = 1, 2
         t2 = conjg(t(jj,nxh2-j))
         t1 = t(jj,j) + t2
         t2 = (t(jj,j) - t2)*t3
         t(jj,j) = ani*(t1 + t2)
         t(jj,nxh2-j) = ani*conjg(t1 - t2)
   90    continue
  100    continue
         ani = 2.*ani
         do 110 jj = 1, 2
         t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
         t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)
     1) - aimag(t(jj,1)))
  110    continue
         sum1 = ani*sum1
         sum2 = ani*sum2
c forward fourier transform
      else if (isign.gt.0) then
         do 130 j = 2, nxhh
         t3 = cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         do 120 jj = 1, 2
         t2 = conjg(t(jj,nxh2-j))
         t1 = t(jj,j) + t2
         t2 = (t(jj,j) - t2)*t3
         t(jj,j) = t1 + t2
         t(jj,nxh2-j) = conjg(t1 - t2)
  120    continue
  130    continue
         do 140 jj = 1, 2
         t(jj,nxhh+1) = 2.0*conjg(t(jj,nxhh+1))
         t(jj,1) = 2.0*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)
     1) - aimag(t(jj,1)))
  140    continue
         sum1 = sum1 + sum1
         sum2 = sum2 + sum2
      endif
c perform recursion and move to real destination
      f(1,1) = real(t(1,1))
      f(2,1) = real(t(2,1))
      f(1,2) = sum1
      f(2,2) = sum2
      do 150 j = 2, nxh
      sum1 = sum1 - aimag(t(1,j))
      sum2 = sum2 - aimag(t(2,j))
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = sum1
      f(2,2*j-1) = real(t(2,j))
      f(2,2*j) = sum2
  150 continue
      f(1,nx+1) = aimag(t(1,1))
      f(2,nx+1) = aimag(t(2,1))
      return
      end
      subroutine FCST1R3(f,t,isign,mixup,sctd,indx,nxd,nxhd)
c this subroutine performs three one dimensional fast real cosine and 
c sine transforms and their inverses, using complex arithmetic
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions.
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = 0, input: all except f, output: mixup, sctd
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 23), where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n) = (1/nx)*(.5*f(1,1) + ((-1)**n)*f(1,nx+1)
c          + sum(f(1,j)*cos(pi*n*j/nx)))
c f(2:3,n) = (1/nx)*sum(f(2:3,j)*sin(pi*n*j/nx))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j) = 2*(.5*f(1,1) + ((-1)**j)*f(1,n+1)
c          + sum(f(1,n)*cos(pi*n*j/nx)))
c f(2:3,j) = 2*sum(f(2:3,n)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxd >= nx + 1
c nxhd >= nx/2
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd
      real f
      complex t, sctd
      dimension f(3,nxd), t(3,nxhd), mixup(nxhd), sctd(nxd)
c local data
      integer nx, nxh, nxhh, nxh2, nx2, indx1, j, k, l, jb, it, lb, ll
      integer j1, nxs, nxs2, km, km2, k1, k2, jj
      real dnx, arg, ani, at1, at2, at3, at4, at5, at6, at7, at8
      real sum1, sum2, sum3
      complex t1, t2, t3, t4
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      nx2 = nx + 2
      indx1 = indx - 1
      if (isign.ne.0) go to 40
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxh
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/real(nx+nx)
      do 30 j = 1, nx
      arg = dnx*real(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
c create auxiliary array and bit-reverse to complex temporary
   40 sum1 = .5*(f(1,1) - f(1,nx+1))
      do 50 j = 2, nxh
      j1 = mixup(j)
      at7 = real(sctd(2*j1-1))
      at5 = aimag(sctd(2*j1-1))
      at1 = (.5 + at5)*f(1,2*j1-1)
      at2 = (.5 - at5)*f(1,nx2-2*j1+1)
      sum1 = sum1 + at7*f(1,2*j1-1)
      at8 = real(sctd(2*j1))
      at6 = aimag(sctd(2*j1))
      at3 = (.5 + at6)*f(1,2*j1)
      at4 = (.5 - at6)*f(1,nx2-2*j1)
      sum1 = sum1 + at8*f(1,2*j1)
      t(1,j) = cmplx(at1+at2,at3+at4)
      at1 = (.5 - at5)*f(2,2*j1-1)
      at2 = -(.5 + at5)*f(2,nx2-2*j1+1)
      at3 = (.5 - at6)*f(2,2*j1)
      at4 = -(.5 + at6)*f(2,nx2-2*j1)
      t(2,j) = cmplx(at1+at2,at3+at4)
      at1 = (.5 - at5)*f(3,2*j1-1)
      at2 = -(.5 + at5)*f(3,nx2-2*j1+1)
      at3 = (.5 - at6)*f(3,2*j1)
      at4 = -(.5 + at6)*f(3,nx2-2*j1)
      t(3,j) = cmplx(at1+at2,at3+at4)
   50 continue
      at8 = real(sctd(2))
      at6 = aimag(sctd(2))
      at3 = (.5 + at6)*f(1,2)
      at4 = (.5 - at6)*f(1,nx)
      sum1 = sum1 + at8*f(1,2)
      t(1,1) = cmplx(.5*(f(1,1)+f(1,nx+1)),at3+at4)
      at3 = (.5 - at6)*f(2,2)
      at4 = -(.5 + at6 )*f(2,nx2-2)
      t(2,1) = cmplx(0.,at3+at4)
      at3 = (.5 - at6)*f(3,2)
      at4 = -(.5 + at6)*f(3,nx2-2)
      t(3,1) = cmplx(0.,at3+at4)
c reduction
      do 80 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 4*km
      do 70 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 60 j = 1, nxs
      t1 = sctd(1+km2*(j-1))
      t2 = t1*t(1,j+k2)
      t3 = t1*t(2,j+k2)
      t4 = t1*t(3,j+k2)
      t(1,j+k2) = t(1,j+k1) - t2
      t(2,j+k2) = t(2,j+k1) - t3
      t(3,j+k2) = t(3,j+k1) - t4
      t(1,j+k1) = t(1,j+k1) + t2
      t(2,j+k1) = t(2,j+k1) + t3
      t(3,j+k1) = t(3,j+k1) + t4
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize result
c inverse fourier transform
      if (isign.lt.0) then
         ani = 1./real(2*nx)
         do 100 j = 2, nxhh
         t3 = cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         do 90 jj = 1, 3
         t2 = conjg(t(jj,nxh2-j))
         t1 = t(jj,j) + t2
         t2 = (t(jj,j) - t2)*t3
         t(jj,j) = ani*(t1 + t2)
         t(jj,nxh2-j) = ani*conjg(t1 - t2)
   90    continue
  100    continue
         ani = 2.*ani
         do 110 jj = 1, 3
         t(jj,nxhh+1) = ani*conjg(t(jj,nxhh+1))
         t(jj,1) = ani*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)
     1) - aimag(t(jj,1)))
  110    continue
         sum1 = ani*sum1
c forward fourier transform
      else if (isign.gt.0) then
         do 130 j = 2, nxhh
         t3 = cmplx(aimag(sctd(2*j-1)),-real(sctd(2*j-1)))
         do 120 jj = 1, 3
         t2 = conjg(t(jj,nxh2-j))
         t1 = t(jj,j) + t2
         t2 = (t(jj,j) - t2)*t3
         t(jj,j) = t1 + t2
         t(jj,nxh2-j) = conjg(t1 - t2)
  120    continue
  130    continue
         do 140 jj = 1, 3
         t(jj,nxhh+1) = 2.0*conjg(t(jj,nxhh+1))
         t(jj,1) = 2.0*cmplx(real(t(jj,1)) + aimag(t(jj,1)),real(t(jj,1)
     1) - aimag(t(jj,1)))
  140    continue
         sum1 = sum1 + sum1
      endif
c perform recursion and move to real destination
      sum2 = .5*real(t(2,1))
      sum3 = .5*real(t(3,1))
      f(1,1) = real(t(1,1))
      f(2,1) = 0.0
      f(3,1) = 0.0
      f(1,2) = sum1
      f(2,2) = sum2
      f(3,2) = sum3
      do 150 j = 2, nxh
      sum1 = sum1 - aimag(t(1,j))
      sum2 = sum2 + real(t(2,j))
      sum3 = sum3 + real(t(3,j))
      f(1,2*j-1) = real(t(1,j))
      f(1,2*j) = sum1
      f(2,2*j-1) = -aimag(t(2,j))
      f(2,2*j) = sum2
      f(3,2*j-1) = -aimag(t(3,j))
      f(3,2*j) = sum3
  150 continue
      f(1,nx+1) = aimag(t(1,1))
      f(2,nx+1) = 0.0
      f(3,nx+1) = 0.0
      return
      end
      subroutine FDST1RX(f,t,isign,mixup,sctd2,indx,nxd,nxhd,nx2d)
c this subroutine performs a one dimensional fast real sine transform
c and its inverse, using complex arithmetic
c algorithm is similar to the one described in Numerical Recipies in
c Fortran, Second Ed., by W. H. Press, B. P. Flannery, S. A. Teukolsky,
c and W. T. Vetterling, [Cambridge Univ. Press, 1992], p. 513.
c for isign = 0, input: all except f, output: mixup, sctd
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22), where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse sine transform DST-III is performed
c f(n) = (1/nx)*(.5*f(nx+1)*(-1)**n + sum(f(j+1)*sin(pi*(n+1/2)*j/nx)))
c if isign = 1, a forward sine transform DST-II is performed
c f(j+1) = 2.0*sum(f(n)*sin(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd2 = sine/cosine table
c nxd >= nx + 1
c nxhd >= nx/2
c nx2d >= 2*nx
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd, nx2d
      real f
      complex t, sctd2
      dimension f(nxd), t(nxhd), mixup(nxhd), sctd2(nx2d)
c local data
      integer nx, nxh, nxhh, nxh2, nx1, nxh1, indx1, j, k, l, lb, ll, jb
      integer it, j1, nxs, nxs2, km, km2, k1, k2
      real dnx, arg, ani, at1, at2, at3, at4, sum1
      complex t1, t2
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      nx1 = nx + 1
      nxh1 = nxh + 1
      indx1 = indx - 1
      if (isign) 50,10,140
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxh
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/real(nx+nx)
      arg = 0.25*6.28318530717959/real(nx)
      at3 = cos(arg)
      at4 = sin(arg)
      do 40 j = 1, nx
      arg = dnx*real(j - 1)
      at1 = cos(arg)
      at2 = sin(arg)
      sctd2(2*j-1) = cmplx(at1,-at2)
      at1 = at1*at4 + at2*at3
      sctd2(2*j) = cmplx(at1,1.0/at1)
   40 continue
      return
c inverse fourier transform
c create auxiliary array for cosine transform to complex temporary
   50 t(1) = cmplx(2.0*f(2),f(nx+1))
      do 60 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd2(2*j1-1))
      at2 = -aimag(sctd2(2*j1-1))
      at4 = f(2*j1) - f(2*j1-2)
      at3 = f(2*j1-1)*at2 + at4*at1
      at4 = at4*at2- f(2*j1-1)*at1
      t(j1) = cmplx(at3,at4)
   60 continue
c scramble coefficients
      do 70 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sctd2(4*j-3)),real(sctd2(4*j-3)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
   70 continue
      t(nxhh+1) = 2.*conjg(t(nxhh+1))
      t(1) = cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1)))
c bit-reverse array elements to real destination
      do 80 j = 1, nxh
      f(2*j-1) = real(t(mixup(j)))
      f(2*j) = aimag(t(mixup(j)))
   80 continue
c move back to complex temporary
      do 90 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
   90 continue
c reduction
      do 120 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 8*km
      do 110 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 100 j = 1, nxs
      t1 = conjg(sctd2(1+km2*(j-1)))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  100 continue
  110 continue
  120 continue
c perform recursion for cosine transform to real destination
      ani = 1./real(4*nx)
      do 130 j = 1, nxh
      at2 = aimag(t(nxh1-j))
      at1 = real(t(j)) + at2
      at2 = real(t(j)) - at2
      at1 = .5*at1*aimag(sctd2(4*j-2))
      f(2*j-1) = ani*(at1 + at2)
      at2 = real(t(nxh1-j))
      at1 = aimag(t(j)) + at2
      at2 = aimag(t(j)) - at2
      at1 = .5*at1*aimag(sctd2(4*j))
      f(2*j) = ani*(at1 + at2)
  130 continue
      f(nx+1) = 0.0
      return
c forward fourier transform
c create auxiliary array and bit-reverse to complex temporary
  140 do 150 j = 1, nxh
      j1 = mixup(j)
      at2 = real(sctd2(4*j1-2))
      at1 = (at2 + .5)*f(2*j1-1)
      at2 = (at2 - .5)*f(nx1-2*j1+1)
      at4 = real(sctd2(4*j1))
      at3 = (at4 + .5)*f(2*j1)
      at4 = (at4 - .5)*f(nx1-2*j1)
      t(j) = cmplx(at1+at2,at3+at4)
  150 continue
c reduction
      do 180 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 8*km
      do 170 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 160 j = 1, nxs
      t1 = sctd2(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  160 continue
  170 continue
  180 continue
c unscramble coefficients and normalize result
      do 190 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sctd2(4*j-3)),-real(sctd2(4*j-3)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
  190 continue
      t(nxhh+1) = 2.0*conjg(t(nxhh+1))
      t(1) = 2.0*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1))
     1)
c perform recursion and move to real destination
      f(nx+1) = aimag(t(1))
      f(2) = 0.5*real(t(1))
      f(1) = 0.0
      sum1 = f(2)
      do 200 j = 2, nxh
      at1 = real(sctd2(2*j-1))
      at2 = -aimag(sctd2(2*j-1))
      at3 = real(t(j))*at2 - aimag(t(j))*at1
      at2 = real(t(j))*at1 + aimag(t(j))*at2
      f(2*j-1) = at3
      sum1 = sum1 + at2
      f(2*j) = sum1
  200 continue
      return
      end
      subroutine FDCT1RX(f,t,isign,mixup,sctd2,indx,nxd,nxhd,nx2d)
c this subroutine performs a one dimensional fast real cosine transform
c and its inverse, using complex arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = 0, input: all except f, output: mixup, sctd
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22), where N = nx/2
c f = input and output data
c t = scratch array
c indx = power of 2 which determines length of transform, nx = 2**indx
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse cosine transform DCT-III is performed
c f(n) = (1/nx)*(.5*f(1) + sum(f(j)*cos(pi*(n+1/2)*j/nx)))
c if isign = 1, a forward cosine transform DCT-II is performed
c f(j) = 2.0*sum(f(n)*cos(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd2 = sine/cosine table
c nxd >= nx + 1
c nxhd >= nx/2
c nx2d >= 2*nx
c written by viktor k. decyk, ucla
c scalar version
      implicit none
      integer isign, mixup, indx, nxd, nxhd, nx2d
      real f
      complex t, sctd2
      dimension f(nxd), t(nxhd), mixup(nxhd), sctd2(nx2d)
c local data
      integer nx, nxh, nxhh, nxh2, nx1, nxh1, indx1, j, k, l, lb, ll, jb
      integer it, j1, nxs, nxs2, km, km2, k1, k2
      real dnx, arg, ani, at1, at2, at3, at4, sum1
      complex t1, t2
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      nx1 = nx + 1
      nxh1 = nxh + 1
      indx1 = indx - 1
      arg = 0.25*6.28318530717959/real(nx)
      if (isign) 50,10,140
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxh
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/real(nx+nx)
      arg = 0.25*6.28318530717959/real(nx)
      at3 = cos(arg)
      at4 = sin(arg)
      do 40 j = 1, nx
      arg = dnx*real(j - 1)
      at1 = cos(arg)
      at2 = sin(arg)
      sctd2(2*j-1) = cmplx(at1,-at2)
      at1 = at1*at4 + at2*at3
      sctd2(2*j) = cmplx(at1,1.0/at1)
   40 continue
      return
c inverse fourier transform
c create auxiliary array for cosine transform to complex temporary
   50 sum1 = f(nx)
      do 60 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd2(2*j1-1))
      at2 = -aimag(sctd2(2*j1-1))
      at4 = f(2*j1) - f(2*j1-2)
      at3 = f(2*j1-1)*at1 + at4*at2
      at4 = f(2*j1-1)*at2 - at4*at1
      t(j1) = cmplx(at3,at4)
   60 continue
      t(1) = cmplx(f(1),-2.0*sum1)
c scramble coefficients
      do 70 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sctd2(4*j-3)),real(sctd2(4*j-3)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
   70 continue
      t(nxhh+1) = 2.*conjg(t(nxhh+1))
      t(1) = cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1)))
c bit-reverse array elements to real destination
      do 80 j = 1, nxh
      f(2*j-1) = real(t(mixup(j)))
      f(2*j) = aimag(t(mixup(j)))
   80 continue
c move back to complex temporary
      do 90 j = 1, nxh
      t(j) = cmplx(f(2*j-1),f(2*j))
   90 continue
c reduction
      do 120 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 8*km
      do 110 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 100 j = 1, nxs
      t1 = conjg(sctd2(1+km2*(j-1)))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  100 continue
  110 continue
  120 continue
c perform recursion for cosine transform to real destination
      ani = 1./real(4*nx)
      do 130 j = 1, nxh
      at2 = aimag(t(nxh1-j))
      at1 = real(t(j)) + at2
      at2 = real(t(j)) - at2
      at2 = .5*at2*aimag(sctd2(4*j-2))
      f(2*j-1) = ani*(at1 - at2)
      at2 = real(t(nxh1-j))
      at1 = aimag(t(j)) + at2
      at2 = aimag(t(j)) - at2
      at2 = .5*at2*aimag(sctd2(4*j))
      f(2*j) = ani*(at1 - at2)
  130 continue
      f(nx+1) = 0.0
      return
c forward fourier transform
c create auxiliary array and bit-reverse to complex temporary
  140 do 150 j = 1, nxh
      j1 = mixup(j)
      at2 = -real(sctd2(4*j1-2))
      at1 = (.5 + at2)*f(2*j1-1)
      at2 = (.5 - at2)*f(nx1-2*j1+1)
      at4 = -real(sctd2(4*j1))
      at3 = (.5 + at4)*f(2*j1)
      at4 = (.5 - at4)*f(nx1-2*j1)
      t(j) = cmplx(at1+at2,at3+at4)
  150 continue
c reduction
      do 180 l = 1, indx1
      nxs = 2**(l - 1)
      nxs2 = nxs + nxs
      km = nxhh/nxs
      km2 = 8*km
      do 170 k = 1, km
      k1 = nxs2*(k - 1)
      k2 = k1 + nxs
      do 160 j = 1, nxs
      t1 = sctd2(1+km2*(j-1))*t(j+k2)
      t(j+k2) = t(j+k1) - t1
      t(j+k1) = t(j+k1) + t1
  160 continue
  170 continue
  180 continue
c unscramble coefficients and normalize result
      do 190 j = 2, nxhh
      t2 = conjg(t(nxh2-j))
      t1 = t(j) + t2
      t2 = (t(j) - t2)*cmplx(aimag(sctd2(4*j-3)),-real(sctd2(4*j-3)))
      t(j) = t1 + t2
      t(nxh2-j) = conjg(t1 - t2)
  190 continue
      t(nxhh+1) = 2.0*conjg(t(nxhh+1))
      t(1) = 2.0*cmplx(real(t(1)) + aimag(t(1)),real(t(1)) - aimag(t(1))
     1)
c perform recursion and move to real destination
      sum1 = -0.5*aimag(t(1))
      do 200 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd2(2*j1-1))
      at2 = -aimag(sctd2(2*j1-1))
      at3 = real(t(j1))*at1 + aimag(t(j1))*at2
      at2 = -real(t(j1))*at2 + aimag(t(j1))*at1
      f(2*j1-1) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(2*j1) = at1
  200 continue
      f(1) = real(t(1))
      f(2) = sum1
      f(nx+1) = 0.0
      return
      end
