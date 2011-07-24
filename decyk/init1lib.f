c-----------------------------------------------------------------------
c 1d PIC library for initialization
c init1lib.f contains procedures to initialize particle
c            co-ordinates:
c DISTR1 initializes x and vx co-ordinates for 1d code.
c DISTR1H initializes x and vx,vy,vz co-ordinates for magnetized
c         1-2/2d codes.
c FDISTR1 initializes x co-ordinate for 1d code, with general
c         distribution in space.
c VDISTR1 initializes vx co-ordinate for 1d code, with maxwellian
c         velocity distribution with drift.
c VDISTR1H initializes vx, vy, vz co-ordinates for 1-2/2d code, with
c          maxwellian velocity distribution with drift.
c GBDISTR1L calculates guiding centers for magnetized 1-2/2d codes.
c ranorm = generates gaussian random numbers.
c randum = generates uniform random numbers.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: july 16, 2011
c-----------------------------------------------------------------------
      subroutine DISTR1(part,vtx,vdx,npx,idimp,nop,nx,ipbc)
c for 1d code, this subroutine calculates initial particle co-ordinate
c and velocity, with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c vtx = thermal velocity of particles in x direction
c vdx = drift velocity of particles x direction
c npx = number of particles distributed in x direction
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, idimp, nop, nx, ipbc
      real part, vtx, vdx
      dimension part(idimp,nop)
c local data
      integer j
      real edgelx, at1, sum1
      double precision dsum1
      double precision ranorm
c set boundary values
      edgelx = 0.0
      at1 = real(nx)/real(npx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
c uniform density profile
      do 10 j = 1, npx
      part(1,j) = edgelx + at1*(real(j) - .5)
   10 continue
c maxwellian velocity distribution
      do 20 j = 1, npx
      part(2,j) = vtx*ranorm()
   20 continue
c add correct drift
      dsum1 = 0.0d0
      do 30 j = 1, npx
      dsum1 = dsum1 + part(2,j)
   30 continue
      sum1 = dsum1
      sum1 = sum1/real(npx) - vdx
      do 40 j = 1, npx
      part(2,j) = part(2,j) - sum1
   40 continue
      return
      end
      subroutine DISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,idimp,nop,nx,i
     1pbc)
c for 1-2/2d code, this subroutine calculates initial particle
c co-ordinate and velocity, with uniform density and maxwellian
c velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx = number of particles distributed in x direction
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, idimp, nop, nx, ipbc
      real part, vtx, vty, vtz, vdx, vdy,v dz
      dimension part(idimp,nop)
c local data
      integer j
      real edgelx, at1, sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
c set boundary values
      edgelx = 0.0
      at1 = real(nx)/real(npx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
c uniform density profile
      do 10 j = 1, npx
      part(1,j) = edgelx + at1*(real(j) - .5)
   10 continue
c maxwellian velocity distribution
      do 20 j = 1, npx
      part(2,j) = vtx*ranorm()
      part(3,j) = vty*ranorm()
      part(4,j) = vtz*ranorm()
   20 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 30 j = 1, npx
      dsum1 = dsum1 + part(2,j)
      dsum2 = dsum2 + part(3,j)
      dsum3 = dsum3 + part(4,j)
   30 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1.0/real(npx)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 40 j = 1, npx
      part(2,j) = part(2,j) - sum1
      part(3,j) = part(3,j) - sum2
      part(4,j) = part(4,j) - sum3
   40 continue
      return
      end
      subroutine FDISTR1(part,fnx,argx1,argx2,argx3,npx,idimp,nop,nx,ipb
     1c,ierr)
c for 1d code, this subroutine calculates initial particle co-ordinates
c with general density profile where density in x is given by
c n(x) = fnx(x,argx1,argx2,argx3,0) and integral of the density is given
c by fnx(x,argx1,argx2,argx3,1)
c part(1,n) = position x of particle n
c fnx = density and density integral function in x direction
c argx1,argx2,argx3 = arguments to fnx
c npx = initial number of particles distributed in x direction
c idimp = size of phase space = 2 or 4
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer npx, idimp, nop, nx, ipbc, ierr
      real argx1, argx2, argx3, part
      dimension part(idimp,nop)
      real fnx
      external fnx
c local data
      integer imax, i, j
      real edgelx, anx, bnx, xt, xt0, x0
      real xn, eps, big, f, fp
      ierr = 0
c eps = convergence criterion
      imax = 20
      eps = 0.0001
      big = 0.5
c set boundary value
      edgelx = 0.
      if (ipbc.eq.2) then
         edgelx = 1.0
      endif
c find normalization for function
      anx = real(nx) - edgelx
      x0 = fnx(edgelx,argx1,argx2,argx3,1)
      bnx = real(npx)/(fnx(anx,argx1,argx2,argx3,1) - x0)
      x0 = bnx*x0 - .5
c density profile in x
      xt0 = edgelx
      xt = xt0 + 0.5/(bnx*fnx(xt0,argx1,argx2,argx3,0))
      do 20 j = 1, npx
      xn = real(j) + x0
c guess next value for xt
      if (j.gt.1) xt = xt + 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
      xt = max(edgelx,min(xt,anx))
      i = 0
   10 f = bnx*fnx(xt,argx1,argx2,argx3,1) - xn
c find improved value for xt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bnx*fnx(xt,argx1,argx2,argx3,0)
            xt0 = xt
            xt = xt - f/fp
            xt = max(edgelx,min(xt,anx))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(xt - xt0)
            xt = xt0 + fp
         else
            fp = xt - xt0
c           xt0 = xt
            xt = xt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 10
         write (2,*) 'newton iteration max exceeded, xt = ', xt
         ierr = ierr + 1
      endif
c store co-ordinate
      part(1,j) = xt
      xt0 = xt
   20 continue
      return
      end
      function FLDISTR1(x,anlx,anxi,shift,intg)
c this function calculates either a density function or its integral
c for a linear density profile.  Used in initializing particle
c coordinates.  The three parameters are redundant, and one can set one
c of them arbitrarily.  A convenient choice is to set  anxi = 1/Lx,
c anlx = NH - NL, shift = (1 - NL)/(NH - NL), where NL is the density
c at the left, and NH at the right compared to the average density
c if intg = 0, n(x) = 1. + anlx*(x*anxi - shift)
c if intg = 1, n(x) = x + .5*anlx*x*(x*anxi - 2.*shift)
      implicit none
      integer intg
      real x, anlx, anxi, shift
c local data
      real FLDISTR1, f
      if (intg.eq.0) then
         f = 1.0 + anlx*(x*anxi - shift)
      else if (intg.eq.1) then
         if (anxi.eq.0.) then
            f = x
         else
            f = x + .5*anlx*x*(x*anxi - 2.*shift)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FLDISTR1 Error: f = ', f
      FLDISTR1 = f
      return
      end
      function FSDISTR1(x,ans,dkx,phase,intg)
c this function calculates either a density function or its integral
c for a sinusoidal density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ans*sin(dkx*x - phase)
c if intg = 1, n(x) = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
      implicit none
      integer intg
      real x, ans, dkx, phase
c local data
      real FSDISTR1, f
      if (intg.eq.0) then
         f = 1.0 + ans*sin(dkx*x - phase)
      else if (intg.eq.1) then
         if (dkx.eq.0.) then
            f = x - ans*sin(phase)*x
         else
            f = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FSDISTR1 Error: f = ', f
      FSDISTR1 = f
      return
      end
      function FGDISTR1(x,ang,wi,x0,intg)
c this function calculates either a density function or its integral
c for a gaussian density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ang*exp(-((x-x0)*wi)**2/2.)
c if intg = 1, n(x) = x + (ang*sqrt(pi/2)/wi)*
c                         (erf((x-x0)*wi/sqrt(2)) + erf(x0*wi/sqrt(2)))
      implicit none
      integer intg
      real x, ang, x0, wi
c local data
      real FGDISTR1, f, sqrt2i, sqtpih, aw, t, erfn
      external erfn
      data sqrt2i, sqtpih /0.7071067811865476,1.253314137397325/
      save sqrt2i, sqtpih
      aw = wi*sqrt2i
      t = (x - x0)*aw
      if (intg.eq.0) then
         if (abs(t).lt.8.) then
            f = 1.0 + ang*exp(-t**2)
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + ang)*x
         else
            f = x + (ang*sqtpih/wi)*(erfn(t) + erfn(x0*aw))
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FGDISTR1 Error: f = ', f
      FGDISTR1 = f
      return
      end
      function FHDISTR1(x,anh,wi,x0,intg)
c this function calculates either a density function or its integral
c for a hyperbolic secant squared density profile.  Used in initializing
c particle coordinates.
c if intg = 0, n(x) = 1.0 + anh*sech((x-x0)*wi)**2
c if intg = 1, n(x) = x + (anh/wi)*(tanh((x-x0)*wi) + tanh(x0*wi))
      implicit none
      integer intg
      real x, anh, x0, wi
c local data
      real FHDISTR1, f, g, t, u
      t = (x - x0)*wi
      if (intg.eq.0) then
         if (abs(t).lt.32.) then
            u = exp(-abs(t))
            f = 1.0 + anh*(2.*u/(1.0 + u*u))**2
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + anh)*x
         else
            if (abs(t).lt.32.) then
               u = exp(-abs(t))**2
               f = (1.0 - u)/(1.0 + u)
            else
               f = 1.0
            endif
            if (t.lt.0.) f = -f
            t = x0*wi
            if (abs(t).lt.32.) then
               u = exp(-abs(t))**2
               g = (1.0 - u)/(1.0 + u)
            else
               g = 1.0
            endif
            if (t.lt.0.) g = -g
            f = x + (anh/wi)*(f + g)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FHDISTR1 Error: f = ', f
      FHDISTR1 = f
      return
      end
      subroutine VDISTR1(part,vtx,vdx,idimp,nop)
c for 1d code, this subroutine calculates initial particle
c velocity with maxwellian velocity with drift
c part(2,n) = velocity vx of particle n
c vtx = thermal velocity of electrons in x direction
c vdx = drift velocity of beam electrons in x direction
c idimp = size of phase space = 2
c nop = number of particles
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer idimp, nop
      real vtx, vdx, part
      dimension part(idimp,nop)
c local data
      integer j
      real sum1, at1
      double precision ranorm
      double precision dsum1
      if (nop.eq.0) return
c maxwellian velocity distribution
      do 10 j = 1, nop
      part(2,j) = vtx*ranorm()
   10 continue
c add correct drift
      dsum1 = 0.0d0
      do 20 j = 1, nop
      dsum1 = dsum1 + part(2,j)
   20 continue
      sum1 = dsum1
      at1 = 1./real(nop)
      sum1 = at1*sum1 - vdx
      do 30 j = 1, nop
      part(2,j) = part(2,j) - sum1
   30 continue
      return
      end
      subroutine VDISTR1H(part,vtx,vty,vtz,vdx,vdy,vdz,idimp,nop)
c for 1-2/2d code, this subroutine calculates initial particle
c velocities with maxwellian velocity with drift
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c idimp = size of phase space = 4
c nop = number of particles
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer idimp, nop
      real vtx, vty, vtz, vdx, vdy, vdz, part
      dimension part(idimp,nop)
c local data
      integer j
      real sum1, sum2, sum3, at1
      double precision ranorm
      double precision dsum1, dsum2, dsum3
      if (nop.eq.0) return
c maxwellian velocity distribution
      do 10 j = 1, nop
      part(2,j) = vtx*ranorm()
      part(3,j) = vty*ranorm()
      part(4,j) = vtz*ranorm()
   10 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 20 j = 1, nop
      dsum1 = dsum1 + part(2,j)
      dsum2 = dsum2 + part(3,j)
      dsum3 = dsum3 + part(4,j)
   20 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1./real(nop)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 30 j = 1, nop
      part(2,j) = part(2,j) - sum1
      part(3,j) = part(3,j) - sum2
      part(4,j) = part(4,j) - sum3
   30 continue
      return
      end
      subroutine GBDISTR1L(part,byz,qbm,idimp,nop,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine reinterprets current particle
c position as position of guiding center, and calculates the actual
c particle position
c in converting from guiding center to actual co-ordinates,
c the following equation is used:
c       x(t) = xg(t) - (vy(t)*omz - vz(t)*omy)/om**2
c where omy = (q/m)*byz(1,xg(t)),
c and   omz = (q/m)*byz(2,xg(t)),
c and the magnetic field components byz(i,x(t)) are approximated
c by interpolation from the nearest grid points:
c byz(i,x) = (1-dx)*byz(i,n)+dx*byz(i,n+1)
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c byz(i,j) = i component of magnetic field at grid j
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxv = first dimension of field arrays, must be >= nx
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, byz, qbm
      dimension part(idimp,nop)
      dimension byz(2,nxv)
c local data
      integer j, nn, n
      real edgelx, edgerx, dxp, amx, omy, omz, at3, omyt, omzt, dx
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c calculate actual position from guiding center
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      dxp = qbm*(part(1,j) - real(nn))
      nn = nn + 1
      amx = qbm - dxp
c find magnetic field
      omy = amx*byz(1,nn) + dxp*byz(1,nn+1)
      omz = amx*byz(2,nn) + dxp*byz(2,nn+1)
      at3 = sqrt(omy*omy + omz*omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3
      omyt = omy*at3
      omzt = omz*at3
c correct position
      dx = part(1,j) - (part(3,j)*omzt - part(4,j)*omyt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         n = abs(dx)/edgerx
         if (dx.lt.edgelx) dx = dx + real(n + 1)*edgerx
         if (dx.ge.edgerx) dx = dx - real(n)*edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
c if x co-ordinate only is out of bounds, try switching vy
            dx = part(1,j) + (part(3,j)*omzt + part(4,j)*omyt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
               part(3,j) = -part(3,j)
            else
c otherwise, try switching both vy and vz
               dx = part(1,j) + (part(4,j)*omzt - part(5,j)*omyt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
                  part(3,j) = -part(3,j)
                  part(4,j) = -part(4,j)
c give up if larmor radius is too large
               else
                 dx = part(1,j)
               endif
            endif
         endif
      endif
c set new position
      part(1,j) = dx
   10 continue
      return
      end
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      integer r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
      function randum()
c this is a version of the random number generator dprandom due to
c c. bingham and the yale computer center, producing numbers
c in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      integer r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
      function erfn(x)
c this function calculates the real error function, according to the
c formulae given in Abramowitz and Stegun, Handbook of Mathematical
c Functions, p. 299.  Error is < 1.5 x 10-7.
      implicit none
      real x
c local data
      real erfn, p, a1, a2, a3, a4, a5, t, f
      data p, a1, a2 /0.3275911,0.254829592,-0.284496736/
      data a3, a4, a5 /1.421413741,-1.453152027,1.061405429/
      save p, a1, a2, a3, a4, a5
      f = abs(x)
      t = 1.0/(1.0 + p*f)
      if (f.le.8.) then
         erfn = 1.0 - t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))*exp(-x*x)
      else
         erfn = 1.0
      endif
      if (x.lt.0.) erfn = -erfn
      return
      end
