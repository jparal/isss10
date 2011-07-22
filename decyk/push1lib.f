c-----------------------------------------------------------------------
c 1d PIC library for pushing particles and depositing charge
c push1lib.f contains procedures to process particles:
c GPOST1 deposits charge density, quadratic interpolation, STANDARD
c        optimization.
c GSPOST1 deposits charge density, quadratic interpolation, LOOKAHEAD
c         optimization
c GSPOST1X deposits charge density, quadratic interpolation, VECTOR
c          optimization.
c GPOST1L deposits charge density, linear interpolation, STANDARD
c         optimization.
c GSPOST1L deposits charge density, linear interpolation, LOOKAHEAD
c          optimization.
c GSPOST1XL deposits charge density, linear interpolation, VECTOR
c           optimization.
c GPUSH1 push particles, quadratic interpolation, STANDARD optimization.
c GSPUSH1 push particles, quadratic interpolation, LOOKAHEAD
c         optimization.
c GPUSH1L push particles, linear interpolation, STANDARD optimization.
c GSPUSH1L push particles, linear interpolation, LOOKAHEAD optimization.
c SORTP1X sort particles by grid, quadratic interpolation, memory
c         conserving algorithm.
c SORTP1XL sort particles by grid, linear interpolation, memory
c          conserving algorithm.
c DSORTP1X sort particles by grid, quadratic interpolation, high
c          performance algorithm.
c DSORTP1XL sort particles by grid, linear interpolation, high
c           performance algorithm.
c RMOVE1 remove particles instead of reflecting at boundary.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: june 29, 2010
c-----------------------------------------------------------------------
      subroutine GPOST1(part,q,qm,nop,idimp,nxv)
c for 1d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c scalar version using guard cells
c 16 flops/particle, 4 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(.75-dx**2), q(n+1)=.5*qm*(.5+dx)**2, q(n-1)=.5*qm*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c q(j) = charge density at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 2
c nxv = first dimension of charge array, must be >= nx+3
      implicit none
      integer nop, idimp, nxv
      real part, q, qm
      dimension part(idimp,nop), q(nxv)
c local data
      integer j, nn
      real qmh, dx
      qmh = .5*qm
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) + .5
      dx = part(1,j) - real(nn)
      nn = nn + 1
c deposit charge
      q(nn) = q(nn) + qmh*(.5 - dx)**2
      q(nn+1) = q(nn+1) + qm*(.75 - dx*dx)
      q(nn+2) = q(nn+2) + qmh*(.5 + dx)**2
   10 continue
      return
      end
      subroutine GSPOST1(part,q,qm,nop,idimp,nxv)
c for 1d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation
c 16 flops/particle, 4 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(.75-dx**2), q(n+1)=.5*qm*(.5+dx)**2, q(n-1)=.5*qm*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c q(j) = charge density at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 2
c nxv = first dimension of charge array, must be >= nx+3
      implicit none
      integer nop, idimp, nxv
      real part, q, qm
      dimension part(idimp,nop), q(nxv)
c local data
      integer nnn, j, nn
      real dxn, qmh, dx
c begin first particle
      nnn = part(1,1) + .5
      dxn = part(1,1) - real(nnn)
      qmh = .5*qm
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      nnn = part(1,j) + .5
      dx = dxn
      dxn = part(1,j) - real(nnn)
c deposit charge
      q(nn) = q(nn) + qmh*(.5 - dx)**2
      q(nn+1) = q(nn+1) + qm*(.75 - dx*dx)
      q(nn+2) = q(nn+2) + qmh*(.5 + dx)**2
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      q(nn) = q(nn) + qmh*(.5 - dxn)**2
      q(nn+1) = q(nn+1) + qm*(.75 - dxn*dxn)
      q(nn+2) = q(nn+2) + qmh*(.5 + dxn)**2
      return
      end
      subroutine GSPOST1X(part,q,qm,nop,idimp,nxv)
c for 1d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with guard cells
c 16 flops/particle, 4 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(.75-dx**2), q(n+1)=.5*qm*(.5+dx)**2, q(n-1)=.5*qm*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c q(j) = charge density at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 2
c nxv = first dimension of charge array, must be >= nx+3
      implicit none
      integer nop, idimp, nxv
      real part, q,  qm
      dimension part(idimp,nop), q(nxv)
c local data
      integer npp, nn, npb, ipp, k, j, je, jb, n, i
      real amx, qmh, dx
      parameter(npp=512)
      dimension nn(3,npp), amx(3,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = real(nop - 1)/real(npb) + 1.
      qmh = .5*qm
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb) + .5
      dx = part(1,j+jb) - real(n)
      n = n + 1
      amx(1,j) = qmh*(.5 - dx)**2
      nn(1,j) = n
      amx(2,j) = qm*(.75 - dx*dx)
      nn(2,j) = n + 1
      amx(3,j) = qmh*(.5 + dx)**2
      nn(3,j) = n + 2
   10 continue
c deposit charge
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 3
      q(nn(i,j)) = q(nn(i,j)) + amx(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GPOST1L(part,q,qm,nop,idimp,nxv)
c for 1d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells
c 7 flops/particle, 3 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(1.-dx) and q(n+1)=qm*dx
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c q(j) = charge density at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 2
c nxv = first dimension of charge array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, q, qm
      dimension part(idimp,nop), q(nxv)
c local data
      integer j, nn
      real dx
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) 
      dx = qm*(part(1,j) - real(nn))
      nn = nn + 1
c deposit charge
      q(nn) = q(nn) + (qm - dx)
      q(nn+1) = q(nn+1) + dx
   10 continue
      return
      end
      subroutine GSPOST1L(part,q,qm,nop,idimp,nxv)
c for 1d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation
c 7 flops/particle, 3 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(1.-dx) and q(n+1)=qm*dx
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c q(j) = charge density at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 2
c nxv = first dimension of charge array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, q, qm
      dimension part(idimp,nop), q(nxv)
c local data
      integer nnn, j, nn
      real dxn, dx
c begin first particle
      nnn = part(1,1)
      dxn = qm*(part(1,1) - real(nnn))
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      nnn = part(1,j)
      dx = dxn
      dxn = qm*(part(1,j) - real(nnn))
c deposit charge
      q(nn) = q(nn) + (qm - dx)
      q(nn+1) = q(nn+1) + dx
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      q(nn) = q(nn) + (qm - dxn)
      q(nn+1) = q(nn+1) + dxn
      return
      end
      subroutine GSPOST1XL(part,q,qm,nop,idimp,nxv)
c for 1d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with guard cells
c 16 flops/particle, 4 loads, 3 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n)=qm*(1.-dx) and q(n+1)=qm*dx
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c q(j) = charge density at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 2
c nxv = first dimension of charge array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, q,  qm
      dimension part(idimp,nop), q(nxv)
c local data
      integer npp, nn, npb, ipp, k, j, je, jb, n, i
      real amx, dx
      parameter(npp=512)
      dimension nn(2,npp), amx(2,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = real(nop - 1)/real(npb) + 1.
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb)
      dx = qm*(part(1,j+jb) - real(n))
      n = n + 1
      amx(1,j) = qm - dx
      nn(1,j) = n
      amx(2,j) = dx
      nn(2,j) = n + 1
   10 continue
c deposit charge
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 2
      q(nn(i,j)) = q(nn(i,j)) + amx(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
c for 1d code, this subroutine updates particle co-ordinate and velocity
c using leap-frog scheme in time and second-order spline interpolation
c in space, with periodic boundary conditions.
c scalar version using guard cells
c 25 flops/particle, 5 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
c and x(t+dt) = x(t) + v(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (.75-dx**2)*fx(n)+.5*(fx(n+1)*(.5+dx)**2+fx(n-1)*(.5-dx)**2)
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fx(j) = force/charge at grid point j, that is convolution of electric
c field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c nxv = dimension of field array, must be >= nx+3
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fx, qbm, dt, ek
      dimension part(idimp,nop), fx(nxv)
c local data
      integer j, nn
      real qtm, edgelx, edgerx, dx
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = float(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = float(nx-1)
      endif
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      dx = part(1,j) - real(nn)
      nn = nn + 1
c find acceleration
      dx = (.75 - dx*dx)*fx(nn+1) + .5*(fx(nn)*(.5 - dx)**2 + fx(nn+2)*(
     1.5 + dx)**2)
c new velocity
      dx = part(2,j) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (part(2,j) + dx)**2
      part(2,j) = dx
c new position
      dx = part(1,j) + dx*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(2,j) = -part(2,j)
         endif
      endif
c set new position
      part(1,j) = dx
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GSPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
c for 1d code, this subroutine updates particle co-ordinate and velocity
c using leap-frog scheme in time and second-order spline interpolation
c in space, with periodic boundary conditions.
c scalar version using guard cells, integer conversion precalculation
c 25 flops/particle, 5 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
c and x(t+dt) = x(t) + v(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (.75-dx**2)*fx(n)+.5*(fx(n+1)*(.5+dx)**2+fx(n-1)*(.5-dx)**2)
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fx(j) = force/charge at grid point j, that is convolution of electric
c field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c nxv = dimension of field array, must be >= nx+3
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fx, qbm, dt, ek
      dimension part(idimp,nop), fx(nxv)
c local data
      integer nnn, nop1, j, nn
      real dxn, qtm, edgelx, edgerx, dx
      double precision sum1
c begin first particle
      nnn = part(1,1) + .5
      dxn = part(1,1) - real(nnn)
      nop1 = nop - 1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = float(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = float(nx-1)
      endif
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      nnn = part(1,j+1) + .5
      dx = dxn
      dxn = part(1,j+1) - real(nnn)
c find acceleration
      dx = (.75 - dx*dx)*fx(nn+1) + .5*(fx(nn)*(.5 - dx)**2 + fx(nn+2)*(
     1.5 + dx)**2)
c new velocity
      dx = part(2,j) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (part(2,j) + dx)**2
      part(2,j) = dx
c new position
      dx = part(1,j) + dx*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(2,j) = -part(2,j)
         endif
      endif
      part(1,j) = dx
   10 continue
c push last particle
      nn = nnn + 1
c find acceleration
      dx = (.75 - dxn*dxn)*fx(nn+1) + .5*(fx(nn)*(.5 - dxn)**2 + fx(nn+2
     1)*(.5 + dxn)**2)
c new velocity
      dx = part(2,nop) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (part(2,nop) + dx)**2
      part(2,nop) = dx
c new position
      dx = part(1,nop) + dx*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(2,nop) = -part(2,nop)
         endif
      endif
      part(1,nop) = dx
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
c for 1d code, this subroutine updates particle co-ordinate and velocity
c using leap-frog scheme in time and first-order linear interpolation
c in space, with periodic boundary conditions.
c scalar version using guard cells
c 16 flops/particle, 4 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
c and x(t+dt) = x(t) + v(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fx(j) = force/charge at grid point j, that is convolution of electric
c field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c nxv = first dimension of charge array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fx, qbm, dt, ek
      dimension part(idimp,nop), fx(nxv)
c local data
      integer j, nn
      real qtm, edgelx, edgerx, dx
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = float(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = float(nx-1)
      endif
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      dx = part(1,j) - real(nn)
      nn = nn + 1
c find acceleration
      dx = (1. - dx)*fx(nn) + dx*fx(nn+1)
c new velocity
      dx = part(2,j) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (part(2,j) + dx)**2
      part(2,j) = dx
c new position
      dx = part(1,j) + dx*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(2,j) = -part(2,j)
         endif
      endif
      part(1,j) = dx
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GSPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ipbc)
c for 1d code, this subroutine updates particle co-ordinate and velocity
c using leap-frog scheme in time and first-order linear interpolation
c in space, with periodic boundary conditions.
c scalar version using guard cells, integer conversion precalculation
c 16 flops/particle, 4 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c v(t+dt/2) = v(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
c and x(t+dt) = x(t) + v(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fx(j) = force/charge at grid point j, that is convolution of electric
c field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((v(t+dt/2)+v(t-dt/2))**2)
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c nxv = first dimension of charge array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fx, qbm, dt, ek
      dimension part(idimp,nop), fx(nxv)
c local data
      integer nnn, nop1, j, nn
      real dxn, qtm, edgelx, edgerx, dx
      double precision sum1
c begin first particle
      nnn = part(1,1)
      dxn = part(1,1) - real(nnn)
      nop1 = nop - 1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = float(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = float(nx-1)
      endif
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      nnn = part(1,j+1)
      dx = dxn
      dxn = part(1,j+1) - real(nnn)
c find acceleration
      dx = (1. - dx)*fx(nn) + dx*fx(nn+1)
c new velocity
      dx = part(2,j) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (part(2,j) + dx)**2
      part(2,j) = dx
c new position
      dx = part(1,j) + dx*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(2,j) = -part(2,j)
         endif
      endif
      part(1,j) = dx
   10 continue
c push last particle
      nn = nnn + 1
c find acceleration
      dx = (1. - dxn)*fx(nn) + dxn*fx(nn+1)
c new velocity
      dx = part(2,nop) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (part(2,nop) + dx)**2
      part(2,nop) = dx
c new position
      dx = part(1,nop) + dx*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(2,nop) = -part(2,nop)
         endif
      endif
      part(1,nop) = dx
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine SORTP1X(part,pt,ip,npic,idimp,nop,nx1)
c this subroutine sorts particles by grid
c quadratic interpolation
c part = particle array
c part(1,n) = position x of particle n
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 2
c nop = number of particles
c nx1 = system length in x direction + 1
      implicit none
      integer ip, npic, idimp, nop, nx1
      real part, pt
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(nx1)
c local data
      integer k, j, n, isum, ist, i
c clear counter array
      do 10 k = 1, nx1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      n = part(1,j) + 0.5
      n = n + 1
      npic(n) = npic(n) + 1
      ip(j) = n
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nx1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      n = ip(j)
      npic(n) = npic(n) + 1
      ip(j) = npic(n)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, nop
      pt(ip(j)) = part(i,j)
   50 continue
      do 60 j = 1, nop
      part(i,j) = pt(j)
   60 continue
   70 continue
      return
      end
      subroutine SORTP1XL(part,pt,ip,npic,idimp,nop,nx1)
c this subroutine sorts particles by grid
c linear interpolation
c part = particle array
c part(1,n) = position x of particle n
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 2
c nop = number of particles
c nx1 = system length in x direction + 1
      implicit none
      integer ip, npic, idimp, nop, nx1
      real part, pt
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(nx1)
c local data
      integer k, j, n, isum, ist, i
c clear counter array
      do 10 k = 1, nx1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      n = part(1,j)
      n = n + 1
      npic(n) = npic(n) + 1
      ip(j) = n
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nx1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      n = ip(j)
      npic(n) = npic(n) + 1
      ip(j) = npic(n)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, nop
      pt(ip(j)) = part(i,j)
   50 continue
      do 60 j = 1, nop
      part(i,j) = pt(j)
   60 continue
   70 continue
      return
      end
      subroutine DSORTP1X(parta,partb,npic,idimp,nop,nx1)
c this subroutine sorts particles by grid
c quadratic interpolation
c parta/partb = input/output particle arrays
c parta(1,n) = position x of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c nx1 = system length in x direction + 1
      implicit none
      integer npic, idimp, nop, nx1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(nx1)
c local data
      integer i, j, k, n, isum, ist, ip
c clear counter array
      do 10 k = 1, nx1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      n = parta(1,j) + 0.5
      n = n + 1
      npic(n) = npic(n) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nx1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      n = parta(1,j) + 0.5
      n = n + 1
      ip = npic(n) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(n) = ip
   50 continue
      return
      end
      subroutine DSORTP1XL(parta,partb,npic,idimp,nop,nx1)
c this subroutine sorts particles by grid
c linear interpolation
c parta/partb = input/output particle arrays
c parta(1,n) = position x of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c nx1 = system length in x direction + 1
      implicit none
      integer npic, idimp, nop, nx1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(nx1)
c local data
      integer i, j, k, n, isum, ist, ip
c clear counter array
      do 10 k = 1, nx1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      n = parta(1,j)
      n = n + 1
      npic(n) = npic(n) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nx1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      n = parta(1,j)
      n = n + 1
      ip = npic(n) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(n) = ip
   50 continue
      return
      end
      subroutine RMOVE1(part,ihole,nx,idimp,nop,ntmax,ipbc)
c this subroutine removes particles which would normally be reflected
c part(1,n) = position x of particle n
c ihole = location of holes left in particle arrays
c nx = system length in x direction
c idimp = size of phase space = 4
c nop = number of particles
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,-1,-2) =
c (none,2d periodic,2d reflecting)
      implicit none
      real part
      integer ihole
      integer nx, idimp, nop, ntmax, ipbc
      dimension part(idimp,nop)
      dimension ihole(ntmax)
c local data
      integer idps
      parameter(idps=2)
      integer i, j, j1, j2, nter, np, jss
      dimension jss(idps)
      real edgelx, edgerx, dx
c set boundary values
      edgelx = 0.
      edgerx = real(nx)
      if (ipbc.eq.(-2)) then
         edgelx = 1.
         edgerx = real(nx-1)
      endif
      nter = 0
      np = nop
c buffer outgoing particles
   10 jss(1) = 0
      jss(2) = 0
      do 20 j = 1, np
      dx = part(1,j)
c periodic boundary conditions
      if (ipbc.eq.(-1)) then
         if (dx.lt.edgelx) part(1,j) = dx + edgerx
         if (dx.ge.edgerx) part(1,j) = part(1,j) - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.(-2)) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (jss(1).lt.ntmax) then
               jss(1) = jss(1) + 1
               ihole(jss(1)) = j
            else
               jss(2) = 1
               go to 30
            endif
         endif
      endif
   20 continue
   30 continue
c fill up holes in particle array with particles from bottom
      do 50 j = 1, jss(1)
      j1 = np - j + 1
      j2 = jss(1) - j + 1
      if (j1.gt.ihole(j2)) then
c move particle only if it is below current hole
         do 40 i = 1, idimp
         part(i,ihole(j2)) = part(i,j1)
   40    continue
      endif
   50 continue
      np = np - jss(1)
c check if buffer overflowed and more particles remain to be checked
      if (jss(2).gt.0) then
         nter = nter + 1
         go to 10
      endif
c information
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, ntmax=', ntmax
      endif
      nop = np
      return
      end
      subroutine DPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxh,nxvh)
c for 1d code, this subroutine calculates particle charge density
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 4*nxh flops/particle
c input: all, output: q
c charge density is calculated from the expression:
c q(n) = sum(qm*exp(-sqrt(-1)*2*n*pi*x/nx))
c part(1,n) = position x of particle n
c q(n) = charge density at fourier grid point n
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx = system length in x direction
c nxh = number of fourier modes kept
c nxvh = dimension of charge array, must be >= nxh
      implicit none
      integer nop, idimp, nx, nxh, nxvh
      real part, qm
      complex q, sctx
      dimension part(idimp,nop), q(nxvh), sctx(nxvh)
c local data
      integer i, j
      real qmn, dnx, dkx
      dnx = 6.28318530717959/real(nx)
      qmn = qm/real(nx)
c find fourier components
      do 30 i = 1, nop
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nx/2
      do 20 j = 2, nxh
      q(j) = q(j) + qmn*sctx(j-1)
   20 continue
      q(1) = cmplx(real(q(1))+qmn,aimag(q(1))+qmn*real(sctx(nxh)))
   30 continue
      return
      end
      subroutine PUSH1GL(part,fx,sctx,qbm,dt,ek,idimp,nop,nx,nxh,nxvh)
c for 1d code, this subroutine updates particle co-ordinate and
c velocity using leap-frog scheme in time using gridless spectral
c version, with periodic boundaries
c scalar version using guard cells
c 2*nxh + 7 flops/particle, 2*nxh loads, 2 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t))*dt, where q/m is charge/mass,
c and x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is calculated from the expression
c fx(x) = sum(fx(n)*exp(sqrt(-1)*2*n*pi*x/nx))
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fx(j) = x component of force/charge at grid (j)
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2)
c idimp = size of phase space = 2
c nop = number of particles
c nx = system length in x direction
c nxh = number of fourier modes kept
c nxvh = dimension of field array, must be >= nxh
      implicit none
      integer nop, idimp, nx, nxh, nxvh
      real part, qbm, dt, ek
      complex fx, sctx
      dimension part(idimp,nop), fx(nxvh), sctx(nxvh)
c local data
      integer i, j
      real zero, anx, dnx, dkx, qtm, dx
      double precision sum1, ex
      zero = 0.
      anx = real(nx)
      dnx = 6.28318530717959/anx
      qtm = qbm*dt
      sum1 = 0.0d0
      do 30 i = 1, nop
c find electric field
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
c mode numbers 0 < kx < nx/2
      do 20 j = 2, nxh
      ex = ex + (real(fx(j)*sctx(j-1)))
   20 continue
      ex = 2.0d0*ex
      ex = ex + real(fx(1))
      dx = ex + aimag(fx(1))*real(sctx(nxh))
c new velocity
      dx = part(2,i) + qtm*dx
c average kinetic energy
      sum1 = sum1 + (dx + part(2,i))**2
      part(2,i) = dx
c new position
      dx = part(1,i) + dx*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
c set new position
      part(1,i) = dx
   30 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
