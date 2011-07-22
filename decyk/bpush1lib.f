c-----------------------------------------------------------------------
c 1d PIC library for pushing particles with magnetic field
c and depositing current
c bpush1lib.f contains procedures to process particles with magnetic
c             fields:
c GJPOST1 deposits current density, quadratic interpolation, STANDARD
c         optimization.
c GSJPOST1 deposits current density, quadratic interpolation, LOOKAHEAD
c          optimization.
c GSJPOST1X deposits current density, quadratic interpolation, VECTOR
c           optimization.
c GJPOST1L deposits current density, linear interpolation, STANDARD
c          optimization.
c GSJPOST1L deposits current density, linear interpolation, LOOKAHEAD
c           optimization.
c GSJPOST1XL deposits current density, linear interpolation, VECTOR
c            optimization.
c GBPUSH13 push particles with magnetic field, quadratic interpolation,
c          STANDARD optimization.
c GSBPUSH13 push particles with magnetic field, quadratic interpolation,
c           LOOKAHEAD optimization.
c GPUSH1L push particles with magnetic field, linear interpolation,
c         STANDARD optimization.
c GSBPUSH13L push particles with magnetic field, linear interpolation,
c            LOOKAHEAD optimization.
c RETARD1 retard particle position a half time-step.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: july 3, 2010
c-----------------------------------------------------------------------
      subroutine GJPOST1(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using second-order spline interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 27 flops/particle, 10 loads, 7 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(.75-dx**2)
c cu(i,n+1)=.5*qci*(.5+dx)**2
c cu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,j+1) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+3
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2,nxv)
c local data
      integer j, nn
      real qmh, edgelx, edgerx, dxp, amx, dxl, vy, vz, dx
      qmh = .5*qm
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) + .5
      dxp = part(1,j) - real(nn)
      nn = nn + 1
      amx = qm*(.75 - dxp*dxp)
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
c deposit current
      vy = part(3,j)
      vz = part(4,j)
      cu(1,nn) = cu(1,nn) + vy*dxl
      cu(2,nn) = cu(2,nn) + vz*dxl
      cu(1,nn+1) = cu(1,nn+1) + vy*amx
      cu(2,nn+1) = cu(2,nn+1) + vz*amx
      cu(1,nn+2) = cu(1,nn+2) + vy*dxp
      cu(2,nn+2) = cu(2,nn+2) + vz*dxp
c advance position half a time-step
      dx = part(1,j) + part(2,j)*dt
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
      return
      end
      subroutine GSJPOST1(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using second-order spline interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation
c 27 flops/particle, 10 loads, 7 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(.75-dx**2)
c cu(i,n+1)=.5*qci*(.5+dx)**2
c cu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,j+1) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+3
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2,nxv)
c local data
      integer nnn, j, nn
      real dxn, qmh, edgelx, edgerx, dxp, amx, dxl, vy, vz, dx
      real dx1, dx2, dx3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      dxn = part(1,1) - real(nnn)
      qmh = .5*qm
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      nnn = part(1,j) + .5
      dxp = dxn
      dxn = part(1,j) - real(nnn)
      amx = qm*(.75 - dxp*dxp)
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
c deposit current
      vy = part(3,j-1)
      vz = part(4,j-1)
      dx1 = cu(1,nn) + vy*dxl
      dxl = cu(2,nn) + vz*dxl
      dx2 = cu(1,nn+1) + vy*amx
      amx = cu(2,nn+1) + vz*amx
      dx3 = cu(1,nn+2) + vy*dxp
      dxp = cu(2,nn+2) + vz*dxp
      cu(1,nn) = dx1
      cu(2,nn) = dxl
      cu(1,nn+1) = dx2
      cu(2,nn+1) = amx
      cu(1,nn+2) = dx3
      cu(2,nn+2) = dxp
c advance position half a time-step
      dx = part(1,j-1) + part(2,j-1)*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(2,j-1) = -part(2,j-1)
         endif
      endif
c set new position
      part(1,j-1) = dx
   10 continue
c deposit current for last particle
      nn = nnn + 1
      amx = qm*(.75 - dxn*dxn)
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
c deposit current
      vy = part(3,nop)
      vz = part(4,nop)
      dx1 = cu(1,nn) + vy*dxl
      dxl = cu(2,nn) + vz*dxl
      dx2 = cu(1,nn+1) + vy*amx
      amx = cu(2,nn+1) + vz*amx
      dx3 = cu(1,nn+2) + vy*dxp
      dxp = cu(2,nn+2) + vz*dxp
      cu(1,nn) = dx1
      cu(2,nn) = dxl
      cu(1,nn+1) = dx2
      cu(2,nn+1) = amx
      cu(1,nn+2) = dx3
      cu(2,nn+2) = dxp
c advance position half a time-step
      dx = part(1,nop) + part(2,nop)*dt
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
c set new position
      part(1,nop) = dx
      return
      end
      subroutine GSJPOST1X(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using second-order spline interpolation,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells
c 27 flops/particle, 10 loads, 7 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(.75-dx**2)
c cu(i,n+1)=.5*qci*(.5+dx)**2
c cu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,j+1) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+3
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2*nxv)
c local data
      integer npp, nn, npb, ipp, k, je, jb, j, n, n2, i
      real dmx, qmh, edgelx, edgerx, dxp, amx, dxl, vy, vz, dx
      parameter(npp=512)
      dimension nn(6,npp), dmx(6,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = real(nop - 1)/real(npb) + 1.
      qmh = .5*qm
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb) + .5
      dxp = part(1,j+jb) - real(n)
      n2 = 2*n + 1
      amx = qm*(.75 - dxp*dxp)
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      nn(1,j) = n2
      nn(2,j) = n2 + 1
      nn(3,j) = n2 + 2
      nn(4,j) = n2 + 3
      nn(5,j) = n2 + 4
      nn(6,j) = n2 + 5
      vy = part(3,j+jb)
      vz = part(4,j+jb)
      dmx(1,j) = vy*dxl
      dmx(2,j) = vz*dxl
      dmx(3,j) = vy*amx
      dmx(4,j) = vz*amx
      dmx(5,j) = vy*dxp
      dmx(6,j) = vz*dxp
c advance position half a time-step
      dx = part(1,j+jb) + part(2,j+jb)*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(2,j+jb) = -part(2,j+jb)
         endif
      endif
c set new position
      part(1,j+jb) = dx
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 6
      cu(nn(i,j)) = cu(nn(i,j)) + dmx(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 14 flops/particle, 8 loads, 5 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,j) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2,nxv)
c local data
      integer j, nn
      real edgelx, edgerx, dxp, amx, vy, vz, dx
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) 
      dxp = qm*(part(1,j) - real(nn))
      nn = nn + 1
      amx = qm - dxp
c deposit current
      vy = part(3,j)
      vz = part(4,j)
      cu(1,nn) = cu(1,nn) + vy*amx
      cu(2,nn) = cu(2,nn) + vz*amx
      cu(1,nn+1) = cu(1,nn+1) + vy*dxp
      cu(2,nn+1) = cu(2,nn+1) + vz*dxp
c advance position half a time-step
      dx = part(1,j) + part(2,j)*dt
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
      return
      end
      subroutine GSJPOST1L(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation
c 14 flops/particle, 8 loads, 5 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,j) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2,nxv)
c local data
      integer nnn, j, nn
      real dxn, edgelx, edgerx, dxp, amx, vy, vz, dx, dx1, dx2
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      dxn = qm*(part(1,1) - real(nnn))
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      nnn = part(1,j)
      dxp = dxn
      dxn = qm*(part(1,j) - real(nnn))
      amx = qm - dxp
c deposit current
      vy = part(3,j-1)
      vz = part(4,j-1)
      dx1 = cu(1,nn) + vy*amx
      amx = cu(2,nn) + vz*amx
      dx2 = cu(1,nn+1) + vy*dxp
      dxp = cu(2,nn+1) + vz*dxp
      cu(1,nn) = dx1
      cu(2,nn) = amx
      cu(1,nn+1) = dx2
      cu(2,nn+1) = dxp
c advance position half a time-step
      dx = part(1,j-1) + part(2,j-1)*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(2,j-1) = -part(2,j-1)
         endif
      endif
c set new position
      part(1,j-1) = dx
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      amx = qm - dxn
c deposit current
      vy = part(3,nop)
      vz = part(4,nop)
      dx1 = cu(1,nn) + vy*amx
      amx = cu(2,nn) + vz*amx
      dx2 = cu(1,nn+1) + vy*dxn
      dxp = cu(2,nn+1) + vz*dxn
      cu(1,nn) = dx1
      cu(2,nn) = amx
      cu(1,nn+1) = dx2
      cu(2,nn+1) = dxp
c advance position half a time-step
      dx = part(1,nop) + part(2,nop)*dt
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
c set new position
      part(1,nop) = dx
      return
      end
      subroutine GSJPOST1XL(part,cu,qm,dt,nop,idimp,nx,nxv,ipbc)
c for 1-2/2d code, this subroutine calculates particle current density
c using first-order linear interpolation,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells
c 14 flops/particle, 8 loads, 5 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vi, where i = y,z
c part(1,n) = position x of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c cu(i,j) = ith component of current density at grid point j
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of current array, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer nop, idimp, nx, nxv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2*nxv)
c local data
      integer npp, nn, npb, ipp, k, je, jb, j, n, n2, i
      real dmx, edgelx, edgerx, dxp, amx, vy, vz, dx
      parameter(npp=512)
      dimension nn(4,npp), dmx(4,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = real(nop - 1)/real(npb) + 1.
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb)
      dxp = qm*(part(1,j+jb) - real(n))
      n2 = 2*n + 1
      amx = qm - dxp
      nn(1,j) = n2
      nn(2,j) = n2 + 1
      nn(3,j) = n2 + 2
      nn(4,j) = n2 + 3
      vy = part(3,j+jb)
      vz = part(4,j+jb)
      dmx(1,j) = vy*amx
      dmx(2,j) = vz*amx
      dmx(3,j) = vy*dxp
      dmx(4,j) = vz*dxp
c advance position half a time-step
      dx = part(1,j+jb) + part(2,j+jb)*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(2,j+jb) = -part(2,j+jb)
         endif
      endif
c set new position
      part(1,j+jb) = dx
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 4
      cu(nn(i,j)) = cu(nn(i,j)) + dmx(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GBPUSH13(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,nx,n
     1xv,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 96 flops/particle, 1 divide, 19 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)),
c position equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (.75-dx**2)*fx(n)+.5*(fx(n+1)*(.5+dx)**2+fx(n-1)*(.5-dx)**2)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxyz(2,j+1,k+1) = y component of force/charge at grid (j,k)
c fxyz(3,j+1,k+1) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c byz(1,j+1,k+1) = y component of magnetic field at grid (j,k)
c byz(2,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxv = second dimension of field arrays, must be >= nx+3
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fxyz, byz, omx, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
c local data
      integer j, nn
      real qtmh, edgelx, edgerx, dxp, amx, dxl, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      dxp = part(1,j) - real(nn)
      nn = nn + 1
      amx = .75 - dxp*dxp
      dxl = .5*(.5 - dxp)**2
      dxp = .5*(.5 + dxp)**2
c find electric field
      dx = amx*fxyz(1,nn+1) + fxyz(1,nn)*dxl + fxyz(1,nn+2)*dxp
      dy = amx*fxyz(2,nn+1) + fxyz(2,nn)*dxl + fxyz(2,nn+2)*dxp
      dz = amx*fxyz(3,nn+1) + fxyz(3,nn)*dxl + fxyz(3,nn+2)*dxp
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn+1) + byz(1,nn)*dxl + byz(1,nn+2)*dxp
      oz = amx*byz(2,nn+1) + byz(2,nn)*dxl + byz(2,nn+2)*dxp
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,j) + dx
      acy = part(3,j) + dy
      acz = part(4,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,j) = dx
      part(3,j) = dy
      part(4,j) = dz
c new position
      dx = part(1,j) + dx*dtc
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
      ek = ek + .5*sum1
      return
      end
      subroutine GSBPUSH13(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,nx,
     1nxv,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 96 flops/particle, 1 divide, 19 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t)), and omz = (q/m)*bz(x(t)),
c position equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (.75-dx**2)*fx(n)+.5*(fx(n+1)*(.5+dx)**2+fx(n-1)*(.5-dx)**2)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxyz(2,j+1,k+1) = y component of force/charge at grid (j,k)
c fxyz(3,j+1,k+1) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c byz(1,j+1,k+1) = y component of magnetic field at grid (j,k)
c byz(2,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxv = second dimension of field arrays, must be >= nx+3
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fxyz, byz, omx, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
c local data
      integer nnn, nop1, j, nn
      real dxn, qtmh, edgelx, edgerx, amx, dxl, dxp, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1) + .5
      dxn = part(1,1) - real(nnn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      nnn = part(1,j+1) + .5
      dxp = dxn
      dxn = part(1,j+1) - real(nnn)
      amx = .75 - dxp*dxp
      dxl = .5*(.5 - dxp)**2
      dxp = .5*(.5 + dxp)**2
c find electric field
      dx = amx*fxyz(1,nn+1) + fxyz(1,nn)*dxl + fxyz(1,nn+2)*dxp
      dy = amx*fxyz(2,nn+1) + fxyz(2,nn)*dxl + fxyz(2,nn+2)*dxp
      dz = amx*fxyz(3,nn+1) + fxyz(3,nn)*dxl + fxyz(3,nn+2)*dxp
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn+1) + byz(1,nn)*dxl + byz(1,nn+2)*dxp
      oz = amx*byz(2,nn+1) + byz(2,nn)*dxl + byz(2,nn+2)*dxp
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,j) + dx
      acy = part(3,j) + dy
      acz = part(4,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,j) = dx
      part(3,j) = dy
      part(4,j) = dz
c new position
      dx = part(1,j) + dx*dtc
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
c push last particle
      nn = nnn + 1
      amx = .75 - dxn*dxn
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
c find electric field
      dx = amx*fxyz(1,nn+1) + fxyz(1,nn)*dxl + fxyz(1,nn+2)*dxp
      dy = amx*fxyz(2,nn+1) + fxyz(2,nn)*dxl + fxyz(2,nn+2)*dxp
      dz = amx*fxyz(3,nn+1) + fxyz(3,nn)*dxl + fxyz(3,nn+2)*dxp
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn+1) + byz(1,nn)*dxl + byz(1,nn+2)*dxp
      oz = amx*byz(2,nn+1) + byz(2,nn)*dxl + byz(2,nn+2)*dxp 
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,nop) + dx
      acy = part(3,nop) + dy
      acz = part(4,nop) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,nop) = dx
      part(3,nop) = dy
      part(4,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
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
c set new position
      part(1,nop) = dx
c normalize kinetic energy
   20 ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,nx,
     1nxv,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 78 flops/particle, 1 divide, 14 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t),y(t)), and omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(1,j,k) = x component of force/charge at grid (j,k)
c fxyz(2,j,k) = y component of force/charge at grid (j,k)
c fxyz(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c byz(1,j,k) = y component of magnetic field at grid (j,k)
c byz(2,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxv = second dimension of field arrays, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fxyz, byz, omx, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
c local data
      integer j, nn
      real qtmh, edgelx, edgerx, dxp, amx, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      dxp = part(1,j) - real(nn)
      nn = nn + 1
      amx = 1.0 - dxp
c find electric field
      dx = amx*fxyz(1,nn) + dxp*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + dxp*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + dxp*fxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn) + dxp*byz(1,nn+1)
      oz = amx*byz(2,nn) + dxp*byz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,j) + dx
      acy = part(3,j) + dy
      acz = part(4,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,j) = dx
      part(3,j) = dy
      part(4,j) = dz
c new position
      dx = part(1,j) + dx*dtc
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
      ek = ek + .5*sum1
      return
      end
      subroutine GSBPUSH13L(part,fxyz,byz,omx,qbm,dt,dtc,ek,idimp,nop,nx
     1,nxv,ipbc)
c for 1-2/2d code, this subroutine updates particle co-ordinate and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 78 flops/particle, 1 divide, 14 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fx(x(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fy(x(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt) + .5*(q/m)*fz(x(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omy = (q/m)*by(x(t),y(t)), and omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt
c fx(x(t)) is approximated by interpolation from the nearest grid points
c fx(x) = (1-dx)*fx(n)+dx*fx(n+1)
c where n = nearest grid point and dx = x-n
c similarly for fy(x), fz(x), by(x), bz(x)
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fxyz(1,j,k) = x component of force/charge at grid (j,k)
c fxyz(2,j,k) = y component of force/charge at grid (j,k)
c fxyz(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c byz(1,j,k) = y component of magnetic field at grid (j,k)
c byz(2,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c nxv = second dimension of field arrays, must be >= nx+1
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, nxv, ipbc
      real part, fxyz, byz, omx, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
c local data
      integer nnn, nop1, j, nn
      real dxn, qtmh, edgelx, edgerx, dxp, amx, dx, dy, dz
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      dxn = part(1,1) - real(nnn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      nnn = part(1,j+1)
      dxp = dxn
      dxn = part(1,j+1) - real(nnn)
      amx = 1.0 - dxp
c find electric field
      dx = amx*fxyz(1,nn) + dxp*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + dxp*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + dxp*fxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn) + dxp*byz(1,nn+1)
      oz = amx*byz(2,nn) + dxp*byz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,j) + dx
      acy = part(3,j) + dy
      acz = part(4,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,j) = dx
      part(3,j) = dy
      part(4,j) = dz
c new position
      dx = part(1,j) + dx*dtc
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
c push last particle
      nn = nnn + 1
      amx = 1.0 - dxn
c find electric field
      dx = amx*fxyz(1,nn) + dxn*fxyz(1,nn+1)
      dy = amx*fxyz(2,nn) + dxn*fxyz(2,nn+1)
      dz = amx*fxyz(3,nn) + dxn*fxyz(3,nn+1)
c find magnetic field
      ox = omx
      oy = amx*byz(1,nn) + dxn*byz(1,nn+1)
      oz = amx*byz(2,nn) + dxn*byz(2,nn+1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(2,nop) + dx
      acy = part(3,nop) + dy
      acz = part(4,nop) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(2,nop) = dx
      part(3,nop) = dy
      part(4,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
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
c set new position
      part(1,nop) = dx
c normalize kinetic energy
   20 ek = ek + .5*sum1
      return
      end
      subroutine RETARD1(part,dtc,idimp,nop,nx,ipbc)
c for 1-2/2d code, particle position is retarded a half time-step
c input: all, output: part
c equation used is:
c x(t+dt) = x(t) - vx(t+dt/2)*dtc
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c dtc = time interval between successive co-ordinate calculations
c idimp = size of phase space = 4
c nop = number of particles
c nx = system length in x direction
c ipbc = particle boundary condition = (0,1,2) =
c (none,2d periodic,2d reflecting)
      implicit none
      integer idimp, nop, nx, ipbc
      real part, dtc
      dimension part(idimp,nop)
c local data
      integer j
      real edgelx, edgerx, dx
c set boundary values
      edgelx = 0.0
      edgerx = real(nx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
c retard position half a time-step for current deposit
      dx = part(1,j) - part(2,j)*dtc
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
      return
      end
