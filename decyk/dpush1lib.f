c-----------------------------------------------------------------------
c 1d PIC library for pushing particles with darwin electric and magnetic
c fields and depositing current and derivative of current
c dpush1lib.f contains procedures to process particles with darwin
c             electric and magnetic fields:
c GMJPOST1 deposits momentum flux, quadratic interpolation, STANDARD
c          optimization.
c GSMJPOST1 deposits momentum flux, quadratic interpolation, LOOKAHEAD
c           optimization.
c GDCJPOST1 deposits momentum flux, acceleration density and current
c           density, quadratic interpolation, STANDARD optimization.
c GSDCJPOST1 deposits momentum flux, acceleration density and current
c            density, quadratic interpolation, LOOKAHEAD optimization.
c GMJPOST1L deposits momentum flux, linear interpolation, STANDARD
c           optimization.
c GSMJPOST1L deposits momentum flux, linear interpolation, LOOKAHEAD
c            optimization.
c GDCJPOST1L deposits momentum flux, acceleration density and current
c            density, linear interpolation, STANDARD optimization.
c GSDCJPOST1L deposits momentum flux, acceleration density and current
c             density, linear interpolation, LOOKAHEAD optimization.
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: august 13, 2009
c-----------------------------------------------------------------------
      subroutine GMJPOST1(part,amu,qm,nop,idimp,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation
c scalar version using guard cells
c 27 flops/particle, 10 loads, 6 stores
c input: all, output: part, amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(.75-dx**2)
c amu(i,n+1)=.5*qci*(.5+dx)**2
c amu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c amu(i,j+1) = ith component of momentum flux at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of flux array, must be >= nx+3
      implicit none
      integer nop, idimp, nxv
      real part, amu, qm
      dimension part(idimp,nop), amu(2,nxv)
      integer j, nn
      real qmh, dxp, amx, dxl, vx, v1, v2
      qmh = .5*qm
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) + .5
      dxp = part(1,j) - float(nn)
      nn = nn + 1
      amx = qm*(.75 - dxp*dxp)
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
c deposit momentum flux
      vx = part(2,j)
      v1 = vx*part(3,j)
      v2 = vx*part(4,j)
      amu(1,nn) = amu(1,nn) + v1*dxl
      amu(2,nn) = amu(2,nn) + v2*dxl
      amu(1,nn+1) = amu(1,nn+1) + v1*amx
      amu(2,nn+1) = amu(2,nn+1) + v2*amx
      amu(1,nn+2) = amu(1,nn+2) + v1*dxp
      amu(2,nn+2) = amu(2,nn+2) + v2*dxp
   10 continue
      return
      end
      subroutine GSMJPOST1(part,amu,qm,nop,idimp,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation
c scalar version using guard cells, integer conversion precalculation
c 27 flops/particle, 10 loads, 6 stores
c input: all, output: part, amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(.75-dx**2)
c amu(i,n+1)=.5*qci*(.5+dx)**2
c amu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c amu(i,j+1) = ith component of momentum flux at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of flux array, must be >= nx+3
      implicit none
      integer nop, idimp, nxv
      real part, amu, qm
      dimension part(idimp,nop), amu(2,nxv)
      integer j, nn, nnn
      real qmh, dxp, amx, dxl, vx, v1, v2
      real dxn, dx1, dx2, dx3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      dxn = part(1,1) - float(nnn)
      qmh = .5*qm
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      nnn = part(1,j) + .5
      dxp = dxn
      dxn = part(1,j) - float(nnn)
      amx = qm*(.75 - dxp*dxp)
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
c deposit momentum flux
      vx = part(2,j-1)
      v1 = vx*part(3,j-1)
      v2 = vx*part(4,j-1)
      dx1 = amu(1,nn) + v1*dxl
      dxl = amu(2,nn) + v2*dxl
      dx2 = amu(1,nn+1) + v1*amx
      amx = amu(2,nn+1) + v2*amx
      dx3 = amu(1,nn+2) + v1*dxp
      dxp = amu(2,nn+2) + v2*dxp
      amu(1,nn) = dx1
      amu(2,nn) = dxl
      amu(1,nn+1) = dx2
      amu(2,nn+1) = amx
      amu(1,nn+2) = dx3
      amu(2,nn+2) = dxp
   10 continue
c deposit momentum flux for last particle
      nn = nnn + 1
      amx = qm*(.75 - dxn*dxn)
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
c deposit momentum flux
      vx = part(2,nop)
      v1 = vx*part(3,nop)
      v2 = vx*part(4,nop)
      dx1 = amu(1,nn) + v1*dxl
      dxl = amu(2,nn) + v2*dxl
      dx2 = amu(1,nn+1) + v1*amx
      amx = amu(2,nn+1) + v2*amx
      dx3 = amu(1,nn+2) + v1*dxp
      dxp = amu(2,nn+2) + v2*dxp
      amu(1,nn) = dx1
      amu(2,nn) = dxl
      amu(1,nn+1) = dx2
      amu(2,nn+1) = amx
      amu(1,nn+2) = dx3
      amu(2,nn+2) = dxp
      return
      end
      subroutine GDCJPOST1(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,idimp,
     1nop,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation.
c scalar version using guard cells
c 139 flops/particle, 1 divide, 22 loads, 18 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(.75-dx**2)
c cu(i,n+1)=.5*qci*(.5+dx)**2
c cu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vj, where j = y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n)=qci*(.75-dx**2)
c dcu(i,n+1)=.5*qci*(.5+dx)**2
c dcu(i,n-1)=.5*qci*(.5-dx)**2
c and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(.75-dx**2)
c amu(i,n+1)=.5*qci*(.5+dx)**2
c amu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c velocity equations at t=t+dt/2 are calculated from:
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
c cu(i,j+1) = ith component of current density
c at grid point j for i = 1, 2
c dcu(i,j+1) = ith component of acceleration density
c at grid point j for i = 1, 2
c amu(i,j+1) = ith component of momentum flux
c at grid point j for i = 1, 2
c omx = magnetic field electron cyclotron frequency in x
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+3
      implicit none
      integer idimp, nop, nxv
      real part, fxyz, byz, cu, dcu, amu, omx, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension cu(2,nxv), dcu(2,nxv), amu(2,nxv)
      integer j, nn
      real qtmh, dti, dxp, amx, dxl, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      dxp = part(1,j) - float(nn)
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
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      amu(1,nn) = amu(1,nn) + v1*dxl
      amu(2,nn) = amu(2,nn) + v2*dxl
      amu(1,nn+1) = amu(1,nn+1) + v1*amx
      amu(2,nn+1) = amu(2,nn+1) + v2*amx
      amu(1,nn+2) = amu(1,nn+2) + v1*dxp
      amu(2,nn+2) = amu(2,nn+2) + v2*dxp
      dcu(1,nn) = dcu(1,nn) + vy*dxl
      dcu(2,nn) = dcu(2,nn) + vz*dxl
      dcu(1,nn+1) = dcu(1,nn+1) + vy*amx
      dcu(2,nn+1) = dcu(2,nn+1) + vz*amx
      dcu(1,nn+2) = dcu(1,nn+2) + vy*dxp
      dcu(2,nn+2) = dcu(2,nn+2) + vz*dxp
      cu(1,nn) = cu(1,nn) + oy*dxl
      cu(2,nn) = cu(2,nn) + oz*dxl
      cu(1,nn+1) = cu(1,nn+1) + oy*amx
      cu(2,nn+1) = cu(2,nn+1) + oz*amx
      cu(1,nn+2) = cu(1,nn+2) + oy*dxp
      cu(2,nn+2) = cu(2,nn+2) + oz*dxp
   10 continue
      return
      end
      subroutine GSDCJPOST1(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,idimp
     1,nop,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 139 flops/particle, 1 divide, 22 loads, 18 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(.75-dx**2)
c cu(i,n+1)=.5*qci*(.5+dx)**2
c cu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vj, where j = y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n)=qci*(.75-dx**2)
c dcu(i,n+1)=.5*qci*(.5+dx)**2
c dcu(i,n-1)=.5*qci*(.5-dx)**2
c and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(.75-dx**2)
c amu(i,n+1)=.5*qci*(.5+dx)**2
c amu(i,n-1)=.5*qci*(.5-dx)**2
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c velocity equations at t=t+dt/2 are calculated from:
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
c cu(i,j+1) = ith component of current density
c at grid point j for i = 1, 2
c dcu(i,j+1) = ith component of acceleration density
c at grid point j for i = 1, 2
c amu(i,j+1) = ith component of momentum flux
c at grid point j for i = 1, 2
c omx = magnetic field electron cyclotron frequency in x
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+3
      implicit none
      integer idimp, nop, nxv
      real part, fxyz, byz, cu, dcu, amu, omx, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension cu(2,nxv), dcu(2,nxv), amu(2,nxv)
      integer nop1, j, nn, nnn
      real qtmh, dti, dxn, dxp, amx, dxl, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, dx1
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      dxn = part(1,1) - float(nnn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      nnn = part(1,j+1) + .5
      dxp = dxn
      dxn = part(1,j+1) - float(nnn)
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
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      dx = amu(1,nn) + v1*dxl
      dy = amu(2,nn) + v2*dxl
      dz = amu(1,nn+1) + v1*amx
      vx = amu(2,nn+1) + v2*amx
      ox = amu(1,nn+2) + v1*dxp
      dx1 = amu(2,nn+2) + v2*dxp
      amu(1,nn) = dx
      amu(2,nn) = dy
      amu(1,nn+1) = dz
      amu(2,nn+1) = vx
      amu(1,nn+2) = ox
      amu(2,nn+2) = dx1
      dx = dcu(1,nn) + vy*dxl
      dy = dcu(2,nn) + vz*dxl
      dz = dcu(1,nn+1) + vy*amx
      vx = dcu(2,nn+1) + vz*amx
      ox = dcu(1,nn+2) + vy*dxp
      dx1 = dcu(2,nn+2) + vz*dxp
      dcu(1,nn) = dx
      dcu(2,nn) = dy
      dcu(1,nn+1) = dz
      dcu(2,nn+1) = vx
      dcu(1,nn+2) = ox
      dcu(2,nn+2) = dx1
      dx = cu(1,nn) + oy*dxl
      dy = cu(2,nn) + oz*dxl
      dz = cu(1,nn+1) + oy*amx
      vx = cu(2,nn+1) + oz*amx
      ox = cu(1,nn+2) + oy*dxp
      dx1 = cu(2,nn+2) + oz*dxp
      cu(1,nn) = dx
      cu(2,nn) = dy
      cu(1,nn+1) = dz
      cu(2,nn+1) = vx
      cu(1,nn+2) = ox
      cu(2,nn+2) = dx1
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
      vx = part(2,nop)
      vy = part(3,nop)
      vz = part(4,nop)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      amu(1,nn) = amu(1,nn) + v1*dxl
      amu(2,nn) = amu(2,nn) + v2*dxl
      amu(1,nn+1) = amu(1,nn+1) + v1*amx
      amu(2,nn+1) = amu(2,nn+1) + v2*amx
      amu(1,nn+2) = amu(1,nn+2) + v1*dxp
      amu(2,nn+2) = amu(2,nn+2) + v2*dxp
      dcu(1,nn) = dcu(1,nn) + vy*dxl
      dcu(2,nn) = dcu(2,nn) + vz*dxl
      dcu(1,nn+1) = dcu(1,nn+1) + vy*amx
      dcu(2,nn+1) = dcu(2,nn+1) + vz*amx
      dcu(1,nn+2) = dcu(1,nn+2) + vy*dxp
      dcu(2,nn+2) = dcu(2,nn+2) + vz*dxp
      cu(1,nn) = cu(1,nn) + oy*dxl
      cu(2,nn) = cu(2,nn) + oz*dxl
      cu(1,nn+1) = cu(1,nn+1) + oy*amx
      cu(2,nn+1) = cu(2,nn+1) + oz*amx
      cu(1,nn+2) = cu(1,nn+2) + oy*dxp
      cu(2,nn+2) = cu(2,nn+2) + oz*dxp
      return
      end
      subroutine GMJPOST1L(part,amu,qm,nop,idimp,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation
c scalar version using guard cells
c 14 flops/particle, 8 loads, 4 stores
c input: all, output: part, amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c amu(i,j) = ith component of momentum flux at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of flux array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, amu, qm
      dimension part(idimp,nop), amu(2,nxv)
      integer j, nn
      real dxp, amx, vx, v1, v2
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) 
      dxp = qm*(part(1,j) - float(nn))
      nn = nn + 1
      amx = qm - dxp
c deposit momentum flux
      vx = part(2,j)
      v1 = vx*part(3,j)
      v2 = vx*part(4,j)
      amu(1,nn) = amu(1,nn) + v1*amx
      amu(2,nn) = amu(2,nn) + v2*amx
      amu(1,nn+1) = amu(1,nn+1) + v1*dxp
      amu(2,nn+1) = amu(2,nn+1) + v2*dxp
   10 continue
      return
      end
      subroutine GSMJPOST1L(part,amu,qm,nop,idimp,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation
c scalar version using guard cells, integer conversion precalculation
c 14 flops/particle, 8 loads, 4 stores
c input: all, output: part, amcu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n) = position x of particle n
c part(2,n) = x velocity of particle n
c part(3,n) = y velocity of particle n
c part(4,n) = z velocity of particle n
c amu(i,j) = ith component of momentum flux at grid point j
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = second dimension of flux array, must be >= nx+1
      implicit none
      integer nop, idimp, nxv
      real part, amu, qm
      dimension part(idimp,nop), amu(2,nxv)
      integer j, nn, nnn
      real dxp, amx, vx, v1, v2
      real dxn, dx1, dx2
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      dxn = qm*(part(1,1) - float(nnn))
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      nnn = part(1,j)
      dxp = dxn
      dxn = qm*(part(1,j) - float(nnn))
      amx = qm - dxp
c deposit momentum flux
      vx = part(2,j-1)
      v1 = vx*part(3,j-1)
      v2 = vx*part(4,j-1)
      dx1 = amu(1,nn) + v1*amx
      amx = amu(2,nn) + v2*amx
      dx2 = amu(1,nn+1) + v1*dxp
      dxp = amu(2,nn+1) + v2*dxp
      amu(1,nn) = dx1
      amu(2,nn) = amx
      amu(1,nn+1) = dx2
      amu(2,nn+1) = dxp
   10 continue
c deposit momentum flux for last particle
      nn = nnn + 1
      amx = qm - dxn
c deposit momentum flux
      vx = part(2,nop)
      v1 = vx*part(3,nop)
      v2 = vx*part(4,nop)
      dx1 = amu(1,nn) + v1*amx
      amx = amu(2,nn) + v2*amx
      dx2 = amu(1,nn+1) + v1*dxn
      dxp = amu(2,nn+1) + v2*dxn
      amu(1,nn) = dx1
      amu(2,nn) = amx
      amu(1,nn+1) = dx2
      amu(2,nn+1) = dxp
      return
      end
      subroutine GDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,idimp
     1,nop,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c scalar version using guard cells
c 108 flops/particle, 1 divide, 16 loads, 12 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj, where j = y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c velocity equations at t=t+dt/2 are calculated from:
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
c cu(i,j) = ith component of current density at grid point j
c for i = 1, 2
c dcu(i,j) = ith component of acceleration density at grid point j
c for i = 1, 2
c amu(i,j) = ith component of momentum flux at grid point j
c for i = 1, 2
c omx = magnetic field electron cyclotron frequency in x
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
      implicit none
      integer idimp, nop, nxv
      real part, fxyz, byz, cu, dcu, amu, omx, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension cu(2,nxv), dcu(2,nxv), amu(2,nxv)
      integer j, nn
      real qtmh, dti, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      dxp = part(1,j) - float(nn)
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
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      amu(1,nn) = amu(1,nn) + v1*amx
      amu(2,nn) = amu(2,nn) + v2*amx
      amu(1,nn+1) = amu(1,nn+1) + v1*dxp
      amu(2,nn+1) = amu(2,nn+1) + v2*dxp
      dcu(1,nn) = dcu(1,nn) + vy*amx
      dcu(2,nn) = dcu(2,nn) + vz*amx
      dcu(1,nn+1) = dcu(1,nn+1) + vy*dxp
      dcu(2,nn+1) = dcu(2,nn+1) + vz*dxp
      cu(1,nn) = cu(1,nn) + oy*amx
      cu(2,nn) = cu(2,nn) + oz*amx
      cu(1,nn+1) = cu(1,nn+1) + oy*dxp
      cu(2,nn+1) = cu(2,nn+1) + oz*dxp
   10 continue
      return
      end
      subroutine GSDCJPOST1L(part,fxyz,byz,cu,dcu,amu,omx,qm,qbm,dt,idim
     1p,nop,nxv)
c for 1-2/2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 108 flops/particle, 1 divide, 16 loads, 12 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n)=qci*(1.-dx) and cu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj, where j = y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n)=qci*(1.-dx) and dcu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*dvj/dt, where j = y,z, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n)=qci*(1.-dx) and amu(i,n+1)=qci*dx
c where n = nearest grid point and dx = x-n
c and qci = qm*vj*vk, where jk = xy,xz for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c velocity equations at t=t+dt/2 are calculated from:
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
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = second dimension of field arrays, must be >= nx+1
      implicit none
      integer idimp, nop, nxv
      real part, fxyz, byz, cu, dcu, amu, omx, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxyz(3,nxv), byz(2,nxv)
      dimension cu(2,nxv), dcu(2,nxv), amu(2,nxv)
      integer nop1, j, nn, nnn
      real qtmh, dti, dxn, dxp, amx, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      dxn = part(1,1) - float(nnn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      nnn = part(1,j+1)
      dxp = dxn
      dxn = part(1,j+1) - float(nnn)
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
      vx = part(2,j)
      vy = part(3,j)
      vz = part(4,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      dx = amu(1,nn) + v1*amx
      dy = amu(2,nn) + v2*amx
      dz = amu(1,nn+1) + v1*dxp
      vx = amu(2,nn+1) + v2*dxp
      amu(1,nn) = dx
      amu(2,nn) = dy
      amu(1,nn+1) = dz
      amu(2,nn+1) = vx
      dx = dcu(1,nn) + vy*amx
      dy = dcu(2,nn) + vz*amx
      dz = dcu(1,nn+1) + vy*dxp
      vx = dcu(2,nn+1) + vz*dxp
      dcu(1,nn) = dx
      dcu(2,nn) = dy
      dcu(1,nn+1) = dz
      dcu(2,nn+1) = vx
      dx = cu(1,nn) + oy*amx
      dy = cu(2,nn) + oz*amx
      dz = cu(1,nn+1) + oy*dxp
      vx = cu(2,nn+1) + oz*dxp
      cu(1,nn) = dx
      cu(2,nn) = dy
      cu(1,nn+1) = dz
      cu(2,nn+1) = vx
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
      vx = part(2,nop)
      vy = part(3,nop)
      vz = part(4,nop)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
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
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxn = qm*dxn
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      v1 = ox*oy
      v2 = ox*oz
      amu(1,nn) = amu(1,nn) + v1*amx
      amu(2,nn) = amu(2,nn) + v2*amx
      amu(1,nn+1) = amu(1,nn+1) + v1*dxn
      amu(2,nn+1) = amu(2,nn+1) + v2*dxn
      dcu(1,nn) = dcu(1,nn) + vy*amx
      dcu(2,nn) = dcu(2,nn) + vz*amx
      dcu(1,nn+1) = dcu(1,nn+1) + vy*dxn
      dcu(2,nn+1) = dcu(2,nn+1) + vz*dxn
      cu(1,nn) = cu(1,nn) + oy*amx
      cu(2,nn) = cu(2,nn) + oz*amx
      cu(1,nn+1) = cu(1,nn+1) + oy*dxn
      cu(2,nn+1) = cu(2,nn+1) + oz*dxn
      return
      end
