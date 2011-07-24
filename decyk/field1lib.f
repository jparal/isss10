c-----------------------------------------------------------------------
c 1d PIC library for solving field equations
c field1lib.f contains procedures to manage guard cells and solve
c             fields equations in fourier space:
c CGUARD1 copy guard cells for 2 component vector array, quadratic
c         interpolation.
c BGUARD1 copy guard cells for 3 component vector array, quadratic
c         interpolation.
c DGUARD1 copy guard cells for scalar array, quadratic interpolation.
c SCGUARD1 initialize field for 2 component vector array, quadratic
c          interpolation.
c SGUARD1 initialize field for scalar array, quadratic interpolation.
c ACGUARD1 add guard cells for 2 component vector array, quadratic
c          interpolation.
c AGUARD1 add guard cells for scalar array, quadratic interpolation.
c CGUARD1L copy guard cells for 2 component vector array, linear
c          interpolation.
c BGUARD1L copy guard cells for 3 component vector array, linear
c          interpolation.
c DGUARD1L copy guard cells for scalar array, linear interpolation.
c SCGUARD1L initialize field for 2 component vector array, linear
c           interpolation.
c SGUARD1L initialize field for scalar array, linear interpolation.
c ACGUARD1L add guard cells for 2 component vector array, linear
c           interpolation.
c AGUARD1L add guard cells for scalar array, linear interpolation.
c POISP1 solve poisson equation for electric force, potential, or
c        smoothing.
c BPOIS13 solve vector poisson equation for magnetic force, vector
c         potential, or smoothing.
c IBPOIS13 solve vector poisson equation for magnetic field.
c MAXWEL1 solve maxwell equation for electric and magnetic fields.
c EMFIELD1 calculate electric force from electric fields given by
c          maxwell and poisson equations.
c BMFIELD1 calculate magnetic force from magnetic field given by maxwell
c          equation.
c EMFIELDR1 calculate electric force from electric fields given by
c           maxwell and poisson equations for real arrays.
c BMFIELDR1 calculate magnetic force from magnetic field given by maxwell
c           equation for real arrays.
c AVPOT13 calculate vector potential from magnetic field.
c AVRPOT13 calculate radiative part of the vector potential
c GTMODES1 extracts selected fourier components from potential array.
c PTMODES1 places selected fourier components into potential array.
c GTVMODES1 extracts selected fourier components from vector potential
c           array.
c PTVMODES1 places selected fourier components into vector potential
c           array.
c SCFGUARD1 initialize 2 component field with scaled vector array,
c           quadratic interpolation.
c SCFGUARD1L initialize 2 component field with scaled vector array,
c            linear interpolation.
c DCUPERP13 calculate transverse derivative of current density from
c           momentum flux.
c ADCUPERP13 calculate transverse derivative of current density from
c            momentum flux and acceleration density.
c EPOIS13 solve vector poisson equation for transverse electric field
c         or force.
c WPMXN1 calculates maximum and minimum plasma frequency.
c ADDQEI1 adds electron and ion densities.
c ADDQEI1X adds electron and ion densities and calculates maximum and
c          minimum plasma frequency.
c BADDEXT1 adds constant to magnetic field in real space for 1-2/2d code
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: july 16, 2011
c-----------------------------------------------------------------------
      subroutine CGUARD1(byz,nx,nxe)
c replicate extended periodic field
c quadratic interpolation
      implicit none
      real byz
      integer nx, nxe
      dimension byz(2,nxe)
      integer i
      do 10 i = 1, 2
      byz(i,1) = byz(i,nx+1)
      byz(i,nx+2) = byz(i,2)
      byz(i,nx+3) = byz(i,3)
   10 continue
      return
      end
      subroutine BGUARD1(fxyz,nx,nxe)
c replicate extended periodic field
c quadratic interpolation
      implicit none
      real fxyz
      integer nx, nxe
      dimension fxyz(3,nxe)
      integer i
      do 10 i = 1, 3
      fxyz(i,1) = fxyz(i,nx+1)
      fxyz(i,nx+2) = fxyz(i,2)
      fxyz(i,nx+3) = fxyz(i,3)
   10 continue
      return
      end
      subroutine DGUARD1(fx,nx,nxe)
c replicate extended periodic field
c quadratic interpolation
      implicit none
      real fx
      integer nx, nxe
      dimension fx(nxe)
      fx(1) = fx(nx+1)
      fx(nx+2) = fx(2)
      fx(nx+3) = fx(3)
      return
      end
      subroutine SCGUARD1(cu,yj0,zj0,nx,nxe)
c initialize extended periodic field
c quadratic interpolation
      implicit none
      real cu, yj0, zj0
      integer nx, nxe
      dimension cu(2,nxe)
      integer i, j
      do 10 j = 1, nx
      cu(1,j+1) = yj0
      cu(2,j+1) = zj0
   10 continue
      do 20 i = 1, 2
      cu(i,1) = 0.
      cu(i,nx+2) = 0.
      cu(i,nx+3) = 0.
   20 continue
      return
      end
      subroutine SGUARD1(q,qi0,nx,nxe)
c initialize extended periodic field
c quadratic interpolation
      implicit none
      real q, qi0
      integer nx, nxe
      dimension q(nxe)
      integer j
      do 10 j = 1, nx
      q(j+1) = qi0
   10 continue
      q(1) = 0.
      q(nx+2) = 0.
      q(nx+3) = 0.
      return
      end
      subroutine ACGUARD1(cu,nx,nxe)
c accumulate extended periodic field
c quadratic interpolation
      implicit none
      real cu
      integer nx, nxe
      dimension cu(2,nxe)
      integer i
      do 10 i = 1, 2
      cu(i,2) = cu(i,2) + cu(i,nx+2)
      cu(i,3) = cu(i,3) + cu(i,nx+3)
      cu(i,nx+1) = cu(i,nx+1) + cu(i,1)
      cu(i,nx+2) = 0.0
      cu(i,nx+3) = 0.0
      cu(i,1) = 0.0
   10 continue
      return
      end
      subroutine AGUARD1(q,nx,nxe)
c accumulate extended periodic field
c quadratic interpolation
      implicit none
      real q
      integer nx, nxe
      dimension q(nxe)
      q(2) = q(2) + q(nx+2)
      q(3) = q(3) + q(nx+3)
      q(nx+1) = q(nx+1) + q(1)
      q(nx+2) = 0.0
      q(nx+3) = 0.0
      q(1) = 0.0
      return
      end
      subroutine CGUARD1L(byz,nx,nxe)
c replicate extended periodic field
c linear interpolation
      implicit none
      real byz
      integer nx, nxe
      dimension byz(2,nxe)
      integer i
      do 10 i = 1, 2
      byz(i,nx+1) = byz(i,1)
   10 continue
      return
      end
      subroutine BGUARD1L(fxyz,nx,nxe)
c replicate extended periodic field
c linear interpolation
      implicit none
      real fxyz
      integer nx, nxe
      dimension fxyz(3,nxe)
      integer i
      do 10 i = 1, 3
      fxyz(i,nx+1) = fxyz(i,1)
   10 continue
      return
      end
      subroutine DGUARD1L(fx,nx,nxe)
c replicate extended periodic field
c linear interpolation
      implicit none
      real fx
      integer nx, nxe
      dimension fx(nxe)
      fx(nx+1) = fx(1)
      return
      end
      subroutine SCGUARD1L(cu,yj0,zj0,nx,nxe)
c initialize extended periodic field
c linear interpolation
      implicit none
      real cu, yj0, zj0
      integer nx, nxe
      dimension cu(2,nxe)
      integer i, j
      do 10 j = 1, nx
      cu(1,j) = yj0
      cu(2,j) = zj0
   10 continue
      do 20 i = 1, 2
      cu(i,nx+1) = 0.
   20 continue
      return
      end
      subroutine SGUARD1L(q,qi0,nx,nxe)
c initialize extended periodic field
c linear interpolation
      implicit none
      real q, qi0
      integer nx, nxe
      dimension q(nxe)
      integer j
      do 10 j = 1, nx
      q(j) = qi0
   10 continue
      q(nx+1) = 0.
      return
      end
      subroutine ACGUARD1L(cu,nx,nxe)
c accumulate extended periodic field
c linear interpolation
      implicit none
      real cu
      integer nx, nxe
      dimension cu(2,nxe)
      integer i
      do 10 i = 1, 2
      cu(i,1) = cu(i,1) + cu(i,nx+1)
      cu(i,nx+1) = 0.0
   10 continue
      return
      end
      subroutine AGUARD1L(q,nx,nxe)
c accumulate extended periodic field
c linear interpolation
      implicit none
      real q
      integer nx, nxe
      dimension q(nxe)
      q(1) = q(1) + q(nx+1)
      q(nx+1) = 0.0
      return
      end
      subroutine POISP1(q,fx,isign,ffc,ax,affp,we,nx)
c this subroutine solves 1d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,affp,nx, output: ffc
c for isign = +-1, input: q,ffc,isign,nx, output: fx,we
c for isign < 0, approximate flop count is: 6*nx
c for isign = 1, approximate flop count is: 3*nx
c for isign = 2, input: q,ffc,isign,nx, output: fx, flop count = nx
c if isign < 0, force/charge is calculated using the equations:
c fx(k) = -sqrt(-1)*k*g(k)*s(k)*q(k), where k = 2pi*j/nx, j=fourier mode,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-(k*ax)**2/2), except for
c fx(k=0) = fx(k=pi) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(k) = g(k)*q(k)*s(k)
c if isign = 2, smoothing is calculated using the equation:
c fx(k) = q(k)*s(k)
c cmplx(q(2*j-1),q(2*j)) = complex charge density for fourier mode j-1
c cmplx(fx(2*j-1),fx(2*j)) = complex force/charge for fourier mode j-1
c if isign = 0, form factor array is prepared
c ffc(2*j) = finite-size particle shape factor s for fourier mode j-1
c ffc(2*j-1) = potential green's function g for fourier mode j-1
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np = number of particles
c electric field energy is also calculated, using
c we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
c nx = system length in x direction
      implicit none
      integer isign, nx
      real ax, affp, we
      real q, fx, ffc
      dimension q(nx), fx(nx), ffc(nx)
c local data
      integer j, nxh
      real dnx, dkx, at1, at2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      ffc(2*j) = exp(-.5*(dkx*ax)**2)
      ffc(2*j-1) = affp*ffc(2*j)/(dkx*dkx)
   10 continue
      ffc(1) = affp
      ffc(2) = 1.0
      return
   20 if (isign.gt.0) go to 40
c calculate force/charge and sum field energy
      wp = 0.0d0
      do 30 j = 2, nxh
      at1 = ffc(2*j-1)*ffc(2*j)
      at2 = dnx*real(j - 1)*at1
      fx(2*j-1) = at2*q(2*j)
      fx(2*j) = -at2*q(2*j-1)
      wp = wp + at1*(q(2*j-1)**2 + q(2*j)**2)
   30 continue
      fx(1) = 0.
      fx(2) = 0.
      we = real(nx)*wp
      return
c calculate potential and sum field energy
   40 if (isign.gt.1) go to 60
      wp = 0.0d0
      do 50 j = 2, nxh
      at2 = ffc(2*j-1)
      at1 = at2*ffc(2*j)
      fx(2*j-1) = at2*q(2*j-1)
      fx(2*j) = at2*q(2*j)
      wp = wp + at1*(q(2*j-1)**2 + q(2*j)**2)
   50 continue
      fx(1) = 0.
      fx(2) = 0.
      we = real(nx)*wp
      return
c calculate smoothing
   60 do 70 j = 2, nxh
      at1 = ffc(2*j)
      fx(2*j-1) = at1*q(2*j-1)
      fx(2*j) = at1*q(2*j)
   70 continue
      fx(1) = ffc(2)*q(1)
      fx(2) = 0.
      end
      subroutine BPOIS13(cu,byz,isign,ffc,ax,affp,ci,wm,nx,nxvh,nxhd)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,affp,nx,nxvh, output: ffc
c for isign = -1, input: cu,ffc,isign,ci,nx,nxvh, output: byz,wm
c approximate flop count is: 30*nxc
c for isign = 1, input: cu,ffc,isign,ci,nx,nxvh, output: byz,wm
c approximate flop count is: 23*nxc
c for isign = 2, input: cu,ffc,isign,nx,nxvh, output: byz
c approximate flop count is: 4*nxc
c where nxc = nx/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx)*s(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx)*s(kx),
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2+)/2), except for
c by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
c if isign = 1, vector potential is calculated using the equation:
c by(kx) = ci*ci*g(kx)*cuy(kx)
c bz(kx) = ci*ci*g(kx)*cuz(kx)
c if isign = 2, smoothing is calculated using the equation:
c by(kx) = cuy(kx)*s(kx)
c bz(kx) = cuz(kx)*s(kx)
c cu(i,j) = complex current density for fourier mode (j-1)
c byz(1,j) = y component of complex magnetic field
c byz(2,j) = z component of complex magnetic field
c all for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2)
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = dimension of form factor array, must be >= nxh
      implicit none
      integer isign, nx, nxvh, nxhd
      real ax, affp, ci, wm
      complex cu, byz, ffc
      dimension cu(2,nxvh), byz(2,nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real dnx, ci2, dkx, at1, at2
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at2 = exp(-.5*(dkx*ax)**2)
      ffc(j) = cmplx(affp*at2/(dkx*dkx),at2)
   10 continue
      ffc(1) = cmplx(affp,1.0)
      return
   20 if (isign.gt.0) go to 40
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 30 j = 2, nxh
      at1 = ci2*real(ffc(j))*aimag(ffc(j))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt2 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
      byz(1,j) = -at2*zt1
      byz(2,j) = at2*zt2
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   30 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = real(nx)*wp
      return
c calculate vector potential and sum field energy
   40 if (isign.gt.1) go to 60
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 50 j = 2, nxh
      at2 = ci2*real(ffc(j))
      at1 = at2*aimag(ffc(j))
      byz(1,j) = at2*cu(1,j)
      byz(2,j) = at2*cu(2,j)
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   50 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = real(nx)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2
   60 do 70 j = 2, nxh
      at1 = aimag(ffc(j))
      byz(1,j) = at1*cu(1,j)
      byz(2,j) = at1*cu(2,j)
   70 continue
      at1 = aimag(ffc(1))
      byz(1,1) = cmplx(at1*real(cu(1,1)),0.)
      byz(2,1) = cmplx(at1*real(cu(2,1)),0.)
      return
      end
      subroutine IBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c magnetic field, with periodic boundary conditions.
c input: cu,ffc,ci,nx,nxv, output: byz,wm
c approximate flop count is: 29*nxc
c where nxc = nx/2 - 1
c the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx),
c where kx = 2pi*j/nx, and j = fourier mode number,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2)/2), except for
c by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
c cu(i,j) = complex current density for fourier mode (j-1)
c byz(i,j) = i component of complex magnetic field
c all for fourier mode (j-1)
c aimag(ffc(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffc(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2), where
c affp = normalization constant = nx/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, wm
      complex cu, byz, ffc
      dimension cu(2,nxvh), byz(2,nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real dnx, ci2, at1, at2
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j))
      zt1 = cmplx(-aimag(cu(2,j)),real(cu(2,j)))
      zt2 = cmplx(-aimag(cu(1,j)),real(cu(1,j)))
      byz(1,j) = -at2*zt1
      byz(2,j) = at2*zt2
      wp = wp + at1*(cu(1,j)*conjg(cu(1,j)) + cu(2,j)*conjg(cu(2,j)))
   10 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = real(nx)*wp
      return
      end
      subroutine MAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, eyz, byz
c approximate flop count is: 87*nxc
c where nxc = nx/2 - 1
c the magnetic field is first updated half a step using the equations:
c by(kx) = by(kx) + .5*dt*sqrt(-1)*kx*ez(kx)
c bz(kx) = bz(kx) - .5*dt*sqrt(-1)*kx*ey(kx)
c the electric field is then updated a whole step using the equations:
c ey(kx) = ey(kx) - c2*dt*sqrt(-1)*kx*bz(kx) - affp*dt*cuy(kx)*s(kx)
c ez(kx) = ez(kx) + c2*dt*sqrt(-1)*kx*by(kx) - affp*dt*cuz(kx)*s(kx)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, c2 = 1./(ci*ci)
c and s(kx) = exp(-((kx*ax)**2)
c j = fourier mode numbers, except for
c ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0) = 0.
c and similarly for by, bz.
c cu(i,j,k) = complex current density
c eyz(i,j,k) = complex transverse electric field
c byz(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1)
c real(ffc(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c aimag(ffc(j)) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*sum((1/affp)*|eyz(kx)|**2)
c magnetic field energy is also calculated, using
c wm = nx*sum((c2/affp)*|byz(kx)|**2)
c nx = system length in x direction
c nxvh = first dimension of field arrays, must be >= nxh
c nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci, dt, wf, wm
      complex eyz, byz, cu, ffc
      dimension eyz(2,nxvh), byz(2,nxvh), cu(2,nxvh), ffc(nxhd)
c local data
      integer j, nxh
      real dnx, dth, c2, cdt, affp, adt, anorm, dkx, afdt
      complex zero, zt1, zt2, zt5, zt6, zt8, zt9
      double precision wp, ws
      if (ci.le.0.) return
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j))
c update magnetic field half time step
      zt1 = cmplx(-aimag(eyz(2,j)),real(eyz(2,j)))
      zt2 = cmplx(-aimag(eyz(1,j)),real(eyz(1,j)))
      zt5 = byz(1,j) + dth*(dkx*zt1)
      zt6 = byz(2,j) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = eyz(1,j) - cdt*(dkx*zt1) - afdt*cu(1,j)
      zt9 = eyz(2,j) + cdt*(dkx*zt2) - afdt*cu(2,j)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      eyz(1,j) = zt8
      eyz(2,j) = zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      byz(1,j) = zt5
      byz(2,j) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
      byz(1,1) = zero
      byz(2,1) = zero
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*ws
      wm = real(nx)*c2*wp
      return
      end
      subroutine EMFIELD1(fxyz,fx,eyz,ffc,nx,nxvh,nxhd)
c this subroutine merges complex vector fields
c includes additional smoothing
      implicit none
      integer nx, nxvh, nxhd
      complex fxyz, fx, eyz, ffc
      dimension fxyz(3,nxvh), fx(nxvh), eyz(2,nxvh)
      dimension ffc(nxhd)
      integer j, nxh
      real at1
      nxh = nx/2
c add the fields
      do 10 j = 1, nxh
      at1 = aimag(ffc(j))
      fxyz(1,j) = fx(j)
      fxyz(2,j) = eyz(1,j)*at1
      fxyz(3,j) = eyz(2,j)*at1
   10 continue
      return
      end
      subroutine BMFIELD1(fyz,eyz,ffc,nx,nxvh,nxhd)
c this subroutine copies complex vector fields
c includes additional smoothing
      implicit none
      integer nx, nxvh, nxhd
      complex fyz, eyz, ffc
      dimension fyz(2,nxvh), eyz(2,nxvh)
      dimension ffc(nxhd)
      integer j, nxh
      real at1
      nxh = nx/2
      do 10 j = 1, nxh
      at1 = aimag(ffc(j))
      fyz(1,j) = eyz(1,j)*at1
      fyz(2,j) = eyz(2,j)*at1
   10 continue
      return
      end
      subroutine EMFIELDR1(fxyz,fx,eyz,ffc,nx,nxe,nxd)
c this subroutine merges real vector fields
c includes additional smoothing
c nxe >= nx+1
      implicit none
      integer nx, nxe, nxd
      real fxyz, fx, eyz
      complex ffc
      dimension fxyz(3,nxe), fx(nxe), eyz(2,nxe)
      dimension ffc(nxd)
      integer j
      real at1
c add the fields
      do 10 j = 1, nx
      at1 = aimag(ffc(j))
      fxyz(1,j) = fx(j)
      fxyz(2,j) = eyz(1,j)*at1
      fxyz(3,j) = eyz(2,j)*at1
   10 continue
      fxyz(1,nx+1) = fx(nx+1)
      fxyz(2,nx+1) = eyz(1,nx+1)
      fxyz(3,nx+1) = eyz(2,nx+1)
      return
      end
      subroutine BMFIELDR1(fyz,eyz,ffc,nx,nxe,nxd)
c this subroutine copies real vector fields
c includes additional smoothing
c nxe >= nx+1
      implicit none
      integer nx, nxe, nxd
      real fyz, eyz
      complex ffc
      dimension fyz(2,nxe), eyz(2,nxe)
      dimension ffc(nxd)
      integer j
      real at1
      do 10 j = 1, nx
      at1 = aimag(ffc(j))
      fyz(1,j) = eyz(1,j)*at1
      fyz(2,j) = eyz(2,j)*at1
   10 continue
      fyz(1,nx+1) = eyz(1,nx+1)
      fyz(2,nx+1) = eyz(2,nx+1)
      return
      end
      subroutine AVPOT13(byz,ayz,nx,nxvh)
c this subroutine calculates 1-2/2d vector potential from magnetic field
c in fourier space with periodic boundary conditions.
c input: byz, nx, nxvh, output: ayz
c approximate flop count is: 10*nxc and nxc divides,
c where nxc = nx/2 - 1
c the vector potential is calculated using the equations:
c ay(kx) = -sqrt(-1)*bz(kx))/kx
c az(kx) = sqrt(-1)*by(kx)/kx
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c ay(kx=pi) = az(kx=pi) = 0, and ay(kx=0) = az(kx=0) = 0.
c byz(i,j,k) = i component of complex magnetic field
c ayz(i,j,k) = i component of complex vector potential
c all for fourier mode (j-1)
c nx = system length in x direction
c nxvh = first dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex byz, ayz
      dimension byz(2,nxvh), ayz(2,nxvh)
c local data
      integer j, nxh
      real dnx, dkx, at2
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.,0.)
c calculate vector potential
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at2 = 1.0/dkx
      zt1 = cmplx(-aimag(byz(2,j)),real(byz(2,j)))
      zt2 = cmplx(-aimag(byz(1,j)),real(byz(1,j)))
      ayz(1,j) = -at2*zt1
      ayz(2,j) = at2*zt2
   10 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      return
      end
      subroutine AVRPOT13(ayz,byz,ffc,ci,nx,nxvh,nxhd)
c this subroutine solves 1-1/2d poisson's equation in fourier space for
c the radiative part of the vector potential
c with periodic boundary conditions.
c input: all, output: ayz
c approximate flop count is: 24*nxc
c where nxc = nx/2 - 1
c the radiative vector potential is updated using the equations:
c ay(kx) = -(sqrt(-1)*kx*bz(kx) + affp*ci2*cuy(kx)*s(kx))/(kx*kx)
c az(kx) = (sqrt(-1)*kx*by(kx) - affp*ci2*cuz(kx)*s(kx))/(kx*kx)
c where kx = 2pi*j/nx, ci2 = ci*ci
c and s(kx) = exp(-((kx*ax)**2)
c j = fourier mode numbers, except for
c ay(kx=pi) = az(kx=pi) = 0, and ay(kx=0) = az(kx=0) = 0.
c ayz(i,j) = on entry, complex current density cu
c ayz(i,j) = on exit, complex current radiative vector potential
c byz(i,j) = complex magnetic field
c for component i, all for fourier mode (j-1)
c real(ffc(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c aimag(ffc()) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c nx = system length in x direction
c nxvh = first dimension of field arrays, must be >= nxh
c nxhd = first dimension of form factor array, must be >= nxh
      implicit none
      integer nx, nxvh, nxhd
      real ci
      complex ayz, byz, ffc
      dimension ayz(2,nxvh), byz(2,nxvh)
      dimension ffc(nxhd)
c local data
      integer nxh, j
      real dnx, afc2, dkx, at1, at2
      complex zero, zt1, zt2
      if (ci.le.0.) return
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      afc2 = real(ffc(1))*ci*ci
      zero = cmplx(0.,0.)
c calculate the radiative vector potential
c mode numbers 0 < kx < nx/2
      do 20 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1.0/(dkx*dkx)
      at2 = afc2*aimag(ffc(j))
c update radiative vector potential
      zt1 = cmplx(-aimag(byz(2,j)),real(byz(2,j)))
      zt2 = cmplx(-aimag(byz(1,j)),real(byz(1,j)))
      ayz(1,j) = -at1*(dkx*zt1 + at2*ayz(1,j))
      ayz(2,j) = at1*(dkx*zt2 - at2*ayz(2,j))
   20 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      return
      end
      subroutine GTMODES1(pot,pott,nx,it,modesx,nxe,nt2,modesxd)
c this subroutine extracts lowest order modes from complex array pot
c and stores them into a location in a time history array pott
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c nxe = dimension of input array pot, nxe >= nx
c nt2 = first dimension of output array pott, nt2 >= 2*it
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, nxe, nt2
      integer modesxd
      real pot, pott
      dimension pot(nxe), pott(nt2,modesxd)
c local data
      integer i1, i2, nxh, jmax, j, j1
      i2 = it + it
      i1 = i2 - 1
      if (i2.gt.nt2) return
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmax
      pott(i1,j) = pot(2*j-1)
      pott(i2,j) = pot(2*j)
   10 continue
c mode numbers kx = 0, nx/2
      pott(i1,1) = pot(1)
      pott(i2,1) = 0.
      if (modesx.gt.nxh) then
         pott(i1,j1) = pot(2)
         pott(i2,j1) = 0.
      endif
      return
      end
      subroutine PTMODES1(pot,pott,nx,it,modesx,nxe,nt2,modesxd)
c this subroutine extracts lowest order modes from a location in a time
c history array pott and stores them into complex array pot
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c nxe = dimension of input array pot, nxe >= nx
c nt2 = first dimension of output array pott, nt2 >= 2*it
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, nxe, nt2
      integer modesxd
      real pot, pott
      dimension pot(nxe), pott(nt2,modesxd)
c local data
      integer i1, i2, nxh, jmax, j,j1
      i2 = it + it
      i1 = i2 - 1
      if (i2.gt.nt2) return
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 
      do 10 j = 2, jmax
      pot(2*j-1) = pott(i1,j)
      pot(2*j) = pott(i2,j)
   10 continue
      do 20 j = jmax+1, nxh
      pot(2*j-1) = 0.
      pot(2*j) = 0.
   20 continue
      pot(1) = pott(i1,1)
      pot(2) = 0.
      if (modesx.gt.nxh) then
         pot(2) = pott(i1,j1)
      endif
      return
      end
      subroutine GTVMODES1(vpot,vpott,nx,it,modesx,ndim,nxvh,nt,modesxd)
c this subroutine extracts lowest order modes from complex vector array
c vpot and stores them into a location in a time history array vpott
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c nt = first dimension of output array vpott, nt >= it
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, ndim, nxvh, nt
      integer modesxd
      complex vpot, vpott
      dimension vpot(ndim,nxvh), vpott(nt,ndim,modesxd)
c local data
      integer nxh, jmax, i, j, j1
      if (it.gt.nt) return
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpott(it,i,j) = vpot(i,j)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 i = 1, ndim
      vpott(it,i,1) = cmplx(real(vpot(i,1)),0.0)
      if (modesx.gt.nxh) then
         vpott(it,i,j1) = cmplx(aimag(vpot(i,1)),0.0)
      endif
   30 continue
      return
      end
      subroutine PTVMODES1(vpot,vpott,nx,it,modesx,ndim,nxvh,nt,modesxd)
c this subroutine extracts lowest order modes from a location in a time
c history array vpott and stores them into complex vector array vpot
c modes stored: kx=(0,1,...,NX/2)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c nt = first dimension of output array vpott, nt >= it
c modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, ndim, nxvh, nt
      integer modesxd
      complex vpot, vpott
      dimension vpot(ndim,nxvh), vpott(nt,ndim,modesxd)
      integer nxh, jmax, i, j, j1
      complex zero
      if (it.gt.nt) return
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmax = min0(modesx,nxh)
      j1 = nxh + 1
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpot(i,j) = vpott(it,i,j)
   10 continue
   20 continue
      do 40 j = jmax+1, nxh
      do 30 i = 1, ndim
      vpot(i,j) = zero
   30 continue
   40 continue
c mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      vpot(i,1) = cmplx(real(vpott(it,i,1)),0.0)
      if (modesx.gt.nxh) then
         vpot(i,1) = cmplx(real(vpot(i,1)),real(vpott(it,i,j1)))
      endif
   50 continue
      return
      end
      subroutine SCFGUARD1(cus,cu,q2m0,nx,nxe)
c initialize extended periodic field with scaled field
c quadratic interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, nxe
      dimension cus(2,nxe), cu(2,nxe)
      integer i, j
c initialize extended field, with zero in the edges
      do 20 j = 1, nx
      do 10 i = 1, 2
      cus(i,j+1) = -q2m0*cu(i,j+1)
   10 continue
   20 continue
      do 30 i = 1, 2
      cus(i,1) = 0.
      cus(i,nx+2) = 0.
      cus(i,nx+3) = 0.
   30 continue
      return
      end
      subroutine SCFGUARD1L(cus,cu,q2m0,nx,nxe)
c initialize extended periodic field with scaled field
c linear interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, nxe
      dimension cus(2,nxe), cu(2,nxe)
      integer i, j
c initialize extended field, with zero in the edges
      do 20 j = 1, nx
      do 10 i = 1, 2
      cus(i,j) = -q2m0*cu(i,j)
   10 continue
   20 continue
      do 30 i = 1, 2
      cus(i,nx+1) = 0.
   30 continue
      return
      end
      subroutine DCUPERP13(dcu,amu,nx,nxvh)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 1-2/2d with periodic boundary conditions.
c the transverse part of the derivative of the current is calculated
c using the equations:
c dcu(1,kx) = -sqrt(-1)*kx*vx*vy
c dcu(2,kx) = -sqrt(-1)*kx*vx*vz
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,kx=0) = 0.
c amu(1,j) = xy component of complex momentum flux
c amu(2,j) = xz component of complex momentum flux
c all for fourier mode (j-1)
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex dcu, amu
      dimension dcu(2,nxvh), amu(2,nxvh)
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(1,j)),-real(amu(1,j)))
      dcu(1,j) = dkx*zt2
      zt1 = cmplx(aimag(amu(2,j)),-real(amu(2,j)))
      dcu(2,j) = dkx*zt1
   10 continue
      dcu(1,1) = zero
      dcu(2,1) = zero
      return
      end
      subroutine ADCUPERP13(dcu,amu,nx,nxvh)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 1-2/2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx) = dcu(1,kx)-sqrt(-1)*kx*vx*vy
c dcu(2,kx) = dcu(2,kx)-sqrt(-1)*kx*vx*vz
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,kx=0) = 0.
c on input:
c dcu(i,j,k) = complex acceleration density for fourier mode (j-1)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1)
c amu(1,j,k) = xy component of complex momentum flux
c amu(2,j,k) = xzx component of complex momentum flux
c all for fourier mode (j-1)
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
      implicit none
      integer nx, nxvh
      complex dcu, amu
      dimension dcu(2,nxvh), amu(2,nxvh)
      integer nxh, j
      real dnx, dkx
      complex zero, zt1, zt2
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      zt2 = cmplx(aimag(amu(1,j)),-real(amu(1,j)))
      dcu(1,j) = dcu(1,j) + dkx*zt2
      zt1 = cmplx(aimag(amu(2,j)),-real(amu(2,j)))
      dcu(2,j) = dcu(2,j) + dkx*zt1
   10 continue
      return
      end
      subroutine EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,nxh
     1d)
c this subroutine solves 1-2/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,affp,wp0,nx,nxvh, output:ffe
c for isign =/ 0, input: dcu,ffe,isign,ci,nx,nxvh,nxhd, output: eyz,wf
c approximate flop count is: 25*nxc
c where nxc = nx/2 - 1
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ey(kx) = -ci*ci*g(kx)*dcuy(kx)*s(kx)
c ez(kx) = -ci*ci*g(kx)*dcuz(kx)*s(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/(kx**2+wp0*ci2*s(kx)**2))*s(kx),
c s(kx) = exp(-((kx*ax)**2+)/2), except for
c ey(kx=pi) = ez(kx=pi) = 0, and ey(kx=0) = ez(kx=0,) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ey(kx) = -ci*ci*g(kx)*dcuy(kx)
c ez(kx) = -ci*ci*g(kx)*dcuz(kx)
c dcu(i,j) = transverse part of complex derivative of current for
c fourier mode (j-1)
c eyz(1,j) = y component of complex transverse electric field
c eyz(2,j) = z component of complex transverse electric field
c all for fourier mode (j-1)
c aimag(ffe(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffe(j)) = potential green's function g
c for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprical of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*sum((affp/(kx**2*ci*ci)**2)*|dcu(kx)*s(kx)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx = system length in x direction
c nxvh = second dimension of field arrays, must be >= nxh
c nxhd = second dimension of form factor array, must be >= nxh
      implicit none
      integer isign, nx, nxvh, nxhd
      real ax, affp, wp0, ci, wf
      complex dcu, eyz, ffe
      dimension dcu(2,nxvh), eyz(2,nxvh)
      dimension ffe(nxhd)
      integer nxh, j
      real dnx, ci2, wpc, dkx, at1, at2
      complex zero
      double precision wp
      nxh = nx/2
      dnx = 6.28318530717959/real(nx)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
      if (isign.ne.0) go to 20
      wpc = wp0*ci2
c prepare form factor array
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at2 = exp(-.5*(dkx*ax)**2)
      ffe(j) = cmplx(affp*at2/(dkx*dkx+ wpc*at2*at2),at2)
   10 continue
      ffe(1) = cmplx(affp,1.0)
      return
c calculate smoothed transverse electric field and sum field energy
   20 if (isign.gt.0) go to 40
      wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 30 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*aimag(ffe(j))
      at2 = at2*at2
      eyz(1,j) = at1*dcu(1,j)
      eyz(2,j) = at1*dcu(2,j)
      wp = wp + at2*(dcu(1,j)*conjg(dcu(1,j)) + dcu(2,j)*conjg(dcu(2,j))
     1)
   30 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   40 wp = 0.0d0
c mode numbers 0 < kx < nx/2
      do 50 j = 2, nxh
      at2 = -ci2*real(ffe(j))
      at1 = at2*at2
      eyz(1,j) = at2*dcu(1,j)
      eyz(2,j) = at2*dcu(2,j)
      wp = wp + at1*(dcu(1,j)*conjg(dcu(1,j)) + dcu(2,j)*conjg(dcu(2,j))
     1)
   50 continue
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = real(nx)*wp/real(ffe(1))
      return
      end
      subroutine WPMXN1(qe,qi0,qbme,wpmax,wpmin,nx,nxe)
c calculates maximum and minimum plasma frequency.
c assumes guard cells have already been added
c qe = charge density for electrons and ions
c qi0 = charge density for ions, assumed constant
c qbme = charge/mass ratio for electrons
c wpmax/wpmin = maximum/minimum plasma frequency
c nx = system length in x direction
c nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      real qe, qi0, qbme, wpmax, wpmin
      integer nx, nxe
      dimension qe(nxe)
      integer j
      real at1
      wpmax = qbme*(qe(1) - qi0)
      wpmin = wpmax
      do 10 j = 1, nx
      at1 = qbme*(qe(j) - qi0)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
      return
      end
      subroutine ADDQEI1(qe,qi,nx,nxe)
c adds electron and ion densities
c assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c nx = system length in x/y direction
c nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      real qe, qi
      integer nx, nxe
      dimension qe(nxe), qi(nxe)
c local data
      integer j
      do 10 j = 1, nx
      qe(j) = qe(j) + qi(j)
   10 continue
      return
      end
      subroutine ADDQEI1X(qe,qi,qbme,qbmi,wpmax,wpmin,nx,nxe)
c adds electron and ion densities, and calculates maximum and minimum
c plasma frequency.  assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c qbme/qbmi = charge/mass ratio for electrons/ions
c wpmax/wpmin = maximum/minimum plasma frequency
c nx = system length in x direction
c nxe = first dimension of charge arrays, nxe must be >= nx
      implicit none
      real qe, qi, qbme, qbmi, wpmax, wpmin
      integer nx, nxe
      dimension qe(nxe), qi(nxe)
c local data
      integer j
      real at1
      wpmax = qbme*qe(1) + qbmi*qi(1)
      wpmin = wpmax
      do 10 j = 1, nx
      at1 = qbme*qe(j) + qbmi*qi(j)
      qe(j) = qe(j) + qi(j)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
      return
      end
      subroutine BADDEXT1(byz,omy,omz,nx,nxe)
c adds constant to magnetic field for 1-2/2d code
c byz = magnetic field
c omy/omz = magnetic field electron cyclotron frequency in y/z 
c nx = system length in x direction
c nxe = second dimension of magnetic field array, nxe must be >= nx
      implicit none
      real byz, omy, omz
      integer nx, nxe
      dimension byz(2,nxe)
      integer j
      do 10 j = 1, nx
      byz(1,j) = byz(1,j) + omy
      byz(2,j) = byz(2,j) + omz
   10 continue
      return
      end
      subroutine VRCOPY1(f,g,nx,ndim,nxv)
c this subroutine copies nx real vector array elements
c to another real array, nxv >= nx
      implicit none
      integer nx, ndim, nxv
      real f, g
      dimension f(ndim,nxv), g(ndim,nxv)
      integer i, j
      do 20 j = 1, nx
      do 10 i = 1, ndim
      g(i,j) = f(i,j)
   10 continue
   20 continue
      return
      end
      subroutine VCCOPY1(f,g,nx,ndim,nxv)
c this subroutine copies nx complex vector array elements
c to another complex array, nxv >= nx
      implicit none
      integer nx, ndim, nxv
      complex f, g
      dimension f(ndim,nxv), g(ndim,nxv)
      integer i, j
      do 20 j = 1, nx
      do 10 i = 1, ndim
      g(i,j) = f(i,j)
   10 continue
   20 continue
      return
      end
