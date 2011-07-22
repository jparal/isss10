c-----------------------------------------------------------------------
c 1d PIC library for solving field equations with dirichlet boundary
c conditions
c dfield2lib.f contains procedures to manage guard cells and solve
c              fields equations in fourier space for dirichlet
c              boundary conditions:
c LCGUARD1 replicates fields for 2 component vector array to replace
c          quadratic with linear interpolation at the edges.
c LBGUARD1 replicates fields for 3 component vector array to replace
c          quadratic with linear interpolation at the edges.
c LDGUARD1 replicates fields for scalar array to replace quadratic with
c          linear interpolation at the edges.
c LSCGUARD1 initialize field for 2 component non-periodic vector array,
c           quadratic interpolation.
c LSGUARD1 initialize field for non-periodic scalar array, quadratic
c          interpolation.
c LACGUARD1 add guard cells for 2 component non-periodic vector array,
c           to replace quadratic with linear interpolation at the edges.
c LAGUARD1 add guard cells for non-periodic scalar array, to replace
c          quadratic with linear interpolation at the edges.
c LSCGUARD1L initialize field for 2 component non-periodic vector array,
c            quadratic interpolation.
c LSGUARD1L initialize field for non-periodic scalar array, quadratic
c           interpolation.
c DBLSIN1A creates doubled array for 2 component vector data to enable
c          various sine/cosine transforms to be perfomed with ffts.
c DBLSIN1D creates doubled array for scalar vector data to enable
c          various sine/cosine transforms to be perfomed with ffts.
c HAFDBL1C extracts data from doubled array for 2 component vector data.
c HAFDBL1B extracts data from doubled array for 3 component vector data.
c HAFDBL1D extracts data from doubled array for scalar data.
c POISDX1 solve 1d poisson equation for electric force, potential, or
c         smoothing with dirichlet boundary conditions, using doubled
c         ffts.
c POISD1 solve 1d poisson equation for electric force, potential, or
c        smoothing with dirichlet boundary conditions, using sine or
c        cosine transforms.
c BPOISDX13 solve 1-1/2d vector poisson equation for magnetic force,
c           vector potential, or smoothing with dirichlet boundary
c           conditions, using doubled ffts.
c BPOISD13 solve 1-1/2d vector poisson equation for magnetic force,
c          vector potential, or smoothing with dirichlet boundary
c          conditions, using sine or cosine transforms.
c IBPOISDX13 solve 1-1/2d vector poisson equation for magnetic field
c            with dirichlet boundary conditions, using doubled ffts.
c IBPOISD13 solve 1-1/2d vector poisson equation for magnetic field
c           with dirichlet boundary conditions, using sine or cosine
c           transforms.
c MAXWELDX1 solve 1d maxwell equation for electric and magnetic fields
c           with dirichlet boundary conditions, using doubled ffts.
c MAXWELD1 solve 1d maxwell equation for electric and magnetic fields
c          with dirichlet boundary conditions, using sine or cosine
c          transforms.
c DMFIELDD1 copies scalar data from doubled fft format to sine format.
c CMFIELDD1 copies 2 component vector data from doubled fft format to
c           sine format.
c EMFIELDD1 combines and smooths 2d electric field in sine or format to
c           doubled fft format.
c BMFIELDD1 smooths 2d magnetic field in sine or format to doubled fft
c           format.
c AVPOTDX13 calculate 1-1/2d vector potential from magnetic field
c           with dirichlet boundary conditions, using doubled ffts.
c AVPOTD13 calculate 1-1/2d vector potential from magnetic field
c          with dirichlet boundary conditions, using sine transforms.
c AVRPOTDX13 calculate 1-1/2d radiative part of the vector potential
c            from current and magnetic field with dirichlet boundary
c            conditions, using doubled ffts.
c AVRPOTD13 calculate 1-1/2d radiative part of the vector potential
c            from current and magnetic field with dirichlet boundary
c            conditions, using sine transforms.
c GTSMODES1 extracts selected fourier sine components from potential
c           array.
c PTSMODES1 places selected fourier sine components into potential
c           array.
c GTVSMODES1 extracts selected fourier sine components from vector
c            potential array.
c PTVSMODES1 places selected fourier sine components into vector
c            potential array.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: july 21, 2010
c-----------------------------------------------------------------------
      subroutine LCGUARD1(byz,nx,nxe)
c this subroutine replicates field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx = system length in x direction
c nxe = first dimension of input array fxy, must be >= nx+3
      implicit none
      real byz
      integer nx, nxe
      dimension byz(2,nxe)
c local data
      integer nx3
      nx3 = nx + 3
      byz(1,1) = 2.*byz(1,2) - byz(1,3)
      byz(2,1) = 2.*byz(2,2) - byz(2,3)
      byz(1,nx3) = 2.*byz(1,nx+2) - byz(1,nx+1)
      byz(2,nx3) = 2.*byz(2,nx+2) - byz(2,nx+1)
      return
      end
      subroutine LBGUARD1(fxyz,nx,nxe)
c this subroutine replicates vector field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx = system length in x direction
c nxe = first dimension of input array bxy, must be >= nx+3
      implicit none
      real fxyz
      integer nx, nxe
      dimension fxyz(3,nxe)
c local data
      integer nx3
      nx3 = nx + 3
      fxyz(1,1) = 2.*fxyz(1,2) - fxyz(1,3)
      fxyz(2,1) = 2.*fxyz(2,2) - fxyz(2,3)
      fxyz(3,1) = 2.*fxyz(3,2) - fxyz(3,3)
      fxyz(1,nx3) = 2.*fxyz(1,nx+2) - fxyz(1,nx+1)
      fxyz(2,nx3) = 2.*fxyz(2,nx+2) - fxyz(2,nx+1)
      fxyz(3,nx3) = 2.*fxyz(3,nx+2) - fxyz(3,nx+1)
      return
      end
      subroutine LDGUARD1(fx,nx,nxe)
c this subroutine replicates scalar field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx = system length in x direction
c nxe = first dimension of input array q, must be >= nx+3
      implicit none
      real fx
      integer nx, nxe
      dimension fx(nxe)
c local data
      integer nx3
      nx3 = nx + 3
      fx(1) = 2.*fx(2) - fx(3)
      fx(nx3) = 2.*fx(nx+2) - fx(nx+1)
      return
      end
      subroutine LSCGUARD1(cu,yj0,zj0,nx,ngx,nxe)
c initialize extended non-periodic field
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real cu, yj0, zj0
      integer nx, ngx, nxe
      dimension cu(2,nxe)
      integer i, j, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 10 j = 2, nxg
      cu(1,j+ngx+1) = yj0
      cu(2,j+ngx+1) = zj0
   10 continue
      do 20 i = 1, 2
      cu(i,1) = 0.
      cu(i,2) = 0.
      cu(i,nx+2) = 0.
      cu(i,nx+3) = 0.
   20 continue
      cu(1,ngx+2) = .5*yj0
      cu(2,ngx+2) = .5*zj0
      cu(1,nx-ngx+2) = .5*yj0
      cu(2,nx-ngx+2) = .5*zj0
      return
      end
      subroutine LSGUARD1(q,qi0,nx,ngx,nxe)
c initialize extended non-periodic scalar field
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real q, qi0
      integer nx, ngx, nxe
      dimension q(nxe)
      integer j, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 10 j = 2, nxg
      q(j+ngx+1) = qi0
   10 continue
      q(1) = 0.
      q(2) = 0.
      q(nx+2) = 0.
      q(nx+3) = 0.
      q(ngx+2) = .5*qi0
      q(nx-ngx+2) = .5*qi0
      return
      end
      subroutine LACGUARD1(cu,nx,nxe)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx = system length in x direction
c nxe = second dimension of input array cu, must be >= nx+3
      implicit none
      real cu
      integer nx, nxe
      dimension cu(2,nxe)
c local data
      integer i, nx1
      nx1 = nx + 1
c add up guard cells
      do 10 i = 1, 2
      cu(i,2) = cu(i,2) + 2.*cu(i,1)
      cu(i,3) = cu(i,3) - cu(i,1)
      cu(i,nx+1) = cu(i,nx+1) - cu(i,nx+3)
      cu(i,nx+2) = cu(i,nx+2) + 2.*cu(i,nx+3)
      cu(i,1) = 0.
      cu(i,nx+3) = 0.
   10 continue
      return
      end
      subroutine LAGUARD1(q,nx,nxe)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c for scalar field
c nx = system length in x direction
c nxe = first dimension of input array q, must be >= nx+3
      implicit none
      real q
      integer nx, nxe
      dimension q(nxe)
c local data
      integer nx1
      nx1 = nx + 1
c add up guard cells
      q(2) = q(2) + 2.*q(1)
      q(3) = q(3) - q(1)
      q(nx+1) = q(nx+1) - q(nx+3)
      q(nx+2) = q(nx+2) + 2.*q(nx+3)
      q(1) = 0.
      q(nx+3) = 0.
      return
      end
      subroutine LSCGUARD1L(cu,yj0,zj0,nx,ngx,nxe)
c initialize extended non-periodic field
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real cu, yj0, zj0
      integer nx, ngx, nxe
      dimension cu(2,nxe)
      integer i, j, nxg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 10 j = 2, nxg
      cu(1,j+ngx) = yj0
      cu(2,j+ngx) = zj0
   10 continue
      do 20 i = 1, 2
      cu(i,1) = 0.
      cu(i,nx+1) = 0.
   20 continue
      cu(1,ngx+1) = .5*yj0
      cu(2,ngx+1) = .5*zj0
      cu(1,nx-ngx+1) = .5*yj0
      cu(2,nx-ngx+1) = .5*zj0
      return
      end
      subroutine LSGUARD1L(q,qi0,nx,ngx,nxe)
c initialize extended non-periodic scalar field
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real q, qi0
      integer nx, ngx, nxe
      dimension q(nxe)
      integer j, nxg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 10 j = 2, nxg
      q(j+ngx) = qi0
   10 continue
      q(1) = 0.
      q(nx+1) = 0.
      q(ngx+1) = .5*qi0
      q(nx-ngx+1) = .5*qi0
      return
      end
      subroutine DBLSIN1A(cu,cu2,nx,nxv,nx2v)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 1d real to complex fft.  y and z components are odd functions in x.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx = system length in x direction
c nxv = second dimension of input array cu, must be >= nx+1
c nx2v = second dimension of output array cu2, must be >= 2*nx
      implicit none
      real cu, cu2
      integer nx, nxv, nx2v
      dimension cu(2,nxv), cu2(2,nx2v)
c local data
      integer j, nxs
c copy to double array
      nxs = nx - 1
      do 10 j = 1, nxs
      cu2(1,j+1) = cu(1,j+1)
      cu2(2,j+1) = cu(2,j+1)
      cu2(1,nx+j+1) = -cu(1,nx-j+1)
      cu2(2,nx+j+1) = -cu(2,nx-j+1)
   10 continue
      cu2(1,1) = 0.
      cu2(2,1) = 0.
      cu2(1,nx+1) = 0.
      cu2(2,nx+1) = 0.
      return
      end
      subroutine DBLSIN1D(q,q2,nx,nxv,nx2v)
c this subroutine creates an odd array q2 from an array q, so that
c a 1d sine transform can be performed with a 1d real to complex fft.
c linear interpolation
c nx = system length in x direction
c nxv = first dimension of input array q, must be >= nx
c nx2v = first dimension of output array q2, must be >= 2*nx
      implicit none
      real q, q2
      integer nx, nxv, nx2v
      dimension q(nxv), q2(nx2v)
c local data
      integer j, nxs
c copy to double array
      nxs = nx - 1
      do 10 j = 1, nxs
      q2(j+1) = q(j+1)
      q2(nx+j+1) = -q(nx-j+1)
   10 continue
      q2(1) = 0.
      q2(nx+1) = 0.
      return
      end
      subroutine DBLCOS1D(q,q2,nx,nxv,nx2v)
c this subroutine creates an even array q2 from an array q, so that
c a 1d cosine transform can be performed with a 1d real to complex fft.
c linear interpolation
c nx = system length in x direction
c nxv = first dimension of input array q, must be >= nx
c nx2v = first dimension of output array q2, must be >= 2*nx
      implicit none
      real q, q2
      integer nx, nxv, nx2v
      dimension q(nxv), q2(nx2v)
c local data
      integer j, nxs
c copy to double array
      nxs = nx - 1
      do 10 j = 1, nxs
      q2(j+1) = q(j+1)
      q2(nx+j+1) = q(nx-j+1)
   10 continue
      q2(1) = q(1)
      q2(nx+1) = q(nx+1)
      return
      end
      subroutine HAFDBL1C(byz,byz2,nx,nxe,nx2v)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx = system length in x direction
c nxe = second dimension of output array fxy, must be >= nx+1
c nx2v = second dimension of input array fxy2, must be >= 2*nx
      implicit none
      real byz, byz2
      integer nx, nxe, nx2v
      dimension byz(2,nxe), byz2(2,nx2v)
c local data
      integer j, nx1
      nx1 = nx + 1
      do 10 j = 1, nx1
      byz(1,j) = byz2(1,j)
      byz(2,j) = byz2(2,j)
   10 continue
      return
      end
      subroutine HAFDBL1B(fxyz,fxyz2,nx,nxe,nx2v)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx = system length in x direction
c nxe = second dimension of output array bxy, must be >= nx+1
c nx2v = second dimension of input array bxy2, must be >= 2*nx
      implicit none
      real fxyz, fxyz2
      integer nx, nxe, nx2v
      dimension fxyz(3,nxe), fxyz2(3,nx2v)
c local data
      integer j, nx1
      nx1 = nx + 1
      do 10 j = 1, nx1
      fxyz(1,j) = fxyz2(1,j)
      fxyz(2,j) = fxyz2(2,j)
      fxyz(3,j) = fxyz2(3,j)
   10 continue
      return
      end
      subroutine HAFDBL1D(q,q2,nx,nxe,nx2v)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c nx = system length in x direction
c nxe = first dimension of output array q, must be >= nx+1
c nx2v = first dimension of input array q2, must be >= 2*nx
      implicit none
      real q, q2
      integer nx, nxe, nx2v
      dimension q(nxe), q2(nx2v)
c local data
      integer j, nx1
      nx1 = nx + 1
      do 10 j = 1, nx1
      q(j) = q2(j)
   10 continue
      return
      end
      subroutine POISDX1(q,fx,isign,ffd,ax,affp,we,nx)
c this subroutine solves 1d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin or cos transform
c for isign = 0, input: isign,ax,affp,nx, output: ffd
c for isign = +-1, input: q,ffd,isign,nx, output: fx,we
c for isign < 0, approximate flop count is: 9*nx
c for isign = 1, approximate flop count is: 5*nx
c for isign = 2, input: q,ffd,isign,nx, output: fx, flop count = nx
c if isign < 0, force/charge is calculated using the equations:
c fx(k) = -sqrt(-1)*k*g(k)*s(k)*q(k), where k = pi*j/nx, j=fourier mode,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-((k*ax)**2), except for
c fx(k=0) = fx(k=pi) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(k) = g(k)*q(k)*s(k)
c if isign = 2, smoothing is calculated using the equation:
c fx(k) = q(kx)*s(kx)
c cmplx(q(2*j-1),q(2*j)) = complex charge density for fourier mode (j-1)
c cmplx(fx(2*j-1),fx(2*j)) = complex force/charge for fourier mode (j-1)
c if isign = 0, form factor array is prepared
c ffd(2*j) = finite-size particle shape factor s
c for fourier mode (j-1)
c ffd(2*j-1) = potential green's function g for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*sum((affp/k**2)*|q(k)*s(k)|**2)
c nx = system length in x direction
      implicit none
      integer isign, nx
      real ax, affp, we
      real q, fx, ffd
      dimension q(2*nx), fx(2*nx), ffd(2*nx)
c local data
      integer j
      real dnx, dkx, at1, at2
      double precision wp
      dnx = 6.28318530717959/float(nx + nx)
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      ffd(2*j) = exp(-.5*(dkx*ax)**2)
      ffd(2*j-1) = affp*ffd(2*j)/(dkx*dkx)
   10 continue
      ffd(1) = affp
      ffd(2) = 1.0
      return
   20 if (isign.gt.0) go to 40
c calculate force/charge and sum field energy
      wp = 0.0d0
      do 30 j = 2, nx
      at1 = ffd(2*j-1)*ffd(2*j)
      at2 = dnx*float(j - 1)*at1
      fx(2*j-1) = at2*q(2*j)
      fx(2*j) = 0.0
      wp = wp + at1*q(2*j)**2
   30 continue
      fx(1) = 0.0
      fx(2) = 0.0
      we = float(nx)*wp
      return
c calculate potential and sum field energy
   40 if (isign.gt.1) go to 60
      wp = 0.0d0
      do 50 j = 2, nx
      at2 = ffd(2*j-1)
      at1 = at2*ffd(2*j)
      fx(2*j-1) = 0.0
      fx(2*j) = at2*q(2*j)
      wp = wp + at1*q(2*j)**2
   50 continue
      fx(1) = 0.0
      fx(2) = 0.0
      we = float(nx)*wp
      return
c calculate smoothing
   60 do 70 j = 2, nx
      at1 = ffd(2*j)
      fx(2*j-1) = 0.0
      fx(2*j) = at1*q(2*j)
   70 continue
      fx(1) = 0.0
      fx(2) = 0.0
      return
      end
      subroutine POISD1(q,fx,isign,ffd,ax,affp,we,nx,nxe,nx2v)
c this subroutine solves 1d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin or cos transforms
c for isign = 0, input: isign,ax,affp,nx,nx2v, output: ffd
c for isign = +-1, input: q,ffd,isign,nx,nx2v, output: fx,we
c for isign < 0, approximate flop count is: 9*nx
c for isign = 1, approximate flop count is: 5*nx
c for isign = 2, input: q,ffd,isign,nx,nx2v, output: fx, flop count = nx
c if isign < 0, force/charge is calculated using the equations:
c fx(k) = -sqrt(-1)*k*g(k)*s(k)*q(k), where k = pi*j/nx, j=fourier mode,
c g(k) = (affp/k**2)*s(k), and s(k) = exp(-((k*ax)**2), except for
c fx(k=0) = fx(k=pi) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(k) = g(k)*q(k)*s(k)
c if isign = 2, smoothing is calculated using the equation:
c fy(k) = q(k)*s(k)
c q(j) = charge density for fourier mode (j-1)
c fx(j) = force/charge for fourier mode (j-1)
c if isign = 0, form factor array is prepared
c ffd(2*j,k) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c ffd(2*j-1) = potential green's function g for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*sum((affp/k**2)*|q(k)*s(k)|**2)
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
c nx2v = first dimension of form factor array, must be >= 2*nx
      implicit none
      integer isign, nx, nxe, nx2v
      real ax, affp, we
      real q, fx, ffd
      dimension q(nxe), fx(nxe), ffd(nx2v)
c local data
      integer j, nx1
      real dnx, dkx, at1, at2
      double precision wp
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      ffd(2*j) = exp(-.5*((dkx*ax)**2))
      ffd(2*j-1) = affp*ffd(2*j)/(dkx*dkx)
   10 continue
      ffd(1) = affp
      ffd(2) = 1.0
      return
   20 if (isign.gt.0) go to 40
c calculate force/charge and sum field energy
      wp = 0.0d0
      do 30 j = 2, nx
      at1 = ffd(2*j-1)*ffd(2*j)
      at2 = dnx*float(j - 1)*at1
      fx(j) = -at2*q(j)
      wp = wp + at1*q(j)**2
   30 continue
      fx(1) = 0.0
      fx(nx+1) = 0.0
      we = float(nx)*wp
      return
c calculate potential and sum field energy
   40 if (isign.gt.1) go to 60
      wp = 0.0d0
      do 50 j = 2, nx
      at2 = ffd(2*j-1)
      at1 = at2*ffd(2*j)
      fx(j) = at2*q(j)
      wp = wp + at1*q(j)**2
   50 continue
      fx(1) = 0.0
      fx(nx+1) = 0.0
      we = float(nx)*wp
      return
c calculate smoothing
   60 do 70 j = 2, nx
      at1 = ffd(2*j)
      fx(j) = at1*q(j)
   70 continue
      fx(1) = 0.0
      fx(nx+1) = 0.0
      return
      end
      subroutine BPOISDX13(cu,byz,isign,ffd,ax,affp,ci,wm,nx,nxv,nxd)
c this subroutine solves 1-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin or cos transform
c for isign = 0, input: isign,ax,affp,nx,nxd, output: ffd
c for isign = -1, input: cu,ffd,isign,ci,nx,nxv,nxd, output: byz,wm
c approximate flop count is: 13*nxc
c for isign = 1, input: cu,ffd,isign,ci,nx,nxv,nxd, output: byz,wm
c approximate flop count is: 9*nxc
c for isign = 2, input: cu,ffd,isign,nx,nxv,nxd, output: byz
c approximate flop count is: 2*nxc
c where nxc = nx - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx)*s(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx)*s(kx),
c where kx = pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2)/2), except for
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
c aimag(ffd(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffd(j)) = potential green's function g
c for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2)
c this expression is valid only if the current is divergence-free
c nx = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nxd = dimension of form factor array, must be >= nx
      implicit none
      integer isign, nx, nxv, nxd
      real ax, affp, ci, wm
      complex cu, byz, ffd
      dimension cu(2,nxv), byz(2,nxv), ffd(nxd)
c local data
      integer j
      real dnx, ci2, dkx, at1, at2
      complex zero
      double precision wp
      dnx = 6.28318530717959/float(nx + nx)
      ci2 = ci*ci
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at2 = exp(-.5*(dkx*ax)**2)
      ffd(j) = cmplx(affp*at2/(dkx*dkx),at2)
   10 continue
      ffd(1) = cmplx(affp,1.0)
      return
   20 if (isign.gt.0) go to 40
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx
      do 30 j = 2, nx
      at1 = ci2*real(ffd(j))*aimag(ffd(j))
      at2 = dnx*float(j - 1)*at1
      byz(1,j) = cmplx(at2*aimag(cu(2,j)),0.0)
      byz(2,j) = cmplx(-at2*aimag(cu(1,j)),0.0)
      wp = wp + at1*(aimag(cu(1,j))**2 + aimag(cu(2,j))**2)
   30 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = float(nx)*wp
      return
c calculate vector potential and sum field energy
   40 if (isign.gt.1) go to 60
      wp = 0.0d0
c mode numbers 0 < kx < nx
      do 50 j = 2, nx
      at2 = ci2*real(ffd(j))
      at1 = at2*aimag(ffd(j))
      byz(1,j) = cmplx(0.,at2*aimag(cu(1,j)))
      byz(2,j) = cmplx(0.,at2*aimag(cu(2,j)))
      wp = wp + at1*(aimag(cu(1,j))**2 + aimag(cu(2,j))**2)
   50 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = float(nx)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx
   60 do 70 j = 2, nx
      at1 = aimag(ffd(j))
      byz(1,j) = cmplx(0.,at1*aimag(cu(1,j)))
      byz(2,j) = cmplx(0.,at1*aimag(cu(2,j)))
   70 continue
      byz(1,1) = zero
      byz(2,1) = zero
      return
      end
c-----------------------------------------------------------------------
      subroutine BPOISD13(cu,byz,isign,ffd,ax,affp,ci,wm,nx,nxe,nxv)
c this subroutine solves 1-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin or cos transforms
c for isign = 0, input: isign,ax,affp,nx,nxv, output: ffd
c for isign = -1, input: cu,ffd,isign,ci,nx,nxe,nxv, output: byz,wm
c approximate flop count is: 13*nxc
c for isign = 1, input: cu,ffd,isign,ci,nx,nxe,nxv, output: byz,wm
c approximate flop count is: 9*nxc
c for isign = 2, input: cu,ffd,isign,nx,nxe,nxv, output: byz
c approximate flop count is: 2*nxc
c where nxc = nx - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx)*s(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx)*s(kx),
c where kx = pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2)/2), except for
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
c aimag(ffd(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffd(j)) = potential green's function g
c for fourier mode (j-1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2)
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxv = dimension of form factor array, must be >= nx
      implicit none
      integer isign, nx, nxe, nxv
      real ax, affp, ci, wm
      real cu, byz
      complex ffd
      dimension cu(2,nxe), byz(2,nxe)
      dimension ffd(nxv)
c local data
      integer j
      real dnx, ci2, dkx, at1, at2
      double precision wp
      dnx = 6.28318530717959/float(nx + nx)
      ci2 = ci*ci
      if (isign.ne.0) go to 20
c prepare form factor array
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at2 = exp(-.5*(dkx*ax)**2)
      ffd(j) = cmplx(affp*at2/(dkx*dkx),at2)
   10 continue
      ffd(1) = cmplx(affp,1.0)
      return
   20 if (isign.gt.0) go to 40
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx
      do 30 j = 2, nx
      at1 = ci2*real(ffd(j))*aimag(ffd(j))
      at2 = dnx*float(j - 1)*at1
      byz(1,j) = -at2*cu(2,j)
      byz(2,j) = at2*cu(1,j)
      wp = wp + at1*(cu(1,j)**2 + cu(2,j)**2)
   30 continue
      byz(1,1) = 0.
      byz(2,1) = 0.
      byz(1,nx+1) = 0.
      byz(2,nx+1) = 0.
      wm = float(nx)*wp
      return
c calculate vector potential and sum field energy
   40 if (isign.gt.1) go to 60
      wp = 0.0d0
c mode numbers 0 < kx < nx
      do 50 j = 2, nx
      at2 = ci2*real(ffd(j))
      at1 = at2*aimag(ffd(j))
      byz(1,j) = at2*cu(1,j)
      byz(2,j) = at2*cu(2,j)
      wp = wp + at1*(cu(1,j)**2 + cu(2,j)**2)
   50 continue
      byz(1,1) = 0.
      byz(2,1) = 0.
      byz(1,nx+1) = 0.
      byz(2,nx+1) = 0.
      wm = float(nx)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx
   60 do 70 j = 2, nx
      at1 = aimag(ffd(j))
      byz(1,j) = at1*cu(1,j)
      byz(2,j) = at1*cu(2,j)
   70 continue
      byz(1,1) = 0.
      byz(2,1) = 0.
      byz(1,nx+1) = 0.
      byz(2,nx+1) = 0.
      return
      end
      subroutine IBPOISDX13(cu,byz,ffd,ci,wm,nx,nxv,nxd)
c this subroutine solves 1-1/2d poisson's equation in fourier space for
c magnetic field, with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin or cos transforms
c input: cu,ffd,ci,nx,nxv,nxd, output: byz,wm
c approximate flop count is: 13*nx
c the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx),
c where kx = pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2)/2), except for
c by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
c cu(i,j) = complex current density for fourier mode (j-1)
c byz(i,j) = i component of complex magnetic field
c all for fourier mode (j-1)
c aimag(ffd(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffd(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2), where
c affp = normalization constant = nx/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxv = first dimension of field arrays, must be >= nx
c nxd = dimension of form factor array, must be >= nx
      implicit none
      integer nx, nxv, nxd
      real ci, wm
      complex cu, byz, ffd
      dimension cu(2,nxv), byz(2,nxv), ffd(nxd)
c local data
      integer j
      real dnx, ci2, at1, at2
      complex zero
      double precision wp
      dnx = 6.28318530717959/float(nx + nx)
      ci2 = ci*ci
      zero = cmplx(0.0,0.0)
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx
      do 10 j = 2, nx
      at1 = ci2*real(ffd(j))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffd(j))
      byz(1,j) = cmplx(at2*aimag(cu(2,j)),0.0)
      byz(2,j) = cmplx(-at2*aimag(cu(1,j)),0.0)
      wp = wp + at1*(aimag(cu(1,j))**2 + aimag(cu(2,j))**2)
   10 continue
      byz(1,1) = zero
      byz(2,1) = zero
      wm = float(nx)*wp
      return
      end
      subroutine IBPOISD13(cu,byz,ffd,ci,wm,nx,nxe,nxv)
c this subroutine solves 1-1/2d poisson's equation in fourier space for
c magnetic field, with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin or cos transforms
c input: cu,ffd,ci,nx,nxe,nxv, output: byz,wm
c approximate flop count is: 13*n
c the magnetic field is calculated using the equations:
c by(kx) = -ci*ci*sqrt(-1)*g(kx)*kx*cuz(kx),
c bz(kx) = ci*ci*sqrt(-1)*g(kx)*kx*cuy(kx),
c where kx = pi*j/nx, and j = fourier mode numbers,
c g(kx) = (affp/kx**2)*s(kx),
c s(kx) = exp(-((kx*ax)**2)/2), except for
c by(kx=pi) = bz(kx=pi) = 0, and by(kx=0) = bz(kx=0) = 0.
c cu(i,j) = complex current density for fourier mode (j-1)
c byz(i,j) = i component of complex magnetic field
c all for fourier mode (j-1)
c aimag(ffd(j)) = finite-size particle shape factor s
c for fourier mode (j-1)
c real(ffd(j)) = potential green's function g
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*sum((affp/kx**2)*ci*ci*|cu(kx)*s(kx)|**2), where
c affp = normalization constant = nx/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxv = dimension of form factor array, must be >= nx
      implicit none
      integer nx, nxe, nxv
      real ci, wm
      real cu, byz
      complex ffd
      dimension cu(2,nxe), byz(2,nxe)
      dimension ffd(nxv)
c local data
      integer j
      real dnx, ci2, at1, at2
      double precision wp
      dnx = 6.28318530717959/float(nx + nx)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx
      do 10 j = 2, nx
      at1 = ci2*real(ffd(j))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffd(j))
      byz(1,j) = -at2*cu(2,j)
      byz(2,j) = at2*cu(1,j)
      wp = wp + at1*(cu(1,j)**2 + cu(2,j)**2)
   10 continue
      byz(1,1) = 0.
      byz(2,1) = 0.
      byz(1,nx+1) = 0.
      byz(2,nx+1) = 0.
      wm = float(nx)*wp
      return
      end
      subroutine MAXWELDX1(eyz,byz,cu,ffd,ci,dt,wf,wm,nx,nxv,nxd)
c this subroutine solves 1d maxwell's equation in fourier space for
c transverse electric and magnetic fields with dirichlet boundary
c conditions (zero potential).
c input: all, output: wf, wm, eyz, byz
c approximate flop count is: 35*nxc, where nxc = nx - 1
c the magnetic field is first updated half a step using the equations:
c by(kx) = by(kx) + .5*dt*sqrt(-1)*kx*ez(kx)
c bz(kx) = bz(kx) - .5*dt*sqrt(-1)*kx*ey(kx)
c the electric field is then updated a whole step using the equations:
c ey(kx) = ey(kx) - c2*dt*sqrt(-1)*kx*bz(kx) - affp*dt*cuy(kx)*s(kx)
c ez(kx) = ez(kx) + c2*dt*sqrt(-1)*kx*by(kx) - affp*dt*cuz(kx)*s(kx)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = pi*j/nx, c2 = 1./(ci*ci)
c and s(kx) = exp(-((kx*ax)**2)**2)
c cu(i,j) = complex current density
c eyz(i,j) = complex transverse electric field
c byz(i,j) = complex magnetic field
c for component i, all for fourier mode (j-1)
c real(ffd(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c aimag(ffd()) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*sum((1/affp)*|eyz(kx)|**2)
c magnetic field energy is also calculated, using
c wm = nx*sum((c2/affp)*|byz(kx)|**2)
c nx = system length in x direction
c nxv = first dimension of field arrays, must be >= nx
c nxd = dimension of form factor array, must be >= nx
      implicit none
      integer nx, nxv, nxd
      real ci, dt, wf, wm
      complex eyz, byz, cu, ffd
      dimension eyz(2,nxv), byz(2,nxv), cu(2,nxv), ffd(nxd)
c local data
      integer j
      real dnx, dth, c2, cdt, affp, adt, anorm, dkx, afdt
      real at5, at6, at8, at9
      complex zero
      double precision wp, ws
      if (ci.le.0.) return
      dnx = 6.28318530717959/float(nx + nx)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffd(1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffd(j))
      at8 = aimag(eyz(1,j))
      at9 = aimag(eyz(2,j))
c update magnetic field half time step, ky > 0
      at5 = real(byz(1,j)) - dth*(dkx*at9)
      at6 = real(byz(2,j)) + dth*(dkx*at8)
c update electric field whole time step
      at8 = at8 - cdt*(dkx*at6) - afdt*aimag(cu(1,j))
      at9 = at9 + cdt*(dkx*at5) - afdt*aimag(cu(2,j))
c update magnetic field half time step and store electric field
      at5 = at5 - dth*(dkx*at9)
      at6 = at6 + dth*(dkx*at8)
      ws = ws + anorm*(at8*at8 + at9*at9)
      wp = wp + anorm*(at5*at5 + at6*at6)
      eyz(1,j) = cmplx(0.0,at8)
      eyz(2,j) = cmplx(0.0,at9)
      byz(1,j) = cmplx(at5,0.0)
      byz(2,j) = cmplx(at6,0.0)
   10 continue
      byz(1,1) = zero
      byz(2,1) = zero
      eyz(1,1) = zero
      eyz(2,1) = zero
      wf = float(nx)*ws
      wm = float(nx)*c2*wp
      return
      end
      subroutine MAXWELD1(eyz,byz,cu,ffd,ci,dt,wf,wm,nx,nxe,nxv)
c this subroutine solves 1d maxwell's equation in fourier space for
c transverse electric and magnetic fields with  dirichlet boundary
c conditions (zero potential).
c input: all, output: wf, wm, eyz, byz
c approximate flop count is: 35*nxc, where nxc = nx - 1
c the magnetic field is first updated half a step using the equations:
c by(kx) = by(kx) + .5*dt*sqrt(-1)*kx*ez(kx)
c bz(kx) = bz(kx) - .5*dt*sqrt(-1)*kx*ey(kx)
c the electric field is then updated a whole step using the equations:
c ey(kx) = ey(kx) - c2*dt*sqrt(-1)*kx*bz(kx) - affp*dt*cuy(kx)*s(kx)
c ez(kx) = ez(kx) + c2*dt*sqrt(-1)*kx*by(kx) - affp*dt*cuz(kx)*s(kx)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = pi*j/nx, c2 = 1./(ci*ci)
c and s(kx) = exp(-((kx*ax)**2)
c cu(i,j) = complex current density
c exy(i,j) = complex transverse electric field
c bxy(i,j) = complex magnetic field
c for component i, all for fourier mode (j-1)
c real(ffd(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c aimag(ffd(j)) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*sum((1/affp)*|eyz(kx)|**2)
c magnetic field energy is also calculated, using
c wm = nx*sum((c2/affp)*|byz(kx)|**2)
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      integer nx, nxe, nxv
      real ci, dt, wf, wm
      real eyz, byz, cu
      complex ffd
      dimension eyz(2,nxe), byz(2,nxe), cu(2,nxe), ffd(nxv)
c local data
      integer j
      real dnx, dth, c2, cdt, affp, adt, anorm, dkx, afdt
      real at5, at6, at8, at9
      double precision wp, ws
      if (ci.le.0.) return
      dnx = 6.28318530717959/float(nx + nx)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffd(1))
      adt = affp*dt
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffd(j))
      at8 = eyz(1,j)
      at9 = eyz(2,j)
c update magnetic field half time step, ky > 0
      at5 = byz(1,j) + dth*(dkx*at9)
      at6 = byz(2,j) - dth*(dkx*at8)
c update electric field whole time step
      at8 = at8 + cdt*(dkx*at6) - afdt*cu(1,j)
      at9 = at9 - cdt*(dkx*at5) - afdt*cu(2,j)
c update magnetic field half time step and store electric field
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 - dth*(dkx*at8)
      ws = ws + anorm*(at8*at8 + at9*at9)
      wp = wp + anorm*(at5*at5 + at6*at6)
      byz(1,j) = at5
      byz(2,j) = at6
      eyz(1,j) = at8
      eyz(2,j) = at9
   10 continue
      byz(1,1) = 0.0
      byz(2,1) = 0.0
      byz(1,nx+1) = 0.0
      byz(2,nx+1) = 0.0
      eyz(1,1) = 0.0
      eyz(2,1) = 0.0
      eyz(1,nx+1) = 0.0
      eyz(2,nx+1) = 0.0
      wf = float(nx)*ws
      wm = float(nx)*c2*wp
      return
      end
      subroutine DMFIELDD1(q2,q,nx,nxv,nxe)
c this subroutine copies the charge density into a smaller array
c which would have been created by a fast sine transform in x
      implicit none
      integer nx, nxv, nxe
      complex q2
      real q
      dimension q2(nxv), q(nxe)
c local data
      integer j
      do 10 j = 1, nx+1
      q(j) = -aimag(q2(j))
   10 continue
      return
      end
      subroutine CMFIELDD1(cu2,cu,nx,nxv,nxe)
c this subroutine copies the current into a smaller array
c which would have been created by fast sine/cosine transforms in x
      implicit none
      integer nx, nxv, nxe
      complex cu2
      real cu
      dimension cu2(2,nxv), cu(2,nxe)
      integer j
      do 10 j = 1, nx+1
      cu(1,j) = -aimag(cu2(1,j))
      cu(2,j) = -aimag(cu2(2,j))
   10 continue
      return
      end
      subroutine EMFIELDD1(fxyz,fx,eyz,ffd,nx,nxv,nxe,nxd)
c this subroutine merges complex vector fields
c in fourier space with dirichlet boundary conditions (zero potential).
c includes additional smoothing
      implicit none
      integer nx, nxv, nxe, nxd
      complex fxyz, fx, ffd
      real eyz
      dimension fxyz(3,nxv), fx(nxv), eyz(2,nxe), ffd(nxd)
c local data
      integer j
      real at1
c add the fields
      do 10 j = 1, nx
      at1 = aimag(ffd(j))
      fxyz(1,j) = fx(j)
      fxyz(2,j) = cmplx(0.0,-eyz(1,j)*at1)
      fxyz(3,j) = cmplx(0.0,-eyz(2,j)*at1)
   10 continue
      return
      end
      subroutine BMFIELDD1(fxyz,eyz,ffd,nx,nxv,nxe,nxd)
c this subroutine copies complex vector fields
c in fourier space with dirichlet boundary conditions (zero potential).
c includes additional smoothing
      implicit none
      integer nx, nxv, nxe, nxd
      complex fxyz, ffd
      real eyz
      dimension fxyz(2,nxv), eyz(2,nxe), ffd(nxd)
c local data
      integer j
      real at1
c copy the magnetic fields
      do 10 j = 1, nx
      at1 = aimag(ffd(j))
      fxyz(1,j) = cmplx(eyz(1,j)*at1,0.0)
      fxyz(2,j) = cmplx(eyz(2,j)*at1,0.0)
   10 continue
      return
      end
      subroutine AVPOTDX13(byz,ayz,nx,nxv)
c this subroutine calculates 1-1/2d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential).
c input: byz, nx, nxv, output: ayz
c approximate flop count is: 5*nxc
c and nxc divides, where nxc = nx - 1,
c the vector potential is calculated using the equations:
c ay(kx) = -sqrt(-1)*bz(kx)/kx
c az(kx) = sqrt(-1)*by(kx)/kx,
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c ay(kx=pi) = az(kx=pi) = 0, and ay(kx=0) = az(kx=0) = 0.
c byz(i,j) = i component of complex magnetic field
c ayz(i,j) = i component of complex vector potential
c all for fourier mode (j-1)
c nx = system length in x direction
c nxv = first dimension of field arrays, must be >= nx
      implicit none
      integer nx, nxv
      complex byz, ayz
      dimension byz(2,nxv), ayz(2,nxv)
c local data
      integer j
      real dnx, dkx, at1, at4, at5, at6
      complex zero
      dnx = 6.28318530717959/float(nx + nx)
      zero = cmplx(0.0,0.0)
c calculate vector potential
c mode numbers 0 < kx < nx
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1.0/dkx
      at4 = real(byz(2,j))
      at5 = real(byz(1,j))
      at6 = at1*at5
      at5 = -at1*at4
      ayz(1,j) = cmplx(0.0,at5)
      ayz(2,j) = cmplx(0.0,at6)
   10 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      return
      end
      subroutine AVPOTD13(byz,ayz,nx,nxe)
c this subroutine calculates 1-1/2d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential).
c input: byz, nx,nxe, output: ayz
c approximate flop count is: 5*nxc
c and nxc divides, where nxc = nx - 1
c the vector potential is calculated using the equations:
c ay(kx) = -sqrt(-1)*bz(kx)/kx
c az(kx) = sqrt(-1)*by(kx)/kx,
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c ay(kx=pi) = az(kx=pi) = 0, and ay(kx=0) = az(kx=0) = 0.
c byz(i,j) = i component of complex magnetic field
c ayz(i,j) = i component of complex vector potential
c all for fourier mode (j-1)
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
      implicit none
      integer nx, nxe
      real byz, ayz
      dimension byz(2,nxe), ayz(2,nxe)
c local data
      integer j
      real dnx, dkx, at1, at4, at5
      dnx = 6.28318530717959/float(nx + nx)
c calculate vector potential
c mode numbers 0 < kx < nx
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1.0/dkx
      at4 = byz(2,j)
      at5 = byz(1,j)
      ayz(1,j) = at1*at4
      ayz(2,j) = -at1*at5
   10 continue
      ayz(1,1) = 0.0
      ayz(2,1) = 0.0
      ayz(1,nx+1) = 0.0
      ayz(2,nx+1) = 0.0
      return
      end
      subroutine AVRPOTDX13(ayz,byz,ffd,ci,nx,nxv,nxd)
c this subroutine solves 1-1/2d poisson's equation in fourier space for
c the radiative part of the vector potential
c with dirichlet boundary conditions (zero potential)
c input: all, output: ayz
c approximate flop count is: 13*nxc,
c and nxc divides, where nxc = nx - 1
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
c aimag(ffd()) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c nx = system length in x direction
c nxv = first dimension of field arrays, must be >= nx
c nxd = first dimension of form factor array, must be >= nx
      implicit none
      integer nx, nxv, nxd
      real ci
      complex ayz, byz, ffd
      dimension ayz(2,nxv), byz(2,nxv)
      dimension ffd(nxd)
c local data
      integer j
      real dnx, afc2, dkx, at1, at2, at4, at5
      complex zero
      if (ci.le.0.0) return
      dnx = 6.28318530717959/float(nx + nx)
      afc2 = real(ffd(1))*ci*ci
      zero = cmplx(0.0,0.0)
c calculate the radiative vector potential
c mode numbers 0 < kx < nx/2
      do 20 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1.0/(dkx*dkx)
      at2 = afc2*aimag(ffd(j))
c update radiative vector potential
      at4 = real(byz(2,j))
      at5 = real(byz(1,j))
      ayz(1,j) = cmplx(0.0,-at1*(dkx*at4 + at2*aimag(ayz(1,j))))
      ayz(2,j) = cmplx(0.0,at1*(dkx*at5 - at2*aimag(ayz(2,j))))
   20 continue
      ayz(1,1) = zero
      ayz(2,1) = zero
      return
      end
      subroutine AVRPOTD13(ayz,byz,ffd,ci,nx,nxe,nxd)
c this subroutine solves 1-1/2d poisson's equation in fourier space for
c the radiative part of the vector potential
c with dirichlet boundary conditions (zero potential)
c input: all, output: ayz
c approximate flop count is: 13*nxc,
c and nxc divides, where nxc = nx - 1
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
c aimag(ffd()) = finite-size particle shape factor s,
c s(kx) = exp(-((kx*ax)**2)/2)
c for fourier mode (j-1)
c ci = reciprical of velocity of light
c nx = system length in x direction
c nxe = first dimension of field arrays, must be >= nx+1
c nxd = first dimension of form factor array, must be >= nx
      implicit none
      integer nx, nxe, nxd
      real ci
      real ayz, byz
      complex ffd
      dimension ayz(2,nxe), byz(2,nxe)
      dimension ffd(nxd)
c local data
      integer j
      real dnx, afc2, dkx, at1, at2, at4, at5
      if (ci.le.0.0) return
      dnx = 6.28318530717959/float(nx + nx)
      afc2 = real(ffd(1))*ci*ci
c calculate the radiative vector potential
c mode numbers 0 < kx < nx/2
      do 20 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1.0/(dkx*dkx)
      at2 = afc2*aimag(ffd(j))
c update radiative vector potential
      at4 = byz(2,j)
      at5 = byz(1,j)
      ayz(1,j) = at1*(dkx*at4 - at2*ayz(1,j))
      ayz(2,j) = -at1*(dkx*at5 + at2*ayz(2,j))
   20 continue
      ayz(1,1) = 0.0
      ayz(2,1) = 0.0
      ayz(1,nx+1) = 0.0
      ayz(2,nx+1) = 0.0
      return
      end
      subroutine GTSMODES1(pot,pott,nx,it,modesx,nxe,nt,modesxd)
c this subroutine extracts lowest order modes from real array pot,
c which contains either sine or cosine coefficients
c and stores them into a location in a time history array pott
c modes stored: kx=(0,1,...,NX)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx+1
c nxe = dimension of input array pot, nxe >= nx+1
c nt = first dimension of output array pott, nt >= it
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, nxe, nt
      integer modesxd
      real pot, pott
      dimension pot(nxe), pott(nt,modesxd)
c local data
      integer jmax, j
      if (it.gt.nt) return
      if ((modesx.le.0).or.(modesx.gt.(nx+1))) return
      jmax = min0(modesx,nx+1)
c mode numbers 0 < kx < nx
      do 10 j = 1, jmax
      pott(it,j) = pot(j)
   10 continue
      return
      end
      subroutine PTSMODES1(pot,pott,nx,it,modesx,nxe,nt,modesxd)
c this subroutine extracts lowest order modes from a location in a time
c history array pott and stores them into real array pot
c which contains either sine or cosine coefficients
c modes stored: kx=(0,1,...,NX)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx+1
c nxe = dimension of input array pot, nxe >= nx+1
c nt = first dimension of output array pott, nt >= it
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, nxe, nt
      integer modesxd
      real pot, pott
      dimension pot(nxe), pott(nt,modesxd)
c local data
      integer jmax, j
      if (it.gt.nt) return
      if ((modesx.le.0).or.(modesx.gt.(nx+1))) return
      jmax = min0(modesx,nx+1)
c mode numbers 0 < kx < nx
      do 10 j = 1, jmax
      pot(j) = pott(it,j)
   10 continue
      do 20 j = jmax+1, nx+1
      pot(j) = 0.0
   20 continue
      return
      end
      subroutine GTVSMODES1(vpot,vpott,nx,it,modesx,ndim,nxv,nt,modesxd)
c this subroutine extracts lowest order modes from real vector array
c which contains either sine or cosine coefficients
c vpot and stores them into a location in a time history array vpott
c modes stored: kx=(0,1,...,NX)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxv >= nx+1
c nt = first dimension of output array vpott, nt >= it
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, ndim, nxv, nt
      integer modesxd
      real vpot, vpott
      dimension vpot(ndim,nxv), vpott(nt,ndim,modesxd)
c local data
      integer jmax, i, j
      if (it.gt.nt) return
      if ((modesx.le.0).or.(modesx.gt.(nx+1))) return
      jmax = min0(modesx,nx+1)
c mode numbers 0 < kx < nx
      do 20 j = 1, jmax
      do 10 i = 1, ndim
      vpott(it,i,j) = vpot(i,j)
   10 continue
   20 continue
      return
      end
      subroutine PTVSMODES1(vpot,vpott,nx,it,modesx,ndim,nxv,nt,modesxd)
c this subroutine extracts lowest order modes from a location in a time
c history array vpott and stores them into real vector array vpot
c which contains either sine or cosine coefficients
c modes stored: kx=(0,1,...,NX)
c nx = system length in x direction
c it = current time
c modesx = number of modes to store in x direction,
c where modesx <= nx+1
c ndim = number of field arrays, must be >= 1
c nxv = second dimension of input array vpot, nxv >= nx+1
c nt = first dimension of output array vpott, nt >= it
c modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, it, modesx, ndim, nxv, nt
      integer modesxd
      real vpot, vpott
      dimension vpot(ndim,nxv), vpott(nt,ndim,modesxd)
      integer jmax, i, j
      if (it.gt.nt) return
      if ((modesx.le.0).or.(modesx.gt.(nx+1))) return
      jmax = min0(modesx,nx+1)
c mode numbers 0 < kx < nx
      do 20 j = 1, jmax
      do 10 i = 1, ndim
      vpot(i,j) = vpott(it,i,j)
   10 continue
   20 continue
      do 40 j = jmax+1, nx+1
      do 30 i = 1, ndim
      vpot(i,j) = 0.0
   30 continue
   40 continue
      return
      end
