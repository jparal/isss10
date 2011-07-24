!-----------------------------------------------------------------------
!
      module dfield1d
!
! Fortran90 interface to 1d PIC Fortran77 library dfield1lib.f
! dfield1mod.f contains procedures to manage guard cells and solve 
!              fields in fourier space with dirichlet boundary
!              conditions:
!              defines module dfield2d
! laguard => ilaguard1 add guard cells for non-periodic scalar array
!            to replace quadratic with linear interpolation at the edges
!            calls LAGUARD1
! laguard => ilacguard1 add guard cells for non-periodic vector array
!            to replace quadratic with linear interpolation at the edges
!            calls LACGUARD1 or LAGUARD1
! lcguard => ilcguard1 copy guard cells for scalar or 2 and 3 component
!            non-periodic vector arrays to replace quadratic with linear
!            interpolation at the edges.
!            calls LDGUARD1, LCGUARD1, or LBGUARD1
! lcguard => ildguard1 copy guard cells for scalar arrays with various
!            interpolations.
!            calls LDGUARD1
! lsguard => ilscguard1 initialize field for scalar or 2 component
!            non-periodic vector array with various interpolations.
!            calls LSCGUARD1, LSGUARD1, LSCGUARD1L, or LSGUARD1L
! lsguard => ilsguard1 initialize field for scalar non-periodic array
!            with various interpolations.
!            calls LSGUARD1 or LSGUARD1L
! dblsin => idblsin1a double array for 1d vector data for dirichlet
!           boundary conditions.
!           calls DBLSIN1A
! dblsin => idblsin1d double array for 1d scalar data for dirichlet
!           boundary conditions.
!           calls DBLSIN1D
! dblcos => idblcos1a double array for 1d vector data for neumann
!           boundary conditions.
!           calls DBLCOS1A
! dblcos => idblcos1d double array for 1d scalar data for neumann
!           boundary conditions.
!           calls DBLCOS1D
! hafdbl => ihafdbl1c extract data from doubled array for 1d vector data
!           calls HAFDBL1D, HAFDBL1C, or HAFDBL1B
! hafdbl => ihafdbl1d extract data from doubled array for 1d scalar data
!           calls HAFDBL1D
! poisd_init => ipoisd1init initializes tables for field solvers, with
!               dirichlet boundary conditions.
!               calls POISDX1
! poisdx => ipoisdx1 solves poisson equation for electric force, potential,
!           or smoothing, with dirichlet boundary conditions, using ffts.
!           calls POISDX1
! poisd => ipoisd1 solves poisson equation for electric force, potential,
!          or smoothing, with dirichlet boundary conditions, using
!          sine/cosine transforms.
!          calls POISD1
! bpoisdx => jbpoisdx13 solves 1-1/2d vector poisson equation for
!            magnetic force with dirichlet boundary conditions, using
!            ffts.
!            calls BPOISDX13
! bpoisd => jbpoisd13 solves 1-1/2d vector poisson equation for magnetic
!           force with dirichlet boundary conditions, using sine/cosine
!           transforms.
!           calls BPOISD13
! spoisdx => ispoisdx13 provides vector smoothing with dirichlet
!            boundary conditions, using ffts.
!            calls BPOISDX13
! spoisd => ispoisd13 provides vector smoothing with dirichlet boundary
!           conditions, sine transform.
!           calls BPOISD13
! apoisdx => iapoisdx13 solves 1-1/2d vector poisson equation for vector
!            potential with dirichlet boundary conditions, using ffts.
!            calls BPOISDX13
! apoisd => iapoisd13 solves 1-1/2d vector poisson equation for vector
!           potential with dirichlet boundary conditions, using
!           sine/cosine transforms.
!           calls BPOISD13
! ibpoisdx => jibpoisdx13 solves vector poisson equation for magnetic
!             field with dirichlet boundary conditions, using ffts.
!             calls IBPOISDX13
! ibpoisd => jibpoisd13 solves vector poisson equation for magnetic
!            field with dirichlet boundary conditions, using sine/cosine
!            transforms.
!            calls IBPOISD13
! maxweldx => imaxweldx1 solves maxwell equation for electric and magnetic
!             fields with dirichlet boundary conditions, using ffts.
!             calls MAXWELDX1
! maxweld => imaxweld1 solves maxwell equation for electric and magnetic
!            fields with dirichlet boundary conditions, using
!            sine/cosine transforms.
!            calls MAXWELD1
! cmfieldd => icmfieldd1 copies vector data from doubled fft to
!             sine format.
!             calls CMFIELDD1
! cmfieldd => idmfieldd1 copies scalar data from doubled fft to
!             sine format.
!             calls DMFIELDD1
! amfieldd => iamfieldd1 copies vector data from doubled fft to
!             cosine format.
!             calls AMFIELDD1
! emfieldd => iemfieldd1 calculates electric forces from fields given by
!             maxwell and poisson equations with dirichlet boundary
!             conditions.
!             calls EMFIELDD1
! emfieldd => ibmfieldd1 calculates magnetic forces from fields given by
!             maxwell equations with dirichlet boundary conditions.
!             calls BMFIELDD1
! avpotdx => iavpotdx13 calculates vector potential from magnetic field
!            with dirichlet boundary conditions, using ffts.
!            calls AVPOTDX13
! avpotd => iavpotd13 calculates vector potential from magnetic field
!           with dirichlet boundary conditions, using sine/cosine
!           transforms.
!           calls AVPOTD13
! avrpotdx => iavrpotdx13 calculates radiative vector potential from
!             magnetic field and current, with dirichlet boundary
!             conditions, using ffts.
!             calling AVRPOTDX13
! avrpotd => iavrpotd13 calculates radiative vector potential from
!             magnetic field and current, with dirichlet boundary
!             conditions, using sine/cosine transforms.
!             calling AVRPOTD13
! gtsmodes => igtsmodes1 extracts selected sine components from
!             potential array.
!             calls GTSMODES1
! gtsmodes => igtvsmodes1 extracts selected sine components from vector
!             potential array.
!             calls GTVSMODES1
! ptsmodes => iptsmodes1 places selected sine components into potential
!             array.
!             calls PTSMODES1
! ptsmodes => iptvsmodes1 places selected sine components into vector
!             potential array.
!             calls PTVSMODES1
! dcuperpdx => idcuperpdx13 calculate transverse derivative of 1-2/2d
!              current density from momentum flux, with dirichlet
!              boundary conditions, using ffts.
!              calls DCUPERPDX13
! dcuperpd => idcuperpd13 calculate transverse derivative of 1-2/2d
!             current density from momentum flux, with dirichlet
!             boundary conditions, using sine/cosine transforms.
!             calls DCUPERPD13
! adcuperpdx => iadcuperpdx13 calculate transverse derivative of 1-2/2d
!               current density from momentum flux and acceleration
!               density, with dirichlet boundary conditions, using ffts.
!               calls ADCUPERPDX13
! adcuperpd => iadcuperpd13 calculate transverse derivative of 1-2/2d
!              current density from momentum flux and acceleration
!              density, with dirichlet boundary conditions, using
!              sine/cosine transforms.
!              calls ADCUPERPD13
! epoisd_init => iepoisd13init initializes tables for darwin field solver,
!                with dirichlet boundary conditions.
!                calls EPOISDX13
! epoisdx => jepoisdx13 solves 1-2/2d vector poisson equation for
!            transverse electric force with dirichlet boundary
!            conditions, using ffts.
!            calls EPOISDX13
! epoisd => jepoisd13 solves 1-2/2d vector poisson equation for
!           transverse electric force with dirichlet boundary
!           conditions, using sine/cosine transforms.
!           calls EPOISD13
! iepoisdx => jiepoisdx13 solves 1-2/2d vector poisson equation
!             for transverse electric field without smoothing, with
!             dirichlet boundary conditions, using ffts.
!             calls EPOISDX13
! iepoisd => jiepoisd13 solves 1-2/2d vector poisson equation
!            for transverse electric field without smoothing, with
!            dirichlet boundary conditions, using sine/cosine transforms
!            calls EPOISD13
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: july 19, 2011
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: LCGUARD1, LBGUARD1, LDGUARD1
      public :: LSCGUARD1, LSGUARD1, LSCGUARD1L, LSGUARD1L
      public :: LACGUARD1, LAGUARD1
      public :: laguard, lcguard, lsguard
      public :: dblsin, dblcos, hafdbl, poisd_init, poisdx, poisd
      public :: bpoisdx, bpoisd, spoisdx, spoisd, apoisdx, apoisd
      public :: ibpoisdx, ibpoisd, maxweldx, maxweld
      public :: cmfieldd, amfieldd, emfieldd
      public :: avpotdx, avpotd, avrpotdx, avrpotd, gtsmodes, ptsmodes
      public :: dcuperpdx, dcuperpd, adcuperpdx, adcuperpd
      public :: epoisd_init, epoisdx, epoisd, iepoisdx, iepoisd
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine LCGUARD1(byz,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(2,nxe) :: byz
         end subroutine
      end interface
      interface
         subroutine LBGUARD1(fxyz,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(3,nxe) :: fxyz
         end subroutine
      end interface
      interface
         subroutine LDGUARD1(fx,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(nxe) :: fx
         end subroutine
      end interface
      interface
         subroutine LSCGUARD1(cu,yj0,zj0,nx,ngx,nxe)
         implicit none
         real :: yj0, zj0
         integer :: nx, ngx, nxe
         real, dimension(2,nxe) :: cu
         end subroutine
      end interface
      interface
         subroutine LSGUARD1(q,qi0,nx,ngx,nxe)
         implicit none
         real :: qi0
         integer :: nx, ngx, nxe
         real, dimension(nxe) :: q
         end subroutine
      end interface
      interface
         subroutine LACGUARD1(cu,nx,nxe)
         implicit none
         integer :: nx, nxe
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine LAGUARD1(q,nx,nxe)
         implicit none
         integer :: nx, nxe
!        real, dimension(*) :: q
         real :: q
         end subroutine
      end interface
      interface
         subroutine LSCGUARD1L(cu,yj0,zj0,nx,ngx,nxe)
         implicit none
         real :: yj0, zj0
         integer :: nx, ngx, nxe
         real, dimension(2,nxe) :: cu
         end subroutine
      end interface
      interface
         subroutine LSGUARD1L(q,qi0,nx,ngx,nxe)
         implicit none
         real :: qi0
         integer :: nx, ngx, nxe
         real, dimension(nxe) :: q
         end subroutine
      end interface
      interface
         subroutine DBLSIN1A(cu,cu2,nx,nxv,nx2v)
         implicit none
         integer :: nx, nxv, nx2v
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(2,nx2v) :: cu2
         end subroutine
      end interface
      interface
         subroutine DBLSIN1D(q,q2,nx,nxv,nx2v)
         implicit none
         integer :: nx, nxv, nx2v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v) :: q2
         end subroutine
      end interface
      interface
         subroutine DBLCOS1A(cu,cu2,nx,nxv,nx2v)
         implicit none
         integer :: nx, nxv, nx2v
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(2,nx2v) :: cu2
         end subroutine
      end interface
      interface
         subroutine DBLCOS1D(q,q2,nx,nxv,nx2v)
         implicit none
         integer :: nx, nxv, nx2v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v) :: q2
         end subroutine
      end interface
      interface
         subroutine HAFDBL1C(byz,byz2,nx,nxe,nx2v)
         implicit none
         integer :: nx, nxe, nx2v
!        real, dimension(*) :: byz
         real :: byz
         real, dimension(2,nx2v) :: byz2
         end subroutine
      end interface
      interface
         subroutine HAFDBL1B(fxyz,fxyz2,nx,nxe,nx2v)
         implicit none
         integer :: nx, nxe, nx2v
!        real, dimension(*) :: fxyz
         real :: fxyz
         real, dimension(3,nx2v) :: fxyz2
         end subroutine
      end interface
      interface
         subroutine HAFDBL1D(q,q2,nx,nxe,nx2v)
         implicit none
         integer :: nx, nxe, nx2v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v) :: q2
         end subroutine
      end interface
      interface
         subroutine POISDX1(q,fx,isign,ffd,ax,affp,we,nx)
         implicit none
         integer :: isign, nx, nxv
         real :: ax, affp, we
         real, dimension(2*nx) :: q
         real, dimension(2*nx) :: fx
         complex, dimension(nx) :: ffd
         end subroutine
      end interface
      interface
         subroutine POISD1(q,fx,isign,ffd,ax,affp,we,nx,nxe,nx2v)
         implicit none
         integer :: isign, nx, nxe, nx2v
         real :: ax, affp, we
!        real, dimension(*) :: q, fx
         real :: q, fx
         complex, dimension(nx2v/2) :: ffd
         end subroutine
      end interface
      interface
         subroutine BPOISDX13(cu,byz,isign,ffd,ax,affp,ci,wm,nx,nxv,nxd)
         implicit none
         integer :: isign, nx, nxv, nxd
         real :: ax, affp, ci, wm
         real, dimension(2,2*nxv) :: cu, byz
         complex, dimension(nxd) :: ffd
         end subroutine
      end interface
      interface
         subroutine BPOISD13(cu,byz,isign,ffd,ax,affp,ci,wm,nx,nxe,nxv)
         implicit none
         integer :: isign, nx, nxe, nxv
         real :: ax, affp, ci, wm
!        real, dimension(*) :: cu, byz
         real :: cu, byz
         complex, dimension(nxv) :: ffd
         end subroutine
      end interface
      interface
         subroutine IBPOISDX13(cu,byz,ffd,ci,wm,nx,nxv,nxd)
         implicit none
         integer :: nx, nxv, nxd
         real :: ci, wm
         real, dimension(2,2*nxv) :: cu
         complex, dimension(2,nxv) :: byz
         complex, dimension(nxd) :: ffd
         end subroutine
      end interface
      interface
         subroutine IBPOISD13(cu,byz,ffd,ci,wm,nx,nxe,nxv)
         implicit none
         integer :: nx, nxe, nxv
         real :: ci, wm
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(2,nxe/2) :: byz
         complex, dimension(nxv) :: ffd
         end subroutine
      end interface
      interface
         subroutine MAXWELDX1(eyz,byz,cu,ffd,ci,dt,wf,wm,nx,nxv,nxd)
         implicit none
         integer :: nx, nxv, nxd
         real :: ci, dt, wf, wm
         complex, dimension(2,nxv) :: eyz, byz
         real, dimension(2,2*nxv) :: cu
         complex, dimension(nxd) :: ffd
         end subroutine
      end interface
      interface
         subroutine MAXWELD1(eyz,byz,cu,ffd,ci,dt,wf,wm,nx,nxe,nxv)
         implicit none
         integer :: nx, nxe, nxv
         real :: ci, dt, wf, wm
         complex, dimension(2,nxe/2) :: eyz, byz
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxv) :: ffd
         end subroutine
      end interface
      interface
         subroutine DMFIELDD1(q2,q,nx,nxv,nxe)
         implicit none
         integer :: nx, nxv, nxe
         real, dimension(2*nxv) :: q2
!        real, dimension(*) :: q
         real :: q
         end subroutine
      end interface
      interface
         subroutine CMFIELDD1(cu2,cu,nx,nxv,nxe)
         implicit none
         integer :: nx, nxv, nxe
         real, dimension(2,2*nxv) :: cu2
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine AMFIELDD1(amu2,amu,nx,nxv,nxe)
         implicit none
         integer :: nx, nxv, nxe
         real, dimension(2,2*nxv) :: amu2
!        real, dimension(*) :: amu
         real :: amu
         end subroutine
      end interface
      interface
         subroutine EMFIELDD1(fxyz,fx,eyz,ffd,nx,nxv,nxe,nxd)
         implicit none
         integer :: nx, nxv, nxe, nxd
         real, dimension(3,2*nxv) :: fxyz
         real, dimension(2*nxv) :: fx
         complex, dimension(2,nxe/2) :: eyz
         complex, dimension(nxd) :: ffd
         end subroutine
      end interface
      interface
         subroutine BMFIELDD1(fxyz,eyz,ffd,nx,nxv,nxe,nxd)
         implicit none
         integer :: nx, nxv, nxe, nxd
         real, dimension(2,2*nxv) :: fxyz
         complex, dimension(2,nxe/2) :: eyz
         complex, dimension(nxd) :: ffd
         end subroutine
      end interface
      interface
         subroutine AVPOTDX13(byz,ayz,nx,nxv)
         implicit none
         integer :: nx, nxv
         complex, dimension(2,nxv) :: byz
         real, dimension(2,2*nxv) :: ayz
         end subroutine
      end interface
      interface
         subroutine AVPOTD13(byz,ayz,nx,nxe)
         implicit none
         integer :: nx, nxe
         complex, dimension(2,nxe/2) :: byz
!        real, dimension(*) :: ayz
         real :: ayz
         end subroutine
      end interface
      interface
         subroutine AVRPOTDX13(ayz,byz,ffd,ci,nx,nxv,nxd)
         implicit none
         integer :: nx, nxv, nxd
         real :: ci
         complex, dimension(2,nxv) :: byz
         real, dimension(2,2*nxv) :: ayz
         complex, dimension(nxd) :: ffd
         end subroutine
      end interface
      interface
         subroutine AVRPOTD13(ayz,byz,ffd,ci,nx,nxe,nxd)
         implicit none
         integer :: nx, nxe, nxd
         real :: ci
         complex, dimension(2,nxe/2) :: byz
!        real, dimension(*) :: ayz
         real :: ayz
         complex, dimension(nxd) :: ffd
         end subroutine
      end interface
      interface
         subroutine GTSMODES1(pot,pott,nx,it,modesx,nxe,nt,modesxd)
         implicit none
         integer :: nx, it, modesx, nxe, nt, modesxd
!        real, dimension(*) :: pot
         real :: pot
         real, dimension(nt,modesxd) :: pott
         end subroutine
      end interface
      interface
         subroutine PTSMODES1(pot,pott,nx,it,modesx,nxe,nt,modesxd)
         implicit none
         integer :: nx, it, modesx, nxe, nt, modesxd
!        real, dimension(*) :: pot
         real :: pot
         real, dimension(nt,modesxd) :: pott
         end subroutine
      end interface
      interface
         subroutine GTVSMODES1(vpot,vpott,nx,it,modesx,ndim,nxv,nt,modes&
     &xd)
         implicit none
         integer :: nx, it, modesx, ndim, nxv, nt, modesxd
!        complex, dimension(*) :: vpot
         real :: vpot
         real, dimension(nt,modesxd) :: vpott
         end subroutine
      end interface
      interface
         subroutine PTVSMODES1(vpot,vpott,nx,it,modesx,ndim,nxv,nt,modes&
     &xd)
         implicit none
         integer :: nx, it, modesx, ndim, nxv, nt, modesxd
!        complex, dimension(*) :: vpot
         real :: vpot
         real, dimension(nt,modesxd) :: vpott
         end subroutine
      end interface
      interface
         subroutine LSCFGUARD1(cus,cu,q2m0,nx,nxe)
         implicit none
         integer :: nx, nxe
         real :: q2m0
         real, dimension(2,nxe) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine LSCFGUARD1L(cus,cu,q2m0,nx,nxe)
         implicit none
         integer :: nx, nxe
         real :: q2m0
         real, dimension(2,nxe) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine DCUPERPDX13(dcu,amu,nx,nxv)
         implicit none
         integer :: nx, nxv
         real, dimension(2,2*nxv) :: dcu, amu
         end subroutine
      end interface
      interface
         subroutine DCUPERPD13(dcu,amu,nx,nxe)
         implicit none
         integer :: nx, nxe
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
      interface
         subroutine ADCUPERPDX23(dcu,amu,nx,nxv)
         implicit none
         integer :: nx, nxv
         real, dimension(2,2*nxv) :: dcu, amu
         end subroutine
      end interface
      interface
        subroutine ADCUPERPD13(dcu,amu,nx,nxe)
         implicit none
         integer :: nx, nxe
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
      interface
         subroutine EPOISDX13(dcu,eyz,isign,fff,ax,affp,wp0,ci,wf,nx,nxv&
     &,nxd)
         implicit none
         integer :: isign, nx, nxv, nxd
         real :: ax, affp, wp0, ci, wf
         real, dimension(2,2*nxv) :: dcu, eyz
         complex, dimension(nxd) :: fff
         end subroutine
      end interface
!     interface
!        subroutine EPOISD13(dcu,eyz,isign,fff,ax,affp,wp0,ci,wf,nx,nxe,&
!    &nxv)
!        implicit none
!        integer :: isign, nx, nxe, nxv
!        real :: ax, affp, wp0, ci, wf
!        real, dimension(*) :: dcu, eyz
!        real :: dcu, eyz
!        complex, dimension(nxv) :: fff
!        end subroutine
!     end interface
!
! define generic interfaces to Fortran90 library
!      
      interface laguard
         module procedure ilaguard1
         module procedure ilacguard1
      end interface
!
      interface lcguard
         module procedure ilcguard1
         module procedure ildguard1
      end interface
!
      interface lsguard
         module procedure ilscguard1
         module procedure ilsguard1
         module procedure ilsfguard1
      end interface
!
      interface dblsin
         module procedure idblsin1a
         module procedure idblsin1d
      end interface
!
      interface dblcos
         module procedure idblcos1a
         module procedure idblcos1d
      end interface
!
      interface hafdbl
         module procedure ihafdbl1c
         module procedure ihafdbl1d
      end interface
!
      interface poisd_init
         module procedure ipoisd1init
      end interface
!
      interface poisdx
         module procedure ipoisdx1
      end interface
!
      interface poisd
         module procedure ipoisd1
      end interface
!
      interface bpoisdx
         module procedure jbpoisdx13
      end interface
!
      interface bpoisd
         module procedure jbpoisd13
      end interface
!
      interface spoisdx
         module procedure ispoisdx13
      end interface
!
      interface spoisd
         module procedure ispoisd13
      end interface
!
      interface apoisdx
         module procedure iapoisdx13
      end interface
!
      interface apoisd
         module procedure iapoisd13
      end interface
!
      interface ibpoisdx
         module procedure jibpoisdx13
      end interface
!
      interface ibpoisd
         module procedure jibpoisd13
      end interface
!
      interface maxweldx
         module procedure imaxweldx1
      end interface
!
      interface maxweld
         module procedure imaxweld1
      end interface
!
      interface cmfieldd
         module procedure icmfieldd1
         module procedure idmfieldd1
      end interface
!
      interface amfieldd
         module procedure iamfieldd1
      end interface
!
      interface emfieldd
         module procedure iemfieldd1
         module procedure ibmfieldd1
      end interface
!
      interface avpotdx
         module procedure iavpotdx13
      end interface
!
      interface avpotd
         module procedure iavpotd13
      end interface
!
      interface avrpotdx
         module procedure iavrpotdx13
      end interface
!
      interface avrpotd
         module procedure iavrpotd13
      end interface
!
      interface gtsmodes
         module procedure igtsmodes1
         module procedure igtvsmodes1
      end interface
!
      interface ptsmodes
         module procedure iptsmodes1
         module procedure iptvsmodes1
      end interface
!
      interface dcuperpdx
         module procedure idcuperpdx13
      end interface
!
      interface dcuperpd
         module procedure idcuperpd13
      end interface
!
      interface adcuperpdx
         module procedure iadcuperpdx13
      end interface
!
      interface adcuperpd
         module procedure iadcuperpd13
      end interface
!
       interface epoisd_init
         module procedure iepoisd13init
      end interface
!
      interface epoisdx
         module procedure jepoisdx13
      end interface
!
      interface epoisd
         module procedure jepoisd13
      end interface
!
      interface iepoisdx
         module procedure jiepoisdx13
      end interface
!
      interface iepoisd
         module procedure jiepoisd13
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ilacguard1(cu,nx,inorder)
! add guard cells for 1d vector data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: cu
! local data
         integer :: nxe, order
         nxe = size(cu,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order /= LINEAR) then
               call LACGUARD1(cu(1,2),nx-2,nxe)
            endif
         case (1)
            if (order /= LINEAR) then
               call LAGUARD1(cu(1,2),nx-2,nxe)
            endif
         end select
         end subroutine ilacguard1
!
         subroutine ilaguard1(q,nx,inorder)
! add guard cells for 1d scalar data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:), pointer :: q
! local data
         integer :: nxe, order
         nxe = size(q,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call LAGUARD1(q(2),nx-2,nxe)
         endif
         end subroutine ilaguard1
!
         subroutine ilcguard1(fxy,nx,inorder)
! copy guard cells for 1d vector data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: fxy
! local data
         integer :: nxe, order
         nxe = size(fxy,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (order /= LINEAR) then
               call LDGUARD1(fxy,nx,nxe)
            endif
         case (2)
            if (order /= LINEAR) then
               call LCGUARD1(fxy,nx,nxe)
            endif
         case (3)
            if (order /= LINEAR) then
               call LBGUARD1(fxy,nx,nxe)
            endif
         end select
         end subroutine ilcguard1
!
         subroutine ildguard1(q,nx,inorder)
! copy guard cells for 1d scalar data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:), pointer :: q
! local data
         integer :: nxe, order
         nxe = size(q,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call LDGUARD1(q,nx,nxe)
         endif
         end subroutine ildguard1
!
         subroutine ilscguard1(cu,yj0,zj0,nx,inorder)
! initialize 1d non-periodic vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: yj0, zj0
         real, dimension(:,:), pointer :: cu
! local data
         integer :: ngx = 1, nxe, order
         nxe = size(cu,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call LSCGUARD1L(cu,yj0,zj0,nx,ngx,nxe)
            else
               call LSCGUARD1(cu,yj0,zj0,nx,ngx,nxe)
            endif
         case (1)
            if (order==LINEAR) then
               call LSGUARD1L(cu,yj0,nx,ngx,nxe)
            else
               call LSGUARD1(cu,yj0,nx,ngx,nxe)
            endif
         end select
         end subroutine ilscguard1
!
         subroutine ilsguard1(q,qi0,nx,inorder)
! initialize non-uniform 1d non-periodic scalar field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qi0
         real, dimension(:), pointer :: q
! local data
         integer :: ngx = 1, nxe, order
         nxe = size(q,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call LSGUARD1L(q,qi0,nx,ngx,nxe)
         else
            call LSGUARD1(q,qi0,nx,ngx,nxe)
         endif
         end subroutine ilsguard1
!
         subroutine idblsin1a(cu,cu2,nx,inorder)
! double array in each dimension for 1d vector data
! for dirichlet boundary conditions
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: cu
         real, dimension(:,:), pointer :: cu2
! local data
         integer :: nxv, nx2v, order
         nxv = size(cu,2)
         nx2v = size(cu2,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call DBLSIN1A(cu(1,1),cu2,nx,nxv,nx2v)
            else
               call DBLSIN1A(cu(1,2),cu2,nx,nxv,nx2v)
            endif
         end select
         end subroutine idblsin1a
!
         subroutine idblsin1d(q,q2,nx,inorder)
! double array in each dimension for 1d scalar data
! for dirichlet boundary conditions
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:), pointer :: q
         real, dimension(:), pointer :: q2
! local data
         integer :: nxv, nx2v, order
         nxv = size(q,1)
         nx2v = size(q2,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call DBLSIN1D(q(1),q2,nx,nxv,nx2v)
         else
            call DBLSIN1D(q(2),q2,nx,nxv,nx2v)
         endif
         end subroutine idblsin1d
!
         subroutine idblcos1a(cu,cu2,nx,inorder)
! double array in each dimension for 1d vector data
! for neumann boundary conditions
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: cu
         real, dimension(:,:), pointer :: cu2
! local data
         integer :: nxv, nx2v, order
         nxv = size(cu,2)
         nx2v = size(cu2,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call DBLCOS1A(cu(1,1),cu2,nx,nxv,nx2v)
            else
               call DBLCOS1A(cu(1,2),cu2,nx,nxv,nx2v)
            endif
         end select
         end subroutine idblcos1a
!
         subroutine idblcos1d(q,q2,nx,inorder)
! double array in each dimension for 1d scalar data
! for neumann boundary conditions
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:), pointer :: q
         real, dimension(:), pointer :: q2
! local data
         integer :: nxv, nx2v, order
         nxv = size(q,1)
         nx2v = size(q2,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call DBLCOS1D(q(1),q2,nx,nxv,nx2v)
         else
            call DBLCOS1D(q(2),q2,nx,nxv,nx2v)
         endif
         end subroutine idblcos1d
!
         subroutine ihafdbl1c(fxy,fxy2,nx,inorder)
! copy from double to normal array in each dimension for 1d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: fxy
         real, dimension(:,:), pointer :: fxy2
! local data
         integer :: ndim, nxe, nx2v, order
         ndim = size(fxy,1)
         nxe = size(fxy,2)
         nx2v = size(fxy2,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (order==LINEAR) then
               call HAFDBL1D(fxy(1,1),fxy2,nx,nxe,nx2v)
            else
               call HAFDBL1D(fxy(1,2),fxy2,nx,nxe,nx2v)
            endif
         case (2)
            if (order==LINEAR) then
               call HAFDBL1C(fxy(1,1),fxy2,nx,nxe,nx2v)
            else
               call HAFDBL1C(fxy(1,2),fxy2,nx,nxe,nx2v)
            endif
         case (3)
            if (order==LINEAR) then
               call HAFDBL1B(fxy(1,1),fxy2,nx,nxe,nx2v)
            else
               call HAFDBL1B(fxy(1,2),fxy2,nx,nxe,nx2v)
            endif
         end select
         end subroutine ihafdbl1c
!
         subroutine ihafdbl1d(q,q2,nx,inorder)
! copy from double to normal array in each dimension for 1d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:), pointer :: q
         real, dimension(:), pointer :: q2
! local data
         integer :: nxe, nx2v, order
         nxe = size(q,1)
         nx2v = size(q2,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call HAFDBL1D(q(1),q2,nx,nxe,nx2v)
         else
            call HAFDBL1D(q(2),q2,nx,nxe,nx2v)
         endif
         end subroutine ihafdbl1d
!
         subroutine ipoisd1init(ffd,ax,affp,nx)
! initialize 1d electric field solver, conducting boundaries
         implicit none
         integer :: nx
         real :: ax, affp
         complex, dimension(:), pointer :: ffd
! local data
         integer :: isign = 0
         real :: we
         real, dimension(1) :: q
         real, dimension(1) :: fx
         call POISDX1(q,fx,isign,ffd,ax,affp,we,nx)
         end subroutine ipoisd1init
!
         subroutine ipoisdx1(q,fx,isign,ffd,we,nx)
! poisson solver for 1d electric field or potential,
! conducting boundaries, using fft
         implicit none
         integer :: isign, nx
         real :: we
         real, dimension(:), pointer :: q, fx
         complex, dimension(:), pointer :: ffd
! local data
         real :: ax, affp
         call POISDX1(q,fx,isign,ffd,ax,affp,we,nx)
         end subroutine ipoisdx1
!
         subroutine ipoisd1(q,fx,isign,ffd,we,nx,order)
! poisson solver for 1d electric field or potential,
! conducting boundaries, using fast sine/cosine transforms
         implicit none
         integer :: isign, nx
         integer, optional :: order
         real :: we
         real, dimension(:), pointer :: q, fx
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxe, nxv, nx2v, inorder
         real :: ax, affp
         nxe = size(q,1); nxv = size(ffd,1)
         nx2v = 2*nxv
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call POISD1(q(1),fx(1),isign,ffd,ax,affp,we,nx,nxe,nx2v)
         else
            call POISD1(q(2),fx(2),isign,ffd,ax,affp,we,nx,nxe,nx2v)
         endif
         end subroutine ipoisd1
!
         subroutine jbpoisdx13(cu,byz,ffd,ci,wm,nx)
! calculates static magnetic field for 1d vector field,
! conducting boundaries, using fft
         implicit none
         integer :: nx
         real :: ci, wm
         real, dimension(:,:), pointer :: cu, byz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: isign = -1, nxe, nxv
         real :: ax, affp
         nxe = size(cu,2)/2; nxv = size(ffd,1)
         select case(size(cu,1))
         case (2)
            call BPOISDX13(cu,byz,isign,ffd,ax,affp,ci,wm,nx,nxe,nxv)
         end select
         end subroutine jbpoisdx13
!
         subroutine jbpoisd13(cu,byz,ffd,ci,wm,nx,order)
! calculates static magnetic field for 1d vector field,
! conducting boundaries, using fast sine/cosine transforms
         implicit none
         integer :: nx
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:), pointer :: cu, byz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: isign = -1, nxe, nxv, inorder
         real :: ax, affp
         nxe = size(cu,2); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(cu,1))
         case (2)
            if (inorder==LINEAR) then
               call BPOISD13(cu(1,1),byz(1,1),isign,ffd,ax,affp,ci,wm,nx&
     &,nxe,nxv)
            else
               call BPOISD13(cu(1,2),byz(1,2),isign,ffd,ax,affp,ci,wm,nx&
     &,nxe,nxv)
            endif
         end select
         end subroutine jbpoisd13
!
         subroutine ispoisdx13(cu,byz,ffd,nx)
! smoother for non-periodic 1d vector field, using fft
         implicit none
         integer :: nx
         real, dimension(:,:), pointer :: cu, byz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: isign = 2, nxe, nxv
         real :: ax, affp, ci, wm
         nxe = size(cu,2)/2; nxv = size(ffd,1)
         select case(size(cu,1))
         case (2)
            call BPOISDX13(cu,byz,isign,ffd,ax,affp,ci,wm,nx,nxe,nxv)
         end select
         end subroutine ispoisdx13
!
         subroutine ispoisd13(cu,ayz,ffd,nx,order)
! smoother for non-periodic 1d vector field, using fast sine transform
         implicit none
         integer :: nx
         integer, optional :: order
         real, dimension(:,:), pointer :: cu, ayz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: isign = 2, nxe, nxv, inorder
         real :: ax, affp, ci, wm
         nxe = size(cu,2); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(cu,1))
         case (2)
            if (inorder==LINEAR) then
               call BPOISD13(cu(1,1),ayz(1,1),isign,ffd,ax,affp,ci,wm,nx&
     &,nxe,nxv)
            else
               call BPOISD13(cu(1,2),ayz(1,2),isign,ffd,ax,affp,ci,wm,nx&
     &,nxe,nxv)
            endif
         end select
         end subroutine ispoisd13
!
         subroutine iapoisdx13(cu,byz,ffd,ci,wm,nx)
! calculates static vector potential for 1d vector field,
! conducting boundaries, using fft
         implicit none
         integer :: nx
         real :: ci, wm
         real, dimension(:,:), pointer :: cu, byz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: isign = 1, nxe, nxv
         real :: ax, affp
         nxe = size(cu,2)/2; nxv = size(ffd,1)
         select case(size(cu,1))
         case (2)
            call BPOISDX13(cu,byz,isign,ffd,ax,affp,ci,wm,nx,nxe,nxv)
         end select
         end subroutine iapoisdx13
!
         subroutine iapoisd13(cu,ayz,ffd,ci,wm,nx,order)
! calculates static vector potential for 1d vector field,
! conducting boundaries
         implicit none
         integer :: nx
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:), pointer :: cu, ayz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: isign = 1, nxe, nxv, inorder
         real :: ax, affp
         nxe = size(cu,2); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(cu,1))
         case (2)
            if (inorder==LINEAR) then
               call BPOISD13(cu(1,1),ayz(1,1),isign,ffd,ax,affp,ci,wm,nx&
     &,nxe,nxv)
            else
               call BPOISD13(cu(1,2),ayz(1,2),isign,ffd,ax,affp,ci,wm,nx&
     &,nxe,nxv)
            endif
         end select
         end subroutine iapoisd13
!
         subroutine jibpoisdx13(cu,byz,ffd,ci,wm,nx)
! calculates static magnetic field for 1d vector field
! conducting boundaries, using fft
         implicit none
         integer :: nx
         real :: ci, wm
         real, dimension(:,:), pointer :: cu
         complex, dimension(:,:), pointer :: byz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxe, nxv
         nxe = size(cu,2)/2; nxv = size(ffd,1)
         call IBPOISDX13(cu,byz,ffd,ci,wm,nx,nxe,nxv)
         end subroutine jibpoisdx13
!
         subroutine jibpoisd13(cu,byz,ffd,ci,wm,nx,order)
! calculates static magnetic field for 1d vector field
! conducting boundaries
         implicit none
         integer :: nx
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:), pointer :: cu
         complex, dimension(:,:), pointer :: byz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxe, nxv, inorder
         nxe = size(cu,2); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call IBPOISD13(cu(1,1),byz,ffd,ci,wm,nx,nxe,nxv)
         else
            call IBPOISD13(cu(1,2),byz,ffd,ci,wm,nx,nxe,nxv)
         endif
         end subroutine jibpoisd13
!
         subroutine imaxweldx1(eyz,byz,cu,ffd,ci,dt,wf,wm,nx)
! calculates maxwell's equation for 1d vector field,
! conducting boundaries, using fft
         implicit none
         integer :: nx
         real :: ci, dt, wf, wm
         complex, dimension(:,:), pointer :: eyz, byz
         real, dimension(:,:), pointer :: cu
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxe, nxv
         nxe = size(cu,2)/2; nxv = size(ffd,1)
         call MAXWELDX1(eyz,byz,cu,ffd,ci,dt,wf,wm,nx,nxe,nxv)
         end subroutine imaxweldx1    
!
         subroutine imaxweld1(eyz,byz,cu,ffd,ci,dt,wf,wm,nx,order)
! calculates maxwell's equation for 1d vector field,
! conducting boundaries, using fast sine/cosine transforms
         implicit none
         integer :: nx
         integer, optional :: order
         real :: ci, dt, wf, wm
         complex, dimension(:,:), pointer :: eyz, byz
         real, dimension(:,:), pointer :: cu
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxe, nxv, inorder
         nxe = size(cu,2); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call MAXWELD1(eyz,byz,cu(1,1),ffd,ci,dt,wf,wm,nx,nxe,nxv)
         else
            call MAXWELD1(eyz,byz,cu(1,2),ffd,ci,dt,wf,wm,nx,nxe,nxv)
         endif
         end subroutine imaxweld1    
!
         subroutine icmfieldd1(cu2,cu,nx,order)
! copies from double to normal array in x dimension for 1d vector data
! for sine transform
         implicit none
         integer :: nx
         integer, optional :: order
         real, dimension(:,:), pointer :: cu2, cu
! local data
         integer :: nxv, nxe, inorder
         nxv = size(cu2,2)/2
         nxe = size(cu,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call CMFIELDD1(cu2,cu(1,1),nx,nxv,nxe)
         else
            call CMFIELDD1(cu2,cu(1,2),nx,nxv,nxe)
         endif
         end subroutine icmfieldd1
!
         subroutine idmfieldd1(q2,q,nx,order)
! copies from double to normal array in x dimension for 1d scalar data
! for sine transform
         implicit none
         integer :: nx
         integer, optional :: order
         real, dimension(:), pointer :: q2, q
! local data
         integer :: nxv, nxe, inorder
         nxv = size(q2,1)/2
         nxe = size(q,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call DMFIELDD1(q2,q(1),nx,nxv,nxe)
         else
            call DMFIELDD1(q2,q(2),nx,nxv,nxe)
         endif
         end subroutine idmfieldd1
!
         subroutine iamfieldd1(amu2,amu,nx,order)
! copies from double to normal array in x dimension for 1d vector data
! for cosine transform
         implicit none
         integer :: nx
         integer, optional :: order
         real, dimension(:,:), pointer :: amu2, amu
! local data
         integer :: nxv, nxe, inorder
         nxv = size(amu2,2)/2
         nxe = size(amu,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call AMFIELDD1(amu2,amu(1,1),nx,nxv,nxe)
         else
            call AMFIELDD1(amu2,amu(1,2),nx,nxv,nxe)
         endif
         end subroutine iamfieldd1
!
         subroutine iemfieldd1(fxyz,fx,eyz,ffd,nx)
! combines and smooths 1d vector fields, conducting boundaries
         implicit none
         integer :: nx
         real, dimension(:,:), pointer :: fxyz
         real, dimension(:), pointer :: fx
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxv, nxe, nxd
         nxv = size(fxyz,2)/2
         nxe = 2*size(eyz,2)
         nxd = size(ffd)
         call EMFIELDD1(fxyz,fx,eyz,ffd,nx,nxv,nxe,nxd)
         end subroutine iemfieldd1
!
         subroutine ibmfieldd1(fxyz,eyz,ffd,nx)
! combines and smooths 1d vector fields, conducting boundaries
         implicit none
         integer :: nx
         real, dimension(:,:), pointer :: fxyz
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxv, nxe, nxd
         nxv = size(fxyz,2)/2
         nxe = 2*size(eyz,2)
         nxd = size(ffd)
         call BMFIELDD1(fxyz,eyz,ffd,nx,nxv,nxe,nxd)
         end subroutine ibmfieldd1
!
         subroutine iavpotdx13(byz,ayz,nx)
! calculates 1d vector potential from magnetic field
! conducting boundaries, using fft
         implicit none
         integer :: nx
         complex, dimension(:,:), pointer :: byz
         real, dimension(:,:), pointer :: ayz
! local data
         integer :: nxe
         nxe = size(ayz,2)
         call AVPOTDX13(byz,ayz,nx,nxe)
         end subroutine iavpotdx13
!
         subroutine iavpotd13(byz,ayz,nx,order)
! calculates 1d vector potential from magnetic field
! conducting boundaries, using fast sine/cosine transforms
         implicit none
         integer :: nx
         integer, optional :: order
         complex, dimension(:,:), pointer :: byz
         real, dimension(:,:), pointer :: ayz
! local data
         integer :: nxe, inorder
         nxe = size(ayz,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call AVPOTD13(byz,ayz(1,1),nx,nxe)
         else
            call AVPOTD13(byz,ayz(1,2),nx,nxe)
         endif
         end subroutine iavpotd13
!
         subroutine iavrpotdx13(ayz,byz,ffd,ci,nx)
! calculates periodic 1d radiative vector potential from magnetic field
! and current, with conducting boundaries, using fft
         implicit none
         integer :: nx
         real :: ci
         complex, dimension(:,:), pointer :: byz
         real, dimension(:,:), pointer :: ayz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxe, nxv
         nxe = size(ayz,2); nxv = size(ffd,1)
         call AVRPOTDX13(ayz,byz,ffd,ci,nx,nxe,nxv)
         end subroutine iavrpotdx13
!
         subroutine iavrpotd13(ayz,byz,ffd,ci,nx,order)
! calculates periodic 1d radiative vector potential from magnetic field
! and current, with conducting boundaries,
! using fast sine/cosine transforms
         implicit none
         integer :: nx
         integer, optional :: order
         real :: ci
         complex, dimension(:,:), pointer :: byz
         real, dimension(:,:), pointer :: ayz
         complex, dimension(:), pointer :: ffd
! local data
         integer :: nxe, nxv, inorder
         nxe = size(ayz,2); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call AVRPOTD13(ayz(1,2),byz,ffd,ci,nx,nxe,nxv)
         else
            call AVRPOTD13(ayz(1,2),byz,ffd,ci,nx,nxe,nxv)
         endif
         end subroutine iavrpotd13
!
         subroutine igtsmodes1(pot,pots,nx,modesx,order)
! extracts lowest order modes from non-periodic 1d scalar field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:), pointer :: pot, pots
! local data
         integer :: nxe, it, nt, modesxd, inorder
         nxe = size(pot,1)
         it = 1; nt = 1
         modesxd = size(pots,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call GTSMODES1(pot(1),pots,nx,it,modesx,nxe,nt,modesxd)
         else
            call GTSMODES1(pot(2),pots,nx,it,modesx,nxe,nt,modesxd)
         endif
         end subroutine igtsmodes1
!
         subroutine iptsmodes1(pot,pots,nx,modesx,order)
! extracts lowest order modes to non-periodic 1d scalar field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:), pointer :: pot, pots
! local data
         integer :: nxe, it, nt, modesxd, inorder
         nxe = size(pot,1)
         it = 1; nt = 1
         modesxd = size(pots,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call PTSMODES1(pot(1),pots,nx,it,modesx,nxe,nt,modesxd)
         else
            call PTSMODES1(pot(2),pots,nx,it,modesx,nxe,nt,modesxd)
         endif
         end subroutine iptsmodes1
!
         subroutine igtvsmodes1(vpot,vpots,nx,modesx,order)
! extracts lowest order modes from non-periodic 1d vector field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:,:), pointer :: vpot, vpots
! local data
         integer :: ndim, nxv, it, nt, modesxd, inorder
         ndim = size(vpot,1); nxv = size(vpot,2)/2
         it = 1; nt = 1
         modesxd = size(vpots,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call GTVSMODES1(vpot(1,1),vpots,nx,it,modesx,ndim,nxv,nt,mod&
     &esxd)
         else
            call GTVSMODES1(vpot(1,2),vpots,nx,it,modesx,ndim,nxv,nt,mod&
     &esxd)
         endif
         end subroutine igtvsmodes1
!
         subroutine iptvsmodes1(vpot,vpots,nx,modesx,order)
! extracts lowest order modes to non-periodic 1d vector field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:,:), pointer :: vpot, vpots
! local data
         integer :: ndim, nxv, it, nt, modesxd, inorder
         ndim = size(vpot,1); nxv = size(vpot,2)
         it = 1; nt = 1
         modesxd = size(vpots,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call PTVSMODES1(vpot(1,1),vpots,nx,it,modesx,ndim,nxv,nt,mod&
     &esxd)
         else
            call PTVSMODES1(vpot(1,2),vpots,nx,it,modesx,ndim,nxv,nt,mod&
     &esxd)
         endif
         end subroutine iptvsmodes1
!
         subroutine ilsfguard1(cus,cu,q2m0,nx,inorder)
! initialize non-periodic 1d vector field with scaled field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: q2m0
         real, dimension(:,:), pointer :: cu, cus
! local data
         integer :: nxe, order
         nxe = size(cus,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cus,1))
         case (2)
            if (order==LINEAR) then
               call LSCFGUARD1L(cus,cu,q2m0,nx,nxe)
            else
               call LSCFGUARD1(cus,cu,q2m0,nx,nxe)
            endif
         end select
         end subroutine ilsfguard1
!
         subroutine idcuperpdx13(dcu,amu,nx)
! calculates the transverse part of non-periodic 1d vector field
! from momentum flux tensor, conducting boundaries, using fft
         implicit none
         integer :: nx
         real, dimension(:,:), pointer :: dcu, amu
! local data
         integer :: nxv
         nxv = size(dcu,2)/2
         select case(size(dcu,1))
         case (2)
            call DCUPERPDX13(dcu,amu,nx,nxv)
         end select
         end subroutine idcuperpdx13
!
         subroutine idcuperpd13(dcu,amu,nx,order)
! calculates the transverse part of non-periodic 1d vector field
! from momentum flux tensor, conducting boundaries, using fast
! sine/cosine transforms
         implicit none
         integer :: nx
         integer, optional :: order
         real, dimension(:,:), pointer :: dcu, amu
! local data
         integer :: nxe, inorder
         nxe = size(dcu,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(dcu,1))
         case (2)
            if (inorder==LINEAR) then
               call DCUPERPD13(dcu(1,1),amu(1,1),nx,nxe)
            else
               call DCUPERPD13(dcu(1,2),amu(1,2),nx,nxe)
            endif
         end select
         end subroutine idcuperpd13
!
         subroutine iadcuperpdx13(dcu,amu,nx)
! calculates the transverse part of non-periodic 1d vector field
! from acceleration vector and momentum flux tensor,
! conducting boundaries, using fft
         implicit none
         integer :: nx
         real, dimension(:,:), pointer :: dcu, amu
! local data
         integer :: nxv
         nxv = size(dcu,2)/2
         select case(size(dcu,1))
         case (2)
            call ADCUPERPDX13(dcu,amu,nx,nxv)
         end select
         end subroutine iadcuperpdx13
!
         subroutine iadcuperpd13(dcu,amu,nx,order)
! calculates the transverse part of non-periodic 1d vector field
! from acceleration vector and momentum flux tensor,
! conducting boundaries using fast sine/cosine transforms
         implicit none
         integer :: nx
         integer, optional :: order
         real, dimension(:,:), pointer :: dcu, amu
! local data
         integer :: nxe, inorder
         nxe = size(dcu,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(dcu,1))
         case (2)
            if (inorder==LINEAR) then
               call ADCUPERPD13(dcu(1,1),amu(1,1),nx,nxe)
            else
               call ADCUPERPD13(dcu(1,2),amu(1,2),nx,nxe)
            endif
         end select
         end subroutine iadcuperpd13
!
         subroutine iepoisd13init(fff,ax,affp,wp0,ci,nx)
! initialize 1d non-periodic transverse electric field solver
         implicit none
         integer :: nx
         real :: ax, affp, wp0, ci
         complex, dimension(:), pointer :: fff
! local data
         integer :: isign = 0, nxe = 1, nxv
         real :: wf
         real, dimension(1,1,1) :: dcu, eyz
         nxv = size(fff,1)
         call EPOISDX13(dcu,eyz,isign,fff,ax,affp,wp0,ci,wf,nx,nxe,nxv)
         end subroutine iepoisd13init
!
         subroutine jepoisdx13(dcu,eyz,fff,ci,wf,nx)
! calculates transverse electric field for periodic 1d vector field
! conducting boundaries, using fft
         implicit none
         integer :: nx
         real :: ci, wf
         real, dimension(:,:), pointer :: dcu, eyz
         complex, dimension(:), pointer :: fff
! local data
         integer :: isign = -1, nxe, nxv
         real :: ax, affp, wp0
         nxe = size(dcu,2)/2; nxv = size(fff,1)
         select case(size(dcu,1))
         case (2)
            call EPOISDX13(dcu,eyz,isign,fff,ax,affp,wp0,ci,wf,nx,nxe,nx&
     &v)
         end select
         end subroutine jepoisdx13
!
         subroutine jepoisd13(dcu,eyz,fff,ci,wf,nx,order)
! calculates transverse electric field for periodic 1d vector field
! conducting boundaries, using fast sine/cosine transforms
         implicit none
         integer :: nx
         integer, optional :: order
         real :: ci, wf
         real, dimension(:,:), pointer :: dcu, eyz
         complex, dimension(:), pointer :: fff
! local data
         integer :: isign = -1, nxe, nxv, inorder
         real :: ax, affp, wp0
         nxe = size(dcu,2); nxv = size(fff,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(dcu,1))
         case (2)
            if (inorder==LINEAR) then
               call EPOISD13(dcu(1,1),eyz(1,1),isign,fff,ax,affp,wp0,ci,&
     &wf,nx,nxe,nxv)
            else
               call EPOISD13(dcu(1,2),eyz(1,2),isign,fff,ax,affp,wp0,ci,&
     &wf,nx,nxe,nxv)
            endif
         end select
         end subroutine jepoisd13
!
         subroutine jiepoisdx13(dcu,eyz,fff,ci,wf,nx)
! calculates transverse electric field for periodic 1d vector field
! without smoothing, conducting boundaries, using fft
         implicit none
         integer :: nx
         real :: ci, wf
         real, dimension(:,:), pointer :: dcu, eyz
         complex, dimension(:), pointer :: fff
! local data
         integer :: isign = 1, nxe, nxv
         real :: ax, affp, wp0
         nxe = size(dcu,2)/2; nxv = size(fff,1)
         select case(size(dcu,1))
         case (2)
            call EPOISDX13(dcu,eyz,isign,fff,ax,affp,wp0,ci,wf,nx,nxe,nx&
     &v)
         end select
         end subroutine jiepoisdx13
!
         subroutine jiepoisd13(dcu,eyz,fff,ci,wf,nx,order)
! calculates transverse electric field for periodic 1d vector field
! without smoothing, conducting boundaries, using fast sine/cosine
! transforms
         implicit none
         integer :: nx
         integer, optional :: order
         real :: ci, wf
         real, dimension(:,:), pointer :: dcu
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: fff
! local data
         integer :: isign = 1, nxe, nxv, inorder
         real :: ax, affp, wp0
         nxe = size(dcu,2); nxv = size(fff,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(dcu,1))
         case (2)
            if (inorder==LINEAR) then
               call EPOISD13(dcu(1,1),eyz,isign,fff,ax,affp,wp0,ci,wf,nx&
     &,nxe,nxv)
            else
               call EPOISD13(dcu(1,2),eyz,isign,fff,ax,affp,wp0,ci,wf,nx&
     &,nxe,nxv)
            endif
         end select
         end subroutine jiepoisd13
!
      end module dfield1d
