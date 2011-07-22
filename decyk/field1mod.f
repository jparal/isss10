!-----------------------------------------------------------------------
!
      module field1d
!
! Fortran90 interface to 1d PIC Fortran77 library field1lib.f
! field1mod.f contains procedures to manage guard cells and solve fields
!             in fourier space:
! cguard => icguard1 copy guard cells for scalar or 2 and 3 component
!           vector arrays with various interpolations.
!           calls CGUARD1, BGUARD1, DGUARD1, CGUARD1L, BGUARD1L,
!           or DGUARD1L
! cguard => idguard1 copy guard cells for scalar arrays with various
!           interpolations.
!           calls DGUARD1, or DGUARD1L
! sguard => isguard1 initialize field for scalar array with various
!           interpolations.
!           calls SGUARD1 or SGUARD1L
! sguard => iscguard1 initialize field for 2 component vector array with
!           various interpolations.
!           calls SCGUARD1 or SCGUARD1L
! sguard => isfguard1 initialize 2 component field with scaled vector
!           array with various interpolations.
!           calls SCFGUARD1 or SCFGUARD1L
! aguard => iaguard1 add guard cells for scalar array with various
!           interpolations.
!           calls AGUARD1 or AGUARD1L
! aguard => iacguard1 add guard cells for 2 component vector array with
!           various interpolations.
!           calls ACGUARD1 or ACGUARD1L
! pois_init => ipois1init initializes tables for field solvers.
!              calls POISP1
! pois => ipois1 solves poisson equation for electric force, potential,
!         or smoothing.
!         calls POISP1
! bpois => jbpois13 solves vector poisson equation for magnetic force.
!          calls BPOIS13
! apois => iapois13 solves vector poisson equation for vector potential.
!          calls BPOIS13
! ibpois => iibpois13 solves vector poisson equation for magnetic field.
!           calls IBPOIS13
! maxwel => imaxwel1 solves maxwell equation for electric and magnetic
!           fields.
!           calls MAXWEL1
! emfield => iemfield1 calculates electric force from electric fields
!            given by maxwell and poisson equations.
!            calls EMFIELD1
! emfield => ibmfield1 calculates magnetic force from magnetic field
!            given by maxwell equation.
!            calls BMFIELD1
! emfieldr => iemfield1 calculates electric force from electric fields
!             given by maxwell and poisson equations for real arrays.
!             calls EMFIELDR1
! emfieldr => ibmfield1 calculates magnetic force from magnetic field
!             given by maxwell equation for real arrays.
!             calls BMFIELDR1
! avpot => iavpot13 calculates vector potential from magnetic field.
!          calls AVPOT13
! avrpot => iavrpot13 calculates radiative vector potential from
!           magnetic field and current
! gtmodes => igtmodes1 extracts selected fourier components from
!            potential array.
!            calls GTMODES1
! gtmodes => igtvmodes1 extracts selected fourier components from vector
!            potential array.
!            calls GTVMODES1
! ptmodes => iptmodes1 places selected fourier components into potential
!            array.
!            calls PTMODES1
! ptmodes => iptvmodes1 places selected fourier components into vector
!            potential array.
!            calls PTVMODES1
! dcuperp => idcuperp13 calculate transverse derivative of current
!            density from momentum flux.
!            calls DCUPERP13
! adcuperp => iadcuperp13 calculate transverse derivative of current
!             density from momentum flux and acceleration density.
!             calls ADCUPERP13
! epois_init => iepois13init initializes tables for darwin field solver.
!               calls EPOIS13
! epois => iepois13 solves vector poisson equation for transverse
!          electric force.
!          calls EPOIS13
! iepois => iiepois13 solves vector poisson equation for transverse
!           electric field.
!           calls EPOIS13
! wpmxn => iwpmxn1 calculates maximum and minimum plasma frequency
!          calls WPMXN1
! baddext => ibaddext1adds constant to magnetic field in real space for
!            1-2/2d code.
!            calls BADDEXT1
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: july 6, 2010
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: cguard, sguard, aguard, pois_init, pois, bpois
      public :: ibpois, maxwel, emfield, emfieldr
      public :: apois, avpot, avrpot, gtmodes, ptmodes
      public :: dcuperp, adcuperp, epois_init, epois, iepois
      public :: wpmxn, baddext
      public :: ivrcopy, ivccopy
!
! define interface to original Fortran77 procedures
      interface
         subroutine CGUARD1(byz,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(2,nxe) :: byz
         end subroutine
      end interface
      interface
         subroutine BGUARD1(fxyz,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(3,nxe) :: fxyz
         end subroutine
      end interface
      interface
         subroutine DGUARD1(fx,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(nxe) :: fx
         end subroutine
      end interface
      interface
         subroutine SCGUARD1(cu,yj0,zj0,nx,nxe)
         implicit none
         real :: yj0, zj0
         integer :: nx, nxe
         real, dimension(2,nxe) :: cu
         end subroutine
      end interface
      interface
         subroutine SGUARD1(q,qi0,nx,nxe)
         implicit none
         real :: qi0
         integer :: nx, nxe
         real, dimension(nxe) :: q
         end subroutine
      end interface
      interface
         subroutine ACGUARD1(cu,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(2,nxe) :: cu
         end subroutine
      end interface
      interface
         subroutine AGUARD1(q,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(nxe) :: q
         end subroutine
      end interface
      interface
         subroutine CGUARD1L(fx,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(nxe) :: fx
         end subroutine
      end interface
      interface
         subroutine BGUARD1L(byz,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(2,nxe) :: byz
         end subroutine
      end interface
      interface
         subroutine SCGUARD1L(cu,yj0,zj0,nx,nxe)
         implicit none
         real :: yj0, zj0
         integer :: nx, nxe
         real, dimension(2,nxe) :: cu
         end subroutine
      end interface
      interface
         subroutine SGUARD1L(q,qi0,nx,nxe)
         implicit none
         real :: qi0
         integer :: nx, nxe
         real, dimension(nxe) :: q
         end subroutine
      end interface
      interface
         subroutine ACGUARD1L(cu,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(2,nxe) :: cu
         end subroutine
      end interface
      interface
         subroutine AGUARD1L(q,nx,nxe)
         implicit none
         integer :: nx, nxe
         real, dimension(nxe) :: q
         end subroutine
      end interface
      interface
         subroutine POISP1(q,fx,isign,ffc,ax,affp,we,nx)
         implicit none
         integer :: isign, nx
         real :: ax, affp, we
         real :: q, fx
         complex, dimension(nx/2) :: ffc
         end subroutine
      end interface
      interface
         subroutine BPOIS13(cu,byz,isign,ffc,ax,affp,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer :: isign, nx, nxvh, nxhd
         real :: ax, affp, ci, wm
!        real, dimension(*) :: cu, byz
         real :: cu, byz
         complex, dimension(nxhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine IBPOIS13(cu,byz,ffc,ci,wm,nx,nxvh,nxhd)
         implicit none
         integer :: nx, nxvh, nxhd
         real :: ci, wm
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(2,nxvh) :: byz
         complex, dimension(nxhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine MAXWEL1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
         implicit none
         integer :: nx, nxvh, nxhd
         real :: ci, dt, wf, wm
         complex, dimension(3,nxvh) :: eyz, byz
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine EMFIELD1(fxyz,fx,eyz,ffc,nx,nxvh,nxhd)
         implicit none
         integer :: nx, nxvh, nxhd
!        real, dimension(*) :: fxyz, fx
         real :: fxyz, fx
         complex, dimension(2,nxvh) :: eyz
         complex, dimension(nxhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine BMFIELD1(fyz,eyz,ffc,nx,nxvh,nxhd)
         implicit none
         integer :: nx, nxvh, nxhd
!        real, dimension(*) :: fyz
         real :: fyz
         complex, dimension(2,nxvh) :: eyz
         complex, dimension(nxhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine EMFIELDR1(fxyz,fx,eyz,ffc,nx,nxe,nxd)
         implicit none
         integer :: nx, nxe, nxd
!        real, dimension(*) :: fxyz, fx
         real :: fxyz, fx
         complex, dimension(2,nxe/2) :: eyz
         complex, dimension(nxd) :: ffc
         end subroutine
      end interface
      interface
         subroutine BMFIELDR1(fyz,eyz,ffc,nx,nxe,nxd)
         implicit none
         integer :: nx, nxe, nxd
!        real, dimension(*) :: fyz
         real :: fyz
         complex, dimension(2,nxe/2) :: eyz
         complex, dimension(nxd) :: ffc
         end subroutine
      end interface
      interface
         subroutine AVPOT13(byz,ayz,nx,nxvh)
         implicit none
         integer :: nx, nxvh
         complex, dimension(2,nxvh) :: byz
!        real, dimension(*) :: ayz
         real :: ayz
         end subroutine
      end interface
      interface
         subroutine AVRPOT13(ayz,byz,ffc,ci,nx,nxvh,nxhd)
         implicit none
         integer :: nx, nxvh, nxhd
         real :: ci
         complex, dimension(2,nxvh) :: byz
!        real, dimension(*) :: ayz
         real :: ayz
         complex, dimension(nxhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine GTMODES1(pot,pott,nx,it,modesx,nxe,nt2,modesxd)
         implicit none
         integer :: nx, it, modesx, nxe, nt2, modesxd
!        real, dimension(*) :: pot
         real :: pot
         complex, dimension(nt2/2,modesxd) :: pott
         end subroutine
      end interface
      interface
         subroutine PTMODES1(pot,pott,nx,it,modesx,nxe,nt2,modesxd)
         implicit none
         integer :: nx, it, modesx, nxe, nt2, modesxd
!        real, dimension(*) :: pot
         real :: pot
         complex, dimension(nt2/2,modesxd) :: pott
         end subroutine
      end interface
      interface
         subroutine GTVMODES1(vpot,vpott,nx,it,modesx,ndim,nxvh,nt,modes&
     &xd)
         implicit none
         integer :: nx, it, modesx, ndim, nxvh, nt, modesxd
!        complex, dimension(*) :: vpot
         real :: vpot
         complex, dimension(nt,modesxd) :: vpott
         end subroutine
      end interface
      interface
         subroutine PTVMODES1(vpot,vpott,nx,it,modesx,ndim,nxvh,nt,modes&
     &xd)
         implicit none
         integer :: nx, it, modesx, ndim, nxvh, nt, modesxd
!        complex, dimension(*) :: vpot
         real :: vpot
         complex, dimension(nt,modesxd) :: vpott
         end subroutine
      end interface
      interface
         subroutine SCFGUARD1(cus,cu,q2m0,nx,nxe)
         implicit none
         real :: q2m0
         integer :: nx, nxe
         real, dimension(2,nxe) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine SCFGUARD1L(cus,cu,q2m0,nx,nxe)
         implicit none
         real :: q2m0
         integer :: nx, nxe
         real, dimension(2,nxe) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine DCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer :: nx, nxvh
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
      interface
         subroutine ADCUPERP13(dcu,amu,nx,nxvh)
         implicit none
         integer :: nx, nxvh
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
!     interface
!        subroutine EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,&
!    &nxhd)
!        implicit none
!        integer :: isign, nx, nxvh, nxhd
!        real :: ax, affp, wp0, ci, wf
!        real, dimension(*) :: dcu, eyz
!        real :: dcu, eyz
!        complex, dimension(nxhd) :: ffe
!        end subroutine
!     end interface
      interface
         subroutine BADDEXT1(byz,omy,omz,nx,nxe)
         implicit none
         integer :: nx, nxe
         real :: omy, omz
!        real, dimension(*,*,*) :: byz
         real :: byz
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface cguard
         module procedure icguard1
         module procedure idguard1
      end interface

      interface sguard
         module procedure isguard1
         module procedure iscguard1
         module procedure isfguard1
      end interface
      
      interface aguard
         module procedure iaguard1
         module procedure iacguard1
      end interface

       interface pois_init
         module procedure ipois1init
      end interface

      interface pois
         module procedure ipois1
      end interface
!
      interface bpois
         module procedure jbpois13
      end interface
!
      interface apois
         module procedure iapois13
      end interface
!
      interface ibpois
         module procedure iibpois13
      end interface
!
      interface maxwel
         module procedure imaxwel1
      end interface
!
      interface emfield
         module procedure iemfield1
         module procedure ibmfield1
      end interface
!
      interface emfieldr
         module procedure iemfieldr1
         module procedure ibmfieldr1
      end interface
!
      interface avpot
         module procedure iavpot13
      end interface
!
      interface avrpot
         module procedure iavrpot13
      end interface
!
      interface gtmodes
         module procedure igtmodes1
         module procedure igtvmodes1
      end interface
!
      interface ptmodes
         module procedure iptmodes1
         module procedure iptvmodes1
      end interface
!
      interface dcuperp
         module procedure idcuperp13
      end interface
!
      interface adcuperp
         module procedure iadcuperp13
      end interface
!
      interface epois_init
         module procedure iepois13init
      end interface
!
      interface epois
         module procedure iepois13
      end interface
!
      interface iepois
         module procedure iiepois13
      end interface
!
      interface wpmxn
         module procedure iwpmxn1
      end interface
!
      interface baddext
         module procedure ibaddext1
      end interface
!
      interface ivrcopy
         module procedure ivrcopy1
      end interface
!
      interface ivccopy
         module procedure ivccopy1
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine icguard1(fxy,nx,inorder)
! copy guard cells for periodic 1d vector data, for 1:2 component vectors
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: fxy
         integer :: nxe, order
         nxe = size(fxy,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (order==LINEAR) then
               call DGUARD1L(fxy,nx,nxe)
            else
               call DGUARD1(fxy,nx,nxe)
            endif
         case (2)
            if (order==LINEAR) then
               call CGUARD1L(fxy,nx,nxe)
            else
               call CGUARD1(fxy,nx,nxe)
            endif
         case (3)
            if (order==LINEAR) then
               call BGUARD1L(fxy,nx,nxe)
            else
               call BGUARD1(fxy,nx,nxe)
            endif
         end select
         end subroutine icguard1
!
         subroutine idguard1(fx,nx,inorder)
! copy guard cells for periodic 1d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:), pointer :: fx
         integer :: nxe, order
         nxe = size(fx)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call DGUARD1L(fx,nx,nxe)
         else
            call DGUARD1(fx,nx,nxe)
         endif
         end subroutine idguard1
!
         subroutine iscguard1(cu,yj0,zj0,nx,inorder)
! initialize periodic 1d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: yj0, zj0
         real, dimension(:,:), pointer :: cu
         integer :: nxe, order
         nxe = size(cu,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call SCGUARD1L(cu,yj0,zj0,nx,nxe)
         else
            call SCGUARD1(cu,yj0,zj0,nx,nxe)
         endif
         end subroutine iscguard1
!
         subroutine isguard1(q,qi0,nx,inorder)
! initialize periodic 1d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qi0
         real, dimension(:), pointer :: q
         integer :: nxe, order
         nxe = size(q)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call SGUARD1L(q,qi0,nx,nxe)
         else
            call SGUARD1(q,qi0,nx,nxe)
         endif
         end subroutine isguard1
!
         subroutine iacguard1(cu,nx,inorder)
! add guard cells for periodic 1d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: cu
         integer :: nxe, order
         nxe = size(cu,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call ACGUARD1L(cu,nx,nxe)
         else
            call ACGUARD1(cu,nx,nxe)
         endif
         end subroutine iacguard1
!
         subroutine iaguard1(q,nx,inorder)
! add guard cells for periodic 1d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:), pointer :: q
         integer :: nxe, order
         nxe = size(q)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call AGUARD1L(q,nx,nxe)
         else
            call AGUARD1(q,nx,nxe)
         endif
         end subroutine iaguard1
!
         subroutine ipois1init(ffc,ax,affp,nx)
! initialize 1d periodic electric field or poisson solver
         implicit none
         integer :: nx
         real :: ax, affp
         complex, dimension(:), pointer :: ffc
         integer :: isign = 0
         real :: we
         real, dimension(1) :: q, fx
         call POISP1(q(1),fx(1),isign,ffc,ax,affp,we,nx)
         end subroutine ipois1init
!
         subroutine ipois1(q,fx,isign,ffc,we,nx,inorder)
! poisson solver for periodic 1d electric field or potential
         implicit none
         integer :: isign, nx
         integer, optional :: inorder
         real :: we
         real, dimension(:), pointer :: q, fx
         complex, dimension(:), pointer :: ffc
         integer :: order
         real :: ax, affp
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call POISP1(q(1),fx(1),isign,ffc,ax,affp,we,nx)
         else
            call POISP1(q(2),fx(2),isign,ffc,ax,affp,we,nx)
         endif
         end subroutine ipois1
!
         subroutine jbpois13(cu,byz,ffc,ci,wm,nx,inorder)
! calculates static magnetic field for periodic 1d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: ci, wm
         real, dimension(:,:), pointer :: cu, byz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: isign = -1, nxvh, nxhd, order
         real :: ax, affp
         nxvh = size(cu,2)/2
         nxhd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call BPOIS13(cu(1,1),byz(1,1),isign,ffc,ax,affp,ci,wm,nx,&
     &nxvh,nxhd)
            else
               call BPOIS13(cu(1,2),byz(1,2),isign,ffc,ax,affp,ci,wm,nx,&
     &nxvh,nxhd)
            endif
         end select
         end subroutine jbpois13
!
         subroutine iapois13(cu,ayz,ffc,ci,wm,nx,inorder)
! calculates static vector potential for periodic 1d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: ci, wm
         real, dimension(:,:), pointer :: cu, ayz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: isign = 1, nxvh, nxhd, order
         real :: ax, affp
         nxvh = size(cu,2)/2
         nxhd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call BPOIS13(cu(1,1),ayz(1,1),isign,ffc,ax,affp,ci,wm,nx,&
     &nxvh,nxhd)
            else
               call BPOIS13(cu(1,2),ayz(1,2),isign,ffc,ax,affp,ci,wm,nx,&
     &nxvh,nxhd)
            endif
         end select
         end subroutine iapois13
!
         subroutine iibpois13(cu,byz,ffc,ci,wm,nx,inorder)
! calculates static magnetic field for periodic 1d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: ci, wm
         real, dimension(:,:), pointer :: cu
         complex, dimension(:,:), pointer :: byz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: nxvh, nxhd, order
         nxvh = size(cu,2)/2
         nxhd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call IBPOIS13(cu(1,1),byz,ffc,ci,wm,nx,nxvh,nxhd)
         else
            call IBPOIS13(cu(1,2),byz,ffc,ci,wm,nx,nxvh,nxhd)
         endif
         end subroutine iibpois13
!
         subroutine imaxwel1(eyz,byz,cu,ffc,ci,dt,wf,wm,nx,inorder)
! calculates maxwell's equation for periodic 1d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: ci, dt, wf, wm
         complex, dimension(:,:), pointer :: eyz, byz
         real, dimension(:,:), pointer :: cu
         complex, dimension(:), pointer :: ffc
! local data
         integer :: nxvh, nxhd, order
         nxvh = size(cu,2)/2
         nxhd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call MAXWEL1(eyz,byz,cu(1,1),ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
         else
            call MAXWEL1(eyz,byz,cu(1,2),ffc,ci,dt,wf,wm,nx,nxvh,nxhd)
         endif
         end subroutine imaxwel1
!
         subroutine iemfield1(fxyz,fx,eyz,ffc,nx,inorder)
! combines and smooths periodic 1d vector fields
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: fxyz
         real, dimension(:), pointer :: fx
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: nxvh, nxhd, order
         nxvh = size(fxyz,2)/2
         nxhd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call EMFIELD1(fxyz(1,1),fx(1),eyz,ffc,nx,nxvh,nxhd)
         else
            call EMFIELD1(fxyz(1,2),fx(2),eyz,ffc,nx,nxvh,nxhd)
         endif
         end subroutine iemfield1
!
         subroutine ibmfield1(fyz,eyz,ffc,nx,inorder)
! copies and smooths periodic 1d vector fields
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: fyz
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: nxvh, nxhd, order
         nxvh = size(fyz,2)/2
         nxhd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call BMFIELD1(fyz(1,1),eyz,ffc,nx,nxvh,nxhd)
         else
            call BMFIELD1(fyz(1,2),eyz,ffc,nx,nxvh,nxhd)
         endif
         end subroutine ibmfield1
!
         subroutine iemfieldr1(fxyz,fx,eyz,ffc,nx,inorder)
! combines and smooths 1d vector fields
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: fxyz
         real, dimension(:), pointer :: fx
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: nxe, nxd, order
         nxe = size(fxyz,2)
         nxd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call EMFIELDR1(fxyz(1,1),fx(1),eyz,ffc,nx,nxe,nxd)
         else
            call EMFIELDR1(fxyz(1,2),fx(2),eyz,ffc,nx,nxe,nxd)
         endif
         end subroutine iemfieldr1
!
         subroutine ibmfieldr1(fyz,eyz,ffc,nx,inorder)
! copies and smooths 1d vector fields
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: fyz
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: nxe, nxd, order
         nxe = size(fyz,2)
         nxd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call BMFIELDR1(fyz(1,1),eyz,ffc,nx,nxe,nxd)
         else
            call BMFIELDR1(fyz(1,2),eyz,ffc,nx,nxe,nxd)
         endif
         end subroutine ibmfieldr1
!
         subroutine iavpot13(byz,ayz,nx,inorder)
! calculates periodic 1d vector potential from magnetic field
         implicit none
         integer :: nx
         integer, optional :: inorder
         complex, dimension(:,:), pointer :: byz
         real, dimension(:,:), pointer :: ayz
! local data
         integer :: nxvh, order
         nxvh = size(byz,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call AVPOT13(byz,ayz(1,1),nx,nxvh)
         else
            call AVPOT13(byz,ayz(1,2),nx,nxvh)
         endif
         end subroutine iavpot13
!
         subroutine iavrpot13(ayz,byz,ffc,ci,nx,inorder)
! calculates periodic 1d radiative vector potential from magnetic field
! and current
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: ci
         complex, dimension(:,:), pointer :: byz
         real, dimension(:,:), pointer :: ayz
         complex, dimension(:), pointer :: ffc
! local data
         integer :: nxvh, nxhd, order
         nxvh = size(ayz,2)/2
         nxhd = size(ffc,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call AVRPOT13(ayz(1,1),byz,ffc,ci,nx,nxvh,nxhd)
         else
            call AVRPOT13(ayz(1,2),byz,ffc,ci,nx,nxvh,nxhd)
         endif
         end subroutine iavrpot13
!
         subroutine igtmodes1(pot,pott,nx,modesx,order)
! extracts lowest order modes from periodic 1d scalar field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:), pointer :: pot
         complex, dimension(:), pointer :: pott
! local data
         integer :: nxe, it, nt2, modesxd, inorder
         nxe = size(pot,1)
         it = 1; nt2 = 2
         modesxd = size(pott,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call GTMODES1(pot(1),pott,nx,it,modesx,nxe,nt2,modesxd)
         else
            call GTMODES1(pot(2),pott,nx,it,modesx,nxe,nt2,modesxd)
         endif
         end subroutine igtmodes1
!
         subroutine iptmodes1(pot,pott,nx,modesx,order)
! extracts lowest order modes to periodic 1d scalar field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:), pointer :: pot
         complex, dimension(:), pointer :: pott
! local data
         integer :: nxe, it, nt2, modesxd, inorder
         nxe = size(pot,1)
         it = 1; nt2 = 2
         modesxd = size(pott,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call PTMODES1(pot(1),pott,nx,it,modesx,nxe,nt2,modesxd)
         else
            call PTMODES1(pot(2),pott,nx,it,modesx,nxe,nt2,modesxd)
         endif
         end subroutine iptmodes1
!
         subroutine igtvmodes1(vpot,vpott,nx,modesx,order)
! extracts lowest order modes from periodic 1d vector field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:,:), pointer :: vpot
         complex, dimension(:,:), pointer :: vpott
! local data
         integer :: ndim, nxvh, it, nt, modesxd, inorder
         ndim = size(vpot,1); nxvh = size(vpot,2)/2
         it = 1; nt = 1
         modesxd = size(vpott,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call GTVMODES1(vpot(1,1),vpott,nx,it,modesx,ndim,nxvh,nt,mod&
     &esxd)
         else
            call GTVMODES1(vpot(1,2),vpott,nx,it,modesx,ndim,nxvh,nt,mod&
     &esxd)
         endif
         end subroutine igtvmodes1
!
         subroutine iptvmodes1(vpot,vpott,nx,modesx,order)
! extracts lowest order modes to periodic 1d vector field
         implicit none
         integer :: nx, modesx
         integer, optional :: order
         real, dimension(:,:), pointer :: vpot
         complex, dimension(:,:), pointer :: vpott
! local data
         integer :: ndim, nxvh, it, nt, modesxd, inorder
         ndim = size(vpot,1); nxvh = size(vpot,2)/2
         it = 1; nt = 1
         modesxd = size(vpott,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call PTVMODES1(vpot(1,1),vpott,nx,it,modesx,ndim,nxvh,nt,mod&
     &esxd)
         else
            call PTVMODES1(vpot(1,2),vpott,nx,it,modesx,ndim,nxvh,nt,mod&
     &esxd)
         endif
         end subroutine iptvmodes1
!
         subroutine isfguard1(cus,cu,q2m0,nx,inorder)
! initialize periodic 1d vector field with scaled field
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
               call SCFGUARD1L(cus,cu,q2m0,nx,nxe)
            else
               call SCFGUARD1(cus,cu,q2m0,nx,nxe)
            endif
         end select
         end subroutine isfguard1
!
         subroutine idcuperp13(dcu,amu,nx,inorder)
! calculates the transverse part of periodic 1d vector field
! from momentum flux tensor
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: dcu, amu
! local data
         integer :: nxvh, order
         nxvh = size(dcu,2)/2
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call DCUPERP13(dcu(1,1),amu(1,1),nx,nxvh)
            else
               call DCUPERP13(dcu(1,2),amu(1,2),nx,nxvh)
            endif
         end select
         end subroutine idcuperp13
!
         subroutine iadcuperp13(dcu,amu,nx,inorder)
! calculates the transverse part of periodic 1d vector field
! from acceleration vector and momentum flux tensor
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:), pointer :: dcu, amu
! local data
         integer :: nxvh, order
         nxvh = size(dcu,2)/2
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call ADCUPERP13(dcu(1,1),amu(1,1),nx,nxvh)
            else
               call ADCUPERP13(dcu(1,2),amu(1,2),nx,nxvh)
            endif
         end select
         end subroutine iadcuperp13
!
         subroutine iepois13init(ffe,ax,affp,wp0,ci,nx)
! initialize 1d periodic transverse electric field solver
         implicit none
         integer :: nx
         real :: ax, affp, wp0, ci
         complex, dimension(:), pointer :: ffe
! local data
         integer :: isign = 0, nxvh = 1, nxhd
         real :: wf, dcu, eyz
         nxhd = size(ffe,1)
         call EPOIS13(dcu,eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,nxvh,nxhd)
         end subroutine iepois13init
!
         subroutine iepois13(dcu,eyz,ffe,ci,wf,nx,inorder)
! calculates transverse electric field for periodic 1d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: ci, wf
         real, dimension(:,:), pointer :: dcu, eyz
         complex, dimension(:), pointer :: ffe
! local data
         integer :: isign = -1, nxvh, nxhd, order
         real :: ax, affp, wp0
         nxvh = size(dcu,2)/2
         nxhd = size(ffe,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call EPOIS13(dcu(1,1),eyz(1,1),isign,ffe,ax,affp,wp0,ci,w&
     &f,nx,nxvh,nxhd)
            else
               call EPOIS13(dcu(1,2),eyz(1,2),isign,ffe,ax,affp,wp0,ci,w&
     &f,nx,nxvh,nxhd)
            endif
         end select
         end subroutine iepois13
!
         subroutine iiepois13(dcu,eyz,ffe,ci,wf,nx,inorder)
! calculates transverse electric field for periodic 1d vector field
! without smoothing
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: ci, wf
         real, dimension(:,:), pointer :: dcu
         complex, dimension(:,:), pointer :: eyz
         complex, dimension(:), pointer :: ffe
! local data
         integer :: isign = 1, nxvh, nxhd, order
         real :: ax, affp, wp0
         nxvh = size(dcu,2)/2
         nxhd = size(ffe,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call EPOIS13(dcu(1,1),eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,&
     &nxvh,nxhd)
            else
               call EPOIS13(dcu(1,2),eyz,isign,ffe,ax,affp,wp0,ci,wf,nx,&
     &nxvh,nxhd)
            endif
         end select
         end subroutine iiepois13
!
         subroutine iwpmxn1(qe,qi0,qbme,wpmax,wpmin,nx,inorder)
! calculates maximum and minimum plasma frequency
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qi0, qbme, wpmax, wpmin
         real, dimension(:), pointer :: qe
! local data
         integer :: nxe, order
         nxe = size(qe,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call WPMXN1(qe(1),qi0,qbme,wpmax,wpmin,nx,nxe)
         else
            call WPMXN1(qe(2),qi0,qbme,wpmax,wpmin,nx,nxe)
         endif
         end subroutine iwpmxn1
!
         subroutine ibaddext1(byz,omy,omz,nx,inorder)
! adds constant to magnetic field
         implicit none
         integer :: nx
         real :: omy, omz
         integer, optional :: inorder
         real, dimension(:,:), pointer :: byz
! local data
         integer :: nxe, order
         nxe = size(byz,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(byz,1))
         case (2)
            if (order==LINEAR) then
               call BADDEXT1(byz(1,1),omy,omz,nx,nxe)
            else
               call BADDEXT1(byz(1,2),omy,omz,nx,nxe)
            endif
         end select
         end subroutine ibaddext1
!
         subroutine ivrcopy1(f,g,nx,inorder)
! copies real vector array elements from f to g
         implicit none
         integer nx
         integer, optional :: inorder
         complex, dimension(:,:) :: f
         real, dimension(:,:) :: g
! local data
         integer :: ndim, nxv, order
         ndim = size(g,1); nxv = size(g,2)
         if (nx > nxv) return
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call VRCOPY1(f,g(1,1),nx,ndim,nxv)
         else
            call VRCOPY1(f,g(1,2),nx,ndim,nxv)
         endif
         end subroutine ivrcopy1
!
         subroutine ivccopy1(f,g,nx,inorder)
! copies complex vector array elements from f to g
         implicit none
         integer nx
         integer, optional :: inorder
         real, dimension(:,:) :: f
         complex, dimension(:,:) :: g
! local data
         integer :: ndim, nxv, order
         ndim = size(g,1); nxv = size(g,2)
         if (nx > nxv) return
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call VCCOPY1(f(1,1),g,nx,ndim,nxv)
         else
            call VCCOPY1(f(1,2),g,nx,ndim,nxv)
         endif
         end subroutine ivccopy1
!
      end module field1d
