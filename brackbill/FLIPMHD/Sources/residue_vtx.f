      subroutine residue_vtx(nvtx,ijkvtx,    &
          srce,    &
          p,residu)
!
!     a routine to calculate the residual error in the solution
!     of poisson's equation
!
      use vast_kind_param, ONLY : double
      use geometry_com_M
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z,   &
        ijkcell,vol,    &
        gradx,grady,gradz,divu
      use cindex_com_M, ONLY : ncells,iwid,jwid,kwid,ibp1,jbp1,kbp1
      use cophys_com_M, ONLY : cdlt,sdlt,dz
      use corgan_com_M, ONLY : strait
      implicit none
!
!
      real(double) ::      &
          p(*),srce(*),     &
          residu(*)
      real(double) :: zero
      integer :: n,ijk,ibp2,jbp2,kbp2
      integer :: nvtx, ijkvtx(*)
!
!     calculate residue
!     NOTE:  if periodic in x or y, the residue is calculated only
!     at vertices i=2,ibp1; j=2,jbp1
!
!     with periodic boundary conditions, BC_SCALAR
!     imposes periodicity on p 
!
      ibp2=ibp1+1
      jbp2=jbp1+1
      kbp2=kbp1+1
      call bc_vertex(ibp2,jbp2,kbp2,     &
          p)
!
      call gradc(ncells,ijkcell,     &
          p,gradx,grady,gradz)
!
      zero=0.0
!
!     if periodic in x and y, TORUSBC sets values of
!     gradpj in i=1,j=1 logical planes, which allows correct 
!     calculation of div at in i=2, j=2 logical planes
!
      call torusbc(ibp1+1,jbp1+1,kbp1+1,    &
          zero,                             &
          gradx,grady,gradz)
!
!     calculate div(grad(p))
!
      call divv(   &
          gradx,grady,gradz,divu)
!
      do 15 n=1,nvtx
      ijk=ijkvtx(n)
      residu(ijk)=srce(ijk)-divu(ijk)
   15 continue
!
      return
      end subroutine residue_vtx
