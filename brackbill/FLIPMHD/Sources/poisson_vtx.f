      subroutine poisson_vtx(nvtx,ijkvtx,refnorm,    &
          PRECON,    &
          itsub,error,rnorm,srce,    &
          bnorm, iter  &
          )
!
!      POISSON_VTX is called by VINIT
!
!      POISSON_VTX calls BC_SCALAR, RESIDUE_VTX, DIAGONAL_VTX,
!                       GMRES_VTX
!
      use vast_kind_param, ONLY : double
      use cindex_com_M, ONLY : ibar,jbar,kbar,ibp1,jbp1,kbp1,    &
        ncells
      use corgan_com_M, ONLY : itdim, strait,cntr
      use cophys_com_M, ONLY : cdlt,sdlt,dz
      use blcom_com_M, ONLY : periodic_x,periodic_y,periodic_z,vol,itmax,   &
         ijkcell,     &
         gradx,grady,gradz
      use numpar_com_M, ONLY : dt
      use gmres_com_M, ONLY : phi,diag,q,residu
      use Timing
!      implicit real*8 (a-h,o-z)
      implicit none
!
      real(double) :: refnorm,            &
               error, rnorm, bnorm, src_norm,  &
               srce(*)
      real(double) :: Tstart, Tfinish
      integer :: itsub, ijk, it, iter, n, m
      integer :: nvtx, ijkvtx(*)
!
      LOGICAL PRECON
      logical preconsav
!
      call cpu_time(Tstart)
    
      preconsav=precon
      precon=.true.
      precon=.false.
!
      call residue_vtx(    &
          srce,    &
          phi,residu)
!
      call L2NormF(srce, src_norm,   &
        1,nvtx,ijkvtx)
!
      call L2NormF(residu, bnorm,   &
        1,nvtx,ijkvtx)
!
      if(bnorm.lt.error*refnorm) then
        write(*,*) 'poisson_vtx: no iteration, bnorm, refnorm=', bnorm, refnorm
        return
      endif
 
!     no preconditioner
         DO 1111 n=1,nvtx
            diag(ijkvtx(n))=1.
 1111 CONTINUE
!
!
!      do it=1,itmax
       it=0
       do while(it.le.itmax.and.rnorm.gt.error*refnorm)
        it=it+1
!
        call residue_vtx(     &
            srce,     &
            phi,residu)
!
        call gmres_vtx(nvtx,ijkvtx,refnorm,     &
            PRECON,     &
            vol,     &
            itsub,iter,error,rnorm,srce,     &
            residu,bnorm,     &
            gradx,grady,gradz)
!
         write(55,*) 'gmres its, res, bnorm ',     &
                  (it-1)*itsub+iter, rnorm,bnorm
         write(8,*) 'gmres its, res, bnorm ',      &
                  (it-1)*itsub+iter, rnorm,bnorm
   !      if(iter.lt.itsub) go to 98
       end do
!
       if(it.eq.itmax.and.rnorm.gt.error*refnorm) then
       write(0,*) 'gmres fails', (it-1)*itsub+iter, rnorm, bnorm
       else
!       go to 99
!  98   continue
       endif
!  99   continue
!
       call cpu_time(Tfinish)
      do n=1,20
         if(RoutineName(n).eq.'poisson_vtx') then
            CPUTime(n)=CPUTime(n)+Tfinish-Tstart
         endif
      enddo
       
      precon=preconsav
!
      return
      end subroutine poisson_vtx
