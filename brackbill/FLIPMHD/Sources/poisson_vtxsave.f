      subroutine poisson_vtx(ncells,ijkcell,nvtx,ijkvtx,
     &     periodic_x,periodic_y,periodic_z,
     &     ITDIM,iwid,jwid,kwid,PRECON,
     &     ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,
     &     vol,
     &     itsub,itmax,error,rnorm,DT,CNTR,srce,
     &     div,residu,Aq,q,bnorm,
     &     gradpx,gradpy,gradpz,p,diag)
c
c      POISSON_VTX is called by VINIT
c
c      POISSON_VTX calls BC_SCALAR, RESIDUE_VTX, DIAGONAL_VTX,
c                       GMRES_VTX
c
      use geometry_com_M
      implicit real*8 (a-h,o-z)
c
      dimension ijkcell(*),ijkvtx(*),
     &          p(*),
     &          Aq(*),q(itdim,20),
     &          srce(*),
     &          vol(*),div(*),residu(*),
     &          gradpx(*),gradpy(*),gradpz(*),diag(*)
c
      LOGICAL PRECON,periodic_x,periodic_y,periodic_z
C
      call residue_vtx(ncells,ijkcell,nvtx,ijkvtx,
     &     periodic_x,periodic_y,periodic_z,
     &     cdlt,sdlt,strait,dz,
     &     iwid,jwid,kwid,ibp1,jbp1,kbp1,
     &     vol,
     &     div,srce,
     &     p,gradpx,gradpy,gradpz,residu)
c
c
      bnorm=0.0
      rl2_norm=0.0
      src_norm=0.0
      do n=1,nvtx
      ijk=ijkvtx(n)
      rl2_norm=rl2_norm+ residu(ijk)*residu(ijk)
      bnorm=bnorm+abs(residu(ijk))
      src_norm = src_norm + abs(srce(ijk))
      enddo
c
      write(*,*)'in poisson_vtx, l2_norm, bnorm, src_norm =', 
     .           sqrt(rl2_norm), bnorm, src_norm
      if(bnorm.lt.error) then
      bnorm=1.0
      end if
c      if(bnorm*error.lt.1.0e-6) then
c       error = 1.0e-6/bnorm
c      end if
 
c     set preconditioner
c
      if(precon) then
c
      call diagonal_vtx(nvtx,ijkvtx,
     &     iwid,jwid,kwid,
     &     dt,vol,diag)
c
      else
c
         DO 1111 n=1,nvtx
            diag(ijkvtx(n))=1.
 1111 CONTINUE
C
      endif

      do n=1,nvtx
         diag(ijkvtx(n))=1.
      end do

c
      do it=1,itmax
      call gmres_vtx(ncells,ijkcell,nvtx,ijkvtx,
     &     periodic_x,periodic_y,periodic_z,
     &     iwid,jwid,kwid,PRECON,
     &     ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,
     &     vol,
     &     itsub,iter,error,rnorm,DT,CNTR,srce,
     &     div,residu,Aq,q,bnorm,
     &     gradpx,gradpy,gradpz,p,diag)
c
       write(*,*) 'gmres its, res, bnorm ',
     &             (it-1)*itsub+iter, rnorm,bnorm
       write(8,*) 'gmres its, res, bnorm ',
     &             (it-1)*itsub+iter, rnorm,bnorm
       if(iter.lt.itsub) go to 98
       end do
c
       write(0,*) 'gmres fails', (it-1)*itsub+iter, rnorm
       go to 99
  98   continue
       write(8,*) 'gmres converges ', (it-1)*itsub+iter, rnorm
  99   continue
c
      return
      end
