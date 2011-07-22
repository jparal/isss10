      subroutine poisson_cg(ncells,ijkcell,nvtx,ijkvtx,     &
          periodic_x,periodic_y,periodic_z,     &
          iwid,jwid,kwid,PRECON,     &
          ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,     &
          vol,vvol,vvolaxis,     &
          q,qtilde,Aqtilde,itmax,error,srce,     &
          residu,Aq,gradpx,gradpy,gradpz,p,diag)
!
!     a routine to solve poisson's equation for the pressure
!     USING A CONJUGATE RESIDUAL SOLVER
!     boundary conditions correspond to a topological torus
!     that is periodic in i and j, and with Dirichlet bc's at k=kbp2
!
      use vast_kind_param, only:  double
      use geometry_com_M
      real(double) :: lambda
      logical :: periodic_x,periodic_y,periodic_z
!
      integer :: ncells,ijkcell(*),nvtx,ijkvtx(*)
      real(double) ::     &
               q(*),qtilde(*),p(*),vvol(*),vvolaxis(*),     &
               Aqtilde(*),Aq(*),srce(*),     &
               vol(*),residu(*),gradpx(*),gradpy(*),gradpz(*),diag(*)
!
!
!     LOGICAL PRECON
!
      write(*,*)'ncells,nvtx,itmax,error=',ncells,nvtx,itmax,error
!
!      if(precon) then
!
      call diagonal_vtx(nvtx,ijkvtx)
!
!      else
!
      DO 1111 N=1,NVTX
      DIAG(IJKvtx(N))=1.
 1111 CONTINUE
!
!      endif
!
      write(*,*) 'done with precond'
!
      do 1112 k=1,kbp1+1
      do 1112 j=1,jbp1+1
      do 1112 i=1,ibp1+1
      ijk=1+(i-1)*iwid+(j-1)*jwid+(k-1)*kwid
      q(ijk)=0.0
      Aq(ijk)=0.0
      qtilde(ijk)=0.0
      residu(ijk)=0.0
 1112 continue
!
      call bc_scalar(ibp1+1,jbp1+1,kbp1+1,   &
          periodic_x,periodic_y,periodic_z,   &
          p)
!
      write(*,*) 'CALCULATE THE INITIAL RESIDUAL ERROR'
!
!     CALCULATE THE INITIAL RESIDUAL ERROR
!
      call residue_vtx(     &
          srce,     &
          p,residu)
!
      write(*,*) 'CALCULATE THE INITIAL RESIDUAL ERROR:after'
!
      rsum=0.0
      bsum=0.0
!
      do 32 n=1,nvtx
      ijk=ijkvtx(n)
      bsum=bsum+abs(srce(ijk))
      rsum=rsum+abs(residu(ijk))
   32 continue
!
      if(rsum.le.error) then
      go to 2
      endif
!
      if(bsum.le.error) bsum=1.0
!
      do 33 n=1,nvtx
      ijk=ijkvtx(n)
      qtilde(ijk)=residu(ijk)/diag(ijk)
      q(ijk)=qtilde(ijk)+p(ijk)
   33 continue
!
      call bc_scalar(ibp1+1,jbp1+1,kbp1+1,     &
          periodic_x,periodic_y,periodic_z,     &
          q)
!
      call residue_vtx(     &
          srce,     &
          q,Aq)
!
      rdot=0.0
      do 34 n=1,nvtx
      ijk=ijkvtx(n)
      q(ijk)=q(ijk)-p(ijk)
      Aq(ijk)=residu(ijk)-Aq(ijk)
      rdot=rdot+residu(ijk)*qtilde(ijk)
   34 continue
!
!     ****************************************************************
!
!     begin conjugate gradient iteration
       write(*,*) 'begin conjugate gradient iteration'
!
!     ****************************************************************
!
      do 1 iter=1,itmax
!
      alfa=0.0
      do 5 n=1,nvtx
      ijk=ijkvtx(n)
      alfa=alfa+q(ijk)*Aq(ijk)
  5   continue
!
      alfa=rdot/alfa
!
      do 14 n=1,nvtx
      ijk=ijkvtx(n)
      p(ijk)=p(ijk)+ALFA*q(ijk)
      residu(ijk)=residu(ijk)-ALFA*Aq(ijk)
   14 continue

!
      rnorm=0.0
      dxnorm=0.0
      xnorm=0.0
      do 13 n=1,nvtx
      ijk=ijkvtx(n)
      rnorm=rnorm+abs(residu(ijk))
      dxnorm=dxnorm+abs(q(ijk))
      xnorm=xnorm+abs(p(ijk))
   13 continue
      write(*,*) 'it,rnorm=',iter,rnorm
      dxnorm=dxnorm*abs(alfa)
!
      if(dxnorm.le.error*xnorm.and.rnorm.le.error*bsum)   &
          go to 2
!
!     apply pre-conditioner
      do 10 n=1,nvtx
      ijk=ijkvtx(n)
      qtilde(ijk)=residu(ijk)/diag(ijk)
   10 continue
!
!     calculate beta
!
!
      alfa=0.0
      do 11 n=1,nvtx
      ijk=ijkvtx(n)
      alfa=alfa+qtilde(ijk)*residu(ijk)
   11 continue
!
      beta=alfa/rdot
      rdot=alfa
!
!dir$ ivdep
      do 12 n=1,nvtx
      ijk=ijkvtx(n)
      q(ijk)=qtilde(ijk)+beta*q(ijk)+p(ijk)
   12 continue
!
      call bc_scalar(ibp1+1,jbp1+1,kbp1+1,     &
          periodic_x,periodic_y,periodic_z,     &
          q)
!
!     calculate Aq
!
!
      call residue_vtx(     &
          srce,     &
          q,     &
          Aq)
!
!
      do 15 n=1,nvtx
      ijk=ijkvtx(n)
      Aq(ijk)=(residu(ijk)-Aq(ijk))
      q(ijk)=q(ijk)-p(ijk)
   15 continue
!
    1 continue
!
      write(6,*)'poisson: iteration fails '
      write(6,*)'dxnorm,xnorm,rnorm,bsum'
      write(6,*) dxnorm,xnorm,rnorm,bsum
      return
!
    2 continue
!     iteration has converged
      write(6,*)'poisson: converges in ',iter
!
      call bc_scalar(ibp1+1,jbp1+1,kbp1+1,    &
          periodic_x,periodic_y,periodic_z,    &
          p)
!
      return
      end subroutine poisson_cg
