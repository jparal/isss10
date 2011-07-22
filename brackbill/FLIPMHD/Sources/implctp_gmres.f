      subroutine diagonal(ncells,ijkcell,NSTART,    &
          iwid,jwid,kwid,    &
          divu,dt,coef,vvol,vol,diag)
!
!     a routine to calculate the diagonal elements of the
!     matrix representing the Laplacian
!
      use vast_kind_param, ONLY : double
      use geometry_com_M
      implicit none
      integer :: ijkcell(*)
      integer :: iwid, jwid, kwid
      real(double) ::  dt,vvol(*),vol(*),diag(*),divu(*),coef(*)
!
      integer :: NSTART,NCELLS,n,ijk
      do 1 n=NSTART,NCELLS
!
      ijk=ijkcell(n)
!
      diag(ijk)=    &
         +((c1x(ijk)**2+c1y(ijk)**2+c1z(ijk)**2)/vvol(ijk+iwid)    &
         +(c2x(ijk)**2+c2y(ijk)**2+c2z(ijk)**2)/vvol(ijk+iwid+jwid)    &
         +(c3x(ijk)**2+c3y(ijk)**2+c3z(ijk)**2)/vvol(ijk+jwid)    &
         +(c4x(ijk)**2+c4y(ijk)**2+c4z(ijk)**2)/vvol(ijk)    &
         +(c5x(ijk)**2+c5y(ijk)**2+c5z(ijk)**2)/vvol(ijk+iwid+kwid)    &
        +(c6x(ijk)**2+c6y(ijk)**2+c6z(ijk)**2)/vvol(ijk+iwid+jwid+kwid)    &
         +(c7x(ijk)**2+c7y(ijk)**2+c7z(ijk)**2)/vvol(ijk+jwid+kwid)    &
         +(c8x(ijk)**2+c8y(ijk)**2+c8z(ijk)**2)/vvol(ijk+kwid))    &
         +vol(ijk)/coef(ijk)
!
    1 continue
!
      do 2 n=1,NSTART-1
!
      ijk=ijkcell(n)
!
      diag(ijk)=    &
         +((c1x(ijk)**2+c1y(ijk)**2+c1z(ijk)**2)/vvol(ijk+iwid)    &
         +(c2x(ijk)**2+c2y(ijk)**2+c2z(ijk)**2)/vvol(ijk+iwid+jwid)    &
         +(c3x(ijk)**2+c3y(ijk)**2+c3z(ijk)**2)/vvol(ijk+jwid)    &
         +(c4x(ijk)**2+c4y(ijk)**2+c4z(ijk)**2)/vvol(ijk)    &
         +(c5x(ijk)**2+c5y(ijk)**2+c5z(ijk)**2)/vvol(ijk+iwid+kwid)    &
        +(c6x(ijk)**2+c6y(ijk)**2+c6z(ijk)**2)/vvol(ijk+iwid+jwid+kwid)    &
         +(c7x(ijk)**2+c7y(ijk)**2+c7z(ijk)**2)/vvol(ijk+jwid+kwid)    &
         +(c8x(ijk)**2+c8y(ijk)**2+c8z(ijk)**2)/vvol(ijk+kwid))    &
         +vol(ijk)/coef(ijk)
!
    2 continue
!
      return
       end subroutine diagonal
!
      subroutine poisson(ncells,ijkcell,nvtxkm,nvtx,ijkvtx,    &
          iwid,jwid,kwid,PRECON,gm1,mv,    &
          ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,    &
          vol,vvol,vvolaxis,    &
          itsub,itmax,error,rnorm,DT,CNTR,srce,rho,divu,csq,    &
          coef,div,residu,Aq,q,bnorm,    &
          gradpx,gradpy,gradpz,p,diag)
!
       use vast_kind_param, ONLY : double
      use geometry_com_M
      implicit none
!
      integer :: ncells,ijkcell(*),nvtx,ijkvtx(*)
      integer :: iwid, jwid, kwid, ibp1, jbp1, kbp1, nvtxkm
      integer :: itsub, itmax
      real(double) :: DT, CNTR, cdlt, sdlt, strait, dz, gm1,   &
          rnorm, bnorm, error, dtsq, dtcntr
      real(double) :: p(*),vvol(*),vvolaxis(*),    &
               Aq(*),q(itdim,20),    &
               srce(*),rho(*),divu(*),csq(*),coef(*),    &
               vol(*),div(*),residu(*),    &
               gradpx(*),gradpy(*),gradpz(*),diag(*)
!
      real(double) ::  mv(*)
      integer :: nstart, ijk, n, iter, it
      LOGICAL :: PRECON
!
!
      nstart=(ibp1-1)*(jbp1-1)+1
      DTSQ=cntr*DT**2
      dtcntr=dt*cntr
!
!     set source term
!
      call bcphi(ibp1,jbp1,kbp1,p)
!
      do 1115 n=1,ncells
      ijk=ijkcell(n)
      coef(ijk)=csq(ijk)*rho(ijk)*dtsq+1.e-10
 1115 continue
!
      do 1116 n=1,ncells
      ijk=ijkcell(n)
      srce(ijk)=p(ijk)/coef(ijk)-divu(ijk)/(dtcntr)
 1116 continue
!
      do  n=1,ncells
      ijk=ijkcell(n)
      diag(ijk)=0.0
      enddo
!
      call bcphi(ibp1,jbp1,kbp1,diag)
!
      call residue(ncells,ijkcell,nvtx,nvtxkm,ijkvtx,    &
          iwid,jwid,kwid,ibp1,jbp1,kbp1,    &
          vol,    &
          srce,coef,mv,gm1,div,divu,    &
          p,gradpx,gradpy,gradpz,vvol,residu)
!
!
      bnorm=0.0
      do n=1,ncells
      ijk=ijkcell(n)
      bnorm=bnorm+abs(residu(ijk))
      enddo
!
      if(bnorm.lt.error) then
      bnorm=1.0
      end if

!     set preconditioner
!
      if(precon) then
!
      call diagonal(ncells,ijkcell,NSTART,    &
          iwid,jwid,kwid,    &
          divu,dt,coef,mv,vol,diag)
!
      else
!
      DO 1111 N=1,NCELLS
      DIAG(IJKCELL(N))=1.
 1111 CONTINUE
!
      endif
!
      do it=1,itmax
      call gmres(ncells,ijkcell,nvtxkm,nvtx,ijkvtx,    &
          ITDIM,iwid,jwid,kwid,PRECON,gm1,mv,    &
          ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,    &
          vol,vvol,vvolaxis,    &
          itsub,iter,error,rnorm,DT,CNTR,srce,rho,divu,csq,    &
          coef,div,residu,Aq,q,bnorm,    &
          gradpx,gradpy,gradpz,p,diag)
!
       if(iter.lt.itsub) go to 98
       end do
!
       write(0,*) 'gmres fails', (it-1)*itsub+iter, rnorm
       go to 99
  98   continue
       write(0,*) 'gmres converges ', (it-1)*itsub+iter, rnorm
  99   continue
!
      return
      end subroutine poisson
!
      subroutine gmres(ncells,ijkcell,nvtxkm,nvtx,ijkvtx,    &
          iwid,jwid,kwid,PRECON,gm1,mv,    &
          ibp1,jbp1,kbp1,cdlt,sdlt,strait,dz,    &
          vol,vvol,vvolaxis,    &
          itsub,iter,error,rnorm,DT,CNTR,srce,rho,divu,csq,    &
          coef,div,residu,Aq,q,bnorm,    &
          gradpx,gradpy,gradpz,p,diag)
!
!
      use vast_kind_param, ONLY : double
      use geometry_com_M
      implicit none
!
      integer :: ncells, ijkcell(*), nvtx, ijkvtx(*)
      integer :: iter, itsub
      integer :: iwid, jwid, kwid, ibp1, jbp1, kbp1, nvtxkm
      real(double) :: gm1, rnorm, DT, CNTR 
      real(double) ::     &
               p(*),vvol(*),vvolaxis(*),     &
               Aq(*),q(itdim,20),     &
               srce(*),rho(*),divu(*),csq(*),coef(*),     &
               vol(*),div(*),residu(*),     &
               gradpx(*),gradpy(*),gradpz(*),diag(*)
!
       real(double) ::  mv(*)
       real(double) :: g(21,2),s(21),h(21,21)
       real(double) :: dot, dtcntr, dtsq, ws, g1, g2, bnorm, error,  &
          sdlt, cdlt, dz, strait
       integer :: nstart, iterp1, ii, jj, i, j, k, ijk, n
!
      LOGICAL PRECON
!
!
      nstart=(ibp1-1)*(jbp1-1)+1
      DTSQ=cntr*DT**2
      dtcntr=dt*cntr
!
!
!     ZERO WORKING ARRAYS
!
      do 1112 n=1,ncells
      ijk=ijkcell(n)
      Aq(ijk)=0.0
 1112 continue
!
      do 1114 ii=1,itsub+1
      s(ii)=0.0
      g(ii,1)=0.0
      g(ii,2)=0.0
      do 1114 jj=1,itsub+1
         h(ii,jj)=0.0
 1114 continue
!
      do 1113 ii=1,itsub
      do 1113 n=1,ncells
      ijk=ijkcell(n)
      q(ijk,ii)=0.0
 1113 continue
!
!
!     CALCULATE THE INITIAL RESIDUAL ERROR
!
      call residue(ncells,ijkcell,nvtx,nvtxkm,ijkvtx,     &
          iwid,jwid,kwid,ibp1,jbp1,kbp1,     &
          vol,     &
          srce,coef,mv,gm1,div,divu,     &
          p,gradpx,gradpy,gradpz,vvol,residu)
!
!
      dot=0.0
!
      do 22 n=1,ncells
      ijk=ijkcell(n)
      dot=dot+residu(ijk)*residu(ijk)
   22 continue
!
      dot=sqrt(dot)
      s(1)=dot
!
      do 25 n=1,ncells
      ijk=ijkcell(n)
      q(ijk,1)=residu(ijk)/dot
  25  continue
!
!     ****************************************************************
!
!     begin gmres
!
!     ****************************************************************
!
      do 1 iter=1,itsub
!     
!     apply preconditioner
!
      do 32 n=1,ncells
      ijk=ijkcell(n)
      Aq(ijk)=-q(ijk,iter)/diag(ijk)
  32  continue
!
!     compute A times preconditioned q 
!
      do 33 n=1,ncells
      ijk=ijkcell(n)
      Aq(ijk)=p(ijk)+Aq(ijk)
   33 continue
      call bcphi(ibp1,jbp1,kbp1,Aq)
!
      call residue(ncells,ijkcell,nvtx,nvtxkm,ijkvtx,     &
          iwid,jwid,kwid,ibp1,jbp1,kbp1,     &
          vol,     &
          srce,coef,mv,gm1,div,divu,     &
          Aq,gradpx,gradpy,gradpz,vvol,Aq)
!
!
!dir$ ivdep
      do 35 n=1,ncells
      ijk=ijkcell(n)
      Aq(ijk)=(residu(ijk)-Aq(ijk))
   35 continue
!
!
!    orthogonalize:
!
      do 13 k=1,iter
      dot=0.0
      do 11 n=1,ncells
      ijk=ijkcell(n)
      dot=dot+Aq(ijk)*q(ijk,k)
   11 continue
      h(k,iter)=dot
!dir$ ivdep
      do 12 n=1,ncells
      ijk=ijkcell(n)
      Aq(ijk)=Aq(ijk)-dot*q(ijk,k)
   12 continue
   13 continue
!
!
      dot=0.0
      do 14 n=1,ncells
      ijk=ijkcell(n)
      dot=dot+Aq(ijk)*Aq(ijk)
   14 continue
      dot=sqrt(dot)
!
      iterp1=iter+1
      h(iterp1,iter)=dot
!
!     apply previous Givens rotations to h (update QR factorization)
!
      do 16 k=1,iter-1
      ws=g(k,1)*h(k,iter)-g(k,2)*h(k+1,iter)
      h(k+1,iter)=g(k,2)*h(k,iter)+g(k,1)*h(k+1,iter)
      h(k,iter)=ws
  16  continue 
!     
!     compute next Givens rotation
      g1=h(iter,iter)
      g2=h(iterp1,iter)
      ws=sqrt(g1*g1+g2*g2)
      g1=g1/ws
      g2=-g2/ws
      g(iter,1)=g1
      g(iter,2)=g2
!
!     apply g to h
      h(iter,iter)=g1*h(iter,iter)-g2*h(iterp1,iter)
      h(iterp1,iter)=0.0
!
!     apply g to s
      ws=g1*s(iter)-g2*s(iterp1)
      s(iterp1)=g2*s(iter)+g1*s(iterp1)
      s(iter)=ws
!
!     |s(iter+1)| is the norm of the current residual
!     check for convergence
      rnorm=abs(s(iterp1))
      if((rnorm.le.error*bnorm).or.(iter.eq.itsub)) go to 2
!
!     normalize next q
!
      do 15 n=1,ncells
      ijk=ijkcell(n)
      q(ijk,iterp1)=Aq(ijk)/dot
  15  continue
!
    1 continue
!
    2 continue
!    
!    update the solution
!    solve h*y=sbar  (sbar contains elements 1,..,iter of s)
!    store y in s
      s(iter)=s(iter)/h(iter,iter)
      do ii=iter-1,1,-1
        ws=0.0
          do jj=ii+1,iter 
          ws=ws+h(ii,jj)*s(jj)
          enddo
        s(ii)=(s(ii)-ws)/h(ii,ii)
      end do
!
!    compute new p
      do 17 n=1,ncells
      ijk=ijkcell(n)
      p(ijk)=-diag(ijk)*p(ijk)
  17  continue
!
      do 18 ii=1,iter
!dir$ ivdep
      do 18 n=1,ncells
      ijk=ijkcell(n)
      p(ijk)=p(ijk)+s(ii)*q(ijk,ii)
  18  continue
!
      do 19 n=1,ncells
      ijk=ijkcell(n)
      p(ijk)=-p(ijk)/diag(ijk)
  19  continue
!
      return
      end subroutine gmres
!
      subroutine residue(ncells,ijkcell,nvtx,nvtxkm,ijkvtx,     &
          iwid,jwid,kwid,ibp1,jbp1,kbp1,     &
          vol,     &
          srce,coef,mv,gm1,div,divu,     &
          p,gradpx,gradpy,gradpz,vvol,residu)
!
!     a routine to calculate the residual error in the solution
!     of poisson's equation
!
      use vast_kind_param
      use geometry_com_M
      implicit none
!
      integer :: ijkcell(*),ijkvtx(*),ncells,nvtx,   &
          ibp1, jbp1, kbp1, iwid, jwid, kwid
      real(double) :: vol(*),vvol(*),p(*),div(*),divu(*),    &
          srce(*),coef(*),mv(*),    &
          gradpx(*),gradpy(*),gradpz(*),residu(*)
      integer :: ijk, n, nvtxkm
      real(double) :: factor, rmv, gm1
        
!
!     calculate residue
!
      call gradf(nvtx,ijkvtx,    &
                       p,gradpx,gradpy,gradpz)
!
      call axisgrad(ibp1,jbp1,    &
          gradpx,gradpy,gradpz)
!
      do 16 n=1,nvtxkm
      ijk=ijkvtx(n)
      rmv=1./mv(ijk)
      gradpx(ijk)=-gradpx(ijk)*rmv
      gradpy(ijk)=-gradpy(ijk)*rmv
      gradpz(ijk)=-gradpz(ijk)*rmv
   16 continue
      factor=0.5*(1.+1./real(2*kbp1))
      do 17 n=nvtxkm+1,nvtx
      ijk=ijkvtx(n)
      rmv=factor/mv(ijk)
      gradpx(ijk)=-gradpx(ijk)*rmv
      gradpy(ijk)=-gradpy(ijk)*rmv
      gradpz(ijk)=-gradpz(ijk)*rmv
   17 continue
!
!
!     calculate div(grad(p)/rho)
!
      call divc(ncells,ijkcell,iwid,jwid,kwid,     &
                        vol,     &
                        GRADPX,GRADPY,GRADPZ,div)
!
!
      do 15 n=1,ncells
      ijk=ijkcell(n)
      residu(ijk)=(srce(ijk)-p(ijk)/coef(ijk)-div(ijk))*vol(ijk)
   15 continue
!
      return
      end subroutine residue
