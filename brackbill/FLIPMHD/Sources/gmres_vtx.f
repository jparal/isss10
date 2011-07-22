      subroutine gmres_vtx(nvtx,ijkvtx,refnorm,    &
          PRECON,    &
          vol,    &
          itsub,iter,error,rnorm,srce,    &
          residu, bnorm,    &
          gradpx,gradpy,gradpz)
!
!
      use vast_kind_param, only:  double
      use corgan_com_M, ONLY : itdim, GMitmax
      use cindex_com_M, ONLY : ncells,  &
         iwid, jwid, kwid, ibp1, jbp1, kbp1
      use blcom_com_M, ONLY : ijkcell,    &
         periodic_x, periodic_y, periodic_z
      use geometry_com_M
      use gmres_com_M, ONLY : q, phi, Jdu,  &
         wKrylov
      use Timing
      implicit none
!
!
      real(double) :: refnorm,        &
               residu(*),             &
               vol(*),srce(*),        &
               gradpx(*),gradpy(*),gradpz(*)
!
       real(double) :: g(21,2),s(21),h(21,21)
!
       real(double) :: error, rnorm, bnorm, dot,    &
               ws, zero, eps, GMRES_tolerance, gamma
       real(double) :: Tstart, Tfinish
       integer :: itsub, iter, i1, ijk, iterp1, itest,   &
               ii, jj, k, l, n, Newton_its, iter_total
       integer :: nvtx, ijkvtx(*)
!
      LOGICAL PRECON
!     
!     solving system A.p=srce
!
!
!     ZERO WORKING ARRAYS, note ensure that **all** of the array is zeroed
!     (use ncells and ijkcell to zero rather than nvtx and ijkvtx); I found that
!     if these arrays are not zeroed completely, non-zero values from the 
!     previous time step will be carried through the calculation and will affect
!     the gradient and divergence calculations. I would assume that correct 
!     implementation of the boundary conditions would prevent any previous
!     non-zero values from leaking thru to the current calculation. Maybe the
!     routine bc_scalar should be checked more closely..

!
!
      call cpu_time(Tstart)
!
      Newton_its=0
      iter_total=0
      do while(rnorm.gt.error*refnorm)
!
         GMRES_tolerance=max(gamma*rnorm,error*refnorm)
         Newton_its=Newton_its+1
!
         do 1112 n=1,nvtx
            ijk=ijkvtx(n)
            wKrylov(ijk)=0.0
 1112    continue

         do 1114 ii=1,itsub+1
            s(ii)=0.0
            g(ii,1)=0.0
            g(ii,2)=0.0
            do jj=1,itsub+1
               h(ii,jj)=0.0
            enddo
 1114    continue
!
      do 1113 ii=1,itsub
      do  n=1,nvtx
         ijk=ijkvtx(n)
         q(ijk,ii)=0.0
      enddo
 1113 continue
!
!
!     CALCULATE THE INITIAL RESIDUAL ERROR, residu=srce-A.p_0 where p_0 is the 
!     initial pressure
!
!
      call residue_vtx(    &
          srce,    &
          phi,residu)
!
      call L2NormF(residu,dot,   &
         1,nvtx,ijkvtx)
!
      s(1)=dot
!
!
!     first Krylov vector..q_1 = (srce-A.p_0)/||srce-A.p_0||

      if(dot.gt.0.0d0) then
      do 25 n=1,nvtx
         ijk=ijkvtx(n)
         q(ijk,1)=-residu(ijk)/dot
  25  continue

      else
!
      do n=1,nvtx
         ijk=ijkvtx(n)
         q(ijk,1)=0.0d0
      enddo
!
      endif
!
!     ****************************************************************
!
!     begin gmres
!
!     ****************************************************************
!
      iter=0
      do while (iter.lt.itsub)
         iter=iter+1
!     
!     apply preconditioner, i.e find y=M^(-1).q
!
    
      call Eps_eval(iter,phi,q,eps,     &
        1,nvtx,ijkvtx)
!
      call OneJacobi(iter,eps,phi,q,residu,Jdu)
!
      do 32 n=1,nvtx
         ijk=ijkvtx(n)
         wKrylov(ijk)=Jdu(ijk)
  32  continue
!
!
      zero=0.0d0
!
!     orthogonalize this with previous
!     Krylov vectors
!
      do 13 k=1,iter
         dot=0.0d0
         do 11 n=1,nvtx
            ijk=ijkvtx(n)
            dot=dot+wKrylov(ijk)*q(ijk,k)
         11 continue
            h(k,iter)=dot
      13 continue
!
      do k=1,iter
         do 12 n=1,nvtx
            ijk=ijkvtx(n)
            wKrylov(ijk)=wKrylov(ijk)-dot*q(ijk,k)
         12 continue
      enddo
!
!
      call L2NormF(wKrylov,dot,   &
        1,nvtx,ijkvtx)
!
      iterp1=iter+1
      h(iterp1,iter)=dot
!
!     normalize next Krylov vector q
!
      do 15 n=1,nvtx
      ijk=ijkvtx(n)
      q(ijk,iterp1)=wKrylov(ijk)/h(iterp1,iter)
  15  continue

!     apply previous Givens rotations to h (update QR factorization)
!
      do 16 k=1,iter-1
         ws=g(k,1)*h(k,iter)-g(k,2)*h(k+1,iter)
         h(k+1,iter)=g(k,2)*h(k,iter)+g(k,1)*h(k+1,iter)
         h(k,iter)=ws
  16  continue 
!     
!     compute next Givens rotation
!
      ws=sqrt(h(iter,iter)**2+h(iterp1,iter)**2)+1.d-50
      g(iter,1)=h(iter,iter)/ws
      g(iter,2)=-h(iterp1,iter)/ws
!
!     apply g to h
      h(iter,iter)=ws
      h(iterp1,iter)=0.0
!
!     apply g to s
      ws=g(iter,1)*s(iter)-g(iter,2)*s(iterp1)
      s(iterp1)=g(iter,2)*s(iter)+g(iter,1)*s(iterp1)
      s(iter)=ws

!SJC  |s(iter+1)| is the norm of the current residual, but these
!     norms are based on the solution of the system Abar.y = srce 
!     where Abar = A(M^(-1))
!
!     |s(iter+1)| is the norm of the current residual
!     check for convergence
     rnorm=abs(s(iterp1))
      if(rnorm.le.GMRES_tolerance) go to 2
!
     enddo
!
    2 continue
!    
!    update the solution
!    solve h*y=sbar  (sbar contains elements 1,..,iter of s)
!    store y in s
      s(iter)=s(iter)/(h(iter,iter)+1.d-50)
      do ii=iter-1,1,-1
        ws=0.0
          do jj=ii+1,iter 
          ws=ws+h(ii,jj)*s(jj)
          enddo
        s(ii)=(s(ii)-ws)/h(ii,ii)
      end do

!     
      do  n=1,nvtx
         ijk=ijkvtx(n)
         ws=0.0
         do  ii=1,iter
            ws=ws+s(ii)*q(ijk,ii)
         enddo
         phi(ijk)=phi(ijk)+ws
      enddo
 
!

!
!  compute norm of actual residual at end of cycle, ||srce-A.p||
!
      call residue_vtx(    &
          srce,    &
          phi,residu)
!
      call L2NormF(residu,rnorm,    &
        1,nvtx,ijkvtx)
!
      iter_total=iter_total+iter
!
      enddo
!
      call cpu_time(Tfinish)
      do l=1,20
         if(RoutineName(l).eq.'gmres_vtx') CPUTime(l)=CPUTime(l)+Tfinish-Tstart
      enddo

      return
      end subroutine gmres_vtx
