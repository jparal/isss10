      subroutine gmres_3dmhd(tol_GMRES,ncells)
!
       USE vast_kind_param, ONLY:  double
       USE nkmhd_com_M
       use geometry_com_M
       use nk_com_M
       use Timing
      implicit none
!
      integer ncells
!
      real(double):: tol_GMRES
      real(double) :: Tstart, Tfinish
!
!	local variables
!
      real(double) ::    &
          ws,g1,g2,dot,sumpq
!
      integer numrstrt, iter, l, n, ii,   &
          iterp1,k,nn,jj
!
!   *************************************************************
!
       call cpu_time(Tstart)
       numrstrt = 0
 1000  continue
!
!     ZERO WORKING ARRAYS
!
      do n=1,num_eq*ncells
         Jq(n)=0.0
      enddo
!
      do ii=1,itsub_mhd+1
         s_local(ii)=0.0
         g(ii,1)=0.0
         g(ii,2)=0.0
         do jj=1,itsub_mhd+1
            h(ii,jj)=0.0
         enddo
      enddo
!
      do ii=1,itsub_mhd
         do n=1,num_eq*ncells
            q(n,ii)=0.0
         enddo
      enddo
!
!     CALCULATE THE INITIAL RESIDUAL ERROR
!     for numrstrt = 0 this is just -(s,sbx,sby,sbz)
!
      if(numrstrt .eq. 0) then
!
      do n=1,ncells
         nn = (n-1)*num_eq +1
         do jj=1,num_eq
            q(nn+jj-1,1) = -s(n,jj) 
            dp(n,jj)=0.0
         enddo
      enddo
!
      else
!
!     input vector is (dp)
!
      call J_times_kvec(dp)
!
!     output vector is (apq)
!
      sumpq=0.0
      rnorm=0.0
      do n=1,ncells
         nn = (n-1)*num_eq +1
         do jj=1,num_eq
            q(nn+jj-1,1) = -s(n,jj) - apq(n,jj)
            sumpq=sumpq    &
             +apq(n,jj)**2
            rnorm=rnorm    &
             +s(n,jj)**2
         enddo
      enddo
      sumpq=sqrt(sumpq)
      rnorm=sqrt(rnorm)
!
      endif
!
      dot=0.0
!
      do n=1,num_eq*ncells
         dot=dot+q(n,1)*q(n,1)
      enddo
!
      dot=sqrt(dot)
      s_local(1)=dot
!
      do n=1,num_eq*ncells
         q(n,1)=q(n,1)/dot
      enddo
!
!     ****************************************************************
!
!     begin gmres
!
!     ****************************************************************
!
      do 1 iter=1,itsub_mhd
!     
!     apply preconditioner
!
      do n=1,ncells
         nn = (n-1)*num_eq + 1
         do jj=1,num_eq
            pq(n,jj)=q(nn+jj-1,iter)
         enddo
      enddo
!
!     compute J times preconditioned q 
!
      call J_times_kvec(pq)
!
!   output vector is (apq)
!
      do n=1,ncells
         nn = (n-1)*num_eq + 1
         do jj=1,num_eq
            Jq(nn+jj-1)=apq(n,jj)
         enddo
      enddo
!
!    orthogonalize:
!
      do k=1,iter
         dot=0.0
!
         do n=1,num_eq*ncells
           dot=dot+Jq(n)*q(n,k)
         enddo
!
         h(k,iter)=dot
!
         do n=1,num_eq*ncells
            Jq(n)=Jq(n)-dot*q(n,k)
         enddo
      enddo
!
!
      dot=0.0
      do n=1,num_eq*ncells
         dot=dot+Jq(n)*Jq(n)
      enddo
!
      dot=sqrt(dot)
!
      iterp1=iter+1
      h(iterp1,iter)=dot
!
!     apply previous Givens rotations to h (update QR factorization)
!
      do k=1,iter-1
         ws=g(k,1)*h(k,iter)-g(k,2)*h(k+1,iter)
         h(k+1,iter)=g(k,2)*h(k,iter)+g(k,1)*h(k+1,iter)
         h(k,iter)=ws
      enddo
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
      ws=g1*s_local(iter)-g2*s_local(iterp1)
      s_local(iterp1)=g2*s_local(iter)+g1*s_local(iterp1)
      s_local(iter)=ws
!
!     |s_local(iter+1)| is the norm of the current residual
!     check for convergence
      rnorm=abs(s_local(iterp1))
!
      if((rnorm.le.tol_GMRES).or.(iter.eq.itsub_mhd)) go to 2
!
!     normalize next q
!
      do n=1,num_eq*ncells
         q(n,iterp1)=Jq(n)/dot
      enddo
!
    1 continue
!
    2 continue
!    
!    update the solution
!    solve h*y=sbar  (sbar contains elements 1,..,iter of s)
!    store y in s
      s_local(iter)=s_local(iter)/h(iter,iter)
      do ii=iter-1,1,-1
        ws=0.0
          do jj=ii+1,iter 
          ws=ws+h(ii,jj)*s_local(jj)
          enddo
        s_local(ii)=(s_local(ii)-ws)/h(ii,ii)
      end do
!
!    compute new preconditioned linear update vector
!
       do n=1,ncells
          do jj=1,num_eq
             pq(n,jj)=0.0
          enddo
       enddo
!
      do ii=1,iter
          do n=1,ncells
             nn=(n-1)*num_eq+1
             do jj=1,num_eq
                pq(n,jj)=pq(n,jj)+s_local(ii)*q(nn+jj-1,ii)
             enddo
          enddo
      enddo
!
!   add linear iteration update vector to initial guess
!   note: inital guess is zero if numrstrt = 0
!
      do n = 1,ncells
        do jj=1,num_eq
           dp(n,jj) = dp(n,jj) + pq(n,jj)
        enddo
      enddo
!
      dxnorm=0.0
      do n=1,ncells
         do jj=1,num_eq
            dxnorm=dxnorm+dp(n,jj)**2
         enddo
      enddo
      if(rnorm .gt. tol_GMRES .and. numrstrt .lt. 4) then
       numrstrt = numrstrt + 1
       goto 1000
      endif
!
!     input vector is (dp)
!
!      call J_times_kvec(dp)
!
!     output vector is (apq)
!
!      do n=1,ncells
!         nn = (n-1)*num_eq +1
!         do jj=1,num_eq
!            q(nn,jj) = -s(n,jj) - apq(n,jj)
!         enddo
!      enddo
!
!
!      dot=0.0
!
!      do n=1,num_eq*ncells
!         dot=dot+q(n,1)*q(n,1)
!      enddo
!
!      dot=sqrt(dot)
!
       call cpu_time(Tfinish)
       do l=1,20
          if(RoutineName(l).eq.'gmres_3dmhd') CPUTime(l)=CPUTime(l)+Tfinish-Tstart
       enddo


      return
      end subroutine gmres_3dmhd




