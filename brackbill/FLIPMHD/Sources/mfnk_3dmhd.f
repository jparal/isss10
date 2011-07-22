      subroutine mfnk_3dmhd(    &
          error,eps,                           &
          numit,itmax)
!
!***BEGIN PROLOGUE MFNK_3dMHD
!***DATE WRITTEN 1/99
!*** constructed from itercg_i3dmhd
!***KEYWORD Matrix-Free Newton Krylov ITERATION, GMRES
!***AUTHORS D.A. Knoll and J.U. Brackbill (LOS ALAMOS NATIONAL LABORATORY)
!***PURPOSE
!     To solve the implicit MHD equations, including the ,adiabatic pressure, changes
!     Faraday's, and momentum equations, J * dx = - residual 
!     where J depends upon ro, Bx, By, Bz
!     residu_3dmhd computes residual error
!     J_times-kvec computes J * dx to begin each Newton iteration
!     gmres_3dmhd solves linear system to specified tolerance each Newton iteration
!
!     num_eq is the number of independent variables per cell
!     xnorm is the L1 norm of the solution vector
!     ynorm is the L1 norm of the RHS
!     dxnorm is the L2 norm of the newton update
!     rnorm is the L2 norm of the nonlinear residual
!     eps is the error tolerance for convergence
!***ROUTINES CALLED SENSITIV LUDECOMP RESIDUM PRECOND
!***END PROLOGUE ITERCG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       USE vast_kind_param, ONLY:  double
       use corgan_com_M
       use cindex_com_M, ONLY : ncells
       use blcom_com_M, ONLY : &
          bxl,byl,bzl,pl,    &
          ijktmp2,           &
          p,csq, ijkcell
       USE nkmhd_com_M
       use nk_com_M, ONLY : solution, s, dp
       use Timing
!
      implicit none

      integer ::  numit,itmax,  &
          l,m,n,ijk, jj
!
      real(double) :: divu_min,ws,zero,      &
          error,eps,Tstart,Tfinish
!
!
      real(double) ::      &
          rnorm0,  rfac, rdamp, tol_GMRES 
!
      call cpu_time(Tstart)
      zero=0.d0
!
! ********************************************************************
!
!     to make indexing compatible with gmres, which is used for other equations,
!     store unknowns in solution array
      do n=1,ncells
         ijk=ijkcell(n)
         solution(n,1)=pl(ijk)
         solution(n,2)=bxl(ijk)
         solution(n,3)=byl(ijk)
         solution(n,4)=bzl(ijk)
      enddo
      
      if(ncells.eq.0) stop
      numit=0
!
!     ********************************************************************
!
!     BEGIN Newton ITERATION
!
!     *******************************************************************
!
   25 continue
      numit=numit+1
!
!     residu_3dmhd evaluates the residual error in the solution
!     with current values of the magnetic field and pressure
!
!
!     residu_3dmhd returns s,sbx,sby,sbz ... residual errors in
!     continuity and Faraday's equations
!
      if(numit .eq. 1)then
!
         call residu_3dmhd(   &
             solution,s)
!
!     compute L2 norm of residue
!
         rnorm = 0.0e0
!
         do n=1,num_eq
            do m=1,ncells
!
            rnorm= rnorm    &
             + s(m,n)**2
!
           enddo
         enddo
!
         rnorm = sqrt(rnorm)
           rnorm0 = rnorm
!
!
           if(rnorm.lt.eps*ynorm) then
              numit=0
              return
           endif
!
      endif
!
!     call preconditioned GMRES
!     returns dp ... Newton updates
!
!     GMRES tolerance will be 0.01 * rnorm (inexact newton method)
!
!
      tol_GMRES = error*rnorm
!
   
      call gmres_3dmhd(tol_GMRES,ncells)
!
!     compute Newton damping scalar, rdamp, based on a maximum
!     allowable change in pressure, pl, of 50 %
!
      rdamp = 1.0e0
!
      do m = 1,ncells
       ijk = ijkcell(m)

      if((solution(m,1)+dp(m,1)) .lt. 0.5e0*solution(m,1) ) then
        rfac = -(0.5e0)*solution(m,1)/dp(m,1)
        if(rfac .lt. rdamp) rdamp = rfac
       elseif((solution(m,1)+dp(m,1)) .gt. 1.5e0*solution(m,1) ) then
        rfac = (0.5e0)*solution(m,1)/dp(m,1)
        if(rfac .lt. rdamp) rdamp = rfac
       endif

      end do
      rdamp=1.0
!
!     increment the iteration variables 
!
      do m=1,ncells
         do jj=1,num_eq
            solution(m,jj)=solution(m,jj)+ rdamp*dp(m,jj)
         enddo
      enddo
!
!     calculate residual with  updated solution
!
      call residu_3dmhd(  &
          solution,s)

!     compute L2 norm of residue
!
      rnorm = 0.0e0
      dxnorm = 0.0e0
!
         do n=1,num_eq
           do m=1,ncells
              ijk = ijkcell(m)
              rnorm= rnorm + s(m,n)**2
              dxnorm= dxnorm + (dp(m,n)/(solution(m,n)+1.d-10))**2
           enddo
         enddo
!
      rnorm = sqrt(rnorm)
      dxnorm = sqrt(dxnorm)/(4.0*ncells)
!
!
!     when both nonlinear residual and newton update are
!     both sufficiently small,
!     iteration has converged
!
!      ratio=rnorm/(ynorm+1.e-50)
!      xratio=dxnorm/(xnorm+1.e-50)
!
         if(rnorm/rnorm0.le.eps) go to 1025
!
      if(numit.gt.itmax) go to 9925
      go to 25
 9925 continue
      write(6,9926) numit,rnorm
!
!    should we cut time step and try again
!
 9926 format(" mfnk failed to converge, numit=",i4,"rnorm=",1pe10.2)
      return
 1025 continue
!
!     transfer of solution back to cell arrays
!     is already done elsewhere
!      do n=1,ncells
!         ijk=ijkcell(n)
!         pl(ijk)=solution(n,1)
!         bxl(ijk)=solution(n,2)
!         byl(ijk)=solution(n,3)
!         bzl(ijk)=solution(n,4)
!      enddo
      call cpu_time(Tfinish)
      do l=1,20
         if(RoutineName(l).eq.'mfnk') CPUTime(l)=CPUTime(l)+Tfinish-Tstart
      enddo
      return
      end
