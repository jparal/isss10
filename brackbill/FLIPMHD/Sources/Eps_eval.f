      subroutine Eps_eval(GMit,u, du, eps,   &
         nstart,nfinish,ijklist)

!     evaluates ueps = u+eps.du

      use vast_kind_param, ONLY: double
      use corgan_com_M, ONLY : itdim
      integer :: GMit
      real(double) :: u(*),du(itdim,*)
      real(double) :: eps, epsfact,sum,sum2
      integer :: n,ij,nstart,nfinish,ijklist(*)
      

      tolerance = 0.00001
      epsfact=1.d1

      sum = 0.0
      sum2 = 0.0
      do n=nstart,nfinish
         ijk=ijklist(n)
         sum = sum + abs(u(ijk))
         sum2 = sum2 + abs(du(ijk,GMit))
      end do

!     if velocity is small, set eps based on maximum shear wave speed..
      if (sum.lt.tolerance) sum = tolerance

!     nfinish-nstart is 3 times as large when each entry doesn't correspond to 3 components
      eps = epsfact*sum/(sum2*(nfinish-nstart+1))
      return
      end

