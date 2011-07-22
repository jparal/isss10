     subroutine L2NormF(F, Fnorm, &
          nstart,nstop,ijklist)

!     calculates the magnitude of a given stored function F = (Fu, Fv)
      use vast_kind_param, ONLY:  double

      integer :: n,nstart,nstop
      integer:: ijklist(*)
      real(double), intent(in) :: F(*)
      real(double) :: Fnorm,sum

      sum = 0.0d0
      do n=nstart,nstop
         ijk=ijklist(n)
         sum = sum + F(ijk)**2
      end do

      Fnorm = dsqrt(sum) 

      return
      end
