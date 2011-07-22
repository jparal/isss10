      function sdot(ncells,a,nstridea,b,nstrideb)
c
c     a routine to compute the inner product of two vectors
c
      dimension a(*),b(*)
c
      sum=0.0
c
      do n=1,ncells
      sum=sum+a(n)*b(n)
      enddo
c
      sdot=sum
c
      return
      end
