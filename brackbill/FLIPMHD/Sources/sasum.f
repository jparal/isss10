      function sasum(ncells,a,nstride1)
c
c     form sum of absolute values
c
      dimension a(*)
c
      sum=0.0
c
      do n=1,ncells
      sum=sum+abs(a(n))
      enddo
c
      sasum=sum
c
      return
      end
