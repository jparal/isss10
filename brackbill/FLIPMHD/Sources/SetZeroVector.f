       subroutine SetZeroVector(nlist,ijklist,varx,vary,varz)
!
      integer ijklist(*)
      real varx(*),vary(*),varz(*)
!
!      set an array to zero
!
      do n=1,nlist
       varx(ijklist(n))=0.0
       vary(ijklist(n))=0.0
       varz(ijklist(n))=0.0
      enddo
!
       return
       end subroutine SetZeroVector
