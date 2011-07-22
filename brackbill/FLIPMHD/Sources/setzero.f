       subroutine setzero(nlist,ijklist,var)
!
      integer ijklist(*)
      real var(*)
!
!      set an array to zero
!
      do n=1,nlist
       var(ijklist(n))=0.0
      enddo
!
       return
       end subroutine setzero
