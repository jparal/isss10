      subroutine vtxindx(iwid,jwid,kwid,ijkv)
!
!     a routine to calculate the indices of cell vertices
!
      implicit real*8 (a-h,o-z)
!
      integer :: ijkv(8)
!
      ijkv(1)=iwid
      ijkv(2)=iwid+jwid
      ijkv(3)=jwid
      ijkv(4)=0
      ijkv(5)=iwid+kwid
      ijkv(6)=iwid+jwid+kwid
      ijkv(7)=jwid+kwid
      ijkv(8)=kwid
!
      return
      end subroutine vtxindx
