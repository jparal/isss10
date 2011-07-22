      subroutine celdex(i,j,k,iwid,jwid,kwid,ijkc)
 !-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, ONLY : double

!...Translated by Pacific-Sierra Research 77to90  4.3E  14:13:36   8/20/02
!...Switches: -yf -x1
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i, j, k
      integer , intent(in) :: iwid
      integer , intent(in) :: jwid
      integer , intent(in) :: kwid
      integer , intent(out) :: ijkc(27)
!-----------------------------------------------
!
      integer :: ijk

!
      ijk=1+iwid*(i-1)+jwid*(j-1)+kwid*(k-1)
!
      ijkc(1)=ijk
      ijkc(2)=ijk+iwid
      ijkc(3)=ijk+iwid+jwid
      ijkc(4)=ijk+jwid
      ijkc(5)=ijk-iwid+jwid
      ijkc(6)=ijk-iwid
      ijkc(7)=ijk-iwid-jwid
      ijkc(8)=ijk-jwid
      ijkc(9)=ijk+iwid-jwid
!
      ijkc(10)=ijk-kwid
      ijkc(11)=ijk+iwid-kwid
      ijkc(12)=ijk+iwid+jwid-kwid
      ijkc(13)=ijk+jwid-kwid
      ijkc(14)=ijk-iwid+jwid-kwid
      ijkc(15)=ijk-iwid-kwid
      ijkc(16)=ijk-iwid-jwid-kwid
      ijkc(17)=ijk-jwid-kwid
      ijkc(18)=ijk+iwid-jwid-kwid
!
      ijkc(19)=ijk+kwid
      ijkc(20)=ijk+iwid+kwid
      ijkc(21)=ijk+iwid+jwid+kwid
      ijkc(22)=ijk+jwid+kwid
      ijkc(23)=ijk-iwid+jwid+kwid
      ijkc(24)=ijk-iwid+kwid
      ijkc(25)=ijk-iwid-jwid+kwid
      ijkc(26)=ijk-jwid+kwid
      ijkc(27)=ijk+iwid-jwid+kwid
!
      return
      end subroutine celdex
