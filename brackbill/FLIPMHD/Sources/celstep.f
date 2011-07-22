      subroutine celstep(iwid, jwid, kwid, ijkstep) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!      use modify_com_M 
 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:13:36   8/20/02  
!...Switches: -yf -x1             
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iwid 
      integer , intent(in) :: jwid 
      integer , intent(in) :: kwid 
      integer , intent(out) :: ijkstep(27) 
!-----------------------------------------------
!
      ijkstep(1) = 0 
      ijkstep(2) = iwid 
      ijkstep(3) = iwid + jwid 
      ijkstep(4) = jwid 
      ijkstep(5) = (-iwid) + jwid 
      ijkstep(6) = -iwid 
      ijkstep(7) = (-iwid) - jwid 
      ijkstep(8) = -jwid 
      ijkstep(9) = iwid - jwid 
!
      ijkstep(10) = -kwid 
      ijkstep(11) = iwid - kwid 
      ijkstep(12) = iwid + jwid - kwid 
      ijkstep(13) = jwid - kwid 
      ijkstep(14) = (-iwid) + jwid - kwid 
      ijkstep(15) = (-iwid) - kwid 
      ijkstep(16) = (-iwid) - jwid - kwid 
      ijkstep(17) = (-jwid) - kwid 
      ijkstep(18) = iwid - jwid - kwid 
!
      ijkstep(19) = kwid 
      ijkstep(20) = iwid + kwid 
      ijkstep(21) = iwid + jwid + kwid 
      ijkstep(22) = jwid + kwid 
      ijkstep(23) = (-iwid) + jwid + kwid 
      ijkstep(24) = (-iwid) + kwid 
      ijkstep(25) = (-iwid) - jwid + kwid 
      ijkstep(26) = (-jwid) + kwid 
      ijkstep(27) = iwid - jwid + kwid 
!
      return  
      end subroutine celstep 
