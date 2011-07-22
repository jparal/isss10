      subroutine bc_wall(ibp1,jbp1,kbp1,iwid,jwid,kwid,   &
                        c5x,c6x,c7x,c8x,                 &
                        c5y,c6y,c7y,c8y,                 &
                        c5z,c6z,c7z,c8z,                 &
                        ul,vl,wl)
!
      USE vast_kind_param, ONLY:  double
      implicit real*8 (a-h,o-z)
!
!     a routine to impose rigid, free-slip wall conditions
!
       real(double) ::                  &
        c5x(*),c6x(*),c7x(*),c8x(*),   &
        c5y(*),c6y(*),c7y(*),c8y(*),   &
        c5z(*),c6z(*),c7z(*),c8z(*),   &
        ul(*),vl(*),wl(*)
!
      real*8 normx,normy,normz
!
      ktop=kbp1+1
      kbot=2
      do 1 j=2,jbp1+1
      do 1 i=2,ibp1+1
          ijk=1+(i-1)*iwid+(j-1)*jwid+(ktop-1)*kwid
!          normx=c8x(ijk-kwid)+c5x(ijk-iwid-kwid)
!     &         +c6x(ijk-iwid-jwid-kwid)+c7x(ijk-jwid-kwid)
!          normy=c8y(ijk-kwid)+c5y(ijk-iwid-kwid)
!     &         +c6y(ijk-iwid-jwid-kwid)+c7y(ijk-jwid-kwid)
!          normz=c8z(ijk-kwid)+c5z(ijk-iwid-kwid)
!     &         +c6z(ijk-iwid-jwid-kwid)+c7z(ijk-jwid-kwid)
!
!          udotn=(ul(ijk)*normx+vl(ijk)*normy+wl(ijk)*normz)
!     &         /(normx**2+normy**2+normz**2+1.e-10)
!
!          ul(ijk)=ul(ijk)-normx*udotn
!          vl(ijk)=vl(ijk)-normy*udotn
!          wl(ijk)=wl(ijk)-normz*udotn
      wl(ijk)=0.0
!
      ijk=1+(i-1)*iwid+(j-1)*jwid+(kbot-1)*kwid
!
      wl(ijk)=0.0
  1   continue
!
      return
      end
 
