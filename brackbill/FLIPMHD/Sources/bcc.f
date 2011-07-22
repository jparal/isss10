      subroutine bcc
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M
      use cindex_com_M
      use numpar_com_M
      use cophys_com_M
!
      implicit none
      integer :: ijkr,ipjkr,ijkb,ijpkb

!
!     initialize quantities in dummy cells
!
!     set quantities in k=kbp2 plane
!
      do 50 j=2,jbp1
      do 50 i=2,ibp1
      ijk=1+(kbp1-1)*kwid+(j-1)*jwid+(i-1)*iwid
      ijkp=ijk+kwid
!
      tsix(ijkp)=tsix(ijk)
      etax(ijkp)=etax(ijk)
      nux(ijkp)=nux(ijk)
!
      tsiy(ijkp)=tsiy(ijk)
      etay(ijkp)=etay(ijk)
      nuy(ijkp)=nuy(ijk)
!
      tsiz(ijkp)=tsiz(ijk)
      etaz(ijkp)=etaz(ijk)
      nuz(ijkp)=nuz(ijk)
!
   50 continue
!     set quantities in i=1 column and in i=ibp2 column
      do 150 k=2,kbp2
         ijk= 1+(k-1)*kwid
      do 110 j=1,jbp2
         ijkr=ijk+iper
         ipjk=ijk+1
         ipjkr=ipjk+iper
         tsix(ijk)=tsix(ijkr)
         tsiy(ijk)=tsiy(ijkr)
         tsiz(ijk)=tsiz(ijkr)
         etax(ijk)=etax(ijkr)
         etay(ijk)=etay(ijkr)
         etaz(ijk)=etaz(ijkr)
         nux(ijk)=nux(ijkr)
         nuy(ijk)=nuy(ijkr)
         nuz(ijk)=nuz(ijkr)
         tsix(ipjkr)=tsix(ipjk)
         tsiy(ipjkr)=tsiy(ipjk)
         tsiz(ipjkr)=tsiz(ipjk)
         etax(ipjkr)=etax(ipjk)
         etay(ipjkr)=etay(ipjk)
         etaz(ipjkr)=etaz(ipjk)
         nux(ipjkr)=nux(ipjk)
         nuy(ipjkr)=nuy(ipjk)
         nuz(ipjkr)=nuz(ipjk)
         ijk=ijk+jwid
  110 continue
!     set quantities in the j=1 and j=jbp2 rows
         ijk= 1+(k-1)*kwid
      do 120 i=1,ibp2
         ijkb=ijk+jper
         ijpk=ijk+jwid
         ijpkb=ijpk+jper
         tsix(ijk) =  tsix(ijkb)*cdlt+tsiy(ijkb)*sdlt
         etax(ijk) =  etax(ijkb)*cdlt+etay(ijkb)*sdlt
         nux(ijk) =  nux(ijkb)*cdlt+nuy(ijkb)*sdlt
         tsiy(ijk) = -tsix(ijkb)*sdlt+tsiy(ijkb)*cdlt
         etay(ijk) = -etax(ijkb)*sdlt+etay(ijkb)*cdlt
         nuy(ijk) = -nux(ijkb)*sdlt+nuy(ijkb)*cdlt
         tsiz(ijk)=tsiz(ijkb)
         etaz(ijk)=etaz(ijkb)
         nuz(ijk)=nuz(ijkb)
         tsix(ijpkb) =  tsix(ijpk)*cdlt-tsiy(ijpk)*sdlt
         etax(ijpkb) =  etax(ijpk)*cdlt-etay(ijpk)*sdlt
         nux(ijpkb) =  nux(ijpk)*cdlt-nuy(ijpk)*sdlt
         tsiy(ijpkb) =  tsix(ijpk)*sdlt+tsiy(ijpk)*cdlt
         etay(ijpkb) =  etax(ijpk)*sdlt+etay(ijpk)*cdlt
         nuy(ijpkb) =  nux(ijpk)*sdlt+nuy(ijpk)*cdlt
         tsiz(ijpkb)=tsiz(ijpk)
         etaz(ijpkb)=etaz(ijpk)
         nuz(ijpkb)=nuz(ijpk)
         ijk=ijk+1
  120 continue
  150 continue
      return
      end
