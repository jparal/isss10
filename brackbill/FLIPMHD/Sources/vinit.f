      subroutine vinit
c
      implicit real*8 (a-h,o-z)
c
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M

      include 'cindex.com'
      include 'numpar.com'
      include 'cophys.com'
c
c     a routine to prepare varibles for the next computation cycle
c     ****************************************************
c
      dummy=0.0
c
       if(rmaj.gt.0.0) then
       call torusbcv(ibp1+1,jbp1+1,kbp1+1,
     &     cdlt,sdlt,DUMMY,dz,
     &     periodic_x,periodic_y,periodic_z,
     &     umom,vmom,wmom)
       else
       call torusbc_scalar(ibp2,jbp2,kbp2,umom)
       call torusbc_scalar(ibp2,jbp2,kbp2,vmom)
       call torusbc_scalar(ibp2,jbp2,kbp2,wmom)
       end if
c
      call torusbc_scalar(ibp1+1,jbp1+1,kbp1+1,   &
          mv)
c
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,
     &     umom)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,
     &     vmom)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,
     &     wmom)
       call axisavg(ibp1,jbp1,iwid,jwid,kwid,
     &     mv)
c
c
      do 6 n=1,nvtx
      ijk=ijkvtx(n)
      factor=mv(ijk)/(mv(ijk)+1.e-10)**2
c
      u(ijk)=umom(ijk)*factor
      v(ijk)=vmom(ijk)*factor
      w(ijk)=wmom(ijk)*factor
c
      ul(ijk)=u(ijk)
      vl(ijk)=v(ijk)
      wl(ijk)=w(ijk)
c
    6 continue
c
c  apply rigid, free-slip wall conditions
c
      call bc_wall(ibp1,jbp1,kbp1,iwid,jwid,kwid,
     &             c5x,c6x,c7x,c8x,
     &             c5y,c6y,c7y,c8y,
     &             c5z,c6z,c7z,c8z,
     &             ul,vl,wl)
 write(*,*) 'vinit: returning from bc_wall'
c
c
c
      do 50 n=1,ncells
      ijk=ijkcell(n)
      rho(ijk)=mc(ijk)/(vol(ijk)+1.e-10)
      bxn(ijk)=bxn(ijk)/(vol(ijk)+1.e-10)
      byn(ijk)=byn(ijk)/(vol(ijk)+1.e-10)
      bzn(ijk)=bzn(ijk)/(vol(ijk)+1.e-10)
      bxl(ijk)=bxn(ijk)
      byl(ijk)=byn(ijk)
      bzl(ijk)=bzn(ijk)
      if(mc(ijk).ne.0.0) then
        sie(ijk)=sie1p(ijk)/(mc(ijk)+1.e-10)
      else
        sie(ijk)=0.0
      end if
  50  continue
c
      if(.not.cartesian) then
c    adjust rho on the k=2 and k=kbp1 boundaries
      fmo=1./(real(kbar)+0.5)
      do 60 i=1,ibp2
cdir$ ivdep
      do 60 j=1,jbp2
      ijkl=1+(i-1)*iwid+(j-1)*jwid
      ijkr=1+(i-1)*iwid+(j-1)*jwid+kbp1*kwid
      rho(ijkl+kwid)=rho(ijkl+kwid)
     &        *(mc(ijkl+kwid)-2.*mc(ijkl))/mc(ijkl+kwid)
      rho(ijkr-kwid)=rho(ijkr-kwid)
     &        *(mc(ijkr-kwid)+mc(ijkr)*fmo)/mc(ijkr-kwid)
   60 continue
c
      endif
c
c    initialize the stress
       do 100 k=1,kbp2
       do 100 j=1,jbp2
       do 100 i=1,ibp2
       ijk=(i-1)*iwid+(j-1)*jwid+(k-1)*kwid+1
       pixx(ijk)=0.0
       pixy(ijk)=0.0
       pixz(ijk)=0.0
       piyy(ijk)=0.0
       piyz(ijk)=0.0
       pizz(ijk)=0.0
 100   continue
c
      return
      end
