      subroutine accel_3dmhd
!
!test      implicit real*8 (a-h,o-z)
!
      use vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M
      use cindex_com_M
      use cophys_com_M
      use numpar_com_M
!
      real*8 dtcntr,rcntr,dummy,factor
      integer n,nvtxtmp
!   a routine to solve the momentum equation
!
      dtcntr=dt*cntr
      rcntr=1./cntr
!
      do 1 n=1,nvtx
         ijk=ijkvtx(n)
         ul(ijk)=0.0
         vl(ijk)=0.0
         wl(ijk)=0.0
  1   continue
!      the stress is set to zero in the ghost cells
!     divpi will include contributions from the interior only
!
       call divpi(nvtx,ijkvtx,        &
               exx,exy,exz,eyy,eyz,ezz,             &
               divpix,divpiy,divpiz)
!
      dummy=0.0d0
!
!     torusbcv assembles values to form periodic vectors
!
      call torusbcv(ibp1+1,jbp1+1,kbp1+1,    &
          cdlt,sdlt,DUMMY,dz,               &
          periodic_x,periodic_y,periodic_z, &
          divpix,divpiy,divpiz)
!
      do 100 n=1,nvtx
      ijk=ijkvtx(n)
      ul(ijk)=dtcntr*divpix(ijk)
      vl(ijk)=dtcntr*divpiy(ijk)
      wl(ijk)=dtcntr*divpiz(ijk)
  100 continue
!
      if(.not.cartesian) then
!
      call axisgrad(ibp1,jbp1,ul,vl,wl)
!
      endif
!
      do 125 n=1,nvtx
      ijk=ijkvtx(n)
      ul(ijk)=ul(ijk)+dtcntr*mv(ijk)*gx
      vl(ijk)=vl(ijk)+dtcntr*mv(ijk)*gy
      wl(ijk)=wl(ijk)+dtcntr*mv(ijk)*gz
  125 continue
!
      do 150 n=1,nvtx
      ijk=ijkvtx(n)
      factor=1./(mv(ijk)+1.e-10)
      ul(ijk)=u(ijk)+ul(ijk)*factor
      vl(ijk)=v(ijk)+vl(ijk)*factor
      wl(ijk)=w(ijk)+wl(ijk)*factor
  150 continue
!
!     apply rigid, free-slip wall conditions
!     at k=2 and k=kbp2 boundaries
!
      call bc_wall(ibp1,jbp1,kbp1,iwid,jwid,kwid,   &
                  c5x,c6x,c7x,c8x,                 &
                  c5y,c6y,c7y,c8y,                 &
                  c5z,c6z,c7z,c8z,                 &
                  ul,vl,wl)           
!
!
!
      if(.not.periodic_y) then
!
!     impose no slip conditions on j=2,j=jbp2 boundaries
!
      call list(2,ibp2,2,2,2,kbp2,iwid,jwid,kwid,    &
          nvtxtmp,ijktmp2)
!
       call bc_noslip(nvtxtmp,ijktmp2,  &
         ul,vl,wl)
!
       call list(2,ibp2,jbp2,jbp2,2,kbp2,iwid,jwid,kwid,   &
         nvtxtmp,ijktmp2)
!
       call bc_noslip(nvtxtmp,ijktmp2,     &
          ul,vl,wl)
!
      endif
!
!
!  calculate the acceleration (velocity increment)
!
      do 2 n=1,nvtx
         ijk=ijkvtx(n)
         ax(ijk)=(ul(ijk)-u(ijk))*rcntr
         ay(ijk)=(vl(ijk)-v(ijk))*rcntr
         az(ijk)=(wl(ijk)-w(ijk))*rcntr
  2   continue
!
!  calculate the work
!
      do 3 n=1,nvtx
         ijk=ijkvtx(n)
         work(ijk)=(ax(ijk)**2+ay(ijk)**2+az(ijk)**2)*cntr
  3   continue
!
      return
      end
