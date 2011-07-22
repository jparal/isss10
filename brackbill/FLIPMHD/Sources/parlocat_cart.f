      subroutine parlocat_cart
!
!     a routine to calculate the logical coordinates 
!     of a particle on a uniform, rectilinear grid
!
      use vast_kind_param, ONLY : double
      use corgan_com_M, ONLY : npart
      use cindex_com_M, ONLY : ncells
      use blcom_com_M, ONLY : ijkcell, iphead,   &
        periodic_x, periodic_y, periodic_z,      &
        xl,xr,yb,yt,ze,zf,                       &
        px, py, pz, pxi, peta, pzta
      use cophys_com_M, ONLY : dx, dy, dz
!
      integer :: xcount, ycount, zcount
!
!     first, require that particle position be in bounds
!     if periodic, return to principal periodic interval
!     if not periodic, reflect particle position in normal
!
      xcount=0
      ycount=0
      zcount=0

      do n=1,ncells
        ijk=ijkcell(n)
        np=iphead(ijk)
!
   21 if(px(np).lt.xl.or.px(np).gt.xr) then
       xcount=xcount+1
!
        if(px(np).lt.xl) then
         if(periodic_x) then
           px(np)=px(np)+xr-xl
         else
           px(np)=xl+2.*(xl-px(np))
         endif
        endif
!
        if(px(np).gt.xr) then
         if(periodic_x) then
           px(np)=px(np)-(xr-xl)
         else
           px(np)=xr-2.*(px(np)-xr)
         endif
        endif
!
        go to 21
        endif
!
   22 if(py(np).lt.yb.or.py(np).gt.yt) then
       ycount=ycount+1
!
        if(py(np).lt.yb) then
         if(periodic_y) then
           py(np)=py(np)+(yt-yb)
         else
           py(np)=py(np)+2.*(yb-py(np))
         endif
        endif
!
        if(py(np).gt.yt) then
         if(periodic_y) then
           py(np)=py(np)-(yt-yb)
         else
           py(np)=yt-2.*(py(np)-yt)
         endif
        endif
        go to 22
        endif
!

!
   23  if(pz(np).lt.ze.or.pz(np).gt.zf) then
        zcount=zcount+1
        if(pz(np).lt.ze) then
         if(periodic_z) then
           pz(np)=pz(np)+(zf-ze)
         else
           pz(np)=ze+2.*(ze-pz(np))
         endif
        endif
!
        if(pz(np).gt.zf) then
         if(periodic_z) then
           pz(np)=pz(np)-(zf-ze)
         else
           pz(np)=zf-2.*(pz(np)-zf)
         endif
        endif
!
!
        go to 23
        endif
        enddo
!
!       calculate the new logical coordinates of the particle
!
        do n=1,ncells
!
         ijk=ijkcell(n)
         np=iphead(ijk)
!
         pxi(np)=2.+(px(np)-xl)/dx
         peta(np)=2.+(py(np)-yb)/dy
         pzta(np)=2.+(pz(np)-ze)/dz
        enddo
!
      return
      end subroutine parlocat_cart
