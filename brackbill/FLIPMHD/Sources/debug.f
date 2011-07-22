       subroutine debug(ncells,ijkcell,iwid,jwid,kwid,   &
         cdlt,sdlt,strait,dz,   &
         ibp1,jbp1,kbp1,nvtx,ijkvtx,xc,yc,zc,   &
         x,y,z,vol,vvol,gradcx,gradcy,gradcz,divu)
!
      use vast_kind_param, only:  double
      use geometry_com_M
!
      integer ::      &
          iwid,jwid,kwid,ibp1,jbp1,kbp1,      &
          ncells,ijkcell(*),nvtx,ijkvtx(*)
      real(double) :: divu(*),xc(*),yc(*),zc(*),vvol(*),      &
          x(*),y(*),z(*),vol(*),gradcx(*),gradcy(*),gradcz(*)
!
      call divc(ncells,ijkcell,iwid,jwid,kwid,      &
          vol,      &
          x,y,z,divu)
!
      call gradc(ncells,ijkcell,      &
          x,gradcx,gradcy,gradcz)
!
      do n=1,ncells
      ijk=ijkcell(n)
      xc(ijk)=0.125*(x(ijk)+x(ijk+iwid)+x(ijk+iwid+jwid)+x(ijk+jwid)      &
           +x(ijk+kwid)+x(ijk+iwid+kwid)      &
           +x(ijk+iwid+jwid+kwid)+x(ijk+jwid+kwid))
      yc(ijk)=0.125*(y(ijk)+x(ijk+iwid)+y(ijk+iwid+jwid)+y(ijk+jwid)      &
           +y(ijk+kwid)+y(ijk+iwid+kwid)      &
           +y(ijk+iwid+jwid+kwid)+y(ijk+jwid+kwid))
      zc(ijk)=0.125*(z(ijk)+z(ijk+iwid)+z(ijk+iwid+jwid)+z(ijk+jwid)      &
           +z(ijk+kwid)+z(ijk+iwid+kwid)      &
           +z(ijk+iwid+jwid+kwid)+z(ijk+jwid+kwid))
      end do
!
!     call torusbc(ibp2,jbp2,kbp2,              &
!         strait,                               &
!         xc,yc,zc)
!
      write(*,*) 'debug:  calling divv,nvtx,iwid,jwid,kwid=',nvtx,iwid,jwid,kwid
      call divv(      &
          xc,yc,zc,divu)
!
      call torusbc_scalar(ibp1+1,jbp1+1,kbp1+1,      &
         periodic_x,periodic_y,periodic_z,      &
         divu)
!
      call axisavg(ibp1,jbp1,iwid,jwid,kwid,      &
          divu)
!
      call bcphi(ibp1,jbp1,kbp1,iwid,jwid,kwid,yc)
!
      call gradf(nvtx,ijkvtx,      &
          yc,gradcx,gradcy,gradcz)
!
      call axisgrad(ibp1,jbp1,gradcx,gradcy,gradcz)
!
      do 10 n=1,nvtx
      ijk=ijkvtx(n)
      gradcx(ijk)=gradcx(ijk)/vvol(ijk)
      gradcy(ijk)=gradcy(ijk)/vvol(ijk)
      gradcz(ijk)=gradcz(ijk)/vvol(ijk)
      divu(ijk)=divu(ijk)/vvol(ijk)
  10  continue
!
      return
      end subroutine debug
