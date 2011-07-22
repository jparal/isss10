      subroutine mshset
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M
      use geometry_com_M
      use cindex_com_M
      use numpar_com_M
      use cophys_com_M
      implicit none
!

      integer :: n, kvacx
      real(double) :: del1,del2,del3,del4,del5,del6,del7,    &
          delt,dtheta,dzstr
!
!
!
!     package for a three-dimensional torus
!     theta increases with i
!     phiang increases with j
!     r increases with k
!     rwall is the minor radius
!     rmaj is the major radius
!     scyllac package
!
!     define rotation, iota
!
!delete      rfkbp2=1./real(kbp2)
!
!
!     define rotation matrix for toroidal geometry
!
      write(*,*)'calling rotate'
      call rotate(ibar,jbar,kbar,    &
          rmaj,dz,strait,iota,istep,    &
          delt,dtheta,dzstr,dphi,    &
          cdlt,sdlt,cdlhf,sdlhf,    &
          cdph,sdph,cdphhf,sdphhf)
!
      write(*,*) 'returning from rotate'
!
      write(*,*)'calling gridinit'
      call gridinit(ibar,jbar,kbar,iwid,jwid,kwid,    &
          delt,dphi,dtheta,dx,dy,dz,rwall,rmaj,dzstr,istep,    &
          del1,del2,del3,del4,del5,del6,del7,    &
          x,y,z,cartesian,xl,xr,yb,yt,ze,zf)
      write(*,*)'returning from gridinit'
!
!
         kvacx= kvac
         kplas= kbar
         rvac = rwall
         kvac= kplas+2
         kvp = kvac+1
         kvm = kvac-1
!
!
!ll   3.  calculate geometric coefficients:
       write(*,*) 'calling geom'
       call geom(ncells,ijkcell,nvtx,ijkvtx,   &
          nvtxkm,   &
          x,y,z,   &
          vol,vvol)
       write(*,*) 'returning from geom,dt,t=',dt,t
   do n=1,nvtx
      ijk=ijkvtx(n)
      if(vvol(ijk).eq.0.0) then
        write(*,*) 'from geom: ijk, vvol=',ijk, vvol(ijk)
        stop
     endif
    enddo

!      call metric(x,y,z,tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz)
!      call metric(ncells,ijkcell,iwid,jwid,kwid,
!     &     x,y,z,
!     &     tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz,sie)
!
      write(*,*) 'mshset: call metricc'
      call metricc(ncells,ijkcell,   &
          vol,   &
          tsix,tsiy,tsiz,etax,etay,etaz,nux,nuy,nuz)
     write(*,*) 'mshset: return from metricc,dt,t=',dt,t
!
!     *********************************
!     baal3 sets boundary conditions for cell-centered variables
!     ************************************
!delete      call bcc
!
      return
      end subroutine mshset
