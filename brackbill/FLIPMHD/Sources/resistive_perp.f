      subroutine resistive_diff
c
c     a routine to calculate diffusion of the magnetic field due to
c     resistance
c
      use blcom_com_M
      use geometry_com_M
      include 'corgan.com'
      include 'cindex.com'
      include 'numpar.com'
      include 'cophys.com'
c
      real*8 exc,eyc,ezc,jxc,jyc,jzc
c
c     calculate the electric field due to resistance
c     to current flow
c
      call curlv(nvtx,ijkvtx,
     &     bxl,byl,bzl,jx,jy,jz)
c
c     impose periodic boundary conditions
c
      zero=0.0
c
      call bc_current(ibp1+1,jbp1+1,kbp1+1,
        zero,                                &
     &     jx,jy,jz)
c
c     calculate the values of the magnetic field
c     at the vertices of the mesh, and store in gradx,grady,gradz
c
      call b_vtx(ncells,ijkcell,iwid,jwid,kwid,
     &    nvtxs,ijkvtx,
     &    bxl,byl,bzl,
     &    vol,upv1,
     &    bxv,byv,bzv)
c
c     calculate the electric field
c
      do n=1,nvtx
c
      ijk=ijkvtx(n)
c
      rvvol=1./vvol(ijk)
      jx(ijk)=jx(ijk)*rvvol
      jy(ijk)=jy(ijk)*rvvol
      jz(ijk)=jz(ijk)*rvvol
c
      jdotb=bxv(ijk)*jx(ijk)+byv(ijk)*jy(ijk)+bzv(ijk)*jz(ijk)
      bsq=bxv(ijk)**2+byv(ijk)**2+bzv(ijk)**2+1.d-20
c
c     the parallel current will not contribute to resistive diffusion
c
      ex(ijk)=resistivity*(jx(ijk)-bxv(ijk)*jdotb/bsq)
      ey(ijk)=resistivity*(jy(ijk)-byv(ijk)*jdotb/bsq)
      ez(ijk)=resistivity*(jz(ijk)-bzv(ijk)*jdotb/bsq)
c
cjub      ex(ijk)=resistivity*jx(ijk)
cjub      ey(ijk)=resistivity*jy(ijk)
cjub      ez(ijk)=resistivity*jz(ijk)
      enddo
   
c
c     calculate the curl of the electric field
c
      call curlc(ncells,ijkcell,
     &     vol,
     &     ex,ey,ez,gradx,grady,gradz)
c
c
      do n=1,ncells
      ijk=ijkcell(n)
c
c     calculate the energy dissipated by resistive diffusion
c
      exc=0.125*(ex(ijk)+ex(ijk+iwid)
     &     +ex(ijk+iwid+jwid)+ex(ijk+jwid)
     &     +ex(ijk+kwid)+ex(ijk+iwid+kwid)
     &     +ex(ijk+iwid+jwid+kwid)+ex(ijk+jwid+kwid))
c
      eyc=0.125*(ey(ijk)+ey(ijk+iwid)
     &     +ey(ijk+iwid+jwid)+ey(ijk+jwid)
     &     +ey(ijk+kwid)+ey(ijk+iwid+kwid)
     &     +ey(ijk+iwid+jwid+kwid)+ey(ijk+jwid+kwid))
c
      ezc=0.125*(ez(ijk)+ez(ijk+iwid)
     &     +ez(ijk+iwid+jwid)+ez(ijk+jwid)
     &     +ez(ijk+kwid)+ez(ijk+iwid+kwid)
     &     +ez(ijk+iwid+jwid+kwid)+ez(ijk+jwid+kwid))
c
      jxc=0.125*(jx(ijk)+jx(ijk+iwid)
     &     +jx(ijk+iwid+jwid)+jx(ijk+jwid)
     &     +jx(ijk+kwid)+jx(ijk+iwid+kwid)
     &     +jx(ijk+iwid+jwid+kwid)+jx(ijk+jwid+kwid))
c
      jyc=0.125*(jy(ijk)+jy(ijk+iwid)
     &     +jy(ijk+iwid+jwid)+jy(ijk+jwid)
     &     +jy(ijk+kwid)+jy(ijk+iwid+kwid)
     &     +jy(ijk+iwid+jwid+kwid)+jy(ijk+jwid+kwid))
c
      jzc=0.125*(jz(ijk)+jz(ijk+iwid)
     &     +jz(ijk+iwid+jwid)+jz(ijk+jwid)
     &     +jz(ijk+kwid)+jz(ijk+iwid+kwid)
     &     +jz(ijk+iwid+jwid+kwid)+jz(ijk+jwid+kwid))
c
      Ohmic_heating(ijk)=(exc*jxc+eyc*jyc+ezc*jzc)
     &   *dt*vol(ijk)

      bxl(ijk)=bxl(ijk)-gradx(ijk)*dt
      byl(ijk)=byl(ijk)-grady(ijk)*dt
      bzl(ijk)=bzl(ijk)-gradz(ijk)*dt
c
      enddo
c
      return
      end
