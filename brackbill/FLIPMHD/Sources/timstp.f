      subroutine timstp
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use logical_com_M
      use blcom_com_M
      use geometry_com_M
      use cindex_com_M
      use numpar_com_M
      use cophys_com_M
      implicit real*8 (a-h,o-z)
!
      real(double) :: wsvis,wsresist,wscourant,wslagrange,rdxsq
      integer ::  n
! +++
! +++ compute the new time step, dt
! +++
      dtvis=1.e20
      dtresist=1.e20
      dtcon=1.e+20
      dtiter=1.e20
!
      if(MAGNETIZED) then
!
      dtiter=dt*(.5+float(itmag)/(float(numit)+float(itmag)))
      else
      dtiter=dtgrow
      endif
!
      do 40 n=1,ncells
      ijk=ijkcell(n)
      rdxsq=       &
                 (c1x(ijk)**2+c1y(ijk)**2+c1z(ijk)**2+       &
                  c2x(ijk)**2+c2y(ijk)**2+c2z(ijk)**2+       &
                  c3x(ijk)**2+c3y(ijk)**2+c3z(ijk)**2+       &
                  c4x(ijk)**2+c4y(ijk)**2+c4z(ijk)**2+       &
                  c5x(ijk)**2+c5y(ijk)**2+c5z(ijk)**2+       &
                  c6x(ijk)**2+c6y(ijk)**2+c6z(ijk)**2+       &
                  c7x(ijk)**2+c7y(ijk)**2+c7z(ijk)**2+       &
                  c8x(ijk)**2+c8y(ijk)**2+c8z(ijk)**2)/vol(ijk)**2
      wate(ijk,1)=(mu+lam)*rdxsq
      wate(ijk,2)=(resistivity)*rdxsq
      wate(ijk,3)=(csq(ijk)+vasq(ijk))*rdxsq
      wate(ijk,4)=exx(ijk)+eyy(ijk)+ezz(ijk)

   40 continue
!
      ijk=ijkcell(1)
      wsvis=wate(ijk,1)
      wsresist=wate(ijk,2) 
      wscourant=wate(ijk,3)
      wslagrange=wate(ijk,4)
      do 45 n=2,ncells
      ijk=ijkcell(n)
      wsvis=max(wsvis,wate(ijk,1))
      wsresist=max(wsresist,wate(ijk,2))
      wscourant=max(wscourant,wate(ijk,3))
      wslagrange=min(wslagrange,wate(ijk,4))
   45 continue
!     
      dtvis=min(dtvis,0.45/(wsvis+1.e-10))
 !     dtresist=min(dtresist,0.45/(wsresist+1.e-1))
      dtresist=1.d20
      courant=sqrt(wscourant)*dt
      lagrange=-wslagrange*dt
!
      dtgrow=1.2*dt
      dtgrow=min(dtgrow,3.*dt/(courant+1.e-20))
!
      if(ncyc.eq.0) dtgrow=dt
!
!     dtpdv- constrains time step so that work done by fluid in timestep
!     cannot exceed energy in cell
!
!     if(dxmax.ne.0.) dtcon=stabl*dtgrow/(dxmax)
      dt=min(dtgrow,dtcon,dtmax,dtvis,dtresist,dtpdv,dtiter)
      if(dt.eq.dtiter) iddt='i'
      if(dt.eq.dtgrow) iddt='g'
      if(dt.eq.dtmax) iddt='m'
      if(dt.eq.dtpdv) iddt='p'
      if(dt.eq.dtvis) iddt='v'
      if(dt.eq.dtresist) iddt='r'
      if(ncyc.eq.1) dtmin=dt*1.e-10
      if(dt.gt.dtmin) return
      write(6,100)  dt,ncyc,iddt
      write(59,100) dt,ncyc,iddt
      twfin=0.
      return
  100 format(4h dt=,1pe12.5,9h at cycle,i5,11h, cause is ,a1)
  200 format ( 7h timstp ,(9i7))
      end subroutine timstp
