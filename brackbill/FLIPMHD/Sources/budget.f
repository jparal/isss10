      subroutine budget
!*********************************************************
!
      USE vast_kind_param, ONLY:  double
      use corgan_com_M
      use blcom_com_M
      use cindex_com_M
      use cophys_com_M
      use numpar_com_M
      use Timing
!
      implicit none
      real(double) :: dummy,cjump,totm,tpmox,tpmoy,tpmoz,sixteenth
      real(double) :: voll, Tstart, Tfinish, zero
      integer :: l,ijkl,ijkr,newcell,newcell_nxt,np,n,   &
          istop,jstop,kr
!
      call cpu_time(Tstart)
      nptotl=0
      toti=0.0
      totmb=0.0d0
      totbb=0.0d0
      totmm=0.0d0
      totk=0.0d0
      totkz=0.0d0
      tote=0.0d0
      totm=0.0d0
      tpmox=0.0
      tpmoy=0.0
      tpmoz=0.0
      do 25 kr=0,nreg
        numtot(kr)=0
   25 continue
!
      zero=0.0d0
!
!    calculate time advanced strain:
      call strain(ncells,ijkcell,iwid,jwid,kwid,   &
          exx,exy,exz,eyy,eyz,ezz,    &
          ul,vl,wl,vol)
!
      do n=1,ncells
      ijk=ijkcell(n)
      divu(ijk)=exx(ijk)+eyy(ijk)+ezz(ijk)
      end do
!
!
      reconnected_flux(ncyc)=0.d0
!
      istop=ibp2
      jstop=jbp2
      if(periodic_x) istop=ibp1
      if(periodic_y) jstop=jbp1
!
        newcell=0
        do 100 n=1,ncells
        ijk=ijkcell(n)
!
!
!     compute total magnetic energy
!
      totmb=totmb     &
          +0.5d0*bmagx(ijk)*bxn(ijk)*vol(ijk)    &
          +0.5d0*bmagy(ijk)*byn(ijk)*vol(ijk)    &
          +0.5d0*bmagz(ijk)*bzn(ijk)*vol(ijk)
!
      totbb=totbb     &
          +0.5d0*(bxn(ijk)*bxn(ijk)   &
                 +byn(ijk)*byn(ijk)   &
                 +bzn(ijk)*bzn(ijk))*vol(ijk)
!
      toti=toti+rho(ijk)*sie(ijk)*vol(ijk)
!
      sixteenth=1./16.
!
      totk=totk+sixteenth*(mv(ijk)*(u(ijk)**2+v(ijk)**2+w(ijk)**2)   &
     + mv(ijk+iwid)*(u(ijk+iwid)**2+v(ijk+iwid)**2+w(ijk+iwid)**2)   &
     + mv(ijk+iwid+jwid)*(u(ijk+iwid+jwid)**2+v(ijk+iwid+jwid)**2+w(ijk+iwid+jwid)**2)   &
     + mv(ijk+jwid)*(u(ijk+jwid)**2+v(ijk+jwid)**2+w(ijk+jwid)**2)   &
     + mv(ijk+kwid)*(u(ijk+kwid)**2+v(ijk+kwid)**2+w(ijk+kwid)**2)   &
     + mv(ijk+iwid+kwid)*(u(ijk+iwid+kwid)**2+v(ijk+iwid+kwid)**2+w(ijk+iwid+kwid)**2)   &
     + mv(ijk+iwid+jwid+kwid)       &
         *(u(ijk+iwid+jwid+kwid)**2+v(ijk+iwid+jwid+kwid)**2+w(ijk+iwid+jwid+kwid)**2)   &
     + mv(ijk+jwid+kwid)*(u(ijk+jwid+kwid)**2+v(ijk+jwid+kwid)**2+w(ijk+jwid+kwid)**2)) 
!
        np=iphead(ijk)
        if(np.ne.0) then
        newcell=newcell+1
        ijkctmp(newcell)=ijk
        iphd2(ijk)=np
        endif
  100   continue
!
       diffbb=totbb-totmb
!
    1 continue
        if(newcell.eq.0) go to 900
!
      do 600 n=1,newcell
      ijk=ijkctmp(n)
!
      np=iphd2(ijk)
      if(np.ne.0) then
      nptotl=nptotl+1
      kr=ico(np)
!
!     calculate the particle energy
!      toti=toti+ep(np)
!      totk=totk+0.5*mass(np)*(up(np)**2+vp(np)**2+wp(np)**2)
      totkz=totkz+0.5*mass(np)*wp(np)**2
      tpmox=tpmox+mass(np)*up(np)
      tpmoy=tpmoy+mass(np)*vp(np)
      tpmoz=tpmoz+mass(np)*wp(np)
!
      totmm=totmm     &
          +0.5d0*(mupx(np)**2+mupy(np)**2+mupz(np)**2)/vol(ijk)
!     save total number in each region	
      numtot(kr)=numtot(kr)+1
!
      totm=totm+mass(np)
      endif
!
  600 continue
!
      newcell_nxt=0
      do 500 n=1,newcell
      ijk=ijkctmp(n)
      np=link(iphd2(ijk))
      if(np.gt.0) then
      newcell_nxt=newcell_nxt+1
      iphd2(ijk)=np
      ijkctmp(newcell_nxt)=ijk
      endif
  500 continue
!
      newcell=newcell_nxt
!
      go to 1
!
  900 continue
!
      ebnrg(nh)=totmb
      eknrg(nh)=totk
      efnrg(nh)=toti
      tote=toti+totk+totmb

!
!    calculate the reconnected flux
!
      call magnetopause(ncells,ijkcell,iwid,jwid,kwid,   &
          vol,     &
          color,gradx,grady,gradz,     &
          bxn,byn,bzn,cjump,     &
          flux,surface_area)
!
      reconnected_flux(ncyc)=flux
!
!    calculate pdv, viscous dissipation 
!    add resistive heating
!
      do 1000 n=1,ncells
      ijk=ijkcell(n)
      pdv(ijk)=(-p(ijk)*divu(ijk)      &
             +exx(ijk)*pixx(ijk)+eyy(ijk)*piyy(ijk)+ezz(ijk)*pizz(ijk)    &
      +2.0*(exy(ijk)*pixy(ijk)+exz(ijk)*pixz(ijk)+eyz(ijk)*piyz(ijk)))    &
         *vol(ijk)*dt/(mc(ijk)+1.e-10)                                    &
          +Ohmic_heating(ijk)/(mc(ijk)+1.e-10)
!
      voll=vol(ijk)*(1.+divu(ijk)*dt)
!
      dbvx(ijk)=(bxl(ijk)*voll-bxn(ijk)*vol(ijk)-divuphix(ijk)*dt)    &
         /(mc(ijk)+1.d-20)
!
      dbvy(ijk)=(byl(ijk)*voll-byn(ijk)*vol(ijk)-divuphiy(ijk)*dt)    &
         /(mc(ijk)+1.d-20)
!
      dbvx(ijk)=(bzl(ijk)*voll-bzn(ijk)*vol(ijk)-divuphiz(ijk)*dt)    &
         /(mc(ijk)+1.d-20)
!!!!         dbvx(ijk)=bmagx(ijk)*(bxl(ijk)-bxn(ijk))*vol(ijk)/(bxn(ijk)*mc(ijk)+1.d-20)
!!!!         dbvy(ijk)=bmagy(ijk)*(byl(ijk)-byn(ijk))*vol(ijk)/(byn(ijk)*mc(ijk)+1.d-20)
!!!!         dbvz(ijk)=bmagz(ijk)*(bzl(ijk)-bzn(ijk))*vol(ijk)/(bzn(ijk)*mc(ijk)+1.d-20)
!
 1000 continue
!
!    periodic boundary conditions in i and j:
!
      dummy=0.0
!
!     no accumulation in this routine, so can call with pdv for
!     all three components
!     must change if grid surfaces are not parallel
!
      call torusbc(ibp2,jbp2,kbp2,    &
          dummy,                      &
          pdv,pdv,pdv)
!
!     reflect in k
!     omf=-1./(real(kbar)+0.5)
      do i=1,ibp2
        do  j=1,jbp2
        ijkl=1+(i-1)*iwid+(j-1)*jwid
        ijkr=1+(i-1)*iwid+(j-1)*jwid+kbp1*kwid
        pdv(ijkl)=pdv(ijkl+kwid)
        pdv(ijkr)=pdv(ijkr-kwid)
        enddo
      enddo
!
 1092 format(" total particle ke=",1pe10.3," total grid ke=",1pe10.3)
 1093 format(" total particle int nrg=",1pe10.3,    &
           " total grid int nrg=",1pe10.3)
 1094 format(" total particle nrg=",1pe10.3," total grid nrg=",1pe10.3)
 1095 format(5e12.4)
 1096 format(6e12.4)
 1097 format(8e12.4)
!
      call cpu_time(Tfinish)
      do l=1,20
        if(RoutineName(l).eq.'budget') CPUTime(l)=CPUTime(l)+Tfinish-Tstart
      enddo
      return
      end
