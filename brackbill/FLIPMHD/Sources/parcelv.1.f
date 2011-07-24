      subroutine parcelv
!
      use vast_kind_param, only:  double
      use corgan_com_M, ONLY : itdim, nspecies
      use cindex_com_M, ONLY :  ncells, nvtx, nsp,   &
        iwid, jwid, kwid
      use numpar_com_M, ONLY : dt
      use blcom_com_M, ONLY: iphead, iphd2, &
        ijkcell, ijkvtx, ijkctmp,           &
        mass, pxi, peta, pzta,              &
        up, vp, wp,link, ico,               &
        wate,                               &
        mv, mv_s, umom, vmom, wmom, numberv, color
      use Scratch_com_M, ONLY : the, zeta, nu,  &
        mv_tmp_s, umom_tmp, vmom_tmp, wmom_tmp, numberv_tmp, color_tmp
      use Timing
      implicit none
!
      real(double) ::     &
               Tstart, Tfinish
!
      integer :: ijkvstep(8),   &
        inew, jnew, knew, is
!
      real(double) ::    &
         wi, wim, wip,   &
         wj, wjm, wjp,   &
         wk, wkm, wkp
!
      real(double) :: parmass, vtxmass
      integer ::     &
          ijk,l,n,MyThread,     &
          NumThreads,           &
          newcell,newcell_nxt,np
      integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM


!      a routine to interpolate particle data to the vertices of grid
!
      call cpu_time(Tstart)
      call vtxindx(iwid,jwid,kwid,ijkvstep)
!
!
        do 11 l=1,8
      do 11 n=1,nvtx
        wate(ijkvtx(n),l)=0.0
   11   continue
!
!
!     **************************************
!
!     set accumulators to zero
!
!     ********************************************
!
      do 25 n=1,nvtx
      ijk=ijkvtx(n)
      mv(ijk)=0.0
      do is=1,nsp
         mv_s(ijk,is)=0.0
      enddo
      umom(ijk)=0.0
      vmom(ijk)=0.0
      wmom(ijk)=0.0
      numberv(ijk)=0.0d0
      color(ijk)=0.0d0
   25 continue
!
!
!
      allocate (the(itdim), zeta(itdim), nu(itdim))
      allocate (mv_tmp_s(itdim,nspecies,8),    &
        umom_tmp(itdim,8), vmom_tmp(itdim,8), wmom_tmp(itdim,8))
      allocate (numberv_tmp(itdim,8), color_tmp(itdim,8))
!
!
!
!
      parmass=0.0
      vtxmass=0.
!
        call OMP_SET_NUM_THREADS(8)
!$omp parallel shared(NumThreads)
      NumThreads=OMP_GET_NUM_THREADS()
      MyThread=OMP_GET_THREAD_NUM() + 1
!$omp do reduction(+:parmass),  &
!$omp private(ijk,np,is,l,inew,jnew,knew),   &
!$omp private(wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp)

      do n=1,ncells
!
      ijk=ijkcell(n)
      np=iphead(ijk)
      do while (np.gt.0)
!
      is=ico(np)
!
      inew=int(pxi(np))
      jnew=int(peta(np))
      knew=int(pzta(np))
!
!
      the(n)=pxi(np)-real(inew)
      zeta(n)=peta(np)-real(jnew)
      nu(n)=pzta(np)-real(knew)
!
!     calculate interpolation weights
!
      wi=1.-the(n)
      wip=the(n)
!
      wj=1.-zeta(n)
      wjp=zeta(n)
!
      wk=1.-nu(n)
      wkp=nu(n)
!
      wate(ijk,1)=wip*wj*wk
      wate(ijk,2)=wip*wjp*wk
      wate(ijk,3)=wi*wjp*wk
      wate(ijk,4)=wi*wj*wk
!
      wate(ijk,5)=wip*wj*wkp
      wate(ijk,6)=wip*wjp*wkp
      wate(ijk,7)=wi*wjp*wkp
      wate(ijk,8)=wi*wj*wkp

      parmass=parmass+mass(np)
!
       do l=1,8

      mv_tmp_s(ijk+ijkvstep(l),1,MyThread)      &
         =mv_tmp_s(ijk+ijkvstep(l),1,MyThread)+mass(np)   &
         *wate(ijk,l)
!
!      mv_tmp_s(ijk+ijkvstep(l),is,MyThread)      &
!         =mv_tmp_s(ijk+ijkvstep(l),is,MyThread)+mass(np)   &
!         *wate(ijk,l)
!
      umom_tmp(ijk+ijkvstep(l),MyThread)=umom_tmp(ijk+ijkvstep(l),MyThread)+   &
            mass(np)*up(np)*wate(ijk,l)
!
      vmom_tmp(ijk+ijkvstep(l),MyThread)=vmom_tmp(ijk+ijkvstep(l),MyThread)+   &
          mass(np)*vp(np)*wate(ijk,l)
!
      wmom_tmp(ijk+ijkvstep(l),MyThread)=wmom_tmp(ijk+ijkvstep(l),MyThread)+   &
           mass(np)*wp(np)*wate(ijk,l)
!
      numberv_tmp(ijk+ijkvstep(l),MyThread)=numberv_tmp(ijk+ijkvstep(l),MyThread)   &
          +wate(ijk,l)
!
      color_tmp(ijk+ijkvstep(l),MyThread)=color_tmp(ijk+ijkvstep(l),MyThread)    &
         +ico(np)*wate(ijk,l)
!
       enddo
       np=link(np)
       enddo
    enddo
!
!$omp end parallel
        call OMP_SET_NUM_THREADS(1)

!        call OMP_SET_NUM_THREADS(8)
!$omp parallel do private(ijk,l,is), reduction(+:vtxmass),  &
!$omp shared(nsp,mv_s,umom,vmom,wmom,numberv,color),    &
!$omp shared(NumThreads)

     do n=1,nvtx
       ijk=ijkvtx(n)
       do l=1,NumThreads
!          do is=1,nsp
             is=1
             mv_s(ijk,is)=mv_s(ijk,is)+mv_tmp_s(ijk,is,l)
!             mv(ijk)=mv(ijk)+mv_tmp_s(ijk,is,l)
!          enddo
          umom(ijk)=umom(ijk)+umom_tmp(ijk,l)
          vmom(ijk)=vmom(ijk)+vmom_tmp(ijk,l)
          wmom(ijk)=wmom(ijk)+wmom_tmp(ijk,l)
          numberv(ijk)=numberv(ijk)+numberv_tmp(ijk,l)
          color(ijk)=color(ijk)+color_tmp(ijk,l)
       enddo
       vtxmass=vtxmass+mv(ijk)
       do is=1,nsp
          vtxmass=vtxmass-mv_s(ijk,is)
       enddo
       if(vtxmass.gt.1.d-10) then
         write(*,*) 'parcelv: vtxmass=',vtxmass
         stop
       endif
     enddo
!
     do n=1,nvtx
       ijk=ijkvtx(n)
       do is=1,nsp
          mv(ijk)=mv(ijk)+mv_s(ijk,is)
       enddo
     enddo
!
      call cpu_time(Tfinish)
      do l=1,20
         if(RoutineName(l).eq.'parcelv') CPUTime(l)=CPUTime(l)+Tfinish-Tstart
      enddo
!
      deallocate (the,zeta,nu)
      deallocate (mv_tmp_s, umom_tmp, vmom_tmp,      &
         wmom_tmp, numberv_tmp, color_tmp)
      return
      end subroutine parcelv
