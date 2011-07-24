      subroutine parcelv(is)
!
      use vast_kind_param, only:  double
      use corgan_com_M, ONLY : itdim
      use cindex_com_M, ONLY :  ncells, nvtx, nsp,   &
        iwid, jwid, kwid
      use numpar_com_M, ONLY : dt
      use blcom_com_M, ONLY:                &
        ijkcell, ijkvtx, ijkctmp,           &
        mass, pxi, peta, pzta,              &
        up, vp, wp,link, ico,               &
        wate,                               &
        mv_s, umom, vmom, wmom, numberv, color
      use Scratch_com_M, ONLY : the, zeta, nu,  &
        mv_tmp, umom_tmp, vmom_tmp, wmom_tmp, numberv_tmp, color_tmp
      use Timing
      use ParticleLists
      implicit none
!
      integer :: is
      real(double) ::     &
               Tstart, Tfinish
!
      integer :: ijkvstep(8),   &
        inew, jnew, knew, npar
!
      real(double) ::    &
         wi, wim, wip,   &
         wj, wjm, wjp,   &
         wk, wkm, wkp
!
      real(double) :: parmass, vtxmass
      integer ::   nptotl,     &
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
      do n=1,nvtx
         ijk=ijkvtx(n)
         mv_s(ijk,is)=0.0
      enddo
!
      npar=0
!
!
      allocate (the(itdim), zeta(itdim), nu(itdim))
      allocate (mv_tmp(itdim,8), umom_tmp(itdim,8), vmom_tmp(itdim,8), wmom_tmp(itdim,8))
      allocate (numberv_tmp(itdim,8), color_tmp(itdim,8))
!
!
      mv_tmp=0.
      umom_tmp=0.
      vmom_tmp=0.
      wmom_tmp=0.
      numberv_tmp=0.
      color_tmp=0.0
!
!
      parmass=0.0
      nptotl=0
      vtxmass=0.
!
        call OMP_SET_NUM_THREADS(8)
!$omp parallel shared(NumThreads,is)
      NumThreads=OMP_GET_NUM_THREADS()
      MyThread=OMP_GET_THREAD_NUM() + 1
!$omp do reduction(+:parmass,nptotl),  &
!$omp private(ijk,np,l,inew,jnew,knew),   &
!$omp private(wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp)

      do n=1,ncells
!
      ijk=ijkcell(n)
      np=ipheadS(ijk,is)
      do while (np.gt.0)
!
      nptotl=nptotl+1
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

      mv_tmp(ijk+ijkvstep(l),MyThread)=mv_tmp(ijk+ijkvstep(l),MyThread)+mass(np)   &
         *wate(ijk,l)
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

   vtxmass=0.0
!    do n=1,nvtx
!       ijk=ijkvtx(n)
!       vtxmass=vtxmass+mv_tmp(ijk,1)
!    enddo
        call OMP_SET_NUM_THREADS(8)
!$omp parallel do private(ijk,l), reduction(+:vtxmass),  &
!$omp shared(mv_s,umom,vmom,wmom,numberv,color,is),    &
!$omp shared(NumThreads)

     do n=1,nvtx
       ijk=ijkvtx(n)
       mv_s(ijk,is)=0.0
       do l=1,NumThreads
          mv_s(ijk,is)=mv_s(ijk,is)+mv_tmp(ijk,l)
          vtxmass=vtxmass+mv_tmp(ijk,l)
          umom(ijk)=umom(ijk)+umom_tmp(ijk,l)
          vmom(ijk)=vmom(ijk)+vmom_tmp(ijk,l)
          wmom(ijk)=wmom(ijk)+wmom_tmp(ijk,l)
          numberv(ijk)=numberv(ijk)+numberv_tmp(ijk,l)
          color(ijk)=color(ijk)+color_tmp(ijk,l)
       enddo
!       vtxmass=vtxmass+mv_s(ijk,is)
     enddo
!
!
      call cpu_time(Tfinish)
      do l=1,20
         if(RoutineName(l).eq.'parcelv') CPUTime(l)=CPUTime(l)+Tfinish-Tstart
      enddo
!
      deallocate (the,zeta,nu)
      deallocate (mv_tmp, umom_tmp, vmom_tmp, wmom_tmp, numberv_tmp, color_tmp)
      return
      end subroutine parcelv
