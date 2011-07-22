      subroutine parcelc(bxn,byn,bzn)
!
      use vast_kind_param, ONLY : double
      use corgan_com_M, ONLY : itdim, npart
      use cindex_com_M, ONLY : ncells, nvtx,     &
          iwid, jwid, kwid,                      &
          ibar, jbar, kbar,                      &
          ibp1, jbp1, kbp1
      use blcom_com_M, ONLY:  ijkcell, ijkvtx, ijkctmp,   &
          periodic_x, periodic_y, periodic_z,    &
          iphead, link,                   &
          mupx, mupy, mupz, number,              &
          ep, mass,                              &
          sie1p, wate, mc,                       &
          pxi, peta, pzta
      use Scratch_com_M, ONLY : the, zeta, nu,   &
          mc_tmp, sie1p_tmp, number_tmp,         &
          bxn_tmp, byn_tmp, bzn_tmp
      use Timing
      implicit none
!
      integer :: MyThread, NumThreads,     &
          ibp2,jbp2,kbp2,    &
          i,j,k,ijk,    &
          l,ll,n,np,ntmp,ntmp_nxt,   &
          nptotl
!
      integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!
      real(double) :: bxn(*), byn(*), bzn(*)
      real(double) :: wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp
      real(double) :: Tstart, Tfinish, TotalMass, sumwate
      real(double) :: numberLocal, massLocal, bxLocal, byLocal, bzLocal, nrgLocal
      integer :: ijkstep(27)
!
!      a routine to interpolate particle data to the grid
!
      allocate (the(itdim), zeta(itdim), nu(itdim))
      allocate (mc_tmp(itdim,8), sie1p_tmp(itdim,8))
      allocate (bxn_tmp(itdim,8), byn_tmp(itdim,8), bzn_tmp(itdim,8))
      allocate (number_tmp(itdim,8))
!
      call cpu_time(Tstart)
      ibp2=ibp1+1
      jbp2=jbp1+1
      kbp2=kbp1+1
!
      call celstep(iwid,jwid,kwid,ijkstep)
!
       ntmp=1+ibp1*iwid+jbp1*jwid+kbp1*kwid

!$omp parallel do private(n,l)
     do n=1,ntmp
       mc(n)=0.0
       sie1p(n)=0.0
       bxn(n)=0.0
       byn(n)=0.0
       bzn(n)=0.0
       number(n)=0.0d0
!
          do l=1,8
            mc_tmp(n,l)=0.0
            sie1p_tmp(n,l)=0.0
            bxn_tmp(n,l)=0.0
            byn_tmp(n,l)=0.0
            bzn_tmp(n,l)=0.0
            number_tmp(n,l)=0.0d0
           enddo
     enddo
!
      call OMP_SET_NUM_THREADS(8)
!
      nptotl=0
!     calculate the mass
!
!$omp parallel shared (ijkstep),                    &
!$omp shared(NumThreads)
      NumThreads=OMP_GET_NUM_THREADS()
      MyThread=OMP_GET_THREAD_NUM()+1
!$omp do reduction(+:nptotl), private(ijk,np,i,j,k,ll),    &
!$omp private(wi,wim,wip,wj,wjm,wjp,wk,wkm,wkp),             &
!$omp private(numberLocal,massLocal,nrgLocal,bxLocal,byLocal,bzLocal)
      
      do n=1,ncells
         ijk=ijkcell(n)
         np=iphead(ijk)
!
         do while (np.gt.0)

           nptotl=nptotl+1
           i=int(pxi(np))
           j=int(peta(np))
           k=int(pzta(np))
!
           the(n)=pxi(np)-real(i)-0.5
           zeta(n)=peta(np)-real(j)-0.5
           nu(n)=pzta(np)-real(k)-0.5
!
           wi=0.75-the(n)**2
           wim=0.5*(0.5-the(n))**2
           wip=0.5*(0.5+the(n))**2
!
           wj=0.75-zeta(n)**2
           wjm=0.5*(0.5-zeta(n))**2
           wjp=0.5*(0.5+zeta(n))**2
!
           wk=0.75-nu(n)**2
           wkm=0.5*(0.5-nu(n))**2
           wkp=0.5*(0.5+nu(n))**2
!
!     k-plane
!
           wate(n,1)=wi*wj*wk
           wate(n,2)=wip*wj*wk
           wate(n,3)=wip*wjp*wk
           wate(n,4)=wi*wjp*wk
           wate(n,5)=wim*wjp*wk
           wate(n,6)=wim*wj*wk
           wate(n,7)=wim*wjm*wk
           wate(n,8)=wi*wjm*wk
           wate(n,9)=wip*wjm*wk
!
!     k-1 - plane
!
           wate(n,10)=wi*wj*wkm
           wate(n,11)=wip*wj*wkm
           wate(n,12)=wip*wjp*wkm
           wate(n,13)=wi*wjp*wkm
           wate(n,14)=wim*wjp*wkm
           wate(n,15)=wim*wj*wkm
           wate(n,16)=wim*wjm*wkm
           wate(n,17)=wi*wjm*wkm
           wate(n,18)=wip*wjm*wkm
!
!     k+1 - plane
!
           wate(n,19)=wi*wj*wkp
           wate(n,20)=wip*wj*wkp
           wate(n,21)=wip*wjp*wkp
           wate(n,22)=wi*wjp*wkp
           wate(n,23)=wim*wjp*wkp
           wate(n,24)=wim*wj*wkp
           wate(n,25)=wim*wjm*wkp
           wate(n,26)=wi*wjm*wkp
           wate(n,27)=wip*wjm*wkp
!
            do ll=1,27
                number_tmp(ijk+ijkstep(ll),MyThread)=number_tmp(ijk+ijkstep(ll),MyThread)   &
                  +wate(n,ll)
                mc_tmp(ijk+ijkstep(ll),MyThread)=mc_tmp(ijk+ijkstep(ll),MyThread)           &
                  +mass(np)*wate(n,ll)
                bxn_tmp(ijk+ijkstep(ll),MyThread)=bxn_tmp(ijk+ijkstep(ll),MyThread)         &
                  +mupx(np)*wate(n,ll)
                byn_tmp(ijk+ijkstep(ll),MyThread)=byn_tmp(ijk+ijkstep(ll),MyThread)         &
                  +mupy(np)*wate(n,ll)
                bzn_tmp(ijk+ijkstep(ll),MyThread)=bzn_tmp(ijk+ijkstep(ll),MyThread)         &
                  +mupz(np)*wate(n,ll)
                sie1p_tmp(ijk+ijkstep(ll),MyThread)=sie1p_tmp(ijk+ijkstep(ll),MyThread)     &
                  +ep(np)*wate(n,ll)

            enddo
            np=link(np)
         enddo
      enddo
!$omp end parallel
!
      call OMP_SET_NUM_THREADS(1)
!
!
!!!!!!      call OMP_SET_NUM_THREADS(8)
!
      nptotl=0
      TotalMass=0.0
!!!!!$omp parallel do private(l), reduction(+:TotalMass, nptotl),   &
!!!!!$omp shared(mc,sie1p,bxn,byn,bzn,number),   &
!!!!!$omp shared(NumThreads)
!
      do n=1,ntmp
!!!!         do l=1,NumThreads
       do l=1,8
            mc(n)=mc(n)+mc_tmp(n,l)
            sie1p(n)=sie1p(n)+sie1p_tmp(n,l)
            bxn(n)=bxn(n)+bxn_tmp(n,l)
            byn(n)=byn(n)+byn_tmp(n,l)
            bzn(n)=bzn(n)+bzn_tmp(n,l)
            number(n)=number(n)+number_tmp(n,l)
         enddo
         TotalMass=TotalMass+mc(n)
         nptotl=nptotl+number(n)
      enddo
!
!`      call OMP_SET_NUM_THREADS(1)

!     impose boundary conditions on particle data
!
      nptotl=0
      do n=1,ncells
         ijk=ijkcell(n)
         nptotl=nptotl+number(ijk)
      enddo
!    bottom and top boundaries
      do i=1,ibp2
         do j=1,jbp2
            ijk=1+(i-1)*ibar+(j-1)*jbar
            nptotl=nptotl+number(ijk)
            ijk=1+(i-1)*ibar+(j-1)*jbar+kbp1*kbar
            nptotl=nptotl+number(ijk)
         enddo
      enddo
!     back and front
      do i=1,ibp2
         do k=2,kbp1
            ijk=1+(i-1)*ibar+(k-1)*kbar
            nptotl=nptotl+number(ijk)
            ijk=1+(i-1)*ibar+jbp1*jbar+(k-1)*kbar
            nptotl=nptotl+number(ijk)
         enddo
      enddo
!right and left
      do j=2,jbp1
         do k=2,kbp1
            ijk=1+(j-1)*jwid+(k-1)*kwid
            nptotl=nptotl+number(ijk)
            ijk=1+ibp1*iwid+(j-1)*jwid+(k-1)*kwid
            nptotl=nptotl+number(ijk)
         enddo
      enddo
      call bc_particles(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          number)
!
      nptotl=0
      do n=1,ncells
         ijk=ijkcell(n)
         nptotl=nptotl+number(ijk)
      enddo
      call bc_particles(ibar,jbar,kbar,iwid,jwid,kwid,   &
          periodic_x,periodic_y,periodic_z,   &
          mc)
!
      call bc_particles(ibar,jbar,kbar,iwid,jwid,kwid,   &
          periodic_x,periodic_y,periodic_z,    &
          bxn)
!
      call bc_particles(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          byn)
!
      call bc_particles(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          bzn)
!
      call bc_particles(ibar,jbar,kbar,iwid,jwid,kwid,    &
          periodic_x,periodic_y,periodic_z,    &
          sie1p)
!
      call cpu_time(Tfinish)
      do l=1,20
        if(RoutineName(l).eq.'parcelc') CPUTime(l)=CPUTime(l)+Tfinish-Tstart
      enddo
!
      deallocate (the, zeta, nu)
      deallocate (mc_tmp, sie1p_tmp)
      deallocate (bxn_tmp, byn_tmp, bzn_tmp)
      deallocate (number_tmp)

      return
      end subroutine parcelc
