      subroutine bc_particles(ibar,jbar,kbar,iwid,jwid,kwid,   &
          periodic_x,periodic_y,periodic_z,   &
          phi)
      use vast_kind_param, ONLY : double
      implicit none
!
      integer, intent(in) :: ibar, jbar, kbar, iwid, jwid, kwid
      logical, intent(in) :: periodic_x, periodic_y, periodic_z
      real(double) :: phi(*)
! 
      integer i,j,k,ibp1,jbp1,kbp1,ijkl,ijkr
!
!     a routine to impose boundary conditions for parcelc
!     data can be either periodic or Neumann in i,j,k
!     for cell-centered data
!
!      called by PARCELC
!
!
      ibp1=ibar+1
      jbp1=jbar+1
      kbp1=kbar+1
!
!     impose boundary conditions in i (logical x)
!
      if(periodic_x) then
!
        do j=1,jbp1+1
        do k=1,kbp1+1
!
        ijkl=(k-1)*kwid+(j-1)*jwid+iwid+1
        ijkr=(k-1)*kwid+(j-1)*jwid+ibar*iwid+1
!
        phi(ijkl)=phi(ijkl)+phi(ijkr+iwid)
        phi(ijkr)=phi(ijkr)+phi(ijkl-iwid)
        phi(ijkr+iwid)=phi(ijkl)
        phi(ijkl-iwid)=phi(ijkr)
!
        enddo
        enddo
!
      else
!
        do j=1,jbp1+1
        do k=1,kbp1+1
!
        ijkl=(k-1)*kwid+(j-1)*jwid+iwid+1
        ijkr=(k-1)*kwid+(j-1)*jwid+ibar*iwid+1
!
        phi(ijkl)=phi(ijkl)+phi(ijkl-iwid)
        phi(ijkr)=phi(ijkr)+phi(ijkr+iwid)
        phi(ijkl-iwid)=phi(ijkl)
        phi(ijkr+iwid)=phi(ijkr)
! 
        enddo
        enddo
!
      endif
!
!
!      impose boundary conditions in j (logical y)
!
      if(periodic_y) then
!
        do i=1,ibp1+1
        do k=1,kbp1+1
!
        ijkl=(k-1)*kwid+jwid+(i-1)*iwid+1
        ijkr=(k-1)*kwid+jbar*jwid+(i-1)*iwid+1
!
        phi(ijkl)=phi(ijkl)+phi(ijkr+jwid)
        phi(ijkr)=phi(ijkr)+phi(ijkl-jwid)
        phi(ijkr+jwid)=phi(ijkl)
        phi(ijkl-jwid)=phi(ijkr)
!
        enddo
        enddo
!
      else
!
        do i=1,ibp1+1
        do k=1,kbp1+1
!
        ijkl=(k-1)*kwid+jwid+(i-1)*iwid+1
        ijkr=(k-1)*kwid+jbar*jwid+(i-1)*iwid+1
!
        phi(ijkl)=phi(ijkl)+phi(ijkl-jwid)
        phi(ijkr)=phi(ijkr)+phi(ijkr+jwid)
        phi(ijkl-jwid)=phi(ijkl)
        phi(ijkr+jwid)=phi(ijkr)
!
        enddo
        enddo
!
      endif
!
!     impose boundary conditions in k (logical z)
!
      if(periodic_z) then
!
        do i=1,ibp1+1
        do j=1,jbp1+1
!
        ijkl=kwid+(j-1)*jwid+(i-1)*iwid+1
        ijkr=kbar*kwid+(j-1)*jwid+(i-1)*iwid+1
!
        phi(ijkl)=phi(ijkl)+phi(ijkr+kwid)
        phi(ijkr)=phi(ijkr)+phi(ijkl-kwid)
        phi(ijkr+kwid)=phi(ijkl)
        phi(ijkl-kwid)=phi(ijkr)
!
        enddo
        enddo

      else
!
        do i=1,ibp1+1
        do j=1,jbp1+1
!
        ijkl=kwid+(j-1)*jwid+(i-1)*iwid+1
        ijkr=kbar*kwid+(j-1)*jwid+(i-1)*iwid+1
!
        phi(ijkl)=phi(ijkl)+phi(ijkl-kwid)
        phi(ijkr)=phi(ijkr)+phi(ijkr+kwid)
        phi(ijkl-kwid)=phi(ijkl)
        phi(ijkr+kwid)=phi(ijkr)
!
        enddo
        enddo
!
      endif
!
      return
      end subroutine bc_particles
