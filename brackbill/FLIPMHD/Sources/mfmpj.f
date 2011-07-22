      subroutine mfmpj(ncells,ijkcell,nvtx,ijkvtx,    &
          vol,    &
          srce,    &
          residu,Aq,Ax,rhs,    &
          gradpx,gradpy,gradpz,p,diag)
!
!
      use vast_kind_param, only:  double
      use geometry_com_M
!
      integer :: ncells,ijkcell(*),nvtx,ijkvtx(*)
      real(double) ::     &
               p(*), Aq(*),Ax(*), rhs(*),     &
               vol(*),residu(*),srce(*),     &
               gradpx(*),gradpy(*),gradpz(*),diag(*)
!
!     
!     
!     apply preconditioner, i.e find y=M^(-1).q
!
      do i = 1,3
      if (i .eq. 1) then
       do n=1,nvtx
        ijk=ijkvtx(n)
        Aq(ijk)=0.0
        Ax(ijk)=0.0
       end do
      endif

      do 32 n=1,nvtx
      ijk=ijkvtx(n)
      Aq(ijk)=(rhs(ijk)-Ax(ijk))/diag(ijk) + Aq(ijk)
!      Aq(ijk)=-rhs(ijk)/diag(ijk) 
  32  continue
!
      zero=0.0d0
!
!     compute A times x (Aq)
!
      do 34 n = 1,ncells
       ijk = ijkcell(n)
       Ax(ijk) = 0.0
   34 continue
!
      do 33 n=1,nvtx
      ijk=ijkvtx(n)
      Ax(ijk)=p(ijk)+Aq(ijk)
   33 continue
!
!
      call residue_vtx(     &
          srce,     &
          Ax,Ax)
!
!
      do 35 n=1,nvtx
      ijk=ijkvtx(n)
      Ax(ijk)=(residu(ijk)-Ax(ijk))
   35 continue
      end do

      return
      end subroutine mfmpj
