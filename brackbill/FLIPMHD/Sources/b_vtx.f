      subroutine b_vtx(ncells,ijkcell,iwid,jwid,kwid,   &
         nvtx,ijkvtx,                                   &
         bx,by,bz,                                      &
         vol,vvol,                                      &
         bxv,byv,bzv)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!DIR$ FREE
      USE vast_kind_param, ONLY:  double
      use blcom_com_M, ONLY :    &
        periodic_x, periodic_y, periodic_z
      use cindex_com_M, ONLY : ibp1, jbp1, kbp1
      use cophys_com_M, ONLY : cdlt, sdlt, dz
!
      implicit none
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ncells, ijkcell(*), iwid, jwid, kwid, nvtx, ijkvtx(*)
      real(double) , intent(in) :: vol(*)
      real(double) , intent(in) :: bx(*)
      real(double) , intent(in) :: by(*)
      real(double) , intent(in) :: bz(*)
      real(double), intent(out) :: vvol(*)
      real(double) , intent(out) :: bxv(*)
      real(double) , intent(out) :: byv(*)
      real(double) , intent(out) :: bzv(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, n, ijk, ijkv(8) ,is
      real(double) :: dummy, zero
      logical :: complement
!-----------------------------------------------

!     a routine to calculate a semi-correct vertex magnetic field
!
     zero=0.0
!
      call vtxindx(iwid,jwid,kwid,ijkv)
!
  do n=1, nvtx
      ijk=ijkvtx(n)
      bxv(ijk)=0.0d0
      byv(ijk)=0.0d0
      bzv(ijk)=0.0d0
  enddo
!
!
      do 1 n=1,ncells
      ijk=ijkcell(n)
!
        do l=1,8
!
          bxv(ijk+ijkv(l))=bxv(ijk+ijkv(l))  &
           +0.125*bx(ijk)*vol(ijk)
!
          byv(ijk+ijkv(l))=byv(ijk+ijkv(l))  &
           +0.125*by(ijk)*vol(ijk)
!
          bzv(ijk+ijkv(l))=bzv(ijk+ijkv(l))  &
           +0.125*bz(ijk)*vol(ijk)
!
       enddo
    1 continue
!
      call bc_vtx(ibp1+1,jbp1+1,kbp1+1,    &
          zero,                                &
          bxv,byv,bzv)
!

      do n=1,nvtx
        ijk=ijkvtx(n)
      bxv(ijk)=bxv(ijk)/vvol(ijk)
      byv(ijk)=byv(ijk)/vvol(ijk)
      bzv(ijk)=bzv(ijk)/vvol(ijk)
!
      enddo
!
      return
      end subroutine b_vtx
