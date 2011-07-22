      subroutine celindex(i1,i2,j1,j2,k1,k2,iwid,jwid,kwid,  &
         ijkcell,ijkctmp,ncells,ijkvtx,nvtxkm,nvtx)
!
      implicit none
      integer, intent(in) :: i1, i2, j1, j2, k1, k2
      integer, intent(in) :: iwid, jwid, kwid
      integer  ijkcell(*), ijkctmp(*), ijkvtx(*)
      integer ncells, nvtxkm, nvtx
!
      integer i, j, k
!
      ncells=0
!
      do 50 k=k1,k2
      do 50 j=j1,j2
      do 50 i=i1,i2
!
      ncells=ncells+1
!
      ijkcell(ncells)=1+(i-1)*iwid+(j-1)*jwid+(k-1)*kwid
      ijkctmp(ncells)=ijkcell(ncells)
!
   50 continue
!
!
      nvtx=0
!
      do 60 k=k1,k2
      do 60 j=j1,j2+1
      do 60 i=i1,i2+1
!
      nvtx=nvtx+1
!
      ijkvtx(nvtx)=1+(i-1)*iwid+(j-1)*jwid+(k-1)*kwid
!
   60 continue
!
      nvtxkm=nvtx
!
      do 70 j=j1,j2+1
      do 70 i=i1,i2+1
!
      nvtx=nvtx+1
!
      ijkvtx(nvtx)=1+(i-1)*iwid+(j-1)*jwid+k2*kwid
!
   70 continue
      return
      end subroutine celindex
