      module ParticleLists
!
!     routines to create and destroy particle lists
!
      use corgan_com_M, ONLY : itdim, nspecies
      integer, dimension(itdim,nspecies) :: ipheadS
!
      contains
!
!     *******************************************************************
!
      subroutine SingleToSpecies
!
!     a routine to construct lists of particles of the same species
!
      use cindex_com_M, ONLY : ncells
      use blcom_com_M, ONLY : ijkcell, iphead, link, ico,  &
         ijkctmp
!
      integer :: n, np, ijk, is, newcell, newcell_nxt
      logical :: nomore
!
      ipheadS(1:itdim,1:nspecies) = 0
      do n=1,ncells
         ijk=ijkcell(n)
         np=iphead(ijk)
         do while (np.gt.0)
            is=ico(np)
            iphead(ijk)=link(np)
            link(np)=ipheadS(ijk,is)
            ipheadS(ijk,is)=np
            np=iphead(ijk)
          enddo
       enddo
!
       return
       end subroutine SingleToSpecies
!
!      *******************************************************************
!
       subroutine SpeciesToSingle
!
!      a routine to consolidate species lists into a single particle list

!
       use cindex_com_M, ONLY : ncells, nsp
       use blcom_com_M, ONLY : ijkcell, iphead, link,   &
          ijkctmp
       logical nomore
!
       integer :: is, n, np, ijk, newcell, newcell_nxt
!
!
!
      do is=1,nsp
         newcell=ncells
         do n=1,ncells
            ijkctmp(n)=ijkcell(n)
         enddo
   1     continue
      nomore=.TRUE.
      newcell_nxt=0
      do n=1,newcell
         if(ipheadS(ijkctmp(n),is).le.0) cycle
            newcell_nxt=newcell_nxt+1
            ijkctmp(newcell_nxt)=ijkctmp(n)
      enddo
!
      newcell=newcell_nxt
      if(newcell.ne.0) then
         nomore=.FALSE.

!
          do n=1,newcell
             ijk=ijkctmp(n)
             if(ipheadS(ijk,is).eq.0) cycle
                np=ipheadS(ijk,is)
                ipheadS(ijk,is)=link(np)
                link(np)=iphead(ijk)
                iphead(ijk)=np
          enddo
!
       endif
!
       if(.not.nomore) go to 1
!
       enddo
!
       return
       end subroutine SpeciesToSingle
!
!     ********************************************************************

       end module ParticleLists

