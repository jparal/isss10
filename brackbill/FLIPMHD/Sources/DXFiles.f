       module DXFiles
!
      use vast_kind_param, ONLY : double
      use corgan_com_M
      use cindex_com_M
      integer, dimension(25,nspecies) :: DXUnit
      character(len=12), dimension(25,nspecies) :: DXFileName
      character(len=9), dimension(25,nspecies) :: DXName
      real(double), dimension(0:999) :: DXTime
!
      contains

      subroutine OpenDXFiles
!
!     a routine to open files for current run
!
      use cplot_com_M, ONLY : iout
      use blcom_com_M, ONLY : SpeciesName
      integer l,is
!
      DXUnit(1:25,1:nspecies) = 0
!
!     grid data
!
!
!     open RSpeciesDensity files
      if(iout(2).gt.0) then
         do is=1,nsp
         DXName(2,is)='R'//SpeciesName(is)
         DXFileName(2,is) =DXName(2,is)//'.dx'
         DXUnit(2,is) = 60+is
         open(DXUnit(2,is),file=DXFileName(2,is),status='replace')
         enddo
       
      endif
!
!     open Current  files
      if(iout(8).gt.0) then
         DXName(8,1)='JCurrent'
         DXFileName(8,1) ='J.dx'
         DXUnit(8,1) = 65
         open(DXUnit(8,1),file=DXFileName(8,1),status='replace')
      endif

!     open Velocity files
      if(iout(1).gt.0) then
         DXName(1,1)='Velocity'
         DXFileName(1,1) ='V.dx'
         DXUnit(1,1) = 70
         open(DXUnit(1,1),file=DXFileName(1,1),status='replace')
      endif
!
!     open Pressure files
      if(iout(3).gt.0) then
         DXName(3,1)='Pressure'
         DXFileName(3,1) ='P.dx'
         DXUnit(3,1) = 85
         open(DXUnit(3,1),file=DXFileName(3,1),status='replace')
      endif

!     open div(V) files
      if(iout(5).gt.0) then
         DXName(5,1)='div(V)'
         DXFileName(5,1) ='div(V).dx'
         DXUnit(5,1) = 75
         open(DXUnit(5,1),file=DXFileName(5,1),status='replace')
      endif
!
!     open internal energy files
      if(iout(6).gt.0) then
         DXName(6,1)='EInrg'
         DXFileName(6,1) ='EInrg.dx'
         DXUnit(6,1) = 80
         open(DXUnit(6,1),file=DXFileName(6,1),status='replace')
      endif
!

!     open  B-field file
!
      if(iout(9).gt.0) then
        DXName(9,1)='BField'
        DXFileName(9,1)='B.dx'
        DXUnit(9,1)=52
        open(DXUnit(9,1), file=DXFileName(9,1),status='replace')
      endif
!
      return
      end subroutine OpenDXFiles
!
      subroutine DXStart
!
!      a routine to write grid and connection data at the beginning of each .dx file
!
      use blcom_com_M, ONLY : x, y, z, ijkcell
      use cindex_com_M, ONLY : ibar, jbar, kbar, nsp
      use cophys_com_M, ONLY : dx, dy, dz
      use cplot_com_M, ONLY : iout
!
      integer :: l, is
!
      do l=1,25
         if(iout(l).gt.0) then
         do is=1,nsp
            if(DXUnit(l,is).gt.0) then
!
            write(DXUnit(l,is),102) kbar+1,jbar+1,ibar+1
            write(DXUnit(l,is),104) z(ijkcell(1)), y(ijkcell(1)), x(ijkcell(1))
            write(DXUnit(l,is),101) dz, zero, zero
            write(DXUnit(l,is),101) zero, dy, zero
            write(DXUnit(l,is),101) zero, zero, dx
!
            write(DXUnit(l,is),103) kbar+1,jbar+1,ibar+1
            write(DXUnit(l,is),*)'attribute "element type" string "cubes"'
!
            endif
         enddo
         endif
       enddo
!
  101 format('delta',3(f12.3))
  102 format('object',2x,' "grid" ',2x,'class gridpositions ',2x,' counts ',i3,2x,i3,2x,i3)
  103 format('object',2x,'"connections"',2x,'class gridconnections counts',1x,i3,2x,i3,2x,i3)
  104 format('origin',3(f12.3))
!
      return
      end subroutine DXStart
!
      subroutine DXOutput
!
      use vast_kind_param, ONLY : double
      use Scratch_com_M
      use cophys_com_M, ONLY : cdlt, sdlt, dx, dy, dz
      use numpar_com_M, ONLY : t
      use blcom_com_M, ONLY :    &
         ijkcell, ijkvtx,     &
          p,             &
          periodic_x, periodic_y, periodic_z,   &
          x, y, z,                       &
          ul, vl, wl,                    &
          bxn, byn, bzn,                 &
          bxv, byv, bzv,                 &
          jx, jy, jz,                    &
          mv, vol, vvol,                 &
          SpeciesName,                   &
          rho
      use cplot_com_M, ONLY : iout, iplot
!
      real(double) :: zero
      zero=0.0d0
      DXTime(iplot)=t
            allocate(plotx(itdim),ploty(itdim),plotz(itdim))
! 
      do l=1,25
         if(l.eq.2.and.iout(2).gt.0) then
            do is=1,nsp
            write(DXUnit(2,is),*)'object ',iplot,' class array type float rank 0 items ',nvtx
            write(DXUnit(2,is),*)'data follows'
               do n=1,nvtx
               write(DXUNit(2,is),101) mv(ijkvtx(n))/vvol(ijkvtx(n))
               enddo
            enddo
         endif
!
         if(l.eq.1.and.iout(1).gt.0) then
!
         do is=1,nsp
            write(DXUNit(1,is),*)'object',iplot,' class array type float rank 1 shape 3',  &
              ' items ',nvtx,' data follows'
               do n=1,nvtx
                  write(DXUnit(1,is),105) ul(ijkvtx(n)), vl(ijkvtx(n)),wl(ijkvtx(n))
               enddo
         enddo
         endif


!
         if(l.eq.8.and.iout(8).gt.0) then
!
      call curlv(nvtx,ijkvtx,      &
          bxn,byn,bzn,jx,jy,jz)
!
!     impose periodic boundary conditions
!
      zero=0.0
!
      call bc_current(ibp1+1,jbp1+1,kbp1+1,    &
          zero,                                &
          jx,jy,jz)

         write(DXUNit(8,1),*)'object',iplot,' class array type float rank 1 shape 3',' items ',nvtx,' data follows'
            do n=1,nvtx
            write(DXUnit(8,1),105) jx(ijkvtx(n)), jy(ijkvtx(n)),jz(ijkvtx(n))
            enddo
!
         endif
!
         if(l.eq.9.and.iout(9).gt.0) then
      call b_vtx(ncells,ijkcell,iwid,jwid,kwid,    &
         nvtx,ijkvtx,    &
         bxn,byn,bzn,    &
         vol,vvol,    &
         bxv,byv,bzv)

         write(DXUnit(9,1),*)'object',iplot,' class array type float rank 1 shape 3',' items ',nvtx
         write(DXUnit(9,1),*) ' data follows'
            do n=1,nvtx
            write(DXUnit(9,1),105) bxv(ijkvtx(n)), byv(ijkvtx(n)),bzv(ijkvtx(n))
            enddo

         endif

!         if(l.eq.6.and.iout(6).gt.0) then

!      call PhiVertex(sie,plotx)

!         write(DXUnit(6,1),*)'object',iplot,'  class array type float rank 0 items ',nvtx
!         write(DXUnit(6,1),*)'data follows'
!            do n=1,nvtx
!            write(DXUnit(6,1),101) plotx(ijkvtx(n))
!            enddo

!         endif


      do is=1,nsp
         if(DXUnit(l,is).gt.0) then
         write(DXUnit(l,is),*)'attribute "dep" string "positions"'
!
         if(iplot.le.9) then
         write(DXUnit(l,is),106) iplot
  106    format('object',' "field',i1,'"',' class field')
         elseif(iplot.le.99) then
         write(DXUnit(l,is),107) iplot
  107    format('object',' "field',i2,'"',' class field')
         elseif(iplot.le.999) then
         write(DXUnit(l,is),109) iplot
  109    format('object',' "field',i3,'"',' class field')
         endif
         write(DXUnit(l,is),*)'component', ' "positions" ',' value ', ' "grid" '
         write(DXUnit(l,is),*)'component', ' "connections" ',' value ', ' "connections"  '
         write(DXUnit(l,is),108) iplot
  108    format('component "data" value ',i3)
         endif
      enddo
   enddo
!
  101 format(f16.9)
  102 format('object',2x,' "grid" ',2x,'class gridpositions ',2x,' counts',i3,2x,i3,2x,i3)
  103 format('object',2x,'"connections"',2x,'class gridconnections counts',1x,i3,2x,i3,2x,i3)
  105 format(3(f12.5))
  115 format(9(f12.5))
!
      deallocate(plotx,ploty,plotz)
      return
      end subroutine DXOutput
!
        subroutine DXFinish
!
!       a routine to write series objects
!
        use cplot_com_M, ONLY : iout, iplot
        use numpar_com_M, ONLY : t
!
        integer l, is, i
!
        do l=1,25
           if(iout(l).gt.0) then
             do is=1,nsp
                if(DXUnit(l,is).gt.0) then
                   write(DXUnit(l,is),101) DXName(l,is)
    101            format('object "',A9,'" class series')
                   do i=0,iplot
                      if(i.le.9) then
                         write(DXUnit(l,is),102) i, DXTime(i), i
    102                  format('member ',i3,' position ',f10.3,' value "field',i1,'"')
                      elseif(i.le.99) then
                         write(DXUnit(l,is),103) i, DXTime(i), i
    103                  format('member ',i3,' position ',f10.3,' value "field',i2,'"')
                      elseif(i.le.999) then
                         write(DXUnit(l,is),104) i, DXTime(i), i
    104                  format('member ',i3,' position ',f10.3,' value "field',i3,'"')
                      endif
                   enddo
                   write(DXUnit(l,is),105) DXName(l,is)
    105            format('attribute "name" string "',A9,'"')
                   write(DXUnit(l,is),*) 'end'
                endif
             enddo
          endif
       enddo
!
       return
       end subroutine DXFinish

      subroutine CloseDXFiles
!
!     a routine to close files for current run
!
      use cplot_com_M, ONLY : iout
!
      integer l, is
!
!     close species density files
      do l=1,25
         if(iout(l).gt.0) then
            do is=1,nsp
               if(DXUnit(l,is).gt.0) then
                  close(DXUnit(l,is),status='keep')
               endif
            enddo
         endif
      enddo
!
      return
      end subroutine CloseDXFiles

      end module DXFiles
