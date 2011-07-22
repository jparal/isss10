      module Timing
      use vast_kind_param, ONLY : double
      character(len=15), dimension(1:20) :: RoutineName
      real(double), dimension(1:20) :: CPUTime
      integer nSamples
      data nSamples/14/
      data RoutineName(1)/'parcelv'/
      data RoutineName(2)/'parmov'/
      data RoutineName(3)/'poisson_vtx'/
      data RoutineName(4)/'mfnk'/
      data RoutineName(5)/'budget'/
      data RoutineName(6)/'parcelc'/
      data RoutineName(7)/'gmres_vtx'/
      data RoutineName(8)/'gmres_3dmhd'/
      data RoutineName(13)/'parset'/

      contains
      subroutine TimingOutput(Tused_flip)
      use vast_kind_param, ONLY : double
      real(double) :: Tused_flip
      integer :: n
      CPUTime=100.*CPUTime/Tused_flip
      write(41,*) 'Tused_flip=',Tused_flip
      write(41,*) 'Routine Name    % Cpu Time'
      do n=1,nSamples
         if(CPUTime(n).gt.0.0) then
            write(41,*) RoutineName(n), CPUTime(n)
         endif
      enddo
      return
      end subroutine TimingOutput

      end module Timing
