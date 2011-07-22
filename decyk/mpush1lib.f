c-----------------------------------------------------------------------
c 1d PIC multi-tasking library for pushing particles and depositing
c charge
c mpush1lib.f contains multi-tasking procedures to process particles:
c MGPOST1 multi-tasking wrapper for GPOST1
c MGSPOST1 multi-tasking wrapper for GSPOST1
c MGSPOST1X multi-tasking wrapper for GSPOST1X
c MGPOST1L multi-tasking wrapper for GPOST1L
c MGSPOST1L multi-tasking wrapper for GSPOST1L
c MGSPOST1XL multi-tasking wrapper for GSPOST1XL
c MGPUSH1 multi-tasking wrapper for GPUSH1
c MGSPUSH1 multi-tasking wrapper for GSPUSH1
c MGPUSH1L multi-tasking wrapper for GPUSH1L
c MGSPUSH1L multi-tasking wrapper for GSPUSH1L
c MSORTP1X multi-tasking wrapper for SORTP1X
c MSORTP1XL multi-tasking wrapper for SORTP1XL
c MDSORTP1X multi-tasking wrapper for DSORTP1X
c MDSORTP1XL multi-tasking wrapper for DSORTP1XL
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: january 8, 2010
c-----------------------------------------------------------------------
      subroutine MGPOST1(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv)
      dimension qp(nxv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GPOST1
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GPOST1,nargs,part(1,npo),qp(1,i),qm,np
     1p,idimp,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GPOST1(part(1,npo),q,qm,npl,idimp,nxv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST1(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv)
      dimension qp(nxv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST1
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST1,nargs,part(1,npo),qp(1,i),qm,n
     1pp,idimp,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST1(part(1,npo),q,qm,npl,idimp,nxv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST1X(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv)
      dimension qp(nxv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST1X
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST1X,nargs,part(1,npo),qp(1,i),qm,
     1npp,idimp,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST1X(part(1,npo),q,qm,npl,idimp,nxv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPOST1L(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv)
      dimension qp(nxv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GPOST1L
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GPOST1L,nargs,part(1,npo),qp(1,i),qm,n
     1pp,idimp,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GPOST1L(part(1,npo),q,qm,npl,idimp,nxv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST1L(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv)
      dimension qp(nxv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST1L
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST1L,nargs,part(1,npo),qp(1,i),qm,
     1npp,idimp,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST1L(part(1,npo),q,qm,npl,idimp,nxv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST1XL(part,q,qm,nop,idimp,nxv,qp,idtask,nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv)
      dimension qp(nxv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST1XL
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST1XL,nargs,part(1,npo),qp(1,i),qm
     1,npp,idimp,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST1XL(part(1,npo),q,qm,npl,idimp,nxv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idtask,n
     1mt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, qbm, dt, ek, ekp
      integer idimp, nop, nx, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GPUSH1
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GPUSH1,nargs,part(1,npo),fx,qbm,dt,ekp
     1(i),idimp,npp,nx,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GPUSH1(part(1,npo),fx,qbm,dt,ek,idimp,npl,nx,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPUSH1(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idtask,
     1nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, qbm, dt, ek, ekp
      integer idimp, nop, nx, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSPUSH1
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSPUSH1,nargs,part(1,npo),fx,qbm,dt,ek
     1p(i),idimp,npp,nx,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSPUSH1(part(1,npo),fx,qbm,dt,ek,idimp,npl,nx,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idtask,
     1nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, qbm, dt, ek, ekp
      integer idimp, nop, nx, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GPUSH1L
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GPUSH1L,nargs,part(1,npo),fx,qbm,dt,ek
     1p(i),idimp,npp,nx,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GPUSH1L(part(1,npo),fx,qbm,dt,ek,idimp,npl,nx,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPUSH1L(part,fx,qbm,dt,ek,idimp,nop,nx,nxv,ekp,idtask
     1,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, qbm, dt, ek, ekp
      integer idimp, nop, nx, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSPUSH1L
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSPUSH1L,nargs,part(1,npo),fx,qbm,dt,e
     1kp(i),idimp,npp,nx,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSPUSH1L(part(1,npo),fx,qbm,dt,ek,idimp,npl,nx,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MSORTP1X(part,pt,ip,npic,idimp,nop,nx1,npicp,idtask,nmt
     1,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npicp
      integer idimp, nop, nx1, idtask, nmt, ierr
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(nx1)
      dimension npicp(nx1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external SORTP1X
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),SORTP1X,nargs,part(1,npo),pt(npo),ip(n
     1po),npicp(1,i),idimp,npp,nx1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call SORTP1X(part(1,npo),pt(npo),ip(npo),npic,idimp,npl,nx1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MSORTP1XL(part,pt,ip,npic,idimp,nop,nx1,npicp,idtask,nm
     1t,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npicp
      integer idimp, nop, nx1, idtask, nmt, ierr
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(nx1)
      dimension npicp(nx1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external SORTP1XL
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),SORTP1XL,nargs,part(1,npo),pt(npo),ip(
     1npo),npicp(1,i),idimp,npp,nx1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call SORTP1XL(part(1,npo),pt(npo),ip(npo),npic,idimp,npl,nx1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDSORTP1X(parta,partb,npic,idimp,nop,nx1,npicp,idtask,n
     1mt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npicp
      integer idimp, nop, nx1, idtask, nmt, ierr
      dimension parta(idimp,nop), partb(idimp,nop), npic(nx1)
      dimension npicp(nx1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external DSORTP1X
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),DSORTP1X,nargs,parta(1,npo),partb(1,np
     1o),npicp(1,i),idimp,npp,nx1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call DSORTP1X(parta(1,npo),partb(1,npo),npic,idimp,npl,nx1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDSORTP1XL(parta,partb,npic,idimp,nop,nx1,npicp,idtask,
     1nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npicp
      integer idimp, nop, nx1, idtask, nmt, ierr
      dimension parta(idimp,nop), partb(idimp,nop), npic(nx1)
      dimension npicp(nx1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external DSORTP1XL
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),DSORTP1XL,nargs,parta(1,npo),partb(1,n
     1po),npicp(1,i),idimp,npp,nx1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call DSORTP1XL(parta(1,npo),partb(1,npo),npic,idimp,npl,nx1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDPOST1GL(part,q,sctx,qm,nop,idimp,nx,nxh,nxvh,qp,sctxp
     1,idtask,nmt,ierr)
c multitasking gridless charge deposition
c qp = charge density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm
      complex q, sctx, qp, sctxp
      integer nop, idimp, nx, nxh, nxvh, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxvh), sctx(nxvh)
      dimension qp(nxvh,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external DPOST1GL
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxvh
      qp(j,i) = cmplx(0.,0.)
   10 continue
      call MP_TASKSTART(idtask(i),DPOST1GL,nargs,part(1,npo),qp(1,i),sct
     1xp(1,i),qm,npp,idimp,nx,nxh,nxvh)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call DPOST1GL(part(1,npo),q,sctx,qm,npl,idimp,nx,nxh,nxvh)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxvh
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPUSH1GL(part,fx,sctx,qbm,dt,ek,idimp,nop,nx,nxh,nxvh,s
     1ctxp,ekp,idtask,nmt,ierr)
c multitasking gridless particle push
c sctxp = scratch arrays for sines and cosines
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, ek, ekp
      complex fx, sctx, sctxp
      integer idimp, nop, nx, nxh, nxvh, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxvh), sctx(nxvh)
      dimension sctxp(nxvh,nmt), ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external PUSH1GL
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PUSH1GL,nargs,part(1,npo),fx,sctxp(1,i
     1),qbm,dt,ekp(i),idimp,npp,nx,nxh,nxvh)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call PUSH1GL(part(1,npo),fx,sctx,qbm,dt,ek,idimp,npl,nx,nxh,nxvh)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
