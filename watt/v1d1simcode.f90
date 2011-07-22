module nrtype
   integer,parameter::I4B = selected_int_kind(9)
   integer,parameter::I2B = selected_int_kind(4)
   integer,parameter::I1B = selected_int_kind(2)
   integer,parameter::SP = kind(1.0)
   integer,parameter::DP = kind(1.0d0)
   integer,parameter::SPC = kind((1.0,1.0))
   integer,parameter::DPC = kind((1.0d0,1.0d0))
   integer,parameter::LGT = kind(.true.)
end module nrtype
!
!********************************************************************
!  version 1.00                             09/11/99
!
module parameters
!
!  To be "included" in driver programs, in order to provide
!  the values of physical constants and the values of unit numbers 
!
!********************************************************************
!
 use nrtype
 implicit none
!
 real(DP),parameter::pi=3.14159265_dp
 real(DP),parameter::ev=1.60217646e-19_dp
 real(DP),parameter::emas=9.10938188e-31_dp
 real(DP),parameter::hmas=1836.0*emas
 real(DP),parameter::hemas=4.0*hmas
 real(DP),parameter::omas=16.0*hmas
 real(DP),parameter::e0=8.85418782e-12_dp
 real(DP),parameter::mu0=4.0e0_dp*pi*1.0e-7_dp
 real(DP),parameter::clight=2.99792458e8_dp
 real(DP),parameter::kb=1.3806503e-23_dp
 real(DP),parameter::rad=pi/180.0e0
 real(DP),parameter::qe=-ev
 real(DP),parameter::qi=ev
 real(DP),parameter::one = 1.0e0_dp
 real(DP),parameter::half = 0.5e0_dp
 real(DP),parameter::zero = 0.0e0_dp
 real(DP),parameter::two = 2.0e0_dp
 real(DP),parameter::Re=6.371d+06
 real(DP),parameter::Mdip=8.05d+22
 real(DP),parameter::Gnewton=6.67428d-11
 real(DP),parameter::Mearth=5.9742d+24
 real(DP),parameter::Omega_earth=7.292115d-05
!
 integer(I4B),parameter::irl=5		        ! read from standard input
 integer(I4B),parameter::iwl=6		        ! write to standard output	
 integer(I4B),parameter::iwout=7		! unit number for .out file
 integer(I4B),parameter::irin=8		        ! unit number for .dat file
 integer(I4B),parameter::iw3=33
 integer(I4B),parameter::iw4=34
 integer(I4B),parameter::iw5=35
 integer(I4B),parameter::iw6=36
 integer(I4B),parameter::iw7=37
 integer(I4B),parameter::iw8=38
 integer(I4B),parameter::iw9=39
 integer(I4B),parameter::iw10=40
 integer(I4B),parameter::iw11=41
 integer(I4B),parameter::iw12=42
 integer(I4B),parameter::iw13=43
 integer(I4B),parameter::iw14=44
 integer(I4B),parameter::iw15=45
 integer(I4B),parameter::iw16=46
 integer(I4B),parameter::iw17=47
 integer(I4B),parameter::iw18=48
 integer(I4B),parameter::iw19=49
 integer(I4B),parameter::iw20=50
 integer(I4B),parameter::iw21=51
 integer(I4B),parameter::iw22=52
 integer(I4B),parameter::iw23=53
!
end module parameters
!
!***********************************************************************
!  Clare E. J. Watt                     26/01/00
!
program iaisim
!
!  DOUBLE-VLASOV simulation code: simulation program
!  Original electron code written by Richard B Horne in f77 
!  Modified by Mervyn P Freeman
!  Written by Clare E J Watt in f90 to include ions  
!
!  The second of two main driver programs which run a Vlasov simulation
!  in one dimension. The original Vlasov code has been modified to 
!  include ions. The first program sets up all the required files, and 
!  the second performs the simulation from the initial data in the 
!  files. This is convenient for continuing from previous simulations
!  which have been stopped.
!
!  This program reads all the required data from the root.out, 
!  root.coord.tmp and root.header.tmp files. Note that it differs from 
!  the original electron Vlasov code in that the grid dimensions are 
!  not given explicitly by the user, they are calculated from the plasma
!  information in the input file. This means that the setup file is used 
!  for setting up the grid, and the initial distribution function etc is 
!  not calculated until the beginning of this program.
!
!  This version of the code separates the Ampere Law into two, using only
!  the perturbed part to evolve the electric field. It is assumed that as
!  the spatially-averaged current evolves, the curl B term changes to 
!  balance it out.
!
!
!***********************************************************************
!
! module containing physical constants and module containing explicit
! interfaces.
!
 use nrtype; use parameters
!
 implicit none
!
! Input parameters
!
 integer(I4B)::ncomps                     ! number of components
 character(len=10),dimension(:),allocatable::form ! form of the distributions 
 integer(I4B),dimension(:),allocatable::species  ! species code for components
 real(DP),dimension(:),allocatable::mass  ! mass of each component
 real(DP),dimension(:),allocatable::nden  ! number density of components
 real(DP),dimension(:),allocatable::temp  ! temperature of the components
 real(DP),dimension(:),allocatable::drift ! drift velocity of the components
 real(DP)::mvth				  ! multiple of thermal speed for vcut
 real(DP)::cf				  ! courant factor
 integer(I4B)::ewrite,fwrite		  ! how often E/f are saved 
 real(DP)::enoise			  ! amplitude of internal E field
!
 integer(I4B)::ilen			  ! length of root
 character(len = 9)::filename	 	  ! filename of output file
 character(len = 8)::fil1
 character(len = 9)::fil2
 character(len = 13)::fil3,fil4,fil5,fil6,fil7,fil8,fil9,fil10
 character(len = 13)::fil11,fil12,fil13,fil14,fil15,fil16,fil17,fil18
 character(len = 16)::fil20
 character(len = 17)::fil21,fil22
 character(len = 12)::fil23
!
! distribution functions and grid values
!
 integer(I4B)::nx			 ! number of spatial grid points
 integer(I4B)::envx			 ! number of e velocity points
 integer(I4B)::invx			 ! number of i velocity points
 integer(I4B)::ntmax			 ! number of timesteps needed
 real(DP)::Lx				 ! length of simulation box
 real(DP)::dx				 ! spatial grid spacing
 real(DP)::dve				 ! e velocity grid spacing
 real(DP)::dvi				 ! i velocity grid spacing
 real(DP)::dt				 ! timestep
 real(DP),dimension(:),allocatable::x	 ! spatial grid points
 real(DP),dimension(:),allocatable::evx	 ! e velocity grid points
 real(DP),dimension(:),allocatable::ivx	 ! i velocity grid points
 real(DP),dimension(:,:),allocatable::fe ! electron distribution function
 real(DP),dimension(:,:),allocatable::fi ! ion distribution function
 real(DP),dimension(:),allocatable::fe0	 ! spatially-averaged e d.f.
 real(DP),dimension(:),allocatable::fi0	 ! spatially-averaged i d.f.
!
! fields and plasma parameters
!
 real(DP),dimension(:),allocatable::ex       ! electric field
 real(DP),dimension(:),allocatable::jx	     ! current density
 real(DP),dimension(:),allocatable::rho	     ! charge density
 real(DP),dimension(:),allocatable::nex	     ! electron number density
 real(DP),dimension(:),allocatable::uex	     ! electron drift velocity
 real(DP),dimension(:),allocatable::tex	     ! electron temperature
 real(DP),dimension(:),allocatable::nix	     ! ion number density
 real(DP),dimension(:),allocatable::uix	     ! ion drift velocity
 real(DP),dimension(:),allocatable::tix	     ! ion temperature
 real(DP),dimension(:),allocatable::exmin1   ! electric field at last timestep
 real(DP),dimension(:),allocatable::anres    ! anomalous resistivity
 real(DP),dimension(:),allocatable::effcol   ! effective collision frequency
 real(DP),dimension(:),allocatable::neve     ! ave n * ave u
 real(DP)::ne0				     ! average nex (for anres)
 real(DP),dimension(:),allocatable::chi	     ! full integral of fe+fi 
 real(DP)::chi0				     ! chi at t=0 
 real(DP),dimension(:),allocatable::chie     ! integral of fe over phase space
 real(DP),dimension(:),allocatable::chii     ! integral of fi over phase space
 real(DP),dimension(:),allocatable::S	     ! entropy
 real(DP),dimension(:),allocatable::aveS     ! "entropy" using f0
 real(DP)::S0
 real(DP)::phystime
 real(DP),dimension(:),allocatable::kine,kini
!
 real(DP)::kmin,kmax			     ! minimum and maximum k
 integer(I4B)::nkmodes			     ! number of growing k modes
 integer(I4B)::i,j,y			     ! count variables
 integer(I4B)::t			     ! timestep
 real(DP)::tnow				     ! current time (in s)
 character(len=4)::a			     ! for finding end of files
 real(DP)::average			     ! average value of things
 real(DP)::nu0,nu0min1
 real(DP)::ttest,btest
 real(DP)::err			
!
! first, set up the filenames
!
 write (iwl,*) ' Please enter input filename: '
 read(irl,'(a)') filename
 ilen = index(filename,'.out')
 if (ilen == 0) then
    write(iwl,*) 'Stopping, data file must have an extension of type .out'
    call zexit
 endif
 fil1(1:ilen)=filename(1:ilen)
 fil1(ilen:ilen+2)='.in'
 fil2(1:ilen)=filename(1:ilen)
 fil2(ilen:ilen+3)='.out'   	      
 fil3(1:ilen)=filename(1:ilen)
 fil3(ilen:ilen+7)='.fex.dat'
 fil4(1:ilen)=filename(1:ilen)
 fil4(ilen:ilen+7)='.fe0.dat'
 fil5(1:ilen)=filename(1:ilen)
 fil5(ilen:ilen+7)='.fix.dat'
 fil6(1:ilen) = filename(1:ilen)
 fil6(ilen:ilen+7) = '.fi0.dat'
 fil7(1:ilen)=filename(1:ilen)
 fil7(ilen:ilen+7) = '.nex.dat'
 fil8(1:ilen)=filename(1:ilen)
 fil8(ilen:ilen+7) = '.nix.dat' 
 fil9(1:ilen) = filename(1:ilen)
 fil9(ilen:ilen+7)='.tex.dat'
 fil10(1:ilen)=filename(1:ilen)
 fil10(ilen:ilen+7) = '.tix.dat'
 fil11(1:ilen)=filename(1:ilen)
 fil11(ilen:ilen+7) = '.uex.dat'
 fil12(1:ilen)=filename(1:ilen)
 fil12(ilen:ilen+7) = '.uix.dat'
 fil13(1:ilen)=filename(1:ilen)
 fil13(ilen:ilen+7) = '.rho.dat'
 fil14(1:ilen)=filename(1:ilen)
 fil14(ilen:ilen+7) = '.jxx.dat'
 fil15(1:ilen)=filename(1:ilen)
 fil15(ilen:ilen+7)='.exx.dat'
 fil16(1:ilen)=filename(1:ilen)
 fil16(ilen:ilen+7) = '.res.dat' 
 fil17(1:ilen)=filename(1:ilen)
 fil17(ilen:ilen+7) = '.dia.dat' 
!
! now all of these files need to be opened ready for reading/writing
!
 open (unit=irin,file=fil1,action='read',form='formatted')
 open (unit=iwout,file=fil2,action='readwrite',form='formatted')
 open (unit=iw3,file=fil3,action='write',form='formatted')
 open (unit=iw4,file=fil4,action='write',form='formatted')
 open (unit=iw5,file=fil5,action='write',form='formatted')
 open (unit=iw6,file=fil6,action='write',form='formatted')
 open (unit=iw7,file=fil7,action='write',form='formatted')
 open (unit=iw8,file=fil8,action='write',form='formatted')
 open (unit=iw9,file=fil9,action='write',form='formatted')
 open (unit=iw10,file=fil10,action='write',form='formatted')
 open (unit=iw11,file=fil11,action='write',form='formatted')
 open (unit=iw12,file=fil12,action='write',form='formatted')
 open (unit=iw13,file=fil13,action='write',form='formatted')
 open (unit=iw14,file=fil14,action='write',form='formatted')
 open (unit=iw15,file=fil15,action='write',form='formatted')
 open (unit=iw16,file=fil16,action='write',form='formatted')
 open (unit=iw17,file=fil17,action='write',form='formatted')
!
! Read in data from input file
!
 read (irin, *) ncomps
!
! Allocate files required for component information
!
 allocate(form(ncomps))
 allocate(species(ncomps))
 allocate(mass(ncomps))
 allocate(nden(ncomps))
 allocate(temp(ncomps))
 allocate(drift(ncomps))!
!
 do y=1,ncomps
    read (irin, *) species(y)
    read (irin, fmt='(a10)') form(y)
    read (irin, *) nden(y)
    read (irin, *) temp(y)
    read (irin, *) drift(y)
 enddo
 do y=1,ncomps
    if (species(y)==1) then
       mass(y) = emas
    else
       mass(y) = hmas
    endif
 enddo
 read (irin, *) mvth
 read (irin, *) cf
 read (irin, *) ntmax
 read (irin, *) fwrite
 read (irin, *) ewrite
 t = 1
!
! Is this an electron-only run, or is there an active proton population too?
!
 if (any(species.eq.2)) then
    read (iwout,'(/////////////// )')
    read (iwout,'(e11.4)') enoise
    read (iwout,'(/////////)')
    read (iwout,'(e11.4)') kmin
    read (iwout,'(e11.4)') kmax
    read (iwout,'(///////////)')
    read (iwout,'(f11.4)') dx
    read (iwout,'(e11.4)') dve
    read (iwout,'(e11.4)') dvi
    read (iwout,'(e11.4)') dt
    read (iwout,'(i6)') nx
    read (iwout,'(i6)') envx
    read (iwout,'(i6)') invx
    read (iwout,'(i7)') ntmax
    read (iwout,'(e11.4)') phystime
    write (iwout, '(i7)') t
    print*, 'E amplitude: ', enoise
    print*, 'min k in box: ', kmin
    print*, 'max k in box: ', kmax
    print*, 'number of spatial points: ', nx
    print*, 'number of electron velocity points: ', envx
    print*, 'number of ion velocity points: ', invx
    nkmodes = nint(kmax/kmin)
    Lx = nx*dx
!
! allocate arrays
!
    allocate(x(0:nx))
    allocate(evx(-envx:envx))
    allocate(ivx(-invx:invx))
    allocate(fe(-1:nx+1,-envx:envx))
    allocate(fe0(-envx:envx))
    allocate(fi(-1:nx+1,-invx:invx))
    allocate(fi0(-invx:invx))
    allocate(ex(0:nx))
    allocate(rho(0:nx))
    allocate(jx(0:nx))
    allocate(nex(0:nx))
    allocate(uex(0:nx))
    allocate(tex(0:nx))
    allocate(nix(0:nx))
    allocate(uix(0:nx))
    allocate(tix(0:nx))
    allocate(exmin1(0:nx))
    allocate(chi(ntmax))
    allocate(chie(ntmax))
    allocate(chii(ntmax))
    allocate(S(ntmax))
    allocate(aveS(ntmax))
    allocate(kine(ntmax))
    allocate(kini(ntmax))
    allocate(neve(ntmax))
    allocate(effcol(ntmax))
    allocate(anres(ntmax))
!
! calculate coords instead of reading them in, more accurate
!
    do i = 0,nx
       x(i) = i*dx
    enddo
    do j = -envx,envx
       evx(j) = j*dve
    enddo
    do j = -invx,invx
       ivx(j) = j*dvi
    enddo 
!
! Assign distribution functions and electrostatic noise at t=1
!
    call initf(ncomps,species,form,nden,temp,drift,nkmodes,enoise,Lx,dx,nx, &
&              envx,invx,fe,fi,x,evx,ivx,nex,uex,ex)
!
! Use distribution functions to get initial plasma parameters by calculating 
! moments.
!
    call inpairs(nx,envx,invx,fe,fi,evx,ivx,dve,dvi,nex,uex,tex, &
&                nix,uix,tix,rho,jx,chie(1),chii(1))
    chi0 = chie(1) + chii(1)
!
! calculate entropy and other simulation checks at t=1
!
    call entropy(nx,envx,invx,fe,fi,evx,ivx,dx,dve,dvi,S(1),aveS(1), &
&                kine(1),kini(1),neve(1),ne0)
    S0 = S(1)
    anres(1) = zero
    effcol(1) = zero
    chi(1) = zero
    S(1) = S(1)-S0
    nu0 = zero
    nu0min1 = sum(nex*(uex-uix))/(nx+1)
!
! write these initial values to the appropriate data files;
! 
    call wrdata(nx,envx,invx,x,evx,ivx,fe,fi,nex,uex,tex,nix,uix,tix, &
&               rho,jx,ex,kine(1),kini(1),chi(1),S(1),effcol(1),anres(1), &
&               zero,1,ewrite,fwrite,dt)
!
! start time integration loop. 
! mccorm is the routine which performs the time integration for 
! both the distributions and the electric field. 
! f1check makes sure that the perturbations at the velocity cutoffs 
! are not too large.
! inpairs calculates the moments of the distribution function
! all other lines of code are mainly concerned with calculating the 
! effective collision freqency and the anomlaous resistivity. 
!
    exmin1 = zero
    do t = 2,ntmax
       tnow = (t-1)*dt
       call mccorm(nx,envx,invx,fe,fi,evx,ivx,dx,dve,dvi,ex,exmin1, &
&                  jx,dt,t)
       call f1check(nx,envx,invx,fe,fi)
       call inpairs(nx,envx,invx,fe,fi,evx,ivx,dve,dvi,nex,uex,tex, &
&                   nix,uix,tix,rho,jx,chie(t),chii(t))
       nu0 = sum(nex*(uex-uix))/(nx+1)    
       chi(t) = chie(t) + chii(t) - chi0
       call entropy(nx,envx,invx,fe,fi,evx,ivx,dx,dve,dvi,S(t),aveS(t), &
&                   kine(t),kini(t),neve(t),ne0)
       effcol(t) = - (nu0-nu0min1)/(nu0*dt)
       anres(t) = emas*effcol(t)/(ne0*ev*ev)
       S(t) = S(t) - S0
       if (mod(t,ewrite)==0) then
          call wrdata(nx,envx,invx,x,evx,ivx,fe,fi,nex,uex,tex,nix,uix,tix, &
&                  rho,jx,ex,kine(t),kini(t),chi(t),S(t),effcol(t),anres(t), &
&                  tnow,t,ewrite,fwrite,dt)
       endif
       if (mod(t,ewrite)==0) then
          print*, t
       endif 
       nu0min1 = nu0
    enddo
!
 else
!
! This option is for an electron only plasma
!
    read (iwout,'(///////////////)')
    read (iwout,'(e11.4)') enoise
    read (iwout,'(/////////)')
    read (iwout,'(e11.4)') kmin
    read (iwout,'(e11.4)') kmax
    read (iwout,'(//////////)')
    read (iwout,'(f11.4)') dx
    read (iwout,'(e11.4)') dve
    read (iwout,'(e11.4)') dt
    read (iwout,'(i6)') nx
    read (iwout,'(i6)') envx
    read (iwout,'(i7)') ntmax
    read (iwout,'(e11.4)') phystime
    write (iwout, '(i7)') t
    print*, 'E amplitude: ', enoise
    print*, 'min k in box: ', kmin
    print*, 'max k in box: ', kmax
    print*, 'number of spatial points: ', nx
    print*, 'number of electron velocity points: ', envx
    nkmodes = nint(kmax/kmin)
    Lx = nx*dx
!
! allocate arrays
!
    allocate(x(0:nx))
    allocate(evx(-envx:envx))
    allocate(fe(-1:nx+1,-envx:envx))
    allocate(fe0(-envx:envx))
    allocate(ex(0:nx))
    allocate(rho(0:nx))
    allocate(jx(0:nx))
    allocate(nex(0:nx))
    allocate(uex(0:nx))
    allocate(uix(0:nx))
    allocate(tex(0:nx))
    allocate(exmin1(0:nx))
    allocate(chi(ntmax))
    allocate(chie(ntmax))
    allocate(S(ntmax))
    allocate(aveS(ntmax))
    allocate(kine(ntmax))
    allocate(neve(ntmax))
    allocate(effcol(ntmax))
    allocate(anres(ntmax))
!
! calculate coords instead of reading them in, more accurate
!
    do i = 0,nx
       x(i) = i*dx
    enddo
    do j = -envx,envx
       evx(j) = j*dve
    enddo
!
! Assign distribution functions and electrostatic noise at t=1
!
    call initf_e(ncomps,species,form,nden,temp,drift,nkmodes,enoise,Lx,dx,nx, &
&                envx,fe,x,evx,nex,uex,ex)
    uix = zero
!
! Use distribution functions to get initial plasma parameters by calculating 
! moments.
!
    call inpairs_e(ncomps,nden,nx,envx,fe,evx,dve,nex,uex,tex, &
&                  rho,jx,chie(1))
    chi0 = chie(1)
!
! calculate entropy and other simulation checks at t=1
!
    call entropy_e(nx,envx,fe,evx,dx,dve,S(1),aveS(1),kine(1),neve(1),ne0)
    S0 = S(1)
    anres(1) = zero
    effcol(1) = zero
    chi(1) = zero
    S(1) = S(1)-S0
    nu0 = zero
    nu0min1 = sum(nex*(uex-uix))/(nx+1)
!
! write these initial values to the appropriate data files;
! 
    call wrdata_e(nx,envx,x,evx,fe,nex,uex,tex,rho,jx,ex,kine(1),chi(1),S(1), &
&                 effcol(1),anres(1),zero,1,ewrite,fwrite,dt)
!
! start time integration loop. 
! mccorm is the routine which performs the time integration for 
! both the distributions and the electric field. 
! f1check makes sure that the perturbations at the velocity cutoffs 
! are not too large.
! inpairs calculates the moments of the distribution function
! all other lines of code are mainly concerned with calculating the 
! effective collision freqency and the anomlaous resistivity. 
!
    exmin1 = zero
    do t = 2,ntmax
       tnow = (t-1)*dt
       call mccorm_e(nx,envx,fe,evx,dx,dve,ex,exmin1,jx,dt,t)
       call f1check_e(nx,envx,fe)
       call inpairs_e(ncomps,nden,nx,envx,fe,evx,dve,nex,uex,tex, &
&                     rho,jx,chie(t))
       nu0 = sum(nex*(uex-uix))/(nx+1)    
       chi(t) = chie(t) - chi0
       call entropy_e(nx,envx,fe,evx,dx,dve,S(t),aveS(t),kine(t),neve(t),ne0)
       effcol(t) = - (nu0-nu0min1)/(nu0*dt)
       anres(t) = emas*effcol(t)/(ne0*ev*ev)
       S(t) = S(t) - S0
       if (mod(t,ewrite)==0) then
          call wrdata_e(nx,envx,x,evx,fe,nex,uex,tex,rho,jx,ex,kine(t), &
&                       chi(t),S(t),effcol(t),anres(t),tnow,t,ewrite,fwrite,dt)
          print*, t
       endif
       nu0min1 = nu0
    enddo
 endif
!
! write to root.out file that run has been successfully completed
!
 write (iwout,*) 'Run has been completed successfully'
!
! Deallocate all allocatable variables
!
 deallocate(x)
 deallocate(evx)
 deallocate(fe)
 deallocate(fe0)
 deallocate(ex)
 deallocate(rho)
 deallocate(jx)
 deallocate(nex)
 deallocate(uex)
 deallocate(tex)
 deallocate(uix)
 deallocate(exmin1)
 deallocate(chi)
 deallocate(chie)
 deallocate(S)
 deallocate(aveS)
 deallocate(kine)
 deallocate(neve)
 deallocate(effcol)
 deallocate(anres)
 deallocate(form)
 deallocate(species)
 deallocate(mass)
 deallocate(nden)
 deallocate(temp)
 deallocate(drift)
!
! close all files
!
 close(iwout) 
 close(iw3); close(iw4); close(iw5); close(iw6); close(iw7)
 close(iw8); close(iw9); close(iw10); close(iw11); close(iw12)
 close(iw13); close(iw14); close(iw15); close(iw16); close(iw17) 
end program iaisim
!
!
!****************************************************************************
!  Clare E. J. Watt                      05/01/00
!
subroutine initf(ncomps,species,form,nden,temp,drift,nk,enoise,Lx,dx,nx, &
&                envx,invx,fe,fi,x,evx,ivx,nex,uex,ex) 
!
!  Routine to initialise the distribution functions at t=0 
!  according to the information in the input file. 
!  Distribution functions are Maxwellian. This routine also forms
!  the noisy electric field at t=0 and modifies the number density
!  and distribution functions accordingly.
!
!****************************************************************************
!
 use nrtype; use parameters
! 
 implicit none
!
! declare input and output arguments
!
 integer(I4B)::ncomps                   
 character(len=10),dimension(ncomps),intent(in)::form 
 integer(I4B),dimension(ncomps),intent(in)::species  
 real(DP),dimension(ncomps),intent(in)::nden
 real(DP),dimension(ncomps),intent(in)::temp
 real(DP),dimension(ncomps),intent(in)::drift  
 integer(I4B),intent(in)::nk	
 real(DP),intent(in)::enoise	
 real(DP),intent(in)::Lx
 real(DP),intent(in)::dx
 integer(I4B),intent(in)::nx
 integer(I4B),intent(in)::envx
 integer(I4B),intent(in)::invx
!
! the electron and ion distribution functions...
! 
 real(DP),dimension(-1:nx+1,-envx:envx),intent(out)::fe
 real(DP),dimension(-1:nx+1,-invx:invx),intent(out)::fi
!
 real(DP),dimension(0:nx),intent(in)::x		! spatial grid points 
 real(DP),dimension(-envx:envx),intent(in)::evx	! electron vel. grid points
 real(DP),dimension(-invx:invx),intent(in)::ivx	! ion vel. grid points
 real(DP),dimension(0:nx),intent(out)::nex	! electron number density
 real(DP),dimension(0:nx),intent(out)::uex	! electron drift velocity
 real(DP),dimension(0:nx),intent(out)::ex	! electric field
!
! variables needed for the subroutine
!
 real(DP),dimension(:),allocatable::alpha       ! thermal velocities
 real(DP),dimension(:),allocatable::mass        ! species mass
 real(DP),dimension(:),allocatable::nu          ! nden fraction
 real(DP)::enden,inden                          ! total ion/electron number den
 real(DP)::edrift,idrift                        ! largest i/e drift
 real(DP)::rtpi					! square root of pi
 integer(I4B)::h,i,j,y				! count variables
 real(DP)::kmin				        ! minimum value of k
 real(DP)::ki,wi				! multiple of k
 integer(I4B)::nmax				! max number of k-values
 real(DP)::phi					! random phase
 real(DP)::rand				        ! seed for the ran routine
 real(DP)::debye2,cs2,wpe2,wpi2			! Debye length, squared
 real(DP)::ampe				        ! amplitude of e spectrum
 integer(I4B)::ierr				! error for fourier subr.
!
! initialise distribution function and electric field by setting null values
!
 allocate(alpha(ncomps))
 allocate(mass(ncomps))
 allocate(nu(ncomps))
 fe = zero
 fi = zero
 ex = zero
 do y=1,ncomps
    if (species(y)==1) then
       mass(y) = emas
    else
       mass(y) = hmas
    endif
 enddo 
 alpha = sqrt(two*temp*ev/mass)
 enden = zero
 inden = zero
 do y=1,ncomps
    if (species(y)==1) then
       enden = enden + nden(y)
    else
       inden = inden + nden(y)
    endif
 enddo
 do y=1,ncomps
    if (species(y)==1) then
       nu(y) = nden(y)/enden
    else
       nu(y) = nden(y)/inden
    endif
 enddo 
 edrift = zero
 idrift = zero
 do y=1,ncomps
    if ((species(y)==1).and.(drift(y).gt.edrift)) then
       edrift = edrift+drift(y)
    endif
    if ((species(y)==2).and.(drift(y).gt.idrift)) then
       idrift = idrift+drift(y)
    endif 
 enddo
 nex = enden
 uex = edrift
 rtpi = sqrt(pi)
!
! calculate noisy electric field, and use to get the number density 
! 
 kmin = (two*pi)/Lx
 nmax = Lx/(two*dx)
 j = 1
 do i = 1,nmax
    ki = i*kmin
    call random_number(rand)
    phi = rand*two*pi
    if (i.lt.nk) then
       ampe = enoise
    else
       ampe = zero
    endif
    ex = ex + ampe*sin(ki*x+phi) 
    nex = nex + e0*ki*ampe*cos(ki*x+phi)/qe
 enddo
!
! In order to avoid discontinuities at the boundary, here at the
! first timestep we explicitly set ex to be the same at each spatial 
! boundary, the distribution function is set to be periodic further on.
! 
 ex(nx) = ex(0)
!
! calculate distribution functions
!
 if (any(form.ne.'maxwellian')) then
    print*, 'STOPPING: Distribution functions not maxwellian'
    call zexit
 endif
 do y=1,ncomps
    if (species(y).eq.1) then
       do i=0,nx-1
          do j=-envx,envx
             fe(i,j) = fe(i,j) + nex(i)*nu(y)* &
&                      exp(-((evx(j)-drift(y))/alpha(y))**2)/(rtpi*alpha(y))
          enddo
       enddo
    else
       do i=0,nx-1
          do j=-invx,invx
             fi(i,j) = fi(i,j) + nden(y)* &
&                      exp(-((ivx(j)-drift(y))/alpha(y))**2)/(rtpi*alpha(y))
          enddo
       enddo
    endif
 enddo
 fi(nx,-invx:invx) = fi(0,-invx:invx) 
 fe(nx,-envx:envx) = fe(0,-envx:envx) 
!
! assign the grid points 'outside the box'.
!
 fe(-1,-envx:envx) = fe(nx-1,-envx:envx)
 fe(nx+1,-envx:envx) = fe(1,-envx:envx)
 fi(-1,-invx:invx) = fi(nx-1,-invx:invx)
 fi(nx+1,-invx:invx) = fi(1,-invx:invx)
!
 deallocate(alpha)
 deallocate(mass)
 deallocate(nu)
!
end subroutine initf                 
!
!
!
!****************************************************************************
!  Clare E. J. Watt                      05/01/00
!
subroutine initf_e(ncomps,species,form,nden,temp,drift,nk,enoise,Lx,dx,nx, &
&                envx,fe,x,evx,nex,uex,ex) 
!
!  Routine to initialise the distribution functions at t=0 
!  according to the information in the input file. 
!  Distribution functions are Maxwellian. This routine also forms
!  the noisy electric field at t=0 and modifies the number density
!  and distribution functions accordingly.
!
!****************************************************************************
!
 use nrtype; use parameters
! 
 implicit none
!
! declare input and output arguments
!
 integer(I4B)::ncomps                   
 character(len=10),dimension(ncomps),intent(in)::form 
 integer(I4B),dimension(ncomps),intent(in)::species  
 real(DP),dimension(ncomps),intent(in)::nden
 real(DP),dimension(ncomps),intent(in)::temp
 real(DP),dimension(ncomps),intent(in)::drift  
 integer(I4B),intent(in)::nk	
 real(DP),intent(in)::enoise	
 real(DP),intent(in)::Lx
 real(DP),intent(in)::dx
 integer(I4B),intent(in)::nx
 integer(I4B),intent(in)::envx
!
! the electron and ion distribution functions...
! 
 real(DP),dimension(-1:nx+1,-envx:envx),intent(out)::fe
!
 real(DP),dimension(0:nx),intent(in)::x		! spatial grid points 
 real(DP),dimension(-envx:envx),intent(in)::evx	! electron vel. grid points
 real(DP),dimension(0:nx),intent(out)::nex	! electron number density
 real(DP),dimension(0:nx),intent(out)::uex	! electron drift velocity
 real(DP),dimension(0:nx),intent(out)::ex	! electric field
!
! variables needed for the subroutine
!
 real(DP),dimension(:),allocatable::alpha       ! thermal velocities
 real(DP),dimension(:),allocatable::mass        ! species mass
 real(DP),dimension(:),allocatable::nu          ! nden fraction
 real(DP)::enden                                ! total electron number den
 real(DP)::edrift                               ! largest e drift
 real(DP)::rtpi					! square root of pi
 integer(I4B)::h,i,j,y				! count variables
 real(DP)::kmin				        ! minimum value of k
 real(DP)::ki,wi				! multiple of k
 integer(I4B)::nmax				! max number of k-values
 real(DP)::phi					! random phase
 real(DP)::rand				        ! seed for the ran routine
 real(DP)::ampe				        ! amplitude of e spectrum
 integer(I4B)::ierr				! error for fourier subr.
!
! initialise distribution function and electric field by setting null values
!
 allocate(alpha(ncomps))
 allocate(mass(ncomps))
 allocate(nu(ncomps))
 fe = zero
 ex = zero
 mass = emas
 alpha = sqrt(two*temp*ev/mass)
 enden = zero
 enden = sum(nden)
 nu = nden/enden
 edrift = sum(drift)
 nex = enden
 uex = edrift
 rtpi = sqrt(pi)
!
! calculate noisy electric field, and use to get the number density 
! 
 kmin = (two*pi)/Lx
 nmax = Lx/(two*dx)
 j = 1
 do i = 1,nmax
    ki = i*kmin
    call random_number(rand)
    phi = rand*two*pi
    if (i.lt.nk) then
       ampe = enoise
    else
       ampe = zero
    endif
    ex = ex + ampe*sin(ki*x+phi) 
    nex = nex + e0*ki*ampe*cos(ki*x+phi)/qe
 enddo
!
! In order to avoid discontinuities at the boundary, here at the
! first timestep we explicitly set ex to be the same at each spatial 
! boundary, the distribution function is set to be periodic further on.
! 
 ex(nx) = ex(0)
!
! calculate distribution functions
!
 if (any(form.ne.'maxwellian')) then
    print*, 'STOPPING: Distribution functions not maxwellian'
    call zexit
 endif
 do y=1,ncomps
    do i=0,nx-1
       do j=-envx,envx
          fe(i,j) = fe(i,j) + nex(i)*nu(y)* &
&                      exp(-((evx(j)-drift(y))/alpha(y))**2)/(rtpi*alpha(y))
       enddo
    enddo
 enddo
 fe(nx,-envx:envx) = fe(0,-envx:envx) 
!
! assign the grid points 'outside the box'.
!
 fe(-1,-envx:envx) = fe(nx-1,-envx:envx)
 fe(nx+1,-envx:envx) = fe(1,-envx:envx)
!
 deallocate(alpha)
 deallocate(mass)
 deallocate(nu)
!
end subroutine initf_e
!
!****************************************************************************
!  Clare E. J. Watt	                    26/01/00
!
subroutine inpairs_e(ncomps,nden,nx,envx,fe,evx,dve,nex,uex,tex, &
&                    rho,jx,chie)
!
!  Routine which calculates the moments of the two  
!  distribution functions in order to get plasma properties. 
!  The velocity integrals are summed in pairs around vx=0.
!
!****************************************************************************
!
 use nrtype; use parameters
!
 implicit none
!
! declare input and ouptut arguments
!
 integer(I4B),intent(in)::ncomps        ! number of components
 real(DP),dimension(ncomps),intent(in)::nden    ! number density
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of electron velocity points
!
! electron distribution functions..
!
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
!
! electron velocity grid points...
!
 real(DP),dimension(-envx:envx),intent(in)::evx
!
 real(DP),intent(in)::dve			! electron velocity grid spacing
 real(DP),dimension(0:nx),intent(out)::nex	! electron number density
 real(DP),dimension(0:nx),intent(out)::uex	! electron drift velocity
 real(DP),dimension(0:nx),intent(out)::tex	! electron temperature
 real(DP),dimension(0:nx),intent(out)::rho	! plasma charge density
 real(DP),dimension(0:nx),intent(out)::jx	! plasma current density
 real(DP),intent(out)::chie		        ! full integral of fe  
!
! interface for inpairsint function
!
 interface
    function inpairsint(nvx,f,dvx)
    use nrtype
    integer(I4B)::nvx    
    real(DP),dimension(-nvx:nvx)::f
    real(DP)::dvx
    real(DP)::inpairsint
    end function inpairsint
 end interface
!
! other variables needed for subroutine
!
 real(DP)::esum1,esum2,esum3			! for electron d.f. integrals
 real(DP),dimension(:),allocatable::fev0,fev,fev2   ! e functions for integrals
 integer(I4B)::i				! count variables
 real(DP)::inden                                ! ion number density
!
! initialise variables.
!
 inden = sum(nden)
 nex = zero
 uex = zero
 tex = zero
 rho = zero
 jx = zero
 chie = zero
!
! allocate functions for integration routines
!
 allocate(fev0(-envx:envx))
 allocate(fev(-envx:envx))
 allocate(fev2(-envx:envx))
!
! perform integration for each spatial grid point
!
 do i = 0,nx
    fev0 = fe(i,-envx:envx)
    fev = evx*fe(i,-envx:envx)
    fev2 = evx*evx*fe(i,-envx:envx)
    esum1 = inpairsint(envx,fev0,dve)
    esum2 = inpairsint(envx,fev,dve)
    esum3 = inpairsint(envx,fev2,dve)
    nex(i) = esum1
    rho(i) = qe*esum1 + qi*inden
    jx(i) = qe*esum2
    uex(i) = esum2/esum1
    tex(i) = (emas/ev)*(esum3/esum1 - uex(i)*uex(i))
    chie = chie + esum1
 enddo
!
!  deallocate arrays
! 
 deallocate(fev0)
 deallocate(fev)
 deallocate(fev2)
end subroutine inpairs_e
!
!
!****************************************************************************
!  Clare E. J. Watt	                    26/01/00
!
subroutine inpairs(nx,envx,invx,fe,fi,evx,ivx,dve,dvi,nex,uex,tex, &
&                  nix,uix,tix,rho,jx,chie,chii)
!
!  Routine which calculates the moments of the two  
!  distribution functions in order to get plasma properties. 
!  The velocity integrals are summed in pairs around vx=0.
!
!****************************************************************************
!
 use nrtype; use parameters
!
 implicit none
!
! declare input and ouptut arguments
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of electron velocity points
 integer(I4B),intent(in)::invx		! number of ion velocity points
!
! electron and ion distribution functions..
!
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
 real(DP),dimension(-1:nx+1,-invx:invx),intent(in)::fi	
!
! electron and ion velocity grid points...
!
 real(DP),dimension(-envx:envx),intent(in)::evx
 real(DP),dimension(-invx:invx),intent(in)::ivx
!
 real(DP),intent(in)::dve			! electron velocity grid spacing
 real(DP),intent(in)::dvi			! ion velocity grid spacing
 real(DP),dimension(0:nx),intent(out)::nex	! electron number density
 real(DP),dimension(0:nx),intent(out)::uex	! electron drift velocity
 real(DP),dimension(0:nx),intent(out)::tex	! electron temperature
 real(DP),dimension(0:nx),intent(out)::nix	! ion number density
 real(DP),dimension(0:nx),intent(out)::uix	! ion drift velocity
 real(DP),dimension(0:nx),intent(out)::tix	! ion temperature
 real(DP),dimension(0:nx),intent(out)::rho	! plasma charge density
 real(DP),dimension(0:nx),intent(out)::jx	! plasma current density
 real(DP),intent(out)::chie		        ! full integral of fe  
 real(DP),intent(out)::chii		        ! full integral of fi 
!
! interface for inpairsint function
!
 interface
    function inpairsint(nvx,f,dvx)
    use nrtype
    integer(I4B)::nvx    
    real(DP),dimension(-nvx:nvx)::f
    real(DP)::dvx
    real(DP)::inpairsint
    end function inpairsint
 end interface
!
! other variables needed for subroutine
!
 real(DP)::esum1,esum2,esum3			! for electron d.f. integrals
 real(DP)::isum1,isum2,isum3			! for ion d.f. integrals
 real(DP),dimension(:),allocatable::fev0,fev,fev2   ! e functions for integrals
 real(DP),dimension(:),allocatable::fiv0,fiv,fiv2   ! i functions for integrals
 integer(I4B)::i				! count variables
!
! initialise variables.
!
 nex = zero
 uex = zero
 tex = zero
 nix = zero
 uix = zero
 tix = zero
 rho = zero
 jx = zero
 chie = zero
 chii = zero
!
! allocate functions for integration routines
!
 allocate(fev0(-envx:envx))
 allocate(fev(-envx:envx))
 allocate(fev2(-envx:envx))
 allocate(fiv0(-invx:invx))
 allocate(fiv(-invx:invx))
 allocate(fiv2(-invx:invx))
!
! perform integration for each spatial grid point
!
 do i = 0,nx
    fev0 = fe(i,-envx:envx)
    fev = evx*fe(i,-envx:envx)
    fev2 = evx*evx*fe(i,-envx:envx)
    fiv0 = fi(i,-invx:invx)
    fiv = ivx*fi(i,-invx:invx)
    fiv2 = ivx*ivx*fi(i,-invx:invx)
    esum1 = inpairsint(envx,fev0,dve)
    esum2 = inpairsint(envx,fev,dve)
    esum3 = inpairsint(envx,fev2,dve)
    isum1 = inpairsint(invx,fiv0,dvi)
    isum2 = inpairsint(invx,fiv,dvi)
    isum3 = inpairsint(invx,fiv2,dvi)
    nex(i) = esum1
    nix(i) = isum1
    rho(i) = qe*esum1 + qi*isum1
    jx(i) = qe*esum2 + qi*isum2
    uex(i) = esum2/esum1
    uix(i) = isum2/isum1
    tex(i) = (emas/ev)*(esum3/esum1 - uex(i)*uex(i))
    tix(i) = (hmas/ev)*(isum3/isum1 - uix(i)*uix(i))
    chie = chie + esum1
    chii = chii + isum1
 enddo
 !
 ! deallocate arrays
 !
 deallocate(fev0)
 deallocate(fev)
 deallocate(fev2)
 deallocate(fiv0)
 deallocate(fiv)
 deallocate(fiv2)
end subroutine inpairs
!
!
!**************************************************************************
!  Clare E. J. Watt                                17/11/99
!
function inpairsint(nvx,f,dvx)
! 
!  Routine to perform the sum required to make the integral. The summing 
!  method is in pairs. This small calculation has been made into a separate 
!  function because it is needed in other routines other than the one used 
!  to calculate the moments. All are velocity integrals so they are 
!  calculated between -nvx and +nvx. 
!
!**************************************************************************
!
 use nrtype; use parameters
 implicit none
!
! declare input and output variables
!
 integer(I4B)::nvx		   	! Number of velocity grid points
 real(DP),dimension(-nvx:nvx)::f  	! function to be integrated
 real(DP)::dvx			        ! step between points in range
 real(DP)::inpairsint		        ! value of integral 
!
! variables needed for subroutine
!
 integer(I4B)::j,jp,jm			! count variables
 real(DP)::int
!
 int = f(0)
 do j = 1,nvx-1
    jp = j
    jm = -j
    int = int + f(jp) + f(jm)
 enddo
 int = int + half*(f(-nvx) + f(nvx))
 inpairsint = int*dvx
end function inpairsint
!
!
!**********************************************************************
!  Clare E. J. Watt		                 05/01/00
!
subroutine mccorm(nx,envx,invx,fe,fi,evx,ivx,dx,dve,dvi,ex,exmin1, &
&                  jx,dt,timestep) 
!
!  Routine to timestep forward the electric field and distribution
!  functions using McCormack's Method (see Computational Fluid Dynamics, 
!  J.F. Wendt, pg128). 
!
!********************************************************************** 
!
 use nrtype; use parameters
!
 implicit none
!
! declare input and output variables
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of e velocity points
 integer(I4B),intent(in)::invx		! number of ion velocity points 
!
!  electron/ion distribution functions, and velocity coordinates
!
 real(DP),dimension(-1:nx+1,-envx:envx),intent(inout)::fe
 real(DP),dimension(-1:nx+1,-invx:invx),intent(inout)::fi
 real(DP),dimension(-envx:envx),intent(in)::evx	
 real(DP),dimension(-invx:invx),intent(in)::ivx
 real(DP),intent(in)::dx			! spatial grid spacing
 real(DP),intent(in)::dve			! electron velocity grid spacing
 real(DP),intent(in)::dvi			! ion velocity grid spacing
!
! electric field and current (functions of x)
!
 real(DP),dimension(0:nx),intent(inout)::ex	! electric field
 real(DP),dimension(0:nx),intent(inout)::exmin1	! electric field at t-1	
 real(DP),dimension(0:nx),intent(in)::jx	! current
 real(DP)::avejx    				! average value of current density
 real(DP),intent(in)::dt			! length of timestep
 integer(I4B),intent(in)::timestep		! number of timestep
!
! variables needed for the subroutine
!
 integer(I4B)::i,j				! count variables
!
! derivatives for McCormack's method
!
 real(DP),dimension(:,:),allocatable::dfedx,dfedv,dfedt,dfedtbar
 real(DP),dimension(:,:),allocatable::febar
 real(DP),dimension(:,:),allocatable::avedfedt
 real(DP),dimension(:,:),allocatable::dfidx,dfidv,dfidt,dfidtbar
 real(DP),dimension(:,:),allocatable::fibar
 real(DP),dimension(:,:),allocatable::avedfidt
 real(DP),dimension(:),allocatable::dEdt,explus1
!
! interface for inpairsint function
!
 interface
    function inpairsint(nvx,f,dvx)
    use nrtype
    integer(I4B)::nvx    
    real(DP),dimension(-nvx:nvx)::f
    real(DP)::dvx
    real(DP)::inpairsint
    end function inpairsint
 end interface
!
! allocate all arrays
!
 allocate(dfedx(0:nx,-envx:envx))
 allocate(dfedv(0:nx,-envx:envx))
 allocate(dfedt(0:nx,-envx:envx))
 allocate(dfedtbar(0:nx,-envx:envx))
 allocate(febar(-1:nx+1,-envx:envx))
 allocate(avedfedt(0:nx,-envx:envx))
 allocate(dfidx(0:nx,-invx:invx))
 allocate(dfidv(0:nx,-invx:invx))
 allocate(dfidt(0:nx,-invx:invx))
 allocate(dfidtbar(0:nx,-invx:invx))
 allocate(fibar(-1:nx+1,-invx:invx))
 allocate(avedfidt(0:nx,-invx:invx))
 allocate(dEdt(0:nx))
 allocate(explus1(0:nx))
!
! step forward electric field. Note that because there is a non-zero
! current in the system, this must be removed in order to isolate 
! the wave-particle interactions
!
 avejx = sum(jx)/(nx+1)
 if (timestep == 2) then
    explus1 = - (clight*clight*mu0)*(jx-avejx)*dt + ex
 else
    explus1 = - two*(clight*clight*mu0)*(jx-avejx)*dt + exmin1
 endif
!
! If timestep is even, use forward difference for the predictor step,
! and rearward difference for the correction. If odd, then vice versa.
! This version of the mccormack algorithm relies on fe/fi having two more
! spatial grid points. This is to ensure that the boundary conditions 
! are dealt with correctly. The two extra grid points are set to their 
! equivalent at the other end of the spatial grid, and are used to 
! calculate the spatial derivatives for the "main" part of the grid.
!
 if (mod(timestep,2).eq.0) then
    do j = -envx,envx-1
       dfedv(0:nx,j) = (fe(0:nx,j+1)-fe(0:nx,j))/dve
    enddo  
    dfedv(0:nx,-envx) =-(25.0_dp*fe(0:nx,-envx) &
&           -48.0_dp*fe(0:nx,-envx+1)+36.0_dp*fe(0:nx,-envx+2) &
&           -16.0_dp*fe(0:nx,-envx+3)+3.0_dp*fe(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*fe(0:nx,envx)-48.0_dp*fe(0:nx,envx-1) &
&                 +36.0_dp*fe(0:nx,envx-2)-16.0_dp*fe(0:nx,envx-3) &
&                 +3.0_dp*fe(0:nx,envx-4))/(12.0_dp*dve)
    do j = -invx,invx-1
       dfidv(0:nx,j) = (fi(0:nx,j+1)-fi(0:nx,j))/dvi
    enddo  
    dfidv(0:nx,-invx) =-(25.0_dp*fi(0:nx,-invx) &
&           -48.0_dp*fi(0:nx,-invx+1)+36.0_dp*fi(0:nx,-invx+2) &
&           -16.0_dp*fi(0:nx,-invx+3)+3.0_dp*fi(0:nx,-invx+4))/(12.0_dp*dvi)
    dfidv(0:nx,invx) = (25.0_dp*fi(0:nx,invx)-48.0_dp*fi(0:nx,invx-1) &
&                 +36.0_dp*fi(0:nx,invx-2)-16.0_dp*fi(0:nx,invx-3) &
&                 +3.0_dp*fi(0:nx,invx-4))/(12.0_dp*dvi)
    do i = 0,nx
       dfedx(i,-envx:envx) = (fe(i+1,-envx:envx)-fe(i,-envx:envx))/dx
    enddo
    do i = 0,nx
       dfidx(i,-invx:invx) = (fi(i+1,-invx:invx)-fi(i,-invx:invx))/dx
    enddo
    do i = 0,nx
       do j = -envx,envx
          dfedt(i,j) = -evx(j)*dfedx(i,j) - (qe/emas)*ex(i)*dfedv(i,j)
          febar(i,j) = fe(i,j) + dfedt(i,j)*dt
       enddo
       do j = -invx,invx
          dfidt(i,j) = -ivx(j)*dfidx(i,j) - (qi/hmas)*ex(i)*dfidv(i,j)
          fibar(i,j) = fi(i,j) + dfidt(i,j)*dt
       enddo
    enddo
    febar(nx,:) = febar(0,:)
    febar(-1,:) = febar(nx-1,:)
    febar(nx+1,:) = febar(1,:)
    fibar(nx,:) = fibar(0,:)
    fibar(-1,:) = fibar(nx-1,:)
    fibar(nx+1,:) = fibar(1,:)
    do j = -envx+1,envx
       dfedv(0:nx,j) = (febar(0:nx,j)-febar(0:nx,j-1))/dve
    enddo 
    dfedv(0:nx,-envx) =-(25.0_dp*febar(0:nx,-envx) &
&       -48.0_dp*febar(0:nx,-envx+1)+36.0_dp*febar(0:nx,-envx+2) &
&       -16.0_dp*febar(0:nx,-envx+3)+3.0_dp*febar(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*febar(0:nx,envx)-48.0_dp*febar(0:nx,envx-1) &
&                 +36.0_dp*febar(0:nx,envx-2)-16.0_dp*febar(0:nx,envx-3) &
&                 +3.0_dp*febar(0:nx,envx-4))/(12.0_dp*dve)
    do j = -invx+1,invx
       dfidv(0:nx,j) = (fibar(0:nx,j)-fibar(0:nx,j-1))/dvi
    enddo 
    dfidv(0:nx,-invx) =-(25.0_dp*fibar(0:nx,-invx) &
&       -48.0_dp*fibar(0:nx,-invx+1)+36.0_dp*fibar(0:nx,-invx+2) &
&       -16.0_dp*fibar(0:nx,-invx+3)+3.0_dp*fibar(0:nx,-invx+4))/(12.0_dp*dvi)
    dfidv(0:nx,invx) = (25.0_dp*fibar(0:nx,invx)-48.0_dp*fibar(0:nx,invx-1) &
&                 +36.0_dp*fibar(0:nx,invx-2)-16.0_dp*fibar(0:nx,invx-3) &
&                 +3.0_dp*fibar(0:nx,invx-4))/(12.0_dp*dvi)
    do i = 0,nx
       dfedx(i,-envx:envx) = (febar(i,-envx:envx)-febar(i-1,-envx:envx))/dx
    enddo 
    do i = 0,nx
       dfidx(i,-invx:invx) = (fibar(i,-invx:invx)-fibar(i-1,-invx:invx))/dx
    enddo 
    do i =0,nx
       do j = -envx,envx
          dfedtbar(i,j) = -evx(j)*dfedx(i,j)-(qe/emas)*explus1(i)*dfedv(i,j)
          avedfedt(i,j) = half*(dfedt(i,j) + dfedtbar(i,j))
          fe(i,j) = fe(i,j)+avedfedt(i,j)*dt
       enddo
       do j = -invx,invx
          dfidtbar(i,j) = -ivx(j)*dfidx(i,j)-(qi/hmas)*explus1(i)*dfidv(i,j)
          avedfidt(i,j) = half*(dfidt(i,j) + dfidtbar(i,j))
          fi(i,j) = fi(i,j)+avedfidt(i,j)*dt
       enddo
    enddo
    fe(0,:) = fe(nx,:)
    fe(-1,:) = fe(nx-1,:)
    fe(nx+1,:) = fe(1,:)
    fi(0,:) = fi(nx,:)
    fi(-1,:) = fi(nx-1,:)
    fi(nx+1,:) = fi(1,:)
!
 else
!
! This is the case of an odd timestep
!
    do j = -envx+1,envx
       dfedv(0:nx,j) = (fe(0:nx,j)-fe(0:nx,j-1))/dve
    enddo
    dfedv(0:nx,-envx) =-(25.0_dp*fe(0:nx,-envx) &
&           -48.0_dp*fe(0:nx,-envx+1)+36.0_dp*fe(0:nx,-envx+2) &
&           -16.0_dp*fe(0:nx,-envx+3)+3.0_dp*fe(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*fe(0:nx,envx)-48.0_dp*fe(0:nx,envx-1) &
&                 +36.0_dp*fe(0:nx,envx-2)-16.0_dp*fe(0:nx,envx-3) &
&                 +3.0_dp*fe(0:nx,envx-4))/(12.0_dp*dve)
    do j = -invx+1,invx
       dfidv(0:nx,j) = (fi(0:nx,j)-fi(0:nx,j-1))/dvi
    enddo
    dfidv(0:nx,-invx) =-(25.0_dp*fi(0:nx,-invx) &
&           -48.0_dp*fi(0:nx,-invx+1)+36.0_dp*fi(0:nx,-invx+2) &
&           -16.0_dp*fi(0:nx,-invx+3)+3.0_dp*fi(0:nx,-invx+4))/(12.0_dp*dvi)
    dfidv(0:nx,invx) = (25.0_dp*fi(0:nx,invx)-48.0_dp*fi(0:nx,invx-1) &
&                 +36.0_dp*fi(0:nx,invx-2)-16.0_dp*fi(0:nx,invx-3) &
&                 +3.0_dp*fi(0:nx,invx-4))/(12.0_dp*dvi)
    do i = 0,nx
       dfedx(i,-envx:envx) = (fe(i,-envx:envx)-fe(i-1,-envx:envx))/dx
    enddo
    do i = 0,nx
       dfidx(i,-invx:invx) = (fi(i,-invx:invx)-fi(i-1,-invx:invx))/dx
    enddo
    do i = 0,nx
       do j = -envx,envx
          dfedt(i,j) = -evx(j)*dfedx(i,j) - (qe/emas)*ex(i)*dfedv(i,j)
          febar(i,j) = fe(i,j) + dfedt(i,j)*dt
       enddo
       do j = -invx,invx
          dfidt(i,j) = -ivx(j)*dfidx(i,j) - (qi/hmas)*ex(i)*dfidv(i,j)
          fibar(i,j) = fi(i,j) + dfidt(i,j)*dt
       enddo
    enddo
    febar(0,:) = febar(nx,:)
    febar(-1,:) = febar(nx-1,:)
    febar(nx+1,:) = febar(1,:)
    fibar(0,:) = fibar(nx,:)
    fibar(-1,:) = fibar(nx-1,:)
    fibar(nx+1,:) = fibar(1,:)
    do j = -envx,envx-1
       dfedv(0:nx,j) = (febar(0:nx,j+1)-febar(0:nx,j))/dve
    enddo
    dfedv(0:nx,-envx) =-(25.0_dp*febar(0:nx,-envx) &
&       -48.0_dp*febar(0:nx,-envx+1)+36.0_dp*febar(0:nx,-envx+2) &
&       -16.0_dp*febar(0:nx,-envx+3)+3.0_dp*febar(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*febar(0:nx,envx)-48.0_dp*febar(0:nx,envx-1) &
&                 +36.0_dp*febar(0:nx,envx-2)-16.0_dp*febar(0:nx,envx-3) &
&                 +3.0_dp*febar(0:nx,envx-4))/(12.0_dp*dve)
    do j = -invx,invx-1
       dfidv(0:nx,j) = (fibar(0:nx,j+1)-fibar(0:nx,j))/dvi
    enddo
    dfidv(0:nx,-invx) =-(25.0_dp*fibar(0:nx,-invx) &
&       -48.0_dp*fibar(0:nx,-invx+1)+36.0_dp*fibar(0:nx,-invx+2) &
&       -16.0_dp*fibar(0:nx,-invx+3)+3.0_dp*fibar(0:nx,-invx+4))/(12.0_dp*dvi)
    dfidv(0:nx,invx) = (25.0_dp*fibar(0:nx,invx)-48.0_dp*fibar(0:nx,invx-1) &
&                 +36.0_dp*fibar(0:nx,invx-2)-16.0_dp*fibar(0:nx,invx-3) &
&                 +3.0_dp*fibar(0:nx,invx-4))/(12.0_dp*dvi)
    do i = 0,nx
       dfedx(i,-envx:envx) = (febar(i+1,-envx:envx)-febar(i,-envx:envx))/dx
    enddo
    do i = 0,nx
       dfidx(i,-invx:invx) = (fibar(i+1,-invx:invx)-fibar(i,-invx:invx))/dx
    enddo
    do i = 0,nx
       do j = -envx,envx
          dfedtbar(i,j) = -evx(j)*dfedx(i,j)-(qe/emas)*explus1(i)*dfedv(i,j)
          avedfedt(i,j) = half*(dfedt(i,j) + dfedtbar(i,j))
          fe(i,j) = fe(i,j) + avedfedt(i,j)*dt
       enddo
       do j = -invx,invx
          dfidtbar(i,j) = -ivx(j)*dfidx(i,j)-(qi/hmas)*explus1(i)*dfidv(i,j)
          avedfidt(i,j) = half*(dfidt(i,j) + dfidtbar(i,j))
          fi(i,j) = fi(i,j) + avedfidt(i,j)*dt
       enddo
    enddo
    fe(nx,:) = fe(0,:)
    fe(-1,:) = fe(nx-1,:)
    fe(nx+1,:) = fe(1,:)
    fi(nx,:) = fi(0,:)
    fi(-1,:) = fi(nx-1,:)
    fi(nx+1,:) = fi(1,:)
 endif
!
! assign exmin1 and ex so that they are ready for next central
! difference calculation
!
 exmin1 = ex
 ex = explus1 
!
 deallocate(dfedx)
 deallocate(dfedv)
 deallocate(dfedt)
 deallocate(dfedtbar)
 deallocate(febar)
 deallocate(avedfedt)
 deallocate(dfidx)
 deallocate(dfidv)
 deallocate(dfidt)
 deallocate(dfidtbar)
 deallocate(fibar)
 deallocate(avedfidt)
 deallocate(dEdt)
 deallocate(explus1)
end subroutine mccorm  
!
!
!**********************************************************************
!  Clare E. J. Watt		                 05/01/00
!
subroutine mccorm_e(nx,envx,fe,evx,dx,dve,ex,exmin1,jx,dt,timestep) 
!
!  Routine to timestep forward the electric field and distribution
!  functions using McCormack's Method (see Computational Fluid Dynamics, 
!  J.F. Wendt, pg128). 
!
!********************************************************************** 
!
 use nrtype; use parameters
!
 implicit none
!
! declare input and output variables
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of e velocity points
!
!  electron/ion distribution functions, and velocity coordinates
!
 real(DP),dimension(-1:nx+1,-envx:envx),intent(inout)::fe
 real(DP),dimension(-envx:envx),intent(in)::evx	
 real(DP),intent(in)::dx			! spatial grid spacing
 real(DP),intent(in)::dve			! electron velocity grid spacing
!
! electric field and current (functions of x)
!
 real(DP),dimension(0:nx),intent(inout)::ex	! electric field
 real(DP),dimension(0:nx),intent(inout)::exmin1	! electric field at t-1	
 real(DP),dimension(0:nx),intent(in)::jx	! current
 real(DP)::avejx    				! average value of current dens.
 real(DP),intent(in)::dt			! length of timestep
 integer(I4B),intent(in)::timestep		! number of timestep
!
! variables needed for the subroutine
!
 integer(I4B)::i,j				! count variables
!
! derivatives for McCormack's method
!
 real(DP),dimension(:,:),allocatable::dfedx,dfedv,dfedt,dfedtbar
 real(DP),dimension(:,:),allocatable::febar
 real(DP),dimension(:,:),allocatable::avedfedt
 real(DP),dimension(:),allocatable::dEdt,explus1
!
! interface for inpairsint function
!
 interface
    function inpairsint(nvx,f,dvx)
    use nrtype
    integer(I4B)::nvx    
    real(DP),dimension(-nvx:nvx)::f
    real(DP)::dvx
    real(DP)::inpairsint
    end function inpairsint
 end interface
!
! allocate all arrays
!
 allocate(dfedx(0:nx,-envx:envx))
 allocate(dfedv(0:nx,-envx:envx))
 allocate(dfedt(0:nx,-envx:envx))
 allocate(dfedtbar(0:nx,-envx:envx))
 allocate(febar(-1:nx+1,-envx:envx))
 allocate(avedfedt(0:nx,-envx:envx))
 allocate(dEdt(0:nx))
 allocate(explus1(0:nx))
!
! step forward electric field. Note that because there is a non-zero
! current in the system, this must be removed in order to isolate 
! the wave-particle interactions
!
 avejx = sum(jx)/(nx+1)
 if (timestep == 2) then
    explus1 = - (clight*clight*mu0)*(jx-avejx)*dt + ex
 else
    explus1 = - two*(clight*clight*mu0)*(jx-avejx)*dt + exmin1
 endif
!
! If timestep is even, use forward difference for the predictor step,
! and rearward difference for the correction. If odd, then vice versa.
! This version of the mccormack algorithm relies on fe/fi having two more
! spatial grid points. This is to ensure that the boundary conditions 
! are dealt with correctly. The two extra grid points are set to their 
! equivalent at the other end of the spatial grid, and are used to 
! calculate the spatial derivatives for the "main" part of the grid.
!
 if (mod(timestep,2).eq.0) then
    do j = -envx,envx-1
       dfedv(0:nx,j) = (fe(0:nx,j+1)-fe(0:nx,j))/dve
    enddo  
    dfedv(0:nx,-envx) =-(25.0_dp*fe(0:nx,-envx) &
&           -48.0_dp*fe(0:nx,-envx+1)+36.0_dp*fe(0:nx,-envx+2) &
&           -16.0_dp*fe(0:nx,-envx+3)+3.0_dp*fe(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*fe(0:nx,envx)-48.0_dp*fe(0:nx,envx-1) &
&                 +36.0_dp*fe(0:nx,envx-2)-16.0_dp*fe(0:nx,envx-3) &
&                 +3.0_dp*fe(0:nx,envx-4))/(12.0_dp*dve)
    do i = 0,nx
       dfedx(i,-envx:envx) = (fe(i+1,-envx:envx)-fe(i,-envx:envx))/dx
    enddo
    do i = 0,nx
       do j = -envx,envx
          dfedt(i,j) = -evx(j)*dfedx(i,j) - (qe/emas)*ex(i)*dfedv(i,j)
          febar(i,j) = fe(i,j) + dfedt(i,j)*dt
       enddo
    enddo
    febar(nx,:) = febar(0,:)
    febar(-1,:) = febar(nx-1,:)
    febar(nx+1,:) = febar(1,:)
    do j = -envx+1,envx
       dfedv(0:nx,j) = (febar(0:nx,j)-febar(0:nx,j-1))/dve
    enddo 
    dfedv(0:nx,-envx) =-(25.0_dp*febar(0:nx,-envx) &
&       -48.0_dp*febar(0:nx,-envx+1)+36.0_dp*febar(0:nx,-envx+2) &
&       -16.0_dp*febar(0:nx,-envx+3)+3.0_dp*febar(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*febar(0:nx,envx)-48.0_dp*febar(0:nx,envx-1) &
&                 +36.0_dp*febar(0:nx,envx-2)-16.0_dp*febar(0:nx,envx-3) &
&                 +3.0_dp*febar(0:nx,envx-4))/(12.0_dp*dve)
    do i = 0,nx
       dfedx(i,-envx:envx) = (febar(i,-envx:envx)-febar(i-1,-envx:envx))/dx
    enddo 
    do i =0,nx
       do j = -envx,envx
          dfedtbar(i,j) = -evx(j)*dfedx(i,j)-(qe/emas)*explus1(i)*dfedv(i,j)
          avedfedt(i,j) = half*(dfedt(i,j) + dfedtbar(i,j))
          fe(i,j) = fe(i,j)+avedfedt(i,j)*dt
       enddo
    enddo
    fe(0,:) = fe(nx,:)
    fe(-1,:) = fe(nx-1,:)
    fe(nx+1,:) = fe(1,:)
!
 else
!
! This is the case of an odd timestep
!
    do j = -envx+1,envx
       dfedv(0:nx,j) = (fe(0:nx,j)-fe(0:nx,j-1))/dve
    enddo
    dfedv(0:nx,-envx) =-(25.0_dp*fe(0:nx,-envx) &
&           -48.0_dp*fe(0:nx,-envx+1)+36.0_dp*fe(0:nx,-envx+2) &
&           -16.0_dp*fe(0:nx,-envx+3)+3.0_dp*fe(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*fe(0:nx,envx)-48.0_dp*fe(0:nx,envx-1) &
&                 +36.0_dp*fe(0:nx,envx-2)-16.0_dp*fe(0:nx,envx-3) &
&                 +3.0_dp*fe(0:nx,envx-4))/(12.0_dp*dve)
    do i = 0,nx
       dfedx(i,-envx:envx) = (fe(i,-envx:envx)-fe(i-1,-envx:envx))/dx
    enddo
    do i = 0,nx
       do j = -envx,envx
          dfedt(i,j) = -evx(j)*dfedx(i,j) - (qe/emas)*ex(i)*dfedv(i,j)
          febar(i,j) = fe(i,j) + dfedt(i,j)*dt
       enddo
    enddo
    febar(0,:) = febar(nx,:)
    febar(-1,:) = febar(nx-1,:)
    febar(nx+1,:) = febar(1,:)
    do j = -envx,envx-1
       dfedv(0:nx,j) = (febar(0:nx,j+1)-febar(0:nx,j))/dve
    enddo
    dfedv(0:nx,-envx) =-(25.0_dp*febar(0:nx,-envx) &
&       -48.0_dp*febar(0:nx,-envx+1)+36.0_dp*febar(0:nx,-envx+2) &
&       -16.0_dp*febar(0:nx,-envx+3)+3.0_dp*febar(0:nx,-envx+4))/(12.0_dp*dve)
    dfedv(0:nx,envx) = (25.0_dp*febar(0:nx,envx)-48.0_dp*febar(0:nx,envx-1) &
&                 +36.0_dp*febar(0:nx,envx-2)-16.0_dp*febar(0:nx,envx-3) &
&                 +3.0_dp*febar(0:nx,envx-4))/(12.0_dp*dve)
    do i = 0,nx
       dfedx(i,-envx:envx) = (febar(i+1,-envx:envx)-febar(i,-envx:envx))/dx
    enddo
    do i = 0,nx
       do j = -envx,envx
          dfedtbar(i,j) = -evx(j)*dfedx(i,j)-(qe/emas)*explus1(i)*dfedv(i,j)
          avedfedt(i,j) = half*(dfedt(i,j) + dfedtbar(i,j))
          fe(i,j) = fe(i,j) + avedfedt(i,j)*dt
       enddo
    enddo
    fe(nx,:) = fe(0,:)
    fe(-1,:) = fe(nx-1,:)
    fe(nx+1,:) = fe(1,:)
 endif
!
! assign exmin1 and ex so that they are ready for next central
! difference calculation
!
 exmin1 = ex
 ex = explus1 
!
 deallocate(dfedx)
 deallocate(dfedv)
 deallocate(dfedt)
 deallocate(dfedtbar)
 deallocate(febar)
 deallocate(avedfedt)
 deallocate(dEdt)
 deallocate(explus1)
end subroutine mccorm_e  
!
!
!*******************************************************************
!   Clare E. J. Watt			14/02/00
!
 subroutine entropy(nx,envx,invx,fe,fi,evx,ivx,dx,dve,dvi,S,aveS, &
&                   kine,kini,neve,ne0)
!
!  Subroutine which calculates the entropy of the system using;
!	
!	S = - sum( int f * ln f dx dv ) 
!
!  where the sum is over all plasma species in the system. This 
!  subroutine also calculates the kinetic energy of the electrons
!  and the ions and the spatially averaged first moments of fe 
!  at each timestep. 
!
!*******************************************************************
!
 use nrtype; use parameters
 implicit none
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of e. velocity grid points
 integer(I4B),intent(in)::invx		! number of ion velocity grid points
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
 real(DP),dimension(-1:nx+1,-invx:invx),intent(in)::fi
 real(DP),dimension(-envx:envx)::evx	! electron velocity grid points
 real(DP),dimension(-invx:invx)::ivx	! ion velocity grid points
 real(DP),intent(in)::dx		! spatial grid spacing
 real(DP),intent(in)::dve		! electron velocity grid spacing
 real(DP),intent(in)::dvi		! ion velocity grid spacing
 real(DP),intent(out)::S		! entropy
 real(DP),intent(out)::aveS		! "entropy" using fe0
 real(DP),intent(out)::kine		! electron kinetic energy
 real(DP),intent(out)::kini		! ion kinetic energy
 real(DP),intent(out)::neve		! spatially-ave 1st moment
 real(DP),intent(out)::ne0		! spatially-ave 0th moment
!
 real(DP),dimension(:),allocatable::felnfe	! for calculating S
 real(DP),dimension(:),allocatable::filnfi	!   "   "    "  
 real(DP)::fesum				!   "   "    "
 real(DP)::fisum
 real(DP),dimension(:),allocatable::fe0,fi0	! spatial average of fe and fi
 real(DP),dimension(:),allocatable::fe0lnfe0    ! for calculating aveS
 real(DP),dimension(:),allocatable::fi0lnfi0    ! for calculating aveS
 real(DP),dimension(:,:),allocatable::fesmall,fismall
 real(DP),dimension(:),allocatable::fev	        ! v * fe
 real(DP),dimension(:),allocatable::fev2,fiv2   ! v^2 * fe, v^2 * fi
 real(DP)::fesum2,fisum2			! for calculating kine/i
 integer(I4B)::i				! counting variable
!
! interfaces
!
 interface
    function inpairsint(nvx,f,dvx)
    use nrtype
    integer(I4B)::nvx    
    real(DP),dimension(-nvx:nvx)::f
    real(DP)::dvx
    real(DP)::inpairsint
    end function inpairsint
 end interface
!
 interface
    function fave(n,nvx,f)
    use nrtype
    integer(I4B)::n,nvx
    real(DP),dimension(0:n,-nvx:nvx)::f
    real(DP),dimension(-nvx:nvx)::fave
    end function fave
 end interface
!
! allocate all internal subroutine arrays
!
 allocate(felnfe(-envx:envx))
 allocate(filnfi(-invx:invx))
 allocate(fe0(-envx:envx))
 allocate(fi0(-invx:invx))
 allocate(fe0lnfe0(-envx:envx))
 allocate(fi0lnfi0(-invx:invx))
 allocate(fesmall(0:nx,-envx:envx))
 allocate(fismall(0:nx,-invx:invx))
 allocate(fev2(-envx:envx))
 allocate(fiv2(-invx:invx))
 allocate(fev(-envx:envx))
!
! calculate the entropy of the system
!
 S = zero
 aveS = zero
 kine = zero
 kini = zero
 do i = 0,nx
    where (fe(i,-envx:envx) .gt. zero)
       felnfe = fe(i,-envx:envx)*log(fe(i,-envx:envx))
    elsewhere
       felnfe = zero
    endwhere
    fesum = inpairsint(envx,felnfe,dve)
    where (fi(i,-invx:invx) .gt. zero)
       filnfi = fi(i,-invx:invx)*log(fi(i,-invx:invx))
    elsewhere
       filnfi = zero
    endwhere
    fisum = inpairsint(invx,filnfi,dvi)
    S = S + fesum + fisum
 enddo
 S = -S*dx
!
! calculate the spatially averaged distribution in order to 
! calculate the "entropy" using fe0
!
 fesmall = fe(0:nx,-envx:envx)
 fe0 = fave(nx,envx,fesmall)
 where (fe0 .gt. zero) 
    fe0lnfe0 = fe0*log(fe0)
 elsewhere
    fe0lnfe0 = zero
 endwhere
 fismall = fi(0:nx,-invx:invx)
 fi0 = fave(nx,invx,fismall)
 where (fi0 .gt. zero)
    fi0lnfi0 = fi0*log(fi0)
 elsewhere
    fi0lnfi0 = zero
 endwhere 
 aveS = inpairsint(envx,fe0lnfe0,dve) + inpairsint(invx,fi0lnfi0,dvi)
 aveS = -aveS
!
! use fe0 to calculate the spatially averaged moments
!
 ne0 = inpairsint(envx,fe0,dve)
 fev = evx*fe0
 neve = inpairsint(envx,fev,dve)
!
! use fe0 to calculate the particle kinetic energy
!
 fev2 = evx*evx*fe0
 fesum2 = inpairsint(envx,fev2,dve)
 kine = half*emas*fesum2
 fiv2 = ivx*ivx*fi0
 fisum2 = inpairsint(invx,fiv2,dvi)
 kini = half*hmas*fisum2
!
! deallocate internal subroutine arrays
! 
 deallocate(felnfe)
 deallocate(filnfi)
 deallocate(fe0)
 deallocate(fi0)
 deallocate(fe0lnfe0)
 deallocate(fi0lnfi0)
 deallocate(fesmall)
 deallocate(fismall)
 deallocate(fev2)
 deallocate(fiv2)
 deallocate(fev)
end subroutine entropy
!
!
!*******************************************************************
!   Clare E. J. Watt			14/02/00
!
 subroutine entropy_e(nx,envx,fe,evx,dx,dve,S,aveS,kine,neve,ne0)
!
!  Subroutine which calculates the entropy of the system using;
!	
!	S = - sum( int f * ln f dx dv ) 
!
!  where the sum is over all plasma species in the system. This 
!  subroutine also calculates the kinetic energy of the electrons
!  and the spatially averaged first moments of fe at each timestep. 
!
!*******************************************************************
!
 use nrtype; use parameters
 implicit none
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of e. velocity grid points
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
 real(DP),dimension(-envx:envx)::evx	! electron velocity grid points
 real(DP),intent(in)::dx		! spatial grid spacing
 real(DP),intent(in)::dve		! electron velocity grid spacing
 real(DP),intent(out)::S		! entropy
 real(DP),intent(out)::aveS		! "entropy" using fe0
 real(DP),intent(out)::kine		! electron kinetic energy
 real(DP),intent(out)::neve		! spatially-ave 1st moment
 real(DP),intent(out)::ne0		! spatially-ave 0th moment
!
 real(DP),dimension(:),allocatable::felnfe	! for calculating S
 real(DP)::fesum				!   "   "    "
 real(DP),dimension(:),allocatable::fe0  	! spatial average of fe
 real(DP),dimension(:),allocatable::fe0lnfe0    ! for calculating aveS
 real(DP),dimension(:,:),allocatable::fesmall
 real(DP),dimension(:),allocatable::fev	        ! v * fe
 real(DP),dimension(:),allocatable::fev2        ! v^2 * fe
 real(DP)::fesum2       			! for calculating kine
 integer(I4B)::i				! counting variable
!
! interfaces
!
 interface
    function inpairsint(nvx,f,dvx)
    use nrtype
    integer(I4B)::nvx    
    real(DP),dimension(-nvx:nvx)::f
    real(DP)::dvx
    real(DP)::inpairsint
    end function inpairsint
 end interface
!
 interface
    function fave(n,nvx,f)
    use nrtype
    integer(I4B)::n,nvx
    real(DP),dimension(0:n,-nvx:nvx)::f
    real(DP),dimension(-nvx:nvx)::fave
    end function fave
 end interface
!
! allocate all internal subroutine arrays
!
 allocate(felnfe(-envx:envx))
 allocate(fe0(-envx:envx))
 allocate(fe0lnfe0(-envx:envx))
 allocate(fesmall(0:nx,-envx:envx))
 allocate(fev2(-envx:envx))
 allocate(fev(-envx:envx))
!
! calculate the entropy of the system
!
 S = zero
 aveS = zero
 kine = zero
 do i = 0,nx
    where (fe(i,-envx:envx) .gt. zero)
       felnfe = fe(i,-envx:envx)*log(fe(i,-envx:envx))
    elsewhere
       felnfe = zero
    endwhere
    fesum = inpairsint(envx,felnfe,dve)
    S = S + fesum
 enddo
 S = -S*dx
!
! calculate the spatially averaged distribution in order to 
! calculate the "entropy" using fe0
!
 fesmall = fe(0:nx,-envx:envx)
 fe0 = fave(nx,envx,fesmall)
 where (fe0 .gt. zero) 
    fe0lnfe0 = fe0*log(fe0)
 elsewhere
    fe0lnfe0 = zero
 endwhere
 aveS = inpairsint(envx,fe0lnfe0,dve)
 aveS = -aveS
!
! use fe0 to calculate the spatially averaged moments
!
 ne0 = inpairsint(envx,fe0,dve)
 fev = evx*fe0
 neve = inpairsint(envx,fev,dve)
!
! use fe0 to calculate the particle kinetic energy
!
 fev2 = evx*evx*fe0
 fesum2 = inpairsint(envx,fev2,dve)
 kine = half*emas*fesum2
!
! deallocate internal subroutine arrays
! 
 deallocate(felnfe)
 deallocate(fe0)
 deallocate(fe0lnfe0)
 deallocate(fesmall)
 deallocate(fev2)
 deallocate(fev)
end subroutine entropy_e
!
!
!*************************************************************************
!  Clare E. J. Watt		          26/01/00
! 
 subroutine wrdata(nx,envx,invx,x,evx,ivx,fe,fi,nex,uex,tex,nix,uix,tix, &
&                  rho,jx,ex,kine,kini,chi,S,effcol,anres,time,stepno, &
&                  ewrite,fwrite,dt)
!
!  Routine which writes out all the relevant data to the correct files.
!  Note that the files are already open, with the unit numbers defined
!  in the parameters module
!
!*************************************************************************
!
 use nrtype; use parameters
!
 implicit none
!
! declare input and output variables
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of e velocity grid points
 integer(I4B),intent(in)::invx		! number of i velocity grid points
 real(DP),dimension(0:nx),intent(in)::x ! x coordinates
 real(DP),dimension(-envx:envx),intent(in)::evx
 real(DP),dimension(-invx:invx),intent(in)::ivx
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
 real(DP),dimension(-1:nx+1,-invx:invx),intent(in)::fi
 real(DP),dimension(0:nx),intent(in)::nex	! electron number density
 real(DP),dimension(0:nx),intent(in)::uex	! electron drift velocity
 real(DP),dimension(0:nx),intent(in)::tex	! electron temperature
 real(DP),dimension(0:nx),intent(in)::nix	! ion number density
 real(DP),dimension(0:nx),intent(in)::uix	! ion drift velocity
 real(DP),dimension(0:nx),intent(in)::tix	! ion temperature
 real(DP),dimension(0:nx),intent(in)::rho	! plasma charge density
 real(DP),dimension(0:nx),intent(in)::jx	! plasma current density
 real(DP),dimension(0:nx),intent(in)::ex	! electric field
 real(DP),intent(in)::kine,kini		        ! kinetic energy values
 real(DP),intent(in)::chi,S			! integral of f, entropy
 real(DP),intent(in)::effcol,anres		! effective collision frequency
 real(DP),intent(in)::time		        ! current time in the sim.
 integer(I4B),intent(in)::stepno		! current timestep
 integer(I4B),intent(in)::ewrite,fwrite	        ! how often E/f are saved
 real(DP),intent(in)::dt
 integer(I4B)::i,j
!
! variables needed for the simulation
!
 real(DP),dimension(:),allocatable::fe0	        ! spatially-averaged e.d.f.
 real(DP),dimension(:),allocatable::fi0	        ! spatially-averaged i.d.f.
 real(DP),dimension(:,:),allocatable::fesmall
 real(DP),dimension(:,:),allocatable::fismall
!
! interface needed for the fave function
!
 interface
    function fave(n,nvx,f)
    use nrtype
    integer(I4B)::n,nvx
    real(DP),dimension(0:n,-nvx:nvx)::f
    real(DP),dimension(-nvx:nvx)::fave
    end function fave
 end interface
!
! backspace outputfile and write the step number to the bottom line
!
 backspace(iwout)
 write(iwout,'(i7)') stepno
!
! allocate the spatially-averaged distribution arrays and the single-
! precision arrays too.
!
 allocate(fesmall(0:nx,-envx:envx)) 
 allocate(fismall(0:nx,-invx:invx))
 allocate(fe0(-envx:envx))
 allocate(fi0(-invx:invx))
!
! get spatially-averaged distributions
!
 fesmall = fe(0:nx,-envx:envx)
 fe0 = fave(nx,envx,fesmall)
 fismall = fi(0:nx,-invx:invx)
 fi0 = fave(nx,invx,fismall)
!
! write the appropriate data to the appropriate data files
!
 write(iw17,'(1p,4(e13.4,2x))') kine,kini,chi,S
 write(iw16,'(e13.4)') anres 
 do i=0,nx
    write(iw7,'(1p,3(e13.4,3x))') x(i), dt*stepno, nex(i)
    write(iw11,'(1p,3(e13.4,3x))') x(i), dt*stepno, uex(i)
    write(iw9,'(1p,3(e13.4,3x))') x(i), dt*stepno, tex(i)
    write(iw8,'(1p,3(e13.4,3x))') x(i), dt*stepno, nix(i)
    write(iw12,'(1p,3(e13.4,3x))') x(i), dt*stepno, uix(i)
    write(iw10,'(1p,3(e13.4,3x))') x(i), dt*stepno, tix(i)
    write(iw13,'(1p,3(e13.4,3x))') x(i), dt*stepno, rho(i)
    write(iw14,'(1p,3(e13.4,3x))') x(i), dt*stepno, jx(i)
    write(iw15,'(1p,3(e13.4,3x))') x(i), dt*stepno, ex(i)
 enddo
 do j=-envx,envx
    write(iw4,'(1p,3(e13.4,3x))') evx(j), dt*stepno, fe0(j)
 enddo
 do j=-invx,invx
    write(iw6,'(1p,3(e13.4,3x))') ivx(j), dt*stepno, fi0(j)
 enddo
!
! Only if it is the correct timestep are the full distribution functions
! written to file
!
 if (mod(stepno,fwrite)==0) then
    do j=-envx,envx
       do i=0,nx
          write (iw3,'(1p,1(e13.4,3x))') fe(i,j)
       enddo
    enddo
    do j=-invx,invx
       do i=0,nx
          write (iw5,'(1p,1(e13.4,3x))') fi(i,j)
       enddo
    enddo
 endif
!
! deallocate arrays
! 
 deallocate(fe0)
 deallocate(fesmall)
 deallocate(fi0)
 deallocate(fismall)
end subroutine wrdata  
!
!
!*************************************************************************
!  Clare E. J. Watt		          26/01/00
! 
 subroutine wrdata_e(nx,envx,x,evx,fe,nex,uex,tex,rho,jx,ex,kine,chi,S, &
&                    effcol,anres,time,stepno,ewrite,fwrite,dt)
!
!  Routine which writes out all the relevant data to the correct files.
!  Note that the files are already open, with the unit numbers defined
!  in the parameters module
!
!*************************************************************************
!
 use nrtype; use parameters
!
 implicit none
!
! declare input and output variables
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of e velocity grid points
 real(DP),dimension(0:nx),intent(in)::x ! x coordinates
 real(DP),dimension(-envx:envx),intent(in)::evx
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
 real(DP),dimension(0:nx),intent(in)::nex	! electron number density
 real(DP),dimension(0:nx),intent(in)::uex	! electron drift velocity
 real(DP),dimension(0:nx),intent(in)::tex	! electron temperature
 real(DP),dimension(0:nx),intent(in)::rho	! plasma charge density
 real(DP),dimension(0:nx),intent(in)::jx	! plasma current density
 real(DP),dimension(0:nx),intent(in)::ex	! electric field
 real(DP),intent(in)::kine		        ! kinetic energy values
 real(DP),intent(in)::chi,S			! integral of f, entropy
 real(DP),intent(in)::effcol,anres		! effective collision frequency
 real(DP),intent(in)::time		        ! current time in the sim.
 integer(I4B),intent(in)::stepno		! current timestep
 integer(I4B),intent(in)::ewrite,fwrite	        ! how often E/f are saved
 real(DP),intent(in)::dt
 integer(I4B)::i,j
!
! variables needed for the simulation
!
 real(DP),dimension(:),allocatable::fe0	        ! spatially-averaged e.d.f.
 real(DP),dimension(:,:),allocatable::fesmall
!
! interface needed for the fave function
!
 interface
    function fave(n,nvx,f)
    use nrtype
    integer(I4B)::n,nvx
    real(DP),dimension(0:n,-nvx:nvx)::f
    real(DP),dimension(-nvx:nvx)::fave
    end function fave
 end interface
!
! backspace outputfile and write the step number to the bottom line
!
 backspace(iwout)
 write(iwout,'(i7)') stepno
!
! allocate the spatially-averaged distribution arrays and the single-
! precision arrays too.
!
 allocate(fesmall(0:nx,-envx:envx)) 
 allocate(fe0(-envx:envx))
!
! get spatially-averaged distributions
!
 fesmall = fe(0:nx,-envx:envx)
 fe0 = fave(nx,envx,fesmall)
!
! write the appropriate data to the appropriate data files
!
 write(iw17,'(1p,3(e13.4,2x))') kine,chi,S
 write(iw16,'(e13.4)') anres 
 do i=0,nx
    write(iw7,'(1p,3(e13.4,3x))') x(i), dt*stepno, nex(i)
    write(iw11,'(1p,3(e13.4,3x))') x(i), dt*stepno, uex(i)
    write(iw9,'(1p,3(e13.4,3x))') x(i), dt*stepno, tex(i)
    write(iw13,'(1p,3(e13.4,3x))') x(i), dt*stepno, rho(i)
    write(iw14,'(1p,3(e13.4,3x))') x(i), dt*stepno, jx(i)
    write(iw15,'(1p,3(e13.4,3x))') x(i), dt*stepno, ex(i)
 enddo
 do j=-envx,envx
    write(iw4,'(1p,3(e13.4,3x))') evx(j), dt*stepno, fe0(j)
 enddo
!
! Only if it is the correct timestep are the full distribution functions
! written to file
!
 if (mod(stepno,fwrite)==0) then
    do j=-envx,envx
       do i=0,nx
          write (iw3,'(1p,1(e13.4,3x))') fe(i,j)
       enddo
    enddo
 endif
!
! deallocate arrays
! 
 deallocate(fe0)
 deallocate(fesmall)
end subroutine wrdata_e
!
!
!
!*************************************************************************
!  Clare E. J. Watt                           19/11/99
!
function fave(n,nvx,f)
!
!  Calculates the spatial average of the distribution function
!
!*************************************************************************
! 
 use nrtype
 implicit none
!
 integer(I4B)::n			! no. of spatial grid points
 integer(I4B)::nvx			! no. of velocity grid points
 real(DP),dimension(0:n,-nvx:nvx)::f	! distribution function
 real(DP),dimension(-nvx:nvx)::fave	! spatially-averaged d.f.
!
 fave = sum(f,dim=1)
 fave = fave/real(n)
end function fave
!
!
!************************************************************************
!   Clare E. J. Watt					06/03/00
!
subroutine f1check(nx,envx,invx,fe,fi)
!
!  This subroutine checks the size of the perturbations at vcut.
!  These can sometimes grow uncontrollably due to vcut being too small. 
!  If this happens, then the simulation stops
!
!************************************************************************
!
 use nrtype; use parameters
 implicit none
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of electron velocity pts
 integer(I4B),intent(in)::invx		! number of ion velocity points
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
 real(DP),dimension(-1:nx+1,-invx:invx),intent(in)::fi
!
 integer(I4B)::i,j				! count variables
 real(DP)::fe0,fi0
 real(DP)::fetest,fitest,fetol,fitol
!
! start with +vcute and +vcuti
!
 fe0 = sum(fe(0:nx,envx))/(nx+1)
 fi0 = sum(fi(0:nx,invx))/(nx+1)
 do i = 0,nx
    fetest = abs(fe(i,envx))-abs(fe0)
    fitest = abs(fi(i,invx))-abs(fi0)
    fetol = maxval(fe(i,:))/1.0e+02_dp
    fitol = maxval(fi(i,:))/1.0e+02_dp
    if (fetest.ge.fetol) then
       write(iwout,*) 'Perturbations at +vcute became too big'
       call zexit
    endif
    if (fitest.ge.fitol) then
       write(iwout,*) 'Perturbations at +vcuti became too big'
       call zexit
    endif
 enddo
!
! now do -vcute and -vcuti
!
 fe0 = sum(fe(0:nx,-envx))/(nx+1)
 fi0 = sum(fi(0:nx,-invx))/(nx+1)
 do i = 0,nx
    fetest = abs(fe(i,-envx))-abs(fe0)
    fitest = abs(fi(i,-invx))-abs(fi0)
    fetol = maxval(fe(i,:))/1.0e+02_dp
    fitol = maxval(fi(i,:))/1.0e+02_dp
    if (fetest.ge.fetol) then
       write(iwout,*) 'Perturbations at -vcute became too big'
       call zexit
    endif
    if (fitest.ge.fitol) then
       write(iwout,*) 'Perturbations at -vcuti became too big'
       call zexit
    endif
 enddo
end subroutine f1check 
!
!************************************************************************
!   Clare E. J. Watt					06/03/00
!
subroutine f1check_e(nx,envx,fe)
!
!  This subroutine checks the size of the perturbations at vcut.
!  These can sometimes grow uncontrollably due to vcut being too small. 
!  If this happens, then the simulation stops
!
!************************************************************************
!
 use nrtype; use parameters
 implicit none
!
 integer(I4B),intent(in)::nx		! number of spatial grid points
 integer(I4B),intent(in)::envx		! number of electron velocity pts
 real(DP),dimension(-1:nx+1,-envx:envx),intent(in)::fe
!
 integer(I4B)::i,j				! count variables
 real(DP)::fe0
 real(DP)::fetest,fetol
!
! start with +vcute
!
 fe0 = sum(fe(0:nx,envx))/(nx+1)
 do i = 0,nx
    fetest = abs(fe(i,envx))-abs(fe0)
    fetol = maxval(fe(i,:))/1.0e+02_dp
    if (fetest.ge.fetol) then
       write(iwout,*) 'Perturbations at +vcute became too big'
       call zexit
    endif
 enddo
!
! now do -vcute
!
 fe0 = sum(fe(0:nx,-envx))/(nx+1)
 do i = 0,nx
    fetest = abs(fe(i,-envx))-abs(fe0)
    fetol = maxval(fe(i,:))/1.0e+02_dp
    if (fetest.ge.fetol) then
       write(iwout,*) 'Perturbations at -vcute became too big'
       call zexit
    endif
 enddo
end subroutine f1check_e 
!          
!********************************************************
!  version 1.00                     15/11/99
!
subroutine zexit
!
!  Routine to halt program nicely if an error occurs
!
!********************************************************
!
 print*, ' Program has stopped prematurely'
 stop
end subroutine zexit           
