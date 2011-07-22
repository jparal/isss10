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
 real(DP),parameter::hmas=1.67262158e-27_dp
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
!  Clare E. J. Watt                                29/02/00 
!
program iaisetup
!
!  Double - VLASOV simulation code: file setup program
!  Original electron code written by Richard B Horne in f77
!  Modified by Mervyn P Freeman
!  Modified by Clare Watt to include ions and to calculate "best
!  grid" for input parameters specified.
!
!  The first of two main driver programs which run a Vlasov simulation
!  in one dimension. The original Vlasov code has been modified to 
!  include ions. The first program sets up the simulation grid
!  and the second performs the simulation from the initial data 
!  given by this program.
!
!  This code is specified to look at unstable ion-acoustic waves, but 
!  minor changes can make it more general.  
!
!***********************************************************************
!
!
 use nrtype; use parameters
 implicit none
!
 integer(I4B)::ilen			! length of root
 character(len = 8)::filename	 	! filename of input file
 character(len = 8)::fil1
 character(len = 9)::fil2
 character(len = 13)::fil3,fil4,fil5,fil6,fil7,fil8,fil9,fil10
 character(len = 13)::fil11,fil12,fil13,fil14,fil15,fil16,fil17,fil18
 character(len = 16)::fil20
 character(len = 17)::fil21,fil22
 character(len = 12)::fil23
!
! Input parameters
!
 integer(I4B)::ncomps                   ! number of components
 character(len=10),dimension(:),allocatable::form ! form of the distributions 
 integer(I4B),dimension(:),allocatable::species  ! species code for components
 real(DP),dimension(:),allocatable::mass! mass of each component
 real(DP),dimension(:),allocatable::nden! number density of components
 real(DP),dimension(:),allocatable::temp! temperature of the components
 real(DP),dimension(:),allocatable::drift  ! drift velocity of the components
 real(DP)::mvth				! multiple of thermal speed for vcut
 real(DP)::cf				! courant factor
 integer(I4B)::ntmax			! total number of timesteps
 integer(I4B)::ewrite,fwrite		! how often E/f are saved
 real(DP)::enoise			! amplitude of internal E field
!
! variables needed for grid calculations
!
 real(DP)::enden,inden                  ! total electron/proton density
 real(DP)::wpe,wpi                      ! plasma frequency
 real(DP)::wpe2,wpi2                    ! plasma frequency squared
 real(DP)::debye2			! debye length, squared
 real(DP)::debye			! debye length
 real(DP)::cs,cs2			! ion acoustic speed (squared)
 real(DP)::etemp                        ! largest electron temperature
 real(DP)::itemp                        ! largest ion temperature
 real(DP)::edrift                       ! largest electron drift velocity
 real(DP)::idrift                       ! largest ion drift velocity
 real(DP)::ae,ai                        ! largest electron and ion thermal sp
 real(DP),dimension(:),allocatable::nu  ! fraction of component nden
 real(DP),dimension(:),allocatable::alpha ! thermal speeds
 real(DP)::kmin,kmax			! min, max values of k
 real(DP)::lmax,lmin			! max,min wavelengths
 real(DP)::dve1,dve2,dvi1,dvi2		! tests for velocity grid spacing
 real(DP)::dt1,dt2			! tests for timesteps	
 real(DP)::dx				! spatial grid spacing
 real(DP)::dve				! electron velocity grid spacing
 real(DP)::dvi				! ion velocity grid spacing
 real(DP)::dt				! timestep
 integer(I4B)::nx			! number of spatial grid points
 integer(I4B)::envx			! no. of electron vel. grid points
 integer(I4B)::invx			! no. of ion vel. grid points
 integer(I4B)::dimrje,dimnut		! no. of data points in rje/i files
 integer(I4B)::dimfe0			! no. of data points in fe0 file
 integer(I4B)::dimfi0			! no. of data points in fi0 file
 integer(I4B)::dimres			! no. of points in resistivity file
 real(DP)::Lx				! length of the simulation box
 integer(I4B),dimension(1)::gloc	! location of max value of gamma
 real(DP)::maxgamma			! maximum value of gamma 
 real(DP)::ttest			! timestep*maxgamma
 real(DP)::vcute,vcuti			! electron/ion cutoff velocity 
 real(DP)::vphmin,vphmax		! min, max phase velocity
 real(DP)::vtest			! phase vel. test
 integer(I4B)::kminloc,kmaxloc		! locations of min/max k
 real(DP)::indication			! indication of calculations 
 real(DP),dimension(:),allocatable::icoord    ! spatial coordinate variable
 real(DP),dimension(:),allocatable::ejcoord   ! electron velocity coord var
 real(DP),dimension(:),allocatable::ijcoord   ! ion velocity coord var
 integer(I4B),parameter::nmax = 200000	! no. of k-values studied
 integer(I4B)::nkmodes		        ! no. of unstable modes
 complex(DPC)::guess			! guess for newton algorithm	
 complex(DPC),dimension(nmax)::k	! all the k-values studied 
 complex(DPC),dimension(nmax)::omega 	! all the omega values 
 real(DP),dimension(nmax)::wr		! the real(DP) part of omega
 real(DP),dimension(nmax)::gamma	! the imaginary part of omega
 real(DP),dimension(nmax)::kr		! the real(DP) part of k
 real(DP)::emax				! expected highest value of efield
 logical::finished
 integer(I4B)::nk
!
! count variables for do-loops
!
 integer(I4B)::i,j,y
!
!*****************************************************************
!
! first, set up the filenames
!
 write (iwl,*) ' Please enter input filename: '
 read(irl,'(a)') filename
 ilen = index(filename,'.in')
 if (ilen ==0) then
    write(iwl,*) 'Stopping, data file must have an extension of type .in'
    call zexit
 endif
 fil1(1:ilen)=filename(1:ilen)
 fil1(ilen:ilen+2)='.in'
 fil2(1:ilen)=filename(1:ilen)
 fil2(ilen:ilen+3)='.out'   	      
 fil20(1:ilen)=filename(1:ilen)
 fil20(ilen:ilen+10)='.xcoord.dat'
 fil21(1:ilen)=filename(1:ilen)
 fil21(ilen:ilen+11)='.evcoord.dat'
 fil22(1:ilen)=filename(1:ilen)
 fil22(ilen:ilen+11)='.ivcoord.dat'
 fil23(1:ilen) = filename(1:ilen)
 fil23(ilen:ilen+6)='.dr.dat'
!
! open files ready for reading/writing
!
 open (unit=irin,file=fil1,form='formatted',action='read')
 open (unit=iwout,file=fil2,form='formatted',action='write')
 open (unit=iw20,file=fil20,form='formatted',action='write')
 open (unit=iw21,file=fil21,form='formatted',action='write')
 open (unit=iw22,file=fil22,form='formatted',action='write')
!
! read input parameters from the input file
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
 allocate(drift(ncomps))
 allocate(alpha(ncomps))
 allocate(nu(ncomps))
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
 read (irin, *) ewrite
 read (irin, *) fwrite
!
! write input data to output file
!
 write (iwout, '(1x,a24)') 'THE INPUT PARAMETERS ARE:'
 do y=1,ncomps
    write (iwout, '(/,a10,1x,a22,i3)') form(y), '     form of component', y
    write (iwout, '(a60)') 'Number density, temperature and drift  &
&                           velocity of component:'
    write (iwout, '(1p,e11.4,3x)') nden(y)
    write (iwout, '(1p,e11.4,3x)') temp(y)
    write (iwout, '(1p,e11.4,3x)') drift(y)
 enddo
 write (iwout, 50) mvth,cf
50 format(/,0p, &
&  f6.3,5x,'     mvth: multiple of thermal speed for vcut',/, &
&  f6.3,5x,'     cf: courant factor')
!
! calculate thermal velocities, plasma frequencies, ion acoustic speed 
! and debye length
!
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
 wpe2 = enden*ev*ev/(e0*emas)
 wpe = sqrt(wpe2)
 wpi2 = enden*ev*ev/(e0*hmas)
 wpi = sqrt(wpi2)
 etemp = zero
 itemp = zero
 ae = zero
 ai = zero
 do y=1,ncomps
    if ((species(y)==1).and.(temp(y).gt.etemp)) then
       etemp = temp(y)
       ae = alpha(y)
    endif
    if ((species(y)==2).and.(temp(y).gt.itemp)) then
       itemp = temp(y)
       ai = alpha(y)
    endif
 enddo
 debye2 = e0*etemp/(enden*ev)
 debye = sqrt(debye2)
 cs2 = ev*(etemp + 3.0_dp*itemp)/hmas
 cs = sqrt(cs2)  
!
! the level of enoise is calculated from the thermal fluctuation level of the 
! plasma: W = kBTe/debye^3 (see Treumann and Baumjohann, Advanced
! Space Plasma Physics, p25)
!
 enoise = sqrt(two*ev*etemp/(e0*debye**3))
!
! write this information to the output file
!
 write (iwout, 55) enoise
55 format(1p,e11.4,'     enoise: Initial e-field perturbation amplitude',/)
 write (iwout, *) 'PLASMA DIAGNOSTICS'
 write (iwout, 60) ae,ai,wpe,wpi,debye,cs
60 format(&
&      1p,e11.4,'     ae: Electron thermal speed (main bulk component)',/, &
&         e11.4,'     ai: Ion thermal speed (main bulk component)',/, &
&         e11.4,'     wpe: Electron plasma frequency',/, & 
&         e11.4,'     wpi: Ion plasma frequency',/, &    
&      0p,f11.4,'     debye: Debye length (using electron temperature)',/, &
&      1p,e11.4,'     cs: Ion acoustic speed',/ )
!
! This code can be run to investigate microinstabilities which involve
! both protons and electrons, and also it can be run in *quick mode* to
! investigate electron-only microinstabilities (e.g. electron bump-on-tail
! instability). If only electron components have been specified, then the
! next bit of code is a little different, so here we have to make a decision
!
 if (inden.gt.zero) then
!
! Look for unstable values of k in an electron/proton plasma.  
! Dispersion relation solver is a newton algorithm which relies on 
! the distribution functions being Maxwellian. For more general solutions, 
! a more general dr solver should be used in place of newton routine below. 
! Note that the *guess* parameter determines which roots of the DR will
! be obtained. By altering this parameter, other important roots can be 
! identified. The type of microinstability should determine a good value for
! the initial guess (see e.g. Theory of Space Plasma Microinstabilities, by
! S. P. Gary for some good approximations).
!
    open (unit = iw23,file=fil23,action='write',form='formatted')
    do j = 1,nmax
       k(j) = cmplx(j*1.0e-06_dp,zero)
       if (j==1) then
          guess = sqrt(((k(j)**2)*cs2)/(1 + (k(j)**2)*debye2))
       else
          guess = omega(j-1)
       endif
       call newton(ncomps,nden,alpha,drift,wpe2,wpi2,nu,species, &
&                  guess,k(j),omega(j))
       write (iw23,'(1p,3(e11.4,2x))') real(k(j)), real(omega(j)), &
&                                     imag(omega(j))
    enddo
    close (iw23)
    wr = real(omega)
    gamma = imag(omega)
    kr = real(k)
    gloc = maxloc(gamma)
    maxgamma = gamma(gloc(1))
    finished = .false.
    kminloc = 2
    vphmin = wr(kminloc)/kr(kminloc)
    vphmax = wr(kminloc)/kr(kminloc)
    do y = 2,nmax
       if ((gamma(y).gt.(5.0e-02_dp*maxgamma)).and.(finished.eqv..false.)) then
          kmin = kr(y)
          kminloc = y
          finished = .true.
       else if ((gamma(y).lt.(5.0e-02_dp*maxgamma)).and.(y.gt.kminloc).and. &
&               (finished.eqv..true.)) then
          kmax = kr(y)
          kmaxloc = y
          exit
       endif  
       vtest = wr(y)/kr(y)
       if ((vtest > vphmax) .and. (gamma(y)>zero)) then
          vphmax = vtest
       endif
       if ((vtest < vphmin) .and. (gamma(y)>zero)) then
          vphmin = vtest
       endif
    enddo
!
! if kmaxloc has not been assigned by this point, then not enough
! k-values have been looked at. Simulation is stopped.
!
    if (kmaxloc == 0) then
       print*, 'still growing solutions at nmax'
       print*, gamma(nmax)
       call zexit
    endif
!
! calculate spatial grid size - dx is a tenth of the shortest 
! growing wavelength, or the Debye length, whichever is shorter. 
! The length of the grid is the longest growing wavelength.
!
    lmin = two*pi/kmax
    dx = lmin/10.0_dp
    if (dx .gt. debye) then
       dx = debye
    endif
    lmax = two*pi/kmin
    Lx = lmax
    nx = ceiling(Lx/dx)
    kmin = kmin/two
!
! calculate velocity grid size. vcute should be large enough to take into 
! account the largest drift velocity, and vcuti should be large enough to 
! cover all the resonant phase velocities (from linear dr for now). 
! dv should be smaller than the smallest resonant phase velocity, and also 
! allow for adequate coverage of the resonant region.   
!
    ae = zero
    ai = zero
    edrift = zero
    idrift = zero
    do y=1,ncomps
       if ((species(y)==1).and.(drift(y).gt.edrift)) then
          edrift = drift(y)
       endif
       if ((species(y)==2).and.(drift(y).gt.idrift)) then
          idrift = drift(y)
       endif
       if ((species(y)==1).and.(alpha(y).gt.ae)) then
          ae = alpha(y)
       endif
       if ((species(y)==2).and.(alpha(y).gt.ai)) then
          ai = alpha(y)
       endif
    enddo
    vcute = mvth*ae + edrift
    if (vcute .gt. 3.0e+08_dp) then
       vcute = 2.95e+08_dp
    endif
    dve1 = vphmin/5.0_dp
    dve2 = wpi/kmax
    if (dve1.lt.dve2) then
       dve = dve1
    else
       dve = 0.9_dp*dve2
    endif
    envx = ceiling(vcute/dve)
    vcute = envx*dve
    dvi1 = 1.0e-01_dp*ai
    dvi2 = vphmin/5.0_dp
    if (dvi1.lt.dvi2) then
       dvi = dvi1
    else
       dvi = dvi2
    endif
    if ((10.0_dp*ai).lt.(10.0_dp*vphmax)) then
       vcuti = 10.0_dp*vphmax
    else
       vcuti = 10.0_dp*ai
    endif
    invx = ceiling(vcuti/dvi)
!
! calculate the size of the timestep, remembering it must satisfy all 
! Courant-type stability conditions, and it must  be smaller than the 
! period of the highest frequency waves in the simulation. 
! 
    emax = 0.2_dp
    if (emas*dve .lt. hmas*dvi) then
       dt1 = cf*emas*dve/(ev*emax)
    else
       dt1 = cf*hmas*dvi/(ev*emax)
    endif
    dt2 = cf*dx/vcute
    if (dt2.gt.dt1) then
       dt = dt1
    else
       dt = dt2
    endif
!
! The timestep should also reflect the size of the growth rates. 
! If the growth rates are high, dt should be shortened accordingly.
! The value of 1.0d-04 is empirical.
!
    ttest = maxgamma*dt
    if (ttest .gt. 1.0e-04_dp) then
       dt = 1.0e-04_dp/maxgamma 
    endif
!
! Now write out all the information gained so far to the output file
!
    write (iwout, *) 'INFORMATION ABOUT UNSTABLE WAVES'
    write (iwout, 70) kmin,kmax,lmin,lmax,vphmin,vphmax,maxgamma
70 format(&
&      1p,e11.4,'     kmin: The smallest value of unstable k',/, &
&         e11.4,'     kmax: The largest value of unstable k',/, &
&         e11.4,'     lmin: The shortest wavelength',/, &
&         e11.4,'     lmax: The longest wavelength',/, &
& 'The length of the box is such that it can hold the longest' ,/, &
& 'wavelength',/, &
&         e11.4,'     vphmin: The phase velocity at kmin',/, &
&         e11.4,'     vphmax: The phase velocity at kmax',/, &
&         e11.4,'     maxgamma: The maximum value of gamma',/)
    write (iwout, *) 'GRID DIMENSIONS'
    write (iwout, 80) Lx,vcute,vcuti,dx,dve,dvi,dt,nx,envx,invx,ntmax, &
&                     ntmax*dt*wpe
80 format(&
&      1p,e11.4,'     Lx: box size in x, in m',/, &
&         e11.4,'     vcute: box size in vx for electrons, in m/s',/, &
&         e11.4,'     vcuti: box size in vx for ions, in m/s',/, &
&      0p,f11.4,'     dx: stepsize in x, in m',/, &
&      1p,e11.4,'     dve: stepsize in vx for electrons, in m/s',/, &
&         e11.4,'     dvi: stepsize in vx for ions, in m/s',/, &
&         e11.4,'     dt: timestep, in s',/, &
&      0p,i6,5x,'     nx: number of spatial grid points',/, &
&         i6,5x,'     envx: number of ve points',/, &
&         i6,5x,'     invx: number of vi points',/, &
&         i7,4x,'     ntmax: number of timesteps',/, &
&      1p,e11.4,'     phystime: total simulation time in 1/wpe',/)
!
! Now write coordinate files 
!
    do i = 0,nx
       write (iw20,'(1p,e13.4)') i*dx
    enddo
    do j = -envx,envx
       write (iw21,'(1p,e13.4)') j*dve
    enddo
    do j = -invx,invx
       write (iw22,'(1p,e13.4)') j*dvi
    enddo  
!
 else
!
! If there are only electron components, then the guess for the newton routine
! is a little different. 
! Note that the *guess* parameter determines which roots of the DR will
! be obtained. By altering this parameter, other important roots can be 
! identified. The type of microinstability should determine a good value for
! the initial guess (see e.g. Theory of Space Plasma Microinstabilities, by
! S. P. Gary for some good approximations).
!
    ae = zero
    ai = zero
    edrift = zero
    idrift = zero
    do y=1,ncomps
       if ((species(y)==1).and.(drift(y).gt.edrift)) then
          edrift = drift(y)
       endif
       if ((species(y)==2).and.(drift(y).gt.idrift)) then
          idrift = drift(y)
       endif
       if ((species(y)==1).and.(alpha(y).gt.ae)) then
          ae = alpha(y)
       endif
       if ((species(y)==2).and.(alpha(y).gt.ai)) then
          ai = alpha(y)
       endif
    enddo
!
    open (unit = iw23,file=fil23,action='write',form='formatted')
    finished = .false.
    nk = 1
    do while (.not.finished)
       k(nk) = cmplx(nk*1.0e-06_dp,zero)
       if (nk==1) then
          guess = k(nk)*edrift
       else
          guess = half*(k(nk)*edrift+omega(nk-1))
       endif
       call newton(ncomps,nden,alpha,drift,wpe2,wpi2,nu,species, &
&                  guess,k(nk),omega(nk))
       write (iw23,'(1p,3(e11.4,2x))') real(k(nk)), real(omega(nk)), &
&                                     imag(omega(nk))
       if (imag(omega(nk)).lt.zero) then
          finished = .true.
       endif
       nk = nk+1
       if (nk.gt.nmax) then
          finished = .true.
       endif
    enddo
    close (iw23)
    wr = real(omega)
    gamma = imag(omega)
    kr = real(k)
    gloc = maxloc(gamma)
    maxgamma = gamma(gloc(1))
    finished = .false.
    kminloc = 2
    vphmin = wr(kminloc)/kr(kminloc)
    vphmax = wr(kminloc)/kr(kminloc)
    do y = 2,nmax
       if ((gamma(y).gt.(1.0e-01_dp*maxgamma)).and.(finished.eqv..false.)) then
          kmin = kr(y)
          kminloc = y
          finished = .true.
       else if ((gamma(y).lt.(1.0e-01_dp*maxgamma)).and.(y.gt.kminloc).and. &
&               (finished.eqv..true.)) then
          kmax = kr(y)
          kmaxloc = y
          exit
       endif  
       vtest = wr(y)/kr(y)
       if ((vtest > vphmax) .and. (gamma(y)>zero)) then
          vphmax = vtest
       endif
       if ((vtest < vphmin) .and. (gamma(y)>zero)) then
          vphmin = vtest
       endif
    enddo
!
! if kmaxloc has not been assigned by this point, then not enough
! k-values have been looked at. Simulation is stopped.
!
    if (kmaxloc == 0) then
       print*, 'still growing solutions at nmax'
       print*, gamma(nmax)
       call zexit
    endif
!
! calculate spatial grid size  
! The length of the grid is the longest growing wavelength.
!
    lmin = two*pi/kmax
    dx = lmin/4.0_dp
    if (dx.gt.debye) then
       dx = debye
    endif
    lmax = two*pi/kmin
    Lx = lmax
    nx = ceiling(Lx/dx)
    kmin = kmin/two
!
! calculate velocity grid size. vcute should be large enough to take into 
! account the largest drift velocity, and vcuti should be large enough to 
! cover all the resonant phase velocities (from linear dr for now). 
! dv should be smaller than the smallest resonant phase velocity, and also 
! allow for adequate coverage of the resonant region.   
!
    vcute = mvth*ae + edrift
    if (vcute .gt. 3.0e+08_dp) then
       vcute = 2.95e+08_dp
    endif
    dve = (vphmax-vphmin)/10.0_dp
    envx = ceiling(vcute/dve)
    vcute = envx*dve
!
! calculate the size of the timestep, remembering it must satisfy all 
! Courant-type stability conditions, and it must  be smaller than the 
! period of the highest frequency waves in the simulation. 
! 
    emax = 0.2_dp
    dt1 = cf*emas*dve/(ev*emax)
    dt2 = cf*dx/vcute
    if (dt2.gt.dt1) then
       dt = dt1
    else
       dt = dt2
    endif
!
! The timestep should also reflect the size of the growth rates. 
! If the growth rates are high, dt should be shortened accordingly.
! The value of 1.0d-04 is empirical.
!
    ttest = maxgamma*dt
    if (ttest .gt. 1.0e-04_dp) then
       dt = 1.0e-04_dp/maxgamma 
    endif
!
! Now write out all the information gained so far to the output file
!
    write (iwout, *) 'INFORMATION ABOUT UNSTABLE WAVES'
    write (iwout, 71) kmin,kmax,lmin,lmax,vphmin,vphmax,maxgamma
71 format(&
&      1p,e11.4,'     kmin: The smallest value of unstable k',/, &
&         e11.4,'     kmax: The largest value of unstable k',/, &
&         e11.4,'     lmin: The shortest wavelength',/, &
&         e11.4,'     lmax: The longest wavelength',/, &
& 'The length of the box is such that it can hold the longest' ,/, &
& 'wavelength',/, &
&         e11.4,'     vphmin: The phase velocity at kmin',/, &
&         e11.4,'     vphmax: The phase velocity at kmax',/, &
&         e11.4,'     maxgamma: The maximum value of gamma',/)
    write (iwout, *) 'GRID DIMENSIONS'
    write (iwout, 81) Lx,vcute,dx,dve,dt,nx,envx,ntmax, &
&                     ntmax*dt*wpe
81 format(&
&      1p,e11.4,'     Lx: box size in x, in m',/, &
&         e11.4,'     vcute: box size in vx for electrons, in m/s',/, &
&      0p,f11.4,'     dx: stepsize in x, in m',/, &
&      1p,e11.4,'     dve: stepsize in vx for electrons, in m/s',/, &
&         e11.4,'     dt: timestep, in s',/, &
&      0p,i6,5x,'     nx: number of spatial grid points',/, &
&         i6,5x,'     envx: number of ve points',/, &
&         i7,4x,'     ntmax: number of timesteps',/, &
&      1p,e11.4,'     phystime: total simulation time in 1/wpe',/)
!
! Now write coordinate files 
!
    do i = 0,nx
       write (iw20,'(1p,e13.4)') i*dx
    enddo
    do j = -envx,envx
       write (iw21,'(1p,e13.4)') j*dve
    enddo
 endif
 close(irin)
 close(iwout)
 close (iw20); close (iw21); close (iw22)
!
! Deallocate all relevant variables
!
 deallocate(form)
 deallocate(species)
 deallocate(mass)
 deallocate(nden)
 deallocate(temp)
 deallocate(drift)
 deallocate(alpha)
 deallocate(nu)
!
end program iaisetup
!
! 
!
!**************************************************************************
! Clare E. J. Watt                                          16/01/00
!
subroutine newton(ncomps,nden,alpha,vd,wpe2,wpi2,nu,species,wguess,k,omega)
!
! This subroutine solves the dispersion relation using a Newton-Raphson 
! numerical algorithm. The dispersion relation is an electrostatic one-
! dimensional equation for Maxwellian components. 
!
!****************************************************************************
!
 use nrtype; use parameters
 implicit none 
!
 integer(I4B),intent(in)::ncomps                 ! number of components
 real(DP),dimension(ncomps),intent(in)::nden     ! number density
 real(DP),dimension(ncomps),intent(in)::alpha    ! thermal speeds
 real(DP),dimension(ncomps),intent(in)::vd       ! drift velocities
 real(DP),intent(in)::wpe2		         ! sq electron plasma frequency
 real(DP),intent(in)::wpi2		         ! sq ion plasma frequency
 real(DP),dimension(ncomps),intent(in)::nu       ! density ratio
 integer(I4B),dimension(ncomps),intent(in)::species  ! 1=electrons,2=protons
 complex(DPC),intent(in)::wguess	         ! input guess at solution
 complex(DPC),intent(in)::k		         ! wavenumber
 complex(DPC),intent(out)::omega	         ! output solution 
 complex(DPC)::zeta			         ! argument for electron p.d.f. 
 complex(DPC)::z			         ! electron p.d.f.
 complex(DPC)::zp			         ! derivative of electron p.d.f.
 complex(DPC)::dr			         ! dispersion relation
 complex(DPC)::dr_prime		                 ! derivative of DR
 complex(DPC)::w1,w2		                 ! used in Newton's method
 logical::finished		                 ! stops loop
 real(DP)::test1,test2		                 ! tests to check convergence 
 integer(I4B)::i,j,counter,y
 complex(DPC),parameter::one_c = (one,zero),two_c = (two,zero)
! 
 counter = 0   
 finished = .false.
 do while (.not. finished)
    if (counter ==0) then
       w1 = wguess
       dr = one_c
       dr_prime = zero
       do y=1,ncomps
          zeta = (w1-vd(y)*k)/(alpha(y)*k)
          call fried(zeta,z,zp)
          if (species(y)==1) then
             dr = dr - (nu(y)*wpe2*zp)/((k**2)*(alpha(y)**2))
             dr_prime = dr_prime + two_c*nu(y)*wpe2*(z + zeta*zp)/ &
&                                  ((k**3)*(alpha(y)**3))
          else
             dr = dr - (nu(y)*wpi2*zp)/((k**2)*(alpha(y)**2))
             dr_prime = dr_prime + two_c*nu(y)*wpi2*(z + zeta*zp)/ &
&                                  ((k**3)*(alpha(y)**3))
          endif
       enddo
    else
       w1 = w2
    endif
    w2 = w1 - (dr/dr_prime)
    counter = counter + 1
    dr = one_c
    dr_prime = zero
    do y=1,ncomps
       zeta = (w2-vd(y)*k)/(alpha(y)*k)
       call fried(zeta,z,zp)
       if (species(y)==1) then
          dr = dr - (nu(y)*wpe2*zp)/((k**2)*(alpha(y)**2))
          dr_prime = dr_prime + two_c*nu(y)*wpe2*(z + zeta*zp)/ &
&                               ((k**3)*(alpha(y)**3))
       else
          dr = dr - (nu(y)*wpi2*zp)/((k**2)*(alpha(y)**2))
          dr_prime = dr_prime + two_c*nu(y)*wpi2*(z + zeta*zp)/ &
&                               ((k**3)*(alpha(y)**3))
       endif
    enddo
    test1 = abs(w2-w1)
    test2 = abs(dr)
    if ((test1 < 1.0e-06_dp) .or. (test2 < 1.0e-10_dp)) then
       finished = .true.
       omega = w2
    else if (counter > 20) then
       finished = .true.
       print*, 'No root found'
       print*, k
       omega = (zero,zero)
       call zexit
    else if (dr_prime == cmplx(zero,zero)) then
       finished = .true.
       omega = (zero,zero)
       print*, 'Turning-point found'
    end if
 end do
end subroutine newton
!
!**************************************************************************
      SUBROUTINE fried(zeta,z,zp)
!**************************************************************************
!
! This subroutine calculates the plasma dispersion function and the 
! derivative of the plasma dispersion function wrt zeta. The routine is 
! borrowed from Richard Horne's hotray code and modified only slightly to fit 
! in with this piece of code.
! 
      use nrtype; use parameters
      IMPLICIT NONE
      INTEGER(I4B),PARAMETER::kc = 10
      COMPLEX(DPC),INTENT(IN)::zeta
      COMPLEX(DPC),INTENT(OUT)::z,zp 
!
      REAL(DP)::x,y,x1,p_r,p_i,xyz,a,y1,t,ar,ai,ppr,ppi
      REAL(DP),PARAMETER::rpi=1.77245385090552_dp
      REAL(DP),PARAMETER::valmax=80.0_dp
!
      x=dreal(zeta)
      y=dimag(zeta)
      x1=abs(x)
      if(y) 11,10,10
 10   call wz1(x1,y,p_r,p_i)
      if(x) 12,13,13
  11  xyz=y*y-x*x
      if(xyz.ge.zero)then
        a=two*exp(xyz)
      else
        a=zero
        if(abs(xyz).lt.valmax)then
          a=two*exp(xyz)
        end if
      end if
      y1=-y
      t=x1*y1+x1*y1
      ar=a*cos(t)
      ai=a*sin(t)
      call wz1(x1,y1,p_r,p_i)
      p_r=-p_r+ar
      p_i= p_i+ai
      if(x) 12,13,13
 12   p_r=p_r
      p_i=-p_i
 13   continue
      a=p_i
      p_i=rpi*p_r
      p_r=-rpi*a
      ppr=-two*(one+x*p_r-y*p_i)
      ppi=-two*(y*p_r+x*p_i)
      z=dcmplx(p_r,p_i)
      zp=dcmplx(ppr,ppi)
      return
      END SUBROUTINE fried
!
!
!**************************************************************************
      SUBROUTINE wz1(x,y,preel,pimag)
!**************************************************************************
!
! This subroutine is called in the above fried routine. It is also 
! borrowed from Richard Horne's hotray code.
! 
      use nrtype; use parameters
      IMPLICIT NONE
      REAL(DP),INTENT(IN)::x,y
      REAL(DP),INTENT(OUT)::preel,pimag 
      INTEGER(I4B),PARAMETER::kc = 10
      INTEGER(I4B)::ierr,icapn,nu,ib,nup1,n,np1
      REAL(DP)::s,h,h2,alamda,r1,r2,s1,s2,t1,t2
      REAL(DP)::c,rich1,rich2,x2
      REAL(DP),PARAMETER::coef=0.112837916709551e01_dp
      REAL(DP),PARAMETER::valmax=80.0_dp
!
      ierr=0
      h2=zero
      alamda=zero
!
!      if((y.ge.0.429d01).or.(x.ge.0.533d01))go to 1
      if((y.lt.0.429e01_dp).and.(x.lt.0.533e01_dp))then
         s=(0.1e01_dp-y/0.429e01_dp)*sqrt(0.1e01_dp-x*x/0.2841e02_dp)
         h=0.16e01_dp*s
         h2=h+h
         icapn=int(0.65e01_dp+0.23e0_dp*s)
         alamda=h2**icapn
         nu=int(0.95e01_dp+0.21e02_dp*s)
!      go to 2
      else
!       continue
!
!  Note that h2 and alamda are not defined here so set a flag
!  so that if they are used due to rounding errors below a 
!  message is printed.
!  Richard horne 16 Oct 91
!
         ierr=1
         h=zero
         icapn=0
         nu=8
      end if
!    continue
      ib=0
      if((h.lt.0.1e-11_dp).or.(alamda.lt.0.1e-11_dp))ib=1
      r1=zero
      r2=zero
      s1=zero
      s2=zero
      nup1=nu+1
      do 3 n=1,nup1
      np1=nup1-n+1
      t1=y+h+dble(np1)*r1
      t2=x-dble(np1)*r2
      c=half/(t1*t1+t2*t2)
      r1=c*t1
      r2=c*t2
!     if((h.le.zero).or.((np1-1).gt.icapn))go to 3
      if((h.gt.zero).and.((np1-1).le.icapn))then
         t1=alamda+s1
         s1=r1*t1-r2*s2
         s2=r2*t1+r1*s2
         alamda=alamda/h2
         if(ierr.eq.1)then
            PRINT*,' wz1:0: rounding error stopping'
            PRINT*,' wz1:0: rounding error stopping'
            call zexit
         end if
      end if
 3    continue
      rich1=dble(ib)
      rich2=dble(1-ib)
      pimag=coef*(rich1*r2+rich2*s2)
      if(y.eq.zero) go to 5
      preel=coef*(dble(ib)*r1+dble(1-ib)*s1)
      go to 999
  5   continue
      preel=zero
      x2=x*x
      if(x2.lt.valmax)then
        preel=exp(-x2)
      end if
 999  continue
      return
      END SUBROUTINE wz1
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
