! $Id: interstellar.f90,v 1.80 2004-03-14 16:07:25 mee Exp $

!  This modules contains the routines for SNe-driven ISM simulations.
!  Still in development. 

module Interstellar

  use Cparam
  use Cdata
  use Density

  implicit none

  real :: x_SN,y_SN,z_SN,rho_SN,lnrho_SN,yH_SN,lnTT_SN,TT_SN,ss_SN,ee_SN
  integer :: l_SN,m_SN,n_SN
  real, dimension(nx) :: dr2_SN     ! Pencil storing radius to SN

  ! Save space for last SNI time
  integer, parameter :: ninterstellarsave=1
  real, dimension(ninterstellarsave) :: interstellarsave
  real :: t_next_SNI=0.0
  real :: t_interval_SNI=impossible

  
  ! Mesh width (in points) of a SNe insertion
  integer :: point_width=4

  ! normalisation factors for 1-d, 2-d, and 3-d profiles like exp(-r^6)
  ! ( 1d: 2    int_0^infty exp(-(r/a)^6)     dr) / a
  !   2d: 2 pi int_0^infty exp(-(r/a)^6) r   dr) / a^2
  !   3d: 4 pi int_0^infty exp(-(r/a)^6) r^2 dr) / a^3 )
  ! ( cf. 3.128289613 -- from where ?!? )
  ! NB: 1d and 2d results just from numerical integration -- calculate
  !      exact integrals at some point...
  real, parameter, dimension(3) :: &
                        cnorm_SN = (/ 1.85544 , 2.80538 , 3.71213666 /) 

!  cp1=1/cp used to convert TT (and ss) into interstellar code units
!  (useful, as many conditions conveniently expressed in terms of TT)
!  code units based on:
!    [length]  = 1kpc  = 3.09 10^21 cm
!    [time]    = 1Gyr  = 3.15 10^16 s             !no, on [u]=1km/s...
!    [rho]     =       = 1.00 10^-24 g/cm^3
!  Lambdaunits converts coolH into interstellar code units.
!   (this should really just be incorporated into coolH coefficients)
!  NB: will start using thermodynamics, and unit_length, etc., imminently...

!  real, parameter :: cp1=27.8   !=R * gamma / (mu * (gamma-1))  27.8 
!  real, parameter :: TTunits=46.6

  double precision :: unit_Lambda

  ! Minimum resulting central temperature of a SN explosion. Move mass to acheive this.
  real :: TT_SN_min_cgs=1.e7
  double precision, parameter :: SNI_area_rate_cgs=1.330982784D-52
  double precision, parameter :: solar_mass_cgs=1.989e33
  real, parameter :: h_SNI_cgs=1.00295e19,h_SNII_cgs=2.7774e18
  real, parameter :: rho_crit_cgs=1.e-24,TT_crit_cgs=4000.
  double precision, parameter :: ampl_SN_cgs=10D51

  ! Minimum resulting central temperature of a SN explosion. Move mass to acheive this.
  real :: TT_SN_min=impossible
  real :: SNI_area_rate=impossible
  real :: h_SNI=impossible,h_SNII=impossible
  real :: solar_mass=impossible, ampl_SN=impossible
  real :: rho_crit=impossible,TT_crit=impossible


  real, parameter :: rhoUV_cgs=0.1
  real, parameter :: TUV_cgs=7000.,T0UV_cgs=12000.,cUV_cgs=5.e-4
  double precision, parameter, dimension(6) ::  &
  coolT_cgs=(/ 300.D0,     2000.D0,    8000.D0,    1.D5,    4.D7,     1.D9 /),  &
  coolH_cgs=(/ 2.2380D-32, 1.0012D-30, 4.6240D-36, 1.7800D-18, 3.2217D-27, 0.D0   /)
  
  real :: rhoUV,TUV,T0UV,cUV
  real, dimension(6) :: coolT, &
    coolB=(/ 2.,       1.5,      2.867,    -.65,    0.5,      0.   /)
  double precision, dimension(6) :: coolH 

  integer :: iproc_SN,ipy_SN,ipz_SN
  logical :: ltestSN = .false.  ! If set .true. SN are only exploded at the
                              ! origin and ONLY the type I scheme is used
                              ! Used in kompaneets test etc.

  ! Should maybe be nondimensionalised
  real :: tau_cloud=2e-2 
  real, parameter :: rho_min=1.e-6
  !real, parameter :: tosolarMkpc3=1.483e7
  real, parameter :: frac_converted=0.02,frac_heavy=0.10,mass_SN=10.

  ! input parameters
  real :: outer_shell_proportion = 2.
  real :: inner_shell_proportion = 1.5
  real :: coolingfunction_scalefactor=1.
  logical :: lnever_move_mass
!tony: disable SNII for debugging 
  logical :: lSNI=.true., lSNII=.false.

! Flag and explosion rates for average interstellar heating
  logical :: laverage_SN_heating = .false.
  logical :: lSN_eth=.true., lSN_ecr=.true.
  real :: r_SNI=3.e+4, r_SNII=4.e+3
  real :: frac_ecr=0.1, frac_eth=0.9

  real :: center_SN_x = impossible,  center_SN_y = impossible, center_SN_z = impossible 

  integer :: dummy 
  namelist /interstellar_init_pars/ dummy

  ! run parameters
  logical:: uniform_zdist_SNI = .false.
  namelist /interstellar_run_pars/ &
      ampl_SN,tau_cloud, &
      uniform_zdist_SNI, ltestSN, lnever_move_mass, &
      lSNI, lSNII, laverage_SN_heating, coolingfunction_scalefactor, &
      point_width, inner_shell_proportion, outer_shell_proportion, &
      center_SN_x, center_SN_y, center_SN_z, &
      frac_ecr, frac_eth, lSN_eth, lSN_ecr, &
      h_SNI, SNI_area_rate

  contains

!***********************************************************************
    subroutine register_interstellar()
!
!  19-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_interstellar called twice')
      first = .false.
!
      linterstellar = .true.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_interstellar: ENTER'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: interstellar.f90,v 1.80 2004-03-14 16:07:25 mee Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_interstellar: nvar > mvar')
      endif
!
    endsubroutine register_interstellar
!***********************************************************************
    subroutine initialize_interstellar(lstarting)
!
!  Perform any post-parameter-read initialization eg. set derived 
!  parameters
!
!  24-nov-02/tony: coded
!
!  read parameters from seed.dat and interstellar.dat
!
      use Cdata
      use General
      use Sub, only: inpui,inpup
      use Mpicomm, only: mpibcast_real, stop_it
      use Ionization, only: getmu
!
      logical, save :: first=.true.
      logical :: lstarting
      logical :: exist
      real :: mu,factor_mysterious

      if (first) then
         if (.not. lstarting) then
            call inpui(trim(directory)//'/seed.dat',seed,nseed)
            if (lroot.and.ip<14) then
               print*, 'initialize_interstellar: reading seed file'
               print*, 'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
            endif
            call random_seed_wrapper(put=seed(1:nseed))
         endif
!
!AB: comment please why interstellar should be read in and what it does
!
         if (lroot) then
            inquire(file=trim(datadir)//'/interstellar.dat',exist=exist)
            if (exist) then 
               if (lroot.and.ip<14) print*, 'initialize_interstellar: read interstellar.dat'
                call inpup(trim(datadir)//'/interstellar.dat',  &
                    interstellarsave,ninterstellarsave)
               if (lroot.and.ip<14) print*, 'initialize_interstellar: t_next_SNI', &
                    interstellarsave(1)
            else
               interstellarsave(1)=t_next_SNI
            endif

         endif
         call mpibcast_real(interstellarsave,1)
         t_next_SNI=interstellarsave(1)
      endif

      if (lroot.and.uniform_zdist_SNI) then
         print*,'initialize_interstellar: using UNIFORM z-distribution of SNI'
      endif
!
!  Calculate unit_Lambda; the following gives the same as
!  unit_Lambda =  unit_velocity**3 / (unit_density*unit_length)
!
      call getmu(mu) 
      if (unit_system=='cgs') then
        unit_Lambda = unit_energy * unit_velocity**3 * unit_time**2 * &
                       mu**2 * m_H**2
!unit_energy / (unit_density*unit_time*unit_mass)
      elseif (unit_system=='SI') then
        call stop_it('initialize_interstellar: SI unit conversions not implemented')
      endif
      if (lroot) print*,'initialize_interstellar: unit_Lambda',unit_Lambda
      coolH = coolH_cgs / unit_Lambda * coolingfunction_scalefactor
      coolT = coolT_cgs / unit_temperature

      if (unit_system=='cgs') then
        TT_SN_min=TT_SN_min_cgs / unit_temperature
        if (SNI_area_rate==impossible) SNI_area_rate=SNI_area_rate_cgs * unit_length**2 * unit_time
        if (h_SNI==impossible)         h_SNI=h_SNI_cgs / unit_length
        h_SNII=h_SNII_cgs / unit_length
        solar_mass=solar_mass_cgs / unit_mass
        rho_crit=rho_crit_cgs / unit_density
        TT_crit=TT_crit_cgs / unit_temperature
        if (ampl_SN == impossible) ampl_SN=ampl_SN_cgs / unit_energy 
      else
        call stop_it('initialize_interstellar: SI unit conversions not implemented')
      endif

      t_interval_SNI = SNI_area_rate * Lxyz(1) * Lxyz(2)


      if (lroot.and.ip<14) then
        print*,'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
        print*,'initialize_interstellar: finished'
      endif

      if (lroot.and. (.not. lstarting)) then
         open(1,file=trim(datadir)//'/sn_series.dat',position='append')
         write(1,'("# ",A)')  &
          '--it-----t----------itype_SN---iproc_SN------x_SN-----------y_SN-----------z_SN-----------rho_SN---------EE_SN-----l_SN--m_SN--n_SN-----'
         close(1)
      endif

      if (ltestSN) then
        t_interval_SNI=1.E10
        t_next_SNI=0.
      endif
!
    endsubroutine initialize_interstellar
!***********************************************************************
    subroutine calc_heat_cool_interstellar(df,rho1,TT1,yH)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  We may want to move it to the entropy module for good, because its use
!  is not restricted to interstellar runs (could be used for solar corona).
!  Also, it doesn't pose an extra load on memory usage or compile time.
!  (We should allow that UV heating can be turned off; so rhoUV should
!  be made an input parameter.)
!
!  19-nov-02/graeme: adapted from calc_heat_cool
!  10-aug-03/axel: TT is used as input
!
      use Cdata
      use Mpicomm
      use Density, only : rho0
      use Sub
      use Ionization
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx), intent(in) :: rho1,TT1,yH
      real, dimension (nx) :: heat,cool,rho,TT
      real :: norm
      integer :: i
!
!  identifier
!
      if(headtt) print*,'calc_heat_cool_interstellar: ENTER'
!
!  rho factor (could perhaps better be calculated in entropy)
!
      rho=1./rho1
      TT=1./TT1
!
!  define T in K, for calculation of both UV heating and radiative cooling
!
!  add T-dept radiative cooling, from Rosen et al., ApJ, 413, 137, 1993
!  cooling is Lambda*rho^2, with (eq 7)
!     Lambda=coolH(i)*TT*coolB(i),   for coolT(i) <= T < coolT(i+1)
!  nb: our coefficients coolH(i) differ from those in Rosen et al. by
!   factor (mu mp)^2, with mu=1.2, since Rosen works in number density, n.
!   (their cooling = Lambda*n^2,  rho=mu mp n.)
!  The factor Lambdaunits converts from cgs units to code units.
!
!  [Currently, coolT(1) is not modified, but this may be necessary
!  to avoid creating gas too cold to resolve.]
!
      cool=0.0
      do i=1,5
        where (coolT(i) <= TT .and. TT < coolT(i+1)) &
          cool=cool+yH*coolH(i)*rho**2*(TT*unit_temperature)**coolB(i)
      enddo
!
!  add UV heating, cf. Wolfire et al., ApJ, 443, 152, 1995
!  with the values above, this gives about 0.012 erg/g/s (T < ~1.e4 K)
!  nb: need rho0 from density_[init/run]_pars, if i want to implement
!      the the arm/interarm scaling.
!
      heat=0.0
!tony: DISABLE UV HEATING -- Requires reformulation
!   heat(:)=rhoUV*(rho0/1.38)**1.4*unit_Lambda*coolH(3)*TUV**coolB(3)*   &
!                               0.5*(1.0+tanh(cUV*(T0UV-TT(:))))
!
!tony: need to do unit_system stuff with scale heights etc.
!  Average SN heating (due to SNI and SNII)
!  The amplitudes of both types is assumed the same (=ampl_SN)
!
      if (laverage_SN_heating) then
        norm=ampl_SN/sqrt(2*pi)
        heat=heat+r_SNI *norm/h_SNI *exp(-(z(n)/h_SNI )**2)
        heat=heat+r_SNII*norm/h_SNII*exp(-(z(n)/h_SNII)**2)
      endif
!
!  For clarity we have constructed the rhs in erg/s/volume [=rho*T*Ds/Dt]
!  so therefore we now need to multiply by rho1*TT1.
!
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+rho1*TT1*(heat-cool)
!
    endsubroutine calc_heat_cool_interstellar
!***********************************************************************
    subroutine check_SN(f)
!
!  Checks for SNe, and implements appropriately:
!   relevant subroutines in entropy.f90
!
    use Cdata
!
    real, dimension(mx,my,mz,mvar+maux) :: f
    logical :: l_SNI=.false.   !only allow SNII if no SNI this step
                               !(may not be worth keeping)
!
    intent(inout) :: f
!
!  identifier  
!
      if(headtt) print*,'check_SN: ENTER'
!
!  Do separately for SNI (simple scheme) and SNII (Boris' scheme)
!
    if (lSNI) call check_SNI (f,l_SNI)
    if (lSNII) call check_SNII(f,l_SNI)
!
    endsubroutine check_SN
!***********************************************************************
    subroutine check_SNI(f,l_SNI)
!
!  If time for next SNI, then implement, and calculate time of subsequent SNI
!
    use Cdata
    use Mpicomm
    use General
!
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(1) :: franSN
    logical :: l_SNI
!
    intent(inout) :: f,l_SNI
!
!  identifier
!
    if(headtt) print*,'check_SNI: ENTER'
!
    l_SNI=.false.
    if (t >= t_next_SNI) then
       call position_SNI(f)
       call explode_SN(f,1)
       !  pre-determine time for next SNI
       if (lroot) then
          if (ip<14) print*,"check_SNI: Old t_next_SNI=",t_next_SNI
          call random_number_wrapper(franSN)   
          t_next_SNI=t + (1.0 + 0.4*(franSN(1)-0.5)) * t_interval_SNI
          if (ip<14) print*,'check_SNI: Next SNI at time = ',t_next_SNI
          interstellarsave(1)=t_next_SNI
          if (ip<14) print*,"check_SNI: New t_next_SNI=",t_next_SNI
       endif
       call mpibcast_real(t_next_SNI,1)
       l_SNI=.true.
    endif
!
    endsubroutine check_SNI
!***********************************************************************
    subroutine check_SNII(f,l_SNI)
!
!  Check for SNII, via self-regulating scheme.
!
    use Cdata
    use General
    use Mpicomm
    use Ionization
! 
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(nx) :: lnrho,rho,rho_cloud,ss,lnTT,TT,yH
!    real :: lnrho,rho,rho_cloud,ss,TT
    real :: mass_cloud,mass_cloud_dim,freq_SNII,prob_SNII,rate_SNII,dv
    real, dimension(1) :: franSN,fsum1,fsum1_tmp,fmpi1
    real, dimension(ncpus) :: mass_cloud_byproc
    integer :: icpu
    logical :: l_SNI
!
    intent(in) :: l_SNI
    intent(inout) :: f
!
!  identifier
!
    if(headtt.and.ip<14) print*,'check_SNII: ENTER'
!
    if (.not. l_SNI) then         ! only do if no SNI this step
!
!  NB: currently no 'nzskip' mechanism to prevent SNII occurring near
!   top or bottom of box.  Non-trivial to implement with nprocz > 1 -- and
!   shouldn't be a real problem, as mass near boundaries should be low.
       mass_cloud=0.0
       do n=n1,n2
           do m=m1,m2
             rho(:)=exp(f(l1:l2,m,n,ilnrho))
             ss(:)=f(l1:l2,m,n,iss)

             call eoscalc(f,yH=yH,lnTT=lnTT)
             TT=exp(lnTT)

             rho_cloud(:)=0.0
             !print*,'min,max TT:  ', min(TT), max(TT)
             !print*,'min,max rho: ', min(rho), max(rho)

             where (rho(:) >= rho_crit .and. TT(:) <= TT_crit)   &
                  rho_cloud(:)=rho(:)
             mass_cloud=mass_cloud+sum(rho_cloud(:))
          enddo
       enddo
       fsum1_tmp=(/ mass_cloud /)
       !print*,'check_SNII, iproc,fsum1_tmp:',iproc,fsum1_tmp(1)
       call mpireduce_sum(fsum1_tmp,fsum1,1) 
       call mpibcast_real(fsum1,1)
       dv=1.
       if (nxgrid/=1) dv=dv*dx
       if (nygrid/=1) dv=dv*dy
       if (nzgrid/=1) dv=dv*dz
       mass_cloud_dim=fsum1(1)*dv/solar_mass
       !print*,'check_SNII: iproc,fsum1:',iproc,fsum1(1)
       ! need convert to dimensional units, for rate/probability calculation only. 
       ! don't overwrite mass_cloud (on individual processors), as it's re-used.
       !if (lroot .and. ip < 14) &
       !     print*,'check_SNII: mass_cloud_dim:',mass_cloud_dim
       !
       freq_SNII=frac_heavy*frac_converted*mass_cloud_dim/mass_SN/tau_cloud
       prob_SNII=freq_SNII*dt
       rate_SNII=freq_SNII*1e-3
       if (Lxyz(1)/=0.) rate_SNII=rate_SNII/Lxyz(1)
       if (Lxyz(2)/=0.) rate_SNII=rate_SNII/Lxyz(2)
       if (lroot) call random_number_wrapper(franSN)   
       call mpibcast_real(franSN,1)
        if (lroot .and. ip < 16) &
             print*,'check_SNII: rate,prob,rnd:',rate_SNII,prob_SNII,franSN(1)
       if (franSN(1) <= prob_SNII) then
          !  position_SNII needs the mass_clouds for each processor;  
          !   communicate and store them here, to avoid recalculation.
          mass_cloud_byproc(:)=0.0
          ! use non-root broadcasts for the communication...
          do icpu=1,ncpus
             fmpi1=mass_cloud
             call mpibcast_real_nonroot(fmpi1,1,icpu-1)
             mass_cloud_byproc(icpu)=fmpi1(1)
          enddo
          ! if (lroot.and.ip<14) print*,'check_SNII: mass_cloud_byproc:',mass_cloud_byproc
          call position_SNII(f,mass_cloud_byproc)
          call explode_SN(f,2)
       endif
    endif
    !
  endsubroutine check_SNII
!***********************************************************************
    subroutine position_SNI(f)
!
!   determine position for next SNI (w/ fixed scale-height)
!
    use Cdata
    use Mpicomm
    use General
!
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
!
    real, dimension(nzgrid) :: cum_prob_SNI
    real :: zn, z00, x00, y00
    real, dimension(3) :: fran3
    integer :: i, l, nzskip=10   !prevent SNI from being too close to boundaries
!
!
!
    if(headtt) print*,'position_SNI: ENTER'
    
    ! Calculate the global (nzgrid) lower z-coordinate
    if (lperi(1)) then; x00=xyz0(1)+.5*dx; else; x00=xyz0(1); endif
    if (lperi(2)) then; y00=xyz0(2)+.5*dy; else; y00=xyz0(2); endif
    if (lperi(3)) then; z00=xyz0(3)+.5*dz; else; z00=xyz0(3); endif
    
    
    !
    !  Pick SN position (l_SN,m_SN,n_SN)
    !
    call random_number(fran3)    ! get 3 random numbers
                                 ! on all processors to keep
                                 ! rnd. generators in sync
    if (lroot) then
       if (ltestSN) then
          if (center_SN_x.eq.impossible) then
            i=int(nxgrid/2)+1
          else
            i=int((center_SN_x-x00)/dx)+1
          endif
       else
          i=int(fran3(1)*nxgrid)+1
       endif
       l_SN=(i-1)+nghost

       if (ltestSN) then
          if (center_SN_y.eq.impossible) then
            i=int(nygrid/2)+1
          else
            i=int((center_SN_y-y00)/dy)+1
          endif
       else
          i=int(fran3(2)*nygrid)+1
       endif
       ipy_SN=(i-1)/ny  ! uses integer division
       m_SN=(i-1)-(ipy_SN*ny)+nghost

       if (ltestSN) then
          if (center_SN_z.eq.impossible) then
            i=int(nzgrid/2)+1
          else
            i=int((center_SN_z-z00)/dz)+1
          endif
          ipz_SN=(i-1)/nz   ! uses integer division
          n_SN=(i-1)-(ipz_SN*nz)+nghost
       elseif (uniform_zdist_SNI) then
          i=int(fran3(3)*nzgrid)+1
          ipz_SN=(i-1)/nz   ! uses integer division
          n_SN=(i-1)-(ipz_SN*nz)+nghost
       else
       !
       !  Cumulative probability function in z currently calculated each time.
       !  It's constant, and could be stored (and calculated in init)
          cum_prob_SNI(1:nzskip)=0.0
          do n=nzskip+1,nzgrid-nzskip
             zn=z00+(n-1)*dz
             cum_prob_SNI(n)=cum_prob_SNI(n-1)+exp(-(zn/h_SNI)**2)
          enddo
          cum_prob_SNI=cum_prob_SNI/cum_prob_SNI(nzgrid-nzskip)
       !  The following should never be needed, but just in case floating point 
       !  errors ever lead to cum_prob_SNI(nzgrid-nzskip) < rnd < 1.
          cum_prob_SNI(nzgrid-nzskip+1:nzgrid)=1.0   
          
          do i=nzskip+1,nzgrid-nzskip
             if (cum_prob_SNI(i-1) <= fran3(3) .and. fran3(3) < cum_prob_SNI(i)) &
                  then
                ipz_SN=(i-1)/nz  ! uses integer division
                n_SN=(i-1)-(ipz_SN*nz)+nghost
                exit
             endif
          enddo
       endif
       iproc_SN=ipz_SN*nprocy + ipy_SN
    endif

    call share_SN_parameters(f)

    endsubroutine position_SNI
!***********************************************************************
    subroutine position_SNII(f,mass_cloud_byproc)
!
!  Determine position for next SNII (using Boris' scheme)
!  It seems impractical to sort all high density points across all processors;
!  instead, we just construct cumulative pdfs that allow us to pick a processor,
!  and then a point on that processor, with probability proportional to rho.
!  As a result, the SN position is *not* independent of ncpus (or of nprocy 
!  and nprocz).  (It is repeatable given fixed nprocy/z though.)
!
    use Cdata
    use General
    use Mpicomm
    use Ionization, only: eoscalc,ilnrho_ss
!
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
    real, intent(in) , dimension(ncpus) :: mass_cloud_byproc
!
    real, dimension(0:ncpus) :: cum_prob_byproc
    real, dimension(1) :: franSN
    real :: mass_cloud,cum_mass,cum_prob_onproc
    real :: lnrho,rho,ss,lnTT,TT,yH
    integer :: icpu,l
!
!
!  identifier
!
      if(lroot.and.ip<20) print*,'position_SNII: ENTER'
!
!  Construct cumulative distribution function, using mass_cloud_byproc.
!  NB: icpu=iproc+1 (iproc in [0,ncpus-1], icpu in [1,ncpus] )
!
    cum_prob_byproc=0.0
    do icpu=1,ncpus
      mass_cloud=mass_cloud_byproc(icpu)
      cum_prob_byproc(icpu)=cum_prob_byproc(icpu-1)+mass_cloud_byproc(icpu)
    enddo
    cum_prob_byproc(:)=cum_prob_byproc(:)/cum_prob_byproc(ncpus)
    if (lroot.and.ip<14) then
      print*,'position_SNII: mass_cloud_byproc=',mass_cloud_byproc
      print*,'position_SNII: cum_prob_byproc=',cum_prob_byproc
      print*,'position_SNII: mass_cloud=',mass_cloud
    endif
!     
!  Use random number to detemine which processor SN is on.
!  (Use root processor for rand, to ensure repeatability.)
!
    if (lroot) call random_number_wrapper(franSN)   
    call mpibcast_real(franSN,1)
    do icpu=1,ncpus
      if (cum_prob_byproc(icpu-1) <= franSN(1) .and.                      &
           franSN(1) < cum_prob_byproc(icpu)) then
        iproc_SN=icpu-1 
        exit
      endif
    enddo
    if (lroot.and.ip<20) &
          print*, 'position_SNII: franSN(1),iproc_SN=',franSN(1),iproc_SN
!
!  Use random number to pick SNII location on the right processor.
!  (No obvious reason to re-use the original random number for this.)
!    franSN(1)=(franSN(1)-cum_prob_byproc(iproc_SN)) /                      &
!              (cum_prob_byproc(iproc_SN+1)-cum_prob_byproc(iproc_SN))
!
    if (lroot) call random_number_wrapper(franSN)   
    call mpibcast_real(franSN,1)
    if (iproc == iproc_SN) then
      cum_mass=0.0
find_SN: do n=n1,n2
        do m=m1,m2
          do l=l1,l2
            lnrho=f(l,m,n,ilnrho)
            rho=exp(lnrho)
            ss=f(l,m,n,iss)
            call eoscalc(ilnrho_ss,lnrho,ss,yH=yH,lnTT=lnTT)
            TT=exp(lnTT)
            if (rho >= rho_crit .and. TT <= TT_crit) then
              cum_mass=cum_mass+rho
              cum_prob_onproc=cum_mass/mass_cloud
              if (franSN(1) <= cum_prob_onproc) then
                l_SN=l; m_SN=m; n_SN=n
                if (ip<20) &
                 print*,'position_SNII: cum_mass,cum_prob_onproc,franSN(1)=', &
                                  cum_mass,cum_prob_onproc,franSN(1)
                exit find_SN
              endif
            endif
          enddo
        enddo
      enddo find_SN
    endif
!
    call share_SN_parameters(f)
!
    endsubroutine position_SNII
!***********************************************************************
    subroutine share_SN_parameters(f)
!   
!   Handle common SN positioning processor communications
!
!
!   27-aug-2003/tony: coded
!    
    use Mpicomm
    use Ionization
      
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
    
    real, dimension(5) :: fmpi5
    integer, dimension(4) :: impi4
!
!  Broadcast position to all processors from root;
!  also broadcast iproc_SN, needed for later broadcast of rho_SN.
!
!
    impi4=(/ iproc_SN, l_SN, m_SN, n_SN /)
    call mpibcast_int(impi4,4)
    iproc_SN=impi4(1)
    l_SN=impi4(2)
    m_SN=impi4(3)
    n_SN=impi4(4)

!
!  With current SN scheme, we need rho at the SN location.
!
      
    if (iproc==iproc_SN) then
      lnrho_SN=f(l_SN,m_SN,n_SN,ilnrho)
      ss_SN=f(l_SN,m_SN,n_SN,iss)
      x_SN=0.; y_SN=0.; z_SN=0.
      if (nxgrid/=1) x_SN=x(l_SN)
      if (nygrid/=1) y_SN=y(m_SN)
      if (nzgrid/=1) z_SN=z(n_SN)
    endif
!
!  Broadcast to all processors.
!
    fmpi5=(/ x_SN, y_SN, z_SN, lnrho_SN, ss_SN /)
    call mpibcast_real_nonroot(fmpi5,5,iproc_SN)

    x_SN=fmpi5(1); y_SN=fmpi5(2); z_SN=fmpi5(3); 
    lnrho_SN=fmpi5(4); ss_SN=fmpi5(5)

    rho_SN=exp(lnrho_SN);

    if (lroot.and.ip<=14) print*, &
 'share_SN_parameters: iproc_SN,x_SN,y_SN,z_SN,l_SN,m_SN,n_SN,rho_SN,ss_SN = ' &
          ,iproc_SN,x_SN,y_SN,z_SN,l_SN,m_SN,n_SN,rho_SN,ss_SN
!
    endsubroutine share_SN_parameters
!***********************************************************************
    subroutine explode_SN(f,itype_SN)
      !
      !  Implement SN (of either type), at pre-calculated position
      !  (This can all be made more efficient, after debugging.)
      !
      !  ??-nov-02/grs : coded from GalaxyCode                        
      !  20-may-03/tony: pencil formulation and broken into subroutines
      !
      use Cdata
      use Mpicomm
      use Ionization
      !
      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
      integer, intent(in) :: itype_SN

      real :: width_SN,width_shell_outer,width_shell_inner,c_SN
      real :: profile_integral, mass_shell, mass_gain
      real :: EE_SN=0.,EE2_SN=0.
      real :: rho_SN_new,lnrho_SN_new,ss_SN_new,yH_SN_new,lnTT_SN_new,ee_SN_new,TT_SN_new,dv
      
      real, dimension(nx) :: deltarho, deltaEE
      real, dimension(1) :: fmpi1, fmpi1_tmp
      real, dimension(2) :: fmpi2, fmpi2_tmp
      real, dimension(nx) ::  lnrho, ss, yH, lnTT, TT, rho_old, ee_old

      logical :: lmove_mass=.false.
      integer :: idim
          
      !
      !  identifier
      !
      if(lroot.and.ip<12) print*,'explode_SN: itype_SN=',itype_SN
      !
      width_SN=point_width*dxmin      
      idim=0                         !allow for 1-d, 2-d and 3-d profiles...
      if (nxgrid /=1) idim=idim+1
      if (nygrid /=1) idim=idim+1
      if (nzgrid /=1) idim=idim+1
      c_SN=ampl_SN/cnorm_SN(idim)/width_SN**idim
      dv=1.
      if (nxgrid/=1) dv=dv*dx
      if (nygrid/=1) dv=dv*dy
      if (nzgrid/=1) dv=dv*dz

      if (lroot.and.ip<=14) print*,'explode_SN: width_SN,c_SN,rho_SN=', width_SN,c_SN,rho_SN
        
      !
      !  Now deal with (if nec.) mass relocation
      !
      call eoscalc(ilnrho_ss,lnrho_SN,ss_SN, &
                              yH=yH_SN,lnTT=lnTT_SN,ee=ee_SN)
      TT_SN=exp(lnTT_SN)

      ee_SN_new = frac_eth*(ee_SN+c_SN/rho_SN)
      call eoscalc(ilnrho_ee,lnrho_SN,ee_SN_new, &
                              ss=ss_SN_new,lnTT=lnTT_SN_new,yH=yH_SN_new)
      TT_SN_new=exp(lnTT_SN_new)

      if(lroot.and.ip<=14) print*, &
         'explode_SN: TT_SN, TT_SN_new, TT_SN_min, ee_SN =', &
                                TT_SN,TT_SN_new,TT_SN_min, ee_SN

      if (TT_SN_new < TT_SN_min) then
         lmove_mass=.not.lnever_move_mass
         ! lmove_mass=.false.  ! use to switch off for debug...

         ! The bit that BREAKS the pencil formulation...
         ! must know the total moved mass BEFORE attempting mass relocation 

         ! ASSUME: SN will fully ionize the gas at its centre
         if (lmove_mass) then
           call getdensity((ee_SN*rho_SN)+c_SN,TT_SN_min,1.,rho_SN_new)
           lnrho_SN_new=log(rho_SN_new)
           ee_SN_new=frac_eth*(ee_SN*rho_SN+c_SN)/rho_SN_new

           call eoscalc(ilnrho_ee,lnrho_SN_new,ee_SN_new, &
                                 ss=ss_SN_new,lnTT=lnTT_SN_new,yH=yH_SN_new)
           TT_SN_new=exp(lnTT_SN_new)

           if(lroot.and.ip<=14) print*, &
              'explode_SN: Relocate mass... TT_SN_new, rho_SN_new=', &
                                                     TT_SN_new,rho_SN_new

           call calcmassprofileintegral_SN(f,width_SN,profile_integral)
           fmpi1_tmp=(/ profile_integral /)
           call mpireduce_sum(fmpi1_tmp,fmpi1,1) 
           call mpibcast_real(fmpi1,1)
           profile_integral=fmpi1(1)*dv
           mass_shell=-(rho_SN_new-rho_SN)*profile_integral
           if (lroot.and.ip<=14) &
             print*, 'explode_SN: mass_shell=',mass_shell
           mass_gain=0.
         endif
      endif
      

      EE_SN=0. ! EE_SN2=0.
      do n=n1,n2
         do m=m1,m2
            
            ! Calculate the distances to the SN origin for all points
            ! in the current pencil and store in the dr2_SN global array
            call proximity_SN()
            
            ! Get the old energy
            lnrho=f(l1:l2,m,n,ilnrho)
            rho_old=exp(lnrho) 
            ss=f(l1:l2,m,n,iss)
            call eoscalc(f,yH=yH,lnTT=lnTT,ee=ee_old)
            TT=exp(lnTT)

            ! Apply perturbations
            call injectenergy_SN(deltaEE,width_SN,c_SN,EE_SN)
            if (lmove_mass) then
              call makecavity_SN(deltarho,width_SN,rho_SN_new-rho_SN, &
                      mass_shell,cnorm_SN(idim),idim,mass_gain)

              lnrho=alog(amax1(rho_old(:)+deltarho(:),rho_min))
            endif
  
            call perturb_energy(lnrho,(ee_old*rho_old+deltaEE*frac_eth)/exp(lnrho),ss,lnTT,yH)
            TT=exp(lnTT)

            if (lcosmicray.and.lSN_ecr) then 
              f(l1:l2,m,n,iecr) = f(l1:l2,m,n,iecr) + (deltaEE * frac_ecr) 
            endif
            !call eoscalc(ilnrho_ee,lnrho,ee_old+(deltaEE/rho_old),ss=ss,lnTT=lnTT,yH=yH)

            ! Save changes 
            if (lSN_eth) then
              f(l1:l2,m,n,ilnrho)=lnrho
              f(l1:l2,m,n,iss)=ss
            endif

            lnTT=log(TT)
            if (ilnTT.ne.0) f(l1:l2,m,n,ilnTT)=lnTT
            if (iyH.ne.0) f(l1:l2,m,n,iyH)=yH
!

       enddo
      enddo

      ! Sum and share diagnostics etc. amongst processors
      fmpi2_tmp=(/ mass_gain, EE_SN /)
      call mpireduce_sum(fmpi2_tmp,fmpi2,4) 
      call mpibcast_real(fmpi2,4)
      mass_gain=fmpi2(1)*dv
      EE_SN=fmpi2(2)*dv
! Extra debug - no longer calculated 
!      EE2_SN=fmpi3(3)*dv; 

      if (lroot.and.ip<=14) print*, &
           'explode_SN: mass_gain=',mass_gain
     
      if (lroot) then
         open(1,file=trim(datadir)//'/sn_series.dat',position='append')
         write(1,'(1i6,1e11.3,2i10,"    ",5e15.7,3i5)')  &
                          it,t,itype_SN,iproc_SN,x_SN,y_SN,z_SN,rho_SN,EE_SN,l_SN,m_SN,n_SN
         close(1)
      endif
      
    endsubroutine explode_SN

!***********************************************************************
    subroutine calcmassprofileintegral_SN(f,width_SN,profile_integral)
!
!  Calculate integral of mass cavity profile  
!
!  22-may-03/tony: coded
!
      use Cdata

      real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
      real, intent(in) :: width_SN
      real, intent(out) :: profile_integral
      real :: dx_SN_in,dx_SN_out_x0,dx_SN_out_x1,dy_SN_in,dy_SN_out_y
      real :: dy_SN_out_x0a,dy_SN_out_x0b,dy_SN_out_x1a,dy_SN_out_x1b
      real :: dz_SN,yshift
      integer :: l,mshift,il,im,in
     
      !
      !  Obtain distance to SN
      !

      profile_integral=0.
      do n=n1,n2
         in=n-nghost
         dz_SN=abs(z(n)-z_SN)
         do m=m1,m2
            im=m-nghost
            !  consider all possible positions in xy plane, to get shortest
            !  can be made more efficient later
            !  y-separations with no boundaries crossed in x
            dy_SN_in=abs(y(m)-y_SN)                      !dyi
            dy_SN_out_y=Ly-dy_SN_in                      !dyii
            !  y-separations across sliding periodic boundary at x=x0
            yshift=y(m)-deltay
            mshift=0
            if (yshift < y0) then 
               mshift=int((y0-yshift)/Ly)+1
               yshift=yshift+mshift*Ly
            elseif (yshift > y0+Ly) then      
               mshift=int((yshift-(y0+Ly))/Ly)+1
               yshift=yshift-mshift*Ly
            endif
            dy_SN_out_x0a=abs(yshift-y_SN)               !dyiii
            dy_SN_out_x0b=Ly-dy_SN_out_x0a               !dyiv
            !  y-separations across sliding periodic boundary at x=x1
            yshift=y(m)+deltay
            mshift=0
            if (yshift < y0) then 
               mshift=int((y0-yshift)/Ly)+1
               yshift=yshift+mshift*Ly
            elseif (yshift > y0+Ly) then
               mshift=int((yshift-(y0+Ly))/Ly)+1
               yshift=yshift-mshift*Ly
            endif
            dy_SN_out_x1a=abs(yshift-y_SN)               !dyv
            dy_SN_out_x1b=Ly-dy_SN_out_x1a               !dyvi
            do l=l1,l2
               il=l-nghost
               !  x-separations associated with each of the above
               dx_SN_in=abs(x(l)-x_SN)                    !dxi=dxii
               dx_SN_out_x0=Lx+(x(l)-x_SN)                !dxiii=dxiv
               dx_SN_out_x1=Lx-(x(l)-x_SN)                !dxv=dxvi
               dr2_SN(il)=min( dx_SN_in**2 + dy_SN_in**2,             &
                    dx_SN_in**2 + dy_SN_out_y**2,          &
                    dx_SN_out_x0**2 + dy_SN_out_x0a**2,    &
                    dx_SN_out_x0**2 + dy_SN_out_x0b**2,    &
                    dx_SN_out_x1**2 + dy_SN_out_x1a**2,    &
                    dx_SN_out_x1**2 + dy_SN_out_x1b**2 )   &
                    + dz_SN**2
            enddo
            profile_integral = profile_integral + sum(exp(-(dr2_SN(:)/width_SN**2)**3))
         enddo
      enddo

    endsubroutine calcmassprofileintegral_SN
!***********************************************************************
    subroutine proximity_SN()
!
!  Calculate pencil of distance to SN explosion site
!
!  20-may-03/tony: extracted from explode_SN code written by grs
!  22-may-03/tony: pencil formulation
!
!
      use Cdata

      real :: dx_SN_in,dx_SN_out_x0,dx_SN_out_x1,dy_SN_in,dy_SN_out_y
      real :: dy_SN_out_x0a,dy_SN_out_x0b,dy_SN_out_x1a,dy_SN_out_x1b
      real :: dz_SN,yshift
      integer :: l,mshift,il,im,in
     
      !
      !  Obtain distance to SN
      !
         im=m-nghost
         in=n-nghost
         dz_SN=abs(z(n)-z_SN)
         !  consider all possible positions in xy plane, to get shortest
         !  can be made more efficient later
         !  y-separations with no boundaries crossed in x
         dy_SN_in=abs(y(m)-y_SN)                      !dyi
         dy_SN_out_y=Ly-dy_SN_in                      !dyii
         !  y-separations across sliding periodic boundary at x=x0
         yshift=y(m)-deltay
         mshift=0
         if (yshift < y0) then 
            mshift=int((y0-yshift)/Ly)+1
            yshift=yshift+mshift*Ly
         elseif (yshift > y0+Ly) then      
            mshift=int((yshift-(y0+Ly))/Ly)+1
            yshift=yshift-mshift*Ly
         endif
         dy_SN_out_x0a=abs(yshift-y_SN)               !dyiii
         dy_SN_out_x0b=Ly-dy_SN_out_x0a               !dyiv
         !  y-separations across sliding periodic boundary at x=x1
         yshift=y(m)+deltay
         mshift=0
         if (yshift < y0) then 
            mshift=int((y0-yshift)/Ly)+1
            yshift=yshift+mshift*Ly
         elseif (yshift > y0+Ly) then
            mshift=int((yshift-(y0+Ly))/Ly)+1
            yshift=yshift-mshift*Ly
         endif
         dy_SN_out_x1a=abs(yshift-y_SN)               !dyv
         dy_SN_out_x1b=Ly-dy_SN_out_x1a               !dyvi
         do l=l1,l2
            il=l-nghost
            !  x-separations associated with each of the above
            dx_SN_in=abs(x(l)-x_SN)                    !dxi=dxii
            dx_SN_out_x0=Lx+(x(l)-x_SN)                !dxiii=dxiv
            dx_SN_out_x1=Lx-(x(l)-x_SN)                !dxv=dxvi
            dr2_SN(il)=min( dx_SN_in**2 + dy_SN_in**2,             &
                 dx_SN_in**2 + dy_SN_out_y**2,          &
                 dx_SN_out_x0**2 + dy_SN_out_x0a**2,    &
                 dx_SN_out_x0**2 + dy_SN_out_x0b**2,    &
                 dx_SN_out_x1**2 + dy_SN_out_x1a**2,    &
                 dx_SN_out_x1**2 + dy_SN_out_x1b**2 )   &
                 + dz_SN**2
         enddo
       
    endsubroutine proximity_SN
!***********************************************************************
    subroutine makecavity_SN(deltarho,width_SN,depth,mass_shell, &
                             cnorm_dim,idim,mass_gain)   
      use Cdata
      !
      real, intent(in) :: width_SN, depth, mass_shell, cnorm_dim
      real, intent(inout) :: mass_gain
      real, intent(out), dimension(nx) :: deltarho
      integer, intent(in) :: idim
      !
      real, dimension(nx) :: profile_shell_outer,profile_shell_inner
      real, dimension(1) :: fsum1,fsum1_tmp
      real :: width_shell_outer, width_shell_inner, c_shell
      ! real, dimension(1) :: fsum1,fsum1_tmp

      width_shell_outer=outer_shell_proportion*width_SN
      width_shell_inner=inner_shell_proportion*width_SN

      deltarho(:) =  depth*exp(-(dr2_SN(:)/width_SN**2)**3)
      
      c_shell=mass_shell /                                  &
           (cnorm_dim*((width_shell_outer**idim)-(width_shell_inner**idim)))

!      if (lroot) print*, &
!           'explode_SN, c_shell:',c_shell
      !  add missing mass back into shell
      
      profile_shell_outer(:)=                              &
           exp(-(dr2_SN(:)/width_shell_outer**2)**3)
      profile_shell_inner(:)=                              &
           exp(-(dr2_SN(:)/width_shell_inner**2)**3)
      
      deltarho(:)=deltarho(:) + c_shell *      &
           (profile_shell_outer(:) - profile_shell_inner(:))
      mass_gain=mass_gain + sum(deltarho(:))  
      
    endsubroutine makecavity_SN

!***********************************************************************
    subroutine injectenergy_SN(deltaEE,width_SN,c_SN,EE_SN)
      use Cdata
      !
      real, intent(in) :: width_SN,c_SN
      real, intent(inout) :: EE_SN
      real, intent(out), dimension(nx) :: deltaEE
      !
      real, dimension(nx) :: profile_SN
      
      ! Whether mass moved or not, inject energy.
      !

      profile_SN=exp(-(dr2_SN(:)/width_SN**2)**3)

      deltaEE(:)=c_SN*profile_SN(:) ! spatial energy density 
      EE_SN=EE_SN+sum(deltaEE(:))   

    endsubroutine injectenergy_SN

endmodule interstellar
!***********************************************************************
