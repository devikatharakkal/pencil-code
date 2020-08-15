!  This module generates the initial conditions for simulating a section of
!  a galaxy disk.
!
!  Assumes that the calculation will be done in the system of units
!  described in Fred Gent's thesis.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
module InitialCondition

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages

  implicit none

  include '../initial_condition.h'
!
! Parameters
  real :: Sigma_g = 10          ! Surface density of gas [Msun pc^-2]
  real :: scaleheight = 0.2    ! scaleheight of the disk [kpc]
  real :: Pmag_to_Ptherm = 0.0 ! Ratio between thermal and magnetic pressure
  real :: Pcr_to_Ptherm = 0.0 ! Ratio between thermal and cosmic ray pressure
  real :: fgas = 0.1           ! Gas fraction (Sigma_g/Sigma_total)
  real :: initial_uu_noise=0.0 ! Amplitude of the noise in uu at the midplane
! Module variables (may be removed later...)
  real :: TTi ! Initial gas temperature ["pencil temperatures"]
  real :: Sigma_t ! Total surface density (stars+DM+gas) [Msun pc^-2]
  real :: pc_Sigma_t ! Total surface density (stars+DM+gas) [Msun pc^-2]
  real :: pc_Sigma_g ! Total surface density (stars+DM+gas) [Msun pc^-2]
! Constants
  real, parameter :: G_pencil_units = 66.42226827645213423329
  real, parameter :: surface_density_pencil_unit = &
                   15.80418771546467448018 ! M_pc/Msun /(1 kpc^2 / pc^2)
  real  :: gammacr=4./3. ! TODO Use sharedvariables instead of this

  namelist /initial_condition_pars/ &
      Sigma_g, fgas, scaleheight, Pmag_to_Ptherm, Pcr_to_Ptherm,  &
      initial_uu_noise

  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
         "$Id: noinitial_condition.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!       use SharedVariables, only: get_shared_variable
      use EquationOfState,only: gamma_m1, gamma1, get_cp1,cs20
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ierr
      real    :: cs20_recommended, TT0, gravz_recommended, TT0_actual

      Sigma_t = Sigma_g/fgas
      pc_Sigma_t = Sigma_t / surface_density_pencil_unit
      pc_Sigma_g = Sigma_g / surface_density_pencil_unit

!     Computes the recommended initial sound speed (to achieve hydrostatic
!      equilibrium)
      cs20_recommended = pi*G_pencil_units*pc_Sigma_t*scaleheight/gamma1  &
                         /(1.0+Pmag_to_Ptherm+Pcr_to_Ptherm)
      ! Computes the TT0 associated with this
      TT0 = gamma1 * cs20_recommended ! Assuming mu=R=1
      ! and what we have
      TT0_actual =  gamma1 * cs20! Assuming mu=R=1

      ! Computes the recommended
      gravz_recommended = 2.0 * pi * G_pencil_units *    pc_Sigma_t
      if (lroot) then
        print *, &
     'Disk section initial condition: The following parameters had been used'
        print *, '   - scaleheight',scaleheight
        print *, '   - fgas',fgas
        print *, '   - pc_Sigma_t',pc_Sigma_t
        print *, 'Disk section initial condition: gravz should be set to ',&
                  gravz_recommended
        print *, 'Disk section initial condition: cs0 should be set to ', &
                  sqrt(cs20_recommended)
        print *, 'Disk section initial condition: present temperature = ', &
                  TT0_actual, ' recomended = ', TT0
        print *, &
        '                           in kelvins:   present temperature= ', &
                  TT0_actual*44.74127*(5./3.)*gamma1/0.62, &
                  'K  recomended = ', TT0*44.74127*(5./3.)*gamma1/0.62
      endif
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!  15-feb-15/MR: optional parameter 'profiles' added
!
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
      if (present(profiles)) then
        call fatal_error('initial_condition_all', &
          'If profiles are asked for, a real initial condition must be specified')
        call keep_compiler_quiet(profiles)
      endif
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field with noise.
!  The amplitude of the superimposed noise is assumed to scale with the
!  density (i.e. those are actually perturbations to the _momentum_).
!
      use General, only: random_number_wrapper
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, tmp, r, p
      integer i,k
!
      if (initial_uu_noise /= 0) then
        if (lroot) print *, &
        'Disk section initial condition: setting uu (adding noise)'
        do n=1,mz
          ampl = initial_uu_noise/sqrt(3.0)*(cosh(z(n)/scaleheight))**(-2.0)
          do m=1,my
            do k=1,mx
              do i=0,2
                call random_number_wrapper(r)
                f(k,m,n,iuu+i)=f(k,m,n,iuu+i)+ampl*(2.0*r-1.0)
              enddo
            enddo
          enddo
        enddo
      endif
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  rho = A*sech^2(z/h) profile -- good for modelling an isothermal disk
!  section where A is computed from the parameters
!
      use Mpicomm,         only: stop_it
      use SharedVariables,  only: put_shared_variable
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: prof, gprof,del2_prof,del6_prof
      real, dimension(mx,my,mz,4) :: imposed_density
      real :: A
      rho_imp =1;grho_imp=2;d2rho_imp=3;d6rho_imp=4

      if (lroot) print *, 'Disk section initial condition: setting lnrho'

      ! Normalization of the profile
      A = pc_Sigma_g/(2.0*scaleheight)

      do n=n1,n2
        !      A * sech2(z/h)
        prof = A*(cosh(z(n)/scaleheight))**(-2.0)
        gprof= -2*(A/scaleheight)*(cosh(z(n)/scaleheight))**(-2.0)* tanh(z(n)/scaleheight)
        del2_prof = 2*(A/scaleheight**2)*((cosh(z(n)/scaleheight))**(-4.0)-&
                            2*(tanh(z(n)/scaleheight)**(2))*((cosh(z(n)/scaleheight))**(-2.0)))
        del6_prof = (-A/scaleheight**(-6.))*(17*(cosh(z(n)/scaleheight))**(-8.0) -&
                     180*(cosh(z(n)/scaleheight))**(-6.0)*tanh(z(n)/scaleheight)**(2) + &
                     114*tanh(z(n)/scaleheight)**(4)*(cosh(z(n)/scaleheight))**(-4.0) -&
                     4 *tanh(z(n)/scaleheight)**(6)*(cosh(z(n)/scaleheight))**(-2.0))
        do m=m1,m2
            f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho) +log(prof)
            imposed_density(l1:l2,m,n,rho_imp) = prof
            imposed_density(l1:l2,m,n,grho_imp) =gprof
            imposed_density(l1:l2,m,n,d2rho_imp)= del2_prof
            imposed_density(l1:l2,m,n,d6rho_imp) = del6_prof
        end do
      end do
                
      call put_shared_variable('imposed_density',imposed_density)

    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!LFSR comment: Is this used at all? I don't think so...
!
!  Initialize entropy calling the subroutine that sets it to isothermal.
!  adapted from
!
      use EquationOfState,only: gamma_m1, gamma1, get_cp1,cs20,get_soundspeed
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: lnrho0 = 0.
      real :: mu = 1! 0.63
      real :: TTi_over_TT0, TT0
      real :: cp1
      real :: cs2target

      if (lroot) print *, 'Disk section initial condition: setting ss'

      call get_cp1(cp1)

      TT0 = cs20*cp1/gamma_m1

      TTi_over_TT0 = pi*G_pencil_units*pc_Sigma_t*scaleheight/gamma1/cs20  &
                         /(1.0+Pmag_to_Ptherm+Pcr_to_Ptherm)

      TTi = TT0 * TTi_over_TT0
      call get_soundspeed(TTi,cs2target)

      if (lroot) print *, &
      'Disk section initial condition: Temperature set to ',TTi
      if (lroot) print *, &
     'Disk section initial condition: cs2bot and cs2top should be set to ',&
               cs2target

      do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) +      &
             (1.0/cp1)*gamma1 * (                    &
             - gamma_m1*(f(l1:l2,m,n,ilnrho)-lnrho0) &
             + log(TTi_over_TT0)         )

        enddo
      enddo
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential with a y-direction field
!  that allows hydrostatic equilibrium.
!
      use Mpicomm,         only: stop_it
      use SharedVariables, only: put_shared_variable
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: prof,prof1,del2_prof
      real, dimension (nx) :: aij_prof_x,aij_prof_y,aij_prof_z
      real, dimension (nx) :: bij_prof_x,bij_prof_y,bij_prof_z
      real :: A
      real,dimension(mx,mz,3,9) :: imposed_b_field
 !     real,dimension(mx,my,mz,3) :: A_imposed=0.


      aa_imp =1; bb_imp =2; del2_imp =3
      aij_impx=4; aij_impy=5; aij_impz=6
      bij_impx=7; bij_impy=8; bij_impz =9
      if (Pmag_to_Ptherm /= 0.0) then
        if (lroot) print *, 'Disk section initial condition: setting aa'
        A = sqrt(pi * G_pencil_units * pc_Sigma_t * pc_Sigma_g &
                  * Pmag_to_Ptherm/(1.0+Pmag_to_Ptherm+Pcr_to_Ptherm))
        do n=n1,n2
          prof = A*2.0*scaleheight*atan(tanh(z(n)/scaleheight/2.0))
            
          del2_prof = A*(tanh(z(n)/scaleheight/2)*cosh(z(n)/scaleheight/2)*&
                       (tanh(z(n)/scaleheight/2)**2 + cosh(z(n)/scaleheight/2)**2 +1))/&
                      (2*scaleheight**2*(tanh(z(n)/scaleheight/2)+1)**2)
          aij_prof_z = A*cosh(z(n)/scaleheight/2)**(-2)/((2*scaleheight)*tanh(z(n)/scaleheight/2) +2*scaleheight)
          bij_prof_z = 2*A*sinh(z(n)/scaleheight)
          prof1 = A*2.0*scaleheight*atan(tanh(z(n)/scaleheight/2.0))
          do m=m1,m2
            f(l1:l2,m,n,iaa)=f(l1:l2,m,n,iaa)+ prof
            imposed_b_field(l1:l2,n,1,aa_imp) = prof1
            imposed_b_field(l1:l2,n,2,bb_imp) = prof
            imposed_b_field(l1:l2,n,1,del2_imp) = del2_prof
            imposed_b_field(l1:l2,n,1,aij_impx) = 0!aij_prof_x
            imposed_b_field(l1:l2,n,2,aij_impx) = 0!aij_prof_y
            imposed_b_field(l1:l2,n,3,aij_impx) = aij_prof_z
            imposed_b_field(l1:l2,n,:,aij_impy) = 0
            imposed_b_field(l1:l2,n,:,aij_impz) = 0
            imposed_b_field(l1:l2,n,1,bij_impx) = 0
            imposed_b_field(l1:l2,n,1,bij_impy) = 0!bij_prof_x
            imposed_b_field(l1:l2,n,2,bij_impy) = 0!bij_prof_y
            imposed_b_field(l1:l2,n,3,bij_impy) = bij_prof_z
            imposed_b_field(l1:l2,n,:,bij_impz) = 0
          end do
        end do
      endif
      call put_shared_variable('imposed_b_field',imposed_b_field)
!
    endsubroutine initial_condition_aa

!***********************************************************************
    subroutine initial_condition_aatest(f)
!
!  Initialize testfield.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aatest
!***********************************************************************
    subroutine initial_condition_uutest(f)
!
!  Initialize testflow.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uutest
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize passive scalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine initial_condition_chiral(f)
!
!  Initialize chiral.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chiral
!***********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chemistry
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!***********************************************************************
    subroutine initial_condition_lnrhon(f)
!
!  Initialize neutral fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrhon
!***********************************************************************
    subroutine initial_condition_ecr(f)
!
!  Initialize the cosmic rays energy density, imposing hydrostatic
!  equilibrium.
!
      use Mpicomm,         only: stop_it

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: prof
      real :: A ! Normalization of the profile
!         pcr=(gammacr-1)*ecr
      if (Pcr_to_Ptherm /= 0.0) then
        if (lroot) print *, 'Disk section initial condition: setting ecr'
        if (lroot) print *, &
        'Disk section initial condition: gammacr should be set to ', gammacr

        A = pi * G_pencil_units* pc_Sigma_t*pc_Sigma_g/2.0 &
             /(gammacr-1.0)*Pcr_to_Ptherm/(1.0+Pmag_to_Ptherm+Pcr_to_Ptherm)
        if (lroot) print *, &
             'Disk section initial condition: Qcr3=',A
        do n=n1,n2
          !      A * sech2(z/h)
          prof = A*(cosh(z(n)/scaleheight))**(-2.0)
          do m=m1,m2
            f(l1:l2,m,n,iecr)=f(l1:l2,m,n,iecr)+prof
          end do
        end do
      endif
!
    endsubroutine initial_condition_ecr
!***********************************************************************
    subroutine initial_condition_fcr(f)
!
!  Initialize cosmic ray flux.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_fcr
!***********************************************************************
    subroutine initial_condition_solid_cells(f)
!
!  Initialize solid cells.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_solid_cells
!***********************************************************************
    subroutine initial_condition_cctest(f)
!
!  Initialize testscalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_cctest
!***********************************************************************
    subroutine initial_condition_xxp(f,fp)
!
!  Initialize particles' positions.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f,fp)
!
!  Initialize particles' velocities.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_vvp
!***********************************************************************
  subroutine read_initial_condition_pars(iostat)

      use File_io, only: parallel_unit

      integer, intent(out) :: iostat
      iostat = 0

      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)

    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)

      integer, intent(in) :: unit

      write(unit, NML=initial_condition_pars)

    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine initial_condition_clean_up
!
!  04-may-11/dhruba: coded
! dummy
!
    endsubroutine initial_condition_clean_up
!***********************************************************************
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
!     include 'initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
