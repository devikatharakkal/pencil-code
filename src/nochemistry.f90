! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Chemistry

  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'chemistry.h'

!!  character, len(50) :: initcustom

! input parameters
!  namelist /special_init_pars/ dummy
!!!eg.    initcustom
! run parameters
!  namelist /special_run_pars/ dummy

!!
!! Declare any index variables necessary for main or
!!
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!!
!! other variables (needs to be consistent with reset list below)
!!
!!   integer :: i_POSSIBLEDIAGNOSTIC=0
!!

  contains

!***********************************************************************
    subroutine register_chemistry()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use Cdata
!
!  Identify CVS/SVN version information.
!
      if (lroot) call cvs_id( &
          "$Id$")
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine initialize_chemistry(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
!!
!!  Initialize any module variables which are parameter dependent
!!
!
! DO NOTHING
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chemistry
!***********************************************************************
    subroutine init_chemistry(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f

!!
!!  SAMPLE IMPLEMENTATION
!!
!!      select case(initspecial)
!!        case('nothing'); if (lroot) print*,'init_special: nothing'
!!        case('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!!        case default
!!          !
!!          !  Catch unknown values
!!          !
!!          if (lroot) print*,'init_special: No such value for initspecial: ', trim(initspecial)
!!          call stop_it("")
!!      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_chemistry
!***********************************************************************
    subroutine pencil_criteria_chemistry()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_chemistry
!***********************************************************************
    subroutine pencil_interdep_chemistry(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      use Sub, only: keep_compiler_quiet
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_chemistry
!***********************************************************************
 subroutine calc_for_chem_mixture(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
 endsubroutine calc_for_chem_mixture
!***********************************************************************
    subroutine dchemistry_dt(f,df,p)
!
!  Dummy routine
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
! Keep compiler quiet by ensuring every parameter is used
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dchemistry_dt
!***********************************************************************
    subroutine read_chemistry_init_pars(unit,iostat)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_chemistry_init_pars
!***********************************************************************
    subroutine write_chemistry_init_pars(unit)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_chemistry_init_pars
!***********************************************************************
    subroutine read_chemistry_run_pars(unit,iostat)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_chemistry_run_pars
!***********************************************************************
    subroutine write_chemistry_run_pars(unit)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_chemistry_run_pars
!***********************************************************************
    subroutine rprint_chemistry(lreset,lwrite)
!
!  reads and registers print parameters relevant to chemistry
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub
!!
!!!   SAMPLE IMPLEMENTATION
!!
!!      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
!!        i_SPECIAL_DIAGNOSTIC=0
      endif
!!
!!      do iname=1,nname
!!        call parse_name(iname,cname(iname),cform(iname),'NAMEOFSPECIALDIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
!!      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        write(3,*) 'i_SPECIAL_DIAGNOSTIC=',i_SPECIAL_DIAGNOSTIC
!!      endif
!!

    endsubroutine rprint_chemistry
!***********************************************************************
    subroutine get_slices_chemistry(f,slices)
!
!  Write slices for animation of chemistry variables.
!
!  26-jun-06/tony: dummy
!
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_chemistry
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine chemistry_boundconds(f,bc)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine chemistry_boundconds
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine calc_cs2x(cs2x,topbot,f)
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (ny,nz) :: cs2x
      character (len=3) :: topbot
      integer :: k
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_cs2x
!********************************************************************
 subroutine calc_cs2y(cs2y,topbot,f)
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (nx,nz) :: cs2y
      character (len=3) :: topbot
      integer :: k
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_cs2y
!********************************************************************
    subroutine get_gammax(gamma0,topbot)
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (ny,nz) :: gamma0
      character (len=3) :: topbot
      integer :: k
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_gammax
!********************************************************************
   subroutine get_gammay(gamma0,topbot)
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (nx,nz) :: gamma0
      character (len=3) :: topbot
      integer :: k
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_gammay
!********************************************************************
    subroutine get_p_infx(p_infx,topbot)
!
      use Mpicomm
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (ny,nz) :: p_infx
      character (len=3) :: topbot
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_p_infx
!********************************************************************
   subroutine get_p_infy(p_infy,topbot)
!
      use Mpicomm
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (nx,nz) :: p_infy
      character (len=3) :: topbot
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_p_infy
!*************************************************************
   subroutine get_rhs_Y(topbot,j,bound_reac_term)

    use Mpicomm
!
!  Dummy routine
!
     real, dimension (ny,nz,nchemspec) :: bound_reac_term
     integer :: j
     character (len=3) :: topbot

     intent(out) :: bound_reac_term
     intent(in) :: j
!
!  set bound_reac_term to keep compiler quiet
!
      bound_reac_term=0.
!
   endsubroutine get_rhs_Y
!*************************************************************
  subroutine get_rhs_T(topbot,j,bound_rhs_T)
 
    use Mpicomm
!
!  Dummy routine
!
     real, dimension (ny,nz) :: bound_rhs_T
     integer :: j
     character (len=3) :: topbot

     intent(out) :: bound_rhs_T
     intent(in) :: j
!
!  set bound_reac_term to keep compiler quiet
!
      bound_rhs_T=0.
!
   endsubroutine get_rhs_T
!*************************************************************
  subroutine get_gamma_full(gamma_full)
!
  use Mpicomm
!
   real, dimension (mx,my,mz) :: gamma_full
!
   intent(out) :: gamma_full
!   
   gamma_full=0.
   call keep_compiler_quiet(gamma_full)
!
  endsubroutine get_gamma_full
!*************************************************************
 subroutine get_cs2_full(cs2_full)

  use Mpicomm

   real, dimension (mx,my,mz) :: cs2_full

   intent(out) :: cs2_full
!
   cs2_full=0.
   call keep_compiler_quiet(cs2_full)
!
  endsubroutine get_cs2_full
!*************************************************************
subroutine chemistry_clean_up()

  use Mpicomm
  real, dimension(mx,my,mz,mfarray) :: f

     call keep_compiler_quiet(f)
!
  endsubroutine chemistry_clean_up
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'special_dummies.inc'
!********************************************************************

endmodule Chemistry

