#include "precompilerdefinitions"
module type_gridsim
use konstanter, only: r8,flyt,lo_huge,lo_hugeint,lo_exitcode_symmetry,lo_exitcode_blaslapack,lo_status,lo_exitcode_param,&
                      lo_exitcode_io,lo_tol,lo_sqtol,lo_volume_A_to_bohr,lo_volume_bohr_to_A,lo_Hartree_to_eV,lo_kb_Hartree,&
                      lo_pi
use gottochblandat, only: walltime,tochar,lo_progressbar_init,lo_progressbar,open_file,&
                   lo_clean_fractional_coordinates,lo_determ,lo_unflatten_2tensor,lo_flattentensor,lo_chop,&
                   lo_looptimer,lo_mean,lo_trapezoid_integration,lo_identitymatrix,lo_linspace,&
                   lo_enforce_linear_constraints,lo_linear_least_squares
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_SUM,MPI_MAX
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_firstorder, only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_jij_secondorder, only: lo_jij_secondorder
use type_forcemap !, only: lo_forcemap
use lo_dielectric_interaction, only: lo_dielectric_tensor
use type_mdsim, only: lo_mdsim
use type_blas_lapack_wrappers, only: lo_gemm,lo_dgels,lo_dgglse,lo_dgesvd,lo_dgelss
use type_equation_of_state, only: lo_eos,lo_eos_1d,lo_eos_2d,lo_eos_birch_murnaghan,lo_eos_vinet,lo_eos_2d_birch_murnaghan
use type_polynomial_interpolation, only: lo_polynomial,lo_grid_interpolation
use ifc_solvers, only: lo_solve_for_borncharges
use helperobjects, only: lo_sparsematrix,reduce_equations
implicit none

private
public :: lo_gridsim

!> keep all the setting in one place to make it less messy, or something
type lo_gridsim_info
    !> settings
    logical :: magnetic_pair_interactions=.false.
    !> what kind of fit are we doing?
    integer :: pairfittype=-lo_hugeint
    !> some stuff that decides how I do the interpolation
    character(len=100), dimension(:), allocatable :: dimension_names
    integer, dimension(:), allocatable :: order_per_dim
    integer :: dim_volume=-lo_hugeint
    integer :: dim_temperature=-lo_hugeint
    integer :: dim_eta=-lo_hugeint
    real(r8) :: temperature_scale=lo_huge
    real(r8) :: distance_scale=lo_huge
    logical :: weighted=.false.
end type

!> simulation data
type lo_gridsim_rawdata
    !> how many different configurations is it?
    integer :: nconf
    !> displacements
    real(r8), dimension(:,:,:), allocatable :: u
    !> forces
    real(r8), dimension(:,:,:), allocatable :: f
    real(r8), dimension(:,:,:), allocatable :: f0
    !> magnetic moments
    real(r8), dimension(:,:,:), allocatable :: m
    !> dielectric tensors
    real(r8), dimension(:,:,:), allocatable :: eps
    !> Born charges
    real(r8), dimension(:,:,:,:), allocatable :: Z
    !> potential energies
    real(r8), dimension(:), allocatable :: e
    real(r8), dimension(:), allocatable :: e_polar
    real(r8), dimension(:), allocatable :: e_fc_pair
    real(r8), dimension(:), allocatable :: e_fc_triplet
    real(r8), dimension(:), allocatable :: e_fc_quartet
    !real(r8), dimension(:), allocatable :: e_magnetic
    !real(r8), dimension(:), allocatable :: e_magnetic_jij
    !real(r8), dimension(:), allocatable :: e_magnetic_tij
    !real(r8), dimension(:), allocatable :: e_magnetic_qij
    !> index on the grid?
    integer, dimension(:), allocatable :: gridind
    !> how many of the total grid-points are relevant?
    integer :: nrelevant_gridpoints
    !> relevant points
    integer, dimension(:), allocatable :: relevant_gridpoints
end type

!> reference data, the pretty small things
type lo_gridsim_refdata
    !> latticevectors
    real(r8), dimension(3,3) :: unitcell_latticevectors
    real(r8), dimension(3,3) :: supercell_latticevectors
    !> reference positions
    real(r8), dimension(:,:), allocatable :: unitcell_positions
    real(r8), dimension(:,:), allocatable :: supercell_positions
    !> reference atomic numbers (maybe later)
    integer, dimension(:), allocatable :: unitcell_atomic_numbers
    integer, dimension(:), allocatable :: supercell_atomic_numbers
    !> dielectric tensor and Born effective charges
    real(r8), dimension(3,3) :: dielectric_tensor
    real(r8), dimension(:,:,:), allocatable :: born_effective_charges
    !real(r8), dimension(:,:,:,:), allocatable :: dipole_forceconstant
    !> Constraints that are different per structure
    integer :: nconstr_pair=-lo_hugeint
    real(r8), dimension(:,:), allocatable :: constr_pair
    !> average forces, energy and magnetic moments?

    !> weight per step for this gridpoint?
    real(r8) :: weight=0.0_r8
end type

!> template that holds polynomial coefficients
type lo_gridsim_coeff
    !> how many variables should it be?
    integer :: nvar=-lo_hugeint
    !> polynomial coefficients
    real(r8), dimension(:,:), allocatable :: coeff
    !> number of constraints on the resulting variables
    integer :: nconstr=-lo_hugeint
    !> constraints
    real(r8), dimension(:,:), allocatable :: constraints
    !> R^2 from the fit
    real(r8), dimension(:), allocatable :: rsquare
    !> polynomical for interpolation
    type(lo_polynomial) :: pl
end type

!> Template that holds the structure
type lo_gridsim_structure
    !> number of internal degrees of freedom
    integer :: ninternal
    !> atomic numbers, for reference
    integer, dimension(:), allocatable :: atomic_numbers
    !> reference lattice vectors
    real(r8), dimension(3,3) :: lv0,ilv0
    !> reference internal positions
    real(r8), dimension(:,:), allocatable :: r0
    !> transformation matrix from irreducible to full
    real(r8), dimension(:,:), allocatable :: coeffM
    !> work-arrays to transform back from irreducible representation
    real(r8), dimension(:), allocatable :: wV,wU
    real(r8), dimension(:,:), allocatable :: wM
    !> interpolation of the internal degrees of freedom
    type(lo_grid_interpolation) :: ipint
    type(lo_grid_interpolation) :: iplv
    type(lo_grid_interpolation) :: ipvol
    contains
        procedure :: interpolate=>interpolate_structure_from_depvar
end type

!> Template that holds delta-energies
type lo_gridsim_energy
    !> xn dependent variables
    real(r8), dimension(:,:), allocatable :: xi,xi_nontr
    !> baseline xi, scaling factor
    real(r8), dimension(:), allocatable :: x0,x0_nontr,ixs,ixs_nontr
    !> characteristic separation
    real(r8) :: xdist
    !> energy to interpolate
    real(r8), dimension(:), allocatable :: energy
    !> polynomical for interpolation
    type(lo_polynomial) :: pl
    !> temperature dimension
    integer :: dim_temperature
    !> temperature scale
    real(r8) :: temperature_scale
    contains
        procedure :: interpolate=>interpolate_energy_from_depvar
end type

!> Template that holds the polar information
type lo_gridsim_polar
    !> how many independent Z
    integer :: nZ=-lo_hugeint
    !> how many indepenedent eps?
    integer :: neps=-lo_hugeint
    !> interpolation for Z
    type(lo_grid_interpolation) :: ipZ
    type(lo_grid_interpolation) :: ipeps
    !> R^2 in forces
    real(r8), dimension(:), allocatable :: rsquare
#ifdef AGRESSIVE_SANITY
    contains
        !> storage size
        procedure :: size_in_mem=>memsize_polar
#endif
end type

!> hold all the pair information
type lo_gridsim_pair
    !> how many forceconstants
    integer :: nfc
    !> the interpolation
    type(lo_grid_interpolation) :: ip
    !> R^2 in forces
    real(r8), dimension(:), allocatable :: rsquare
#ifdef AGRESSIVE_SANITY
    contains
        !> storage size
        procedure :: size_in_mem=>memsize_pair
#endif
end type

!> hold all grid interpolation information
type lo_gridsim_gf
    !> how many forceconstants
    integer :: nfc_pair
    integer :: nfc_triplet
    integer :: nfc_quartet
    !> the interpolation
    type(lo_grid_interpolation) :: ip2
    type(lo_grid_interpolation) :: ip3
    type(lo_grid_interpolation) :: ip4
end type

!> hold all the pair information
type lo_gridsim_magpair
    !> how many forceconstants
    integer :: nfc,nfc_T
    !> the interpolation
    type(lo_grid_interpolation) :: ipM0
    type(lo_grid_interpolation) :: ipMD
    type(lo_grid_interpolation) :: ipqij
    type(lo_grid_interpolation) :: ipJ
    type(lo_grid_interpolation) :: ipT
    !> R^2 in forces
    real(r8), dimension(:), allocatable :: rsquare_force
    real(r8), dimension(:), allocatable :: rsquare_energy_jij
    real(r8), dimension(:), allocatable :: rsquare_energy_tij
    real(r8), dimension(:), allocatable :: rsquare_energy_qij
    !> coefficient matrix
    real(r8), dimension(:,:), allocatable :: coeff
    real(r8), dimension(:,:), allocatable :: coeff_T
end type

!> grid of simulations
type lo_gridsim
    !> how many atoms in the supercell?
    integer :: na_ss=-lo_hugeint
    !> how many atoms in the unitcell
    integer :: na_uc=-lo_hugeint
    !> how many simulations, in total
    integer :: nsim=-lo_hugeint
    !> number of dimensions of the grid
    integer :: ndim=-lo_hugeint
    !> x,y,z-values for the simulation grid.
    real(r8), dimension(:,:), allocatable :: grid_coordinates
    !> polynomial coefficients
    type(lo_gridsim_coeff) :: pair
    type(lo_gridsim_coeff) :: triplet
    type(lo_gridsim_coeff) :: quartet
    type(lo_gridsim_coeff) :: eps_global
    type(lo_gridsim_coeff) :: eps_singlet
    type(lo_gridsim_coeff) :: eps_pair
    type(lo_gridsim_coeff) :: Z_singlet
    type(lo_gridsim_coeff) :: Z_pair
    type(lo_gridsim_coeff) :: Z_triplet
    !> settings
    type(lo_gridsim_info) :: info
    !> force-displacement-magnetic-energy data
    type(lo_gridsim_rawdata) :: raw
    !> reference data for every gridpoint
    type(lo_gridsim_refdata), dimension(:), allocatable :: ref
    !> data used for structure interpolation
    type(lo_gridsim_structure) :: structure
    !> data used for energy interpolation
    type(lo_gridsim_energy) :: energy
    !> data used for polar interpolation
    type(lo_gridsim_polar) :: polar
    !> grid-fitted forceconstants
    type(lo_gridsim_gf) :: gf
    !>
    type(lo_gridsim_magpair) :: magpair
    !> polynomial
    type(lo_polynomial) :: poly
    !> zero Kelvin equation of state, used for transformations and things
    class(lo_eos), allocatable :: eos
    !> new attempt at pair force constants
    type(lo_gridsim_pair) :: gridpair
    contains
        !> set up the basic things
        procedure :: init=>initialize_gridsim
        !> grab all simulation data from file
        procedure :: read_simulations_from_file
        !> get the coefficients for the structure
        procedure :: create_structure_interpolation
        !> get the coefficients for the energy
        procedure :: create_energy_interpolation
        !> get the coefficients for the polar part
        procedure :: create_polar_interpolation
        !> get the forceconstants
        procedure :: solve_secondorder
        procedure :: solve_secondorder_gridfit
        procedure :: solve_thirdorder
        procedure :: solve_fourthorder
        !procedure :: solve_magnetic_onsite
        !procedure :: solve_magnetic_pair_polyfit
        !procedure :: solve_magnetic_pair_gridfit
        procedure :: solve_gridfit
        procedure :: solve_dielectric
        !> remove forces
        procedure :: subtract_polar_forces
        procedure :: subtract_secondorder_forces
        procedure :: subtract_thirdorder_forces
        procedure :: subtract_fourthorder_forces
        !procedure :: subtract_magnetic_forces
        !> evaluate everything at a certain point
        procedure :: eval=>evaluate_irreducible
        !> transform coordinates
        procedure :: coordinate_transformation
        !> dump diagnostic data
        !procedure :: dump_diagnostic_data
#ifdef AGRESSIVE_SANITY
        procedure :: memreport
#endif
end type

!> some default parameters I have to fix later
integer, parameter :: maxdim=3
integer, parameter :: maxorder=4
integer, parameter :: maxcoeff=640
!> how often (in seconds) to report progress to stdout for the time-consuming parts
real(r8), parameter :: timereport=15.0_r8
!> max size of arrays, in MiB, that I will try to allocate.
real(r8), parameter :: maxmem=250.0_r8

contains

#include "type_gridsim_solve_grfit.f90"
#include "type_gridsim_solve_polar.f90"
#include "type_gridsim_solve_secondorder.f90"
#include "type_gridsim_solve_thirdorder.f90"
#include "type_gridsim_solve_fourthorder.f90"
#include "type_gridsim_polynomial.f90"
#include "type_gridsim_structure.f90"
#include "type_gridsim_energy.f90"
#include "type_gridsim_io.f90"

!> initialize things and set all heuristics
subroutine initialize_gridsim(gs,filename,map,order,temperaturescale,distancescale,pairfittype,weighted,verbosity,mw)
    !> simulation grid
    class(lo_gridsim), intent(out) :: gs
    !> filename with settings
    character(len=*), intent(in) :: filename
    !> forcemap to steal settings from
    type(lo_forcemap), intent(in) :: map
    !> order of polynomial fit
    integer, intent(in) :: order
    !> quantum scale temperatures in a clever way
    real(r8), intent(in) :: temperaturescale
    !> scale distances in a clever way
    real(r8), intent(in) :: distancescale
    !> what type of fit for the second order
    integer, intent(in) :: pairfittype
    !> weight timesteps
    logical, intent(in) :: weighted
    !> how much to talk
    integer, intent(in) :: verbosity
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    if ( verbosity .gt. 0 ) then
        write(*,*) ''
        write(*,*) 'SETTING UP SOLUTION FOR A GRID OF SIMULATIONS'
    endif

    ! First grab the polynomial coefficients and set that up. Future Olle will
    ! add the option of using different order polynomical for different things.
    checkfile: block
        real(r8) :: f1,f2,f3,f4,f5,f6,f7,f8,f9
        character(len=1000) :: dum,eosname
        integer :: u,i

        ! briefly open the input files and grab some data
        u=open_file('in',trim(filename))
            read(u,*) gs%ndim                  ! number of dimensions
            read(u,*) gs%nsim                  ! number of simulations
            lo_allocate(gs%info%dimension_names(gs%ndim))
            lo_allocate(gs%info%order_per_dim(gs%ndim))
            read(u,*) gs%info%dimension_names  ! what are the dimensions called
            read(u,*) gs%info%order_per_dim    ! polynomial order per dimension
            read(u,*) eosname
            ! figure out which equation of state I want to use.
            if ( eosname(1:5) .eq. 'Birch' .or. eosname(1:5) .eq. 'birch' ) then
                ! Birch-murnaghan equation of state! Read four parameters
                read(u,*) f1,f2,f3,f4
                ! generate the equation of state.
                allocate(lo_eos_birch_murnaghan::gs%eos)
                select type(eos=>gs%eos); type is(lo_eos_birch_murnaghan)
                    call eos%generate( f1,f2,f3,f4 )
                end select
                if ( verbosity .gt. 0 ) write(*,*) '... Got Birch-Murnaghan equation of state.'
            elseif ( eosname(1:5) .eq. 'Vinet' .or. eosname(1:5) .eq. 'vinet' ) then
                ! Vinet equation of state! Read four parameters
                read(u,*) f1,f2,f3,f4
                ! generate the equation of state.
                allocate(lo_eos_vinet::gs%eos)
                select type(eos=>gs%eos); type is(lo_eos_vinet)
                    call eos%generate( f1,f2,f3,f4 )
                end select
                if ( verbosity .gt. 0 ) write(*,*) '... Got Vinet equation of state.'
            elseif ( eosname(1:5) .eq. '2D-Bi' ) then
                ! 2D Birch-Murnaghan EOS
                read(u,*) f1,f2,f3,f4,f5,f6,f7,f8,f9
                ! generate the equation of state.
                allocate(lo_eos_2d_birch_murnaghan::gs%eos)
                select type(eos=>gs%eos); type is(lo_eos_2d_birch_murnaghan)
                    call eos%generate_2d( f1,f2,f3,f4,f5,f6,f7,f8,f9 )
                end select
                if ( verbosity .gt. 0 ) write(*,*) '... Got 2D Birch-Murnaghan equation of state.'
            elseif ( eosname(1:4).eq. 'none' .or. eosname(1:4) .eq. 'null' ) then
                ! I'll make it a null type. Should be reasonable. Or not.
                ! No equation of state at all. Not sure what to do actually.
                !call lo_stop_gracefully(['No equation of state provided. Have not decided how to handle that.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            else
                call lo_stop_gracefully(['No equation of state provided. Try "null" if you really do not want one.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
                ! we have no equation of state, specify that it's the null type.
            endif

            ! Ok got the dimensions of the grid, and perhaps an equation of state. Now read in the grid coordinates.
            lo_allocate(gs%grid_coordinates(gs%ndim,gs%nsim))
            gs%grid_coordinates=0.0_r8
            do i=1,gs%nsim
                read(u,*) gs%grid_coordinates(:,i),dum
            enddo
        close(u)

        ! Figure out if we have temperature,volume or eta.
        gs%info%dim_temperature=-1
        gs%info%dim_volume=-1
        gs%info%dim_eta=-1
        do i=1,gs%ndim
            if ( trim( gs%info%dimension_names(i) ) .eq. 'T' )   gs%info%dim_temperature=i
            if ( trim( gs%info%dimension_names(i) ) .eq. 'V' )   gs%info%dim_volume=i
            if ( trim( gs%info%dimension_names(i) ) .eq. 'eta' ) gs%info%dim_eta=i
        enddo

        ! Convert volumes to atomic units
        if ( gs%info%dim_volume .gt. 0 ) gs%grid_coordinates(gs%info%dim_volume,:)=gs%grid_coordinates(gs%info%dim_volume,:)*lo_volume_A_to_bohr

        ! Also stop now if I don't have dimensions that make sense.
        select type(eos=>gs%eos)
        class is(lo_eos_1d)
            if ( gs%info%dim_volume .eq. -1 ) then
                call lo_stop_gracefully(['Equation of state provided, but no dimension is labelled volume (indicated by "V")'],&
                                        lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        class is(lo_eos_2d)
            if ( gs%info%dim_volume .eq. -1 ) then
                call lo_stop_gracefully(['Equation of state provided, but no dimension is labelled volume (indicated by "V")'],&
                                        lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
            if ( gs%info%dim_eta .eq. -1 ) then
                call lo_stop_gracefully(['Equation of state provided, but no dimension is labelled eta (indicated by "eta")'],&
                                        lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        end select

        ! Scale the temperatures cleverly. But to do that one of the dimensions has to be temperature.
        if ( temperaturescale .gt. 0.0_r8 ) then
            gs%info%temperature_scale=temperaturescale
            if ( gs%info%dim_temperature .eq. -1 ) then
                call lo_stop_gracefully(['Temperature scale provided, but no dimension is labelled temperature (indicated by "T")'],&
                                        lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif
        else
            gs%info%temperature_scale=-1.0_r8
        endif
        ! and some information regarding the second order fits
        gs%info%distance_scale=distancescale
        gs%info%pairfittype=pairfittype
        gs%info%weighted=weighted
    end block checkfile

    ! get the skeleton for the polynomials
    buildpoly: block
        real(r8), dimension(:,:), allocatable :: transformed_coordinates
        integer :: i

        ! transform the coordinates
        lo_allocate(transformed_coordinates(gs%ndim,gs%nsim))
        do i=1,gs%nsim
            call coordinate_transformation(gs,gs%grid_coordinates(:,i),transformed_coordinates(:,i))
        enddo

        ! Now particularly meaningful block at the moment, but I plan to add different polynomials for the different orders soon.
        call gs%poly%init( order,gs%ndim,transformed_coordinates,gs%info%dimension_names, gs%info%order_per_dim )
        lo_deallocate(transformed_coordinates)
    end block buildpoly

    ! Set what do do, and how many variables are needed to do so.
    heur: block
        integer :: i,j,k,l
        character(len=2000) :: ds

        ! Number of atoms?
        gs%na_ss=map%n_atom_ss
        gs%na_uc=map%n_atom_uc
        ! l=gs%info%nphi_singlet +&
        !   gs%info%nphi_pair    +&
        !   gs%info%nphi_triplet +&
        !   gs%info%nphi_quartet +&
        !   gs%info%nphi_eps     +&
        !   gs%info%nphi_Z
        l=0
        if ( verbosity .gt. 0 ) then
            write(*,*) '          total number of coefficients: ',tochar(l)
            write(*,*) '     using a ',tochar(gs%poly%ndim),'-D polynomial of order ',tochar(gs%poly%order),': '
            do i=1,gs%poly%ncoeff
                ds=''
                do j=1,gs%poly%ndim
                    k=gs%poly%exponents(j,i)
                    ds=trim(ds)//trim(gs%poly%variablename(j))//'^'//tochar(k)
                enddo
                write(*,*) '               coeff ',tochar(i,-2), ': ',trim(ds)
            enddo
        endif
    end block heur
end subroutine

!> evaluate the irreducible representation at a certain point
subroutine evaluate_irreducible(gs,map,gridcoord,pairconstr)
    !> simulation grid
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> point on grid
    real(r8), dimension(:), intent(in) :: gridcoord
    !> pair constraints
    real(r8), dimension(:,:), intent(in), optional :: pairconstr

    real(r8), dimension(gs%ndim) :: transformed_gridcoord
    real(r8) :: f0,f1
    integer :: i

    ! Apply my transformation
    call coordinate_transformation(gs,gridcoord,transformed_gridcoord)

    select case(map%polar)
    case(1)
        ! Just the Born charges and dielectric constant
        do i=1,gs%polar%neps
            map%xuc%x_eps_global(i)=gs%polar%ipeps%eval( i,gridcoord )
        enddo
        do i=1,gs%polar%nZ
            map%xuc%x_Z_singlet(i)=gs%polar%ipZ%eval( i,gridcoord )
        enddo
    case(2)
        ! All the fancy dielectric stuff.
        do i=1,map%xuc%nx_eps_global
            f0=gs%polar%ipeps%eval( i,gridcoord )
            f1=gs%poly%eval(transformed_gridcoord,gs%eps_global%coeff(:,i))
            map%xuc%x_eps_global(i)=f0+f1 !gs%poly%eval(transformed_gridcoord,gs%eps_global%coeff(:,i))
        enddo
        do i=1,map%xuc%nx_eps_singlet
            map%xuc%x_eps_singlet(i)=gs%poly%eval(transformed_gridcoord,gs%eps_singlet%coeff(:,i))
        enddo
        do i=1,map%xuc%nx_eps_pair
            map%xuc%x_eps_pair(i)=gs%poly%eval(transformed_gridcoord,gs%eps_pair%coeff(:,i))
        enddo

        do i=1,map%xuc%nx_Z_singlet
            map%xuc%x_Z_singlet(i)=gs%poly%eval(transformed_gridcoord,gs%Z_singlet%coeff(:,i))
        enddo
        do i=1,map%xuc%nx_Z_pair
            map%xuc%x_Z_pair(i)=gs%poly%eval(transformed_gridcoord,gs%Z_pair%coeff(:,i))
        enddo
        if ( map%have_Z_triplet ) then
            do i=1,map%xuc%nx_Z_triplet
                map%xuc%x_Z_triplet(i)=gs%poly%eval(transformed_gridcoord,gs%Z_triplet%coeff(:,i))
            enddo
            if ( map%constraints%neqz3 .gt. 0 ) then
                call lo_enforce_linear_constraints(map%constraints%eqz3,map%xuc%x_Z_triplet)
            endif
        endif

    end select

    if ( gs%info%magnetic_pair_interactions ) then
        write(*,*) 'FIXME MAGNETIC PAIR',__FILE__
        ! Local Polynomial. Almost always bad.
        !do i=1,map%ntheta_magpair
        !    map%theta_magpair(i)=gs%magpair%ipJ%eval(i, transformed_gridcoord, rexp=1)
        !enddo
        ! Local polynomial for crossterm
        !do i=1,map%ntheta_magpair_T
        !    map%theta_magpair_T(i)=gs%magpair%ipT%eval(i, transformed_gridcoord, rexp=1)
        !enddo
        ! do i=1,map%n_atom_uc
        !     j=map%uc(i)%J%irreducible_shell
        !     map%theta_magsinglet_mi(j)=gs%magpair%ipM0%eval(j, transformed_gridcoord, rexp=1)
        !     map%theta_magsinglet_si(j)=gs%magpair%ipMD%eval(j, transformed_gridcoord, rexp=1)
        !     map%uc(i)%J%m0=map%theta_magsinglet_mi(j)
        !     map%uc(i)%J%m0dev=map%theta_magsinglet_si(j)
        ! enddo
        ! do j=1,map%ntheta_magpair_qij
        !     map%theta_magpair_qij(j)=gs%magpair%ipqij%eval(j, transformed_gridcoord, rexp=1)
        ! enddo
        ! ! Global polynomial
        ! do i=1,gs%magpair%nfc
        !     map%theta_magpair_jij(i)=gs%poly%eval( transformed_gridcoord,gs%magpair%coeff(:,i))
        ! enddo
        ! ! Global polynomial
        ! do i=1,gs%magpair%nfc_T
        !     map%theta_magpair_tij(i)=gs%poly%eval( transformed_gridcoord,gs%magpair%coeff_T(:,i))
        ! enddo
    endif

    select case(gs%info%pairfittype)
    case(1)
        ! Global polynomial
        if ( map%have_fc_pair ) then
            do i=1,gs%pair%nvar
                map%xuc%x_fc_pair(i)=gs%poly%eval(transformed_gridcoord,gs%pair%coeff(:,i))
            enddo
            if ( present(pairconstr) ) then
                !call enforce_constraints( map%xuc%x_fc_pair,pairconstr )
                call lo_enforce_linear_constraints(pairconstr,map%xuc%x_fc_pair)
            endif
        endif
        if ( map%have_fc_triplet ) then
            do i=1,gs%triplet%nvar
                map%xuc%x_fc_triplet(i)=gs%poly%eval(transformed_gridcoord,gs%triplet%coeff(:,i))
            enddo
            if ( gs%triplet%nconstr .gt. 0 ) then
                !call enforce_constraints( map%xuc%x_fc_triplet,gs%triplet%constraints )
                call lo_enforce_linear_constraints(gs%triplet%constraints,map%xuc%x_fc_triplet)
            endif
        endif
        if ( map%have_fc_quartet ) then
            do i=1,gs%quartet%nvar
                map%xuc%x_fc_quartet(i)=gs%poly%eval(transformed_gridcoord,gs%quartet%coeff(:,i))
            enddo
            if ( gs%quartet%nconstr .gt. 0 ) then
                !call enforce_constraints( map%xuc%x_fc_quartet,gs%quartet%constraints )
                call lo_enforce_linear_constraints(gs%quartet%constraints,map%xuc%x_fc_quartet)
            endif
        endif
    case(2)
        ! Local polynomial
        if ( map%have_fc_pair ) then
            do i=1,map%xuc%nx_fc_pair
                map%xuc%x_fc_pair(i)=gs%gf%ip2%eval(i, gridcoord, rexp=2)
            enddo
            if ( present(pairconstr) ) then
                !call enforce_constraints( map%xuc%x_fc_pair,pairconstr )
                call lo_enforce_linear_constraints(pairconstr,map%xuc%x_fc_pair)
            endif
        endif
        if ( map%have_fc_triplet ) then
            do i=1,map%xuc%nx_fc_triplet
                map%xuc%x_fc_triplet(i)=gs%gf%ip3%eval(i, gridcoord, rexp=2)
            enddo
            if ( gs%triplet%nconstr .gt. 0 ) then
                !call enforce_constraints( map%xuc%x_fc_triplet,gs%triplet%constraints )
                call lo_enforce_linear_constraints(gs%triplet%constraints,map%xuc%x_fc_triplet)
            endif
        endif
        if ( map%have_fc_quartet ) then
            do i=1,map%xuc%nx_fc_quartet
                map%xuc%x_fc_quartet(i)=gs%gf%ip4%eval(i, gridcoord, rexp=2)
            enddo
            if ( gs%quartet%nconstr .gt. 0 ) then
                !call enforce_constraints( map%xuc%x_fc_quartet,gs%quartet%constraints )
                call lo_enforce_linear_constraints(gs%quartet%constraints,map%xuc%x_fc_quartet)
            endif
        endif
    end select
end subroutine

end module
