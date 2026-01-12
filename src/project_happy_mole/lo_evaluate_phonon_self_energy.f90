module lo_evaluate_phonon_self_energy
use konstanter, only: r8,i8,lo_huge,lo_hugeint,lo_pi,lo_twopi,lo_exitcode_param,lo_freqtol,lo_imag,lo_frequency_Hartree_to_THz
use gottochblandat, only: walltime,lo_progressbar,lo_progressbar_init,lo_trueNtimes,lo_linspace,&
                          lo_gauss,lo_planck,lo_mean,tochar,lo_trapezoid_integration
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh,lo_qpoint
use type_blas_lapack_wrappers, only: lo_gemm
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix
use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
!use lo_brents_method, only: lo_brent_helper
use quadratures_stencils, only: lo_centraldifference, lo_gaussianquadrature
use lo_distributed_phonon_dispersion_relations, only: lo_distributed_phonon_dispersions,lo_distributed_phonon_dispersions_qpoint
use lo_spectralfunction_convolution, only: lo_convolution_handle
use lo_spectralfunction_helpers, only: lo_evaluate_spectral_function, lo_gaussian_smear_spectral_function,&
    lo_tapering_function,lo_interpolated_spectral_function,lo_find_spectral_function_max_and_fwhm
implicit none

private
public :: lo_phonon_selfenergy

!> place to hold all the self-energy stuff, for a single q-point
type lo_phonon_selfenergy
    ! how many frequencies
    integer :: n_energy = -lo_hugeint
    ! number of bands
    integer :: n_mode = -lo_hugeint
    ! energy axis for the self-energy
    real(r8), dimension(:), allocatable :: energy_axis
    ! Imaginary part of self energy (energy,mode)
    real(r8), dimension(:, :), allocatable :: im
    ! Real part of self energy (energy,mode)
    real(r8), dimension(:, :), allocatable :: re
    !> Direction of probe
    real(r8), dimension(3) :: qdir = -lo_huge
    ! Smearing parameters
    real(r8) :: smearing_prefactor = -lo_huge
    ! integrationtype
    integer :: integrationtype = -lo_hugeint
    ! some options
    logical :: isotope_scattering = .false.
    logical :: thirdorder_scattering = .false.
    logical :: fourthorder_scattering = .false.
    ! auxiliary information about the resulting lineshapes:
    !> peak location:
    real(r8), dimension(:), allocatable :: xmid
    !> left fwhm:
    real(r8), dimension(:), allocatable :: xlo
    !> right fwhm:
    real(r8), dimension(:), allocatable :: xhi
    !> before normalization, how far are we from being correct
    real(r8), dimension(:), allocatable :: normalization_residual
    contains
        procedure :: generate
        procedure :: size_in_mem=>se_size_in_mem
        procedure :: destroy=>destroy_se
end type

! Prefactors, to make sure they are the same all the time.
real(r8), parameter :: threephonon_prefactor = lo_pi/2.0_r8
real(r8), parameter :: isotope_prefactor = lo_pi/4.0_r8

interface
    module subroutine generate(se, qpoint, qdir, uc, fc, fct, fcf, ise, qp, ddr,&
        temperature, max_energy, n_energy, integrationtype, sigma,&
        use_isotope, use_thirdorder, use_fourthorder, &
        tmr, mw, mem, verbosity)
        class(lo_phonon_selfenergy), intent(out) :: se
        class(lo_qpoint), intent(in) :: qpoint
        real(r8), dimension(3), intent(in) :: qdir
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_forceconstant_thirdorder), intent(in) :: fct
        type(lo_forceconstant_fourthorder), intent(in) :: fcf
        type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_distributed_phonon_dispersions), intent(in) :: ddr
        real(r8), intent(in) :: temperature
        real(r8), intent(in) :: max_energy
        integer, intent(in) :: n_energy
        integer, intent(in) :: integrationtype
        real(r8), intent(in) :: sigma
        logical, intent(in) :: use_isotope
        logical, intent(in) :: use_thirdorder
        logical, intent(in) :: use_fourthorder
        type(lo_timer), intent(inout) :: tmr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface

interface
    module subroutine threephonon_imaginary_selfenergy(se,qpoint,ompoint,qp,dr,ise,p,fc2,fc3,temperature,mw,mem,verbosity)
        type(lo_phonon_selfenergy), intent(inout) :: se
        class(lo_qpoint), intent(in) :: qpoint
        type(lo_distributed_phonon_dispersions_qpoint), intent(in) :: ompoint
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_distributed_phonon_dispersions), intent(in) :: dr
        type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
        type(lo_crystalstructure), intent(in) :: p
        type(lo_forceconstant_secondorder), intent(inout) :: fc2
        type(lo_forceconstant_thirdorder), intent(in) :: fc3
        real(r8), intent(in) :: temperature
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface

interface ! to helpers
    module subroutine add_quadratic_to_fix_normalization(omega,sigmare,sigmaim,x,initial_residual,mw)
        real(r8), dimension(:), intent(in) :: omega
        real(r8), dimension(:,:), intent(inout) :: sigmare
        real(r8), dimension(:,:), intent(inout) :: sigmaim
        real(r8), dimension(:), intent(in) :: x
        real(r8), dimension(:), intent(out) :: initial_residual
        type(lo_mpi_helper), intent(inout) :: mw
    end subroutine
    module pure function isotope_scattering_strength(p, egv1, egv2, om1, om2) result(lambda)
        type(lo_crystalstructure), intent(in) :: p
        real(r8), intent(in) :: om1
        real(r8), intent(in) :: om2
        complex(r8), dimension(:), intent(in) :: egv1
        complex(r8), dimension(:), intent(in) :: egv2
        real(r8) :: lambda
    end function
end interface

contains

!> size in memory, in bytes
function se_size_in_mem(se) result(mem)
    !> atom
    class(lo_phonon_selfenergy), intent(in) :: se
    !> size in memory, bytes
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(se)

    if (allocated(se%energy_axis)) mem = mem + storage_size(se%energy_axis)*size(se%energy_axis)
    if (allocated(se%im)) mem = mem + storage_size(se%im)*size(se%im)
    if (allocated(se%re)) mem = mem + storage_size(se%re)*size(se%re)
    if (allocated(se%xmid)) mem = mem + storage_size(se%xmid)*size(se%xmid)
    if (allocated(se%xlo)) mem = mem + storage_size(se%xlo)*size(se%xlo)
    if (allocated(se%xhi)) mem = mem + storage_size(se%xhi)*size(se%xhi)
    if (allocated(se%normalization_residual)) mem = mem + storage_size(se%normalization_residual)*size(se%normalization_residual)
    mem = mem/8
end function

subroutine destroy_se(self)
    class(lo_phonon_selfenergy), intent(inout) :: self

    self%n_energy = -lo_hugeint
    self%n_mode = -lo_hugeint
    self%qdir=-lo_huge
    self%smearing_prefactor = -lo_huge
    self%integrationtype = -lo_hugeint
    self%isotope_scattering = .false.
    self%thirdorder_scattering = .false.
    self%fourthorder_scattering = .false.

    if ( allocated(self%energy_axis           ) ) deallocate(self%energy_axis           )
    if ( allocated(self%im                    ) ) deallocate(self%im                    )
    if ( allocated(self%re                    ) ) deallocate(self%re                    )
    if ( allocated(self%xmid                  ) ) deallocate(self%xmid                  )
    if ( allocated(self%xlo                   ) ) deallocate(self%xlo                   )
    if ( allocated(self%xhi                   ) ) deallocate(self%xhi                   )
    if ( allocated(self%normalization_residual) ) deallocate(self%normalization_residual)
end subroutine

end module
