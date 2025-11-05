module lo_evaluate_phonon_self_energy
use konstanter, only: r8, i8, lo_iou, lo_sqtol, lo_twopi, lo_pi, lo_tol, lo_huge, lo_hugeint, lo_freqtol, lo_imag, lo_degenvector, lo_status, &
                      lo_frequency_Hartree_to_THz, lo_exitcode_symmetry, lo_exitcode_param, lo_kappa_au_to_SI, lo_kb_hartree, &
                      lo_A_to_bohr, lo_frequency_THz_to_Hartree
use gottochblandat, only: walltime, lo_progressbar_init, lo_progressbar, lo_lorentz, lo_trapezoid_integration, &
                          lo_linear_interpolation, lo_trueNtimes, lo_gauss, tochar, lo_planck, lo_linspace, &
                          lo_put_function_on_new_axis, lo_looptimer, lo_chop, lo_mean, lo_sqnorm, lo_outerproduct, &
                          lo_harmonic_oscillator_cv, open_file, lo_flattentensor, lo_unflatten_2tensor, &
                          lo_clean_fractional_coordinates
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_operate_on_vector, lo_expandoperation_pair
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_qpoint, lo_LV_tetrahedron_weights, lo_get_small_group_of_qpoint
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use hdf5_wrappers, only: lo_hdf5_helper
use fftw_wrappers, only: lo_convolution, lo_fft, lo_ifft, lo_abcd_convolution

use scatteringrates, only: lo_listofscatteringrates
use lineshape_helper, only: evaluate_spectral_function, taperfn_im, gaussian_smear_spectral_function, &
                            lo_spectralfunction_helper, index_on_grid, find_spectral_function_max_and_fwhm, integrate_spectral_function, &
                            lo_convolution_helper, gaussian_smear_distribution
use type_blas_lapack_wrappers, only: lo_gemm
use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lo_brents_method, only: lo_brent_helper
use quadratures_stencils, only: lo_centraldifference, lo_gaussianquadrature
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
    real(r8), dimension(:, :), allocatable :: im_3ph, im_iso
    ! Real part of self energy (energy,mode)
    real(r8), dimension(:, :), allocatable :: re !re_3ph, re_4ph
    !> Direction of probe
    real(r8), dimension(3) :: qdir
    ! Smearing parameters
    real(r8) :: smearing_prefactor = -lo_huge
    ! integrationtype
    integer :: integrationtype = -lo_hugeint
    ! some options
    logical :: isotope_scattering = .false.
    logical :: thirdorder_scattering = .false.
    logical :: fourthorder_scattering = .false.
    logical :: diagonal = .false.
    logical :: skipsym = .false.
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
end type

! Prefactors, to make sure they are the same all the time.
real(r8), parameter :: threephonon_prefactor = lo_pi/16.0_r8
real(r8), parameter :: fourphonon_prefactor = 1.0_r8/8.0_r8
real(r8), parameter :: isotope_prefactor = lo_pi/4.0_r8

interface
    module subroutine generate(se, qpoint, qdir, wp, uc, fc, fct, fcf, ise, qp, dr, &
        temperature, max_energy, n_energy, integrationtype, sigma,&
        use_isotope, use_thirdorder, use_fourthorder, offdiagonal, &
        tmr, mw, mem, verbosity)
        !> self energy
        class(lo_phonon_selfenergy), intent(out) :: se
        !> qpoint of interest
        class(lo_qpoint), intent(in) :: qpoint
        !> direction of probe?
        real(r8), dimension(3), intent(in) :: qdir
        !> harmonic properties at this q-point
        type(lo_phonon_dispersions_qpoint), intent(in) :: wp
        !> crystal structure
        type(lo_crystalstructure), intent(in) :: uc
        !> second order force constant
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        !> third order force constant
        type(lo_forceconstant_thirdorder), intent(in) :: fct
        !> fourth order force constant
        type(lo_forceconstant_fourthorder), intent(in) :: fcf
        !> tabulated and smeared spectral functions
        type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
        !> q-point mesh
        class(lo_qpoint_mesh), intent(in) :: qp
        !> harmonic properties on this mesh
        type(lo_phonon_dispersions), intent(in) :: dr
        !> temperature
        real(r8), intent(in) :: temperature
        !> max energy to consider
        real(r8), intent(in) :: max_energy
        !> number of energies
        integer, intent(in) :: n_energy
        !> how are we going to integrate
        integer, intent(in) :: integrationtype
        !> adjustment of smearing parameter
        real(r8), intent(in) :: sigma
        !> include isotope scattering
        logical, intent(in) :: use_isotope
        !> include thirdorder scattering
        logical, intent(in) :: use_thirdorder
        !> include fourthorder scattering
        logical, intent(in) :: use_fourthorder
        !> include off-diagonal terms in the self-energy
        logical, intent(in) :: offdiagonal
        !> timer (start and stop outside this routine)
        type(lo_timer), intent(inout) :: tmr
        !> MPI communicator
        type(lo_mpi_helper), intent(inout) :: mw
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> talk a lot
        integer, intent(in) :: verbosity
    end subroutine
end interface

interface
    module subroutine isotope_imaginary_selfenergy(wp,se,qp,dr,sr,mw,mem,verbosity)
        type(lo_phonon_dispersions_qpoint), intent(in) :: wp
        type(lo_phonon_selfenergy), intent(inout) :: se
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_listofscatteringrates), intent(in) :: sr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine threephonon_imaginary_selfenergy(wp,se,qp,dr,sr,ise,p,temperature,mw,mem,verbosity)
        type(lo_phonon_dispersions_qpoint), intent(in) :: wp
        type(lo_phonon_selfenergy), intent(inout) :: se
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_listofscatteringrates), intent(in) :: sr
        type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
        type(lo_crystalstructure), intent(in) :: p
        real(r8), intent(in) :: temperature
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine fourphonon_real_selfenergy(qpoint, wp, gp, qp, uc, temperature, dr, fcf, delta, skipsym, mw, mem, verbosity)
        type(lo_qpoint), intent(in) :: qpoint
        type(lo_phonon_dispersions_qpoint), intent(in) :: wp
        type(lo_phonon_dispersions_qpoint), intent(in) :: gp
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_forceconstant_fourthorder), intent(in) :: fcf
        real(r8), intent(in) :: temperature
        real(r8), dimension(:), intent(out) :: delta
        logical, intent(in) :: skipsym
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface

interface ! to helpers
    module subroutine tapering_function(x,y)
        real(r8), dimension(:), intent(in) :: x
        real(r8), dimension(:), intent(out) :: y
    end subroutine
    module subroutine add_quadratic_to_fix_normalization(omega,sigmare,sigmaim,x,initial_residual,mw)
        !> harmonic frequencies
        real(r8), dimension(:), intent(in) :: omega
        !> real part of self-energy
        real(r8), dimension(:,:), intent(inout) :: sigmare
        !> imaginary part of self-energy
        real(r8), dimension(:,:), intent(inout) :: sigmaim
        !> energy axis
        real(r8), dimension(:), intent(in) :: x
        !> spectral function residual
        real(r8), dimension(:), intent(out) :: initial_residual
        !> mpi helper
        type(lo_mpi_helper), intent(inout) :: mw
    end subroutine
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
    if (allocated(se%im_3ph)) mem = mem + storage_size(se%im_3ph)*size(se%im_3ph)
    if (allocated(se%im_iso)) mem = mem + storage_size(se%im_iso)*size(se%im_iso)
    if (allocated(se%re)) mem = mem + storage_size(se%re)*size(se%re)
    !if (allocated(se%re_3ph)) mem = mem + storage_size(se%re_3ph)*size(se%re_3ph)
    !if (allocated(se%re_4ph)) mem = mem + storage_size(se%re_4ph)*size(se%re_4ph)
    if (allocated(se%xmid)) mem = mem + storage_size(se%xmid)*size(se%xmid)
    if (allocated(se%xlo)) mem = mem + storage_size(se%xlo)*size(se%xlo)
    if (allocated(se%xhi)) mem = mem + storage_size(se%xhi)*size(se%xhi)
    if (allocated(se%normalization_residual)) mem = mem + storage_size(se%normalization_residual)*size(se%normalization_residual)
    mem = mem/8
end function

end module
