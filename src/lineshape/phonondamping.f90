module phonondamping
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
use type_qpointmesh, only: lo_qpoint_mesh, lo_qpoint, lo_LV_tetrahedron_weights, lo_get_small_group_of_qpoint
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_phonon_dos, only: lo_phonon_dos
use hdf5_wrappers, only: lo_hdf5_helper
use fftw_wrappers, only: lo_convolution, lo_fft, lo_ifft, lo_abcd_convolution

use options, only: lo_opts
use scatteringrates, only: lo_listofscatteringrates
use lo_realspace_selfenergy, only: lo_interpolated_selfenergy
use lo_thermal_transport, only: lo_thermal_conductivity
use lineshape_helper, only: evaluate_spectral_function, taperfn_im, gaussian_smear_spectral_function, &
                            lo_spectralfunction_helper, index_on_grid, find_spectral_function_max_and_fwhm, integrate_spectral_function, &
                            lo_convolution_helper, gaussian_smear_distribution
use type_blas_lapack_wrappers, only: lo_gemm
implicit none

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
    real(r8), dimension(:, :), allocatable :: re_3ph, re_4ph
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
    ! auxiliary information about lineshapes:
    !> peak location:
    real(r8), dimension(:), allocatable :: xmid
    !> left fwhm:
    real(r8), dimension(:), allocatable :: xlo
    !> right fwhm:
    real(r8), dimension(:), allocatable :: xhi
    !> scaling factor to get proper normalization
    real(r8), dimension(:), allocatable :: scalingfactor
contains
    procedure :: generate
    procedure :: generate_interp
    procedure :: size_in_mem => se_size_in_mem
end type

! Prefactors, to make sure they are the same all the time.
real(r8), parameter :: threephonon_prefactor = lo_pi/16.0_r8
real(r8), parameter :: fourphonon_prefactor = 1.0_r8/8.0_r8
real(r8), parameter :: isotope_prefactor = lo_pi/4.0_r8

! dos interfaces
interface
    module subroutine get_intensity_as_dos(pd, qpd, drd, uc, fc, fct, fcf, ise, sf, qp, dr, opts, tmr, mw, mem, verbosity)
        type(lo_phonon_dos), intent(out) :: pd
        class(lo_qpoint_mesh), intent(in) :: qpd
        type(lo_phonon_dispersions), intent(inout) :: drd
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_forceconstant_thirdorder), intent(in) :: fct
        type(lo_forceconstant_fourthorder), intent(in) :: fcf
        type(lo_interpolated_selfenergy), intent(in) :: ise
        type(lo_spectralfunction_helper), intent(out) :: sf
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_opts), intent(in) :: opts
        type(lo_timer), intent(inout) :: tmr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine get_intensity_as_dos_interp(pd, tc, qp, dr, uc, ise, opts, mw, mem, verbosity)
        type(lo_phonon_dos), intent(out) :: pd
        type(lo_thermal_conductivity), intent(out) :: tc
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(inout) :: dr
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_interpolated_selfenergy), intent(in) :: ise
        type(lo_opts), intent(in) :: opts
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface

! grid interfaces
interface
    module subroutine get_selfenergy_on_closed_grid(sf, tc, pd, qp, dr, uc, fc, fct, fcf, ise, opts, tmr, mw, mem, verbosity)
        type(lo_spectralfunction_helper), intent(out) :: sf
        type(lo_thermal_conductivity), intent(out) :: tc
        type(lo_phonon_dos), intent(out) :: pd
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(inout) :: dr
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_forceconstant_thirdorder), intent(in) :: fct
        type(lo_forceconstant_fourthorder), intent(in) :: fcf
        type(lo_interpolated_selfenergy), intent(in) :: ise
        type(lo_opts), intent(in) :: opts
        type(lo_timer), intent(inout) :: tmr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine get_selfenergy_on_points(qvec, wp, qp, dr, uc, fc, fct, fcf, ise, opts, sigmaRe, sigmaIm, tmr, mw, mem, verbosity)
        real(r8), dimension(:, :), intent(in) :: qvec
        type(lo_phonon_dispersions_qpoint), dimension(:), intent(inout) :: wp
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(inout) :: dr
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_forceconstant_thirdorder), intent(in) :: fct
        type(lo_forceconstant_fourthorder), intent(in) :: fcf
        type(lo_interpolated_selfenergy), intent(in) :: ise
        type(lo_opts), intent(in) :: opts
        real(r8), dimension(:, :, :), intent(inout) :: sigmaRe
        real(r8), dimension(:, :, :), intent(inout) :: sigmaIm
        type(lo_timer), intent(inout) :: tmr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface

! path interfaces
interface
    module subroutine spectral_function_along_path(bs, uc, fc, fct, fcf, ise, qp, dr, opts, tmr, mw, mem)
        type(lo_phonon_bandstructure), intent(inout) :: bs
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_forceconstant_thirdorder), intent(in) :: fct
        type(lo_forceconstant_fourthorder), intent(in) :: fcf
        type(lo_interpolated_selfenergy), intent(in) :: ise
        class(lo_qpoint_mesh), intent(in) :: qp
        type(lo_phonon_dispersions), intent(in) :: dr
        type(lo_opts), intent(in) :: opts
        type(lo_timer), intent(inout) :: tmr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine spectral_function_along_path_interp(bs, uc, ise, opts, mw, mem)
        type(lo_phonon_bandstructure), intent(inout) :: bs
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_interpolated_selfenergy), intent(in) :: ise
        type(lo_opts), intent(in) :: opts
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
end interface

contains

#include "phonondamping_threephonon.f90"
#include "phonondamping_generation.f90"
#include "phonondamping_gaussian.f90"
#include "phonondamping_tetrahedron.f90"
#include "phonondamping_fourthorder.f90"

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
    if (allocated(se%re_3ph)) mem = mem + storage_size(se%re_3ph)*size(se%re_3ph)
    if (allocated(se%re_4ph)) mem = mem + storage_size(se%re_4ph)*size(se%re_4ph)
    if (allocated(se%xmid)) mem = mem + storage_size(se%xmid)*size(se%xmid)
    if (allocated(se%xlo)) mem = mem + storage_size(se%xlo)*size(se%xlo)
    if (allocated(se%xhi)) mem = mem + storage_size(se%xhi)*size(se%xhi)
    if (allocated(se%scalingfactor)) mem = mem + storage_size(se%scalingfactor)*size(se%scalingfactor)
    mem = mem/8
end function

end module
