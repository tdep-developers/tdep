module lo_collision_matrix
!! Container for things related to thermal transport. Not how it is actually calculated,
!! that is generated elsewhere, but how it is stored and written to file.
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_freqtol, lo_pi, lo_exitcode_param, &
                      lo_groupvel_Hartreebohr_to_ms, lo_kb_Hartree, lo_kappa_au_to_SI, lo_frequency_Hartree_to_THz, &
                      lo_frequency_Hartree_to_icm, lo_frequency_Hartree_to_meV, lo_time_au_to_s, lo_bohr_to_A
use gottochblandat, only: lo_trapezoid_integration
! use gottochblandat, only: tochar, walltime, lo_chop, &
!                           lo_real_nullspace_coefficient_matrix, lo_transpositionmatrix, lo_real_pseudoinverse, &
!                           lo_flattentensor, lo_trapezoid_integration, lo_linspace, lo_outerproduct, lo_planck, &
!                           lo_harmonic_oscillator_cv, lo_unflatten_2tensor, lo_linear_interpolation
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh,lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions,lo_phonon_dispersions_qpoint
use type_symmetryoperation, only: lo_expandoperation_pair, lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm
!use lineshape_helper, only: evaluate_spectral_function, taperfn_im, integrate_spectral_function, integrate_two_spectral_functions
use quadratures_stencils, only: lo_centraldifference

use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lo_spectralfunction_convolution, only: lo_convolution_handle
use lo_spectralfunction_helpers, only: lo_evaluate_spectral_function,lo_gaussian_smear_spectral_function
implicit none

private
public :: lo_scattering_matrix

type lo_scattering_matrix
    contains
        procedure :: generate=>create_scattering_matrix
end type

type lo_scattering_matrix_convolutions
end type

contains

subroutine create_scattering_matrix(scm,p,fc2,fc3,ise,temperature,qp,mw,mem,verbosity)
    !> scattering matrix
    class(lo_scattering_matrix), intent(inout) :: scm
    !> structure
    type(lo_crystalstructure), intent(inout) :: p
    !> second order forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc2
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fc3
    !> interpolated self-energy
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> temperature
    real(r8), intent(in) :: temperature
    !> q-point grid
    class(lo_qpoint_mesh), intent(inout) :: qp
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_phonon_dispersions) :: dr
    type(lo_convolution_handle) :: ch

    init: block
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'Creating scattering matrix'
        endif

        ! Going to need harmonic dispersions
        call dr%generate(qp,fc2,p,mw,mem,verbosity=-1)

        ! And a convolution helper
        call ch%generate(ise%omega,temperature,dr%n_mode)
    end block init

    !call ise%evaluate()

end subroutine

!> calculate a specific entry in the scattering matrix
subroutine scattering_matrix_entry(iq,jq,p,qp,dr,ise,ch,adaptive_prefactor,mem)
    !> index of first q-point, irreducible
    integer, intent(in) :: iq
    !> index of second q-point, full
    integer, intent(in) :: jq
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point grid
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> interpolated self-energy
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> convolution handle
    type(lo_convolution_handle), intent(inout) :: ch
    !> adaptive smearing prefactor
    real(r8), intent(in) :: adaptive_prefactor
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(3) :: v0
    integer, dimension(3) :: gi
    integer :: kq

    init: block

        ! First step, I guess, is figure out the third q-point?
        select type(qp)
        type is(lo_fft_mesh)
            v0=qp%ip(iq)%r + qp%ap(jq)%r
            v0=matmul(p%inv_reciprocal_latticevectors,v0)
            gi=qp%index_from_coordinate(v0)
            ! Third q-point index, good enough I guess.
            kq=qp%gridind2ind(gi(1),gi(2),gi(3))
        end select
    end block init

    interpolate: block
        real(r8), dimension(:,:), allocatable :: buf_re,buf_im,buf_j
        real(r8) :: sigma,f0
        integer :: imode

        allocate(buf_re(ise%n_energy,dr%n_mode))
        allocate(buf_im(ise%n_energy,dr%n_mode))
        allocate(buf_j(ise%n_energy,dr%n_mode))
        buf_re=0.0_r8
        buf_im=0.0_r8
        buf_j=0.0_r8

        ! This is the first q-point, fetch self-energies
        call ise%evaluate(p,qp%ip(iq)%r,dr%iq(iq)%omega,dr%iq(iq)%egv,buf_re,buf_im,mem)
        ! Turn self-energy into smeared normalized spectral functions
        do imode=1,dr%n_mode
            ! skip acoustic
            if ( dr%iq(iq)%omega(imode) .lt. lo_freqtol ) cycle
            ! raw spectral function
            call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),dr%iq(iq)%omega(imode),buf_j(:,imode))
            ! smear the spectral function
            sigma=qp%smearingparameter( dr%iq(iq)%vel(:,imode), dr%default_smearing(imode), adaptive_prefactor)
            call lo_gaussian_smear_spectral_function(ise%omega,sigma,buf_j(:,imode))
            ! and normalize? why not
            f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
        enddo

        ! Pre-buffer in convolution helper
        !call ch%buffer_and_transform_J(buf_j,1)

        ! Repeat for second q-point
        call ise%evaluate(p,qp%ap(jq)%r,dr%aq(jq)%omega,dr%aq(jq)%egv,buf_re,buf_im,mem)
        do imode=1,dr%n_mode
            if ( dr%aq(jq)%omega(imode) .lt. lo_freqtol ) cycle
            call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),dr%aq(jq)%omega(imode),buf_j(:,imode))
            sigma=qp%smearingparameter( dr%aq(jq)%vel(:,imode), dr%default_smearing(imode), adaptive_prefactor)
            call lo_gaussian_smear_spectral_function(ise%omega,sigma,buf_j(:,imode))
            f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
        enddo
        !call ch%buffer_and_transform_J(buf_j,2)

        ! And for the third q-point
        call ise%evaluate(p,qp%ap(kq)%r,dr%aq(kq)%omega,dr%aq(kq)%egv,buf_re,buf_im,mem)
        do imode=1,dr%n_mode
            if ( dr%aq(kq)%omega(imode) .lt. lo_freqtol ) cycle
            call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),dr%aq(kq)%omega(imode),buf_j(:,imode))
            sigma=qp%smearingparameter( dr%aq(kq)%vel(:,imode), dr%default_smearing(imode), adaptive_prefactor)
            call lo_gaussian_smear_spectral_function(ise%omega,sigma,buf_j(:,imode))
            f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
        enddo
        !call ch%buffer_and_transform_J(buf_j,3)
    end block interpolate

    matrixelement: block
        ! Then I have to evaluate the matrix element. Hmmm.
    end block matrixelement

end subroutine

end module
