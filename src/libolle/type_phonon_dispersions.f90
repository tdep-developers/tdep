module type_phonon_dispersions
!!
!! Handles phonon dispersion relations on a q-point mesh
!!
use konstanter, only: r8, i8, lo_huge, lo_hugeint, lo_degenvector, lo_status, lo_freqtol, lo_tol, &
                      lo_sqtol, lo_temperaturetol, &
                      lo_status, lo_hbar_hartree, lo_frequency_Hartree_to_Hz, lo_groupvel_HartreeBohr_to_ms, &
                      lo_exitcode_param, lo_hartree_to_eV, lo_bohr_to_m, lo_kappa_au_to_SI, lo_kb_hartree, &
                      lo_Hartree_to_Joule, lo_exitcode_physical, lo_imag, lo_time_s_to_au
use gottochblandat, only: open_file, walltime, tochar, lo_trueNtimes, lo_classical_harmonic_oscillator_free_energy, &
                          lo_harmonic_oscillator_cv, lo_harmonic_oscillator_entropy, lo_harmonic_oscillator_free_energy, &
                          lo_planck, lo_sqnorm, lo_progressbar_init, lo_progressbar, qsort, lo_outerproduct, lo_chop, &
                          lo_symmetric_eigensystem_3x3matrix, lo_flattentensor, lo_unflatten_2tensor, &
                          lo_linear_least_squares, lo_nullspace_coefficient_matrix, lo_identitymatrix, &
                          lo_real_pseudoinverse, lo_determ
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_blas_lapack_wrappers, only: lo_gemm, lo_gemv
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint, lo_qpoint_mesh, lo_monkhorst_pack_mesh, lo_wedge_mesh, lo_fft_mesh, &
                           lo_get_small_group_of_qpoint
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix, lo_operate_on_vector, lo_expandoperation_pair
use hdf5_wrappers, only: lo_hdf5_helper

implicit none
private
public :: lo_phonon_dispersions_qpoint
public :: lo_phonon_dispersions

!> A q-point in dispersion relations
type lo_phonon_dispersions_qpoint
    !> frequencies ( in angular frequency in Hz)
    real(r8), dimension(:), allocatable :: omega
    !> group velocities
    real(r8), dimension(:, :), allocatable :: vel
    !> mode eigenvectors
    complex(r8), dimension(:, :), allocatable :: egv
    !> mode gruneisen parameter
    real(r8), dimension(:), allocatable :: gruneisen
    !> linewidth
    real(r8), dimension(:), allocatable :: linewidth
    !> anharmonic frequency shift
    real(r8), dimension(:), allocatable :: shift3, shift4
    !> plus scattering rate
    real(r8), dimension(:), allocatable :: p_plus
    !> minus scattering rate
    real(r8), dimension(:), allocatable :: p_minus
    !> isotope scattering rate
    real(r8), dimension(:), allocatable :: p_iso
    !> QS parameter for thermal conductivity
    real(r8), dimension(:), allocatable :: qs
    !> Helper for thermal conductivity
    real(r8), dimension(:, :), allocatable :: F0
    !> Helper for thermal conductivity
    real(r8), dimension(:, :), allocatable :: Fn
    !> mode decomposed thermal conductivity
    real(r8), dimension(:, :, :), allocatable :: kappa
    !> harmonic free energy
    real(r8), dimension(:), allocatable :: F
    !> harmonic entropy
    real(r8), dimension(:), allocatable :: S
    !> anharmonic free energy
    real(r8), dimension(:), allocatable :: deltaF3, deltaF4
    !> anharmonic entropy
    real(r8), dimension(:), allocatable :: deltaS3, deltaS4
    !> electron-phonon phasespace
    real(r8), dimension(:), allocatable :: electronphononphasespace
    !> frequencies from elastic dispersions
    real(r8), dimension(:), allocatable :: omega_elastic
    !> three-phonon phasespace
    real(r8), dimension(:), allocatable :: threephononphasespace
    !> phonon mean free path
    real(r8), dimension(:, :), allocatable :: mfp
    !> scalar phonon mean free path
    real(r8), dimension(:), allocatable :: scalar_mfp
    !> thermal pre-factor for neutron scattering
    real(r8), dimension(:), allocatable :: thermal_prefactor
    !> pyroelectric amplitude
    real(r8), dimension(:), allocatable :: pyroelectric_amplitude

    !> how degenerate is the mode?
    integer, dimension(:), allocatable :: degeneracy
    !> with what mode is it degenerate?
    integer, dimension(:, :), allocatable :: degenmode

    !> for vectorial quantities, we have a matrix that can fix things
    real(r8), dimension(3, 3) :: invariance_fixer = -lo_huge

    !> Hessian in q-space, derivative of group velocity
    !real(r8), dimension(:,:,:), allocatable :: qhessian
contains
    !> calculate all normal phonon things for a single point
    procedure :: generate => harmonic_things_for_single_point
    !> count how large the point becomes once packed into a character buffer
    procedure :: size_packed => phonon_dispersions_qpoint_size_packed
    !> pack this point to a buffer
    procedure :: pack_to_buf => pack_phonon_dispersions_qpoint
    !> unpack this point from a buffer
    procedure :: unpack_from_buf => unpack_phonon_dispersions_qpoint
    !> size in memory, in bytes
    procedure :: size_in_mem => phonon_dispersions_qpoint_size_in_mem
end type

!> Phonon dispersion relations in the full BZ
type lo_phonon_dispersions
    !> how many phonon modes
    integer :: n_mode = -lo_hugeint
    !> number of points in the irreducible part
    integer :: n_irr_qpoint = -lo_hugeint
    !> total number of points
    integer :: n_full_qpoint = -lo_hugeint
    !> points in the irreducible part
    type(lo_phonon_dispersions_qpoint), dimension(:), allocatable :: iq
    !> points in the full zone
    type(lo_phonon_dispersions_qpoint), dimension(:), allocatable :: aq
    !> largest frequency
    real(r8) :: omega_max = -lo_huge
    !> smallest nonzero frequency
    real(r8) :: omega_min = -lo_huge
    !> sensible default smearing per band
    real(r8), dimension(:), allocatable :: default_smearing
contains
    !> create the full dispersions
    procedure :: generate
    !> calculate gruneisen parameters
    !procedure :: gruneisen
    !> phonon free energy
    procedure :: phonon_free_energy
    !> classical phonon free energy
    procedure :: phonon_free_energy_classical
    !> phonon entropy
    procedure :: phonon_entropy
    !> phonon heat capacity
    procedure :: phonon_cv
    !> phonon angular momentum
    procedure :: phonon_angular_momentum_matrix
    !> thermal displacement covariance matrix
    procedure :: thermal_displacement_matrix
    !> store the full grid in a hdf file
    procedure :: write_to_hdf5
    !> store the irreducible grid in a hdf5 file
    procedure :: write_irreducible_to_hdf5
    !> make sure things that are supposed to be degenerate really are
    procedure :: enforce_degeneracy
    !> communication helper
    procedure :: allgather_irreducible
    !> communication helper
    procedure :: allgather_fullmesh
    !> measure size in memory, in bytes
    procedure :: size_in_mem => phonon_dispersions_size_in_mem
    !> destroy
    procedure :: destroy
end type

!> interfaces to phonon_dispersion_relations_generation
interface
    module subroutine generate(dr, qp, fc, p, mw, mem, verbosity)
        class(lo_phonon_dispersions), intent(out) :: dr
        class(lo_qpoint_mesh), intent(inout) :: qp
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_crystalstructure), intent(inout) :: p
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    ! module subroutine gruneisen(dr,qp,fct,p,verbosity,mpi_communicator)
    !     class(lo_phonon_dispersions), intent(inout) :: dr
    !     class(lo_qpoint_mesh), intent(inout) :: qp
    !     type(lo_forceconstant_thirdorder), intent(in) :: fct
    !     type(lo_crystalstructure), intent(in) :: p
    !     integer, intent(in), optional :: verbosity
    !     integer, intent(in), optional :: mpi_communicator
    ! end subroutine
    module subroutine allgather_irreducible(dr, mw, mem)
        class(lo_phonon_dispersions), intent(inout) :: dr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine allgather_fullmesh(dr, mw, mem)
        class(lo_phonon_dispersions), intent(inout) :: dr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
end interface

contains

! !> quick test to check wether the dispersions are stable at all q-points.
! function is_unstable(dr,qp) result(unstable)
!     !> dispersions
!     class(lo_phonon_dispersions), intent(in) :: dr
!     !> q-mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> unstable?
!     logical :: unstable
!
!     integer :: i,j,ctr
!
!     unstable=.false.
!     ctr=0
!     do i=1,dr%n_full_qpoint
!     do j=1,dr%n_mode
!         if ( dr%aq(i)%omega(j) .lt. lo_freqtol ) then
!             ctr=ctr+1
!             if ( lo_sqnorm(qp%ap(i)%r) .gt. lo_sqtol ) unstable=.true.
!         endif
!     enddo
!     enddo
!     do i=1,dr%n_irr_qpoint
!     do j=1,dr%n_mode
!         if ( dr%iq(i)%omega(j) .lt. lo_freqtol ) then
!             ctr=ctr+1
!             if ( lo_sqnorm(qp%ip(i)%r) .gt. lo_sqtol ) unstable=.true.
!         endif
!     enddo
!     enddo
!     ! 6 or 8. Why? 1-D materials have 4 acoustic branches.
!     if ( ctr .gt. 8 ) unstable=.true.
! end function

!> make sure degenerate things are as degenerate as they should be. Not all quantities right now, but most. Also only works on the irreducible wedge.
subroutine enforce_degeneracy(dr)
    !> dispersions
    class(lo_phonon_dispersions), intent(inout) :: dr

    integer :: i, j, b1, b2
    integer, dimension(dr%n_mode) :: di
    real(r8) :: p_plus, p_minus, p_iso, avgfactor
    real(r8), dimension(3) :: vel, F0, Fn
    logical :: do_velocity, do_scatteringrate, do_F0, do_Fn

    ! average velocities?
    if (allocated(dr%iq(1)%vel)) then
        do_velocity = .true.
    else
        do_velocity = .false.
    end if
    ! average scatteringrates?
    if (allocated(dr%iq(1)%p_plus)) then
        do_scatteringrate = .true.
    else
        do_scatteringrate = .false.
    end if
    if (allocated(dr%iq(1)%F0)) then
        do_F0 = .true.
    else
        do_F0 = .false.
    end if
    if (allocated(dr%iq(1)%Fn)) then
        do_Fn = .true.
    else
        do_Fn = .false.
    end if

    do i = 1, dr%n_irr_qpoint
        ! what is the max level of degeneracy here?
        di = 1
        do b1 = 1, dr%n_mode
            ! skip if already taken care of
            if (di(b1) .eq. 0) cycle
            ! now average and spread over this mode, if necessary
            if (dr%iq(i)%degeneracy(b1) .gt. 1) then
                p_plus = 0.0_r8
                p_minus = 0.0_r8
                p_iso = 0.0_r8
                vel = 0.0_r8
                F0 = 0.0_r8
                Fn = 0.0_r8
                avgfactor = 1.0_r8/(dr%iq(i)%degeneracy(b1))
                ! average it
                do j = 1, dr%iq(i)%degeneracy(b1)
                    b2 = dr%iq(i)%degenmode(j, b1)
                    if (do_velocity) then
                        vel = vel + dr%iq(i)%vel(:, b2)*avgfactor
                    end if
                    if (do_scatteringrate) then
                        p_plus = p_plus + dr%iq(i)%p_plus(b2)*avgfactor
                        p_minus = p_minus + dr%iq(i)%p_minus(b2)*avgfactor
                        p_iso = p_iso + dr%iq(i)%p_iso(b2)*avgfactor
                    end if
                    if (do_F0) then
                        F0 = F0 + dr%iq(i)%F0(:, b2)*avgfactor
                    end if
                    if (do_Fn) then
                        Fn = Fn + dr%iq(i)%Fn(:, b2)*avgfactor
                    end if
                end do
                ! spread  it out
                do j = 1, dr%iq(i)%degeneracy(b1)
                    b2 = dr%iq(i)%degenmode(j, b1)
                    if (do_velocity) then
                        dr%iq(i)%vel(:, b2) = vel
                    end if
                    if (do_scatteringrate) then
                        dr%iq(i)%p_plus(b2) = p_plus
                        dr%iq(i)%p_minus(b2) = p_minus
                        dr%iq(i)%p_iso(b2) = p_iso
                    end if
                    if (do_F0) then
                        dr%iq(i)%F0(:, b2) = F0
                    end if
                    if (do_Fn) then
                        dr%iq(i)%Fn(:, b2) = Fn
                    end if
                    ! don't bother with this mode anymore
                    di(b2) = 0
                end do
            end if
        end do
    end do
end subroutine

!> write everything on the mesh to an hdf5 file
subroutine write_to_hdf5(dr, qp, uc, filename, mem, temperature)
    !> phonon dispersions
    class(lo_phonon_dispersions), intent(in) :: dr
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> filename
    character(len=*), intent(in) :: filename
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> optional temperature
    real(r8), intent(in), optional :: temperature

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:, :, :, :), allocatable :: dddd
    real(r8), dimension(:, :, :), allocatable :: ddd
    real(r8), dimension(:, :), allocatable :: dd
    real(r8), dimension(3) :: v0, v1
    real(r8) :: n1, n, f0, omega, cv, tau
    integer :: i, j, k, l
    character(len=1000) :: dname

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))
    ! Some general attributes first:
    call h5%store_attribute(dr%n_mode/3, h5%file_id, 'number_of_atoms', lo_status)
    call h5%store_attribute(dr%n_mode, h5%file_id, 'number_of_bands', lo_status)
    call h5%store_attribute(dr%n_full_qpoint, h5%file_id, 'number_of_qpoints', lo_status)
    select type (qp)
    type is (lo_monkhorst_pack_mesh)
        call h5%store_attribute('Monkhorst-Pack', h5%file_id, 'meshtype', lo_status)
    type is (lo_fft_mesh)
        call h5%store_attribute('FFT', h5%file_id, 'meshtype', lo_status)
    type is (lo_wedge_mesh)
        call h5%store_attribute('Wedge', h5%file_id, 'meshtype', lo_status)
    end select

    ! q-points
    dname = 'qpoints'
    call mem%allocate(dd, [3, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    do i = 1, qp%n_full_point
        dd(:, i) = qp%ap(i)%r
    end do
    call h5%store_data(dd, h5%file_id, trim(dname), enhet='1/A', dimensions='q-vector,xyz')
    call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! frequencies
    if (allocated(dr%aq(1)%omega)) then
        dname = 'frequencies'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            dd(:, i) = dr%aq(i)%omega*lo_frequency_hartree_to_Hz
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz (angular)', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Group velocities
    if (allocated(dr%aq(1)%vel)) then
        dname = 'group_velocities'
        call mem%allocate(ddd, [3, dr%n_mode, qp%n_full_point], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            ddd(:, :, i) = dr%aq(i)%vel*lo_groupvel_hartreebohr_to_ms
        end do
        call h5%store_data(ddd, h5%file_id, trim(dname), enhet='m/s', dimensions='q-vector,mode,xyz')
        call mem%deallocate(ddd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Eigenvectors
    if (allocated(dr%aq(1)%egv)) then
        dname = 'eigenvectors_re'
        call mem%allocate(ddd, [dr%n_mode, dr%n_mode, qp%n_full_point], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            ddd(:, :, i) = real(dr%aq(i)%egv)
        end do
        call h5%store_data(ddd, h5%file_id, trim(dname), enhet='dimensionless', dimensions='q-vector,mode,atom-xyz')
        dname = 'eigenvectors_im'
        do i = 1, qp%n_full_point
            ddd(:, :, i) = aimag(dr%aq(i)%egv)
        end do
        call h5%store_data(ddd, h5%file_id, trim(dname), enhet='dimensionless', dimensions='q-vector,mode,atom-xyz')
        call mem%deallocate(ddd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Gruneisen parameters
    if (allocated(dr%aq(1)%gruneisen)) then
        dname = 'gruneisen_parameters'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            dd(:, i) = dr%aq(i)%gruneisen
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='dimensionless', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Lifetimes
    if (allocated(dr%iq(1)%linewidth)) then
        dname = 'lifetimes'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            do j = 1, dr%n_mode
                if (dr%iq(qp%ap(i)%irreducible_index)%linewidth(j) .gt. lo_freqtol/10) then
                    dd(j, i) = 0.5_r8/(dr%iq(qp%ap(i)%irreducible_index)%linewidth(j)*lo_frequency_hartree_to_Hz)
                else
                    dd(j, i) = 0.0_r8
                end if
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='s', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! unit cell vectors
    dname = 'lattice_vectors'
    call mem%allocate(dd, [3, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ! TODO: check transposition convention for lattice vectors
    dd = lo_bohr_to_m*1.d10*transpose(uc%latticevectors)
    call h5%store_data(dd, h5%file_id, trim(dname), enhet='A', dimensions='xyz,xyz')
    call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Linewidths
    if (allocated(dr%iq(1)%linewidth)) then
        dname = 'linewidths'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            do j = 1, dr%n_mode
                dd(j, i) = dr%iq(qp%ap(i)%irreducible_index)%linewidth(j)*lo_frequency_hartree_to_Hz
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Shifts
    if (allocated(dr%iq(1)%shift3)) then
        dname = 'shifts_thirdorder'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            do j = 1, dr%n_mode
                dd(j, i) = dr%iq(qp%ap(i)%irreducible_index)%shift3(j)*lo_frequency_hartree_to_Hz
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if
    if (allocated(dr%iq(1)%shift4)) then
        dname = 'shifts_fourthorder'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            do j = 1, dr%n_mode
                dd(j, i) = dr%iq(qp%ap(i)%irreducible_index)%shift4(j)*lo_frequency_hartree_to_Hz
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Anharmonic free energy
    if (allocated(dr%iq(1)%deltaF3)) then
        dname = 'anharmonic_free_energy_thirdorder'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            do j = 1, dr%n_mode
                dd(j, i) = dr%iq(qp%ap(i)%irreducible_index)%deltaF3(j)*lo_hartree_to_eV
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='eV', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if
    if (allocated(dr%iq(1)%deltaF4)) then
        dname = 'anharmonic_free_energy_fourthorder'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            do j = 1, dr%n_mode
                dd(j, i) = dr%iq(qp%ap(i)%irreducible_index)%deltaF4(j)*lo_hartree_to_eV
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='eV', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Mean free path
    if (allocated(dr%iq(1)%linewidth)) then
        dname = 'mean_free_paths'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            do j = 1, dr%n_mode
                if (dr%iq(qp%ap(i)%irreducible_index)%linewidth(j) .gt. lo_freqtol) then
                    dd(j, i) = (norm2(dr%aq(i)%vel(:, j))*0.5_r8)/(dr%iq(qp%ap(i)%irreducible_index)%linewidth(j))
                    if (dd(j, i) .gt. lo_tol) then
                        dd(j, i) = dd(j, i)*lo_Bohr_to_m
                    else
                        dd(j, i) = 0.0_r8
                    end if
                else
                    dd(j, i) = 0.0_r8
                end if
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='m', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Thermal conductivity stuff
    if (allocated(dr%iq(1)%p_plus) .and. present(temperature)) then
        dname = 'scattering_rates_plus'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
            k = qp%ap(i)%irreducible_index
            do j = 1, dr%n_mode
                if (dr%aq(i)%omega(j) .lt. dr%omega_min*0.5_r8) then
                    dd(j, i) = 0.0_r8
                else
                    n1 = lo_planck(temperature, dr%aq(i)%omega(j))
                    dd(j, i) = dr%iq(k)%p_plus(j)*lo_frequency_hartree_to_Hz/(n1*(n1 + 1.0_r8))
                end if
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    if (allocated(dr%iq(1)%p_minus) .and. present(temperature)) then
        dname = 'scattering_rates_minus'
        do i = 1, qp%n_full_point
            k = qp%ap(i)%irreducible_index
            do j = 1, dr%n_mode
                if (dr%aq(i)%omega(j) .lt. dr%omega_min*0.5_r8) then
                    dd(j, i) = 0.0_r8
                else
                    n1 = lo_planck(temperature, dr%aq(i)%omega(j))
                    dd(j, i) = dr%iq(k)%p_minus(j)*lo_frequency_hartree_to_Hz/(n1*(n1 + 1.0_r8))
                end if
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz', dimensions='q-vector,mode')
        dname = 'scattering_rates_isotope'
        do i = 1, qp%n_full_point
            k = qp%ap(i)%irreducible_index
            do j = 1, dr%n_mode
                if (dr%aq(i)%omega(j) .lt. dr%omega_min*0.5_r8) then
                    dd(j, i) = 0.0_r8
                else
                    n1 = lo_planck(temperature, dr%aq(i)%omega(j))
                    dd(j, i) = dr%iq(k)%p_iso(j)*lo_frequency_hartree_to_Hz/(n1*(n1 + 1.0_r8))
                end if
            end do
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    if (allocated(dr%iq(1)%Fn) .and. present(temperature)) then
        call mem%allocate(dddd, [3, 3, dr%n_mode, qp%n_full_point], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        dname = 'thermal_conductivity'
        dddd = 0.0_r8
        do i = 1, qp%n_full_point
            l = qp%ap(i)%irreducible_index
            k = qp%ap(i)%operation_from_irreducible
            do j = 1, dr%n_mode
                if (dr%aq(i)%omega(j) .lt. dr%omega_min*0.5_r8) cycle
                if (k .gt. 0) then
                    v0 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%Fn(:, j), reciprocal=.true.)
                    v1 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%vel(:, j), reciprocal=.true.)
                else
                    v0 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%Fn(:, j), reciprocal=.true.)
                    v1 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%vel(:, j), reciprocal=.true.)
                end if
                ! Get kappa for this q-point
                omega = dr%iq(l)%omega(j)
                n = lo_planck(temperature, omega)
                f0 = omega*(n + 1)*n
                dddd(:, :, j, i) = f0*lo_outerproduct(v0, v1)/(uc%volume*lo_kb_hartree*temperature)
            end do
        end do
        dddd = dddd*lo_kappa_au_to_SI
        call h5%store_data(dddd, h5%file_id, trim(dname), enhet='W/mK', dimensions='q-vector,mode,xyz,xyz')
        call mem%deallocate(dddd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    else if (allocated(dr%iq(1)%linewidth) .and. present(temperature)) then
        do i = 1, qp%n_full_point
            l = qp%ap(i)%irreducible_index
            k = qp%ap(i)%operation_from_irreducible
            do j = 1, dr%n_mode
                if (dr%aq(i)%omega(j) .lt. dr%omega_min*0.5_r8) cycle
                if (k .gt. 0) then
                    v0 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%vel(:, j), reciprocal=.true.)
                else
                    v0 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%vel(:, j), reciprocal=.true.)
                end if
                ! Get kappa for this q-point
                omega = dr%iq(l)%omega(j)
                cv = lo_harmonic_oscillator_cv(temperature, omega)
                if (dr%iq(l)%linewidth(j) .gt. lo_freqtol) then
                    tau = 1.0_r8/(2.0_r8*dr%iq(l)%linewidth(j))
                    dddd(:, :, j, i) = lo_outerproduct(v0, v0)*cv*tau/uc%volume
                end if
            end do
        end do
        dddd = dddd*lo_kappa_au_to_SI
        call h5%store_data(dddd, h5%file_id, trim(dname), enhet='W/mK', dimensions='q-vector,mode,xyz,xyz')
        call mem%deallocate(dddd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! If a temperature is specified, some more stuff
    if (present(temperature)) then
        call h5%store_attribute(temperature, h5%file_id, 'temperature', lo_status)
        ! Might as well dump heat capacities while I'm here
        dname = 'harmonic_heat_capacity'
        call mem%allocate(dd, [dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_full_point
        do j = 1, dr%n_mode
            dd(j, i) = lo_harmonic_oscillator_cv(temperature, dr%aq(i)%omega(j))
        end do
        end do
        dd = dd*lo_Hartree_to_Joule
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='J/K', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> write the dispersion on the irreducible mesh to a hdf5 file
subroutine write_irreducible_to_hdf5(dr, qp, uc, filename, mem)
    !> phonon dispersions
    class(lo_phonon_dispersions), intent(in) :: dr
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> filename
    character(len=*), intent(in) :: filename
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:, :, :, :), allocatable :: dddd
    real(r8), dimension(:, :, :), allocatable :: ddd
    real(r8), dimension(:, :), allocatable :: dd
    real(r8), dimension(:), allocatable :: di
    real(r8), dimension(3) :: v0, v1
    real(r8) :: n1, n, f0, omega, cv, tau
    integer :: i, j, k, l
    character(len=1000) :: dname

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Some general attributes first:
    call h5%store_attribute(dr%n_mode/3, h5%file_id, 'number_of_atoms', lo_status)
    call h5%store_attribute(dr%n_mode, h5%file_id, 'number_of_bands', lo_status)
    call h5%store_attribute(dr%n_irr_qpoint, h5%file_id, 'number_of_qpoints', lo_status)
    select type (qp)
    type is (lo_monkhorst_pack_mesh)
        call h5%store_attribute('Monkhorst-Pack', h5%file_id, 'meshtype', lo_status)
    type is (lo_fft_mesh)
        call h5%store_attribute('FFT', h5%file_id, 'meshtype', lo_status)
    type is (lo_wedge_mesh)
        call h5%store_attribute('Wedge', h5%file_id, 'meshtype', lo_status)
    end select

    ! q-points
    dname = 'qpoints'
    call mem%allocate(dd, [3, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    do i = 1, qp%n_irr_point
        dd(:, i) = qp%ip(i)%r
    end do
    call h5%store_data(dd, h5%file_id, trim(dname), enhet='1/A', dimensions='q-vector,xyz')
    call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! frequencies
    if (allocated(dr%iq(1)%omega)) then
        dname = 'frequencies'
        call mem%allocate(dd, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_irr_point
            dd(:, i) = dr%iq(i)%omega*lo_frequency_hartree_to_Hz
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='Hz (angular)', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Group velocities
    if (allocated(dr%iq(1)%vel)) then
        dname = 'group_velocities'
        call mem%allocate(ddd, [3, dr%n_mode, qp%n_irr_point], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_irr_point
            ddd(:, :, i) = dr%iq(i)%vel*lo_groupvel_hartreebohr_to_ms
        end do
        call h5%store_data(ddd, h5%file_id, trim(dname), enhet='m/s', dimensions='q-vector,mode,xyz')
        call mem%deallocate(ddd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Gruneisen parameters
    if (allocated(dr%iq(1)%gruneisen)) then
        dname = 'gruneisen_parameters'
        call mem%allocate(dd, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do i = 1, qp%n_irr_point
            dd(:, i) = dr%iq(i)%gruneisen
        end do
        call h5%store_data(dd, h5%file_id, trim(dname), enhet='dimensionless', dimensions='q-vector,mode')
        call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! unit cell vectors
    dname = 'lattice_vectors'
    call mem%allocate(dd, [3, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ! TODO: check transposition convention for lattice vectors
    dd = lo_bohr_to_m*1.d10*transpose(uc%latticevectors)
    call h5%store_data(dd, h5%file_id, trim(dname), enhet='A', dimensions='xyz,xyz')
    call mem%deallocate(dd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! integration weights
    dname = 'integration_weights'
    call mem%allocate(di, qp%n_irr_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    do i = 1, qp%n_irr_point
        di(i) = qp%ip(i)%integration_weight
    end do
    call h5%store_data(di, h5%file_id, trim(dname), enhet='1', dimensions='q-vector')
    call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> calculate the phonon free energy as a direct sum
pure function phonon_free_energy(dr, temperature) result(f)
    !> the phonon dispersions
    class(lo_phonon_dispersions), intent(in) :: dr
    !> the temperature
    real(r8), intent(in) :: temperature
    !> the free energy
    real(r8) :: f

    integer :: i, j
    ! Return a stupid number if there are imaginary modes:
    if (dr%omega_min .lt. 0.0_r8) then
        f = 123456789.0_r8
        return
    end if
    ! If not, calculate it the normal way
    f = 0.0_r8
    do i = 1, dr%n_full_qpoint
        do j = 1, dr%n_mode
            f = f + lo_harmonic_oscillator_free_energy(temperature, dr%aq(i)%omega(j))
        end do
    end do
    f = f/dr%n_full_qpoint/(dr%n_mode/3)
end function

!> calculate the phonon free energy as a direct sum
pure function phonon_free_energy_classical(dr, temperature) result(f)
    !> the phonon dispersions
    class(lo_phonon_dispersions), intent(in) :: dr
    !> the temperature
    real(r8), intent(in) :: temperature
    !> the free energy
    real(r8) :: f

    integer :: i, j
    ! Return a stupid number if there are imaginary modes:
    if (dr%omega_min .lt. 0.0_r8) then
        f = 123456789.0_r8
        return
    end if
    ! If not, calculate it the normal way
    f = 0.0_r8
    do i = 1, dr%n_full_qpoint
        do j = 1, dr%n_mode
            f = f + lo_classical_harmonic_oscillator_free_energy(temperature, dr%aq(i)%omega(j))
        end do
    end do
    f = f/dr%n_full_qpoint/(dr%n_mode/3)
end function

!> calculate the phonon heat capacity as a direct sum
pure function phonon_cv(dr, temperature) result(cv)
    !> the phonon dispersions
    class(lo_phonon_dispersions), intent(in) :: dr
    !> the temperature
    real(r8), intent(in) :: temperature
    !> the heat capacity
    real(r8) :: cv

    integer :: i, j
    cv = 0.0_r8
    do i = 1, dr%n_full_qpoint
        do j = 1, dr%n_mode
            cv = cv + lo_harmonic_oscillator_cv(temperature, dr%aq(i)%omega(j))
        end do
    end do
    cv = cv/dr%n_full_qpoint/(dr%n_mode/3)
end function

!> calculate the phonon entropy as a direct sum
pure function phonon_entropy(dr, temperature, modenum, sitenum) result(s)
    !> the phonon dispersions
    class(lo_phonon_dispersions), intent(in) :: dr
    !> the temperature
    real(r8), intent(in) :: temperature
    !> the entropy in eV/K/atom (entropy for the modenum-th mode if modenum is specified)
    real(r8) :: s
    !> calculate vibrational entropy for a specific mode
    integer, intent(in), optional :: modenum
    !> project entropy onto a specific site
    integer, intent(in), optional :: sitenum

    complex(r8), dimension(3) :: cv0
    real(r8) :: w
    integer :: i, j

    ! Calculate the entropy of only one mode if modenum is specified
    if (present(modenum)) then
        j = modenum
        s = 0.0_r8
        do i = 1, dr%n_full_qpoint
            s = s + lo_harmonic_oscillator_entropy(temperature, dr%aq(i)%omega(j))
        end do
        s = s/dr%n_full_qpoint/(dr%n_mode/3)
        return
    end if

    ! If modenum or sitenum was not specified, calculate the total entropy
    if (present(sitenum)) then
        s = 0.0_r8
        do i = 1, dr%n_full_qpoint
            do j = 1, dr%n_mode
                cv0 = dr%aq(i)%egv((sitenum - 1)*3 + 1:sitenum*3, j)
                w = abs(dot_product(cv0, cv0))
                s = s + w*lo_harmonic_oscillator_entropy(temperature, dr%aq(i)%omega(j))
            end do
        end do
        s = s/dr%n_full_qpoint/(dr%n_mode/3)
        return
    end if

    ! If modenum or sitenum was not specified, calculate the total entropy
    s = 0.0_r8
    do i = 1, dr%n_full_qpoint
        do j = 1, dr%n_mode
            s = s + lo_harmonic_oscillator_entropy(temperature, dr%aq(i)%omega(j))
        end do
    end do
    s = s/dr%n_full_qpoint/(dr%n_mode/3)
end function

!> calculate phonon angular momentum matrix
subroutine phonon_angular_momentum_matrix(dr, qp, uc, temperature, alpha, mw)
    !> dispersion relations
    class(lo_phonon_dispersions), intent(in) :: dr
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> temperature
    real(r8), intent(in) :: temperature
    !> displacement covariance matrix
    real(r8), dimension(3, 3), intent(out) :: alpha
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    real(r8), parameter :: faketau = 10E-12_r8*lo_time_s_to_au ! lifetime of 10 ps
    complex(r8), dimension(3, 3) :: Mx, My, Mz
    logical :: havetau

    ! Some initial things
    init: block
        ! Representation of the angular momentum operator
        Mx = 0.0_r8
        My = 0.0_r8
        Mz = 0.0_r8
        Mx(2, 3) = 1
        Mx(3, 2) = -1
        My(1, 3) = -1
        My(3, 1) = 1
        Mz(1, 2) = 1
        Mz(2, 1) = -1
        Mx = -Mx*lo_imag
        My = -My*lo_imag
        Mz = -Mz*lo_imag

        ! Check if I have real linewidths or just use a fake one
        if (allocated(dr%iq(1)%linewidth)) then
            havetau = .true.
        else
            havetau = .false.
        end if
    end block init

    ! Calculate actual angular momentum thingy
    calc: block
        complex(r8), dimension(3) :: cv0, cv1, cv2
        real(r8), dimension(3) :: v0, v1, w0, w1
        real(r8) :: f0
        integer :: i, j, k, l, o

        alpha = 0.0_r8
        l = 0
        do i = 1, qp%n_full_point
        do j = 1, dr%n_mode
            l = l + 1
            if (mod(l, mw%n) .ne. mw%r) cycle
            ! Skip acoustic
            if (dr%aq(i)%omega(j) .lt. lo_freqtol) cycle
            ! Get the weird rotation guy
            cv0 = 0.0_r8
            do k = 1, uc%na
                cv1 = dr%aq(i)%egv((k - 1)*3 + 1:k*3, j)
                cv2 = matmul(Mx, cv1)
                cv0(1) = cv0(1) + dot_product(cv1, cv2)
                cv2 = matmul(My, cv1)
                cv0(2) = cv0(2) + dot_product(cv1, cv2)
                cv2 = matmul(Mz, cv1)
                cv0(3) = cv0(3) + dot_product(cv1, cv2)
            end do
            ! Now average over the small group. Seems sensible? Yes no maybe.
            v0 = 0.0_r8
            w0 = 0.0_r8
            v1 = real(cv0)
            w1 = dr%aq(i)%vel(:, j)
            do k = 1, qp%ap(i)%n_invariant_operation
                o = qp%ap(i)%invariant_operation(k)
                v0 = v0 + matmul(uc%sym%op(o)%m, v1)
                w0 = w0 + matmul(uc%sym%op(o)%m, w1)
            end do
            v0 = v0/real(qp%ap(i)%n_invariant_operation, r8)
            w0 = w0/real(qp%ap(i)%n_invariant_operation, r8)

            ! dn/dT
            f0 = lo_harmonic_oscillator_cv(temperature, dr%aq(i)%omega(j))/dr%aq(i)%omega(j)

            f0 = f0*qp%ap(i)%integration_weight
            if (havetau) then
                ! This is divided by tau
                k = qp%ap(i)%irreducible_index
                f0 = f0/(2.0_r8*dr%iq(k)%linewidth(j))
            else
                ! This is divided by random constant number
                f0 = f0*faketau
            end if
            alpha = alpha + lo_outerproduct(v0, w0)*f0
        end do
        end do
        call mw%allreduce('sum', alpha)
        ! And finally scale with volume. I should really do a
        ! dimensionality analysis on this to figure out the unit.
        alpha = alpha/uc%volume
        f0 = norm2(alpha)
        alpha = lo_chop(alpha, f0*1E-10_r8)
    end block calc
end subroutine

!> Calculate the thermal displacement covariance matrix
subroutine thermal_displacement_matrix(dr, qp, uc, temperature, sigma, mw, mem)
    !> dispersion relations
    class(lo_phonon_dispersions), intent(in) :: dr
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> temperature
    real(r8), intent(in) :: temperature
    !> displacement covariance matrix
    real(r8), dimension(:, :, :), intent(out) :: sigma
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    complex(r8), dimension(:, :, :), allocatable :: cs
    complex(r8), dimension(3) :: cv0
    real(r8), dimension(:, :), allocatable :: coeffM, rotM, invarM, IM
    real(r8), dimension(:), allocatable :: wV, vX
    real(r8), dimension(9, 9) :: m9, T
    real(r8) :: f0
    integer :: q, j, a1, a2, ctr, ne, nx, o

    ! Some parameters and temporary space
    ne = uc%na*9
    call mem%allocate(cs, [3, 3, uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(rotM, [ne, ne], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(invarM, [ne, ne], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(IM, [ne, ne], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(wV, ne, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    cs = 0.0_r8
    rotM = 0.0_r8
    invarM = 0.0_r8
    IM = 0.0_r8
    wV = 0.0_r8

    ! Add up the thermal displacement matrix
    sigma = 0.0_r8
    cs = 0.0_r8
    ctr = 0
    do q = 1, qp%n_full_point
    do j = 1, dr%n_mode
        ctr = ctr + 1
        if (mod(ctr, mw%n) .ne. mw%r) cycle
        if (dr%aq(q)%omega(j) .lt. lo_freqtol) cycle
        f0 = 0.5_r8*(1.0_r8 + 2.0_r8*lo_planck(temperature, dr%aq(q)%omega(j)))/dr%aq(q)%omega(j)
        do a1 = 1, uc%na
            cv0 = dr%aq(q)%egv((a1 - 1)*3 + 1:a1*3, j)
            cs(:, :, a1) = cs(:, :, a1) + lo_outerproduct(cv0, cv0)*f0*qp%ap(q)%integration_weight
        end do
    end do
    end do
    call mw%allreduce('sum', cs)

    ! Now make sure the symmetry is proper.
    T(:, 1) = [1, 0, 0, 0, 0, 0, 0, 0, 0]
    T(:, 2) = [0, 0, 0, 1, 0, 0, 0, 0, 0]
    T(:, 3) = [0, 0, 0, 0, 0, 0, 1, 0, 0]
    T(:, 4) = [0, 1, 0, 0, 0, 0, 0, 0, 0]
    T(:, 5) = [0, 0, 0, 0, 1, 0, 0, 0, 0]
    T(:, 6) = [0, 0, 0, 0, 0, 0, 0, 1, 0]
    T(:, 7) = [0, 0, 1, 0, 0, 0, 0, 0, 0]
    T(:, 8) = [0, 0, 0, 0, 0, 1, 0, 0, 0]
    T(:, 9) = [0, 0, 0, 0, 0, 0, 0, 0, 1]

    invarm = 0.0_r8
    rotM = 0.0_r8
    call lo_identitymatrix(IM)
    do o = 1, uc%sym%n
        ! Expand the operation
        m9 = lo_expandoperation_pair(uc%sym%op(o)%m)
        rotM = 0.0_r8
        do a1 = 1, uc%na
            a2 = uc%sym%op(o)%fmap(a1)
            rotM((a2 - 1)*9 + 1:a2*9, (a1 - 1)*9 + 1:a1*9) = m9
        end do
        invarM = invarM + rotM - IM
    end do
    ! And the transposition
    rotM = 0.0_r8
    do a1 = 1, uc%na
        rotM((a1 - 1)*9 + 1:a1*9, (a1 - 1)*9 + 1:a1*9) = T
    end do
    invarM = invarM + rotM - IM
    call lo_nullspace_coefficient_matrix(invarM, coeffM, nx)
    call mem%allocate(vX, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    vX = 0.0_r8

    do a1 = 1, uc%na
        wV((a1 - 1)*9 + 1:a1*9) = lo_flattentensor(real(cs(:, :, a1)))/uc%mass(a1)
    end do
    call lo_linear_least_squares(coeffM, wV, vX)
    call lo_gemv(coeffM, vX, wV)
    do a1 = 1, uc%na
        sigma(:, :, a1) = lo_unflatten_2tensor(wV((a1 - 1)*9 + 1:a1*9))
    end do

    ! And some cleanup
    call mem%deallocate(cs, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(rotM, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(invarM, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(IM, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(wV, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(vX, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Calculate all the harmonic things for a q-point
subroutine harmonic_things_for_single_point(ompoint, fc, p, mem, qpoint, qvec, qdirection)
    !> point in the dispersions
    class(lo_phonon_dispersions_qpoint), intent(inout) :: ompoint
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> q-point
    class(lo_qpoint), intent(in), optional :: qpoint
    !> q-point as just a vector
    real(r8), dimension(3), intent(in), optional :: qvec
    !> q-direction?
    real(r8), dimension(3), intent(in), optional :: qdirection

    real(r8), dimension(3), parameter :: default_direction = [1.0_r8, 0.0_r8, 0.0_r8]
    type(lo_qpoint) :: dumqpoint
    real(r8), dimension(3) :: qdir
    integer :: nb
    logical :: skipna

    ! So this is slightly too complicated for something this simple, one would think.
    ! But it's really not that bad. Just makes sense to make the thing flexible.
    ! Sort out some simple things first:
    init: block
        ! Really basic sanity tests. Will make more when I think of it.
        if (p%na .ne. fc%na) then
            call lo_stop_gracefully(['Different number of atoms in structure and forceconstant.'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if

        ! Get the number of modes
        nb = fc%na*3

        ! Check if there is space, if not, make some
        if (.not. allocated(ompoint%omega)) allocate (ompoint%omega(nb))
        if (.not. allocated(ompoint%egv)) allocate (ompoint%egv(nb, nb))
        if (.not. allocated(ompoint%vel)) allocate (ompoint%vel(3, nb))
        if (.not. allocated(ompoint%degeneracy)) allocate (ompoint%degeneracy(nb))
        if (.not. allocated(ompoint%degenmode)) allocate (ompoint%degenmode(nb, nb))
        ! if ( .not. allocated(ompoint%qhessian)   ) allocate(ompoint%qhessian(3,3,nb))
        ompoint%omega = 0.0_r8
        ompoint%egv = 0.0_r8
        ompoint%vel = 0.0_r8
        ompoint%degeneracy = 0
        ompoint%degenmode = 0
        ! ompoint%qhessian=0.0_r8

        ! Figure out what to do with the q-point
        if (present(qpoint) .and. present(qvec)) then
            call lo_stop_gracefully(['Either provide a q-point, or a q-vector, not both'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if

        ! Is there a direction specified?
        if (present(qdirection)) then
            if (lo_sqnorm(qdirection) .gt. lo_sqtol) then
                qdir = qdirection
                skipna = .false.
            else
                qdir = 0.0_r8
                skipna = .true.
            end if
        else
            qdir = default_direction
            skipna = .false.
        end if
        qdir = qdir/norm2(qdir)

        ! If there is no q-point supplied, construct one from the input vector
        if (present(qvec)) then
            dumqpoint%r = qvec - p%bz%gshift(qvec + lo_degenvector)
            call lo_get_small_group_of_qpoint(dumqpoint, p)
        end if
    end block init

    ! Solve for harmonic things
    slv: block
        ! complex(r8), dimension(:,:,:), allocatable :: Dqq
        complex(r8), dimension(:, :, :), allocatable :: Dq
        complex(r8), dimension(:, :), allocatable :: D
        real(r8), dimension(:, :), allocatable :: coeffM, pinv_coeffM
        real(r8), dimension(3, 3) :: m0, I3
        integer, dimension(:, :), allocatable :: di
        integer :: i, l, b1, b2, iop, nvar

        ! Get the frequencies
        ! call mem%allocate(Dqq,[nb,nb,6],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(Dq, [nb, nb, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(D, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(di, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        !Dqq=0.0_r8
        Dq = 0.0_r8
        D = 0.0_r8
        di = 0
        if (present(qpoint)) then
            call fc%dynamicalmatrix( &
                p, qpoint, D, mem, Dq, qdirection=qdir, skipnonanalytical=skipna)
            call fc%frequencies_eigenvectors_groupvelocities( &
                D, ompoint%omega, mem, Dq, eigenvectors=ompoint%egv, groupvelocities=ompoint%vel, qpoint=qpoint)
        else
            call fc%dynamicalmatrix( &
                p, dumqpoint, D, mem, Dq, qdirection=qdir, skipnonanalytical=skipna)
            call fc%frequencies_eigenvectors_groupvelocities( &
                D, ompoint%omega, mem, Dq, eigenvectors=ompoint%egv, groupvelocities=ompoint%vel, qpoint=dumqpoint)
        end if
        ! Sort out degeneracies
        di = 0
        do b1 = 1, nb
        do b2 = b1, nb
            if (abs(ompoint%omega(b1) - ompoint%omega(b2)) .lt. lo_freqtol) then
                di(b1, b2) = 1
                di(b2, b1) = 1
            end if
        end do
        end do
        ! store
        do b1 = 1, nb
            ompoint%degeneracy(b1) = sum(di(:, b1))
        end do
        ! some space to store the degeneracy
        ompoint%degenmode = 0
        do b1 = 1, nb
            l = 0
            do b2 = 1, nb
                if (di(b1, b2) .eq. 1) then
                    l = l + 1
                    ompoint%degenmode(l, b1) = b2
                end if
            end do
        end do

        ! Then sort out how group velocities should transform, perhaps?
        I3 = 0.0_r8
        do i = 1, 3
            I3(i, i) = 1.0_r8
        end do
        m0 = 0.0_r8
        if (present(qvec)) then
            do i = 1, dumqpoint%n_invariant_operation
                iop = dumqpoint%invariant_operation(i)
                if (iop .gt. 0) then
                    m0 = m0 + p%sym%op(iop)%m - I3
                else
                    m0 = m0 - p%sym%op(iop)%m - I3
                end if
            end do
        else
            do i = 1, qpoint%n_invariant_operation
                iop = qpoint%invariant_operation(i)
                if (iop .lt. 0) then
                    m0 = m0 + p%sym%op(iop)%m - I3
                else
                    m0 = m0 - p%sym%op(iop)%m - I3
                end if
            end do
        end if
        ! Get the invariance thing?
        if (norm2(m0) .gt. lo_sqtol) then
            call lo_nullspace_coefficient_matrix(m0, coeffM, nvar)
        else
            nvar = 0
        end if
        if (nvar .gt. 0) then
            allocate (pinv_coeffM(nvar, 3))
            call lo_real_pseudoinverse(coeffM, pinv_coeffM)
            ! Get the magic invariance matrix thing
            ompoint%invariance_fixer = lo_chop(matmul(coeffM, pinv_coeffM), 1E-12_r8)
            deallocate (coeffM)
            deallocate (pinv_coeffM)
        else
            ompoint%invariance_fixer = I3
        end if

        ! Cleanup
        ! call mem%deallocate(Dqq,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(Dq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(D, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (present(qvec)) then
            deallocate (dumqpoint%invariant_operation)
        end if
    end block slv
end subroutine

!> measure size in memory, in bytes
function phonon_dispersions_size_in_mem(dr) result(mem)
    !> dispersions
    class(lo_phonon_dispersions), intent(in) :: dr
    !> memory in bytes
    integer :: mem, i

    mem = 0
    ! easy things
    mem = mem + storage_size(dr)
    if (allocated(dr%default_smearing)) mem = mem + storage_size(dr%default_smearing)*size(dr%default_smearing)
    mem = mem/8
    ! more annoying things
    if (allocated(dr%iq)) then
        do i = 1, size(dr%iq)
            mem = mem + dr%iq(i)%size_in_mem()
        end do
    end if
    if (allocated(dr%aq)) then
        do i = 1, size(dr%aq)
            mem = mem + dr%aq(i)%size_in_mem()
        end do
    end if
end function

!> measure size in memory, in bytes
function phonon_dispersions_qpoint_size_in_mem(p) result(mem)
    !> one qpoint from the dispersions
    class(lo_phonon_dispersions_qpoint), intent(in) :: p
    !> memory in bytes
    integer :: mem

    mem = 0
    mem = mem + storage_size(p)
    if (allocated(p%omega)) mem = mem + storage_size(p%omega)*size(p%omega)
    if (allocated(p%vel)) mem = mem + storage_size(p%vel)*size(p%vel)
    if (allocated(p%egv)) mem = mem + storage_size(p%egv)*size(p%egv)
    if (allocated(p%gruneisen)) mem = mem + storage_size(p%gruneisen)*size(p%gruneisen)
    if (allocated(p%linewidth)) mem = mem + storage_size(p%linewidth)*size(p%linewidth)
    if (allocated(p%shift3)) mem = mem + storage_size(p%shift3)*size(p%shift3)
    if (allocated(p%shift4)) mem = mem + storage_size(p%shift4)*size(p%shift4)
    if (allocated(p%p_plus)) mem = mem + storage_size(p%p_plus)*size(p%p_plus)
    if (allocated(p%p_minus)) mem = mem + storage_size(p%p_minus)*size(p%p_minus)
    if (allocated(p%p_iso)) mem = mem + storage_size(p%p_iso)*size(p%p_iso)
    if (allocated(p%qs)) mem = mem + storage_size(p%qs)*size(p%qs)
    if (allocated(p%F0)) mem = mem + storage_size(p%F0)*size(p%F0)
    if (allocated(p%Fn)) mem = mem + storage_size(p%Fn)*size(p%Fn)
    if (allocated(p%kappa)) mem = mem + storage_size(p%kappa)*size(p%kappa)
    if (allocated(p%F)) mem = mem + storage_size(p%F)*size(p%F)
    if (allocated(p%S)) mem = mem + storage_size(p%S)*size(p%S)
    if (allocated(p%deltaF3)) mem = mem + storage_size(p%deltaF3)*size(p%deltaF3)
    if (allocated(p%deltaF4)) mem = mem + storage_size(p%deltaF4)*size(p%deltaF4)
    if (allocated(p%deltaS3)) mem = mem + storage_size(p%deltaS3)*size(p%deltaS3)
    if (allocated(p%deltaS4)) mem = mem + storage_size(p%deltaS4)*size(p%deltaS4)
    if (allocated(p%degeneracy)) mem = mem + storage_size(p%degeneracy)*size(p%degeneracy)
    if (allocated(p%degenmode)) mem = mem + storage_size(p%degenmode)*size(p%degenmode)
    if (allocated(p%pyroelectric_amplitude)) mem = mem + storage_size(p%pyroelectric_amplitude)*size(p%pyroelectric_amplitude)
    mem = mem/8
end function

!> measure size once packed into a character buffer.
pure function phonon_dispersions_qpoint_size_packed(p) result(mem)
    !> one qpoint from the dispersions
    class(lo_phonon_dispersions_qpoint), intent(in) :: p
    !> memory in bytes
    integer :: mem

    logical, dimension(23) :: fld
    fld = .false.

    mem = 0
    mem = mem + storage_size(fld)*size(fld)/8
    if (allocated(p%omega)) mem = mem + storage_size(p%omega)*size(p%omega)/8
    if (allocated(p%vel)) mem = mem + storage_size(p%vel)*size(p%vel)/8
    if (allocated(p%egv)) mem = mem + storage_size(p%egv)*size(p%egv)/8
    if (allocated(p%gruneisen)) mem = mem + storage_size(p%gruneisen)*size(p%gruneisen)/8
    if (allocated(p%linewidth)) mem = mem + storage_size(p%linewidth)*size(p%linewidth)/8
    if (allocated(p%shift3)) mem = mem + storage_size(p%shift3)*size(p%shift3)/8
    if (allocated(p%shift4)) mem = mem + storage_size(p%shift4)*size(p%shift4)/8
    if (allocated(p%p_plus)) mem = mem + storage_size(p%p_plus)*size(p%p_plus)/8
    if (allocated(p%p_minus)) mem = mem + storage_size(p%p_minus)*size(p%p_minus)/8
    if (allocated(p%p_iso)) mem = mem + storage_size(p%p_iso)*size(p%p_iso)/8
    if (allocated(p%qs)) mem = mem + storage_size(p%qs)*size(p%qs)/8
    if (allocated(p%F0)) mem = mem + storage_size(p%F0)*size(p%F0)/8
    if (allocated(p%Fn)) mem = mem + storage_size(p%Fn)*size(p%Fn)/8
    if (allocated(p%kappa)) mem = mem + storage_size(p%kappa)*size(p%kappa)/8
    if (allocated(p%F)) mem = mem + storage_size(p%F)*size(p%F)/8
    if (allocated(p%S)) mem = mem + storage_size(p%S)*size(p%S)/8
    if (allocated(p%deltaF3)) mem = mem + storage_size(p%deltaF3)*size(p%deltaF3)/8
    if (allocated(p%deltaF4)) mem = mem + storage_size(p%deltaF4)*size(p%deltaF4)/8
    if (allocated(p%deltaS3)) mem = mem + storage_size(p%deltaS3)*size(p%deltaS3)/8
    if (allocated(p%deltaS4)) mem = mem + storage_size(p%deltaS4)*size(p%deltaS4)/8
    if (allocated(p%degeneracy)) mem = mem + storage_size(p%degeneracy)*size(p%degeneracy)/8
    if (allocated(p%degenmode)) mem = mem + storage_size(p%degenmode)*size(p%degenmode)/8
    if (allocated(p%pyroelectric_amplitude)) mem = mem + storage_size(p%pyroelectric_amplitude)*size(p%pyroelectric_amplitude)/8
end function

!> pack a phonon dispersion q-point to mpi character buffer
subroutine pack_phonon_dispersions_qpoint(p, buf, pos, mw)
    !> one qpoint from the dispersions
    class(lo_phonon_dispersions_qpoint), intent(in) :: p
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    logical, dimension(23) :: fld

    ! Figure out the relevant fields?
    fld = .false.
    if (allocated(p%omega)) fld(1) = .true.
    if (allocated(p%vel)) fld(2) = .true.
    if (allocated(p%egv)) fld(3) = .true.
    if (allocated(p%gruneisen)) fld(4) = .true.
    if (allocated(p%linewidth)) fld(5) = .true.
    if (allocated(p%shift3)) fld(6) = .true.
    if (allocated(p%shift4)) fld(7) = .true.
    if (allocated(p%p_plus)) fld(8) = .true.
    if (allocated(p%p_minus)) fld(9) = .true.
    if (allocated(p%p_iso)) fld(10) = .true.
    if (allocated(p%qs)) fld(11) = .true.
    if (allocated(p%F0)) fld(12) = .true.
    if (allocated(p%Fn)) fld(13) = .true.
    if (allocated(p%kappa)) fld(14) = .true.
    if (allocated(p%F)) fld(15) = .true.
    if (allocated(p%S)) fld(16) = .true.
    if (allocated(p%deltaF3)) fld(17) = .true.
    if (allocated(p%deltaF4)) fld(18) = .true.
    if (allocated(p%deltaS3)) fld(19) = .true.
    if (allocated(p%deltaS4)) fld(20) = .true.
    if (allocated(p%degeneracy)) fld(21) = .true.
    if (allocated(p%degenmode)) fld(22) = .true.
    if (allocated(p%pyroelectric_amplitude)) fld(23) = .true.

    ! Pack the list of fields
    call mw%pack(fld, buf, pos)

    ! Now start packing everything else.
    if (fld(1)) call mw%pack(p%omega, buf, pos)
    if (fld(2)) call mw%pack(p%vel, buf, pos)
    if (fld(3)) call mw%pack(p%egv, buf, pos)
    if (fld(4)) call mw%pack(p%gruneisen, buf, pos)
    if (fld(5)) call mw%pack(p%linewidth, buf, pos)
    if (fld(6)) call mw%pack(p%shift3, buf, pos)
    if (fld(7)) call mw%pack(p%shift4, buf, pos)
    if (fld(8)) call mw%pack(p%p_plus, buf, pos)
    if (fld(9)) call mw%pack(p%p_minus, buf, pos)
    if (fld(10)) call mw%pack(p%p_iso, buf, pos)
    if (fld(11)) call mw%pack(p%qs, buf, pos)
    if (fld(12)) call mw%pack(p%F0, buf, pos)
    if (fld(13)) call mw%pack(p%Fn, buf, pos)
    if (fld(14)) call mw%pack(p%kappa, buf, pos)
    if (fld(15)) call mw%pack(p%F, buf, pos)
    if (fld(16)) call mw%pack(p%S, buf, pos)
    if (fld(17)) call mw%pack(p%deltaF3, buf, pos)
    if (fld(18)) call mw%pack(p%deltaF4, buf, pos)
    if (fld(19)) call mw%pack(p%deltaS3, buf, pos)
    if (fld(20)) call mw%pack(p%deltaS4, buf, pos)
    if (fld(21)) call mw%pack(p%degeneracy, buf, pos)
    if (fld(22)) call mw%pack(p%degenmode, buf, pos)
    if (fld(23)) call mw%pack(p%pyroelectric_amplitude, buf, pos)
end subroutine

!> pack a phonon dispersion q-point to mpi character buffer
subroutine unpack_phonon_dispersions_qpoint(p, buf, pos, mw)
    !> one qpoint from the dispersions
    class(lo_phonon_dispersions_qpoint), intent(inout) :: p
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    logical, dimension(23) :: fld
    ! Unpack the list of fields
    call mw%unpack(fld, buf, pos)
    ! Then unpack the relevant fields
    if (fld(1)) call mw%unpack(p%omega, buf, pos)
    if (fld(2)) call mw%unpack(p%vel, buf, pos)
    if (fld(3)) call mw%unpack(p%egv, buf, pos)
    if (fld(4)) call mw%unpack(p%gruneisen, buf, pos)
    if (fld(5)) call mw%unpack(p%linewidth, buf, pos)
    if (fld(6)) call mw%unpack(p%shift3, buf, pos)
    if (fld(7)) call mw%unpack(p%shift4, buf, pos)
    if (fld(8)) call mw%unpack(p%p_plus, buf, pos)
    if (fld(9)) call mw%unpack(p%p_minus, buf, pos)
    if (fld(10)) call mw%unpack(p%p_iso, buf, pos)
    if (fld(11)) call mw%unpack(p%qs, buf, pos)
    if (fld(12)) call mw%unpack(p%F0, buf, pos)
    if (fld(13)) call mw%unpack(p%Fn, buf, pos)
    if (fld(14)) call mw%unpack(p%kappa, buf, pos)
    if (fld(15)) call mw%unpack(p%F, buf, pos)
    if (fld(16)) call mw%unpack(p%S, buf, pos)
    if (fld(17)) call mw%unpack(p%deltaF3, buf, pos)
    if (fld(18)) call mw%unpack(p%deltaF4, buf, pos)
    if (fld(19)) call mw%unpack(p%deltaS3, buf, pos)
    if (fld(20)) call mw%unpack(p%deltaS4, buf, pos)
    if (fld(21)) call mw%unpack(p%degeneracy, buf, pos)
    if (fld(22)) call mw%unpack(p%degenmode, buf, pos)
    if (fld(23)) call mw%unpack(p%pyroelectric_amplitude, buf, pos)
end subroutine

!> destroy
subroutine destroy(dr)
    !> dispersions
    class(lo_phonon_dispersions), intent(inout) :: dr

    if (allocated(dr%iq)) deallocate (dr%iq)
    if (allocated(dr%aq)) deallocate (dr%aq)
    if (allocated(dr%default_smearing)) deallocate (dr%default_smearing)
    dr%n_mode = -lo_hugeint
    dr%n_irr_qpoint = -lo_hugeint
    dr%n_full_qpoint = -lo_hugeint
    dr%omega_max = -lo_huge
    dr%omega_min = -lo_huge
end subroutine

end module
