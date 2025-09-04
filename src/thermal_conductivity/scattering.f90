#include "precompilerdefinitions"
module scattering
use konstanter, only: r8, i8, lo_freqtol, lo_twopi, lo_exitcode_param, lo_hugeint, lo_pi, lo_tol, &
                      lo_phonongroupveltol, lo_tol, lo_frequency_THz_to_Hartree, lo_kb_hartree, lo_huge, &
                      lo_frequency_Hartree_to_meV
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_gauss, lo_planck, lo_return_unique
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_symmetryoperation, only: lo_operate_on_vector
use lo_randomnumbers, only: lo_mersennetwister
use lo_fftgrid_helper, only: lo_montecarlo_grid, singlet_to_triplet, triplet_to_singlet, &
                             fft_third_grid_index, fft_fourth_grid_index
use lo_timetracker, only: lo_timer

use options, only: lo_opts

use lo_sorting, only: lo_qsort

implicit none
private
public :: lo_scattering_rates

real(r8), parameter :: isotope_prefactor = lo_pi/4.0_r8
real(r8), parameter :: threephonon_prefactor = lo_pi/16.0_r8
real(r8), parameter :: fourphonon_prefactor = lo_pi/96.0_r8

! Container for scattering rates
type lo_scattering_rates
    !> The number of qpoint/mode on this rank
    integer :: my_nqpoints
    !> The list of qpoint and modes for this rank
    integer, dimension(:), allocatable :: my_qpoints, my_modes
    !> Bose-Einstein and squared smearing for each mode on irreducible q-point
    real(r8), dimension(:, :), allocatable :: be, sigsq
    !> The scattering matrix
    real(r8), dimension(:, :), allocatable :: Xi
    !> The reciprocal lattice vectors scaled by the qgrid, for adaptive broadening
    real(r8), dimension(3, 3) :: reclat = -lo_huge
    !> The smearing parameter
    real(r8) :: sigma, thresh_sigma

contains
    !> Generate the scattering amplitudes
    procedure :: generate
    !> destroy the scattering amplitues
    procedure :: destroy => sr_destroy
    !> Measure size in memory, in bytes
    procedure :: size_in_mem => sr_size_in_mem
end type

contains
subroutine generate(sr, qp, dr, uc, fct, fcf, opts, tmr, mw, mem)
    !> The scattering rate
    class(lo_scattering_rates), intent(out) :: sr
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The third order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> The fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> The options
    type(lo_opts), intent(in) :: opts
    !> Timer
    type(lo_timer), intent(inout) :: tmr
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> The random number generator
    type(lo_mersennetwister) :: rng
    !> The grids for monte-carlo integration
    type(lo_montecarlo_grid) :: mcg3, mcg4

    init: block
        !> The average broadening factor for each branch
        real(r8), dimension(:), allocatable :: sigavg
        !> A buffer for a vector
        real(r8), dimension(3) :: v0
        !> The calculated sigma for a mode
        real(r8) :: sigma
        !> To initialize the random number generator and timing
        real(r8) :: rseed
        !> The q-point grid dimension
        integer, dimension(3) :: dims
        !> Some integers for the do loop/indices
        integer :: q1, q2, b1, il, j, my_nqpoints, ctr

        ! grid dimensions
        select type (qp)
        class is (lo_fft_mesh)
            dims = qp%griddensity
        class default
            call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
        end select

        if (opts%seed .gt. 0) then
            rseed = 1.0 / real(opts%seed, r8)
        else
            rseed = walltime()
            if (mw%talk) write(*, *) '... walltime() used to generate random state'
        end if

        ! Initialize the random number generator
        call rng%init(iseed=mw%r, rseed=rseed)

        if (mw%talk) write (*, *) '... creating Monte-Carlo grid'
        ! Initialize the monte-carlo grid
        if (opts%thirdorder) then
            call mcg3%initialize(dims, opts%qg3ph)
        end if
        if (opts%fourthorder) then
            call mcg4%initialize(dims, opts%qg4ph)
        end if

        ! We can start some precomputation
        allocate (sr%be(qp%n_irr_point, dr%n_mode))
        allocate (sr%sigsq(qp%n_irr_point, dr%n_mode))
        call mem%allocate(sigavg, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        sigavg = 0.0_r8

        sr%be = 0.0_r8
        sr%sigsq = 0.0_r8
        ! This could be useful if integrationtype=1 -> fixed gaussian smearing
        sr%sigma = opts%sigma * lo_frequency_THz_to_Hartree
        ! This could be useful if integrationtype=6 -> Adapative smearing from groupvel diff
        do j=1, 3
            sr%reclat(:, j) = uc%reciprocal_latticevectors(:, j) / dims(j)
        end do

        ! And we get the values in the array
        do q1 = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                ! Skip gamma point
                if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle

                if (opts%classical) then
                    sr%be(q1, b1) = lo_kb_hartree*opts%temperature/dr%iq(q1)%omega(b1)
                else
                    sr%be(q1, b1) = lo_planck(opts%temperature, dr%iq(q1)%omega(b1))
                end if

                v0 = matmul(abs(dr%iq(q1)%vel(:, b1)), sr%reclat)**2
                ! This allows to work around problems with 2D materials
                sigma = maxval(v0*0.5_r8)
                sr%sigsq(q1, b1) = sigma

                ! Let's accumulate the average broadening factor for each mode
                sigavg(b1) = sigavg(b1) + sigma * qp%ip(q1)%integration_weight
            end do
        end do

        ! We add a little baseline, to avoid pathological case when group velocities are near zero
        sr%sigsq = sr%sigsq + maxval(sigavg) * 1e-4_r8
        ! This is to have a sanity check in memory, for integrationtype 6
        sr%thresh_sigma = maxval(sigavg) * 1e-4_r8
        call mem%deallocate(sigavg, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        if (mw%talk) write (*, *) '... distributing q-point/modes on MPI ranks'
        ctr = 0
        my_nqpoints = 0
        do q1 = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                ! We skip the acoustic mode at Gamma
                if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
                ctr = ctr + 1

                ! MPI thing
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                my_nqpoints = my_nqpoints + 1
            end do
        end do

        ! We can allocate all we need
        sr%my_nqpoints = my_nqpoints
        allocate (sr%my_qpoints(my_nqpoints))
        allocate (sr%my_modes(my_nqpoints))
        allocate (sr%Xi(my_nqpoints, qp%n_full_point*dr%n_mode))
        sr%my_qpoints = -lo_hugeint
        sr%my_modes = -lo_hugeint
        sr%Xi = 0.0_r8

        ! Let's attribute the q1/b1 indices to the ranks
        il = 0
        ctr = 0
        do q1 = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                ! We skip the acoustic mode at Gamma
                if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
                ctr = ctr + 1

                ! MPI thing
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                il = il + 1
                sr%my_qpoints(il) = q1
                sr%my_modes(il) = b1
            end do
        end do
        if (mw%talk) write (*, *) '... everything is ready, starting scattering computation'
        call tmr%tock('initialization')
    end block init

    scatt: block
        !> Buffer to contains the linewidth
        real(r8), dimension(:, :), allocatable :: buf_lw
        !> Buffer for the linewidth of the local point
        real(r8) :: buf, f0, velnorm, t0
        !> Some integers for the loops
        integer :: il, j, q1, b1, b2

        call mem%allocate(buf_lw, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_lw = 0.0_r8

        t0 = walltime()
        if (mw%talk) call lo_progressbar_init()
        do il = 1, sr%my_nqpoints
            buf = 0.0_r8
            if (opts%isotopescattering) then
                call compute_isotope_scattering(il, sr, qp, dr, uc, opts%temperature, buf, &
                                                opts%integrationtype, opts%sigma, mw, mem)
                call tmr%tock('isotope scattering')
            end if
            if (opts%thirdorder) then
                call compute_threephonon_scattering(il, sr, qp, dr, uc, fct, mcg3, rng, &
                                                    buf, opts%integrationtype, opts%sigma, mw, mem)
                call tmr%tock('threephonon scattering')
            end if
            if (opts%fourthorder) then
                call compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, mcg4, rng, &
                                                   buf, opts%integrationtype, opts%sigma, mw, mem)
                call tmr%tock('fourphonon scattering')
            end if
            ! We end with the boundary scattering
            if (opts%mfp_max .gt. 0.0_r8) then
                velnorm = norm2(dr%iq(sr%my_qpoints(il))%vel(:, sr%my_modes(il)))
                if (velnorm .gt. lo_phonongroupveltol) then
                    buf = buf + velnorm/opts%mfp_max
                end if
            end if
            ! Now we can update the linewidth for this mode
            buf_lw(sr%my_qpoints(il), sr%my_modes(il)) = buf

            if (mw%talk .and. lo_trueNtimes(il, 127, sr%my_nqpoints)) then
                call lo_progressbar(' ... computing scattering amplitude', il, sr%my_nqpoints, walltime() - t0)
            end if
        end do
        if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', sr%my_nqpoints, sr%my_nqpoints, walltime() - t0)

        ! Reduce the linewidth
        call mw%allreduce('sum', buf_lw)

        ! Distribute it after fixing the degeneracies
        do q1 = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
                ! First we fix the degeneracy
                f0 = 0.0_r8
                do j = 1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    f0 = f0 + buf_lw(q1, b2)
                end do
                f0 = f0/real(dr%iq(q1)%degeneracy(b1), r8)
                do j = 1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    buf_lw(q1, b2) = f0
                end do
                ! Now we can set the linewidth
                dr%iq(q1)%linewidth(b1) = buf_lw(q1, b1)

                ! While we are at it, we can set other things
                dr%iq(q1)%qs(b1) = 2.0_r8*dr%iq(q1)%linewidth(b1)
                velnorm = norm2(dr%iq(q1)%vel(:, b1))
!               if (velnorm .gt. lo_phonongroupveltol) then
                if (dr%iq(q1)%linewidth(b1) .gt. lo_freqtol) then
                    dr%iq(q1)%mfp(:, b1) = dr%iq(q1)%vel(:, b1)/dr%iq(q1)%qs(b1)
                    dr%iq(q1)%scalar_mfp(b1) = velnorm/dr%iq(q1)%qs(b1)
                    dr%iq(q1)%F0(:, b1) = dr%iq(q1)%mfp(:, b1)
                    dr%iq(q1)%Fn(:, b1) = dr%iq(q1)%F0(:, b1)
                end if
            end do
        end do
        call mem%deallocate(buf_lw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block scatt

    symmetrize_scatmat: block
        !> To hold the q-point in reduced coordinates
        real(r8), dimension(3) :: qv2, qv2p
        !> To get the index of the new triplet on the fft_grid
        integer, dimension(3) :: gi
        !> Some buffer
        real(r8), dimension(:), allocatable :: buf
        !> Integer for do loops and so on
        integer :: il, jl, q1, j, k, q2, b2, q2p, n

        call tmr%tick()
        if (mw%talk) write (*, *) '... symmetrizing scattering matrix'

        ! We use the relation Xi_{R*q, R*q'} = Xi_{q, q'''} to enforce the symmetry of Xi
        ! TODO look if these irreducible pair could reduce the cost
        call mem%allocate(buf, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do il = 1, sr%my_nqpoints
            q1 = sr%my_qpoints(il)
            allq2: do q2 = 1, qp%n_full_point
                buf = 0.0_r8
                n = 0
                do j = 1, qp%ip(q1)%n_invariant_operation
                    k = qp%ip(q1)%invariant_operation(j)
                    ! Here, we generate q''=R*q from the symmetry that leaves q invariant
                    select type (qp); type is (lo_fft_mesh)
                        qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
                        qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
                        if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle
                        gi = qp%index_from_coordinate(qv2p)
                        q2p = qp%gridind2ind(gi(1), gi(2), gi(3)) ! this is R*q'
                    end select
                    if (q2p .lt. q2) cycle allq2  ! If q2p < q2, we already did this guy
                    n = n + 1
                    ! We accumulate the values for each bands
                    do b2 = 1, dr%n_mode
                        jl = (q2p - 1)*dr%n_mode + b2
                        buf(b2) = buf(b2) + sr%Xi(il, jl)
                    end do
                end do
                if (n .eq. 0) cycle
                buf = buf/real(n, r8)
                ! And now we distribute
                do j = 1, qp%ip(q1)%n_invariant_operation
                    k = qp%ip(q1)%invariant_operation(j)
                    ! Here, we generate q''=R*q from the symmetry that leaves q invariant
                    select type (qp); type is (lo_fft_mesh)
                        qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
                        qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
                        if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle
                        gi = qp%index_from_coordinate(qv2p)
                        q2p = qp%gridind2ind(gi(1), gi(2), gi(3)) ! this is R*q'
                    end select
                    do b2 = 1, dr%n_mode
                        jl = (q2p - 1)*dr%n_mode + b2
                        sr%Xi(il, jl) = buf(b2)
                    end do
                end do
            end do allq2
        end do
        call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call tmr%tock('scattering matrix symmetrization')
    end block symmetrize_scatmat
end subroutine

#include "scattering_isotope.f90"
#include "scattering_threephonon.f90"
#include "scattering_fourphonon.f90"

subroutine sr_destroy(sr)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr

    integer :: il

    if (allocated(sr%Xi)) deallocate (sr%Xi)
    if (allocated(sr%my_qpoints)) deallocate (sr%my_qpoints)
    if (allocated(sr%my_modes)) deallocate (sr%my_modes)
    if (allocated(sr%be)) deallocate (sr%be)
    if (allocated(sr%sigsq)) deallocate (sr%sigsq)
    sr%my_nqpoints = -lo_hugeint
end subroutine

! Function to measure the size of the memory
function sr_size_in_mem(sr) result(mem)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr
    !> size in memory, bytes
    integer(i8) :: mem

    mem = storage_size(sr)
    if (allocated(sr%my_qpoints)) mem = mem + storage_size(sr%my_qpoints)*size(sr%my_qpoints)
    if (allocated(sr%my_modes)) mem = mem + storage_size(sr%my_modes)*size(sr%my_modes)
    if (allocated(sr%be)) mem = mem + storage_size(sr%be)*size(sr%be)
    if (allocated(sr%Xi)) mem = mem + storage_size(sr%Xi)*size(sr%Xi)
    if (allocated(sr%sigsq)) mem = mem + storage_size(sr%sigsq)*size(sr%sigsq)
    mem = mem/8
end function
end module
