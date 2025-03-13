#include "precompilerdefinitions"
module kappa
use konstanter, only: r8, lo_sqtol, lo_kb_hartree, lo_freqtol, lo_kappa_au_to_SI, &
                      lo_groupvel_Hartreebohr_to_ms, lo_twopi
use gottochblandat, only: lo_sqnorm, lo_planck, lo_outerproduct, lo_chop, lo_harmonic_oscillator_cv
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_symmetryoperation, only: lo_operate_on_vector, lo_eigenvector_transformation_matrix, &
                                  lo_operate_on_secondorder_tensor
use type_blas_lapack_wrappers, only: lo_gemm, lo_gemv
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder

use scattering, only: lo_scattering_rates

implicit none

private
public :: get_kappa
public :: get_kappa_offdiag
public :: iterative_solution
public :: symmetrize_kappa
public :: get_viscosity
contains

!> Calculate the thermal conductivity
subroutine get_kappa(dr, qp, uc, temperature, classical, kappa)
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> temperature
    real(r8), intent(in) :: temperature
    !> Are we in the classical limit ?
    logical, intent(in) :: classical
    !> thermal conductivity tensor
    real(r8), dimension(3, 3), intent(out) :: kappa

    real(r8), dimension(3) :: v0, v1
    real(r8) :: om1, cv
    integer :: j

    integer :: q1, b1
    real(r8), dimension(3, 3) :: v2, buf

    kappa = 0.0_r8
    do q1 = 1, qp%n_irr_point
        dr%iq(q1)%kappa = 0.0_r8
        do b1 = 1, dr%n_mode
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle

            if (classical) then
                cv = lo_kb_hartree
            else
                cv = lo_harmonic_oscillator_cv(temperature, om1)
            end if

            ! To ensure the symmetry, we average over the symmetry operation of the crystal
            v2 = 0.0_r8
            do j = 1, uc%sym%n
                v0 = lo_operate_on_vector(uc%sym%op(j), dr%iq(q1)%Fn(:, b1), reciprocal=.true.)
                v1 = lo_operate_on_vector(uc%sym%op(j), dr%iq(q1)%vel(:, b1), reciprocal=.true.)
                v2 = v2 + lo_outerproduct(v0, v1)/uc%sym%n
            end do
            buf = cv*v2/uc%volume
            dr%iq(q1)%kappa(:, :, b1) = buf
            kappa = kappa + buf*qp%ip(q1)%integration_weight
        end do
    end do
    kappa = lo_chop(kappa, sum(abs(kappa))*1e-6_r8)
end subroutine

subroutine get_kappa_offdiag(dr, qp, uc, fc, temperature, classical, mem, mw, kappa_offdiag)
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> temperature
    real(r8), intent(in) :: temperature
    !> Are we in the classical limit ?
    logical, intent(in) :: classical
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> thermal conductivity tensor
    real(r8), dimension(3, 3), intent(out) :: kappa_offdiag

    !> The off diagonal group velocity
    real(r8), dimension(:, :, :), allocatable :: buf_vel
    !> The off diagonal group velocity, squared
    real(r8), dimension(:, :, :, :), allocatable :: buf_velsq
    !> The qpoint
    integer :: iq

    call mem%allocate(buf_vel, [3, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_velsq, [3, 3, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    kappa_offdiag = 0.0_r8

    do iq = 1, qp%n_irr_point
        if (mod(iq, mw%n) .ne. mw%r) cycle
        buf_vel = 0.0_r8
        buf_velsq = 0.0_r8

        ! Calculate the off-diagonal group velocity.
        groupvel: block
            complex(r8), dimension(:, :, :), allocatable :: buf_grad_dynmat
            complex(r8), dimension(:, :), allocatable :: kronegv, buf_egv, buf_egw, buf_cm0, buf_cm1, buf_cm2
            complex(r8), dimension(3) :: cv0
            real(r8), dimension(3) :: v0, v1
            integer :: a1, a2, ia, ib, ic, ix, iy, iz, k, iop, i, ii, j, jj

            ! Some buffers
            call mem%allocate(buf_grad_dynmat, [dr%n_mode, dr%n_mode, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm0, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm1, [dr%n_mode**2, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm2, [dr%n_mode**2, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_egv, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_egw, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(kronegv, [dr%n_mode**2, dr%n_mode**2], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            buf_grad_dynmat = 0.0_r8
            buf_cm0 = 0.0_r8
            buf_cm1 = 0.0_r8
            buf_cm2 = 0.0_r8
            buf_egv = 0.0_r8
            buf_egw = 0.0_r8
            kronegv = 0.0_r8

            ! Dynamical matrix and derivatives
            call fc%dynamicalmatrix(uc, qp%ip(iq), buf_cm0, mem, buf_grad_dynmat, qdirection=[1.0_r8, 0.0_r8, 0.0_r8])

            ! Flatten gradient of dynamical matrix
            do iz = 1, 3
                do a1 = 1, uc%na !ise%n_atom
                do a2 = 1, uc%na !ise%n_atom
                do ix = 1, 3
                do iy = 1, 3
                    ib = (a1 - 1)*3 + ix
                    ia = (a2 - 1)*3 + iy
                    ic = flattenind(a1, a2, ix, iy, dr%n_mode)
                    buf_cm1(ic, iz) = buf_grad_dynmat(ia, ib, iz)/(uc%invsqrtmass(a1)*uc%invsqrtmass(a2))
                end do
                end do
                end do
                end do
            end do

            ! Average over all operations
            kronegv = 0.0_r8
            do k = 1, qp%ip(iq)%n_invariant_operation
                iop = qp%ip(iq)%invariant_operation(k)
                ! Rotate eigenvectors
                call lo_eigenvector_transformation_matrix(buf_cm0, uc%rcart, qp%ip(iq)%r, uc%sym%op(abs(iop)))
                if (iop .lt. 0) then
                    call lo_gemm(buf_cm0, dr%iq(iq)%egv, buf_egw)
                    buf_egw = conjg(buf_egw)
                else
                    call lo_gemm(buf_cm0, dr%iq(iq)%egv, buf_egw)
                end if

                do i = 1, dr%n_mode
                    if (dr%iq(iq)%omega(i) .gt. lo_freqtol) then
                        do a1 = 1, uc%na
                        do ix = 1, 3
                            ib = (a1 - 1)*3 + ix
                            buf_egw(ib, i) = buf_egw(ib, i)*uc%invsqrtmass(a1)/sqrt(dr%iq(iq)%omega(i)*2.0_r8)
                        end do
                        end do
                    else
                        buf_egw(:, i) = 0.0_r8
                    end if
                end do

                do i = 1, dr%n_mode
                do j = 1, dr%n_mode
                    buf_egv = buf_egw*conjg(buf_egw(j, i))
                    do ii = 1, dr%n_mode
                    do jj = 1, dr%n_mode
                        ia = (i - 1)*dr%n_mode + ii
                        ib = (j - 1)*dr%n_mode + jj
                        kronegv(ia, ib) = kronegv(ia, ib) + buf_egv(jj, ii)
                    end do
                    end do
                end do
                end do
            end do
            kronegv = kronegv/real(qp%ip(iq)%n_invariant_operation, r8)
            ! this means sandwich with eigenvectors,  frequencies,
            ! prefactors and masses are already in there.
            call lo_gemm(kronegv, buf_cm1, buf_cm2)

            ! Keep the group velocities?
            do i = 1, dr%n_mode
            do j = 1, dr%n_mode
                ii = (i - 1)*dr%n_mode + j
                jj = (j - 1)*dr%n_mode + i
                cv0 = buf_cm2(ii, :)
                ! remove tiny numbers.
                cv0 = lo_chop(cv0, 1E-10/(lo_groupvel_Hartreebohr_to_ms/1000))
                ! I can take the real part since at the end we sum over
                ! both modes and the imaginary components disappear.
                buf_vel(:, i, j) = real(cv0, r8)
            end do
            end do

            ! Be very careful with degeneracies. Will give very wrong contribution
            ! to thermal transport if not taken care of properly, I'm pretty sure.
            ! Feels like there could be significant double-counting otherwise.
            ! This seems ok for now, but have to keep an eye out.

            ! Rotate out the group velocities to get the symmetry averaged ones
            buf_velsq = 0.0_r8
            do i = 1, dr%n_mode
            do j = 1, dr%n_mode
                v0 = buf_vel(:, i, j)
                do k = 1, qp%ip(iq)%n_full_point
                    iop = qp%ip(iq)%operation_full_point(k)
                    v1 = matmul(uc%sym%op(abs(iop))%m, v0)
                    buf_velsq(:, :, i, j) = buf_velsq(:, :, i, j) + lo_outerproduct(v1, v1)
                end do
            end do
            end do
            buf_velsq = buf_velsq/real(qp%ip(iq)%n_full_point, r8)

            ! cleanup
            call mem%deallocate(buf_grad_dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_egw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(kronegv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block groupvel

        ! Now we can compute the off diagonal contribution
        offdiag: block
            !> prefactor and buffers
            real(r8) :: pref, f0, tau, om1, om2, tau1, tau2
            !> Some integers for do loop on so on
            integer :: jmode, kmode

            pref = qp%ip(iq)%integration_weight/uc%volume

            do jmode = 1, dr%n_mode
                ! Skip gamma for acoustic branches
                if (dr%iq(iq)%omega(jmode) .lt. lo_freqtol) cycle
                om1 = dr%iq(iq)%omega(jmode)
                tau1 = dr%iq(iq)%linewidth(jmode)
                do kmode = 1, dr%n_mode
                    if (jmode .eq. kmode) cycle  ! We only want the off diagonal contribution
                    ! Skip gamma for acoustic branches
                    om2 = dr%iq(iq)%omega(kmode)
                    if (om2 .lt. lo_freqtol) cycle

                    tau2 = dr%iq(iq)%linewidth(kmode)

                    ! This is consistent with the paper, but a bit different from QHGK
                    ! This comes from the fact that we don't work with creation/annihilation but
                    ! directly with displacement/momentup operator
                    if (classical) then
                        f0 = lo_kb_hartree
                    else
                        f0 = 0.5_r8*(lo_harmonic_oscillator_cv(temperature, om1) + &
                                     lo_harmonic_oscillator_cv(temperature, om2))
                    end if

                    tau = (tau1 + tau2)/((tau1 + tau2)**2 + (om1 - om2)**2)
                    kappa_offdiag(:, :) = kappa_offdiag(:, :) + buf_velsq(:, :, jmode, kmode)*tau*f0*pref
                end do ! k mode
            end do ! j mode
        end block offdiag
    end do ! i qpt
    call mem%deallocate(buf_vel, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_velsq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mw%allreduce('sum', kappa_offdiag)

    kappa_offdiag = lo_chop(kappa_offdiag, sum(abs(kappa_offdiag))*1E-6_r8)

contains
    ! Consistent index flattening? Impossibru to get consistent.
    function flattenind(a1, a2, ix, iy, nb) result(i)
        integer, intent(in) :: a1, a2, ix, iy, nb
        integer :: i

        integer :: ia, ib

        ia = (a1 - 1)*3 + ix
        ib = (a2 - 1)*3 + iy
        i = (ib - 1)*nb + ia
    end function
end subroutine

subroutine symmetrize_kappa(kappa, uc)
    !> The kappa to symmetrize
    real(r8), dimension(3, 3), intent(inout) :: kappa
    !> The unit cell
    type(lo_crystalstructure), intent(in) :: uc

    real(r8), dimension(3, 3) :: tmp
    integer :: iop

    tmp = 0.0_r8
    do iop = 1, uc%sym%n
        tmp = tmp + lo_operate_on_secondorder_tensor(uc%sym%op(iop), kappa)
    end do

    kappa = tmp/uc%sym%n
    kappa = lo_chop(kappa, sum(abs(kappa))*1e-6_r8)
end subroutine

subroutine iterative_solution(sr, dr, qp, uc, temperature, niter, tol, classical, mw, mem)
    !> integration weights
    type(lo_scattering_rates), intent(inout) :: sr
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Temperature
    real(r8), intent(in) :: temperature
    !> Max number of iterations
    integer, intent(in) :: niter
    !> Tolerance
    real(r8), intent(in) :: tol
    !> Are we in the classical limit
    logical, intent(in) :: classical
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> The F vector on cart, mode, irr-qpoint
    real(r8), dimension(:, :, :), allocatable :: Fnb
    !> The F vector on full q-point, flattened along mode/qpoint for BLAS vector-matrix multiplication
    real(r8), dimension(:, :), allocatable :: Fbb, buf
    !> The thermal conductivity
    real(r8), dimension(3, 3) :: kappa
    !> Buffer to check convergence
    real(r8), dimension(niter) :: scfcheck
    !> The mixing parameter between iterations
    real(r8) :: mixingparameter
    !> Integer for the do loop
    integer :: iter

    ! set some things and make space
    init: block
        real(r8), dimension(3, 3) :: m0
        mixingparameter = 0.95_r8
        call mem%allocate(Fnb, [3, dr%n_mode, qp%n_irr_point], persistent=.false., &
                          scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(Fbb, [3, dr%n_mode*qp%n_full_point], persistent=.false., &
                          scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf, [3, sr%my_nqpoints], persistent=.false., &
                          scalable=.false., file=__FILE__, line=__LINE__)
        Fnb = 0.0_r8
        Fbb = 0.0_r8
        scfcheck = 0.0_r8
        ! Get the first kappa-value
        call get_kappa(dr, qp, uc, temperature, classical, kappa)
        m0 = kappa*lo_kappa_au_to_SI
        if (mw%talk) write (*, "(1X,I4,6(1X,F14.4))") 0, m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
    end block init

    scfloop: do iter = 1, niter
        ! get the Fn values all across the BZ
        foldout: block
            real(r8), dimension(3) :: v
            integer :: iq, jq, iop, b1, i2

            Fnb = 0.0_r8
            Fbb = 0.0_r8
            do iq = 1, qp%n_full_point
                if (mod(iq, mw%n) .ne. mw%r) cycle
                iop = qp%ap(iq)%operation_from_irreducible
                jq = qp%ap(iq)%irreducible_index
                do b1 = 1, dr%n_mode
                    i2 = (iq - 1)*dr%n_mode + b1
                    if (iop .gt. 0) then
                        v = lo_operate_on_vector(uc%sym%op(iop), dr%iq(jq)%Fn(:, b1), reciprocal=.true.)
                        Fbb(:, i2) = v
                    else
                        v = -lo_operate_on_vector(uc%sym%op(abs(iop)), dr%iq(jq)%Fn(:, b1), reciprocal=.true.)
                        Fbb(:, i2) = v
                    end if
                end do
            end do
            call mw%allreduce('sum', Fbb)
        end block foldout

        ! Multiply the off-diagonal part of Xi with F
        applyXi: block
            integer :: il, b1, q1, a

            ! We use BLAS for this, and we have to do it for each cartesian direction
            do a = 1, 3
                call lo_gemv(sr%Xi, Fbb(a, :), buf(a, :))
            end do
            ! And now we distribute the results on the irreducible qpoints
            do il = 1, sr%my_nqpoints
                q1 = sr%my_qpoints(il)
                b1 = sr%my_modes(il)
                Fnb(:, b1, q1) = -buf(:, il)/dr%iq(q1)%qs(b1)
            end do
            call mw%allreduce('sum', Fnb)
        end block applyXi

        ! make sure degeneracies are satisfied properly and add the previous thing
        distributeF: block
            real(r8), dimension(3) :: v0
            integer :: q1, b1, b2, j
            do q1 = 1, qp%n_irr_point
                do b1 = 1, dr%n_mode
                    v0 = 0.0_r8
                    do j = 1, dr%iq(q1)%degeneracy(b1)
                        b2 = dr%iq(q1)%degenmode(j, b1)
                        v0 = v0 + Fnb(:, b2, q1)
                    end do
                    v0 = v0/real(dr%iq(q1)%degeneracy(b1), r8)
                    do j = 1, dr%iq(q1)%degeneracy(b1)
                        b2 = dr%iq(q1)%degenmode(j, b1)
                        Fnb(:, b2, q1) = v0
                    end do
                    ! Add the previous thing
                    Fnb(:, b1, q1) = dr%iq(q1)%F0(:, b1) + Fnb(:, b1, q1)
                end do
            end do
        end block distributeF

        ! Add everything together and check convergency
        addandcheck: block
            real(r8), dimension(3, 3) :: m0
            real(r8) :: g0, g1, g2
            integer :: i, j

            g0 = 0.0_r8
            do i = 1, dr%n_irr_qpoint
            do j = 1, dr%n_mode
                g1 = lo_sqnorm(dr%iq(i)%Fn(:, j) - Fnb(:, j, i))
                g2 = lo_sqnorm(dr%iq(i)%Fn(:, j))
                if (g2 .gt. lo_sqtol) then
                    g0 = g0 + g1/g2
                end if
                dr%iq(i)%Fn(:, j) = dr%iq(i)%Fn(:, j)*(1.0_r8 - mixingparameter) + &
                                    mixingparameter*(Fnb(:, j, i))
            end do
            end do
            scfcheck(iter) = g0/qp%n_irr_point/dr%n_mode

            ! Get the current kappa, to print to stdout.
            call get_kappa(dr, qp, uc, temperature, classical, kappa)
            m0 = kappa*lo_kappa_au_to_SI
            if (mw%r .eq. 0) write (*, "(1X,I4,6(1X,F14.4),2X,ES10.3)") &
                iter, m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3), scfcheck(iter)

            ! Check for convergence. The criterion is that the relative difference between the new
            ! and old Fn is to be one part in 1E-5, for two consecutive iterations.
            if (iter .ge. 3) then
                g0 = sum(scfcheck(iter - 2:iter))
                if (g0 .lt. tol) then
                    exit scfloop
                end if
            end if

            ! We are not converged if we made it here.
            ! If we had too many iterations I want to adjust the mixing a little
            if (iter .gt. 15 .and. mixingparameter .gt. 0.50) then
                mixingparameter = mixingparameter*0.98_r8
            end if
        end block addandcheck
    end do scfloop
    call mem%deallocate(Fnb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(Fbb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Calculate the thermal conductivity
subroutine get_viscosity(dr, qp, uc, temperature, classical, mu)
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> temperature
    real(r8), intent(in) :: temperature
    !> Are we in the classical limit ?
    logical, intent(in) :: classical
    !> thermal conductivity tensor
    real(r8), dimension(3, 3, 3, 3), intent(out) :: mu

    real(r8), dimension(3) :: v0, v1, kr
    real(r8) :: om1, cv, pref, n
    integer :: j

    integer :: q1, b1, iq, a, b, c, d, iop
    real(r8), dimension(3, 3) :: v2, buf

    mu = 0.0_r8
    do q1=1, qp%n_full_point
        iq = qp%ap(q1)%irreducible_index
        kr = qp%ap(q1)%r * lo_twopi
        iop = qp%ap(q1)%operation_from_irreducible
        do b1=1, dr%n_mode
            om1 = dr%aq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle
            n = lo_planck(temperature, om1)
            pref = n * (n + 1.0_r8) / uc%volume / lo_kb_hartree / temperature
            v2 = 0.0_r8
            if (iop .gt. 0) then
                v0 = lo_operate_on_vector(uc%sym%op(iop), dr%iq(iq)%Fn(:, b1))
                v1 = lo_operate_on_vector(uc%sym%op(iop), dr%iq(iq)%vel(:, b1))
            else
                v0 = -lo_operate_on_vector(uc%sym%op(abs(iop)), dr%iq(iq)%Fn(:, b1))
                v1 = -lo_operate_on_vector(uc%sym%op(abs(iop)), dr%iq(iq)%vel(:, b1))
            end if
            v2 = lo_outerproduct(v0, v1)
            do a=1, 3
            do b=1, 3
            do c=1, 3
            do d=1, 3
                mu(a, b, c, d) = mu(a, b, c, d) + pref * kr(a) * v0(b) * kr(c) * v1(d) / qp%n_full_point
!               mu(a, b, c, d) = mu(a, b, c, d) + pref * v2(a, b) * kr(c) * kr(d) / qp%n_full_point
            end do
            end do
            end do
            end do
        end do
    end do
    mu = lo_chop(mu, sum(abs(mu))*1e-6_r8)
end subroutine
end module
