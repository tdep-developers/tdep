#include "precompilerdefinitions"
module pbe
use konstanter, only: r8, lo_tol, lo_sqtol, lo_pi, lo_kb_hartree, lo_freqtol, lo_huge, lo_kappa_au_to_SI, &
                      lo_phonongroupveltol
use gottochblandat, only: lo_sqnorm, lo_planck, lo_outerproduct, lo_chop
use mpi_wrappers, only: lo_mpi_helper, MPI_SUM, MPI_DOUBLE_PRECISION, MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_symmetryoperation, only: lo_operate_on_vector
use type_blas_lapack_wrappers, only: lo_dgelss

! local
use phononevents
implicit none
private
public :: get_kappa, calculate_qs, get_selfconsistent_solution
contains

!> Calculate the thermal conductivity
subroutine get_kappa(dr, qp, uc, temperature, kappa)
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> temperature
    real(r8), intent(in) :: temperature
    !> thermal conductivity tensor
    real(r8), dimension(3, 3), intent(out) :: kappa

    real(r8), dimension(3) :: v0, v1
    real(r8) :: n, f0, omega, omthres, prefactor
    integer :: i, j, k, l
    !integer :: iop,ifull

    omthres = dr%omega_min*0.5_r8
    prefactor = 1.0_r8/(uc%volume*lo_kb_hartree*temperature)
    do i = 1, qp%n_full_point
        dr%aq(i)%kappa = 0.0_r8
        k = qp%ap(i)%operation_from_irreducible
        do j = 1, dr%n_mode
            ! Skip gamma for acoustic branches
            if (dr%aq(i)%omega(j) .lt. omthres) cycle
            ! Which operation takes this point from the wedge to here
            l = qp%ap(i)%irreducible_index
            ! Rotate things to this points. Negative is the time reversal thingy, but does not really matter here.
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
            dr%aq(i)%kappa(:, :, j) = prefactor*f0*lo_outerproduct(v1, v0)
        end do
    end do

    ! Sum it ip!
    kappa = 0.0_r8
    do i = 1, qp%n_full_point
        do j = 1, dr%n_mode
            kappa = kappa + dr%aq(i)%kappa(:, :, j)/qp%n_full_point
        end do
    end do
    f0 = sum(abs(kappa))
    kappa = lo_chop(kappa, f0*1E-6_r8)
end subroutine

!> Calculate the QS term
subroutine calculate_qs(qp, sc, dr, temperature, mw, mem)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !>  scattering rates and strengths
    type(lo_threephononevents), intent(inout) :: sc
    !> temperature
    real(r8), intent(in) :: temperature
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), parameter :: threephonon_prefactor = lo_pi/4.0_r8
    real(r8), parameter :: isotope_prefactor = lo_pi/2.0_r8
    real(r8), dimension(:, :), allocatable :: buf_plus, buf_minus, buf_iso
    real(r8) :: om1, om2, om3, n1, n2, n3, f0, f1, f2, omthres, qs_boundary, velnorm
    integer :: gi1, gi2, gi3, b1, b2, b3, iq, lqp, j

    ! Threshold for omega to be zero
    omthres = dr%omega_min*0.2_r8
    ! Set things to zero
    do iq = 1, qp%n_irr_point
        dr%iq(iq)%linewidth = 0.0_r8
        dr%iq(iq)%p_plus = 0.0_r8
        dr%iq(iq)%p_minus = 0.0_r8
        dr%iq(iq)%p_iso = 0.0_r8
        dr%iq(iq)%qs = 0.0_r8
        dr%iq(iq)%F0 = 0.0_r8
        dr%iq(iq)%Fn = 0.0_r8
        dr%iq(iq)%mfp = 0.0_r8
        dr%iq(iq)%scalar_mfp = 0.0_r8
    end do

    ! Some temporary buffers
    call mem%allocate(buf_plus, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., &
                      file=__FILE__, line=__LINE__)
    call mem%allocate(buf_minus, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., &
                      file=__FILE__, line=__LINE__)
    call mem%allocate(buf_iso, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., &
                      file=__FILE__, line=__LINE__)
    buf_plus = 0.0_r8
    buf_minus = 0.0_r8
    buf_iso = 0.0_r8

    ! fetch all the scattering rates for this rank
    lqp = 0
    qploop1: do iq = 1, qp%n_irr_point
        ! Make it parallel
        if (mod(iq, mw%n) .ne. mw%r) cycle
        ! Increment counter
        lqp = lqp + 1
        b1loop: do b1 = 1, dr%n_mode
            ! Three-phonon part
            do b2 = 1, dr%n_mode
            do b3 = 1, dr%n_mode
                ! First the plus-events
                do j = 1, sc%q(lqp)%plus(b1, b2, b3)%n
                    ! grid-indices
                    gi1 = sc%q(lqp)%gi1
                    gi2 = sc%q(lqp)%plus(b1, b2, b3)%e(j)%gi2
                    gi3 = sc%q(lqp)%plus(b1, b2, b3)%e(j)%gi3
                    ! frequencies
                    om1 = dr%aq(gi1)%omega(b1)
                    om2 = dr%aq(gi2)%omega(b2)
                    om3 = dr%aq(gi3)%omega(b3)
                    ! init scattering rate to zero
                    sc%q(lqp)%plus(b1, b2, b3)%e(j)%W = 0.0_r8
                    ! skip things with gamma
                    if (minval([om1, om2, om3]) .lt. omthres) cycle
                    ! distribution functions
                    n1 = lo_planck(temperature, om1)
                    n2 = lo_planck(temperature, om2)
                    n3 = lo_planck(temperature, om3)

                    f0 = n1*n2*(n3 + 1)*threephonon_prefactor
                    f0 = f0*sc%q(lqp)%plus(b1, b2, b3)%e(j)%psisquare*sc%q(lqp)%plus(b1, b2, b3)%e(j)%deltafunction
                    ! Store intermediate scattering strength
                    sc%q(lqp)%plus(b1, b2, b3)%e(j)%W = f0
                    ! Store actual scattering strength
                    buf_plus(b1, iq) = buf_plus(b1, iq) + f0
                end do
                ! Then the minus-events
                do j = 1, sc%q(lqp)%minus(b1, b2, b3)%n
                    ! grid-indices
                    gi1 = sc%q(lqp)%gi1
                    gi2 = sc%q(lqp)%minus(b1, b2, b3)%e(j)%gi2
                    gi3 = sc%q(lqp)%minus(b1, b2, b3)%e(j)%gi3
                    ! frequencies
                    om1 = dr%aq(gi1)%omega(b1)
                    om2 = dr%aq(gi2)%omega(b2)
                    om3 = dr%aq(gi3)%omega(b3)
                    sc%q(lqp)%minus(b1, b2, b3)%e(j)%W = 0.0_r8
                    if (minval([om1, om2, om3]) .lt. omthres) cycle
                    ! distribution functions
                    n1 = lo_planck(temperature, om1)
                    n2 = lo_planck(temperature, om2)
                    n3 = lo_planck(temperature, om3)
                    !
                    f0 = n1*(n2 + 1)*(n3 + 1)*threephonon_prefactor
                    f0 = f0*sc%q(lqp)%minus(b1, b2, b3)%e(j)%psisquare* &
                         sc%q(lqp)%minus(b1, b2, b3)%e(j)%deltafunction
                    sc%q(lqp)%minus(b1, b2, b3)%e(j)%W = f0
                    ! accumulate actual scattering strength
                    buf_minus(b1, iq) = buf_minus(b1, iq) + f0
                end do
            end do
            end do

            ! Isotope part
            do b2 = 1, dr%n_mode
                do j = 1, sc%iq(lqp)%band(b1, b2)%n
                    ! grid-indices
                    gi1 = sc%iq(lqp)%gi1
                    gi2 = sc%iq(lqp)%band(b1, b2)%e(j)%gi2
                    ! frequencies
                    om1 = dr%aq(gi1)%omega(b1)
                    om2 = dr%aq(gi2)%omega(b2)
                    ! distribution functions
                    n1 = lo_planck(temperature, om1)
                    n2 = lo_planck(temperature, om2)
                    ! sanity check
                    sc%iq(lqp)%band(b1, b2)%e(j)%W = 0.0_r8
                    if (minval([om1, om2]) .lt. omthres) cycle
                    f0 = isotope_prefactor*n1*(n2 + 1.0_r8)*om1*om2* &
                         sc%iq(lqp)%band(b1, b2)%e(j)%scatterstrength* &
                         sc%iq(lqp)%band(b1, b2)%e(j)%deltafunction
                    sc%iq(lqp)%band(b1, b2)%e(j)%W = f0
                    ! Just add it straight
                    buf_iso(b1, iq) = buf_iso(b1, iq) + f0
                end do
            end do
        end do b1loop

        ! Fix degeneracies? Always a sensible idea
        if (sc%correctionlevel .ge. 2) then
            ! make sure degeneracies are satisfied properly.
            do b1 = 1, dr%n_mode
                f0 = 0.0_r8
                f1 = 0.0_r8
                f2 = 0.0_r8
                do j = 1, dr%iq(iq)%degeneracy(b1)
                    b2 = dr%iq(iq)%degenmode(j, b1)
                    f0 = f0 + buf_plus(b2, iq)
                    f1 = f1 + buf_minus(b2, iq)
                    f2 = f2 + buf_iso(b2, iq)
                end do
                f0 = f0/real(dr%iq(iq)%degeneracy(b1), r8)
                f1 = f1/real(dr%iq(iq)%degeneracy(b1), r8)
                f2 = f2/real(dr%iq(iq)%degeneracy(b1), r8)
                do j = 1, dr%iq(iq)%degeneracy(b1)
                    b2 = dr%iq(iq)%degenmode(j, b1)
                    buf_plus(b2, iq) = f0
                    buf_minus(b2, iq) = f1
                    buf_iso(b2, iq) = f2
                end do
            end do
        end if
    end do qploop1

    ! reduce over MPI
    call mw%allreduce('sum', buf_plus)
    call mw%allreduce('sum', buf_minus)
    call mw%allreduce('sum', buf_iso)

    ! ! Sum stuff up
    do iq = 1, qp%n_irr_point
    do b1 = 1, dr%n_mode
        ! Store the p+ p- piso thingies
        dr%iq(iq)%p_plus(b1) = buf_plus(b1, iq)
        dr%iq(iq)%p_minus(b1) = buf_minus(b1, iq)
        dr%iq(iq)%p_iso(b1) = buf_iso(b1, iq)

        n1 = lo_planck(temperature, dr%iq(iq)%omega(b1))
        velnorm = norm2(dr%iq(iq)%vel(:, b1))
        ! Boundary scattering.
        if (sc%mfpmax .gt. 0.0_r8) then
            ! Per mode. Say that the phonon will be scattered after travelling lmax distance.
            if (velnorm .gt. lo_phonongroupveltol) then
                f1 = sc%mfpmax/velnorm
                qs_boundary = (n1*(n1 + 1.0_r8))/f1
            else
                ! very large distance that should never matter.
                qs_boundary = 0.0_r8
            end if
        else
            qs_boundary = 0.0_r8
        end if

        ! Get qs, linewidth and mean free path
        if (dr%iq(iq)%omega(b1) .gt. omthres) then
            dr%iq(iq)%qs(b1) = dr%iq(iq)%p_plus(b1) + &
                               dr%iq(iq)%p_minus(b1)*0.5_r8 + &
                               dr%iq(iq)%p_iso(b1) + &
                               qs_boundary
            dr%iq(iq)%linewidth(b1) = 0.5_r8*dr%iq(iq)%qs(b1)/(n1*(n1 + 1.0_r8))

            if (velnorm .gt. lo_phonongroupveltol) then
                dr%iq(iq)%mfp(:, b1) = dr%iq(iq)%vel(:, b1)*0.5_r8/dr%iq(iq)%linewidth(b1)
                dr%iq(iq)%scalar_mfp(b1) = velnorm*0.5_r8/dr%iq(iq)%linewidth(b1)
                dr%iq(iq)%F0(:, b1) = dr%iq(iq)%mfp(:, b1)*dr%iq(iq)%omega(b1)/temperature
                dr%iq(iq)%Fn(:, b1) = dr%iq(iq)%F0(:, b1)
            end if
        end if
    end do
    end do

    call mem%deallocate(buf_plus, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_minus, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_iso, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Self-consistently solve everything
subroutine get_selfconsistent_solution(sc, dr, qp, uc, temperature, niter, tol, mw, mem)
    !> integration weights
    type(lo_threephononevents), intent(inout) :: sc
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
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:, :, :), allocatable :: Fnb, Fbb
    real(r8), dimension(3, 3) :: kappa
    real(r8), dimension(niter) :: scfcheck
    real(r8) :: mixingparameter
    integer :: iter

    ! set some things and make space
    init: block
        real(r8), dimension(3, 3) :: m0
        mixingparameter = 0.95_r8
        allocate (Fnb(3, dr%n_mode, dr%n_irr_qpoint))
        allocate (Fbb(3, dr%n_mode, dr%n_full_qpoint))
        Fnb = 0.0_r8
        Fbb = 0.0_r8
        scfcheck = 0.0_r8
        ! Get the first kappa-value
        call get_kappa(dr, qp, uc, temperature, kappa)
        m0 = kappa*lo_kappa_au_to_SI
        if (mw%talk) write (*, "(1X,I4,6(1X,F14.4))") 0, m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
    end block init

    scfloop: do iter = 1, niter
        ! get the Fn values all across the BZ
        foldout: block
            real(r8), dimension(3) :: v
            integer :: iq, jq, iop, b1

            Fbb = 0.0_r8
            do iq = 1, qp%n_full_point
                if (mod(iq, mw%n) .ne. mw%r) cycle
                iop = qp%ap(iq)%operation_from_irreducible
                jq = qp%ap(iq)%irreducible_index
                do b1 = 1, dr%n_mode
                    if (iop .gt. 0) then
                        v = lo_operate_on_vector(uc%sym%op(iop), dr%iq(jq)%Fn(:, b1), reciprocal=.true.)
                        Fbb(:, b1, iq) = v
                    else
                        v = -lo_operate_on_vector(uc%sym%op(abs(iop)), dr%iq(jq)%Fn(:, b1), reciprocal=.true.)
                        Fbb(:, b1, iq) = v
                    end if
                end do
            end do
            call mw%allreduce('sum', Fbb)
        end block foldout

        ! update F to new values
        updateF: block
            real(r8), dimension(3) :: Fp, Fpp, v0
            real(r8) :: iQs, W
            integer :: lqp, i, j, b1, b2, b3, ii, jj
            integer :: iq

            ! Some temporary space

            Fnb = 0.0_r8
            lqp = 0
            do iq = 1, qp%n_irr_point
                if (mod(iq, mw%n) .ne. mw%r) cycle
                lqp = lqp + 1
                b1loop: do b1 = 1, dr%n_mode
                    ! prefetch some stuff
                    if (dr%iq(iq)%linewidth(b1) .gt. lo_freqtol) then
                        iQS = 1.0_r8/dr%iq(iq)%qs(b1)
                    else
                        cycle b1loop
                    end if
                    ! three-phonon stuff
                    v0 = 0.0_r8
                    do b2 = 1, dr%n_mode
                    do b3 = 1, dr%n_mode
                        do j = 1, sc%q(lqp)%plus(b1, b2, b3)%n
                            ! Scatterstrength
                            W = sc%q(lqp)%plus(b1, b2, b3)%e(j)%W
                            ii = sc%q(lqp)%plus(b1, b2, b3)%e(j)%gi2
                            jj = sc%q(lqp)%plus(b1, b2, b3)%e(j)%gi3
                            Fp = Fbb(:, b2, ii)
                            Fpp = Fbb(:, b3, jj)
                            ! F-difference
                            v0 = v0 + (Fp + Fpp)*W*iQs
                        end do
                        do j = 1, sc%q(lqp)%minus(b1, b2, b3)%n
                            ! Scatterstrength
                            W = sc%q(lqp)%minus(b1, b2, b3)%e(j)%W
                            ii = sc%q(lqp)%minus(b1, b2, b3)%e(j)%gi2
                            jj = sc%q(lqp)%minus(b1, b2, b3)%e(j)%gi3
                            Fp = Fbb(:, b2, ii)
                            Fpp = Fbb(:, b3, jj)
                            v0 = v0 + (Fp + Fpp)*W*iQs*0.5_r8
                        end do
                    end do
                    end do
                    ! isotope stuff
                    do b2 = 1, dr%n_mode
                        do j = 1, sc%iq(lqp)%band(b1, b2)%n
                            W = sc%iq(lqp)%band(b1, b2)%e(j)%W
                            ii = sc%iq(lqp)%band(b1, b2)%e(j)%gi2
                            Fp = Fbb(:, b2, ii)
                            v0 = v0 + W*Fp*iQs
                        end do
                    end do
                    ! Update Fn
                    Fnb(:, b1, iq) = Fnb(:, b1, iq) - v0
                end do b1loop

                if (sc%correctionlevel .ge. 2) then
                    ! make sure degeneracies are satisfied properly.
                    do b1 = 1, dr%n_mode
                        v0 = 0.0_r8
                        do i = 1, dr%iq(iq)%degeneracy(b1)
                            b2 = dr%iq(iq)%degenmode(i, b1)
                            v0 = v0 + Fnb(:, b2, iq)
                        end do
                        v0 = v0/real(dr%iq(iq)%degeneracy(b1), r8)
                        do i = 1, dr%iq(iq)%degeneracy(b1)
                            b2 = dr%iq(iq)%degenmode(i, b1)
                            Fnb(:, b2, iq) = v0
                        end do
                    end do
                end if

                if (sc%correctionlevel .ge. 3) then
                end if
            end do
            call mw%allreduce('sum', Fnb)

            ! Add the previous thing
            do iq = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                Fnb(:, b1, iq) = dr%iq(iq)%F0(:, b1) + Fnb(:, b1, iq)
            end do
            end do

        end block updateF

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
                dr%iq(i)%Fn(:, j) = dr%iq(i)%Fn(:, j)*(1.0_r8 - mixingparameter) + mixingparameter*(Fnb(:, j, i))
            end do
            end do
            scfcheck(iter) = g0/qp%n_irr_point/dr%n_mode

            ! Check for convergence. The criterion is that the relative difference between the new
            ! and old Fn is to be one part in 1E-5, for two consecutive iterations.
            if (iter .ge. 3) then
                g0 = sum(scfcheck(iter - 2:iter))
                if (g0 .lt. tol) then
                    exit scfloop
                end if
            end if

            ! We are not converged if we made it here. Get the current kappa, to print to stdout.
            call get_kappa(dr, qp, uc, temperature, kappa)
            m0 = kappa*lo_kappa_au_to_SI
            if (mw%r .eq. 0) write (*, "(1X,I4,6(1X,F14.4),2X,E10.3)") &
                iter, m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3), scfcheck(iter)

            ! If we had too many iterations I want to adjust the mixing a little
            if (iter .gt. 15) then
                mixingparameter = mixingparameter*0.98_r8
            end if
        end block addandcheck
    end do scfloop

end subroutine

! !> Self-consistently solve everything
! subroutine get_selfconsistent_solution(sc,dr,qp,uc,mw,temperature,niter,tol)
!     !> integration weights
!     type(lo_threephononevents), intent(inout) :: sc
!     !> dispersions
!     type(lo_phonon_dispersions), intent(inout) :: dr
!     !> q-mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> structure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> Temperature
!     real(r8), intent(in) :: temperature
!     !> Max number of iterations
!     integer, intent(in) :: niter
!     !> Tolerance
!     real(r8), intent(in) :: tol
!
!     real(r8), dimension(:,:,:), allocatable :: Fnb,Fbb
!     real(r8), dimension(3,3) :: kappa
!     real(r8), dimension(niter) :: scfcheck
!     real(r8) :: mixingparameter
!     integer :: iter
!
!     ! set some things and make space
!     init: block
!         real(r8), dimension(3,3) :: m0
!         mixingparameter=0.95_r8
!         lo_allocate(Fnb(3,dr%n_mode,dr%n_irr_qpoint))
!         lo_allocate(Fbb(3,dr%n_mode,dr%n_full_qpoint))
!         Fnb=0.0_r8
!         Fbb=0.0_r8
!         scfcheck=0.0_r8
!         ! Get the first kappa-value
!         call get_kappa(dr,qp,uc,temperature,kappa)
!         m0=kappa*lo_kappa_au_to_SI
!         if ( mw%r .eq. 0 ) write(*,"(1X,I4,6(1X,F14.4))") 0,m0(1,1),m0(2,2),m0(3,3),m0(1,2),m0(1,3),m0(2,3)
!     end block init
!
! scfloop: do iter=1,niter
!     ! get the Fn values all across the BZ
!     foldout: block
!         real(r8), dimension(3) :: v
!         integer :: lqp,i,b1,ii,jj
!
!         Fbb=0.0_r8
!         lqp=0
!         do i=1,qp%n_irr_point
!             ! Make it parallel
!             if ( mod(i,mw%n) .ne. mw%r ) cycle
!             ! Increment counter
!             lqp=lqp+1
!             ii=qp%ap(i)%operation_from_irreducible
!             jj=qp%ap(i)%irreducible_index
!             ! Take care about time-reversal symmetries
!             do b1=1,dr%n_mode
!                 if ( ii .gt. 0 ) then
!                     v=lo_operate_on_vector(uc%sym%op(ii),dr%iq(jj)%Fn(:,b1),reciprocal=.true.)
!                     Fbb(:,b1,i)=v
!                 else
!                     v=-lo_operate_on_vector(uc%sym%op(abs(ii)),dr%iq(jj)%Fn(:,b1),reciprocal=.true.)
!                     Fbb(:,b1,i)=v
!                 endif
!             enddo
!         enddo
!         call mw%allreduce('sum',Fbb)
!     end block foldout
!
!     ! update F to new values
!     updateF: block
!         real(r8), dimension(3) :: Fp,Fpp,v0
!         real(r8) :: iQs,W
!         integer :: lqp,i,j,b1,b2,b3,ii,jj
!         Fnb=0.0_r8
!         ! Do the whole SCF thing
!         lqp=0
!         do i=1,qp%n_irr_point
!             ! Make it parallel
!             if ( mod(i,mw%n) .ne. mw%r ) cycle
!             ! Increment counter
!             lqp=lqp+1
!             do b1=1,dr%n_mode
!                 ! prefetch some stuff
!                 if ( dr%iq(i)%linewidth(b1) .gt. lo_freqtol ) then
!                     iQS=1.0_r8/dr%iq(i)%qs(b1)
!                 else
!                     cycle
!                 endif
!                 ! three-phonon stuff
!                 v0=0.0_r8
!                 do b2=1,dr%n_mode
!                 do b3=1,dr%n_mode
!                     do j=1,sc%q(lqp)%plus(b1,b2,b3)%n
!                         ! Scatterstrength
!                         W=sc%q(lqp)%plus(b1,b2,b3)%e(j)%W
!                         ii=sc%q(lqp)%plus(b1,b2,b3)%e(j)%gi2
!                         jj=sc%q(lqp)%plus(b1,b2,b3)%e(j)%gi3
!                         Fp=Fbb(:,b2,ii)
!                         Fpp=Fbb(:,b3,jj)
!                         ! F-difference
!                         v0=v0+(Fp+Fpp)*W*iQs
!                     enddo
!                     do j=1,sc%q(lqp)%minus(b1,b2,b3)%n
!                         ! Scatterstrength
!                         W=sc%q(lqp)%minus(b1,b2,b3)%e(j)%W
!                         ii=sc%q(lqp)%minus(b1,b2,b3)%e(j)%gi2
!                         jj=sc%q(lqp)%minus(b1,b2,b3)%e(j)%gi3
!                         Fp=Fbb(:,b2,ii)
!                         Fpp=Fbb(:,b3,jj)
!                         v0=v0+(Fp+Fpp)*W*iQs*0.5_r8
!                     enddo
!                 enddo
!                 enddo
!                 ! isotope stuff
!                 do b2=1,dr%n_mode
!                     do j=1,sc%iq(lqp)%band(b1,b2)%n
!                         W=sc%iq(lqp)%band(b1,b2)%e(j)%W
!                         ii=sc%iq(lqp)%band(b1,b2)%e(j)%gi2
!                         Fp=Fbb(:,b2,ii)
!                         v0=v0+W*Fp*iQs
!                     enddo
!                 enddo
!                 ! Update Fn
!                 Fnb(:,b1,i)=Fnb(:,b1,i)-v0/qp%n_full_point
!             enddo
!         enddo
!     end block updateF
!
!     if ( sc%correctionlevel .ge. 5 ) then
!         !! make sure the update is nice and symmetric
!         !symandstore: block
!         !    real(r8), dimension(3) :: v0,v1
!         !    integer :: lqp,i,b1
!
!         !    lqp=0
!         !    do i=1,qp%n_irr_point
!         !        ! Make it parallel
!         !        if ( mod(i,mw%n) .ne. mw%r ) cycle
!         !        ! Increment counter
!         !        lqp=lqp+1
!         !        ! coefficient matrix thingy
!         !        if ( sc%q(lqp)%nvar .eq. 0 ) then
!         !            ! If symmetry does not help, just plain update.
!         !            Fnb(:,:,i)=Fnb(:,:,i)+dr%iq(i)%F0
!         !            cycle
!         !        else
!         !            ! symmetry probably helps
!         !            do b1=1,dr%n_mode
!         !                ! updated F-thing
!         !                v0=Fnb(:,b1,i)
!         !                ! lsq solver
!         !                v1(1:sc%q(lqp)%nvar)=matmul(transpose(sc%q(lqp)%invariant_operator),v0)
!         !                v1=matmul(sc%q(lqp)%invariant_operator,v1(1:sc%q(lqp)%nvar))
!         !                ! Update with the symmetrized thing
!         !                Fnb(:,b1,i)=v1+dr%iq(i)%F0(:,b1)
!         !            enddo
!         !        endif
!         !    enddo
!         !end block symandstore
!     else
!         ! just store it with no corrections
!         onlystore: block
!             integer :: i,b1
!             do i=1,qp%n_irr_point
!                 if ( mod(i,mw%n) .ne. mw%r ) cycle
!                 do b1=1,dr%n_mode
!                     Fnb(:,b1,i)=dr%iq(i)%F0(:,b1)+Fnb(:,b1,i)
!                 enddo
!             enddo
!         end block onlystore
!     endif
!
!     ! perhaps fix degeneracies?
!     if ( sc%correctionlevel .ge. 3 ) then
!     fixdegen: block
!         real(r8), dimension(3) :: v0
!         real(r8) :: avgfactor
!         integer, dimension(dr%n_mode) :: di
!         integer :: lqp,i,j,b1,b2
!
!         lqp=0
!         do i=1,qp%n_irr_point
!             ! Make it parallel
!             if ( mod(i,mw%n) .ne. mw%r ) cycle
!             ! Increment counter
!             lqp=lqp+1
!             di=1
!             do b1=1,dr%n_mode
!                 ! skip if already taken care of
!                 if ( di(b1) .eq. 0 ) cycle
!                 ! average and distribute
!                 if ( dr%iq(i)%degeneracy(b1) .gt. 1 ) then
!                     v0=0.0_r8
!                     avgfactor=1.0_r8/(dr%iq(i)%degeneracy(b1))
!                     do j=1,dr%iq(i)%degeneracy(b1)
!                         b2=dr%iq(i)%degenmode(j,b1)
!                         v0=v0+Fnb(:,b2,i)*avgfactor
!                     enddo
!                     do j=1,dr%iq(i)%degeneracy(b1)
!                         b2=dr%iq(i)%degenmode(j,b1)
!                         Fnb(:,b2,i)=v0
!                         di(b2)=0
!                     enddo
!                 endif
!             enddo
!         enddo
!     end block fixdegen
!     endif
!
!     ! Add everything together and check convergency
!     addandcheck: block
!         real(r8), dimension(3,3) :: m0
!         real(r8) :: g0,g1,g2
!         integer :: i,j
!         ! sum up all the contributions
!         call mpi_allreduce(MPI_IN_PLACE,Fnb,qp%n_irr_point*dr%n_mode*3,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!         !
!         g0=0.0_r8
!         do i=1,dr%n_irr_qpoint
!         do j=1,dr%n_mode
!             g1=lo_sqnorm(dr%iq(i)%Fn(:,j)-Fnb(:,j,i))
!             g2=lo_sqnorm(dr%iq(i)%Fn(:,j))
!             if ( g2 .gt. lo_sqtol ) then
!                 g0=g0+g1/g2
!             endif
!             dr%iq(i)%Fn(:,j)=dr%iq(i)%Fn(:,j)*(1.0_r8-mixingparameter)+mixingparameter*(Fnb(:,j,i))
!         enddo
!         enddo
!         scfcheck(iter)=g0/qp%n_irr_point/dr%n_mode
!
!         ! Check for convergence. The criterion is that the relative difference between the new
!         ! and old Fn is to be one part in 1E-5, for two consecutive iterations.
!         if ( iter .ge. 3 ) then
!             g0=sum(scfcheck(iter-2:iter))
!             if ( g0 .lt. tol ) then
!                 exit scfloop
!             endif
!         endif
!
!         ! We are not converged if we made it here. Get the current kappa, to print to stdout.
!         call get_kappa(dr,qp,uc,temperature,kappa)
!         m0=kappa*lo_kappa_au_to_SI
!         if ( mw%r .eq. 0 ) write(*,"(1X,I4,6(1X,F14.4),2X,E10.3)") iter,m0(1,1),m0(2,2),m0(3,3),m0(1,2),m0(1,3),m0(2,3),scfcheck(iter)
!
!         ! If we had too many iterations I want to adjust the mixing a little
!         if ( iter .gt. 15 ) then
!             mixingparameter=mixingparameter*0.98_r8
!         endif
!     end block addandcheck
! enddo scfloop
!
! end subroutine

end module
