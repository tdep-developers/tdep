module fourthorder
!! Compute the fourth order contribution to the anharmonic free energy
use konstanter, only: r8, lo_freqtol, lo_twopi, lo_imag, lo_kb_Hartree, lo_exitcode_param
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_planck, &
                          lo_planck_deriv, lo_planck_secondderiv
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_timetracker, only: lo_timer
use lo_fftgrid_helper, only: fft_fourth_grid_index

implicit none
private

public :: free_energy_fourthorder
public :: free_energy_fourthorder_secondorder

contains
!> Calculate the fourth order free energy
subroutine free_energy_fourthorder(uc, fcf, qp, dr, temperature, df4, s4, cv4, quantum, mw, mem)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8) :: temperature
    !> free energies, entropy and heat capacity
    real(r8), intent(out) :: df4, s4, cv4
    !> what to compute
    logical, intent(in) :: quantum
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3, egv4
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf
    !> Outer products of the eigenvectors, one column per mode, and the two
    !> intermediates of the contraction above.
    complex(r8), dimension(:, :), allocatable :: V1, V2, Wm, Sm
    !> Frequencies, bose-einstein occupation and scattering strength and some other buffer
    real(r8) :: sig, om1, om2, om3, om4, n1, n2, n3, n4, psisq, f0, f1, f2, t0, prefactor
    !>
    real(r8) :: s1, s2, s3, mult, dn1, dn2, ddn1, ddn2, df0, ddf0
    !> The complex threephonon matrix element
    complex(r8) :: c0
    !> Integers for do loops and counting
    integer :: qi, q1, q2, q3, q4, b1, b2, b3, b4, i, ctr

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(V1, [dr%n_mode**2, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(V2, [dr%n_mode**2, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(Wm, [dr%n_mode**2, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(Sm, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    t0 = walltime()

    df4 = 0.0_r8
    s4 = 0.0_r8
    cv4 = 0.0_r8

    ctr = 0
    do q1=1, qp%n_irr_point
    do q2=1, qp%n_full_point
        ctr = ctr + 1
        if (mod(ctr, mw%n) .ne. mw%r) cycle

        prefactor = qp%ip(q1)%integration_weight*qp%ap(q2)%integration_weight/uc%na
        ! pre-transform the matrix element
        call pretransform_phi4_first(fcf, qp%ip(q1)%r, qp%ap(q2)%r, ptf)
        ! Both bands at once: psisq = Re( V2^T ptf V1 ).
        !
        ! The outer products the band loop kept rebuilding are just the columns of
        ! V1 and V2, so they only need forming once per q-pair. ptf goes with
        ! its first index slowest, so it can be used here as a column-major
        ! n_mode^2 by n_mode^2 matrix without a copy.
        V1 = 0.0_r8
        V2 = 0.0_r8
        do b1 = 1, dr%n_mode
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle
            egv1 = dr%iq(q1)%egv(:, b1)/sqrt(om1)
            call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv1, 1, egv1, 1, V1(:, b1), dr%n_mode)
        end do
        do b2 = 1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle
            egv2 = dr%aq(q2)%egv(:, b2)/sqrt(om2)
            call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv2, 1, V2(:, b2), dr%n_mode)
        end do
        call zgemm('N', 'N', dr%n_mode**2, dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), &
                   ptf, dr%n_mode**2, V1, dr%n_mode**2, (0.0_r8, 0.0_r8), Wm, dr%n_mode**2)
        call zgemm('T', 'N', dr%n_mode, dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), &
                   V2, dr%n_mode**2, Wm, dr%n_mode**2, (0.0_r8, 0.0_r8), Sm, dr%n_mode)

        do b1=1, dr%n_mode
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle
            do b2=1, dr%n_mode
                om2 = dr%aq(q2)%omega(b2)
                if (om2 .lt. lo_freqtol) cycle

                psisq = real(Sm(b2, b1), r8)

                if (quantum) then
                    ! Phonon occupation
                    n1 = lo_planck(temperature, om1)
                    n2 = lo_planck(temperature, om2)
                    ! First derivative with respect to temperature
                    dn1 = lo_planck_deriv(temperature, om1)
                    dn2 = lo_planck_deriv(temperature, om2)
                    ! Second derivative with respect to temperature
                    ddn1 = lo_planck_secondderiv(temperature, om1)
                    ddn2 = lo_planck_secondderiv(temperature, om2)
                    ! The free energy
                    f0 = (2.0_r8 * n1 + 1.0_r8) * (2.0_r8 * n2 + 1.0_r8) / 32.0_r8
                    ! The entropy
                    df0 = 2.0_r8 * dn1 * (2.0_r8 * n2 + 1.0_r8) + (2.0_r8 * n1 + 1.0_r8) * 2.0_r8 * dn2
                    df0 = df0 / 32.0_r8
                    ddf0 = 2.0_r8 * ddn1 * (2.0_r8 * n2 + 1.0_r8) + 2.0_r8 * ddn2 * (2.0_r8 * n1 + 1.0_r8) + &
                         8.0_r8 * dn1 * dn2
                    ddf0 = ddf0 * temperature / 32.0_r8
                else
                    ! Much simpler in the classical case
                    ! Free energy
                    f0 = (lo_kb_Hartree*temperature)**2 / (om1*om2) / 8.0_r8
                    ! Entropy
                    df0 = lo_kb_Hartree**2 * temperature / (om1*om2) / 4.0_r8
                    ! Heat capacity
                    ddf0 = lo_kb_Hartree**2 * temperature / (om1*om2) / 4.0_r8
                end if
                ! And we accumulate
                df4 = df4 + f0 * psisq * prefactor
                s4 = s4 - df0 * psisq * prefactor
                cv4 = cv4 - ddf0 * psisq * prefactor
            end do
        end do
    end do
    if (mw%talk .and. lo_trueNtimes(q1, 127, qp%n_irr_point)) then
        call lo_progressbar(' ... fourth order', q1, qp%n_irr_point, walltime() - t0)
    end if
    end do
    ! For a nice progressbar
    if (mw%talk) then
        call lo_progressbar(' ... fourth order', qp%n_irr_point, qp%n_irr_point, walltime() - t0)
    end if
    ! Reduce on all ranks
    call mw%allreduce('sum', df4)
    call mw%allreduce('sum', s4)
    call mw%allreduce('sum', cv4)

    ! And deallocate
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(V1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(V2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(Wm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(Sm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Calculate the fourth order free energy
subroutine free_energy_fourthorder_secondorder(uc, fcf, qp, dr, temperature, fe4, s4, cv4, quantum, mw, mem)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8) :: temperature
    !> free energies, entropy and heat capacity
    real(r8), intent(out) :: fe4, s4, cv4
    !> what to compute
    logical, intent(in) :: quantum
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> For the smearing parameters
    real(r8), dimension(:, :), allocatable :: sigsq
    !> Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3, egv4
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2, evp3
    !> Phonon occupations
    real(r8), dimension(4) :: n, dn, ddn
    !> Frequencies, bose-einstein occupation and scattering strength and some other buffer
    real(r8) :: sigma, om1, om2, om3, om4, psisq, f0, f1, f2, f3, t0, prefactor
    !>
    real(r8) :: s1, s2, s3, mult, df0, ddf0, sig, sig1, sig2, sig3, sig4
    !> The phonon occupation plus 1
    real(r8) :: np1, np2, np3, np4
    !> The multiplicities
    real(r8) :: m0, m1, m2, m3
    !> The complex threephonon matrix element
    complex(r8) :: c0
    !> Integers for do loops and counting
    integer :: qi, q1, q2, q3, q4, b1, b2, b3, b4, i, ctr
    !> The dimension of the q-grid
    integer, dimension(3) :: dims

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv4, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(sigsq, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    t0 = walltime()

    select type (qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    fe4 = 0.0_r8
    s4 = 0.0_r8
    cv4 = 0.0_r8

    ! First we compute the broadening parameter for each modes
    sigsq = 0.0_r8
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            sigsq(q1, b1) = qp%smearingparameter(dr%iq(q1)%vel(:, b1), dr%default_smearing(b1), 1.0_r8)**2
        end do
    end do

    n = 0.0_r8
    dn = 0.0_r8
    ddn = 0.0_r8

    ! Same omission as free_energy_thirdorder, though nothing calls this one yet.
    ctr = 0
    do q1=1, qp%n_irr_point
    do q2=1, qp%n_full_point
    do q3=1, qp%n_full_point
        ctr = ctr + 1
        if (mod(ctr, mw%n) .ne. mw%r) cycle
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)

        if (q3 .lt. q2) cycle
        if (q4 .lt. q3) cycle

        ! Prefactors to account for equivalent points
        if (q2 .eq. q3 .and. q3 .eq. q4) then
            m0 = 1.0_r8
            m1 = 1.0_r8
            m2 = 1.0_r8
            m3 = 1.0_r8
        else if ((q2 .ne. q3 .and. q3 .eq. q4) .or. &
                 (q3 .ne. q2 .and. q2 .eq. q4) .or. &
                 (q4 .ne. q2 .and. q2 .eq. q3)) then
            m0 = 3.0_r8
            m1 = 1.0_r8
            m2 = 1.0_r8
            m3 = 1.0_r8
        else
            m0 = 6.0_r8
            m1 = 2.0_r8
            m2 = 2.0_r8
            m3 = 2.0_r8
        end if

        prefactor = qp%ip(q1)%integration_weight*qp%ap(q2)%integration_weight*qp%ap(q2)%integration_weight/uc%na
        ! pre-transform the matrix element
        call pretransform_phi4(fcf, qp%ap(q2)%r, qp%ap(q3)%r, qp%ap(q4)%r, ptf)

        do b1=1, dr%n_mode
            ! Get first phonon
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle
            n(1) = lo_planck(temperature, om1)
            dn(1) = lo_planck_deriv(temperature, om1)
            ddn(1) = lo_planck_secondderiv(temperature, om1)
            egv1 = dr%iq(q1)%egv(:, b1)/sqrt(om1)
            do b2=1, dr%n_mode
                ! Get second phonon
                om2 = dr%aq(q2)%omega(b2)
                if (om2 .lt. lo_freqtol) cycle
                n(2) = lo_planck(temperature, om2)
                dn(2) = lo_planck_deriv(temperature, om2)
                ddn(2) = lo_planck_secondderiv(temperature, om2)
                egv2 = dr%aq(q2)%egv(:, b2)/sqrt(om2)

                ! Multiply first and second phonon
                evp1 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                do b3=1, dr%n_mode
                    ! Get third phonon
                    om3 = dr%aq(q3)%omega(b3)
                    if (om3 .lt. lo_freqtol) cycle
                    n(3) = lo_planck(temperature, om3)
                    dn(3) = lo_planck_deriv(temperature, om3)
                    ddn(3) = lo_planck_secondderiv(temperature, om3)
                    egv3 = dr%aq(q3)%egv(:, b3)/sqrt(om3)

                    ! Project on third phonon
                    evp2 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                    do b4=1, dr%n_mode
                        ! Get fourth phonon
                        om4 = dr%aq(q4)%omega(b4)
                        if (om4 .lt. lo_freqtol) cycle
                        n(4) = lo_planck(temperature, om4)
                        dn(4) = lo_planck_deriv(temperature, om4)
                        ddn(4) = lo_planck_secondderiv(temperature, om4)
                        egv4 = dr%aq(q4)%egv(:, b4)/sqrt(om4)

                        ! Project on fourth phonon
                        evp3 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                        evp3 = conjg(evp3)

                        ! Compute the scattering matrix element
                        c0 = dot_product(evp3, ptf)
                        psisq = abs(c0*conjg(c0))*prefactor

                        if (quantum) then
                            ! Phonon occupation
                            np1 = n(1) + 1.0_r8
                            np2 = n(2) + 1.0_r8
                            np3 = n(3) + 1.0_r8
                            np4 = n(4) + 1.0_r8

                            ! We get the broadening parameter
                            sig1 = sigsq(q1, b1)
                            sig2 = sigsq(qp%ap(q2)%irreducible_index, b2)
                            sig3 = sigsq(qp%ap(q3)%irreducible_index, b3)
                            sig4 = sigsq(qp%ap(q4)%irreducible_index, b3)
                            sig = sqrt(sig1 + sig2 + sig3 + sig4)

                            f1 = (np1*np2*np3*np4 - n(1)*n(2)*n(3)*n(4)) / real(1.0_r8/(om1+om2+om3+om4+lo_imag*sig), r8)
                            f2 = 4.0_r8*(n(1)*np2*np3*np4 - np1*n(2)*n(3)*n(4)) * real(1.0/(-om1+om2+om3+om4 + lo_imag*sig), r8)
                            f3 = 3.0_r8*(n(1)*n(2)*(n(3)+n(4)+1.0_r8) - n(3)*n(4)*(n(1)+n(2)+1.0_r8))*real(1.0/(om1+om2-om3-om4 + lo_imag*sig), r8)

                            f0 = (f1 + f2 + f3) / 768.0_r8
                        else
                            ! Much simpler in the classical case
                            ! Free energy
                            f0 = (lo_kb_Hartree*temperature)**3 / (om1*om2*om3*om4) / 48.0_r8
                            ! Entropy
                            df0 = lo_kb_Hartree**3*temperature**2 / (om1*om2*om3*om4) / 16.0_r8
                            ! Heat capacity
                            ddf0 = lo_kb_Hartree**3*temperature**2 / (om1*om2*om3*om4) / 8.0_r8
                        end if
                        ! And we accumulate
                        fe4 = fe4 - f0 * psisq
                        s4 = s4 + df0 * psisq
                        cv4 = cv4 + ddf0 * psisq
                    end do
                end do
            end do
        end do
    end do
    end do
    if (mw%talk .and. lo_trueNtimes(q1, 127, qp%n_irr_point)) then
        call lo_progressbar(' ... fourth order', q1, qp%n_irr_point, walltime() - t0)
    end if
    end do
    ! For a nice progressbar
    if (mw%talk) then
        call lo_progressbar(' ... fourth order', qp%n_irr_point, qp%n_irr_point, walltime() - t0)
    end if
    ! Reduce on all ranks
    call mw%allreduce('sum', fe4)
    call mw%allreduce('sum', s4)
    call mw%allreduce('sum', cv4)

    ! And deallocate memory
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Get the Fourier transform of the fourth order matrix element
subroutine pretransform_phi4_first(fcf, q1, q2, ptf)
    !> fourth order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q1, q2
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l, m

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3, rv4
    real(r8) :: iqr
    integer :: a1, a2, a3, a4, ia, ib, ic, id, q, nb

    nb = fcf%na*3
    ptf = 0.0_r8
    do a1 = 1, fcf%na
    do q = 1, fcf%atom(a1)%n
        a2 = fcf%atom(a1)%quartet(q)%i2
        a3 = fcf%atom(a1)%quartet(q)%i3
        a4 = fcf%atom(a1)%quartet(q)%i4

        rv2 = fcf%atom(a1)%quartet(q)%lv2
        rv3 = fcf%atom(a1)%quartet(q)%lv3
        rv4 = fcf%atom(a1)%quartet(q)%lv4

        iqr = -dot_product(q1, rv2) + dot_product(q2, rv3) - dot_product(q2, rv4)
        iqr = iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do l = 1, 3
        do k = 1, 3
        do j = 1, 3
        do i = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            id = (a4 - 1)*3 + l
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            m = (ia - 1)*nb*nb*nb + (ib - 1)*nb*nb + (ic - 1)*nb + id
            ptf(m) = ptf(m) + fcf%atom(a1)%quartet(q)%mwm(i, j, k, l)*expiqr
        end do
        end do
        end do
        end do
    end do
    end do
end subroutine

!> Do half the transform to q-space
pure subroutine pretransform_phi4(fcf, q2, q3, q4, ptf)
    !> Fourth order forceconstant
    class(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> The q-vectors, without the 2pi
    real(r8), dimension(3), intent(in) :: q2, q3, q4
    !> Flattened, pre-transformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l, m

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3, rv4
    real(r8) :: iqr
    integer :: a1, a2, a3, a4, ia, ib, ic, id, q, nb

    nb = fcf%na*3
    ptf = 0.0_r8
    do a1 = 1, fcf%na
    do q = 1, fcf%atom(a1)%n
        a2 = fcf%atom(a1)%quartet(q)%i2
        a3 = fcf%atom(a1)%quartet(q)%i3
        a4 = fcf%atom(a1)%quartet(q)%i4

        rv2 = fcf%atom(a1)%quartet(q)%lv2
        rv3 = fcf%atom(a1)%quartet(q)%lv3
        rv4 = fcf%atom(a1)%quartet(q)%lv4

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3) + dot_product(q4, rv4)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do l = 1, 3
        do k = 1, 3
        do j = 1, 3
        do i = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            id = (a4 - 1)*3 + l
            ! Now for the grand flattening scheme, consistent with the zgeru operations.
            m = (ia - 1)*nb*nb*nb + (ib - 1)*nb*nb + (ic - 1)*nb + id
            ptf(m) = ptf(m) + fcf%atom(a1)%quartet(q)%mwm(i, j, k, l)*expiqr
        end do
        end do
        end do
        end do
    end do
    end do
end subroutine
end module
