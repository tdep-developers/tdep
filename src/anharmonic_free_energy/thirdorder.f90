module thirdorder
!! Compute the third order contribution to the anharmonic free energy
use konstanter, only: r8, lo_freqtol, lo_twopi, lo_imag, lo_kb_Hartree, lo_exitcode_param
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_planck, &
                          lo_planck_deriv, lo_planck_secondderiv, lo_invert_real_matrix, &
                          lo_harmonic_oscillator_cv, lo_harmonic_oscillator_internal_energy
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_symmetryoperation, only: lo_operate_on_secondorder_tensor
use lo_timetracker, only: lo_timer
use lo_fftgrid_helper, only: fft_third_grid_index

use thermodynamic_helpers, only: full_to_voigt

implicit none
private

public :: elastic_thirdorder
public :: free_energy_thirdorder

contains
subroutine elastic_thirdorder(uc, fc, fct, qp, dr, temperature, p3ph, alpha, quantum, mw, mem)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order force constants
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8), intent(in) :: temperature
    !> free energies and heat capacity
    real(r8), dimension(3, 3), intent(out) :: p3ph, alpha
    !> what to compute
    logical, intent(in) :: quantum
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> The eigenvector for the current mode
    complex(r8), dimension(:), allocatable :: egv
    !> Complex coefficients for gruneisen
    complex(r8), dimension(3, 3) :: c0
    !> stiffness and compliance tensors, in voigt notation
    real(r8), dimension(6, 6) :: cij, sij
    !> Buffers for Gruneisen and thermal expansion, with full anisotropy
    real(r8), dimension(3, 3) :: grun, bufalpha, grunbuf
    !> q-points and lattice distances
    real(r8), dimension(3) :: qv, rv2, rv3
    !> Phase factor for Fourier transform and product of eigenvectors
    complex(r8) :: expiqr, evprod
    !> Harmonic stuffs
    real(r8) :: om, iqr, invsqrtmass, u0, cv0
    !> Integers for do loops
    integer :: q1, b1, a1, t, a2, a3, i, j, k, l, ivoigt, jvoigt, ctr, ifull

    ! We start by allocating space for the eigenvector
    call mem%allocate(egv, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Set everything to zero
    p3ph = 0.0_r8
    bufalpha = 0.0_r8

    ! And here we go
    ctr = 0
    do q1=1, qp%n_irr_point
        qv = qp%ip(q1)%r * lo_twopi
        do b1=1, dr%n_mode
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! We start by computing the mode Gruneisen parameter, with full anisotropy
            om = dr%iq(q1)%omega(b1)
            if (om .lt. lo_freqtol) cycle
            c0 = 0.0_r8
            egv = dr%iq(q1)%egv(:, b1)
            do a1=1, fct%na
            do t=1, fct%atom(a1)%n
                a2 = fct%atom(a1)%triplet(t)%i2
                a3 = fct%atom(a1)%triplet(t)%i3
                rv2 = fct%atom(a1)%triplet(t)%lv2
                rv3 = fct%atom(a1)%triplet(t)%rv3
                invsqrtmass = uc%invsqrtmass(a1) * uc%invsqrtmass(a2)
                iqr = -dot_product(qv, rv2)
                expiqr = cmplx(cos(iqr), sin(iqr), r8)
                do i=1, 3
                do j=1, 3
                    evprod = conjg(egv((a1-1) * 3 + i)) * egv((a2-1) * 3 + j) * expiqr
                    do k=1, 3
                    do l=1, 3
                        c0(k, l) = c0(k, l) + fct%atom(a1)%triplet(t)%m(i, j, k) * evprod * rv3(l) * invsqrtmass
                    end do
                    end do
                end do
                end do
            end do
            end do
            grun = -real(c0, r8) / (6.0_r8 * om**2)

            ! Here we sum over all symmetry equivalent q-point, taking into account the symop
            grunbuf = 0.0_r8
            do ifull=1, qp%ip(q1)%n_full_point
                grunbuf = grunbuf + lo_operate_on_secondorder_tensor(uc%sym%op(ifull), grun)
            end do
            grun = grunbuf

            ! Now we can accumulate things using the Gruneisen parameter
            if (quantum) then
                u0 = lo_harmonic_oscillator_internal_energy(temperature, om)
                cv0 = lo_harmonic_oscillator_cv(temperature, om)
            else
                u0 = lo_kb_Hartree * temperature
                cv0 = lo_kb_Hartree
            end if

            ! And we accumulate
            p3ph(:, :) = p3ph(:, :) - grun(:, :) * u0 / uc%volume / qp%n_full_point
            bufalpha(:, :) = bufalpha(:, :) + cv0 * grun(:, :) / qp%n_full_point
        end do
    end do
    ! Reduce the calculation
    call mw%allreduce('sum', p3ph)
    call mw%allreduce('sum', bufalpha)

    ! We need the compliance tensor to get the thermal expansion
    call fc%get_elastic_constants(uc, mw, -1)
    cij = fc%elastic_constants_voigt + fc%elastic_constants_voigt_longrange
    call lo_invert_real_matrix(cij, sij)

    ! And we just have to add the compliance tensor part to get the thermal expansion
    alpha = 0.0_r8
    do i=1, 3
    do j=1, 3
    do k=1, 3
    do l=1, 3
        ivoigt = full_to_voigt(i, j)
        jvoigt = full_to_voigt(k, l)

        ! That's simply the formula
        alpha(i, j) = alpha(i, j) + sij(ivoigt, jvoigt) * bufalpha(k, l) / uc%volume
    end do
    end do
    end do
    end do
    call mem%deallocate(egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Calculate the third order free energy
subroutine free_energy_thirdorder(uc, fct, qp, dr, temperature, fe3, s3, cv3, quantum, mw, mem)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8) :: temperature
    !> free energies and heat capacity
    real(r8), intent(out) :: fe3, s3, cv3
    !> what to compute
    logical, intent(in) :: quantum
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> For the smearing parameters
    real(r8), dimension(:, :), allocatable :: sigsq
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf
    !> Frequency-scaled eigenvectors of the three q-points, one column per mode,
    !> the two intermediates of the contraction, and every matrix element of
    !> this q-pair
    complex(r8), dimension(:, :), allocatable :: E1, E2, E3, T1, T2, Cm
    complex(r8), dimension(:, :, :), allocatable :: c0all
    !> Frequencies, bose-einstein occupation and scattering strength and some other buffer
    real(r8) :: sigma, om1, om2, om3, n2, n3, psisq, f0, f1, plf0, plf1, perm, f2, n1, pref, t0, prefactor
    !>
    real(r8) :: sig1, sig2, sig3, dn1, dn2, dn3, ddn1, ddn2, ddn3, df1, df2, ddf1, ddf2, df0, ddf0
    !> The complex threephonon matrix element
    complex(r8) :: c0
    !> Integers for do loops and counting
    integer :: qi, q1, q2, q3, q2p, q3p, b1, b2, b3, i, ctr
    !> The dimension of the q-grid
    integer, dimension(3) :: dims

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(E1, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(E2, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(E3, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(T1, [dr%n_mode**2, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(T2, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(Cm, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(c0all, [dr%n_mode, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(sigsq, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    t0 = walltime()

    select type (qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    fe3 = 0.0_r8
    s3 = 0.0_r8
    cv3 = 0.0_r8

    ! First we compute the broadening parameter for each modes
    sigsq = 0.0_r8
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            sigsq(q1, b1) = qp%smearingparameter(dr%iq(q1)%vel(:, b1), dr%default_smearing(b1), 1.0_r8)**2
        end do
    end do

    if (mw%talk) call lo_progressbar_init()

    ! mod(ctr, mw%n) == mw%r only hands each q-pair to one rank if every rank is
    ! counting from the same place. Uninitialised they are not: one 64 rank run
    ! here started them anywhere between 32764 and 32767, which leaves some pairs
    ! summed twice and others not at all. Extremely small differences.
    ctr = 0
    do q1=1, qp%n_irr_point
    do q2=1, qp%n_full_point
        ctr = ctr + 1
        if (mod(ctr, mw%n) .ne. mw%r) cycle
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)
       

        prefactor = qp%ip(q1)%integration_weight*qp%ap(q2)%integration_weight/uc%na
        ! pre-transform the matrix element
        call pretransform_phi3(fct, qp%ap(q2)%r, qp%ap(q3)%r, ptf)

        ! All three bands at once, peeling one index off ptf at a time.
        !
        ! Two things to be careful of. There is no conjugation left in the sum:
        ! the conjg on evp2 and the one dot_product applies to its first argument
        ! cancel, so putting either back changes the physics. And ptf comes out
        ! with its first index slowest, which is what lets it be used directly as
        ! a column-major n_mode^2 by n_mode matrix here.
        E1 = 0.0_r8
        E2 = 0.0_r8
        E3 = 0.0_r8
        do b1 = 1, dr%n_mode
            if (dr%iq(q1)%omega(b1) .gt. lo_freqtol) &
                E1(:, b1) = dr%iq(q1)%egv(:, b1)/sqrt(dr%iq(q1)%omega(b1))
            if (dr%aq(q2)%omega(b1) .gt. lo_freqtol) &
                E2(:, b1) = dr%aq(q2)%egv(:, b1)/sqrt(dr%aq(q2)%omega(b1))
            if (dr%aq(q3)%omega(b1) .gt. lo_freqtol) &
                E3(:, b1) = dr%aq(q3)%egv(:, b1)/sqrt(dr%aq(q3)%omega(b1))
        end do
        call zgemm('N', 'N', dr%n_mode**2, dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), &
                   ptf, dr%n_mode**2, E1, dr%n_mode, (0.0_r8, 0.0_r8), T1, dr%n_mode**2)
        do b1 = 1, dr%n_mode
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
            call zgemm('N', 'N', dr%n_mode, dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), &
                       T1(:, b1), dr%n_mode, E2, dr%n_mode, (0.0_r8, 0.0_r8), T2, dr%n_mode)
            call zgemm('T', 'N', dr%n_mode, dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), &
                       E3, dr%n_mode, T2, dr%n_mode, (0.0_r8, 0.0_r8), Cm, dr%n_mode)
            c0all(:, :, b1) = Cm
        end do

        do b1=1, dr%n_mode
            ! Get first phonon
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle
            n1 = lo_planck(temperature, om1)
            dn1 = lo_planck_deriv(temperature, om1)
            ddn1 = lo_planck_secondderiv(temperature, om1)
            do b2=1, dr%n_mode
                ! Get second phonon
                om2 = dr%aq(q2)%omega(b2)
                if (om2 .lt. lo_freqtol) cycle
                n2 = lo_planck(temperature, om2)
                dn2 = lo_planck_deriv(temperature, om2)
                ddn2 = lo_planck_secondderiv(temperature, om2)

                do b3=1, dr%n_mode
                    ! Get third phonon
                    om3 = dr%aq(q3)%omega(b3)
                    if (om3 .lt. lo_freqtol) cycle
                    n3 = lo_planck(temperature, om3)
                    dn3 = lo_planck_deriv(temperature, om3)
                    ddn3 = lo_planck_secondderiv(temperature, om3)

                    c0 = c0all(b3, b2, b1)
                    psisq = real(conjg(c0)*c0, r8) * prefactor

                    if (quantum) then
                        ! Get a smearing parameter
                        sig1 = sigsq(q1, b1)
                        sig2 = sigsq(qp%ap(q2)%irreducible_index, b2)
                        sig3 = sigsq(qp%ap(q3)%irreducible_index, b3)
                        sigma = sqrt(sig1 + sig2 + sig3)

                        ! This are the formula to get the free energy
                        f1 = (n1 + 1.0_r8)*(n2 + n3 + 1.0_r8) + n2*n3
                        f2 = n3*(n1 + n2 + 1.0_r8) - n1*n2
                        ! Here is for the entropy
                        df1 = dn1 * (n2 + n3 + 1.0_r8) + (n1 + 1.0_r8) * (dn2 + dn3) + (dn2 * n3) + (n2 * dn3)
                        df2 = dn3 * (n1 + n2 + 1.0_r8) + (dn1 * (n3 - n2)) + (dn2 * (n3 - n1))
                        ! And the heat capacity
                        ddf1 = ddn1 * (n2 + n3 + 1.0_r8) + 2.0_r8 * dn1 * (dn2 + dn3) + (n1 + 1.0_r8) * (ddn2 + ddn3) + &
                               ddn2 * n3 + n2 * ddn3 + 2.0_r8 * dn2 * dn3
                        ddf2 = (ddn3 * (n1 + n2 + 1.0_r8)) + (dn3 * (dn1 + dn2)) + &
                               (ddn1 * (n3 - n2)) + (dn1 * (dn3 - dn2)) + &
                               (ddn2 * (n3 - n1)) + (dn2 * (dn3 - dn1))

                        ! Compute the principal values
                        f1 = f1*real(1.0_r8/(om1 + om2 + om3 + lo_imag*sigma), r8)
                        df1 = df1*real(1.0_r8/(om1 + om2 + om3 + lo_imag*sigma), r8)
                        ddf1 = ddf1*real(1.0_r8/(om1 + om2 + om3 + lo_imag*sigma), r8)
                        f2 = 3.0_r8*f2*real(1.0/(om1 + om2 - om3 + lo_imag*sigma), r8)
                        df2 = 3.0_r8 * df2 * real(1.0/(om1 + om2 - om3 + lo_imag*sigma), r8)
                        ddf2 = 3.0_r8 * ddf2 * real(1.0/(om1 + om2 - om3 + lo_imag*sigma), r8)

                        ! Get everything together
                        f0 = (f1 + f2) / 48.0_r8
                        df0 = (df1 + df2) / 48.0_r8
                        ddf0 = (ddf1 + ddf2) * temperature / 48.0_r8
                    else
                        ! It's much simpler in the classical case
                        ! Free energy
                        f0 = (lo_kb_Hartree*temperature)**2 / (om1*om2*om3) / 12.0_r8
                        ! Entropy
                        df0 = lo_kb_Hartree**2 * temperature / (om1*om2*om3) / 6.0_r8
                        ! Heat capacity
                        ddf0 = lo_kb_Hartree**2 * temperature / (om1*om2*om3) / 6.0_r8
                    end if
                    ! And accumulate
                    fe3 = fe3 - f0 * psisq
                    s3 = s3 + df0 * psisq
                    cv3 = cv3 + ddf0 * psisq
                end do
            end do
        end do
    end do
    if (mw%talk .and. lo_trueNtimes(q1, 127, qp%n_irr_point)) then
        call lo_progressbar(' ... third order', q1, qp%n_irr_point, walltime() - t0)
    end if
    end do
    if (mw%talk) then
        call lo_progressbar(' ... third order', qp%n_irr_point, qp%n_irr_point, walltime() - t0)
    end if
    ! Reduce on all ranks
    call mw%allreduce('sum', fe3)
    call mw%allreduce('sum', s3)
    call mw%allreduce('sum', cv3)

    ! And we can deallocate everything
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(E1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(E2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(E3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(T1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(T2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(Cm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(c0all, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(sigsq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Get the Fourier transform of the third order matrix element
subroutine pretransform_phi3(fct, q2, q3, ptf)
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2, q3
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3
    real(r8) :: iqr
    integer :: a1, a2, a3, ia, ib, ic, t, nb

    nb = fct%na*3
    ptf = 0.0_r8
    do a1 = 1, fct%na
    do t = 1, fct%atom(a1)%n
        a2 = fct%atom(a1)%triplet(t)%i2
        a3 = fct%atom(a1)%triplet(t)%i3

        rv2 = fct%atom(a1)%triplet(t)%lv2
        rv3 = fct%atom(a1)%triplet(t)%lv3

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            l = (ia - 1)*nb*nb + (ib - 1)*nb + ic
            ptf(l) = ptf(l) + fct%atom(a1)%triplet(t)%mwm(i, j, k)*expiqr
        end do
        end do
        end do
    end do
    end do
end subroutine
end module
