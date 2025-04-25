submodule(lo_evaluate_phonon_self_energy) lo_evaluate_phonon_self_energy_isotope
implicit none
contains

!> evaluate the isotope self-energy
module subroutine isotope_imaginary_selfenergy(wp,se,qp,dr,sr,mw,mem,verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> scattering rates
    type(lo_listofscatteringrates), intent(in) :: sr
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    select case (se%integrationtype)
    case (1:2)
        call isotope_imaginary_selfenergy_gaussian(wp, qp, dr, se, sr, mw, mem, verbosity)
    case (3)
        call isotope_imaginary_selfenergy_tetrahedron(wp, qp, dr, se, sr, mw, mem, verbosity)
    case (4:5)
        ! Nothing, taken care of with the other integrals
    case default
        call lo_stop_gracefully(['Unknown integration type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end select

end subroutine

!> Gaussian integration of the isotope self energy
subroutine isotope_imaginary_selfenergy_gaussian(wp, qp, dr, se, sr, mw, mem, verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> scattering rates
    type(lo_listofscatteringrates), intent(in) :: sr
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:), allocatable :: buf, sigma
    real(r8), dimension(3) :: v0
    real(r8) :: prefactor, invf, t0
    integer :: q, b1, b2, i, ii, jj, ctr

    ! Some constants
    t0 = walltime()
    invf = se%n_energy/se%energy_axis(se%n_energy)
    ! Some space
    call mem%allocate(sigma, se%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (verbosity .gt. 0) call lo_progressbar_init()
    ! start integrating
    ctr = 0
    se%im_iso = 0.0_r8
    if (sr%atgamma) then
        do q = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! MPI thing
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! Get the smearing parameter, if adaptive
            select case (se%integrationtype)
            case (1)
                do b2 = 1, dr%n_mode
                    sigma(b2) = dr%default_smearing(b2)*se%smearing_prefactor
                end do
            case (2)
                do b2 = 1, dr%n_mode
                    v0 = dr%iq(q)%vel(:, b2)
                    sigma(b2) = qp%smearingparameter(v0, dr%default_smearing(b2), se%smearing_prefactor)
                end do
            end select
            ! Get the prefactor
            prefactor = isotope_prefactor*qp%ip(q)%integration_weight

            ! sum things up
            do b2 = 1, dr%n_mode
                ! Add it to the self-energy
                ii = max(floor((dr%iq(q)%omega(b2) - 4*sigma(b2))*invf), 1)
                jj = min(ceiling((dr%iq(q)%omega(b2) + 4*sigma(b2))*invf), se%n_energy)
                do i = ii, jj
                    se%im_iso(i, b1) = se%im_iso(i, b1) + prefactor*lo_gauss(se%energy_axis(i), dr%iq(q)%omega(b2), sigma(b2))*sr%psi_iso(b1, b2, q)
                end do
            end do

            ! Report?
            if (verbosity .gt. 0) then
                if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) call lo_progressbar(' ... isotope imaginary selfenergy', ctr, dr%n_mode*qp%n_irr_point)
            end if
        end do
        end do
    else
        do q = 1, qp%n_full_point
        do b1 = 1, dr%n_mode
            ! MPI thing
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! Get the smearing parameter, if adaptive
            select case (se%integrationtype)
            case (1)
                do b2 = 1, dr%n_mode
                    sigma(b2) = dr%default_smearing(b2)*se%smearing_prefactor
                end do
            case (2)
                do b2 = 1, dr%n_mode
                    v0 = dr%aq(q)%vel(:, b2)
                    sigma(b2) = qp%smearingparameter(v0, dr%default_smearing(b2), se%smearing_prefactor)
                end do
            end select
            ! Get the prefactor
            prefactor = isotope_prefactor*qp%ap(q)%integration_weight

            ! sum things up
            do b2 = 1, dr%n_mode
                ! Add it to the self-energy
                ii = max(floor((dr%aq(q)%omega(b2) - 4*sigma(b2))*invf), 1)
                jj = min(ceiling((dr%aq(q)%omega(b2) + 4*sigma(b2))*invf), se%n_energy)
                do i = ii, jj
                    se%im_iso(i, b1) = se%im_iso(i, b1) + prefactor*lo_gauss(se%energy_axis(i), dr%aq(q)%omega(b2), sigma(b2))*sr%psi_iso(b1, b2, q)
                end do
            end do

            ! Report?
            if (verbosity .gt. 0) then
                if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) call lo_progressbar(' ... isotope imaginary selfenergy', ctr, dr%n_mode*qp%n_full_point)
            end if
        end do
        end do
    end if
    ! Sum it up
    call mw%allreduce('sum', se%im_iso)

    ! Fix degeneracies
    call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf = 0.0_r8
    do b1 = 1, dr%n_mode
        buf = 0.0_r8
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            buf = buf + se%im_iso(:, b2)
        end do
        buf = buf/real(wp%degeneracy(b1), r8)
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            se%im_iso(:, b2) = buf
        end do
    end do
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(sigma, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (verbosity .gt. 0) call lo_progressbar(' ... isotope imaginary selfenergy', dr%n_mode*qp%n_full_point, dr%n_mode*qp%n_full_point, walltime() - t0)
end subroutine

!> Integrate the isotope part with the tetrahedron method
subroutine isotope_imaginary_selfenergy_tetrahedron(wp, qp, dr, se, sr, mw, mem, verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> the q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> the dipersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> the resulting self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> the scattering rates
    type(lo_listofscatteringrates), intent(in) :: sr
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(4) :: bvals, cvals, wts1, wts2
    real(r8), dimension(:, :, :), allocatable :: psisq
    real(r8), dimension(:, :), allocatable :: omr2
    real(r8), dimension(:), allocatable :: buf
    real(r8) :: minc, maxc, omthres, delta, int1, int2, bigOM, invf, t0, sigma
    integer :: i, t, q, ii, jj, b1, b2

    ! Start timer
    t0 = walltime()

    ! Make a little space
    call mem%allocate(psisq, [4, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr2, [4, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    psisq = 0.0_r8
    omr2 = 0.0_r8

    ! Threshold for locating acoustic modes
    omthres = dr%omega_min*0.5_r8
    ! Miniscule smearing in the integration
    delta = (se%energy_axis(2) - se%energy_axis(1))*0.25_r8
    ! Precalculated thing to make integrations faster
    invf = se%n_energy/se%energy_axis(se%n_energy)

    ! Start integrating
    if (verbosity .gt. 0) call lo_progressbar_init()
    se%im_iso = 0.0_r8
    do t = 1, qp%n_full_tet
        ! mpi stuff
        if (mod(t, mw%n) .ne. mw%r) cycle
        ! grab frequencies and matrix elements for this tetrahedron
        do ii = 1, 4
            q = qp%at(t)%full_index(ii)
            omr2(ii, :) = dr%aq(q)%omega
            do b1 = 1, dr%n_mode
            do b2 = 1, dr%n_mode
                psisq(ii, b1, b2) = sr%psi_iso(b1, b2, q)
            end do
            end do
        end do
        ! Add the prefactors directly to psisq
        psisq = psisq*qp%at(t)%integration_weight*isotope_prefactor
        ! loop over frequencies and bands
        do b2 = 1, dr%n_mode
            ! smearing for degenerate tetrahedrons

            ! get tetrahedron values
            minc = lo_huge
            maxc = -lo_huge
            do ii = 1, 4
                cvals(ii) = omr2(ii, b2)
                minc = min(minc, cvals(ii))
                maxc = max(maxc, cvals(ii))
            end do
            minc = minc - 2*delta
            maxc = maxc + 2*delta
            ii = max(floor(minc*invf), 1)
            jj = min(ceiling(maxc*invf), se%n_energy)
            do b1 = 1, dr%n_mode
                sigma = dr%default_smearing(b1)
                do i = ii, jj
                    bigOM = se%energy_axis(i)
                    ! is it relevant?
                    if (bigOM .gt. minc .and. bigOM .lt. maxc) then
                        bvals = psisq(:, b1, b2)
                        wts1 = lo_LV_tetrahedron_weights(cvals, bigOM - delta, lo_freqtol, sigma)
                        wts2 = lo_LV_tetrahedron_weights(cvals, bigOM + delta, lo_freqtol, sigma)
                        int1 = sum(wts1*bvals)
                        int2 = sum(wts2*bvals)
                        se%im_iso(i, b1) = se%im_iso(i, b1) + (int1 + int2)*0.5_r8
                    end if
                end do
            end do
        end do
        if (verbosity .gt. 0) then
            if (lo_trueNtimes(t, 63, qp%n_full_tet)) call lo_progressbar(' ... isotope imaginary selfenergy', t, qp%n_full_tet)
        end if
    end do
    ! Clean a little
    call mem%deallocate(psisq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ! Sync
    call mw%allreduce('sum', se%im_iso)
    ! Fix degeneracies
    call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf = 0.0_r8
    ! Fix degeneracies
    do b1 = 1, dr%n_mode
        buf = 0.0_r8
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            buf = buf + se%im_iso(:, b2)
        end do
        buf = buf/real(wp%degeneracy(b1), r8)
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            se%im_iso(:, b2) = buf
        end do
    end do
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Report that we are done!
    if (verbosity .gt. 0) call lo_progressbar(' ... isotope imaginary selfenergy', se%n_energy, se%n_energy, walltime() - t0)
end subroutine

end submodule