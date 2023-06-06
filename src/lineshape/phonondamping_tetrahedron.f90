
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

    real(r8), parameter :: isotope_prefactor = lo_pi/4.0_r8
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

!> The three-phonon imaginary selfenergy with the tetrahedron method
subroutine threephonon_imaginary_selfenergy_tetrahedron(wp, qp, sr, dr, temperature, se, mw, mem, verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> the q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> the phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> the scattering rates
    type(lo_listofscatteringrates), intent(inout) :: sr
    !> the self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> the temperature
    real(r8), intent(in) :: temperature
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:, :, :, :), allocatable :: psisq1, psisq2
    real(r8), dimension(:, :), allocatable :: omr2, omr3
    real(r8), dimension(:), allocatable :: buf
    real(r8), dimension(4) :: wts1, wts2, bvals, cvals
    real(r8) :: minc, maxc, sigma, omthres, f0, f1, delta
    real(r8) :: bigOM, invf, t0
    integer :: i, j, ii, jj, q
    integer :: b1, b2, b3

    ! Fix the third order stuff
    t0 = walltime()
    delta = (se%energy_axis(2) - se%energy_axis(1))*0.1_r8
    invf = se%n_energy/se%energy_axis(se%n_energy)
    omthres = dr%omega_min*0.5_r8

    call mem%allocate(omr2, [4, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr3, [4, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(psisq1, [4, dr%n_mode, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(psisq2, [4, dr%n_mode, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    omr2 = 0.0_r8
    omr3 = 0.0_r8
    psisq1 = 0.0_r8
    psisq2 = 0.0_r8

    se%im_3ph = 0.0_r8
    if (verbosity .gt. 0) call lo_progressbar_init()
    do j = 1, qp%n_full_tet
        ! mpi stuff
        if (mod(j, mw%n) .ne. mw%r) cycle
        ! fetch frequencies and matrix elements for this tetrahedron
        do ii = 1, 4
            q = qp%at(j)%full_index(ii)
            ! get frequencies in tetrahedron
            omr2(ii, :) = dr%aq(q)%omega
            omr3(ii, :) = sr%omega3(:, q)
            ! matrix elements, weighted with occupation functions
            do b1 = 1, dr%n_mode
            do b2 = 1, dr%n_mode
            do b3 = 1, dr%n_mode
                if (omr2(ii, b2) .gt. omthres .and. omr3(ii, b3) .gt. omthres) then
                    f1 = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))
                    f0 = 1.0_r8 + lo_planck(temperature, omr2(ii, b2)) + lo_planck(temperature, omr3(ii, b3))
                    psisq1(ii, b1, b2, b3) = f1*f0
                    f0 = lo_planck(temperature, omr3(ii, b3)) - lo_planck(temperature, omr2(ii, b2))
                    psisq2(ii, b1, b2, b3) = f1*f0
                else
                    psisq1(ii, b1, b2, b3) = 0.0_r8
                    psisq2(ii, b1, b2, b3) = 0.0_r8
                end if
            end do
            end do
            end do

            ! Adjusted version
            do b1 = 1, dr%n_mode
            do b2 = 1, dr%n_mode
            do b3 = 1, dr%n_mode
                if (b2 .eq. b3) then
                    if (omr2(ii, b2) .gt. omthres .and. omr3(ii, b3) .gt. omthres) then
                        f1 = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))
                        f0 = 1.0_r8 + lo_planck(temperature, omr2(ii, b2)) + lo_planck(temperature, omr3(ii, b3))
                        psisq1(ii, b1, b2, b3) = f1*f0
                        f0 = lo_planck(temperature, omr3(ii, b3)) - lo_planck(temperature, omr2(ii, b2))
                        psisq2(ii, b1, b2, b3) = f1*f0
                    else
                        psisq1(ii, b1, b2, b3) = 0.0_r8
                        psisq2(ii, b1, b2, b3) = 0.0_r8
                    end if
                else
                    if (omr2(ii, b2) .gt. omthres .and. omr3(ii, b3) .gt. omthres) then
                        f1 = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q))) + abs(sr%psi_3ph(b1, b3, b2, q)*conjg(sr%psi_3ph(b1, b3, b2, q)))
                        f0 = 1.0_r8 + lo_planck(temperature, omr2(ii, b2)) + lo_planck(temperature, omr3(ii, b3))
                        psisq1(ii, b1, b2, b3) = f1*f0
                        f0 = lo_planck(temperature, omr3(ii, b3)) - lo_planck(temperature, omr2(ii, b2))
                        psisq2(ii, b1, b2, b3) = f1*f0
                    else
                        psisq1(ii, b1, b2, b3) = 0.0_r8
                        psisq2(ii, b1, b2, b3) = 0.0_r8
                    end if
                end if
            end do
            end do
            end do
        end do
        ! And add the prefactor
        psisq1 = psisq1*threephonon_prefactor*qp%at(j)%integration_weight
        psisq2 = psisq2*threephonon_prefactor*qp%at(j)%integration_weight

        ! Sum over bands and frequencies
        do b2 = 1, dr%n_mode
            !do b3=1,dr%n_mode
            do b3 = b2, dr%n_mode
                ! delta(bigOM-om2-om3)
                minc = lo_huge
                maxc = -lo_huge
                do ii = 1, 4
                    cvals(ii) = omr2(ii, b2) + omr3(ii, b3)
                    minc = min(minc, cvals(ii))
                    maxc = max(maxc, cvals(ii))
                end do
                ii = max(floor(minc*invf), 1)
                jj = min(ceiling(maxc*invf), se%n_energy)
                do b1 = 1, dr%n_mode
                    sigma = dr%default_smearing(b1)
                    bvals = psisq1(:, b1, b2, b3)
                    do i = ii, jj
                        bigOM = se%energy_axis(i)
                        wts1 = lo_LV_tetrahedron_weights(cvals, bigOM - delta, lo_freqtol, sigma)
                        wts2 = lo_LV_tetrahedron_weights(cvals, bigOM + delta, lo_freqtol, sigma)
                        se%im_3ph(i, b1) = se%im_3ph(i, b1) + sum((wts1 + wts2)*bvals)*0.5_r8
                    end do
                end do

                ! delta(bigOM+om2+om3)
                minc = lo_huge
                maxc = -lo_huge
                do ii = 1, 4
                    cvals(ii) = -omr2(ii, b2) - omr3(ii, b3)
                    minc = min(minc, cvals(ii))
                    maxc = max(maxc, cvals(ii))
                end do
                ii = max(floor(minc*invf), 1)
                jj = min(ceiling(maxc*invf), se%n_energy)
                do b1 = 1, dr%n_mode
                    sigma = dr%default_smearing(b1)
                    bvals = psisq1(:, b1, b2, b3)
                    do i = ii, jj
                        bigOM = se%energy_axis(i)
                        wts1 = lo_LV_tetrahedron_weights(cvals, bigOM - delta, lo_freqtol, sigma)
                        wts2 = lo_LV_tetrahedron_weights(cvals, bigOM + delta, lo_freqtol, sigma)
                        se%im_3ph(i, b1) = se%im_3ph(i, b1) - sum((wts1 + wts2)*bvals)*0.5_r8
                    end do
                end do

                ! delta(bigOM-om2+om3)
                minc = lo_huge
                maxc = -lo_huge
                do ii = 1, 4
                    cvals(ii) = omr2(ii, b2) - omr3(ii, b3)
                    minc = min(minc, cvals(ii))
                    maxc = max(maxc, cvals(ii))
                end do
                ii = max(floor(minc*invf), 1)
                jj = min(ceiling(maxc*invf), se%n_energy)
                do b1 = 1, dr%n_mode
                    sigma = dr%default_smearing(b1)
                    bvals = psisq2(:, b1, b2, b3)
                    do i = ii, jj
                        bigOM = se%energy_axis(i)
                        wts1 = lo_LV_tetrahedron_weights(cvals, bigOM - delta, lo_freqtol, sigma)
                        wts2 = lo_LV_tetrahedron_weights(cvals, bigOM + delta, lo_freqtol, sigma)
                        se%im_3ph(i, b1) = se%im_3ph(i, b1) + sum((wts1 + wts2)*bvals)*0.5_r8
                    end do
                end do

                ! delta(bigOM+om2-om3)
                minc = lo_huge
                maxc = -lo_huge
                do ii = 1, 4
                    cvals(ii) = -omr2(ii, b2) + omr3(ii, b3)
                    minc = min(minc, cvals(ii))
                    maxc = max(maxc, cvals(ii))
                end do
                ii = max(floor(minc*invf), 1)
                jj = min(ceiling(maxc*invf), se%n_energy)
                do b1 = 1, dr%n_mode
                    sigma = dr%default_smearing(b1)
                    bvals = psisq2(:, b1, b2, b3)
                    do i = ii, jj
                        bigOM = se%energy_axis(i)
                        wts1 = lo_LV_tetrahedron_weights(cvals, bigOM - delta, lo_freqtol, sigma)
                        wts2 = lo_LV_tetrahedron_weights(cvals, bigOM + delta, lo_freqtol, sigma)
                        se%im_3ph(i, b1) = se%im_3ph(i, b1) - sum((wts1 + wts2)*bvals)*0.5_r8
                    end do
                end do
            end do
        end do

        ! Dump a bit
        if (verbosity .gt. 0) then
            if (lo_trueNtimes(j, 127, qp%n_full_tet)) call lo_progressbar(' ... threephonon imaginary selfenergy', j, qp%n_full_tet)
        end if
    end do
    ! Add it up
    call mw%allreduce('sum', se%im_3ph)

    call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ! Fix degeneracies
    do b1 = 1, dr%n_mode
        buf = 0.0_r8
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            buf = buf + se%im_3ph(:, b2)
        end do
        buf = buf/real(wp%degeneracy(b1), r8)
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            se%im_3ph(:, b2) = buf
        end do
    end do
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (verbosity .gt. 0) call lo_progressbar(' ... threephonon imaginary selfenergy', qp%n_full_tet, qp%n_full_tet, walltime() - t0)

    call mem%deallocate(omr2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(psisq1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(psisq2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine
