submodule(lo_evaluate_phonon_self_energy) lo_evaluate_phonon_self_energy_threephonon
implicit none
contains

!> evaluate the isotope self-energy
module subroutine threephonon_imaginary_selfenergy(wp,se,qp,dr,sr,ise,temperature,mw,mem,verbosity)
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
    !> self-energy interpolation
    type(lo_interpolated_selfenergy_grid), intent(in) :: ise
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    select case (se%integrationtype)
    case (1:2)
        call threephonon_imaginary_selfenergy_gaussian(wp, se, sr, qp, dr, temperature, mw, mem, verbosity)
    case (3)
        call threephonon_imaginary_selfenergy_tetrahedron(wp, qp, sr, dr, temperature, se, mw, mem, verbosity)
    case (4)
        call threephonon_imaginary_selfenergy_convolution(se, wp, qp, dr, sr, ise, temperature, mw, mem, verbosity)
    case default
        call lo_stop_gracefully(['Unknown integration type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end select

end subroutine

!> threephonon self energy gaussian
subroutine threephonon_imaginary_selfenergy_gaussian(wp, se, sr, qp, dr, temperature, mw, mem, verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> scattering rates
    type(lo_listofscatteringrates), intent(in) :: sr
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:), allocatable :: buf
    real(r8) :: psisquare, psisq23, psisq32
    real(r8) :: sigma, om1, om2, om3, invf, s2, s3
    real(r8) :: plf1, plf2, t0, pref
    integer :: q, ctr, i, ii, jj, b1, b2, b3, ilo, ihi

    ! set the unit
    t0 = walltime()

    invf = se%n_energy/se%energy_axis(se%n_energy)

    ctr = 0
    se%im_3ph = 0.0_r8
    call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf = 0.0_r8

    if (sr%atgamma) then
        ! If we are doing the Gamma-point, we can use symmetry
        qloopirr: do q = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                ! Make it parallel
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                ! Prefactor for this q-point
                pref = threephonon_prefactor*qp%ip(q)%integration_weight
                ! Slightly smarter version that should be a little faster.
                do b2 = 1, dr%n_mode
                do b3 = b2, dr%n_mode
                    ! reset buffer
                    buf = 0.0_r8
                    ilo = lo_hugeint
                    ihi = -lo_hugeint
                    ! Get the smearing parameter, in case it's adaptive
                    select case (se%integrationtype)
                    case (1)
                        !sigma=(dr%default_smearing(b2)+dr%default_smearing(b3))*0.5_r8*se%smearing_prefactor
                        s2 = dr%default_smearing(b2)
                        s3 = dr%default_smearing(b3)
                        sigma = sqrt(s2**2 + s3**2)
                    case (2)
                        s2 = qp%adaptive_sigma(qp%ip(q)%radius, dr%iq(q)%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                        s3 = qp%adaptive_sigma(qp%ip(q)%radius, dr%iq(q)%vel(:, b3), dr%default_smearing(b3), se%smearing_prefactor)
                        sigma = sqrt(s2**2 + s3**2)
                        !v0=dr%iq(q)%vel(:,b2)
                        !v1=dr%iq(q)%vel(:,b3)
                        !sigma=qp%smearingparameter(v0-v1,(dr%default_smearing(b2)+dr%default_smearing(b3))*0.5_r8,se%smearing_prefactor)
                    end select

                    ! fetch frequencies
                    om1 = wp%omega(b1)
                    om2 = dr%iq(q)%omega(b2)
                    om3 = dr%iq(q)%omega(b3)
                    ! Get the S-function
                    plf1 = 1.0_r8 + lo_planck(temperature, om2) + lo_planck(temperature, om3)
                    plf2 = lo_planck(temperature, om3) - lo_planck(temperature, om2)

                    ! delta(bigOM-om2-om3)
                    ii = max(floor((om2 + om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((om2 + om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) + plf1*lo_gauss(se%energy_axis(i), om2 + om3, sigma)
                    end do
                    ! delta(bigOM+om2+om3)
                    ii = max(floor((-om2 - om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((-om2 - om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) - plf1*lo_gauss(se%energy_axis(i), -om2 - om3, sigma)
                    end do
                    ! delta(bigOM-om2+om3)
                    ii = max(floor((om2 - om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((om2 - om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) + plf2*lo_gauss(se%energy_axis(i), om2 - om3, sigma)
                    end do
                    ! delta(bigOM+om2-om3)
                    ii = max(floor((-om2 + om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((-om2 + om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) - plf2*lo_gauss(se%energy_axis(i), -om2 + om3, sigma)
                    end do

                    ! Increment the self-energy
                    if (ilo .lt. ihi) then
                        if (b2 .eq. b3) then
                            psisquare = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))*pref
                            se%im_3ph(ilo:ihi, b1) = se%im_3ph(ilo:ihi, b1) + buf(ilo:ihi)*psisquare
                        else
                            psisq23 = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))*pref
                            psisq32 = abs(sr%psi_3ph(b1, b3, b2, q)*conjg(sr%psi_3ph(b1, b3, b2, q)))*pref
                            se%im_3ph(ilo:ihi, b1) = se%im_3ph(ilo:ihi, b1) + buf(ilo:ihi)*(psisq23 + psisq32)
                        end if
                    end if
                end do
                end do

                if (verbosity .gt. 0) then
                    if (lo_trueNtimes(ctr, 127, qp%n_irr_point*dr%n_mode)) call lo_progressbar(' ... threephonon imaginary selfenergy', ctr, dr%n_mode*qp%n_irr_point)
                end if
            end do
        end do qloopirr
    else
        ! If not Gamma we have to do all of it.
        qloopfull: do q = 1, qp%n_full_point
            do b1 = 1, dr%n_mode
                ! Make it parallel
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                ! Prefactor for this q-point
                pref = threephonon_prefactor*qp%ap(q)%integration_weight
                ! Slightly smarter version that should be a little faster.
                do b2 = 1, dr%n_mode
                    !do b3=b2,dr%n_mode
                    do b3 = 1, dr%n_mode
                        ! reset buffer
                        buf = 0.0_r8
                        ilo = lo_hugeint
                        ihi = -lo_hugeint
                        ! Get the smearing parameter, in case it's adaptive
                        select case (se%integrationtype)
                        case (1)
                            ! s2=dr%default_smearing(b2)
                            ! s3=dr%default_smearing(b3)
                            ! sigma=sqrt(s2**2 + s3**2)
                            sigma = 1.0_r8*lo_frequency_THz_to_Hartree
                        case (2)
                            s2 = qp%adaptive_sigma(qp%ap(q)%radius, dr%aq(q)%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                            s3 = qp%adaptive_sigma(qp%ap(q)%radius, sr%vel3(:, b3, q), dr%default_smearing(b3), se%smearing_prefactor)
                            sigma = sqrt(s2**2 + s3**2)
                        end select
                        ! fetch frequencies
                        om1 = wp%omega(b1)
                        om2 = dr%aq(q)%omega(b2)
                        om3 = sr%omega3(b3, q)
                        ! Get the S-function
                        plf1 = 1.0_r8 + lo_planck(temperature, om2) + lo_planck(temperature, om3)
                        plf2 = lo_planck(temperature, om3) - lo_planck(temperature, om2)

                        ! delta(bigOM-om2-om3)
                        ii = max(floor((om2 + om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((om2 + om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) + plf1*lo_gauss(se%energy_axis(i), om2 + om3, sigma)
                        end do
                        ! delta(bigOM+om2+om3)
                        ii = max(floor((-om2 - om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((-om2 - om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) - plf1*lo_gauss(se%energy_axis(i), -om2 - om3, sigma)
                        end do
                        ! delta(bigOM-om2+om3)
                        ii = max(floor((om2 - om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((om2 - om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) + plf2*lo_gauss(se%energy_axis(i), om2 - om3, sigma)
                        end do
                        ! delta(bigOM+om2-om3)
                        ii = max(floor((-om2 + om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((-om2 + om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) - plf2*lo_gauss(se%energy_axis(i), -om2 + om3, sigma)
                        end do

                        ! Increment the self-energy
                        if (ilo .lt. ihi) then
                            !if ( b2 .eq. b3 ) then
                            psisquare = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))*pref
                            se%im_3ph(ilo:ihi, b1) = se%im_3ph(ilo:ihi, b1) + buf(ilo:ihi)*psisquare
                            !else
                            !    psisq23=abs( sr%psi_3ph(b1,b2,b3,q)*conjg(sr%psi_3ph(b1,b2,b3,q)) )*pref
                            !    psisq32=abs( sr%psi_3ph(b1,b3,b2,q)*conjg(sr%psi_3ph(b1,b3,b2,q)) )*pref
                            !    se%im_3ph(ilo:ihi,b1)=se%im_3ph(ilo:ihi,b1)+buf(ilo:ihi)*(psisq23+psisq32)
                            !endif
                        end if
                    end do
                end do

                if (verbosity .gt. 0) then
                    if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) call lo_progressbar(' ... threephonon imaginary selfenergy', ctr, dr%n_mode*qp%n_full_point)
                end if
            end do
        end do qloopfull
    end if

    ! Sum it up
    call mw%allreduce('sum', se%im_3ph)

    ! Fix degeneracies?
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
    if (verbosity .gt. 0) call lo_progressbar(' ... threephonon imaginary selfenergy', dr%n_mode*qp%n_full_point, dr%n_mode*qp%n_full_point, walltime() - t0)
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
    type(lo_listofscatteringrates), intent(in) :: sr
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

!> Get the three-phonon self-energy via convolutions of the spectral function.
subroutine threephonon_imaginary_selfenergy_convolution(se, wp, qp, dr, sr, ise, temperature, mw, mem, verbosity)
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> scattering rates
    type(lo_listofscatteringrates), intent(in) :: sr
    !> tabulated spectral functions
    type(lo_interpolated_selfenergy_grid), intent(in) :: ise
    !> structure
    !type(lo_crystalstructure), intent(in) :: p
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity


    real(r8) :: timer, t0, t1

    timer = walltime()
    t0 = timer
    t1 = timer

    ! Might add some more pre-processing here.
    init: block
        if (verbosity .gt. 0) then
            write (lo_iou, *) '... convolution based integration'
        end if
        ! Reset values
        se%im_iso = 0.0_r8
        se%im_3ph = 0.0_r8
    end block init


    ! Try this again, this time with more feeling. The next attempt will
    ! make do with much less Fourier transforms, this one is to get it right.
    newattempt: block
        type(lo_convolution_helper) :: ch
        real(r8), dimension(:), allocatable :: buf_sfun, buf_sigma2, buf_sigma3
        real(r8), dimension(3) :: v0
        real(r8) :: pref, sigma2, sigma3, sigma, psisq
        integer :: imode, iq, jq, kq, mode1, mode2, mode3

        ! Set up the helper guy for the convolutions
        call ch%generate(se%energy_axis, temperature, dr%n_mode)

        ! Some temporary buffers
        allocate (buf_sigma2(dr%n_mode))
        allocate (buf_sigma3(dr%n_mode))
        buf_sigma2 = 0.0_r8
        buf_sigma3 = 0.0_r8
        allocate (buf_sfun(se%n_energy))
        buf_sfun = 0.0_r8

        if (verbosity .gt. 0) call lo_progressbar_init()
        qploopfull: do iq = 1, qp%n_full_point
            ! Make it MPI parallel
            if (mod(iq, mw%n) .ne. mw%r) cycle

            ! Pre-fetch and massage spectral functions before the mode loop
            fixspectralfunctions: block

                ! So, I will need the spectral functions for mode 2 and mode 3. Think it
                ! makese sense to evaluate them on-the-fly, as to not dramatically increase
                ! memory usage.
                !call ise%evaluate(qp,qp%ap(iq)%r,dr%aq(iq),sigma_Re,sigma_Im)
                !call ise%evaluate(qp,sr%qvec3(:,iq),sr%wp3(iq),sigma_Re,sigma_Im)

            end block fixspectralfunctions

!                     ! This is the closed-grid or Gamma-point case, where everything has been prepared beforehand.
!                     ! First sort out the q-indices. This is q'
!                     jq = qp%ap(iq)%irreducible_index
!                     ! Sort out q''
!                     v0 = matmul(p%inv_reciprocal_latticevectors, sr%qvec3(:, iq))
!                     kq = index_on_grid(qp, v0)
!                     kq = qp%ap(kq)%irreducible_index

!                     ! Pre-fetch the spectral functions
!                     call ch%buffer_spectral_functions(isf%spectralfunction(:, :, jq), isf%spectralfunction(:, :, kq))

!                     ! ! Pre-fetch the spectral functions
!                     ! do imode=1,dr%n_mode
!                     !     do ie=1,n
!                     !         buf_sf2_planck(ie,imode)=isf%spectralfunction(ie+1,imode,jq)
!                     !         buf_sf2_planck(-ie,imode)=-buf_sf2_planck(ie,imode)
!                     !         buf_sf3_planck(ie,imode)=isf%spectralfunction(ie+1,imode,kq)
!                     !         buf_sf3_planck(-ie,imode)=-buf_sf3_planck(ie,imode)
!                     !     enddo
!                     !     buf_sf2_planck(:,imode)=buf_sf2_planck(:,imode)*buf_planck_plus_one
!                     !     buf_sf3_planck(:,imode)=buf_sf3_planck(:,imode)*buf_planck_plus_one
!                     !     ! Insert FFT here
!                     ! enddo

!                     ! ! Set them to sharp values instead, for debugging
!                     ! buf_sf2_planck=0.0_r8
!                     ! do imode=1,dr%n_mode
!                     !     if ( dr%iq(iq)%omega(imode) .lt. lo_freqtol ) cycle
!                     !     pref=lo_huge
!                     !     j=-1
!                     !     do ie=1,n
!                     !         if ( abs(buf_x(ie)-dr%iq(iq)%omega(imode)) .lt. pref ) then
!                     !             pref=abs(buf_x(ie)-dr%iq(iq)%omega(imode))
!                     !             j=ie
!                     !         endif
!                     !     enddo
!                     !     buf_sf2_planck(j,imode)=1.0_r8/deltax
!                     !     buf_sf2_planck(-j,imode)=-1.0_r8/deltax
!                     !     buf_sf2_planck(:,imode)=buf_sf2_planck(:,imode)*buf_planck_plus_one
!                     ! enddo
!                 case default
!                     call lo_stop_gracefully(['Bad integration type'], lo_exitcode_param, __FILE__, __LINE__)
!                 end select

!                 ! Prefactor for this q-point
!                 pref = threephonon_prefactor*qp%ap(iq)%integration_weight

!                 ! Excellent. Now start convoluting every which way.
!                 mode2loop2: do mode2 = 1, dr%n_mode
!                 mode3loop2: do mode3 = 1, dr%n_mode
!                     ! Skip acoustic modes?
!                     if (dr%iq(jq)%omega(mode2) .lt. lo_freqtol) cycle
!                     if (dr%iq(kq)%omega(mode3) .lt. lo_freqtol) cycle
!                     ! Smear?
!                     sigma2 = qp%adaptive_sigma(qp%ip(jq)%radius, dr%iq(jq)%vel(:, mode2), dr%default_smearing(mode2), se%smearing_prefactor)
!                     sigma3 = qp%adaptive_sigma(qp%ip(kq)%radius, dr%iq(kq)%vel(:, mode3), dr%default_smearing(mode3), se%smearing_prefactor)
!                     sigma = sqrt(sigma2**2 + sigma3**2)

!                     ! ! Do the convolution over modes
!                     ! call lo_convolution(buf_sf2_planck(:,mode2),buf_sf2_planck(:,mode3),deltax,.false.,buf_conv)
!                     ! ! Adjust with planck factors
!                     ! buf_conv=buf_conv*buf_inv_planck_plus_one
!                     ! buf_sfun=0.0_r8
!                     ! do ie=1,n
!                     !     buf_sfun(ie+1)=( buf_conv(ie)-buf_conv(-ie) )*0.5_r8
!                     ! enddo
!                     ! call gaussian_smear_spectral_function(se%energy_axis,sigma,buf_sfun)

!                     call ch%convolute_to_sfun(mode2, mode3, sigma, buf_sfun)

!                     ! Excellent, everything was convoluted properly. Start accumulating self-energy
!                     mode1loop2: do mode1 = 1, dr%n_mode
!                         psisq = abs(sr%psi_3ph(mode1, mode2, mode3, iq)*conjg(sr%psi_3ph(mode1, mode2, mode3, iq)))
!                         se%im_3ph(:, mode1) = se%im_3ph(:, mode1) + buf_sfun*pref*psisq
!                     end do mode1loop2
!                 end do mode3loop2
!                 end do mode2loop2

!                 ! Since I have the spectral functions available I might as well sort
!                 ! out the isotope scattering here.
!                 pref = isotope_prefactor*qp%ap(iq)%integration_weight
!                 do mode2 = 1, dr%n_mode
!                     buf_sfun = isf%spectralfunction(:, mode2, jq)
!                     sigma2 = qp%adaptive_sigma(qp%ip(jq)%radius, dr%iq(jq)%vel(:, mode2), dr%default_smearing(mode2), se%smearing_prefactor)
!                     call gaussian_smear_spectral_function(se%energy_axis, sigma2, buf_sfun)
!                     do mode1 = 1, dr%n_mode
!                         se%im_iso(:, mode1) = se%im_iso(:, mode1) + buf_sfun*sr%psi_iso(mode1, mode2, iq)*pref
!                     end do
!                 end do

!                 if (verbosity .gt. 0 .and. iq .lt. qp%n_full_point) then
!                     t1 = walltime()
!                     call lo_progressbar(' ... imaginary selfenergy', iq, qp%n_full_point, t1 - t0)
!                 end if
            end do qploopfull
        ! Cleanup
        call ch%destroy()

        ! Add together
        call mw%allreduce('sum', se%im_3ph)
        call mw%allreduce('sum', se%im_iso)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... imaginary selfenergy', qp%n_full_point, qp%n_full_point, t1 - t0)
            t0 = t1
        end if

        ! if ( mw%talk ) then
        !         call h5%init(__FILE__,__LINE__)
        !         call h5%open_file('write','sbuf.hdf5')
        !
        !         call h5%store_data(allsbuf,h5%file_id,'sfun')
        !
        !         call h5%close_file()
        !         call h5%destroy(__FILE__,__LINE__)
        ! endif
    end block newattempt

!     ! actualcalculation: block
!     !     real(r8), dimension(:,:), allocatable :: bufRe,bufIm,buf_spectral2,buf_spectral3
!     !     real(r8), dimension(:), allocatable :: sfun,taperfn
!     !     real(r8), dimension(:), allocatable :: buf_fpre,buf_om,buf_conv
!     !     real(r8), dimension(:), allocatable :: buf_odd_sf2a,buf_odd_sf2b
!     !     real(r8), dimension(:), allocatable :: buf_odd_sf3a,buf_odd_sf3b
!     !     real(r8), dimension(3) :: v0
!     !     real(r8) :: psisq,pref,minIm,sigma
!     !     integer :: mode1,mode2,mode3,iq,jq,ctr
!     !     integer :: ie,n
!     !
!     !     ! Oddly allocated buffers:
!     !     n=se%n_energy
!     !     allocate(buf_fpre(-n+1:n))
!     !     allocate(buf_om(-n+1:n))
!     !     allocate(buf_odd_sf2a(-n+1:n))
!     !     allocate(buf_odd_sf2b(-n+1:n))
!     !     allocate(buf_odd_sf3a(-n+1:n))
!     !     allocate(buf_odd_sf3b(-n+1:n))
!     !     allocate(buf_conv(-n+1:n))
!     !     buf_fpre=0.0_r8
!     !     buf_om=0.0_r8
!     !     buf_odd_sf2a=0.0_r8
!     !     buf_odd_sf2b=0.0_r8
!     !     buf_odd_sf3a=0.0_r8
!     !     buf_odd_sf3b=0.0_r8
!     !     buf_conv=0.0_r8
!     !
!     !     ! get the planck prefactor
!     !     buf_fpre=0.0_r8
!     !     do ie=2,se%n_energy
!     !         buf_fpre(ie)=lo_planck(temperature,se%energy_axis(ie))+0.5_r8
!     !         buf_fpre(1-ie)=-buf_fpre(ie)
!     !     enddo
!     !
!     !     ! Normally allocated buffers:
!     !     call mem%allocate(sfun,se%n_energy,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%allocate(taperfn,se%n_energy,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%allocate(bufRe,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%allocate(bufIm,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%allocate(buf_spectral2,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%allocate(buf_spectral3,[se%n_energy,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     sfun=0.0_r8
!     !     taperfn=0.0_r8
!     !     bufRe=0.0_r8
!     !     bufIm=0.0_r8
!     !     buf_spectral2=0.0_r8
!     !     buf_spectral3=0.0_r8
!     !     minIm=( se%energy_axis(2)-se%energy_axis(1) )
!     !     call taperfn_im(se%energy_axis,dr%omega_max,dr%omega_min,taperfn)
!     !     ! Stop the init timer
!     !
!     !     if ( sr%atgamma ) then
!     !         ! Use symmetry
!     !         if ( verbosity .gt. 0 ) call lo_progressbar_init()
!     !         qploopirr: do iq=1,qp%n_irr_point
!     !             ! Make it MPI parallel
!     !             if ( mod(iq,mw%n) .ne. mw%r ) cycle
!     !
!     !             ! First thing, I need the spectral functions at q'
!     !
!     !             ! We need the spectral functions at q', which is also the one at q''. Two different ways to go about it.
!     !             select case(se%integrationtype)
!     !             case(4)
!     !                 ! Here I calculate the spectral functions from the interpolated self-energy, and then smear them appropriately.
!     !
!     !                 ! First thing, I need the spectral functions at q'
!     !                 call ise%diagonal_selfenergy(p,qp%ip(iq)%r,dr%iq(iq)%omega,dr%iq(iq)%egv,bufRe,bufIm)
!     !                 bufIm=max(bufIm,minIm)
!     !                 do mode2=1,se%n_mode
!     !                     if ( dr%aq(iq)%omega(mode2) .gt. lo_freqtol ) then
!     !                         bufIm(:,mode2)=bufIm(:,mode2)*taperfn
!     !                         call evaluate_spectral_function(ise%energy,bufIm(:,mode2),bufRe(:,mode2),dr%iq(iq)%omega(mode2),buf_spectral2(:,mode2))
!     !                         sigma=qp%adaptive_sigma(qp%ip(iq)%radius,dr%iq(iq)%vel(:,mode2),dr%default_smearing(mode2),se%smearing_prefactor)
!     !                         call gaussian_smear_spectral_function(ise%energy,sigma,buf_spectral2(:,mode2))
!     !                         buf_spectral2(:,mode2)=buf_spectral2(:,mode2)/lo_trapezoid_integration(ise%energy,buf_spectral2(:,mode2))
!     !                     else
!     !                         buf_spectral2(:,mode2)=0.0_r8
!     !                     endif
!     !                 enddo
!     !             case(5)
!     !                 ! This is the closed-grid or Gamma-point case, where everything has been prepared beforehand.
!     !                 ! Sort out q' first
!     !                 buf_spectral2=isf%spectralfunction(:,:,iq)
!     !             end select
!     !
!     !             ! Prefactor for this q-point
!     !             pref=threephonon_prefactor*qp%ip(iq)%integration_weight
!     !
!     !             ! Excellent. Now start convoluting every which way.
!     !             mode2loop1: do mode2=1,dr%n_mode
!     !             mode3loop1: do mode3=mode2,dr%n_mode
!     !                 ! Prepare for convolution
!     !                 buf_odd_sf2a=0.0_r8
!     !                 buf_odd_sf2b=0.0_r8
!     !                 buf_odd_sf3a=0.0_r8
!     !                 buf_odd_sf3b=0.0_r8
!     !                 do ie=2,se%n_energy
!     !                     buf_odd_sf2a(ie)=buf_spectral2(ie,mode2)
!     !                     buf_odd_sf3a(ie)=buf_spectral2(ie,mode3)
!     !                     buf_odd_sf2a(1-ie)=-buf_odd_sf2a(ie)
!     !                     buf_odd_sf3a(1-ie)=-buf_odd_sf3a(ie)
!     !                 enddo
!     !                 buf_odd_sf2b=buf_odd_sf2a*buf_fpre
!     !                 buf_odd_sf3b=buf_odd_sf3a*buf_fpre
!     !
!     !                 call lo_abcd_convolution(buf_odd_sf2a,buf_odd_sf3b,buf_odd_sf3a,buf_odd_sf2b,se%energy_axis(2)-se%energy_axis(1),buf_conv)
!     !                 ! Post-convolution the values are FFT-shifted oddly so collect to the right place.
!     !                 sfun=0.0_r8
!     !                 do ie=1,se%n_energy
!     !                     sfun(ie)=buf_conv(se%n_energy+1-ie)
!     !                 enddo
!     !
!     !                 ! Excellent, everything was convoluted properly. Start accumulating self-energy
!     !                 mode1loop1: do mode1=1,dr%n_mode
!     !                     if ( mode2 .eq. mode3 ) then
!     !                         psisq=abs( sr%psi_3ph(mode1,mode2,mode3,iq) * conjg( sr%psi_3ph(mode1,mode2,mode3,iq) ) )
!     !                         se%im_3ph(:,mode1)=se%im_3ph(:,mode1)-sfun*pref*psisq
!     !                     else
!     !                         psisq=abs( sr%psi_3ph(mode1,mode2,mode3,iq) * conjg( sr%psi_3ph(mode1,mode2,mode3,iq) ) )+&
!     !                               abs( sr%psi_3ph(mode1,mode3,mode2,iq) * conjg( sr%psi_3ph(mode1,mode3,mode2,iq) ) )
!     !                         se%im_3ph(:,mode1)=se%im_3ph(:,mode1)-sfun*pref*psisq
!     !                     endif
!     !                 enddo mode1loop1
!     !
!     !             enddo mode3loop1
!     !             enddo mode2loop1
!     !
!     !             ! Since I have the spectral functions available I might as well sort
!     !             ! out the isotope scattering here.
!     !             pref=isotope_prefactor*qp%ip(iq)%integration_weight
!     !             do mode2=1,dr%n_mode
!     !                 sfun=buf_spectral2(:,mode2)*pref
!     !                 do mode1=1,dr%n_mode
!     !                     se%im_iso(:,mode1)=se%im_iso(:,mode1)+sfun*sr%psi_iso(mode1,mode2,iq)
!     !                 enddo
!     !             enddo
!     !
!     !             if ( verbosity .gt. 0 .and. iq .lt. qp%n_irr_point ) then
!     !                 t1=walltime()
!     !                 call lo_progressbar(' ... imaginary selfenergy',iq,qp%n_irr_point,t1-t0)
!     !             endif
!     !         enddo qploopirr
!     !     else
!     !         ! Do not use symmetry
!     !         if ( verbosity .gt. 0 ) call lo_progressbar_init()
!     !         ctr=0
!     !         qploopfull: do iq=1,qp%n_full_point
!     !             ! Make it MPI parallel
!     !             if ( mod(iq,mw%n) .ne. mw%r ) cycle
!     !
!     !             ! We need the spectral functions at q' and q''. Two different ways to go about it.
!     !             select case(se%integrationtype)
!     !             case(4)
!     !                 ! Here I calculate the spectral functions from the interpolated self-energy, and then smear them appropriately.
!     !
!     !                 ! First thing, I need the spectral functions at q'
!     !                 call ise%diagonal_selfenergy(p,qp%ap(iq)%r,dr%aq(iq)%omega,dr%aq(iq)%egv,bufRe,bufIm)
!     !                 bufIm=max(bufIm,minIm)
!     !                 do mode2=1,se%n_mode
!     !                     if ( dr%aq(iq)%omega(mode2) .gt. lo_freqtol ) then
!     !                         bufIm(:,mode2)=bufIm(:,mode2)*taperfn
!     !                         call evaluate_spectral_function(ise%energy,bufIm(:,mode2),bufRe(:,mode2),dr%aq(iq)%omega(mode2),buf_spectral2(:,mode2))
!     !                         sigma=qp%adaptive_sigma(qp%ap(iq)%radius,dr%aq(iq)%vel(:,mode2),dr%default_smearing(mode2),se%smearing_prefactor)
!     !                         call gaussian_smear_spectral_function(ise%energy,sigma,buf_spectral2(:,mode2))
!     !                         buf_spectral2(:,mode2)=buf_spectral2(:,mode2)/lo_trapezoid_integration(ise%energy,buf_spectral2(:,mode2))
!     !                     else
!     !                         buf_spectral2(:,mode2)=0.0_r8
!     !                     endif
!     !                 enddo
!     !
!     !                 ! Then the spectral functions at q''
!     !                 call ise%diagonal_selfenergy(p,sr%qvec3(:,iq),sr%omega3(:,iq),sr%egv3(:,:,iq),bufRe,bufIm)
!     !                 bufIm=max(bufIm,minIm)
!     !                 do mode3=1,se%n_mode
!     !                     if ( sr%omega3(mode3,iq) .gt. lo_freqtol ) then
!     !                         bufIm(:,mode3)=bufIm(:,mode3)*taperfn
!     !                         call evaluate_spectral_function(ise%energy,bufIm(:,mode3),bufRe(:,mode3),sr%omega3(mode3,iq),buf_spectral3(:,mode3))
!     !                         sigma=qp%adaptive_sigma(qp%ap(iq)%radius,sr%vel3(:,mode3,iq),dr%default_smearing(mode3),se%smearing_prefactor)
!     !                         call gaussian_smear_spectral_function(ise%energy,sigma,buf_spectral3(:,mode3))
!     !                         buf_spectral3(:,mode3)=buf_spectral3(:,mode3)/lo_trapezoid_integration(ise%energy,buf_spectral3(:,mode3))
!     !                     else
!     !                         buf_spectral3(:,mode3)=0.0_r8
!     !                     endif
!     !                 enddo
!     !             case(5)
!     !                 ! This is the closed-grid case, where everything has been prepared beforehand.
!     !                 ! Sort out q' first
!     !                 jq=qp%ap(iq)%irreducible_index
!     !                 buf_spectral2=isf%spectralfunction(:,:,jq)
!     !
!     !                 ! Sort out q''
!     !                 v0=matmul(p%inv_reciprocal_latticevectors,sr%qvec3(:,iq))
!     !                 jq=index_on_grid(qp,v0)
!     !                 jq=qp%ap(jq)%irreducible_index
!     !                 buf_spectral3=isf%spectralfunction(:,:,jq)
!     !             end select
!     !
!     !             ! Prefactor for this q-point
!     !             pref=threephonon_prefactor*qp%ap(iq)%integration_weight
!     !             ! Excellent. Now start convoluting every which way.
!     !             mode2loop2: do mode2=1,dr%n_mode
!     !             mode3loop2: do mode3=1,dr%n_mode
!     !                 ! ! Prepare for convolution
!     !                 buf_odd_sf2a=0.0_r8
!     !                 buf_odd_sf2b=0.0_r8
!     !                 buf_odd_sf3a=0.0_r8
!     !                 buf_odd_sf3b=0.0_r8
!     !                 do ie=2,se%n_energy
!     !                     buf_odd_sf2a(ie)=buf_spectral2(ie,mode2)
!     !                     buf_odd_sf3a(ie)=buf_spectral3(ie,mode3)
!     !                     buf_odd_sf2a(1-ie)=-buf_odd_sf2a(ie)
!     !                     buf_odd_sf3a(1-ie)=-buf_odd_sf3a(ie)
!     !                 enddo
!     !                 buf_odd_sf2b=buf_odd_sf2a*buf_fpre
!     !                 buf_odd_sf3b=buf_odd_sf3a*buf_fpre
!     !
!     !                 call lo_abcd_convolution(buf_odd_sf2a,buf_odd_sf3b,buf_odd_sf3a,buf_odd_sf2b,se%energy_axis(2)-se%energy_axis(1),buf_conv)
!     !                 ! Post-convolution the values are FFT-shifted oddly so collect to the right place.
!     !                 sfun=0.0_r8
!     !                 do ie=1,se%n_energy
!     !                     sfun(ie)=buf_conv(se%n_energy+1-ie)
!     !                 enddo
!     !
!     !                 ! Excellent, everything was convoluted properly. Start accumulating self-energy
!     !                 mode1loop2: do mode1=1,dr%n_mode
!     !                     psisq=abs( sr%psi_3ph(mode1,mode2,mode3,iq) * conjg( sr%psi_3ph(mode1,mode2,mode3,iq) ) )
!     !                     se%im_3ph(:,mode1)=se%im_3ph(:,mode1)-sfun*pref*psisq
!     !                 enddo mode1loop2
!     !
!     !             enddo mode3loop2
!     !             enddo mode2loop2
!     !
!     !             ! Since I have the spectral functions available I might as well sort
!     !             ! out the isotope scattering here.
!     !             pref=isotope_prefactor*qp%ap(iq)%integration_weight
!     !             do mode2=1,dr%n_mode
!     !                 sfun=buf_spectral2(:,mode2)*pref
!     !                 do mode1=1,dr%n_mode
!     !                     se%im_iso(:,mode1)=se%im_iso(:,mode1)+sfun*sr%psi_iso(mode1,mode2,iq)
!     !                 enddo
!     !             enddo
!     !
!     !             if ( verbosity .gt. 0 .and. iq .lt. qp%n_full_point ) then
!     !                 t1=walltime()
!     !                 call lo_progressbar(' ... imaginary selfenergy',iq,qp%n_full_point,t1-t0)
!     !             endif
!     !         enddo qploopfull
!     !     endif ! use symmetry
!     !
!     !     call mem%deallocate(sfun,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%deallocate(taperfn,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%deallocate(bufRe,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%deallocate(bufIm,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%deallocate(buf_spectral2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !     call mem%deallocate(buf_spectral3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     !
!     !     ! Add together
!     !     call mw%allreduce('sum',se%im_3ph)
!     !     call mw%allreduce('sum',se%im_iso)
!     !
!     !     if ( verbosity .gt. 0 ) then
!     !         t1=walltime()
!     !         call lo_progressbar(' ... imaginary selfenergy',qp%n_full_point,qp%n_full_point,t1-t0)
!     !         t0=t1
!     !     endif
!     ! end block actualcalculation

!     ! Make sure the whole thing makes sense in the end.
!     massage: block
!         real(r8), dimension(:), allocatable :: taper, buf1, buf2
!         integer :: mode1, mode2, i

!         call mem%allocate(taper, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf1, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%allocate(buf2, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         taper = 0.0_r8
!         buf1 = 0.0_r8
!         buf2 = 0.0_r8

!         ! Sort out degeneracies
!         do mode1 = 1, dr%n_mode
!             buf1 = 0.0_r8
!             buf2 = 0.0_r8
!             do i = 1, wp%degeneracy(mode1)
!                 mode2 = wp%degenmode(i, mode1)
!                 buf1 = buf1 + se%im_iso(:, mode2)
!                 buf2 = buf2 + se%im_3ph(:, mode2)
!             end do
!             buf1 = buf1/real(wp%degeneracy(mode1), r8)
!             buf2 = buf2/real(wp%degeneracy(mode1), r8)
!             do i = 1, wp%degeneracy(mode1)
!                 mode2 = wp%degenmode(i, mode1)
!                 se%im_iso(:, mode2) = buf1
!                 se%im_3ph(:, mode2) = buf2
!             end do
!         end do
!         ! Get the tapering function
!         call taperfn_im(se%energy_axis, dr%omega_max, dr%omega_min, taper)
!         ! Taper the self-energies so that they are zero where they should be.
!         do mode1 = 1, dr%n_mode
!             se%im_3ph(:, mode1) = max(se%im_3ph(:, mode1), 0.0_r8)
!             se%im_iso(:, mode1) = max(se%im_iso(:, mode1), 0.0_r8)
!             se%im_3ph(:, mode1) = se%im_3ph(:, mode1)*taper
!             se%im_iso(:, mode1) = se%im_iso(:, mode1)*taper
!         end do

!         call mem%deallocate(taper, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!         call mem%deallocate(buf2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!     end block massage

    if ( mw%talk) write(*,*) 'done here for now ',__FILE__,__LINE__
    call mw%destroy()
    stop


end subroutine

end submodule