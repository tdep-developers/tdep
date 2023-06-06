
!> count three-phonons scattering events with gaussian smearing
subroutine threephonon_gaussian_oneqp(qp, dr, scq, gi1, thres, smearing_prefactor, integrationtype, mem)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> scattering phase space
    type(lo_3phqp2), intent(inout) :: scq
    !> index to q1
    integer, intent(in) :: gi1
    !> gaussian threshold
    real(r8), intent(in) :: thres
    !> smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    !> how to integrate
    integer, intent(in) :: integrationtype
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:, :), allocatable :: vel1, vel2, vel3
    real(r8), dimension(:), allocatable :: omr1, omr2, omr3
    real(r8) :: deltafunction, om1, om2, om3, sigma, omthres
    real(r8) :: sig1, sig2, sig3
    integer, dimension(:, :, :), allocatable :: bc1, bc2
    integer, dimension(3) :: dims
    integer :: gi2, gi3, b1, b2, b3, i

    ! Threshold for small frequencies
    omthres = dr%omega_min*0.2_r8

    ! Some temporary space
    call mem%allocate(bc1, [dr%n_mode, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(bc2, [dr%n_mode, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(vel1, [3, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(vel2, [3, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(vel3, [3, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    bc1 = 0
    bc2 = 0
    omr1 = 0.0_r8
    omr2 = 0.0_r8
    omr3 = 0.0_r8
    vel2 = 0.0_r8
    vel3 = 0.0_r8

    ! grid dimensions
    select type (qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    ! Set the q-index and omega for q, and reset counters
    scq%gi1 = gi1
    bc1 = 0
    bc2 = 0
    omr1 = dr%aq(gi1)%omega
    vel1 = dr%aq(gi1)%vel

    ! Do the actual counting
    do i = 1, qp%n_full_point
        ! Get q'', This is q1+q2+q3=G
        gi2 = i
        gi3 = fft_third_grid_index(gi1, gi2, dims)

        omr2 = dr%aq(gi2)%omega
        omr3 = dr%aq(gi3)%omega
        vel2 = dr%aq(gi2)%vel
        vel3 = dr%aq(gi3)%vel

        ! Count events
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
        do b3 = 1, dr%n_mode
            om1 = omr1(b1)
            om2 = omr2(b2)
            om3 = omr3(b3)
            ! plus-events, first get sigma
            select case (integrationtype)
            case (1)
                ! sig1=dr%default_smearing(b1)
                ! sig2=dr%default_smearing(b2)
                ! sig3=dr%default_smearing(b3)
                ! sigma=sqrt(sig1**2 + sig2**2 + sig3**2)
                sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_prefactor
            case (2)
                !sigma=qp%smearingparameter(vel2(:,b2)-vel3(:,b3),(dr%default_smearing(b2)+dr%default_smearing(b3))*0.5_r8,smearing_prefactor)
                sig1 = qp%adaptive_sigma(qp%ap(gi1)%radius, vel1(:, b1), dr%default_smearing(b1), smearing_prefactor)
                sig2 = qp%adaptive_sigma(qp%ap(gi2)%radius, vel2(:, b2), dr%default_smearing(b2), smearing_prefactor)
                sig3 = qp%adaptive_sigma(qp%ap(gi3)%radius, vel3(:, b3), dr%default_smearing(b3), smearing_prefactor)
                sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                !sigma=sqrt(sig2**2 + sig3**2)
            end select

            if (abs(om1 - om3 + om2) .lt. thres*sigma) then
            if (om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres) then
                bc1(b1, b2, b3) = bc1(b1, b2, b3) + 1
            end if
            end if

            if (abs(om1 - om2 - om3) .lt. thres*sigma) then
            if (om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres) then
                bc2(b1, b2, b3) = bc2(b1, b2, b3) + 1
            end if
            end if
        end do
        end do
        end do
    end do

    ! Allocate the storage
    do b1 = 1, dr%n_mode
    do b2 = 1, dr%n_mode
    do b3 = 1, dr%n_mode
        scq%plus(b1, b2, b3)%n = bc1(b1, b2, b3)
        scq%minus(b1, b2, b3)%n = bc2(b1, b2, b3)
        if (scq%plus(b1, b2, b3)%n .gt. 0) then
            allocate (scq%plus(b1, b2, b3)%e(bc1(b1, b2, b3)))
        end if
        if (scq%minus(b1, b2, b3)%n .gt. 0) then
            allocate (scq%minus(b1, b2, b3)%e(bc2(b1, b2, b3)))
        end if
    end do
    end do
    end do

    ! Count again and store things
    bc1 = 0
    bc2 = 0
    do i = 1, qp%n_full_point
        ! This is q1+q2+q3=G
        gi2 = i
        gi3 = fft_third_grid_index(gi1, gi2, dims)
        omr2 = dr%aq(gi2)%omega
        omr3 = dr%aq(gi3)%omega
        vel2 = dr%aq(gi2)%vel
        vel3 = dr%aq(gi3)%vel
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
        do b3 = 1, dr%n_mode
            om1 = omr1(b1)
            om2 = omr2(b2)
            om3 = omr3(b3)
            ! First calculate the adaptive sigma in a consistent way
            ! plus-events, first get sigma
            select case (integrationtype)
            case (1)
                ! sig1=dr%default_smearing(b1)
                ! sig2=dr%default_smearing(b2)
                ! sig3=dr%default_smearing(b3)
                ! sigma=sqrt(sig1**2 + sig2**2 + sig3**2)
                sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_prefactor
            case (2)
                sig1 = qp%adaptive_sigma(qp%ap(gi1)%radius, vel1(:, b1), dr%default_smearing(b1), smearing_prefactor)
                sig2 = qp%adaptive_sigma(qp%ap(gi2)%radius, vel2(:, b2), dr%default_smearing(b2), smearing_prefactor)
                sig3 = qp%adaptive_sigma(qp%ap(gi3)%radius, vel3(:, b3), dr%default_smearing(b3), smearing_prefactor)
                sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                !sigma=sqrt(sig2**2 + sig3**2)
            end select
            if (abs(om1 - om3 + om2) .lt. thres*sigma) then
            if (om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres) then
                deltafunction = lo_gauss(om1, om3 - om2, sigma)
                bc1(b1, b2, b3) = bc1(b1, b2, b3) + 1
                scq%plus(b1, b2, b3)%e(bc1(b1, b2, b3))%gi2 = gi2
                scq%plus(b1, b2, b3)%e(bc1(b1, b2, b3))%gi3 = gi3
                scq%plus(b1, b2, b3)%e(bc1(b1, b2, b3))%deltafunction = deltafunction*qp%ap(i)%integration_weight
            end if
            end if

            if (abs(om1 - om2 - om3) .lt. thres*sigma) then
            if (om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres) then
                bc2(b1, b2, b3) = bc2(b1, b2, b3) + 1
                deltafunction = lo_gauss(om1, om2 + om3, sigma)
                scq%minus(b1, b2, b3)%e(bc2(b1, b2, b3))%gi2 = gi2
                scq%minus(b1, b2, b3)%e(bc2(b1, b2, b3))%gi3 = gi3
                scq%minus(b1, b2, b3)%e(bc2(b1, b2, b3))%deltafunction = deltafunction*qp%ap(i)%integration_weight
            end if
            end if
        end do
        end do
        end do
    end do
    ! Cleanup
    call mem%deallocate(bc1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(bc2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(vel1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(vel2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(vel3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Gaussian integration weights for isotope scattering, from one q-point
subroutine iso_gaussian_oneqp(qp, dr, scq, gi1, thres, smearing_adjustment, integrationtype, mem)
    !> q-point mesh
    !type(lo_fft_mesh), intent(in) :: qp
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> qpoint in question
    type(lo_iso2), intent(inout) :: scq
    !> gridindex in question
    integer, intent(in) :: gi1
    !> threshold to cut of gaussian
    real(r8), intent(in) :: thres
    !> smearing adjustment
    real(r8), intent(in) :: smearing_adjustment
    !> how to integrate
    integer, intent(in) :: integrationtype
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer :: i, ii, b1, b2, gi2
    integer, dimension(:, :), allocatable :: bandcounter
    real(r8) :: deltafunction, om1, om2, omthres, sigma

    ! Values for q
    omthres = dr%omega_min*0.2_r8
    scq%gi1 = gi1

    ! count first
    call mem%allocate(bandcounter, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    bandcounter = 0
    do i = 1, qp%n_full_point
        gi2 = i
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
            om1 = dr%aq(gi1)%omega(b1)
            om2 = dr%aq(gi2)%omega(b2)
            if (om1 .gt. omthres .and. om2 .gt. omthres) then
                select case (integrationtype)
                case (1)
                    sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_adjustment
                case (2)
                    sigma = qp%smearingparameter(dr%aq(gi2)%vel(:, b2), dr%default_smearing(b2), smearing_adjustment)
                end select

                if (abs(om1 - om2) .lt. thres*sigma) then
                    bandcounter(b1, b2) = bandcounter(b1, b2) + 1
                end if
            end if
        end do
        end do
    end do

    ! Allocate storage
    do b1 = 1, dr%n_mode
    do b2 = 1, dr%n_mode
        scq%band(b1, b2)%n = bandcounter(b1, b2)
        if (scq%band(b1, b2)%n .gt. 0) then
            allocate (scq%band(b1, b2)%e(scq%band(b1, b2)%n))
        end if
    end do
    end do

    ! Count again and store
    bandcounter = 0
    do i = 1, qp%n_full_point
        gi2 = i
        ! It should not bounce to itself, perhaps
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
            om1 = dr%aq(gi1)%omega(b1)
            om2 = dr%aq(gi2)%omega(b2)
            if (om1 .gt. omthres .and. om2 .gt. omthres) then
                select case (integrationtype)
                case (1)
                    sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_adjustment
                case (2)
                    sigma = qp%smearingparameter(dr%aq(gi2)%vel(:, b2), dr%default_smearing(b2), smearing_adjustment)
                end select

                if (abs(om1 - om2) .lt. thres*sigma) then
                    bandcounter(b1, b2) = bandcounter(b1, b2) + 1
                    ii = bandcounter(b1, b2)
                    deltafunction = lo_gauss(om1, om2, sigma)
                    scq%band(b1, b2)%e(ii)%deltafunction = deltafunction/qp%n_full_point
                    scq%band(b1, b2)%e(ii)%gi2 = gi2
                end if
            end if
        end do
        end do
    end do
    call mem%deallocate(bandcounter, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine
