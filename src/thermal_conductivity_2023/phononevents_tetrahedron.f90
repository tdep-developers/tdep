
!> Tetrahedron integration weights for threephonon scattering, from one q-point
subroutine threephonon_tetrahedron_oneqp(qp, dr, scq, gi1, mem)
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> q-point in question
    type(lo_3phqp2), intent(inout) :: scq
    !> grid-index in question
    integer, intent(in) :: gi1
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer :: i, j, ii, jj, b1, b2, b3
    integer :: gi2, gi3, i_gamma
    integer, dimension(3) :: dims
    integer, dimension(4) :: tetqpoints
    real(r8), dimension(:, :, :, :), allocatable :: qpwp, qpwm
    real(r8), dimension(:, :), allocatable :: omr2, omr3
    real(r8), dimension(:), allocatable :: omr1
    real(r8), dimension(4) :: cval_p, cval_m, wts1, wts2
    real(r8) :: om1, omthres, minc_m, minc_p, maxc_m, maxc_p, wthres, sigmatol
    real(r8) :: sigma

    ! Some temporary space
    call mem%allocate(omr1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr2, [4, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr3, [4, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qpwp, [dr%n_mode, dr%n_mode, dr%n_mode, qp%n_full_point], persistent=.false., &
                      scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qpwm, [dr%n_mode, dr%n_mode, dr%n_mode, qp%n_full_point], persistent=.false., &
                      scalable=.false., file=__FILE__, line=__LINE__)
    omr1 = 0.0_r8
    omr2 = 0.0_r8
    omr3 = 0.0_r8
    qpwp = 0.0_r8
    qpwm = 0.0_r8

    sigma = 1.0_r8*lo_frequency_THz_to_Hartree

    ! smallest frequency allowed, should allow everything except acoustic modes at Gamma
    omthres = dr%omega_min*0.5_r8
    ! smallest weight to care about
    wthres = 1E-30_r8
    ! Dimensions of q-point grid
    select type (qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        write (*, *) 'Really need an fft grid for this'
        stop
    end select
    ! The reference q-point, q1 in q1+q2+q3=G
    scq%gi1 = gi1
    ! Some tolerances
    sigmatol = maxval(dr%default_smearing)*1E-4_r8

    qpwp = 0.0_r8
    qpwm = 0.0_r8
    omr1 = dr%aq(gi1)%omega
    tetloop: do i = 1, qp%n_full_tet
        ! Grab frequencies
        do j = 1, 4
            gi2 = qp%at(i)%full_index(j)
            gi3 = fft_third_grid_index(gi1, gi2, dims)
            omr2(j, :) = dr%aq(gi2)%omega
            omr3(j, :) = dr%aq(gi3)%omega
            tetqpoints(j) = gi2
        end do
        ! Add integration weights
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
        do b3 = 1, dr%n_mode
            minc_p = lo_huge
            maxc_p = -lo_huge
            minc_m = lo_huge
            maxc_m = -lo_huge
            do j = 1, 4
                cval_p(j) = omr3(j, b3) - omr2(j, b2)   ! plus events
                cval_m(j) = omr3(j, b3) + omr2(j, b2)   ! minus events
                minc_p = min(minc_p, cval_p(j))
                minc_m = min(minc_m, cval_m(j))
                maxc_p = max(maxc_p, cval_p(j))
                maxc_m = max(maxc_m, cval_m(j))
            end do
            ! add a small delta for safety, to avid pathological cases
            minc_p = minc_p - sigmatol
            minc_m = minc_m - sigmatol
            maxc_p = maxc_p + sigmatol
            maxc_m = maxc_m + sigmatol
            ! check if it's relevant
            om1 = omr1(b1)
            if (om1 .gt. minc_p .and. om1 .lt. maxc_p) then !.and. om1 .gt. omthres ) then
                wts1 = lo_LV_tetrahedron_weights(cval_p, om1, lo_freqtol, dr%default_smearing(b1))
                !wts1=lo_LV_tetrahedron_weights(cval_p,om1,lo_freqtol,sigma)
                qpwp(b1, b2, b3, tetqpoints) = qpwp(b1, b2, b3, tetqpoints) + wts1*qp%at(i)%integration_weight
            end if
            if (om1 .gt. minc_m .and. om1 .lt. maxc_m) then !.and. om1 .gt. omthres ) then
                wts1 = lo_LV_tetrahedron_weights(cval_m, om1, lo_freqtol, dr%default_smearing(b1))
                !wts1=lo_LV_tetrahedron_weights(cval_m,om1,lo_freqtol,sigma)
                qpwm(b1, b2, b3, tetqpoints) = qpwm(b1, b2, b3, tetqpoints) + wts1*qp%at(i)%integration_weight
            end if
        end do
        end do
        end do
    end do tetloop

    ! Probably a good idea to locate gamma
    i_gamma = -1
    do i = 1, qp%n_full_point
        if (lo_sqnorm(qp%ap(i)%r) .lt. lo_sqtol) then
            i_gamma = i
            exit
        end if
    end do
    if (i_gamma .lt. 0) call lo_stop_gracefully(['FFT mesh does not contain gamma'], lo_exitcode_symmetry, __FILE__, __LINE__)

    ! There should be no scattering allowed with gamma
    do b3 = 1, dr%n_mode
    do b2 = 1, dr%n_mode
    do b1 = 1, dr%n_mode
        if (dr%aq(i_gamma)%omega(b1) .lt. lo_freqtol .or. &
            dr%aq(i_gamma)%omega(b2) .lt. lo_freqtol .or. &
            dr%aq(i_gamma)%omega(b3) .lt. lo_freqtol) then
            qpwp(b1, b2, b3, i_gamma) = 0.0_r8
            qpwm(b1, b2, b3, i_gamma) = 0.0_r8
        end if
    end do
    end do
    end do
    qpwp(1:3, 1:3, 1:3, i_gamma) = 0.0_r8
    qpwm(1:3, 1:3, 1:3, i_gamma) = 0.0_r8

    ! Store weights per q-point
    do b1 = 1, dr%n_mode
    do b2 = 1, dr%n_mode
    do b3 = 1, dr%n_mode
        ! count q-points
        ii = 0
        jj = 0
        do i = 1, qp%n_full_point
            if (abs(qpwp(b1, b2, b3, i)) .gt. wthres) then
                ii = ii + 1
            end if
            if (abs(qpwm(b1, b2, b3, i)) .gt. wthres) then
                jj = jj + 1
            end if
        end do
        ! make some space
        scq%plus(b1, b2, b3)%n = ii
        scq%minus(b1, b2, b3)%n = jj
        if (scq%plus(b1, b2, b3)%n .gt. 0) then
            allocate (scq%plus(b1, b2, b3)%e(ii))
        end if
        if (scq%minus(b1, b2, b3)%n .gt. 0) then
            allocate (scq%minus(b1, b2, b3)%e(jj))
        end if
        ! store weights and indices
        ii = 0
        jj = 0
        do i = 1, qp%n_full_point
            if (abs(qpwp(b1, b2, b3, i)) .gt. wthres) then
                ii = ii + 1
                gi2 = i
                scq%plus(b1, b2, b3)%e(ii)%gi2 = gi2
                scq%plus(b1, b2, b3)%e(ii)%gi3 = fft_third_grid_index(gi1, gi2, dims)
                scq%plus(b1, b2, b3)%e(ii)%deltafunction = qpwp(b1, b2, b3, i)
            end if
            if (abs(qpwm(b1, b2, b3, i)) .gt. wthres) then
                jj = jj + 1
                gi2 = i
                scq%minus(b1, b2, b3)%e(jj)%gi2 = gi2
                scq%minus(b1, b2, b3)%e(jj)%gi3 = fft_third_grid_index(gi1, gi2, dims)
                scq%minus(b1, b2, b3)%e(jj)%deltafunction = qpwm(b1, b2, b3, i)
            end if
        end do
    end do
    end do
    end do

    ! Cleanup
    call mem%deallocate(omr1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qpwp, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qpwm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Tetrahedron integration weights for isotope scattering, from one q-point
subroutine iso_tetrahedron_oneqp(qp, dr, scq, gi1, mem)
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersoins
    type(lo_phonon_dispersions), intent(in) :: dr
    !> current q-point
    type(lo_iso2), intent(inout) :: scq
    !> current grid-index
    integer, intent(in) :: gi1
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), parameter :: thres_weight = lo_tiny
    integer, dimension(4) :: tetqpoints
    integer :: gi2
    integer :: i, j, ii, b1, b2, i_gamma
    real(r8), dimension(4) :: cval, wts
    real(r8) :: om1, omthres, minc, maxc, sigmatol
    real(r8), dimension(:, :, :), allocatable :: qpw
    real(r8), dimension(:, :), allocatable :: omr2
    real(r8), dimension(:), allocatable :: omr1

    call mem%allocate(omr1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(omr2, [4, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qpw, [dr%n_mode, dr%n_mode, qp%n_full_point], persistent=.false., &
                      scalable=.false., file=__FILE__, line=__LINE__)
    omr1 = 0.0_r8
    omr2 = 0.0_r8
    qpw = 0.0_r8
    omthres = dr%omega_min*0.2_r8
    sigmatol = maxval(dr%default_smearing)/1000.0_r8

    ! Calculate all integration weights for the relevant tetrahedra
    scq%gi1 = gi1
    qpw = 0.0_r8
    omr1 = dr%aq(gi1)%omega
    tetloop: do i = 1, qp%n_full_tet
        ! Grab frequencies
        do j = 1, 4
            gi2 = qp%at(i)%full_index(j)
            omr2(j, :) = dr%aq(gi2)%omega
            tetqpoints(j) = gi2
        end do

        ! Add integration weights
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
            minc = lo_huge
            maxc = -lo_huge
            do j = 1, 4
                cval(j) = omr2(j, b2)
                minc = min(minc, cval(j))
                maxc = max(maxc, cval(j))
            end do
            ! add a small delta for safety
            minc = minc - 3*sigmatol
            maxc = maxc + 3*sigmatol
            ! check if it's relevant
            om1 = omr1(b1)
            if (om1 .gt. minc .and. om1 .lt. maxc .and. om1 .gt. omthres) then
                ! sum up integration weights
                wts = lo_LV_tetrahedron_weights(cval, om1, lo_freqtol, dr%default_smearing(b1))
                qpw(b1, b2, tetqpoints) = qpw(b1, b2, tetqpoints) + wts*qp%at(i)%integration_weight
            end if
        end do
        end do
    end do tetloop

    ! Set the integration weights for the acoustic branches to zero
    i_gamma = -1
    do i = 1, qp%n_full_point
        if (lo_sqnorm(qp%ap(i)%r) .lt. lo_sqtol) then
            i_gamma = i
            exit
        end if
    end do
    if (i_gamma .lt. 0) call lo_stop_gracefully(['FFT mesh does not contain gamma'], lo_exitcode_symmetry, __FILE__, __LINE__)

    do b2 = 1, dr%n_mode
    do b1 = 1, dr%n_mode
        if (dr%aq(i_gamma)%omega(b1) .lt. lo_freqtol .or. dr%aq(i_gamma)%omega(b2) .lt. lo_freqtol) then
            qpw(b1, b2, i_gamma) = 0.0_r8
        end if
    end do
    end do

    ! Store weights per q-point, and only the relevant q-points
    do b1 = 1, dr%n_mode
    do b2 = 1, dr%n_mode
        ! count q-points
        ii = 0
        do i = 1, qp%n_full_point
            if (abs(qpw(b1, b2, i)) .gt. thres_weight) then
                ii = ii + 1
            end if
        end do
        ! make some space
        scq%band(b1, b2)%n = ii
        if (scq%band(b1, b2)%n .gt. 0) then
            allocate (scq%band(b1, b2)%e(ii))
        end if
        ! store weights and indices
        ii = 0
        do i = 1, qp%n_full_point
            if (abs(qpw(b1, b2, i)) .gt. thres_weight) then
                ii = ii + 1
                scq%band(b1, b2)%e(ii)%gi2 = i
                scq%band(b1, b2)%e(ii)%deltafunction = qpw(b1, b2, i)
            end if
        end do
    end do
    end do

    call mem%deallocate(omr1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(omr2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qpw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine
