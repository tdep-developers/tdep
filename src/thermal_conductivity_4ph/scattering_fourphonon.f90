
subroutine compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, mcg, rng, &
                                         g0, integrationtype, smearing, mw, mem)
    !> The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> The monte-carlo grid
    type(lo_montecarlo_grid), intent(in) :: mcg
    !> The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    !> The linewidth for this mode
    real(r8), intent(inout) :: g0
    !> what kind of integration are we doing
    integer, intent(in) :: integrationtype
    !> The smearing width
    real(r8), intent(in) :: smearing
    !> Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3, egv4
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2, evp3
    !> The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull1, qgridfull2
    !> The qpoints in cartesian coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4
    !> The complex scattering amplitude
    complex(r8) :: c0
    !> Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, psisq
    ! The gaussian integration width
    real(r8) :: sigma
    !> Stuff for the linewidths
    real(r8) :: n2, n3, n4, n2p, n3p, n4p, plf0, plf1, plf2, plf3
    !> Integers for do loops
    integer :: q1, q2, q3, q4, q2p, q3p, q4p, b1, b2, b3, b4, qi, qj, i
    !> Is the quartet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult0, mult1, mult2, mult3
    !> All the prefactors for the scattering
    real(r8) :: f0, f1, f2, f3, f4, f5, f6, f7, fall
    !> The reducible triplet corresponding to the currently computed quartet
    integer, dimension(:, :), allocatable :: red_quartet

    real(r8), dimension(:, :), allocatable :: od_terms

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv4, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(od_terms, [qp%n_full_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    od_terms = 0.0_r8

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1)/sqrt(om1)

    ! Prepare the grid for the monte-carlo average
    call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mcg%generate_grid(qgridfull1, rng)
    call mcg%generate_grid(qgridfull2, rng)

    compute_loop: do qi = 1, mcg%npoints
    do qj = 1, mcg%npoints
        q2 = qgridfull1(qi)
        q3 = qgridfull2(qj)
        if (q3 .lt. q2) cycle
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, mcg%full_dims)
        if (q4 .lt. q3) cycle

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, red_quartet, mw, mem)
        if (isred) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        qv4 = qp%ap(q4)%r
        call pretransform_phi4(fcf, qv2, qv3, qv4, ptf)

        ! We can already take care of the multiplicity caused by permutation
        if (q2 .eq. q3 .and. q3 .eq. q4) then
            mult0 = 1.0_r8
            mult1 = 1.0_r8
            mult2 = 1.0_r8
            mult3 = 1.0_r8
        else if ((q2 .ne. q3 .and. q3 .eq. q4) .or. &
                 (q3 .ne. q2 .and. q2 .eq. q4) .or. &
                 (q4 .ne. q2 .and. q2 .eq. q3)) then
            mult0 = 3.0_r8
            mult1 = 1.0_r8
            mult2 = 1.0_r8
            mult3 = 1.0_r8
        else
            mult0 = 6.0_r8
            mult1 = 2.0_r8
            mult2 = 2.0_r8
            mult3 = 2.0_r8
        end if

        do b2 = 1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
            n2p = n2 + 1.0_r8
            egv2 = dr%aq(q2)%egv(:, b2)/sqrt(om2)

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3 = 1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                n3 = sr%be(qp%ap(q3)%irreducible_index, b3)
                n3p = n3 + 1.0_r8
                egv3 = dr%aq(q3)%egv(:, b3)/sqrt(om3)

                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                do b4 = 1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    n4 = sr%be(qp%ap(q4)%irreducible_index, b4)
                    n4p = n4 + 1.0_r8
                    egv4 = dr%aq(q4)%egv(:, b4)/sqrt(om4)

                    select case (integrationtype)
                    case (1)
                        sigma = smearing*lo_frequency_THz_to_Hartree
                    case (2)
                        sigma = sqrt(sr%sigsq(q1, b1) + &
                                     sr%sigsq(qp%ap(q2)%irreducible_index, b2) + &
                                     sr%sigsq(qp%ap(q3)%irreducible_index, b3) + &
                                     sr%sigsq(qp%ap(q4)%irreducible_index, b4))
                    end select

                    evp3 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                    evp3 = conjg(evp3)
                    c0 = dot_product(evp3, ptf)
                    psisq = fourphonon_prefactor*abs(c0*conjg(c0))*mcg%weight**2

                    ! Prefactors, only the Bose-Einstein distributions
                    plf0 = n2p*n3p*n4p - n2*n3*n4
                    plf1 = 3.0_r8*n2*n3p*n4p - n2p*n3*n4
                    plf2 = 3.0_r8*n3*n2p*n4p - n3p*n2*n4
                    plf3 = 3.0_r8*n4*n3p*n2p - n4p*n3*n2

                    ! Prefactors, including the matrix elements and dirac
                    f0 = mult0*psisq*plf0*(lo_gauss(om1, om2 + om3 + om4, sigma) - lo_gauss(om1, -om2 - om3 - om4, sigma))
                    f1 = mult1*psisq*plf1*(lo_gauss(om1, -om2 + om3 + om4, sigma) - lo_gauss(om1, om2 - om3 - om4, sigma))
                    f2 = mult2*psisq*plf2*(lo_gauss(om1, -om3 + om2 + om4, sigma) - lo_gauss(om1, om3 - om2 - om4, sigma))
                    f3 = mult3*psisq*plf3*(lo_gauss(om1, -om4 + om3 + om2, sigma) - lo_gauss(om1, om4 - om3 - om2, sigma))

                    fall = f0 + f1 + f2 + f3

                    ! Add everything to the linewidth
                    do i = 1, size(red_quartet, 2)
                        g0 = g0 + fall

                        q2p = red_quartet(1, i)
                        q3p = red_quartet(2, i)
                        q4p = red_quartet(3, i)

                        od_terms(q2p, b2) = od_terms(q2p, b2) + 2.0_r8*fall*om2/om1
                        od_terms(q3p, b3) = od_terms(q3p, b3) + 2.0_r8*fall*om3/om1
                        od_terms(q4p, b4) = od_terms(q4p, b4) + 2.0_r8*fall*om4/om1
                    end do
                end do
            end do
        end do
    end do
    end do compute_loop

    ! Now we can symmetrize the off-diagonal contribution
    ! This can be done in a way to put a value to for mode that have been skipped by the Monte-Carlo !
    symmetrize_and_distribute: block
        integer, dimension(dr%n_mode) :: nn
        !> To hold the q-point in reduced coordinates
        real(r8), dimension(3) :: qv2p
        !> To get the index of the new triplet on the fft_grid
        integer, dimension(3) :: gi
        !> Some buffer
        real(r8), dimension(dr%n_mode) :: buf_xi
        !> Some integers for the do loop
        integer :: j, k, i2

        ! Let's average the off diagonal term
        allq2: do q2 = 1, qp%n_full_point
            buf_xi = 0.0_r8
            nn = 0
            do j = 1, qp%ip(q1)%n_invariant_operation
                k = qp%ip(q1)%invariant_operation(j)
                select type (qp); type is (lo_fft_mesh)
                    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
                    qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
                    if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle
                    gi = qp%index_from_coordinate(qv2p)
                    q2p = qp%gridind2ind(gi(1), gi(2), gi(3)) ! this is R*q'
                end select
                if (q2p .lt. q2) cycle allq2  ! If q2p < q2, we already did this guy
                do b2 = 1, dr%n_mode
                    if (od_terms(q2p, b2) .gt. 0.0_r8) then
                        nn(b2) = nn(b2) + 1
                        buf_xi(b2) = buf_xi(b2) + od_terms(q2p, b2)
                    end if
                end do
            end do
            do j = 1, qp%ip(q1)%n_invariant_operation
                k = qp%ip(q1)%invariant_operation(j)
                select type (qp); type is (lo_fft_mesh)
                    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
                    qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
                    if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle
                    gi = qp%index_from_coordinate(qv2p)
                    q2p = qp%gridind2ind(gi(1), gi(2), gi(3)) ! this is R*q'
                end select
                do b2 = 1, dr%n_mode
                    if (nn(b2) .eq. 0) cycle
                    od_terms(q2p, b2) = buf_xi(b2)/real(nn(b2), r8)
                end do
            end do
        end do allq2
        ! And now we distribute the off-diagonal terms on the scattering matrix
        do q2 = 1, qp%n_full_point
            do b2 = 1, dr%n_mode
                i2 = (q2 - 1)*dr%n_mode + b2
                sr%Xi(il, i2) = sr%Xi(il, i2) + od_terms(q2, b2)
            end do
        end do
    end block symmetrize_and_distribute

    ! And we can deallocate everything
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(od_terms, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

subroutine quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, red_quartet, mw, mem)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2, q3 and q4 are full
    integer, intent(in) :: q1, q2, q3, q4
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> The equivalent triplet
    integer, dimension(:, :), allocatable, intent(out) :: red_quartet
    !> Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The new-qpoints and the temporary invariant triplet
    integer, dimension(:, :), allocatable :: newqp, newqp_sort
    !> The number of invariant operation for q1
    integer :: n

    n = qp%ip(q1)%n_invariant_operation

    ! We will need this to keep equivalent point
    call mem%allocate(newqp, [3, n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(newqp_sort, [3, n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    newqp = -lo_hugeint
    newqp_sort = -lo_hugeint

    ! The first part is to generate all qpoint invariant in little star of q1
    get_equivalent: block
        !> To hold the q-point in reduced coordinates
        real(r8), dimension(3) :: qv2, qv3, qv4, qv2p, qv3p, qv4p
        !> To get the index of the new triplet on the fft_grid
        integer, dimension(3) :: gi
        !> The new triplet after the operation
        integer, dimension(3) :: qpp
        !> Integers for the do loops
        integer :: j, k

        ! First get the reciprocal lattice vectors, in reduce coordinates
        qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
        qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)
        qv4 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q4)%r)

        isred = .false.
        ! Let's try all operations that leaves q1 invariant
        do j = 1, n
            k = qp%ip(q1)%invariant_operation(j)
            qpp = -lo_hugeint
            select type (qp); type is (lo_fft_mesh)
                ! Transform q2 and check if it's on the grid
                qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
                if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle

                ! Transform q3 and check if it's on the grid
                qv3p = lo_operate_on_vector(uc%sym%op(k), qv3, reciprocal=.true., fractional=.true.)
                if (qp%is_point_on_grid(qv3p) .eqv. .false.) cycle

                ! Transform q4 and check if it's on the grid
                qv4p = lo_operate_on_vector(uc%sym%op(k), qv4, reciprocal=.true., fractional=.true.)
                if (qp%is_point_on_grid(qv4p) .eqv. .false.) cycle

                ! If everything is on the grid, get the location of each point
                gi = qp%index_from_coordinate(qv2p)
                qpp(1) = qp%gridind2ind(gi(1), gi(2), gi(3))
                gi = qp%index_from_coordinate(qv3p)
                qpp(2) = qp%gridind2ind(gi(1), gi(2), gi(3))
                gi = qp%index_from_coordinate(qv4p)
                qpp(3) = qp%gridind2ind(gi(1), gi(2), gi(3))
            end select
            ! We need to keep the symmetry equivalent point in the order they came in
            newqp(:, j) = qpp
            ! We also need them sorted, to check for redundancy and reducibility
            call lo_qsort(qpp)
            newqp_sort(:, j) = qpp
            if (qpp(1) .gt. q2 .or. qpp(2) .gt. q3) isred = .true.
        end do

        if (minval(newqp) .lt. 0) then
            do j = 1, size(newqp, 2)
                if (any(newqp(:, j) .lt. 0)) newqp(:, j) = [q2, q3, q4]
                if (any(newqp(:, j) .lt. 0)) newqp_sort(:, j) = [q2, q3, q4]
            end do
        end if
    end block get_equivalent

    ! Then we just get rid of redundant point in the list
    sort_reducible: block
        !> Integers for the loops
        integer :: j, k, ctr, ctr2

        ! Now we have the same problem of permutation as in the third order
        ! So first, we count the quartet equivalent by permutation, a little bit harder
        ctr = 0
        do j = 1, n
            ctr2 = 0
            do k = j, n
                if (k .eq. j) cycle
                if (all(newqp_sort(:, j) .eq. newqp_sort(:, k))) ctr2 = ctr2 + 1
            end do
            if (ctr2 .eq. 0) ctr = ctr + 1
        end do

        ! Now we can create the list of equivalent quartet with permutation removed
        allocate (red_quartet(3, ctr))
        ctr = 0
        do j = 1, n
            ctr2 = 0
            do k = j, n
                if (k .eq. j) cycle
                if (all(newqp_sort(:, j) .eq. newqp_sort(:, k))) ctr2 = ctr2 + 1
            end do
            if (ctr2 .eq. 0) then
                ctr = ctr + 1
                red_quartet(1, ctr) = newqp(1, j)
                red_quartet(2, ctr) = newqp(2, j)
                red_quartet(3, ctr) = newqp(3, j)
            end if
        end do

    end block sort_reducible
    call mem%deallocate(newqp, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(newqp_sort, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
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
