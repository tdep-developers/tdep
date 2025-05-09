
subroutine compute_threephonon_scattering(il, sr, qp, dr, uc, fct, mcg, rng, &
                                          g0, integrationtype, smearing, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Third order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
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
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2
    !> The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull
    !> The qpoints and the dimension of the qgrid
    real(r8), dimension(3) :: qv2, qv3
    !> Frequencies, bose-einstein occupation and scattering strength and some other buffer
    real(r8) :: sigma, om1, om2, om3, n2, n3, psisq, f0, f1, plf0, plf1, perm
    !> The complex threephonon matrix element
    complex(r8) :: c0
    !> Integers for do loops and counting
    integer :: qi, q1, q2, q3, q2p, q3p, b1, b2, b3, i
    !> Is the triplet irreducible ?
    logical :: isred
    !> The reducible triplet corresponding to the currently computed triplet
    integer, dimension(:, :), allocatable :: red_triplet
    !> buff to keep the off diagonal scattering matrix elements
    real(r8), dimension(:, :), allocatable :: od_terms
    !> The approximated dirac
    real(r8) :: d0, d1

    real(r8), dimension(3, 3) :: reclat, wig
    real(r8), dimension(3) :: w, sigsig

    logical :: isok

    do i=1, 3
        reclat(:, i) = uc%reciprocal_latticevectors(:, i) / mcg%full_dims(i)
    end do

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(od_terms, [qp%n_full_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    od_terms = 0.0_r8

    ! Already set some values for mode (q1, b1)
    q1 = sr%my_qpoints(il)
    b1 = sr%my_modes(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1)/sqrt(om1)

    call mem%allocate(qgridfull, mcg%npoints, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mcg%generate_grid(qgridfull, rng)

    compute_loop: do qi = 1, mcg%npoints
        q2 = qgridfull(qi)
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, mcg%full_dims)
        if (q3 .lt. q2) cycle
        call triplet_is_irreducible(qp, uc, q1, q2, q3, isred, red_triplet, mw, mem)
        if (isred) cycle

        ! Let's compute the multiplicity due to permutation now
        if (q2 .eq. q3) then
            perm = 1.0_r8
        else
            perm = 2.0_r8
        end if

        ! This get the ifc3 in Fourier space, but not on phonons
        call pretransform_phi3(fct, qp%ap(q2)%r, qp%ap(q3)%r, ptf)
        do b2 = 1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2)/sqrt(om2)

            ! This is the multiplication of eigv of phonons 1 and 2
            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3 = 1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                egv3 = dr%aq(q3)%egv(:, b3)/sqrt(om3)

                call get_dirac(sr, qp, dr, q1, q2, q3, b1, b2, b3, integrationtype, d0, d1, isok)

                if (.not. isok) cycle

                ! This is the multiplication of eigv of phonons 1 and 2 and now 3
                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                evp2 = conjg(evp2)
                ! And with this, we have the scattering matrix element coming in
                c0 = dot_product(evp2, ptf)
                psisq = threephonon_prefactor*abs(c0*conjg(c0))*mcg%weight

                ! Let's get the Bose-Einstein distributions
                n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
                n3 = sr%be(qp%ap(q3)%irreducible_index, b3)

                ! The prefactor for the scattering
                plf0 = n2 - n3
                plf1 = n2 + n3 + 1.0_r8
!               f0 = perm*psisq*plf0*(lo_gauss(om1, -om2 + om3, sigma) - lo_gauss(om1, om2 - om3, sigma))
!               f1 = perm*psisq*plf1*(lo_gauss(om1, om2 + om3, sigma) - lo_gauss(om1, -om2 - om3, sigma))
                f0 = perm*psisq*plf0*d0
                f1 = perm*psisq*plf1*d1

                ! And we add everything for each triplet equivalent to the one we are actually computing
                do i = 1, size(red_triplet, 2)
                    ! We accumulate the linewidth
                    g0 = g0 + f0 + f1

                    ! And also the off-diagonal part of the scattering matrix
                    q2p = red_triplet(1, i)
                    q3p = red_triplet(2, i)
                    od_terms(q2p, b2) = od_terms(q2p, b2) + 2.0_r8*(f0 + f1)*om2/om1
                    od_terms(q3p, b3) = od_terms(q3p, b3) + 2.0_r8*(f0 + f1)*om3/om1
                end do
            end do
        end do
    end do compute_loop

    ! Now we can symmetrize the off-diagonal contribution
    ! For this, we only compute the average values from q-point on the Monte-Carlo grid.
    ! But we distribute averaged entries of the scattering matrix for every equivalent points
    symmetrize_and_distribute: block
        !> To keep track of the number of mode we actually add
        integer, dimension(dr%n_mode) :: nn
        !> To hold the q-point in reduced coordinates
        real(r8), dimension(3) :: qv2p
        !> To get the index of the new triplet on the fft_grid
        integer, dimension(3) :: gi
        !> Some buffer
        real(r8), dimension(dr%n_mode) :: buf_xi
        !> Some integers for the do loops and stuff
        integer :: j, k, i2

        ! Let's average the off diagonal term
        allq2: do qi = 1, mcg%npoints
            q2 = qgridfull(qi)
            buf_xi = 0.0_r8
            nn = 0
            ! First we get the average value
            do j = 1, qp%ip(q1)%n_invariant_operation
                k = qp%ip(q1)%invariant_operation(j)
                ! First we generate q2'
                select type (qp); type is (lo_fft_mesh)
                    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
                    qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
                    if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle
                    gi = qp%index_from_coordinate(qv2p)
                    q2p = qp%gridind2ind(gi(1), gi(2), gi(3)) ! this is R*q'
                end select
                ! If q2p < q2, we already did this guy
                if (q2p .lt. q2) cycle allq2
                ! Accumulate values for each equivalent term, but avoid those which where not computed
                do b2 = 1, dr%n_mode
                    if (od_terms(q2p, b2) .gt. 0.0_r8) then
                        nn(b2) = nn(b2) + 1
                        buf_xi(b2) = buf_xi(b2) + od_terms(q2p, b2)
                    end if
                end do
            end do
            ! And now we can distribute it
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

        ! And now we add things, with the normalization
        od_terms = od_terms
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
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(od_terms, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    deallocate (red_triplet)

    contains
subroutine get_dirac(sr, qp, dr, q1, q2, q3, b1, b2, b3, integrationtype, d0, d1, isok)
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(in) :: sr
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> Harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> The qpoints
    integer, intent(in) :: q1, q2, q3
    !> The modes
    integer, intent(in) :: b1, b2, b3
    !> THe integration type
    integer, intent(in) :: integrationtype
    !> The approximated dirac
    real(r8), intent(out) :: d0, d1

    logical, intent(out) :: isok

    !> The frequencies
    real(r8) :: om1, om2, om3
    ! Buffer to contains all the sigmas
    real(r8), dimension(3, 3) :: allsig
    !> An integer for the do loop
    integer :: i

    integer :: j

    select case (integrationtype)
    ! Gaussian smearing with fixed width
    case (1)
        sigma = smearing * lo_frequency_THz_to_Hartree
    ! Adaptive Gaussian smearing with smeared frequency approach
    case (2)
        sigma = sqrt(sr%sigsq(q1, b1) + &
                     sr%sigsq(qp%ap(q2)%irreducible_index, b2) + &
                     sr%sigsq(qp%ap(q3)%irreducible_index, b3))
    ! Adaptive Gaussian smearing with velocity difference approach
    case (6)
        allsig = 0.0_r8
        allsig(:, 1) = matmul(dr%iq(q1)%vel(:, b1) - dr%aq(q2)%vel(:, b2), sr%reclat)**2
        allsig(:, 2) = matmul(dr%iq(q1)%vel(:, b1) - dr%aq(q3)%vel(:, b3), sr%reclat)**2
        allsig(:, 3) = matmul(dr%aq(q2)%vel(:, b2) - dr%aq(q3)%vel(:, b3), sr%reclat)**2

        sigma = 0.0_r8
        do i=1, 3
            sigma = sigma + sqrt(maxval(allsig(:, i)) * 0.5_r8) / 3.0_r8
        end do

!       sigma = sqrt(maxval(matmul(dr%aq(q2)%vel(:, b2) - dr%aq(q3)%vel(:, b3), sr%reclat)**2) * 0.5_r8)
    end select

    om1 = dr%iq(q1)%omega(b1)
    om2 = dr%aq(q2)%omega(b2)
    om3 = dr%aq(q3)%omega(b3)
    d0 = 0.0_r8
    d1 = 0.0_r8
    j = 0
    isok = .false.
    ! We cannot have interaction with same mode (or degenerate)
    if (abs(om1 - om2) .gt. lo_freqtol .and. &
        abs(om1 - om3) .gt. lo_freqtol .and. &
        abs(om2 - om3) .gt. lo_freqtol) then
!       d0 = lo_gauss(om1, -om2 + om3, sigma) - lo_gauss(om1, om2 - om3, sigma)
!       d1 = lo_gauss(om1, om2 + om3, sigma) - lo_gauss(om1, -om2 - om3, sigma)

        if (abs(om1 + om2 - om3) .lt. 4.0_r8 * sigma) then
            d0 = d0 + lo_gauss(om1, -om2 + om3, sigma)
            j = j + 1
        end if
        if (abs(om1 - om2 + om3) .lt. 4.0_r8 * sigma) then
            d0 = d0 - lo_gauss(om1, om2 - om3, sigma)
            j = j + 1
        end if
        if (abs(om1 - om2 - om3) .lt. 4.0_r8 * sigma) then
            d1 = d1 + lo_gauss(om1, om2 + om3, sigma)
            j = j + 1
        end if
        if (abs(om1 + om2 + om3) .lt. 4.0_r8 * sigma) then
            d1 = d1 - lo_gauss(om1, -om2 - om3, sigma)
            j = j + 1
        end if
    end if
    if (j .gt. 0) isok = .true.
end subroutine
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

subroutine triplet_is_irreducible(qp, uc, q1, q2, q3, isred, red_triplet, mw, mem)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2 and q3 are full
    integer, intent(in) :: q1, q2, q3
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> The equivalent triplet
    integer, dimension(:, :), allocatable, intent(out) :: red_triplet
    !> Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The new-qpoints and the temporary invariant triplet
    integer, dimension(:, :), allocatable :: newqp, newqp_sort
    !> An integer to keep size of thing
    integer :: n

    ! Little explanation of what is happening here
    ! First I generate all triplet point that are equivalent to the input in the little star of q1
    ! Also, I return the fact that it's irreducible if possible by permutation
    ! Then I remove the doublet (including permutation) in order to return the list of equivalent triplet
    ! With this, I can generate equivalent triplet even if the grid does not respect the space group of the crystal
    ! You should make your grid respect it, but you never now what people do

    n = qp%ip(q1)%n_invariant_operation

    ! We will need this to keep equivalent point
    call mem%allocate(newqp, [2, n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(newqp_sort, [2, n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    newqp = -lo_hugeint
    newqp_sort = -lo_hugeint

    ! The first part is to generate all qpoint invariant in little star of q1
    get_equivalent: block
        !> To hold the q-point in reduced coordinates
        real(r8), dimension(3) :: qv2, qv3, qv2p, qv3p
        !> To get the index of the new triplet on the fft_grid
        integer, dimension(3) :: gi
        !> The new triplet after the operation
        integer, dimension(2) :: qpp
        !> Integers for the do loops
        integer :: j, k

        ! First get the q-points in reduce coordinates
        qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
        qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)

        isred = .false.
        ! Let's try all operations that leaves q1 invariant
        do j = 1, n
            k = qp%ip(q1)%invariant_operation(j)
            qpp = -lo_hugeint
            select type (qp); type is (lo_fft_mesh)
                ! Rotate q2 and look if it's the on grid
                qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
                if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle

                ! Rotate q3 and look if it's the on grid
                qv3p = lo_operate_on_vector(uc%sym%op(k), qv3, reciprocal=.true., fractional=.true.)
                if (qp%is_point_on_grid(qv3p) .eqv. .false.) cycle

                ! If everything is on the grid, get the index of each point
                gi = qp%index_from_coordinate(qv2p)
                qpp(1) = qp%gridind2ind(gi(1), gi(2), gi(3))
                gi = qp%index_from_coordinate(qv3p)
                qpp(2) = qp%gridind2ind(gi(1), gi(2), gi(3))
            end select
            ! If we got here, it means that the new triplet is on the grid
            ! We need to keep the symmetry equivalent point in the order they came in
            newqp(:, j) = qpp
            ! We also need them sorted, to check for redundancy and reducibility
            call lo_qsort(qpp)
            newqp_sort(:, j) = qpp
            ! Now I check if I can reduce this triplet
            if (qpp(1) .gt. q2) isred = .true.
        end do

        ! For stability, replace points not on the grid with starting q-point
        ! Should not be needed if the grid respect the symmetries of the lattice though
        if (minval(newqp) .lt. 0) then
            do j = 1, n
                if (any(newqp(:, j) .lt. 0)) then
                    newqp(:, j) = [q2, q3]
                    newqp_sort(:, j) = [q2, q3]
                end if
            end do
        end if
    end block get_equivalent

    ! Then we just get rid of redundant point in the list
    sort_reducible: block
        !> Integers for the loops
        integer :: j, k, ctr, ctr2

        ctr = 0
        do j = 1, n
            ctr2 = 0
            do k = j, n
                if (k .eq. j) cycle
                if (all(newqp_sort(:, j) .eq. newqp_sort(:, k))) ctr2 = ctr2 + 1
            end do
            ! If ctr2 is 0, it means that this guy is unique by permutation
            if (ctr2 .eq. 0) ctr = ctr + 1
        end do
        allocate (red_triplet(2, ctr))
        ctr = 0
        do j = 1, n
            ctr2 = 0
            do k = j, n
                if (k .eq. j) cycle
                if (all(newqp_sort(:, j) .eq. newqp_sort(:, k))) ctr2 = ctr2 + 1
            end do
            ! If ctr2 is 0, it means that this guy is unique by permutation
            if (ctr2 .eq. 0) then
                ctr = ctr + 1
                red_triplet(1, ctr) = newqp(1, j)
                red_triplet(2, ctr) = newqp(2, j)
            end if
        end do
    end block sort_reducible

    call mem%deallocate(newqp, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(newqp_sort, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine
