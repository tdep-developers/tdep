submodule(ifc_solvers) ifc_solvers_trainpred
!!
!! Helper that selects slices from the full simulation for crossvalidation.
!!
implicit none

contains

!> Divide the simulation across ranks and sort out division
module subroutine divide_simulation(tp, sim, wts, mw)
    !> division helper
    type(lo_trainpred_helper), intent(out) :: tp
    !> MD simulation
    type(lo_mdsim), intent(in) :: sim
    !> Weight per timestep
    real(r8), dimension(:), intent(in) :: wts
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    integer :: i, j, k, l

    ! Count steps on this rank
    tp%nt = 0
    do i = 1, sim%nt
        if (mod(i, mw%n) .eq. mw%r) tp%nt = tp%nt + 1
    end do

    ! Map out the division thingies in a logical way, and count how
    ! many steps I have for prediction and training.
    if (tp%nt .gt. 0) then
        ! This is the case where this rank got some timesteps
        allocate (tp%timestep_ind(tp%nt))
        allocate (tp%weight(tp%nt))
        allocate (tp%train_ts(tp%nt, n_division))
        allocate (tp%pred_ts(tp%nt, n_division))
        tp%timestep_ind = 0
        tp%weight = 0.0_r8
        tp%train_ts = .false.
        tp%pred_ts = .false.
        l = 0
        do i = 1, sim%nt
            if (mod(i, mw%n) .ne. mw%r) cycle
            l = l + 1
            ! make a note which timestep
            tp%timestep_ind(l) = i
            ! store the weight
            tp%weight(l) = sqrt(wts(i))
            j = mod(i, n_division) + 1
            do k = 1, n_division
                if (j .eq. k) then
                    tp%pred_ts(l, k) = .true.
                else
                    tp%train_ts(l, k) = .true.
                end if
            end do
        end do
        ! Now there might be an idiot case where I have no steps for training or
        ! prediction in any rank. Things will get wonky, but it should not crash.
        if (sim%nt .lt. 3) then
            l = 0
            do i = 1, sim%nt
                if (mod(i, mw%n) .ne. mw%r) cycle
                l = l + 1
                tp%pred_ts(l, :) = .true.
                tp%train_ts(l, :) = .true.
            end do
        end if
        ! Now count
        tp%nstep_pred = 0
        tp%nstep_train = 0
        l = 0
        do i = 1, sim%nt
            if (mod(i, mw%n) .ne. mw%r) cycle
            l = l + 1
            do k = 1, n_division
                if (tp%pred_ts(l, k)) tp%nstep_pred(k) = tp%nstep_pred(k) + 1
                if (tp%train_ts(l, k)) tp%nstep_train(k) = tp%nstep_train(k) + 1
            end do
        end do
    else
        ! This rank got no timesteps. It happens.
        allocate (tp%train_ts(1, n_division))
        allocate (tp%pred_ts(1, n_division))
        allocate (tp%timestep_ind(1))
        allocate (tp%weight(1))
        tp%timestep_ind = -lo_hugeint
        tp%weight = -lo_huge
        tp%train_ts = .false.
        tp%pred_ts = .false.
        tp%nstep_pred = 0
        tp%nstep_train = 0
    end if
end subroutine

!> collect subsets of data for either training or test of predictions
module subroutine collect_training_prediction( &
    tp, map, interaction, div, training, nrow, mem, dh, ih, coeff, values, displacements, energies)
    !> division helper
    class(lo_trainpred_helper), intent(in) :: tp
    !> map over interactions
    type(lo_forcemap), intent(in) :: map
    !> which interaction to collect
    integer, intent(in) :: interaction
    !> which division are we working on
    integer, intent(in) :: div
    !> are we collecting training data or response data
    logical, intent(in) :: training
    !> how many rows did we get (on this rank)
    integer, intent(out) :: nrow
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> coefficient matrices for born charges
    type(lo_dielectric_helper), intent(in), optional :: dh
    !> coefficient matrices for force constants
    type(lo_ifc_helper), intent(in), optional :: ih
    !> coefficient matrix?
    real(r8), dimension(:, :), allocatable, intent(out), optional :: coeff
    !> values to predict?
    real(r8), dimension(:), allocatable, intent(out), optional :: values
    !> displacements corresponding to values
    real(r8), dimension(:), allocatable, intent(out), optional :: displacements
    !> energies for certain configurations
    real(r8), dimension(:), allocatable, intent(out), optional :: energies

    integer :: nf, nx, ne, nt
    integer :: i, ctr, i1, i2, j1, j2

    ! Which interaction are we working with?
    select case (interaction)
    case (-3)
        ! fetch raw Z
        nf = 9*map%n_atom_ss
        nx = -1
    case (-2)
        ! fetch raw eps
        nf = 9
        nx = -1
    case (-1)
        ! fetch polar data
        nf = 3*map%n_atom_ss
        nx = -1
    case (0)
        ! fetch raw force/displacement/energy data
        nf = 3*map%n_atom_ss
        nx = -1
    case (1)
        ! first order force constants
        nf = 3*map%n_atom_ss
        nx = map%xuc%nx_fc_singlet
    case (2)
        ! second order force constants
        nf = 3*map%n_atom_ss
        nx = map%xuc%nx_fc_pair
    case (3)
        ! third order force constants
        nf = 3*map%n_atom_ss
        nx = map%xuc%nx_fc_triplet
    case (4)
        ! third order force constants
        nf = 3*map%n_atom_ss
        nx = map%xuc%nx_fc_quartet
    case (10)
        ! global dielectric constant
        nf = 9
        nx = map%xuc%nx_eps_global
    case (11)
        ! dielectric singlet
        nf = 9
        nx = map%xuc%nx_eps_singlet
    case (12)
        ! dielectric pair
        nf = 9
        nx = map%xuc%nx_eps_pair
    case (21)
        ! Born charge singlets
        nf = 9*map%n_atom_ss
        nx = map%xuc%nx_Z_singlet
    case (22)
        ! Born charge pairs
        nf = 9*map%n_atom_ss
        nx = map%xuc%nx_Z_pair
    case (23)
        ! Born charge triplets
        nf = 9*map%n_atom_ss
        nx = map%xuc%nx_Z_triplet
    case default
        call lo_stop_gracefully(['Unknown interaction'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    ! How many equations/timesteps?
    if (training) then
        ne = tp%nstep_train(div)*nf
        nt = tp%nstep_train(div)
    else
        ne = tp%nstep_pred(div)*nf
        nt = tp%nstep_pred(div)
    end if

    if (nt .eq. 0) then
        ! Nothing on this rank
        nrow = 0
        if (present(coeff)) then
            call mem%allocate(coeff, [1, 1], persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            coeff = 0.0_r8
        end if
        if (present(values)) then
            call mem%allocate(values, 1, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            values = 0.0_r8
        end if
        if (present(displacements)) then
            call mem%allocate(displacements, 1, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            displacements = 0.0_r8
        end if
        if (present(energies)) then
            call mem%allocate(energies, 1, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            energies = 0.0_r8
        end if
        ! and exit early
        return
    end if

    ! space for the things to return
    if (present(coeff)) then
        call mem%allocate(coeff, [ne, nx], persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
        coeff = 0.0_r8
    end if
    if (present(values)) then
        call mem%allocate(values, ne, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
        values = 0.0_r8
    end if
    if (present(displacements)) then
        call mem%allocate(displacements, ne, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
        displacements = 0.0_r8
    end if
    if (present(energies)) then
        call mem%allocate(energies, nt, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
        energies = 0.0_r8
    end if

    ! collect values
    ctr = 0
    do i = 1, tp%nt
        ! Figure out if I should skip the step
        if (training) then
            if (tp%pred_ts(i, div)) cycle
        else
            if (tp%train_ts(i, div)) cycle
        end if
        ! get the ranges
        ctr = ctr + 1
        i1 = (ctr - 1)*nf + 1
        j1 = ctr*nf
        i2 = (i - 1)*nf + 1
        j2 = i*nf
        select case (interaction)
        case (-3)
            ! collect raw Z
            if (present(values)) values(i1:j1) = dh%Z0(i2:j2)
        case (-2)
            ! collect raw eps
            if (present(values)) values(i1:j1) = dh%epso(i2:j2)
        case (-1)
            ! collect polar energies and forces
            if (present(energies)) energies(ctr) = ih%epol(i)
            if (present(values)) values(i1:j1) = ih%fp(i2:j2)
        case (0)
            ! collect raw values
            if (present(energies)) energies(ctr) = ih%epot(i)
            if (present(displacements)) displacements(i1:j1) = ih%u0(i2:j2)
            if (present(values)) values(i1:j1) = ih%f0(i2:j2) + ih%fp(i2:j2)
        case (1)
            if (present(coeff)) coeff(i1:j1, :) = ih%cm1(i2:j2, :)
            if (present(values)) values(i1:j1) = ih%f1(i2:j2)
        case (2)
            if (present(coeff)) coeff(i1:j1, :) = ih%cm2(i2:j2, :)
            if (present(values)) values(i1:j1) = ih%f2(i2:j2)
        case (3)
            if (present(coeff)) coeff(i1:j1, :) = ih%cm3(i2:j2, :)
            if (present(values)) values(i1:j1) = ih%f3(i2:j2)
        case (4)
            if (present(coeff)) coeff(i1:j1, :) = ih%cm4(i2:j2, :)
            if (present(values)) values(i1:j1) = ih%f4(i2:j2)
        case (10)
            if (present(coeff)) coeff(i1:j1, :) = dh%cm_eps0(i2:j2, :)
            if (present(values)) values(i1:j1) = dh%eps0(i2:j2)
        case (11)
            if (present(coeff)) coeff(i1:j1, :) = dh%cm_eps1(i2:j2, :)
            if (present(values)) values(i1:j1) = dh%eps1(i2:j2)
        case (12)
            if (present(coeff)) coeff(i1:j1, :) = dh%cm_eps2(i2:j2, :)
            if (present(values)) values(i1:j1) = dh%eps2(i2:j2)
        case (21)
            if (present(coeff)) coeff(i1:j1, :) = dh%cm_Z1(i2:j2, :)
            if (present(values)) values(i1:j1) = dh%Z1(i2:j2)
        case (22)
            if (present(coeff)) coeff(i1:j1, :) = dh%cm_Z2(i2:j2, :)
            if (present(values)) values(i1:j1) = dh%Z2(i2:j2)
        case (23)
            if (present(coeff)) coeff(i1:j1, :) = dh%cm_Z3(i2:j2, :)
            if (present(values)) values(i1:j1) = dh%Z3(i2:j2)
        end select
    end do

    ! And return the number of rows we got
    nrow = ne
end subroutine

end submodule
