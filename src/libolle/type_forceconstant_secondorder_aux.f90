submodule(type_forceconstant_secondorder) type_forceconstant_secondorder_aux
use konstanter, only: lo_exitcode_blaslapack
use lo_randomnumbers, only: lo_mersennetwister
use gottochblandat, only: lo_planck, lo_harmonic_oscillator_free_energy, lo_negsqrt
use type_distancetable, only: lo_distancetable
use type_blas_lapack_wrappers, only: lo_dsyevr, lo_dsyevd
use lo_sorting, only: lo_qsort
implicit none
contains

!> calculate the frequencies of the commensurate modes
module subroutine commensurate_modes(fcss, ss, uc, fc, mw)
    !> forceconstant for the supercell
    class(lo_forceconstant_secondorder), intent(inout) :: fcss
    !> reference supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> reference unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> reference unitcell forceconstant
    type(lo_forceconstant_secondorder), intent(in) :: fc
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    real(r8), dimension(:, :, :, :), allocatable :: Dlr, Dna
    real(r8), dimension(:, :), allocatable :: D
    real(r8), dimension(:), allocatable :: dr
    real(r8) :: f0
    integer, dimension(:), allocatable :: ind
    integer :: nb, a1, a2, i, j
    integer :: solrnk

    nb = ss%na*3
    if (.not. allocated(fcss%omega)) allocate (fcss%omega(nb))
    if (.not. allocated(fcss%amplitudes)) allocate (fcss%amplitudes(nb))
    if (.not. allocated(fcss%eigenvectors)) allocate (fcss%eigenvectors(nb, nb))
    fcss%omega = 0.0_r8
    fcss%amplitudes = 0.0_r8
    fcss%eigenvectors = 0.0_r8
    solrnk = mw%n - 1

    ! Idiot sanity test. Should think of some more just in case.
    f0 = real(ss%na, r8)/real(uc%na, r8)
    if (abs(f0 - anint(f0)) .gt. lo_sqtol) then
        call lo_stop_gracefully(['Unit and supercell not commensurate'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
    end if

    if (mw%r .eq. solrnk) then
        ! Get the dynamical matrix
        allocate (D(nb, nb))
        D = 0.0_r8
        do a1 = 1, fcss%na
            do i = 1, fcss%atom(a1)%n
                a2 = fcss%atom(a1)%pair(i)%i2
                f0 = ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
                D((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) = D((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) + fcss%atom(a1)%pair(i)%m*f0
            end do
        end do
        ! add the polar part
        if (fc%polar .and. fc%loto%correctiontype .eq. 3) then
            allocate (Dlr(3, 3, ss%na, ss%na))
            allocate (Dna(3, 3, ss%na, ss%na))
            Dlr = 0.0_r8
            Dna = 0.0_r8
            call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, Dlr, 1E-10_r8)
            call nonanalytical_dynamical_matrix(fcss, ss, [1.0_r8, 0.0_r8, 0.0_r8], Dna)
            do a1 = 1, fcss%na
                do a2 = 1, fcss%na
                    f0 = ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
                    D((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) = D((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) + (Dlr(:, :, a1, a2) + Dna(:, :, a1, a2))*f0
                end do
            end do
            deallocate (Dlr)
            deallocate (Dna)
        end if

        ! Make sure it is symmetric?
        do i = 1, nb
        do j = i + 1, nb
            f0 = (D(i, j) + D(j, i))*0.5_r8
            D(i, j) = f0
            D(j, i) = f0
        end do
        end do

        ! Diagonalize
        call lo_dsyevr(D, fcss%omega, z=fcss%eigenvectors, info=i)
        if (i .ne. 0) then
            call lo_stop_gracefully(['dsyevr exit code'//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__, mw%comm)
        end if

        ! Fix acoustic modes, just set the three smallest to zero
        allocate (dr(nb))
        allocate (ind(nb))
        dr = abs(fcss%omega)
        call lo_qsort(dr, ind)
        do i = 1, 3
            fcss%omega(ind(i)) = 0.0_r8
        end do
        ! Convert to frequencies
        fcss%omega = lo_negsqrt(fcss%omega)

        deallocate (dr)
        deallocate (ind)
        deallocate (D)
    end if

    ! Now distribute over ranks.
    call mw%bcast(fcss%omega, solrnk, __FILE__, __LINE__)
    call mw%bcast(fcss%eigenvectors, solrnk, __FILE__, __LINE__)
end subroutine

!> create a fake forceconstant from a Debye temperature
module subroutine fake_forceconstant(fc, uc, ss, debye_temperature, maximum_frequency, verbosity)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(out) :: fc
    !> unit cell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> a debye temperature to match
    real(r8), intent(in), optional :: debye_temperature
    !> a maximum frequency to match
    real(r8), intent(in), optional :: maximum_frequency
    !> how much to talk
    integer, intent(in), optional :: verbosity

    type(lo_forceconstant_secondorder) :: fcss
    real(r8) :: t0
    integer :: verb, opttype

    ! Figure out what to do
    init: block
        type(lo_distancetable) :: dt
        integer :: a1, i
        t0 = walltime()
        if (present(verbosity)) then
            verb = verbosity
        else
            verb = 0
        end if

        ! figure out what to do
        opttype = 0
        if (present(debye_temperature) .and. present(maximum_frequency)) then
            write (*, *) 'Decide wether you want to match to Debye temperature or max frequency. Not both.'
            stop
        end if
        if (present(debye_temperature)) opttype = 1
        if (present(maximum_frequency)) opttype = 2

        ! say something?
        if (verb .gt. 0) then
            write (*, *) ''
            select case (opttype)
            case (0)
                write (*, *) 'You have to optimize with respect to something.'
                stop
            case (1)
                write (*, *) 'Creating a fake forceconstant, corresponding to a Debye temperature of ', &
                    tochar(debye_temperature)
            case (2)
                write (*, *) 'Creating a fake forceconstant, corresponding to a max frequency of ', &
                    tochar(maximum_frequency*lo_frequency_Hartree_to_THz), 'THz'
            end select
        end if

        ! create an empty forceconstant
        call dt%generate(uc%r, uc%latticevectors, ss%maxcutoff()*1.3_r8, verbosity=0)
        fc%na = uc%na
        fc%cutoff = dt%cutoff
        allocate (fc%atom(fc%na))
        do a1 = 1, fc%na
            fc%atom(a1)%n = dt%particle(a1)%n
            allocate (fc%atom(a1)%pair(fc%atom(a1)%n))
            do i = 1, fc%atom(a1)%n
                fc%atom(a1)%pair(i)%m = 0.0_r8
                fc%atom(a1)%pair(i)%r = dt%particle(a1)%v(:, i)
                fc%atom(a1)%pair(i)%lv1 = dt%particle(a1)%v(:, i)
                fc%atom(a1)%pair(i)%lv2 = dt%particle(a1)%lv(:, i)
                fc%atom(a1)%pair(i)%i1 = a1
                fc%atom(a1)%pair(i)%i2 = dt%particle(a1)%ind(i)
            end do
        end do
        if (verb .gt. 0) write (*, *) '... built skeleton for forceconstant'

        ! Map it to the supercell
        call fc%remap(uc, ss, fcss)
        if (verb .gt. 0) write (*, *) '... mapped unit and supercell forceconstants to each other'
    end block init

    ! Now optimize
    select case (opttype)
    case (1)
        debyefit: block
            ! Fit the phonons to a Debye temperature.
            real(r8) :: target_zpe, zpe, alpha, step
            integer :: i

            ! the target zero-point energy (per atom)
            target_zpe = 9.0_r8*lo_kb_hartree*debye_temperature/8.0_r8

            ! choose something too small for my crappy solver
            alpha = 1.0_r8
            i = 0
            do
                i = i + 1
                call set_values(alpha, fcss)
                zpe = get_zpe(ss, fcss)
                if (zpe .gt. target_zpe) then
                    alpha = alpha*0.5_r8
                else
                    exit
                end if
            end do

            ! and fit it properly. this is reasonably fast.
            step = 2.0_r8
            do i = 1, 2000 ! or something
                call set_values(alpha, fcss)
                zpe = get_zpe(ss, fcss)
                if (zpe .lt. target_zpe) then
                    alpha = alpha*step
                else
                    alpha = alpha/step
                    step = (step + 1.0_r8)*0.5_r8
                end if
                if (abs(zpe - target_zpe) .lt. lo_tol) exit
            end do
            call set_values(alpha, fc)
        end block debyefit
    case (2)
        freqfit: block
            ! Fit the phonons to a maximum temperature
            real(r8) :: target_maxf, alpha, maxom, f0

            ! To Hz
            target_maxf = maximum_frequency
            ! guess something
            alpha = 1.0_r8
            call set_values(alpha, fcss)
            maxom = maxomega(ss, fcss)
            ! omega goes as sqrt(alpha), so adjust alpha
            f0 = (target_maxf/maxom)**2
            alpha = alpha*f0
            call set_values(alpha, fc)
        end block freqfit
    end select

    if (verb .gt. 0) write (*, *) '... got a fake forceconstant in ', tochar(walltime() - t0), 's'
contains

    !> returns the supercell zero point energy
    function maxomega(ss, fcss) result(maxom)
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_forceconstant_secondorder), intent(in) :: fcss
        real(r8) :: maxom

        real(r8), dimension(:, :), allocatable :: dynmat
        real(r8), dimension(:), allocatable :: om
        integer :: a1, a2, i

        allocate (dynmat(ss%na*3, ss%na*3))
        allocate (om(ss%na*3))
        dynmat = 0.0_r8
        om = 0.0_r8
        do a1 = 1, fcss%na
        do i = 1, fcss%atom(a1)%n
            a2 = fcss%atom(a1)%pair(i)%i2
            dynmat((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) = &
                dynmat((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) + fcss%atom(a1)%pair(i)%m*ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
        end do
        end do
        call lo_dsyevd(dynmat, om, 'V')
        maxom = sqrt(abs(maxval(om)))
        deallocate (dynmat)
        deallocate (om)
    end function

    !> returns the supercell zero point energy
    function get_zpe(ss, fcss) result(zpe)
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_forceconstant_secondorder), intent(in) :: fcss
        real(r8) :: zpe

        real(r8), dimension(:, :), allocatable :: dynmat
        real(r8), dimension(:), allocatable :: om
        real(r8) :: f0
        integer :: a1, a2, i

        allocate (dynmat(ss%na*3, ss%na*3))
        allocate (om(ss%na*3))
        dynmat = 0.0_r8
        om = 0.0_r8
        do a1 = 1, fcss%na
        do i = 1, fcss%atom(a1)%n
            a2 = fcss%atom(a1)%pair(i)%i2
            dynmat((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) = &
                dynmat((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) + fcss%atom(a1)%pair(i)%m*ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
        end do
        end do
        call lo_dsyevd(dynmat, om, 'V')
        ! get the zero-point energy
        zpe = 0.0_r8
        do i = 1, ss%na*3
            if (om(i) .lt. 0.0_r8) then
                if (abs(om(i)) .gt. lo_freqtol) then
                    zpe = zpe - 1E50_r8
                end if
            else
                f0 = sqrt(abs(om(i)))
                if (f0 .gt. lo_freqtol) then
                    zpe = zpe + lo_harmonic_oscillator_free_energy(0.0_r8, f0)
                end if
            end if
        end do
        ! energy per atom
        zpe = zpe/ss%na
        deallocate (dynmat)
        deallocate (om)
    end function

    !> return a forceconstant for a specific scaling factor
    subroutine set_values(alpha, fc)
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        real(r8), intent(in) :: alpha
        !
        integer :: a1, i
        real(r8), dimension(3, 3) :: m
        real(r8) :: r, rx, ry, rz, beta, fpp

        ! set values for the forceconstant
        beta = 6.0_r8
        do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            r = norm2(fc%atom(a1)%pair(i)%r)
            if (r .lt. lo_tol) then
                fpp = 0.0_r8
                fc%atom(a1)%pair(i)%m = 0.0_r8
                cycle
            else
                fpp = alpha/(r**beta)
            end if
            rx = fc%atom(a1)%pair(i)%r(1)
            ry = fc%atom(a1)%pair(i)%r(2)
            rz = fc%atom(a1)%pair(i)%r(3)
            m(1, 1) = -rx**2*r*fpp
            m(2, 2) = -ry**2*r*fpp
            m(3, 3) = -rz**2*r*fpp
            m(1, 2) = rx*ry*(-r*fpp)
            m(1, 3) = rx*rz*(-r*fpp)
            m(2, 3) = ry*rz*(-r*fpp)
            m(2, 1) = m(1, 2)
            m(3, 1) = m(1, 3)
            m(3, 2) = m(2, 3)
            m = m/(r**2)
            m = lo_chop(m, 1E-10_r8)
            fc%atom(a1)%pair(i)%m = m
        end do
        end do
        call fc%setsumtozero()

    end subroutine
end subroutine

!> use the harmonic model to initialize a cell
module subroutine initialize_cell(fcss, ss, uc, fc, temperature, quantum, exact, closest_distance, mw, nosync)
    !> force constant for this (super) cell
    class(lo_forceconstant_secondorder), intent(inout) :: fcss
    !> supercell to be thermally populated
    type(lo_crystalstructure), intent(inout) :: ss
    !> reference unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> reference unitcell forceconstant
    type(lo_forceconstant_secondorder), intent(in) :: fc
    !> the temperature
    real(r8), intent(in) :: temperature
    !> classical or quantum statistics?
    logical, intent(in) :: quantum
    !> set the temperature exactly to what we want?
    logical, intent(in) :: exact
    !> make sure atoms are not too close to each another
    real(r8), intent(in) :: closest_distance
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> No syncing, different distribution on all ranks
    logical, intent(in), optional :: nosync

    ! Not sure about save attribute here.
    type(lo_mersennetwister), save :: tw
    integer :: solrnk
    logical :: sync

    init: block
        ! Seed rng if needed
        if (tw%initialized .eqv. .false.) then
            ! make sure to seed it with the rank to make it different on all MPI ranks.
            call tw%init(iseed=mw%r, rseed=walltime())
        end if

        ! Decide on syncinc thingies
        if (present(nosync)) then
            sync = .not. nosync
        else
            sync = .true.
        end if

        ! Which rank to solve on
        if (sync) then
            ! Only solve on one rank
            solrnk = mw%n - 1
        else
            ! Solve on all ranks
            solrnk = mw%r
        end if

        ! If it's the first call, we need to calculate the harmonic stuff
        if (allocated(fcss%omega) .eqv. .false.) then
            call fcss%commensurate_modes(ss, uc, fc, mw)
        end if
    end block init

    if (solrnk .eq. mw%r) then
        ! Set the amplitudes
        setamplitude: block
            integer :: i
            do i = 1, fcss%na*3
                if (fcss%omega(i) .gt. lo_freqtol) then
                    ! Choose quantum or classical statistics
                    if (quantum) then
                        fcss%amplitudes(i) = sqrt((2*lo_planck(temperature, fcss%omega(i)) + 1)*0.5_r8/fcss%omega(i))
                    else
                        fcss%amplitudes(i) = sqrt(lo_kb_hartree*temperature)/fcss%omega(i)
                    end if
                else
                    ! set to zero for acoustic modes
                    fcss%amplitudes(i) = 0.0_r8
                end if
            end do
        end block setamplitude

        builddisplacements: block
            integer, parameter :: maxiter = 100
            type(lo_distancetable) :: dt
            real(r8), dimension(3) :: ctot, ptot, v0
            real(r8) :: x1, x2, f0, f1, nndist
            integer :: a1, i, j, l, iter

            ! Build the displacements and stuff
            ss%u = 0.0_r8
            ss%v = 0.0_r8
            if (closest_distance .gt. 0.0_r8) then
                ! Build a distance table to be used for distance checking
                nndist = ss%nearest_neighbour_distance()
                call dt%generate(ss%r, ss%latticevectors, nndist*2.0_r8, verbosity=0)
            end if

            ! Iterate to make sure I get a configuration with no atoms too close to each other?
            dstloop: do iter = 1, maxiter
                ! Generate displacements
                modeloop: do i = 1, ss%na*3
                    call tw%rnd_boxmuller_pair(1.0_r8, 0.0_r8, x1, x2)
                    l = 0
                    do a1 = 1, ss%na
                    do j = 1, 3
                        l = l + 1
                        ss%u(j, a1) = ss%u(j, a1) + fcss%amplitudes(i)*ss%invsqrtmass(a1)*x1*fcss%eigenvectors(l, i)
                        ss%v(j, a1) = ss%v(j, a1) - fcss%amplitudes(i)*ss%invsqrtmass(a1)*x2*fcss%eigenvectors(l, i)*fcss%omega(i)
                    end do
                    end do
                end do modeloop

                if (closest_distance .gt. 0.0_r8) then
                    ! Check distances
                    f0 = lo_huge
                    do i = 1, dt%np
                    do j = 2, dt%particle(i)%n
                        v0 = dt%particle(i)%v(:, j) + ss%u(:, dt%particle(i)%ind(j)) - ss%u(:, i)
                        f0 = min(f0, norm2(v0))
                    end do
                    end do
                    f0 = f0/nndist
                    if (f0 .lt. closest_distance) then
                        ! It's ok, configuration is fine
                        exit dstloop
                    end if
                else
                    ! No check, just skip iterating
                    exit dstloop
                end if

                ! If we went through too many times and it did not help, I will manually
                ! force the displacements small enough.
                if (iter .eq. maxiter) then
                    do i = 1, ss%na
                        f0 = norm2(ss%u(:, i))
                        if (f0/nndist .gt. closest_distance*0.5_r8) then
                            ss%u(:, i) = ss%u(:, i)*closest_distance*0.5_r8/f0
                        end if
                    end do
                end if
            end do dstloop

            ! Remove all drift, make sure that momentum add up to zero, and center of mass does not move.
            ctot = 0.0_r8
            ptot = 0.0_r8
            do i = 1, ss%na
                ctot = ctot + ss%u(:, i)*ss%mass(i)
                ptot = ptot + ss%v(:, i)*ss%mass(i)
            end do
            ctot = ctot/sum(ss%mass)
            ptot = ptot/sum(ss%mass)
            do i = 1, ss%na
                ss%u(:, i) = ss%u(:, i) - ctot
                ss%v(:, i) = ss%v(:, i) - ptot
            end do

            ! Maybe we want to adjust it exactly to the desired temperature
            if (exact) then
                f0 = fcss%potential_energy(ss%u)
                f1 = ss%kinetic_energy()
                f0 = f0/(ss%na - 1)
                f1 = f1/(ss%na - 1)
                ! temperature now
                f0 = (f0 + f1)/(3*lo_kb_hartree)
                ! factor that fixes it
                f0 = sqrt(temperature/f0)
                ss%v = ss%v*f0
                ss%u = ss%u*f0
            end if
            ! Add the displacements to the positions
            do a1 = 1, ss%na
                ss%r(:, a1) = lo_clean_fractional_coordinates(ss%r(:, a1) + matmul(ss%inv_latticevectors, ss%u(:, a1)))
                ss%rcart(:, a1) = ss%fractional_to_cartesian(ss%r(:, a1))
            end do
        end block builddisplacements
    end if

    ! Then sync to all ranks.
    if (sync) then
        call mw%bcast(ss%r, solrnk, __FILE__, __LINE__)
        call mw%bcast(ss%rcart, solrnk, __FILE__, __LINE__)
        call mw%bcast(ss%u, solrnk, __FILE__, __LINE__)
        call mw%bcast(ss%v, solrnk, __FILE__, __LINE__)
    end if
end subroutine

!> Ensure translational invariance is obeyed
module subroutine setsumtozero(fc)
    !> the force constant
    class(lo_forceconstant_secondorder), intent(inout) :: fc

    real(r8), dimension(3, 3) :: m
    integer :: a1, vec, l

    do a1 = 1, fc%na
        ! find the force constant corresponding to the self term
        l = 0
        do vec = 1, fc%atom(a1)%n
            if (lo_sqnorm(fc%atom(a1)%pair(vec)%r) .lt. lo_sqtol) then
                l = vec
                exit
            end if
        end do
        ! sanity check
        if (l .eq. 0) then
            call lo_stop_gracefully(['Could not find the self-term when enforcing secondorder translational invariance'], &
                                    lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
        ! add all the ifcs together
        m = 0.0_r8
        do vec = 1, fc%atom(a1)%n
            if (vec .ne. l) then
                m = m + fc%atom(a1)%pair(vec)%m(:, :)
            end if
        end do
        ! make sure sum is zero
        fc%atom(a1)%pair(l)%m = lo_chop(-m, lo_sqtol)
    end do
end subroutine

end submodule
