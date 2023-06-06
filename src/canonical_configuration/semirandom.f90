module semirandom
!!
!! Generate configurations with decent condition number
!!
use konstanter, only: r8, lo_iou, lo_huge, lo_hugeint, lo_exitcode_param
use gottochblandat, only: walltime, tochar, lo_progressbar, lo_progressbar_init, lo_real_singular_value_decomposition
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use lo_randomnumbers, only: lo_mersennetwister
use type_forcemap, only: lo_forcemap, lo_coeffmatrix_eps_singlet, lo_coeffmatrix_eps_pair
!use type_mdsim, only: lo_mdsim
implicit none
private
public :: generate_semirandom_configurations

contains

!> create configurations that are not quite random but still quite random
subroutine generate_semirandom_configurations(uc, ss, fc, fcss, temperature, quantum, dc2, dc3, output_format, nconf, mw, mem, verbosity)
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss
    !> unitcell forceconstant
    type(lo_forceconstant_secondorder), intent(in) :: fc
    !> supercell forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fcss
    !> temperature
    real(r8), intent(in) :: temperature
    !> quantum statistics
    logical, intent(in) :: quantum
    !> second order dielectric cutoff
    real(r8), intent(in) :: dc2
    !> third order dielectric cutoff
    real(r8), intent(in) :: dc3
    !> which format to write in
    integer, intent(in) :: output_format
    !> number of configurations to generate, approx
    integer, intent(in) :: nconf
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! need a symmetry list thing
    type(lo_forcemap) :: map
    !type(lo_mdsim), dimension(:), allocatable :: set
    type(lo_mersennetwister) :: tw
    real(r8) :: t0, t1
    integer :: nconf_per_set, nx
    integer :: nset, iset

    ! start timers
    t0 = walltime()
    t1 = t0

    ! Set some basics
    init: block
        type(lo_interaction_tensors) :: sl

        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'GENERATING SEMIRANDOM CONFIGURATIONS'
        end if

        ! Seed random numbers?
        call tw%init(iseed=mw%r, rseed=walltime())

        ! Work out the symmetry?
        call sl%generate(uc, ss, dc2, -1.0_r8, -1.0_r8, .true., mw, mem, verbosity + 1, &
                         dielcutoff2=dc2, dielcutoff3=dc3)
        call map%generate(uc, ss, polarcorrectiontype=3, st=sl, mw=mw, mem=mem, verbosity=verbosity)

        ! So, the purpose of this routine is to get decent configurations that can handle
        ! extraction of polarizability tensor to high orders. If you are unlucky it can be a very
        ! ill-conditioned problem, and I'm trying to fix that. First thing to check is how many
        ! independent we are looking for?
        nx = 0
        if (map%xuc%nx_eps_singlet .gt. 0) nx = nx + map%xuc%nx_eps_singlet
        if (map%xuc%nx_eps_pair .gt. 0) nx = nx + map%xuc%nx_eps_pair
        ! check that the input is not too weird.
        if (nx .eq. 0) then
            call lo_stop_gracefully(['No dielectric interactions, this is strange'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
        ! My logic is that we meaure 6 values per configuration, and we have nx configurations. I want
        ! to overdetermine by a tiny bit and as such choose
        !nconf_per_set=ceiling( 2.0_r8*nx/6.0_r8 )*50
        ! to be the number of configurations that give a solveable set of configurations. Then I generate many of these.
        !nset=ceiling( real(nconf,r8)/nconf_per_set )
        nconf_per_set = nconf
        nset = 1

        if (verbosity .gt. 0) then
            write (lo_iou, *) '  number of unknown to solve for:', nx
            write (lo_iou, *) '          configurations per set:', nconf_per_set
            write (lo_iou, *) '                  number of sets:', nset
        end if

    end block init

    do iset = 1, nset
        buildset: block
            integer, parameter :: nouter = 200     ! picked at random, kind of
            integer, parameter :: ninner = 20      ! picked at random, kind of
            type(lo_crystalstructure) :: p
            real(r8), dimension(:, :, :), allocatable :: old_UM, old_VM, old_RM
            real(r8), dimension(:, :, :), allocatable :: new_UM, new_VM, new_RM
            real(r8), dimension(:, :, :), allocatable :: best_UM, best_VM, best_RM
            real(r8), dimension(:, :), allocatable :: best_CM, old_CM, new_CM
            real(r8), dimension(:, :), allocatable :: CMA1, CMA2
            real(r8), dimension(:), allocatable :: dr
            real(r8) :: mctemp, f0, f1
            real(r8) :: old_cnd, new_cnd, best_cnd
            integer :: i, j, k, l, ii, jj
            integer :: ilo, ihi, jlo, jhi
            integer :: outiter, initer, ctr_accept, ctr_best, ctr_slowiter
            logical :: accept
            character(len=200) :: fname

            ! Some temporary space
            call mem%allocate(best_UM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(best_VM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(best_RM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(old_UM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(old_VM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(old_RM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(new_UM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(new_VM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(new_RM, [3, ss%na, nconf_per_set], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            best_UM = 0.0_r8
            best_VM = 0.0_r8
            best_RM = 0.0_r8
            old_UM = 0.0_r8
            old_VM = 0.0_r8
            old_RM = 0.0_r8
            new_UM = 0.0_r8
            new_VM = 0.0_r8
            new_RM = 0.0_r8

            call mem%allocate(dr, mw%n, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dr = 0.0_r8

            if (map%xuc%nx_eps_singlet .gt. 0) then
                call mem%allocate(CMA1, [9, map%xuc%nx_eps_singlet], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                CMA1 = 0.0_r8
                ilo = 1
                ihi = map%xuc%nx_eps_singlet
            else
                ilo = -1
                ihi = -1
            end if
            if (map%xuc%nx_eps_pair .gt. 0) then
                call mem%allocate(CMA2, [9, map%xuc%nx_eps_pair], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                CMA2 = 0.0_r8
                jlo = 1
                jhi = map%xuc%nx_eps_pair
                if (map%have_eps_singlet) then
                    jlo = jlo + map%xuc%nx_eps_singlet
                    jhi = jhi + map%xuc%nx_eps_singlet
                end if
            else
                jlo = -1
                jhi = -1
            end if

            ! The large coefficient matrix
            call mem%allocate(old_CM, [9*nconf_per_set, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(new_CM, [9*nconf_per_set, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(best_CM, [9*nconf_per_set, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            old_CM = 0.0_r8
            new_CM = 0.0_r8
            best_CM = 0.0_r8

            ! Copy of structure to work with
            p = ss

            ! Generate a starting set of configurations
            t0 = walltime()
            if (verbosity .gt. 0) call lo_progressbar_init()
            do i = 1, nconf_per_set
                p%r = ss%r
                p%rcart = ss%rcart
                p%v = 0.0_r8
                p%u = 0.0_r8
                call fcss%initialize_cell(p, uc, fc, temperature, quantum, .false., -1.0_r8, mw, nosync=.true.)
                old_UM(:, :, i) = p%u
                old_VM(:, :, i) = p%v
                old_RM(:, :, i) = p%r
                ! while we are at it get the coefficient matrices
                if (map%have_eps_singlet) then
                    call lo_coeffmatrix_eps_singlet(map, p%u, CMA1)
                    ii = (i - 1)*9 + 1
                    jj = i*9
                    old_CM(ii:jj, ilo:ihi) = CMA1
                end if
                if (map%have_eps_pair) then
                    call lo_coeffmatrix_eps_pair(map, p%u, CMA2)
                    ii = (i - 1)*9 + 1
                    jj = i*9
                    old_CM(ii:jj, jlo:jhi) = CMA2
                end if

                if (verbosity .gt. 0) then
                    t1 = walltime()
                    call lo_progressbar(' ... generating initial guess', i, nconf_per_set, t1 - t0)
                end if
            end do
            t0 = t1

            ! Starting condition number
            old_cnd = condition_number(old_CM)
            mctemp = 1E-3_r8 ! Or whatever?
            best_cnd = 0.0_r8
            ctr_slowiter = 0

            outloop: do outiter = 1, nouter
                ctr_accept = 0
                ctr_best = 0
                f0 = 1.0_r8
                inloop: do initer = 1, ninner

                    ! reset
                    new_CM = old_CM
                    new_UM = old_UM
                    new_VM = old_VM
                    new_RM = old_RM
                    ! How many configs do I want to change?
                    if (outiter .lt. nouter/10) then
                        k = 1 !tw%rnd_int(nconf_per_set/2)
                    else
                        k = 1
                    end if
                    ! change some configurations
                    confloop: do j = 1, k
                        ! pick a configuration to change
                        i = tw%rnd_int(nconf_per_set)
                        ! get a new config
                        p%r = ss%r
                        p%rcart = ss%rcart
                        p%v = 0.0_r8
                        p%u = 0.0_r8
                        call fcss%initialize_cell(p, uc, fc, temperature, quantum, .false., -1.0_r8, mw, nosync=.true.)
                        new_UM(:, :, i) = p%u
                        new_VM(:, :, i) = p%v
                        new_RM(:, :, i) = p%r

                        if (map%have_eps_singlet) then
                            call lo_coeffmatrix_eps_singlet(map, p%u, CMA1)
                            ii = (i - 1)*9 + 1
                            jj = i*9
                            new_CM(ii:jj, ilo:ihi) = CMA1
                        end if
                        if (map%have_eps_pair) then
                            call lo_coeffmatrix_eps_pair(map, p%u, CMA2)
                            ii = (i - 1)*9 + 1
                            jj = i*9
                            new_CM(ii:jj, jlo:jhi) = CMA2
                        end if
                    end do confloop
                    ! new condition number
                    new_cnd = condition_number(new_CM)

                    ! do something about this
                    if (new_cnd .gt. old_cnd) then
                        accept = .true.
                    else
                        f0 = exp((new_cnd - old_cnd)/mctemp)
                        f1 = tw%rnd_real()
                        if (f0 .gt. f1) then
                            accept = .true.
                        else
                            accept = .false.
                        end if
                    end if

                    ! Keep the best?
                    if (new_cnd .gt. best_cnd) then
                        ctr_best = ctr_best + 1
                        best_cnd = new_cnd
                        best_UM = new_UM
                        best_VM = new_VM
                        best_RM = new_RM
                        best_CM = new_CM
                    end if

                    if (accept) then
                        old_cnd = new_cnd
                        old_UM = new_UM
                        old_VM = new_VM
                        old_RM = new_RM
                        old_CM = new_CM
                        ctr_accept = ctr_accept + 1
                    end if

                    !if ( mw%talk ) write(*,*) outiter,initer,old_cnd,new_cnd,f0
                end do inloop

                if (mw%talk) then
                    write (*, *) outiter, new_cnd, best_cnd, mctemp, ctr_accept
                end if

                ! maybe lower temperature?
                if (ctr_accept .gt. ninner/2) then
                    mctemp = mctemp*0.25_r8
                elseif (ctr_accept .le. 1) then
                    mctemp = mctemp*1.82_r8
                end if

                ! compare results across ranks and continue from the best one
                dr = 0.0_r8
                dr(mw%r + 1) = best_cnd
                call mw%allreduce('max', dr)
                f0 = -lo_huge
                j = 0
                do i = 1, mw%n
                    if (dr(i) .gt. f0) then
                        f0 = dr(i)
                        j = i
                    end if
                end do

                ! make sure everyone agrees what is best.
                call mw%bcast(best_UM, from=j - 1)
                call mw%bcast(best_VM, from=j - 1)
                call mw%bcast(best_RM, from=j - 1)
                call mw%bcast(best_CM, from=j - 1)
                call mw%bcast(best_cnd, from=j - 1)

                ! Replace current with the best
                old_UM = best_UM
                old_VM = best_VM
                old_RM = best_RM
                old_CM = best_CM
                old_cnd = best_cnd

                ! check how we are doing
                call mw%allreduce('max', ctr_best)
                if (ctr_best .eq. 0) then
                    ctr_slowiter = ctr_slowiter + 1
                else
                    ctr_slowiter = 0
                end if

                if (ctr_slowiter .ge. 3) then
                    ! we have reached something
                    exit outloop
                end if
            end do outloop

            ! Now write these configurations before moving on?
            if (mw%talk) then
                do i = 1, nconf_per_set
                    j = (iset - 1)*nconf_per_set + i
                    p = ss
                    p%r = best_RM(:, :, i)
                    p%u = best_UM(:, :, i)
                    p%v = best_VM(:, :, i)
                    select case (output_format)
                    case (1) ! vasp output
                        fname = 'contcar_conf'//tochar(j, 4)
                        call p%writetofile(trim(fname), output_format, write_velocities=.true.)
                    case (2) ! abinit output
                        fname = 'abinput_conf'//tochar(j, 4)
                        call p%writetofile(trim(fname), output_format, write_velocities=.true.)
                    case (3) ! LAMMPS output
                        fname = 'lammps_conf'//tochar(j, 4)
                        call p%writetofile(trim(fname), output_format, write_velocities=.true.)
                    case (4) ! AIMS output
                        fname = 'aims_conf'//tochar(j, 4)
                        call p%writetofile(trim(fname), output_format, write_velocities=.true.)
                    case (5) ! Siesta output
                        fname = 'siesta_conf'//tochar(j, 4)
                        call p%writetofile(trim(fname), output_format, write_velocities=.true.)
                    end select
                end do
            end if

        end block buildset
    end do

    call mw%destroy()
    stop

end subroutine

!> get the condition number of a matrx
function condition_number(A) result(cnd)
    !> matrix
    real(r8), dimension(:, :), intent(in) :: A
    !> condition number
    real(r8) :: cnd

    real(r8), dimension(:), allocatable :: S

    call lo_real_singular_value_decomposition(A, S)
    cnd = minval(S)/maxval(S)
    deallocate (S)
end function

end module
