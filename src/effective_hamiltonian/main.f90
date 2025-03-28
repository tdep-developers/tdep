#include "precompilerdefinitions"
program effective_hamiltonian
!!{!src/effective_hamiltonian/manual.md!}
use konstanter, only: r8, lo_tol, lo_kb_hartree, lo_bohr_to_A, lo_frequency_Hartree_to_THz, lo_eV_to_Hartree, &
                      lo_exitcode_physical, lo_exitcode_param, lo_temperaturetol, lo_Hartree_to_eV
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: tochar, walltime, lo_stop_gracefully, open_file, lo_progressbar_init, &
                            lo_progressbar, lo_does_file_exist, lo_trueNtimes

use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_mdsim, only: lo_mdsim

use lo_epot, only: lo_energy_differences

use type_crystalstructure, only: lo_crystalstructure
use options, only: lo_opts

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: ss, uc
type(lo_energy_differences) :: pot
type(lo_mdsim) :: sim

type(lo_forceconstant_secondorder) :: fc2
type(lo_forceconstant_thirdorder) :: fc3
type(lo_forceconstant_fourthorder) :: fc4
! POLAR IFCS?

type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
real(r8), dimension(:, :), allocatable :: ebuf

logical :: generate_configs = .false.

init: block

    integer :: f, i, j, l, readrank
    logical :: readonthisrank, mpiparallel
    real(r8) :: t0

    ! Get CLI options
    call opts%parse()

    if (opts%nconf .gt. 0) generate_configs = .true.

    if (mw%talk) then
        write (*, *) 'Recap of the parameters governing the calculation'
        write (*, '(1X,A40,L3)') 'Thirdorder contribution                 ', opts%thirdorder
        write (*, '(1X,A40,L3)') 'Fourthorder contribution                ', opts%fourthorder
        write (*, '(1X,A40,I5)') 'Stride                                  ', opts%stride
        write(*, '(1X,A40,L3)') 'Generate canonical configs               ', generate_configs
        write (*, '(1X,A40,L3)') 'Quantum configurations                  ', opts%quantum
        write (*, '(1X,A40,F20.12)') 'Temperature                         ', opts%temperature
    end if


    call mw%init()
    call mem%init()
    if (.not. mw%talk) opts%verbosity = -100

    ! Read structures
    if (mw%talk) write (*, *) '... reading infiles'
    call ss%readfromfile('infile.ssposcar', verbosity=opts%verbosity)
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)

    ! Match the supercell to the unitcell, always a good idea
    call uc%classify('spacegroup', timereversal=.true.)
    call ss%classify('supercell', uc)

    call fc2%readfromfile(uc, 'infile.forceconstant', mem, verbosity=-1)
    if (opts%thirdorder) call fc3%readfromfile(uc, 'infile.forceconstant_thirdorder')
    if (opts%fourthorder) call fc4%readfromfile(uc, 'infile.forceconstant_fourthorder')
    if (mw%talk) write (*, *) '... read forceconstants'

    if (generate_configs) then
        if ((opts%quantum .eqv. .false.) .and. (opts%temperature .lt. lo_temperaturetol)) then
            call lo_stop_gracefully(['For classical statistics temperature has to be nonzero'], lo_exitcode_physical, __FILE__, __LINE__)
        end if
        if (mw%talk .and. opts%quantum) write (*, *) '... will generate quantum configurations'
        if (mw%talk .and. opts%quantum) write (*, *) '... will generate classical configurations'
    end if

    call pot%setup(uc, ss, fc2, fc3, fc4, mw, opts%verbosity)
    if (mw%talk) write (*, *) '... setup potential energy calculator'

    if (.not. generate_configs) then
        ! If packed simulation is there might as well read it
        if (lo_does_file_exist('infile.sim.hdf5')) then
            call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=opts%stride)
        else 
            ! Just read the infile.positions/ infile.meta dont need forces or stat files
            ! Will need to populate some fields of sim manually

            ! should I just do sim%readfromfile? would let me calculate things like U0 and other cumulants
            ! to get infile.positions you need to run MD anyway and might as well dump the others

            sim%crystalstructure = ss ! pretty sure this copies, but oh well its not a lot of data


            ! We are reading on one rank, then broadcasting
            readrank = 0
            if (mw%r .eq. readrank) then
                readonthisrank = .true.
            else
                readonthisrank = .false.
            end if


            if (readonthisrank) then
                f = open_file('in', 'infile.meta')
                read (f, *) sim%na
                read (f, *) sim%nt
                close (f)
            end if

            lo_allocate(sim%r(3, sim%na, sim%nt))
            sim%r = 0.0_r8
            
            if (readonthisrank) then
                t0 = walltime()
                if (opts%verbosity .gt. 0) call lo_progressbar_init()
                f = open_file('in', 'infile.positions')
                i = 0
                do l = 1, sim%nt

                    i = i + 1
                    do j = 1, sim%na
                        read (f, *) sim%r(:, j, i)
                    end do

                    if (opts%verbosity .gt. 0) then
                    if (lo_trueNtimes(l, 25, sim%nt)) then
                        call lo_progressbar('... reading positions', l, sim%nt, walltime() - t0)
                    end if
                    end if
                end do
                close (f)
            end if

            call mw%bcast(sim%r, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%na, readrank, __FILE__, __LINE__)
            call mw%bcast(sim%nt, readrank, __FILE__, __LINE__)

            ! get non-pbc positions
            call sim%get_nonpbc_positions()
            if (opts%verbosity .gt. 0) write (*, *) '... unwrapped positions'
            ! get displacements
            call sim%get_displacements()
            if (opts%verbosity .gt. 0) write (*, *) '... got displacements'

        end if
        if (mw%talk) write (*, *) '... parsed simulation data'
    end if

end block init




energy : block

    integer :: ctr, i
    real(r8), dimension(:, :), allocatable :: f2, f3, f4, fp
    real(r8) :: e2, e3, e4, ep

    call mem%allocate(ebuf, [opts%nconf, 4], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ebuf = 0.0_r8

    if (generate_configs) then
        call pot%statistical_sampling(uc, ss, fc2, opts%nconf, opts%temperature, ebuf, mw, mem, opts%verbosity)
    else

        ! Dummy space for force
        call mem%allocate(f2, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(f3, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(f4, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(fp, [3, ss%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        f2 = 0.0_r8
        f3 = 0.0_r8
        f4 = 0.0_r8
        fp = 0.0_r8
        
        ! how do I add a progress bar here without race condition?
        ctr = 0
        do i = 1, sim%nt  !! TODO SKIP BASED ON STRIDE??

            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle

            ! Calculate the energy, e2/e3/e4/ep are zeroed inside of this call
            call pot%energies_and_forces(sim%u(:,:,i), e2, e3, e4, ep, f2, f3, f4, fp)

            ebuf(i, 1) = e2
            ebuf(i, 2) = e3
            ebuf(i, 3) = e4
            ebuf(i, 4) = ep
        end do
    end if
    if (mw%talk) write (*, *) '... calculated energies'


end block energy


io : block

    integer :: i, u
    real(r8) :: total_energy

    ! file for the energies
    u = open_file('out', 'outfile.energies')
    write (u, '(A,A)') '# Unit:      ', 'eV/atom'
    write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
    write (u, "(A)") '#  conf    Etotal                  Epolar                &
        &Epair                 Etriplet              Equartet'

    do i = 1, sim%nt
        total_energy = sum(ebuf(i,:))
        write (u, "(1X,I5,5(2X,E20.12))") i, total_energy*lo_Hartree_to_eV, ebuf(i, 4)*lo_Hartree_to_eV, &
                                            ebuf(i, 1)*lo_Hartree_to_eV, ebuf(i, 2)*lo_Hartree_to_eV, &
                                            ebuf(i, 3)*lo_Hartree_to_eV
    end do

    ! close outfile.energies
    write (*, '(A)') ' ... energies writen to `outfile.energies`'
    close (u)

end block io

call mw%destroy()

end program