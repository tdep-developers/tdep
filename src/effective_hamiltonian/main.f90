#include "precompilerdefinitions"
program effective_hamiltonian
!!{!src/effective_hamiltonian/manual.md!}
use konstanter, only: r8, lo_tol, lo_kb_hartree, lo_bohr_to_A, lo_frequency_Hartree_to_THz, lo_eV_to_Hartree, &
                      lo_exitcode_physical, lo_exitcode_param, lo_temperaturetol, lo_Hartree_to_eV
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: tochar, walltime, lo_stop_gracefully, open_file, lo_progressbar_init, &
                            lo_progressbar, lo_does_file_exist, lo_trueNtimes, lo_mean

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
real(r8), dimension(:, :), allocatable :: pebuf
! real(r8), dimension(:), allocatable :: kebuf

logical :: generate_configs = .false.


call opts%parse()
call mw%init()
call mem%init()

init: block

    integer :: f, i, j, l, readrank
    logical :: readonthisrank, mpiparallel
    real(r8) :: t0


    if (.not. mw%talk) opts%verbosity = -100

    if (opts%nconf .gt. 0) generate_configs = .true.

    if (mw%talk) then
        write (*, *) 'Recap of the parameters governing the calculation'
        write (*, '(1X,A40,L3)') 'Thirdorder contribution                 ', opts%thirdorder
        write (*, '(1X,A40,L3)') 'Fourthorder contribution                ', opts%fourthorder
        write (*, '(1X,A40,I5)') 'Stride                                  ', opts%stride
        write(*, '(1X,A40,L3)') 'Generate canonical configs               ', generate_configs
        if(generate_configs) then
            write(*, '(1X,A40,I5)') 'Num Configs                          ', opts%nconf
            write (*, '(1X,A40,L3)') 'Quantum configurations              ', opts%quantum
            write (*, '(1X,A40,F20.12)') 'Temperature                     ', opts%temperature
        end if
    end if

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
    
        if (opts%temperature .lt. lo_temperaturetol) then
            call lo_stop_gracefully(['If nconf is passed, the temperature flag must be set'], lo_exitcode_param, __FILE__, __LINE__)
        end if

    end if

    call pot%setup(uc, ss, fc2, fc3, fc4, mw, opts%verbosity + 1)
    if (mw%talk) write (*, *) '... setup potential energy calculator'

    if (.not. generate_configs) then
        ! If packed simulation is there might as well read it
        if (lo_does_file_exist('infile.sim.hdf5')) then
            call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=opts%stride)
        else 
            call sim%read_from_file(verbosity = opts%verbosity + 2, stride = opts%stride, &
                                     magnetic=.false., dielectric=.false., nrand=-1, mw=mw)
        end if
        if (mw%talk) write (*, *) '... parsed simulation data'
    end if

end block init




energy : block

    integer :: i, u 
    real(r8), dimension(:, :), allocatable :: f2, f3, f4, fp
    real(r8) :: e2, e3, e4, ep, total_energy, to_ev_per_atom


    if (generate_configs) then

        if (mw%talk) write (*, *) '... generating canonical configurations'
        call mem%allocate(pebuf, [opts%nconf, 4], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        pebuf = 0.0_r8

        call pot%statistical_sampling(uc, ss, fc2, opts%nconf, opts%temperature, opts%quantum, pebuf, mw, mem, opts%verbosity)

    else

        call mem%allocate(pebuf, [sim%nt, 4], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        pebuf = 0.0_r8

        ! call mem%allocate(kebuf, [sim%nt], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ! kebuf = 0.0_r8

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
        do i = 1, sim%nt 

            if (mod(i, mw%n) .ne. mw%r) cycle

            ! Calculate the energy, e2/e3/e4/ep are zeroed inside of this call
            call pot%energies_and_forces(sim%u(:,:,i), e2, e3, e4, ep, f2, f3, f4, fp)

            pebuf(i, 1) = e2
            pebuf(i, 2) = e3
            pebuf(i, 3) = e4
            pebuf(i, 4) = ep
        end do

        call mw%allreduce('sum', pebuf)

    end if

    if (mw%talk) write (*, *) '... calculated energies'

end block energy

! Modified from anharmonic_free_energy
epotthings: block
    real(r8), dimension(3, 5) :: cumulant
    real(r8), dimension(:, :), allocatable :: ediff
    real(r8) :: inverse_kbt, T_actual, U0, total_energy, hartree_to_mev, to_mev_per_atom
    integer :: i, u

    to_mev_per_atom = 1000*lo_Hartree_to_eV / real(ss%na, r8)

    if (generate_configs) then
        T_actual = opts%temperature
    else
        T_actual = sim%temperature_thermostat
    end if

    ! Don't have sim%stat%potential energy when generating canonical_configs
    ! So we cannot get U0 or cumulants
    if (.not. generate_configs) then
        ! Calculate the baseline energy
        allocate (ediff(sim%nt, 5))
        ediff = 0.0_r8

        do i = 1, sim%nt
            if (mod(i, mw%n) .ne. mw%r) cycle
            ediff(i, 1) = sim%stat%potential_energy(i)
            ediff(i, 2) = sim%stat%potential_energy(i) -  pebuf(i, 1)
            ediff(i, 3) = sim%stat%potential_energy(i) -  pebuf(i, 1) - pebuf(i, 4)
            ediff(i, 4) = sim%stat%potential_energy(i) -  pebuf(i, 1) - pebuf(i, 4) - pebuf(i, 2)
            ediff(i, 5) = sim%stat%potential_energy(i) -  pebuf(i, 1) - pebuf(i, 4) - pebuf(i, 2) - pebuf(i, 3)
        end do
        call mw%allreduce('sum', ediff)

        if (T_actual .gt. 1E-5_r8) then
            inverse_kbt = 1.0_r8/lo_kb_Hartree/T_actual
        else
            inverse_kbt = 0.0_r8
        end if

        ! Compute the first and second order cumulants
        do i = 1, 5
            cumulant(1, i) = lo_mean(ediff(:, i))
            cumulant(2, i) = lo_mean((ediff(:, i) - cumulant(1, i))**2)
            cumulant(2, i) = cumulant(2, i)*inverse_kbt*0.5_r8
            cumulant(3, i) = lo_mean((ediff(:, i) - cumulant(1, i))**3)
            cumulant(3, i) = cumulant(3, i)*inverse_kbt**2/6.0_r8
        end do

        ! And normalize it to be per atom
        cumulant = cumulant*to_mev_per_atom
        if (mw%talk) then
            u = open_file('out', 'outfile.cumulants')

            write (u, *) '# Temperature (K) : ', T_actual
            write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
            write (u, *) '# Potential energy [eV / atom]:'
            write (u, "(1X,A,E21.14,1X,A,F21.14)") '                  mean(E): ', cumulant(1, 1), ' upper bound:', cumulant(2, 1)
            write (u, "(1X,A,E21.14,1X,A,F21.14)") '                  mean(E-E2): ', cumulant(1, 2), ' upper bound:', cumulant(2, 2)
            write (u, "(1X,A,E21.14,1X,A,F21.14)") '            mean(E-E2-Epolar): ', cumulant(1, 3), ' upper bound:', cumulant(2, 3)
            if (opts%thirdorder) then
                write (u, "(1X,A,E21.14,1X,A,F21.14)") '        mean(E-E2-Epolar-E3): ', cumulant(1, 4), ' upper bound:', cumulant(2, 4)
            end if
            if (opts%fourthorder) then
                write (u, "(1X,A,E21.14,1X,A,F21.14)") '    mean(E-E2-Epolar-E3-E4): ', cumulant(1, 5), ' upper bound:', cumulant(2, 5)
            end if

            close(u)
            write (*, '(A)') ' ... cumulants writen to `outfile.cumulants`'
        end if

        if(opts%thirdorder .and. (.not. opts%fourthorder)) then
            U0 = cumulant(1,4)
        else if(opts%fourthorder) then
            U0 = cumulant(1,5)
        else
            U0 = cumulant(1,3)
        end if

        ! Now that we have U0 we can get E_total_tdep
        if (mw%talk) then
            u = open_file('out', 'outfile.energies')
            write (u, '(A,A)') '# Unit:      ', 'meV/atom'
            write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
            write (u, *) '# U0 [meV/atom]: ', U0
            write (u, "(A)") '#  conf      Etotal_actual            Etotal_tdep                Epolar              &
                &Epair               Etriplet            Equartet'
            
            do i = 1, sim%nt
                total_energy = sum(pebuf(i,:))*to_mev_per_atom + U0
                write (u, "(1X,I5,6(2X,E20.12))") i, sim%stat%potential_energy(i)*to_mev_per_atom, total_energy, pebuf(i, 4)*to_mev_per_atom, &
                                                    pebuf(i, 1)*to_mev_per_atom, pebuf(i, 2)*to_mev_per_atom, &
                                                    pebuf(i, 3)*to_mev_per_atom
            end do

            ! close outfile.energies
            close (u)

            write (*, '(A)') ' ... energies writen to `outfile.energies`'
        end if

    else
        if (mw%talk) then
            u = open_file('out', 'outfile.energies')
            write (u, '(A,A)') '# Unit:      ', 'meV/atom'
            write (u, *) '# Temperature (K) : ', T_actual
            write (u, '(A,A)') '# no. atoms: ', tochar(ss%na)
            write (u, "(A)") '# True potential energy unknown when using canonical configs, cannot calculate U0 or other cumulants.'
            write (u, "(A)") '#  conf      Epolar              &
                &Epair               Etriplet            Equartet'
            
            do i = 1, opts%nconf
                write (u, "(1X,I5,4(2X,E20.12))") i, pebuf(i, 4)*to_mev_per_atom, pebuf(i, 1)*to_mev_per_atom, &
                                                    pebuf(i, 2)*to_mev_per_atom, pebuf(i, 3)*to_mev_per_atom
            end do

            close (u)

            write (*, '(A)') ' ... energies writen to `outfile.energies`'
        end if
    end if
    

end block epotthings

call mw%destroy()

end program