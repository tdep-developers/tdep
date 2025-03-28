#include "precompilerdefinitions"
program effective_hamiltonian
!!{!src/effective_hamiltonian/manual.md!}
use konstanter, only: r8, lo_tol, lo_kb_hartree, lo_bohr_to_A, lo_frequency_Hartree_to_THz, lo_eV_to_Hartree, &
                      lo_exitcode_physical, lo_exitcode_param, lo_temperaturetol
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: tochar, walltime, lo_stop_gracefully, open_file

use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_mdsim, only: lo_mdsim

use lo_epot, only: statistical_sampling, setup_potential_energy_differences

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

logical :: generate_configs = .false.

init: block
    ! Get CLI options
    call opts%parse()

    if (opts%nconf .gt. 0) generate_configs = .true.

    if (mw%talk) then
        write (*, *) 'Recap of the parameters governing the calculation'
        write (*, '(1X,A40,L3)') 'Thirdorder scattering                   ', opts%thirdorder
        write (*, '(1X,A40,L3)') 'Fourthorder scattering                  ', opts%fourthorder
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

    call fc%readfromfile(uc, 'infile.forceconstant', mem, verbosity=-1)
    if (opts%thirdorder) call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
    if (opts%fourthorder) call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
    if (mw%talk) write (*, *) '... read forceconstants'

    if (generate_configs) then
        if (opts%temperature .lt. 0.0_r8) then
            call lo_stop_gracefully(['You need to provide a positive temperature to generate configurations'], &
                                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
        if (mw%talk .and. opts%quantum) write (*, *) '... will generate quantum configurations'
        if (mw%talk .and. opts%quantum) write (*, *) '... will generate classical configurations'
    end if

    call pot%setup(uc, ss, fc2, fc3, fc4, mw, verbosity)
    if (mw%talk) write (*, *) '... setup potential energy calculator'

    if (.not. generate_configs) then
        readmag = .false.
        readdiel = .false.
        call sim%read_from_file(verbosity=opts%verbosity + 2, stride=opts%stride, mw=mw)
        if (mw%talk) write (*, *) '... parsed infile.positions'
    end if

end block init


energy : block

    real(r8), dimension(:, :), allocatable :: ebuf

    if (generate_configs) then
        call pot%statistical_sampling(uc, ss, fc2, otps%nconf, opts%temperature, mw, mem, opts%verbosity)
    else
        ! use parsed sim to evalulate energies
    end if


end block energy