#include "precompilerdefinitions"
program canonical_configuration
!!{!src/canonical_configuration/manual.md!}
use konstanter, only: r8, lo_tol, lo_kb_hartree, lo_bohr_to_A, lo_frequency_Hartree_to_THz, lo_eV_to_Hartree, &
                      lo_exitcode_physical, lo_temperaturetol, lo_exitcode_param
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: tochar, walltime, lo_stop_gracefully, open_file

use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_jij_secondorder, only: lo_jij_secondorder
use type_mdsim, only: lo_mdsim
use lo_randomnumbers, only: lo_mersennetwister

use options, only: lo_opts
use semirandom, only: generate_semirandom_configurations

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: ss, uc, sqsid, sqsrl
type(lo_forceconstant_secondorder) :: fc, fcss
type(lo_jij_secondorder) :: jij
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
integer, dimension(:), allocatable :: alloy_permutation
type(lo_mersennetwister), save :: tw

! Get the necessary things first
init: block
    real(r8) :: seed
    ! Get CLI options
    call opts%parse()
    call mw%init()
    call mem%init()
    if (.not. mw%talk) opts%verbosity = -100

    ! initialize random numbers
    if (opts%seed < 0) then
        ! fully random
        seed = walltime()
        if (mw%talk) print *, '... walltime() used to initialize random state'
    else
        ! use inverse of seed because integer part is ignored in lo_randomnumbers.f90
        if (mw%talk) print *, '... RANDOM SEED: ', opts%seed
        if (opts%seed == 0) then
            seed = 0.0_r8
        else
            seed = 1.0_r8/float(opts%seed)
        end if
    end if
    call tw%init(iseed=mw%r, rseed=seed)

    ! Just make sure no clever person starts running this in parallel.
    if (mw%n .gt. 1) then
        call lo_stop_gracefully(['Do not run this in parallel.'], lo_exitcode_param, __FILE__, __LINE__)
    end if

    ! Read structures
    if (mw%talk) write (*, *) '... reading infiles'
    call ss%readfromfile('infile.ssposcar', verbosity=opts%verbosity)
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    ! Match the supercell to the unitcell, always a good idea
    call uc%classify('spacegroup', timereversal=.true.)
    call ss%classify('supercell', uc)

    ! get a forceconstant somehow
    if (opts%debye_temperature .gt. 0.0_r8) then
        call fc%fake_forceconstant(uc, ss, debye_temperature=opts%debye_temperature, verbosity=opts%verbosity)
        write (*, *) '... constructed fake forceconstant corresponding to Td = ', tochar(opts%debye_temperature), 'K'
        call fc%writetofile(uc, 'outfile.fakeforceconstant')
        write (*, *) '... wrote it to "outfile.fakeforceconstant", check that the frequency range is reasonable'
    elseif (opts%maximum_frequency .gt. 0.0_r8) then
        call fc%fake_forceconstant(uc, ss, maximum_frequency=opts%maximum_frequency, verbosity=opts%verbosity)
        write (*, *) '... constructed fake forceconstant corresponding to max(omega) = ', &
            tochar(opts%maximum_frequency*lo_frequency_Hartree_to_THz), ' THz'
        call fc%writetofile(uc, 'outfile.fakeforceconstant')
        write (*, *) '... wrote it to "outfile.fakeforceconstant", check that dispersions look reasonable'
    else
        ! read the forceconstant from file
        call fc%readfromfile(uc, 'infile.forceconstant', mem, opts%verbosity)
    end if

    ! Create fake magnetic exchange interactions
    if (abs(opts%exchange_J) .gt. lo_tol) then
        ! Convert from the meV in the input to Hartree
        opts%exchange_J = opts%exchange_J*lo_eV_to_Hartree/1000
        opts%exchange_J = opts%exchange_J/(opts%mean_moment**2)
        call jij%fake_jij(uc, opts%exchange_J, opts%mean_moment)
        call jij%writetofile(uc, 'outfile.fakejij')
    end if

    ! Remap force constant to supercell
    if (uc%info%alloy) then
        write (*, *) '... alloy detected'
        ! We need the ideal SQS, and the relaxed SQS. These are represented as normal
        ! structures, without any alloy stuff.
        call sqsid%readfromfile('infile.sqs_ideal', verbosity=opts%verbosity)
        call sqsrl%readfromfile('infile.sqs_relaxed', verbosity=opts%verbosity)
        ! Figure out how it is supposed to be permuted
        allocate (alloy_permutation(ss%na))
        alloy_permutation = 0
        call ss%alloy_site_permutation(sqsid, alloy_permutation)
        ! Permute the supercell
        call ss%permute_positions(alloy_permutation, forward=.false.)
        ! Then do the magical remapping
        call fc%remap(uc, ss, fcss)
        ! And permute back, just to be on the safe side
        call ss%permute_positions(alloy_permutation, forward=.true.)
    else
        ! in normal case just remap it
        call fc%remap(uc, ss, fcss)
    end if
    if (mw%talk) write (*, *) '... remapped fc'

    if ((opts%zpm .eqv. .false.) .and. (opts%temperature .lt. lo_temperaturetol)) then
        call lo_stop_gracefully(['For classical statistics temperature has to be nonzero'], lo_exitcode_physical, __FILE__, __LINE__)
    end if
end block init

if (opts%semirandom) then
    ! Generate configurations the fancy way
    if (opts%dielcutoff2 < lo_tol) then
        call lo_stop_gracefully(['diel. cutoff (-dc2) has to be nonzero when using --semirandom'], &
                                lo_exitcode_param, __FILE__, __LINE__)
    end if
    call generate_semirandom_configurations(uc, ss, fc, fcss, opts%temperature, opts%zpm, opts%dielcutoff2, opts%dielcutoff3, opts%output_format, opts%nconf, mw, mem, opts%verbosity)
end if

! fkdev: write each mode displacement separately
if (opts%modes) then
    opts%nconf = 6*fcss%na
    if (mw%talk) write (*, *) '!!! generate 2 (+/-) displacements for each mode'
    if (mw%talk) write (*, '(1X,A,I5,A)') '--> ', opts%nconf, ' structures'
end if

if (opts%modes) then
    write (*, "(1X,A)") '  no.   +/-      ek(K)          ep(K)          <ek/ep>         T(K)          <T>(K)      <msd>(A)         frequency(THz)'
else
    write (*, "(1X,A)") '  no.   +/-      ek(K)          ep(K)          <ek/ep>         T(K)          <T>(K)      <msd>(A)'
end if
dumpconf: block
    type(lo_mdsim) :: sim
    type(lo_crystalstructure) :: p
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(:, :, :, :), allocatable :: polar_fc
    real(r8) :: ep, ek, temp, avgtemp, avgmsd, msd, ratio, rek, rep
    integer :: i, a1, a2, iconf, imode
    logical :: invert
    character(len=1000) :: fname
    character(len=1) :: invert_char

    ! Get a copy of the structure
    if (uc%info%alloy) then
        ! Copy of the relaxed SQS
        p = sqsrl
    else
        ! Just a normal copy
        p = ss
    end if

    ! Get a copy of the polar forceconstant to evaluate potential energy
    if (fc%polar) then
        ! Can't do polar and alloys at the same time, for now.
        if (uc%info%alloy) then
            call lo_stop_gracefully(['Can not do polar alloys at the moment. Will fix.'], lo_exitcode_param, __FILE__, __LINE__)
        end if

        lo_allocate(polar_fc(3, 3, p%na, p%na))
        polar_fc = 0.0_r8
        call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, polar_fc, 1E-10_r8)
    end if

    ! Maybe create a fake dataset?
    if (opts%fakesim) then
        call sim%init_empty(uc, ss, opts%nconf, opts%temperature, &
                            timestep=-1.0_r8, &
                            magnetic_moments=.false., &
                            variable_lattice=.false.)
    end if

    ! Now start generating configurations
    ratio = 0.0_r8
    avgtemp = 0.0_r8
    avgmsd = 0.0_r8
    rek = 0.0_r8
    rep = 0.0_r8
    imode = 1
    invert = .false.

    do iconf = 1, opts%nconf
        ! reset the structure
        if (uc%info%alloy) then
            ! Copy of the relaxed SQS
            p%r = sqsrl%r
            p%rcart = sqsrl%rcart
        else
            ! Just a normal copy
            p%r = ss%r
            p%rcart = ss%rcart
        end if
        p%v = 0.0_r8
        p%u = 0.0_r8

        ! create the thermal configuration
        if (opts%modes) then
            ! create +/- displacement for each mode
            imode = (iconf - 1)/2 + 1
            if (mod(iconf, 2) == 0) then
                invert = .true.
            else
                invert = .false.
            end if
            call fcss%initialize_cell(p, uc, fc, opts%temperature, opts%zpm, .false., opts%mindist, mw, imode=imode, invert=invert, tw=tw)
        else
            call fcss%initialize_cell(p, uc, fc, opts%temperature, opts%zpm, .false., opts%mindist, mw, tw=tw)
        end if

        ! dump to file
        select case (opts%output_format)
        case (1) ! vasp output
            fname = 'contcar_conf'//tochar(iconf, 4)
            call p%writetofile(trim(fname), opts%output_format, write_velocities=.true.)
        case (2) ! abinit output
            fname = 'abinput_conf'//tochar(iconf, 4)
            call p%writetofile(trim(fname), opts%output_format, write_velocities=.true.)
        case (3) ! LAMMPS output
            fname = 'lammps_conf'//tochar(iconf, 4)
            call p%writetofile(trim(fname), opts%output_format, write_velocities=.true.)
        case (4) ! AIMS output
            fname = 'aims_conf'//tochar(iconf, 4)
            call p%writetofile(trim(fname), opts%output_format, write_velocities=.true.)
        case (5) ! Siesta output
            fname = 'siesta_conf'//tochar(iconf, 4)
            call p%writetofile(trim(fname), opts%output_format, write_velocities=.true.)
        end select

        ! just measure some stuff, for no good reason.
        ! First kinetic energy:
        ek = p%kinetic_energy()/(p%na)
        ! Then forces
        p%f = 0.0_r8
        do a1 = 1, p%na
            do i = 1, fcss%atom(a1)%n
                a2 = fcss%atom(a1)%pair(i)%i2
                p%f(:, a1) = p%f(:, a1) - matmul(fcss%atom(a1)%pair(i)%m, p%u(:, a2))
            end do
        end do
        if (fc%polar) then
            do a1 = 1, p%na
            do a2 = 1, p%na
                p%f(:, a1) = p%f(:, a1) - matmul(polar_fc(:, :, a1, a2), p%u(:, a2))
            end do
            end do
        end if

        ! Then potential energy
        ep = 0.0_r8
        do a1 = 1, p%na
            ep = ep - dot_product(p%u(:, a1), p%f(:, a1))*0.5_r8/real(p%na, r8)
        end do

        ! Now, I might want to create a fake dataset to play around with.
        if (opts%fakesim) then
            m0 = 0.0_r8 ! No stress
            if (sim%alloy) then
                ! In case it's an alloy. Hmmm. Probably permute the things.
                call p%permute_positions(alloy_permutation, forward=.true.)
                call sqsid%permute_positions(alloy_permutation, forward=.true.)
                ! Store the timestep
                call sim%add_timestep(p%r, p%f, ep*p%na, ek*p%na, ek/(1.5_r8*lo_kb_hartree), m0, atomic_numbers=p%atomic_number, reference_positions=sqsid%r)
                ! Permute back
                call p%permute_positions(alloy_permutation, forward=.false.)
                call sqsid%permute_positions(alloy_permutation, forward=.false.)
            else
                ! Just store normally
                call sim%add_timestep(p%r, p%f, ep*p%na, ek*p%na, ek/(1.5_r8*lo_kb_hartree), m0)
            end if
        end if

        ! Accumulate averages
        rek = rek + ek
        rep = rep + ep
        ratio = rek/rep
        temp = (ek + ep)/(3*lo_kb_hartree)
        msd = 0.0_r8
        do a1 = 1, p%na
            msd = msd + norm2(p%u(:, a1))/p%na
        end do

        avgmsd = avgmsd + msd
        avgtemp = avgtemp + temp
        if (invert) then
            invert_char = '-'
        else
            invert_char = '+'
        end if
        if (opts%modes) then
            write (*, "(1X,I4,5X,A,5(2X,F13.5),2X,F15.8,2X,F10.3)") iconf, invert_char, &
                ek/lo_kb_hartree/1.5_r8, ep/lo_kb_hartree/1.5_r8, ratio, temp, &
                avgtemp/iconf, avgmsd*lo_bohr_to_A/iconf, fcss%omega(imode)*lo_frequency_Hartree_to_THz
        else
            write (*, "(1X,I4,5X,A,5(2X,F13.5),2X,F15.8)") iconf, invert_char, &
                ek/lo_kb_hartree/1.5_r8, ep/lo_kb_hartree/1.5_r8, ratio, temp, &
                avgtemp/iconf, avgmsd*lo_bohr_to_A/iconf
        end if
    end do

    ! And, at the end, maybe dump the fake simulation
    if (opts%fakesim) then
        if (fc%polar) then
            ! Make sure to store the polar stuff
            call sim%write_to_hdf5(uc, ss, 'outfile.fakesim.hdf5', opts%verbosity, eps=fc%loto%eps, Z=fc%loto%born_effective_charges)
        else
            ! Just write normally
            call sim%write_to_hdf5(uc, ss, 'outfile.fakesim.hdf5', opts%verbosity)
        end if
    end if
end block dumpconf

end program
