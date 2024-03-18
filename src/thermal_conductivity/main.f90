#include "precompilerdefinitions"
program thermal_conductivity
!!{!src/thermal_conductivity/manual.md!}
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol
use gottochblandat, only: walltime, tochar, lo_linspace, lo_logspace, lo_mean
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use dump_data, only: lo_dump_gnuplot_2d_real

! unique
use options, only: lo_opts
use scatteringstrengths, only: calculate_scattering_amplitudes
use pbe, only: get_kappa, get_kappa_offdiag, calculate_qs, get_selfconsistent_solution
use phononevents, only: lo_threephononevents, lo_find_all_scattering_events
use mfp, only: lo_mfp, get_cumulative_plots, write_cumulative_plots

implicit none
! Standard, from libolle
type(lo_opts) :: opts
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_phonon_dispersions) :: dr
type(lo_phonon_dos) :: pd
type(lo_crystalstructure) :: uc
class(lo_qpoint_mesh), allocatable :: qp
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
! Unique
type(lo_threephononevents) :: sc
type(lo_mfp) :: mf
! Small stuff
real(r8), dimension(:, :), allocatable :: thermal_cond
real(r8), dimension(:), allocatable :: temperatures
! timers
real(r8) :: timer_init, timer_count, timer_matrixelements, timer_scf
real(r8) :: timer_kappa, timer_qs, timer_cumulative, tt0

! Set up all harmonic properties. That involves reading all the input file,
! creating grids, getting the harmonic properties on those grids.
initharmonic: block
    integer :: i, j
    ! Start MPI and timers
    tt0 = walltime()
    timer_init = tt0
    timer_qs = 0.0_r8
    timer_kappa = 0.0_r8
    timer_scf = 0.0_r8
    timer_cumulative = 0.0_r8
    call mw%init()
    ! Get options
    call opts%parse()
    if (mw%r .ne. 0) opts%verbosity = -100
    ! Init memory tracker
    call mem%init()

    if (mw%talk) write (*, *) '... using ', tochar(mw%n), ' MPI ranks'
    ! There is a bunch of stuff that all ranks need, first the unit cell:
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=opts%timereversal)
    if (mw%talk) write (*, *) '... read unitcell poscar'

    ! Perhaps non-natural isotope distribution
    ! write (*, *) 'FIXME OUTPUT UNITS'
    if (opts%readiso) then
        if (mw%talk) write (*, *) '... reading isotope distribution from file'
        call uc%readisotopefromfile()
        if (mw%talk) then
            do i = 1, uc%na
                do j = 1, uc%isotope(i)%n
                    write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                        j, uc%isotope(i)%conc(j), uc%isotope(i)%mass(j)
                end do
                write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                    trim(uc%atomic_symbol(uc%species(i))), uc%isotope(i)%mean_mass, &
                    uc%isotope(i)%disorderparameter
            end do
        end if
    elseif (opts%verbosity .gt. 0) then
        do i = 1, uc%na
            do j = 1, uc%isotope(i)%n
                write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                    j, uc%isotope(i)%conc(j), uc%isotope(i)%mass(j)
            end do
            write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                trim(uc%atomic_symbol(uc%species(i))), uc%isotope(i)%mean_mass, &
                uc%isotope(i)%disorderparameter
        end do
    end if

    ! Read the force constants
    call fc%readfromfile(uc, 'infile.forceconstant', mem, -1)
    if (mw%talk) write (*, *) '... read second order forceconstant'
    call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
    if (mw%talk) write (*, *) '... read third order forceconstant'

    ! Get q-point mesh
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, &
                           headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity, nosym=.not. opts%qpsymmetry)

    ! Get frequencies in the whole BZ
    if (mw%talk) then
        write (*, *) '... getting the full dispersion relations'
    end if
    call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)
    ! Also the phonon DOS, for diagnostics
    call pd%generate(dr, qp, uc, mw, mem, verbosity=opts%verbosity, &
                     sigma=opts%sigma, n_dos_point=opts%mfppts*2, integrationtype=opts%integrationtype)

    ! Make sure it's stable, no point in going further if it is unstable.
    if (dr%omega_min .lt. -lo_freqtol) then
        write (*, *) ''
        write (*, *) 'FOUND UNSTABLE MODES. WILL STOP NOW.'
        call mpi_barrier(mw%comm, mw%error)
        call mpi_finalize(lo_status)
        stop
    end if

    ! now I have all harmonic things, stop the init timer
    timer_init = walltime() - timer_init
end block initharmonic

! Get the integration weights and matrix elements
weights_elements: block
    real(r8) :: t0
    t0 = walltime()
    timer_count = walltime()
    call lo_find_all_scattering_events(sc, qp, dr, uc, mw, mem, opts%sigma, opts%thres, opts%integrationtype, &
                                       opts%correctionlevel, opts%mfp_max, opts%isotopescattering)
    call mpi_barrier(mw%comm, mw%error)

    ! stop counting timer, start matrixelement timer
    timer_count = walltime() - timer_count
    timer_matrixelements = walltime()

    ! Calculate scattering amplitudes
    t0 = walltime()
    call calculate_scattering_amplitudes(uc, qp, sc, dr, fct, mw)
    call mpi_barrier(mw%comm, mw%error)
    ! stop matrix element timer, start some other timer
    timer_matrixelements = walltime() - timer_matrixelements
    if (mw%talk) write (*, *) 'Counted and got scattering amplitudes in ', tochar(walltime() - t0)
end block weights_elements

! Make space and initialize everything to calculate thermal conductivity
initkappa: block
    integer :: i

    ! space to store the actual thermal conductivity
    allocate (thermal_cond(10, opts%trangenpts))
    thermal_cond = 0.0_r8

    ! temperature axis
    allocate (temperatures(opts%trangenpts))
    if (opts%logtempaxis) then
        call lo_logspace(opts%trangemin, opts%trangemax, temperatures)
    else
        call lo_linspace(opts%trangemin, opts%trangemax, temperatures)
    end if
    ! Setup the mean-free-path vs kappa plots.
    ! how many points on the x-axis?
    mf%np = opts%mfppts
    ! how many temperatures?
    mf%nt = opts%trangenpts
    ! one plot for each temperature
    allocate (mf%temp(mf%nt))

    ! Make some space to keep intermediate values
    do i = 1, qp%n_irr_point
        allocate (dr%iq(i)%p_plus(dr%n_mode))
        allocate (dr%iq(i)%p_minus(dr%n_mode))
        allocate (dr%iq(i)%p_iso(dr%n_mode))
        allocate (dr%iq(i)%qs(dr%n_mode))
        allocate (dr%iq(i)%linewidth(dr%n_mode))
        allocate (dr%iq(i)%F0(3, dr%n_mode))
        allocate (dr%iq(i)%Fn(3, dr%n_mode))
        allocate (dr%iq(i)%mfp(3, dr%n_mode))
        allocate (dr%iq(i)%scalar_mfp(dr%n_mode))
        dr%iq(i)%linewidth = 0.0_r8
        dr%iq(i)%p_plus = 0.0_r8
        dr%iq(i)%p_minus = 0.0_r8
        dr%iq(i)%p_iso = 0.0_r8
        dr%iq(i)%qs = 0.0_r8
        dr%iq(i)%F0 = 0.0_r8
        dr%iq(i)%Fn = 0.0_r8
        dr%iq(i)%mfp = 0.0_r8
        dr%iq(i)%scalar_mfp = 0.0_r8
    end do
    do i = 1, qp%n_full_point
        allocate (dr%aq(i)%kappa(3, 3, dr%n_mode))
        dr%aq(i)%kappa = 0.0_r8
    end do
end block initkappa

! Iteratively solve the BTE for each temperature. Additionally calculate
! mean free path plots and things like that.
getkappa: block
    real(r8), dimension(3, 3) :: kappa, kappa_offdiag, m0
    real(r8) :: t0
    integer :: i

    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'THERMAL CONDUCTIVITY'
        if (opts%scfiterations .eq. 0) then
            write (*, "(3X,A11,6(1X,A14))") 'Temperature', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        end if
    end if
    ! Main loop over temperatures to solve the BTE
    do i = 1, opts%trangenpts

        ! I might get a silly tiny temperature, then things will break.
        if (temperatures(i) .lt. lo_temperaturetol) then
            kappa = 0.0_r8
            thermal_cond(1, i) = temperatures(i)
            thermal_cond(2:10, i) = 0.0_r8
            cycle
        end if

        ! Scattering rates
        t0 = walltime()
        call calculate_qs(qp, sc, dr, temperatures(i), mw, mem)
        timer_qs = timer_qs + walltime() - t0

        call get_kappa_offdiag(dr, qp, uc, temperatures(i), fc, mem, mw, kappa_offdiag)

        ! Get the self-consistent solution
        call mpi_barrier(mw%comm, mw%error)
        if (opts%scfiterations .gt. 0) then
            if (mw%talk) then
                write (*, *) ''
                write (*, *) 'Temperature: ', tochar(temperatures(i))
                write (*, "(1X,A4,6(1X,A14),2X,A10)") 'iter', &
                    'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   ', 'DeltaF/F'
            end if
            t0 = walltime()
            call get_selfconsistent_solution(sc, dr, qp, uc, temperatures(i), opts%scfiterations, opts%scftol, mw, mem)
            !call get_selfconsistent_solution(sc,dr,qp,uc,mw,temperatures(i),opts%scfiterations,opts%scftol)
            timer_scf = timer_scf + walltime() - t0
            call get_kappa(dr, qp, uc, temperatures(i), kappa)
            m0 = kappa*lo_kappa_au_to_SI
            if (mw%talk) write (*, "(5X,6(1X,F14.4),2X,E10.3)") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        else
            call get_kappa(dr, qp, uc, temperatures(i), kappa)
            m0 = kappa*lo_kappa_au_to_SI
            if (mw%talk) write (*, "(1X,F12.3,6(1X,F14.4),2X,E10.3)") &
                temperatures(i), m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        end if

        if (mw%talk) then
            m0 = kappa_offdiag*lo_kappa_au_to_SI
            write (*, "(1X,A25)") 'Off diagonal contribution'
            write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
            write (*, "(5X,6(1X,F14.4),2X,E10.3)") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
            m0 = (kappa+kappa_offdiag)*lo_kappa_au_to_SI
            write (*, "(1X,A19)") 'Total kappa (W/m/K)'
            write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
            write (*, "(5X,6(1X,F14.4),2X,E10.3)") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        end if

        ! Store thermal conductivity tensor
        thermal_cond(1, i) = temperatures(i)
        thermal_cond(2, i) = (kappa(1, 1) + kappa_offdiag(1, 1))*lo_kappa_au_to_SI
        thermal_cond(3, i) = (kappa(2, 2) + kappa_offdiag(2, 2))*lo_kappa_au_to_SI
        thermal_cond(4, i) = (kappa(3, 3) + kappa_offdiag(3, 3))*lo_kappa_au_to_SI
        thermal_cond(5, i) = (kappa(1, 3) + kappa_offdiag(1, 3))*lo_kappa_au_to_SI
        thermal_cond(6, i) = (kappa(2, 3) + kappa_offdiag(2, 3))*lo_kappa_au_to_SI
        thermal_cond(7, i) = (kappa(1, 2) + kappa_offdiag(1, 2))*lo_kappa_au_to_SI
        thermal_cond(8, i) = (kappa(3, 1) + kappa_offdiag(3, 1))*lo_kappa_au_to_SI
        thermal_cond(9, i) = (kappa(3, 2) + kappa_offdiag(3, 2))*lo_kappa_au_to_SI
        thermal_cond(10, i) = (kappa(2, 1) + kappa_offdiag(2, 1))*lo_kappa_au_to_SI

        ! Calculate the cumulative plots
        t0 = walltime()
        call mpi_barrier(mw%comm, mw%error)
        call get_cumulative_plots(mf%temp(i), qp, dr, pd, uc, opts%mfppts, temperatures(i), opts%sigma, kappa, mw, mem)

        timer_cumulative = timer_cumulative + walltime() - t0
        call mpi_barrier(mw%comm, mw%error)
    end do
end block getkappa

! dump things to file and print timings
finalize_and_write: block
    real(r8) :: t0

    ! Write thermal conductivity to file
    if (mw%talk) call lo_dump_gnuplot_2d_real(thermal_cond, 'outfile.thermal_conductivity', &
                                              ylabel='\kappa W/mK', xlabel='Temperature (K)')

    ! Write the cumulative kappa
    if (mw%talk) call write_cumulative_plots(mf, pd, uc, 'thz', 'outfile.cumulative_kappa.hdf5', opts%verbosity)

    ! Maybe dump data on a grid
    if (mw%talk .and. opts%dumpgrid .and. opts%trangenpts .eq. 1) then
        write (*, *) '... dumping auxiliary data to files:'
        call dr%write_to_hdf5(qp, uc, 'outfile.grid_thermal_conductivity.hdf5', mem, temperatures(1))
    end if

    ! sum up the total time
    if (mw%talk) tt0 = walltime() - tt0

    ! Print timings
    if (mw%talk) then
        t0 = timer_init + timer_count + timer_matrixelements + timer_qs + timer_kappa + timer_scf + timer_cumulative
        write (*, *) ' '
        write (*, *) 'Timings:'
        write (*, "(A,F12.3,A,F7.3,A)") '            initialization:', timer_init, ' s, ', real(timer_init*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '       integration weights:', timer_count, ' s, ', real(timer_count*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '           matrix elements:', timer_matrixelements, &
            ' s, ', real(timer_matrixelements*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '            QS calculation:', timer_qs, ' s, ', real(timer_qs*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '                     kappa:', timer_kappa, ' s, ', real(timer_kappa*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '          self consistency:', timer_scf, ' s, ', real(timer_scf*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '          cumulative plots:', timer_cumulative, &
            ' s, ', real(timer_cumulative*100/tt0), '%'
        write (*, "(A,F12.3,A)") '                     total:', tt0, ' seconds'
    end if
end block finalize_and_write

! And we are done!
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)

end program
