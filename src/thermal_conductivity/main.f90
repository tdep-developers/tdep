#include "precompilerdefinitions"
program thermal_conductivity
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol, lo_m_to_Bohr, lo_emu_to_amu
use gottochblandat, only: walltime, tochar, open_file
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use lo_timetracker, only: lo_timer

use options, only: lo_opts
use kappa, only: get_kappa, get_kappa_offdiag, iterative_solution, symmetrize_kappa
use scattering, only: lo_scattering_rates
use cumulative_kappa, only: lo_cumulative_kappa

implicit none

! Standard from libolle
type(lo_opts) :: opts
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_forceconstant_fourthorder) :: fcf
type(lo_phonon_dispersions) :: dr
type(lo_phonon_dos) :: pd
type(lo_crystalstructure) :: uc
class(lo_qpoint_mesh), allocatable :: qp
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_timer) :: tmr_init, tmr_scat, tmr_kappa, tmr_tot
! The scattering rates
type(lo_scattering_rates) :: sr
type(lo_cumulative_kappa) :: mf
real(r8) :: t0

! Set up all harmonic properties. That involves reading all the input file,
! creating grids, getting the harmonic properties on those grids.
initharmonic: block
    integer :: i, j, q1
    ! Start MPI and timers
    call mw%init()
    t0 = walltime()
    ! Start the initialization timer
    call tmr_tot%start()
    call tmr_init%start()
    ! Get options
    call opts%parse()
    if (.not. mw%talk) opts%verbosity = -100
    ! Init memory tracker
    call mem%init()

    if (mw%talk) then
        write (*, *) 'Recap of the parameters governing the calculation'
        write (*, '(1X,A40,F20.12)') 'Temperature                             ', opts%temperature
        write (*, '(1X,A40,L3)') 'Thirdorder scattering                   ', opts%thirdorder
        write (*, '(1X,A40,L3)') 'Fourthorder scattering                  ', opts%fourthorder
        write (*, '(1X,A40,L3)') 'Isotope scattering                      ', opts%isotopescattering
        write (*, '(1X,A40,L3)') 'Classical limit                         ', opts%classical
        write (*, '(1X,A40,I4,I4,I4)') 'full q-point grid                       ', opts%qgrid
        write (*, '(1X,A40,I4,I4,I4)') 'Monte-Carlo 3rd order q-point grid      ', opts%qg3ph
        write (*, '(1X,A40,I4,I4,I4)') 'Monte-Carlo 4th order q-point grid      ', opts%qg4ph
        write (*, '(1X,A40,I5)') 'Max number of iteration                 ', opts%itermaxsteps
        write (*, '(1X,A40,E20.12)') 'Max mean free path (in m)               ', opts%mfp_max/lo_m_to_Bohr
        write (*, '(1X,A40,E20.12)') 'Tolerance for the iterative solution    ', opts%itertol
        select case (opts%integrationtype)
        case (1)
            write (*, '(1X,A40,2X,A)') 'Integration type                        ', 'Gaussian with fixed broadening'
        case (2)
            write (*, '(1X,A40,2X,A)') 'Integration type                        ', 'Adaptive Gaussian'
        write (*, '(1X,A40,E20.12)') 'Sigma factor for gaussian smearing      ', opts%sigma
        end select
        write (*, '(1X,A40,2X,I4)') 'Number of mfp point for cumulative kappa ', opts%mfppts
        write (*, '(1X,A40,2X,I4)') 'Number of freq point for spectral kappa ', opts%freqpts
        select case (opts%dosintegrationtype)
        case(1)
            write (*, '(1X,A40,2X,A)') 'Integration type for spectral kappa       ', 'Gaussian with fixed broadening'
        case(2)
            write (*, '(1X,A40,2X,A)') 'Integration type for spectral kappa       ', 'Adaptive Gaussian'
        case(3)
            write (*, '(1X,A40,2X,A)') 'Integration type for spectral kappa       ', 'Tetrahedron'
        end select
        write (*, '(1X,A40,I10)') 'Number of MPI ranks                     ', mw%n
        if (opts%seed .gt. 0) write(*, '(1X,A40,I10)') 'Random seed                             ', opts%seed
        write (*, *) ''
    end if

    if (mw%talk) write (*, *) 'Initialize calculation'
    ! There is a bunch of stuff that all ranks need, first the unit cell:
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=opts%timereversal)
    if (mw%talk) write (*, *) '... read unitcell poscar'

    ! Perhaps non-natural isotope distribution
    if (opts%readiso) then
        if (mw%talk) write (*, *) '... reading isotope distribution from file'
        call uc%readisotopefromfile()
        if (mw%talk) then
            do i = 1, uc%na
                do j = 1, uc%isotope(i)%n
                    write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                        j, uc%isotope(i)%conc(j), uc%isotope(i)%mass(j) * lo_emu_to_amu
                end do
                write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                    trim(uc%atomic_symbol(uc%species(i))), uc%isotope(i)%mean_mass * lo_emu_to_amu, &
                    uc%isotope(i)%disorderparameter
            end do
        end if
    elseif (mw%talk .and. opts%verbosity .gt. 0) then
        do i = 1, uc%na
            do j = 1, uc%isotope(i)%n
                write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                    j, uc%isotope(i)%conc(j), uc%isotope(i)%mass(j) * lo_emu_to_amu
            end do
            write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                trim(uc%atomic_symbol(uc%species(i))), uc%isotope(i)%mean_mass * lo_emu_to_amu, &
                uc%isotope(i)%disorderparameter
        end do
    end if

    ! Read the force constants
    call fc%readfromfile(uc, 'infile.forceconstant', mem, -1)
    if (mw%talk) write (*, *) '... read second order forceconstant'
    call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
    if (mw%talk) write (*, *) '... read third order forceconstant'
    if (opts%fourthorder) then
        call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
        if (mw%talk) write (*, *) '... read fourth order forceconstant'
    end if

    call tmr_init%tock('read input files')

    if (mw%talk) write (*, *) '... generating q-point mesh'
    ! Get q-point mesh
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, &
                           headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity, nosym=.not. opts%qpsymmetry)

    call tmr_init%tock('generated q-mesh')

    ! Get frequencies in the whole BZ
    if (mw%talk) write (*, *) '... generating harmonic properties on the q-point mesh'
    call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

    ! Make sure it's stable, no point in going further if it is unstable.
    if (dr%omega_min .lt. -lo_freqtol) then
        write (*, *) ''
        write (*, *) 'FOUND UNSTABLE MODES. WILL STOP NOW.'
        call mpi_barrier(mw%comm, mw%error)
        call mpi_finalize(lo_status)
        stop
    end if

    ! Make some space to keep intermediate values
    do q1 = 1, qp%n_irr_point
        allocate (dr%iq(q1)%linewidth(dr%n_mode))
        allocate (dr%iq(q1)%F0(3, dr%n_mode))
        allocate (dr%iq(q1)%Fn(3, dr%n_mode))
        allocate (dr%iq(q1)%qs(dr%n_mode))
        allocate (dr%iq(q1)%mfp(3, dr%n_mode))
        allocate (dr%iq(q1)%scalar_mfp(dr%n_mode))
        allocate (dr%iq(q1)%kappa(3, 3, dr%n_mode))
        dr%iq(q1)%linewidth = 0.0_r8
        dr%iq(q1)%F0 = 0.0_r8
        dr%iq(q1)%Fn = 0.0_r8
        dr%iq(q1)%qs = 0.0_r8
        dr%iq(q1)%mfp = 0.0_r8
        dr%iq(q1)%scalar_mfp = 0.0_r8
        dr%iq(q1)%kappa = 0.0_r8
    end do
    call tmr_init%tock('harmonic dispersions')
    call tmr_tot%tock('initialization')

    ! now I have all harmonic things, stop the init timer
    t0 = walltime() - t0
    if (mw%talk) write (*, "(1X,A,F12.3,A)") '... done in ', t0, ' s'
    call tmr_init%stop()
end block initharmonic

get_scattering_rates: block
    call tmr_scat%start()
    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'Calculating scattering events'
    end if

    t0 = walltime()
    call sr%generate(qp, dr, uc, fct, fcf, opts, tmr_scat, mw, mem)
    t0 = walltime() - t0

    call tmr_scat%stop()
    call tmr_tot%tock('scattering computation')

    if (mw%talk) write (*, "(1X,A,F12.3,A)") '... done in ', t0, ' s'
end block get_scattering_rates

blockkappa: block
    real(r8), dimension(3, 3) :: kappa_iter, kappa_offdiag, kappa_sma, m0
    real(r8) :: t0
    integer :: i, u, q1, b1

    call tmr_kappa%start()
    t0 = walltime()

    ! I might get a silly tiny temperature, then things will break.
    if (opts%temperature .lt. lo_temperaturetol) then
        kappa_iter = 0.0_r8
        kappa_sma = 0.0_r8
        kappa_offdiag = 0.0_r8
    end if

    if (mw%talk) write (*, *) ''
    if (mw%talk) write (*, *) 'Thermal conductivity calculation'

    if (mw%talk) write (*, *) '... computing kappa in the single mode approximation'
    call get_kappa(dr, qp, uc, opts%temperature, opts%classical, kappa_sma)
    call tmr_kappa%tock('single mode approximation')
    if (mw%talk) write (*, *) '... computing off diagonal (coherence) contribution'
    call get_kappa_offdiag(dr, qp, uc, fc, opts%temperature, opts%classical, mem, mw, kappa_offdiag)
    call tmr_kappa%tock('off-diagonal contribution')
    if (opts%itermaxsteps .gt. 0) then
        if (mw%talk) then
            write (*, *) '... solving iteratively the collective contribution'
            write (*, "(1X,A4,6(1X,A14),2X,A10)") 'iter', &
                'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   ', 'DeltaF/F'
        end if
        t0 = walltime()
        call iterative_solution(sr, dr, qp, uc, opts%temperature, opts%itermaxsteps, opts%itertol, opts%classical, mw, mem)
        t0 = walltime() - t0
        if (mw%talk) write (*, "(1X,A,F12.3,A)") '... done in ', t0, ' s'
        call tmr_kappa%tock('collective contribution')
    end if
    call get_kappa(dr, qp, uc, opts%temperature, opts%classical, kappa_iter)
    if (mw%talk) write (*, *) ''
    if (mw%talk) write (*, *) '... symmetrizing the thermal conductivity tensors'
    call symmetrize_kappa(kappa_iter, uc)
    call symmetrize_kappa(kappa_offdiag, uc)
    call symmetrize_kappa(kappa_sma, uc)
    call tmr_kappa%tock('symmetrization')
    if (mw%talk) then
        ! First we write in the standard output
        u = open_file('out', 'outfile.thermal_conductivity')
        write (u, '(A2,A5,15X,A)') '# ', 'Unit:', 'W/m/K'
        write (u, '(A2,A12,8X,E20.12)') '# ', 'Temperature:', opts%temperature

        write (*, *) ''
        write (*, "(1X,A52)") 'Decomposition of the thermal conductivity (in W/m/K)'
        m0 = kappa_sma*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A)") 'Single mode approximation (SMA)'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A)") '# Single mode approximation'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)

        m0 = (kappa_iter - kappa_sma)*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A)") 'Correction to include collective contribution via iterative procedure'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A)") '# Collective contribution'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)

        m0 = kappa_offdiag*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A)") 'Off diagonal (coherence) contribution'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A)") '# Off diagonal (coherence) contribution'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)

        m0 = (kappa_iter + kappa_offdiag)*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A26)") 'Total thermal conductivity'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A28)") '# Total thermal conductivity'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)

        close (u)
        write(*, *) ''
    end if

    ! First we get the cumulative kappa with the mean free path
    if (mw%talk) write(*, *) '... computing cumulative kappa'
    call mf%get_cumulative_kappa(qp, dr, uc, opts%mfppts, opts%temperature, opts%sigma, mw, mem)
    call tmr_kappa%tock('cumulative kappa')

    ! Then the spectral kappa
    if (mw%talk) write(*, *) '... computing spectral kappa'
    call pd%generate(dr, qp, uc, mw, mem, verbosity=opts%verbosity, &
                     sigma=opts%sigma, n_dos_point=opts%freqpts, integrationtype=opts%dosintegrationtype)
    call mf%get_spectral_kappa(uc, qp, dr, pd, mw, mem)
    call tmr_kappa%tock('spectral kappa')

    ! In a last step, we have to add a prefactor to Fn, to have the right kappa per mode in the outfile
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
            dr%iq(q1)%Fn(:, b1) = dr%iq(q1)%Fn(:, b1)*dr%iq(q1)%omega(b1)/opts%temperature
        end do
    end do

    call tmr_tot%tock('thermal conductivity computation')
    call tmr_kappa%stop()
    t0 = walltime() - t0
end block blockkappa

finalize_and_write: block
    if (mw%talk) then
        write (*, *) ''
        write (*, *) '... dumping auxiliary data to files'
        call dr%write_to_hdf5(qp, uc, 'outfile.thermal_conductivity_grid.hdf5', mem, opts%temperature)
        call mf%write_to_hdf5(pd, uc, opts%unit, 'outfile.cumulative_thermal_conductivity.hdf5', mem)

        write (*, *) ''
        write (*, '(A,A)') 'Scattering rates can be found in                                ', 'outfile.thermal_conductivity_grid.hdf5'
        write (*, '(A,A)') 'Thermal conductivity tensor can be found in                     ', 'outfile.thermal_conductivity'
        write (*, '(A,A)') 'Cumulative and spectral thermal conductivity can be found in    ', 'outfile.cumulative_thermal_conductivity.hdf5'

        ! Print timings
        write (*, *) ''
        write (*, '(1X,A)') 'SUGGESTED CITATIONS:'
        write (*, '(1X,A41,A)') 'Software: ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'
        write (*, '(1X,A41,A)') 'Method: ', 'D. A. Broido et al., Appl Phys Lett 91, 231922 (2007)'
        write (*, '(1X,A41,A)') 'Iterative Boltzmann transport equation: ', 'M. Omini et al., Phys Rev B 53, 9064 (1996)'
        write (*, '(1X,A41,A)') 'Algorithm: ', 'A. H. Romero et al., Phys Rev B 91, 214310 (2015)'
        write (*, '(1X,A41,A)') 'Off-diagonal (coherences) contribution: ', 'M. Simoncelli et al., Nat Phys 15 809-813  (2019)'
        write (*, '(1X,A41,A)') '                                         ', 'L. Isaeva et al., Nat Commun 10 3853 (2019)'
        write (*, '(1X,A41,A)') '                                         ', 'A. Fiorentino et al., Phys Rev B 107, 054311  (2023)'
        write (*, '(1X,A41,A)') 'Theory : ', 'A. Castellano et al, J. Chem. Phys. 159 (23) (2023)'
        write (*, '(1X,A41,A)') 'Theory and algorithm : ', 'A. Castellano et al, ArXiv:2411.14949 (2024)'
    end if
    call tmr_tot%tock('io')

    call tmr_tot%stop()
    if (mw%talk) write (*, *) ''
    call tmr_init%dump(mw, 'Initialization timings:')
    call tmr_scat%dump(mw, 'Scattering timings:')
    call tmr_kappa%dump(mw, 'Thermal conductivity timings:')
    call tmr_tot%dump(mw, 'Total timings:')
end block finalize_and_write

! And we are done!
call sr%destroy()
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)
end program
