program lineshape
!!{!src/lineshape/manual.md!}
use konstanter, only: r8, lo_exitcode_param, lo_pi
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use lo_dielectric_interaction, only: lo_dielectric_tensor
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh, lo_read_qmesh_from_file, lo_get_small_group_of_qpoint
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use gottochblandat, only: open_file, walltime, lo_chop, lo_points_on_sphere
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer

use phonondamping
use dielscatter, only: lo_dielectric_response
use options, only: lo_opts
use io, only: write_lineshape_to_hdf5
use lo_realspace_selfenergy, only: lo_interpolated_selfenergy
use lo_thermal_transport, only: lo_thermal_conductivity
use lineshape_helper, only: lo_spectralfunction_helper

implicit none
type(lo_opts) :: opts
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_timer) :: tmr_init, tmr_calc, tmr_diel

type(lo_crystalstructure) :: uc
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_forceconstant_fourthorder) :: fcf
type(lo_interpolated_selfenergy) :: ise
type(lo_dielectric_tensor) :: di
type(lo_phonon_dispersions) :: dr, drd
type(lo_phonon_dos) :: pd
class(lo_qpoint_mesh), allocatable :: qp, qpd

real(r8) :: timer_init, timer_total

timer_total = walltime()
timer_init = walltime()

! Read information from file and work out the heuristics.
init: block
    ! Init MPI!
    call mw%init()
    ! Start the initialization timer
    call tmr_init%start()
    ! some options
    call opts%parse()
    ! only be verbose on the first rank
    if (.not. mw%talk) opts%verbosity = -100

    ! Read structure
    call uc%readfromfile('infile.ucposcar')
    call uc%classify('wedge', timereversal=.true.)
    if (mw%talk) write (*, *) '... using ', tochar(mw%n), ' MPI ranks'
    if (mw%talk) write (*, *) '... read structure'
    if (opts%readiso) then
        if (mw%talk) write (*, *) '... reading isotope distribution from file'
        call uc%readisotopefromfile()
    end if

    ! Read forceconstants
    call fc%readfromfile(uc, 'infile.forceconstant', mem, opts%verbosity)
    if (mw%talk) write (*, *) '... read second order forceconstant'
    if (opts%thirdorder) then
        call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
        if (mw%talk) write (*, *) '... read third order forceconstant'
    end if
    if (opts%fourthorder) then
        call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
        if (mw%talk) write (*, *) '... read fourth order forceconstant'
    end if
    if (opts%dielectric) then
        call di%readfromfile(uc, 'infile.dielectric_interaction')
        if (mw%talk) write (*, *) '... read dielectric interactions'
    end if
    if (opts%fancyinterp .or. opts%integrationtype .eq. 4) then
        call ise%read_from_hdf5('infile.interpolated_selfenergy.hdf5', mw)
        if (abs(ise%temperature - opts%temperature) .gt. lo_tol) then
            call lo_stop_gracefully(['Interpolated selfenergy is for a different temperature than input temperature.'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
        if (mw%talk) write (*, *) '... read fancy interpolation'
    end if
    call tmr_init%tock('read input files')

    ! Get a q-point mesh
    if (opts%readqmesh) then
        call lo_read_qmesh_from_file(qp, uc, 'infile.qgrid.hdf5', mem, verbosity=opts%verbosity)
    else
        select case (opts%meshtype)
        case (1)
            call lo_generate_qmesh(qp, uc, opts%qgrid, 'monkhorst', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity, nosym=opts%nosym)
        case (2)
            call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case (3)
            call lo_generate_qmesh(qp, uc, opts%qgrid, 'wedge', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case default
            call lo_stop_gracefully(['Unknown mesh type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end select
    end if

    call tmr_init%tock('generated q-mesh')

    ! Dispersions for everyone!
    call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)
    if (mw%talk) write (*, *) '... got the full dispersion relations'
    timer_init = walltime() - timer_init

    ! Now I can decide what the maximum frequency on the self-energy axis will be.
    ! This is quite a generous margin.
    opts%maxf = 2*(dr%omega_max*1.1 + maxval(dr%default_smearing)*3)

    call tmr_init%tock('harmonic dispersions')
    call tmr_init%stop()

    !if ( opts%integrationtype .eq. 4 ) then
    ! With this integration type I need to obtain spectral functions somehow.
    !call sfb%generate(qp,dr,opts%nf,opts%maxf,.false.,mw,mem,opts%verbosity+1)
    !endif

    call tmr_init%dump(mw, 'Initialization timings:')
end block init

! Now everything is read and set, time to calculate something. Three main modes: single point,
! bandstructure or dos. Only do one of them, gets annoying and confusing otherwise, I think.
if (opts%oneqpoint) then
    onepoint: block
        type(lo_phonon_dispersions_qpoint) :: wp
        type(lo_qpoint) :: qpoint
        type(lo_phonon_selfenergy) :: se
        type(lo_spectralfunction_helper) :: isf
        type(lo_dielectric_response) :: dir
        real(r8), dimension(3) :: v0, qin
        real(r8) :: f0

        call tmr_calc%start()

        if (mw%talk) then
            write (*, *) ''
            write (*, *) 'Calculating lineshapes for a single q-point'
            opts%verbosity = opts%verbosity + 1
        end if

        ! coordinates for this single point
        if (trim(opts%highsymmetrypoint) .ne. 'none') then
            v0 = uc%coordinate_from_high_symmetry_point_label(opts%highsymmetrypoint)
        else
            v0 = uc%fractional_to_cartesian(opts%qpoint, reciprocal=.true.)
        end if
        qpoint%r = lo_chop(v0, lo_sqtol)
        call lo_get_small_group_of_qpoint(qpoint, uc)

        ! From what direction are we approaching Gamma?
        qin = opts%q_in
        f0 = norm2(qin)
        if (f0 .lt. lo_sqtol) then
            call lo_stop_gracefully(['Q-direction needs to be specified!'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        else
            ! Define the Cartesian probing q-direction
            qin = qin/f0
        end if

        ! harmonic properties at this point
        call wp%generate(fc, uc, mem, qpoint, qdirection=qin)

        call tmr_calc%tock('harmonic properties at q')

        ! If applicable get spectralfunctions from file
        select case (opts%integrationtype)
        case (5)
            ! Read from file
            call isf%read_from_hdf5('infile.grid_spectral_function.hdf5', qp, opts%maxf, opts%temperature, opts%nf, mw, opts%verbosity)
            ! Smear the self-energy
            !call isf%smear(qp,dr,opts%maxf,opts%sigma,mw,mem)

            call tmr_calc%tock('reading and smearing input spectral function')
        end select

        ! get the self-energy
        if (opts%fancyinterp) then
            call se%generate_interp(qpoint, qin, wp, uc, ise, opts, mw, mem)
        else
            call se%generate(qpoint, qin, wp, uc, fc, fct, fcf, ise, isf, qp, dr, opts, tmr_calc, mw, mem, opts%verbosity)
        end if

        call tmr_calc%stop()
        call tmr_calc%dump(mw, 'Phonon selfenergy timings:')

        ! it could be the case that I have higher order Born charges and dielectric tensors!
        if (opts%dielectric) then
            call tmr_diel%start()
            call dir%generate(wp, di, uc, qp, dr, fc, fct, se, isf, opts, tmr_diel, mw, mem, opts%verbosity)
            call tmr_diel%stop()
            call tmr_diel%dump(mw, 'Dielectric response timings:')
        end if

        ! dump it to file
        if (mw%talk) then
            if (opts%dielectric) then
                call write_lineshape_to_hdf5(se, uc, qpoint, wp, dir, di, qp, opts%enhet, opts%temperature, 'outfile.dielectric_response.hdf5', mem, opts%verbosity, dielectric=.true.)
            else
                call write_lineshape_to_hdf5(se, uc, qpoint, wp, dir, di, qp, opts%enhet, opts%temperature, 'outfile.phonon_self_energy.hdf5', mem, opts%verbosity, dielectric=.false.)
            end if
        end if
    end block onepoint
end if

if (opts%qpointpath) then
    qppath: block
        type(lo_phonon_bandstructure) :: bs
        integer :: npts, i

        ! Start the calculation timer
        call tmr_calc%start()

        if (mw%talk) then
            write (*, *) ''
            write (*, *) 'Calculating lineshapes on a path in the BZ.'
        end if

        ! Generate the path, and the harmonic things on said path
        npts = 0
        do i = 1, opts%nq_on_path
            npts = npts + opts%stride
            if (npts .ge. opts%nq_on_path - opts%stride) exit
        end do
        npts = npts + 1
        call bs%generate(uc, fc, timereversal=opts%timereversal, mw=mw, mem=mem, verbosity=opts%verbosity, npts=npts, readpathfromfile=opts%readpathfromfile)

        ! Calculate the actual thingy. Quite long, so moved to it's own routine.
        if (opts%fancyinterp) then
            call spectral_function_along_path_interp(bs, uc, ise, opts, mw, mem)
        else
            call spectral_function_along_path(bs, uc, fc, fct, fcf, ise, qp, dr, opts, tmr_calc, mw, mem)
        end if

        ! Dump this
        if (mw%talk) then
            write (*, *) '... writing output'
            call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations.hdf5', mem)
            call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function.hdf5')
        end if

        ! Report timings
        call tmr_calc%stop()
        call tmr_calc%dump(mw, 'Timings:')
    end block qppath
end if

! Calculate spectral function on a (closed) grid
if (opts%grid) then
    makegrid: block
        type(lo_spectralfunction_helper) :: sf
        type(lo_thermal_conductivity) :: tc

        ! Start the calculation timer
        call tmr_calc%start()

        if (mw%talk) then
            write (*, *) ''
            write (*, *) 'Calculating the spectral functions on a grid'
        end if

        if (opts%meshtype .ne. 2) then
            call lo_stop_gracefully(['Grid calculations need FFT meshes'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if

        ! Calculate the self-energy on a closed grid
        call get_selfenergy_on_closed_grid(sf, tc, pd, qp, dr, uc, fc, fct, fcf, ise, opts, tmr_calc, mw, mem, opts%verbosity + 1)

        call tmr_calc%tick()
        ! Store the calculated spectral function to file
        if (mw%r .eq. mw%n - 1) then
            call sf%write_to_hdf5('outfile.grid_spectral_function.hdf5')
            ! Destroy sf here!
        end if

        ! Dump thermal transport to file
        call tc%write_to_hdf5(qp, dr, uc, 'outfile.thermal_conductivity.hdf5', opts%enhet, mw, mem)

        ! Dump phonon dos to file
        if (mw%talk) then
            call pd%write_to_hdf5(uc, opts%enhet, 'outfile.phonon_spectralfunction_dos.hdf5', mem)
        end if

        call tmr_calc%tock('io')

        ! Report timings
        call tmr_calc%stop()
        call tmr_calc%dump(mw, 'Timings:')
    end block makegrid
end if

if (opts%phonondos) then
    makedos: block
        type(lo_spectralfunction_helper) :: sf
        type(lo_thermal_conductivity) :: tc

        call tmr_calc%start()

        if (mw%talk) write (*, *) ''
        if (mw%talk) write (*, *) 'Calculating the phonon dos'
        select case (opts%meshtype)
        case (1)
            call lo_generate_qmesh(qpd, uc, opts%qgrid_dos, 'monkhorst', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case (2)
            call lo_generate_qmesh(qpd, uc, opts%qgrid_dos, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case (3)
            call lo_generate_qmesh(qpd, uc, opts%qgrid_dos, 'wedge', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case default
            call lo_stop_gracefully(['Unknown mesh type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end select

        ! and dispersions on this grid
        call drd%generate(qpd, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

        ! Adjust the max frequency, should be with respect to the tight grid!
        opts%maxf = 2*(drd%omega_max*1.1 + maxval(drd%default_smearing)*3)

        ! Get the actual dos
        if (opts%fancyinterp) then
            call get_intensity_as_dos_interp(pd, tc, qpd, drd, uc, ise, opts, mw, mem, opts%verbosity + 1)
        else
            call get_intensity_as_dos(pd, qpd, drd, uc, fc, fct, fcf, ise, sf, qp, dr, opts, tmr_calc, mw, mem, opts%verbosity + 1)
        end if

        ! Store the calculated spectral function to file
        if (mw%r .eq. mw%n - 1) then
            call sf%write_to_hdf5('outfile.grid_spectral_function.hdf5')
            ! Destroy sf here!
        end if

        ! Massage and dump thermal transport stuff to file.
        ! if ( mw%talk ) then
        !     write(*,*) '... writing thermal conductivity to file'
        ! endif
        ! call tc%write_to_hdf5(qpd,uc,'outfile.spectral_thermal_conductivity.hdf5',opts%enhet,mw,mem)

        ! Dump the dos to file
        if (mw%talk) then
            ! Dump thermal conductivity here

            call pd%write_to_file(uc, opts%enhet, 'outfile.phonon_spectralfunction_dos')
            call pd%write_to_hdf5(uc, opts%enhet, 'outfile.phonon_spectralfunction_dos.hdf5', mem)
        end if

        ! Dump all the data on the dos grid to file
        if (mw%talk .and. opts%dumpgrid) then
            call drd%write_to_hdf5(qpd, uc, 'outfile.grid_dispersions.hdf5', mem)
        end if

        call tmr_calc%stop()
        call tmr_calc%dump(mw, 'Timings:')
    end block makedos
end if

! Generate an interpolation thing. Fancy stuff.
if (opts%geninterp) then
    buildinterp: block
        real(r8), parameter :: maxrad = 0.25_r8
        integer, parameter :: nextra = 20
        type(lo_phonon_dispersions_qpoint), dimension(:), allocatable :: wp
        real(r8), dimension(:, :, :), allocatable :: sigmaRe, sigmaIm
        real(r8), dimension(:, :), allocatable :: qvecs
        real(r8), dimension(3, nextra) :: extraq
        real(r8) :: f0
        integer :: i, nq

        call tmr_calc%start()

        if (mw%talk) write (*, *) ''
        if (mw%talk) write (*, *) 'Building magical interpolation'
        select case (2)
        case (1)
            call lo_generate_qmesh(qpd, uc, opts%qgrid_dos, 'monkhorst', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case (2)
            call lo_generate_qmesh(qpd, uc, opts%qgrid_dos, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case (3)
            call lo_generate_qmesh(qpd, uc, opts%qgrid_dos, 'wedge', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        case default
            call lo_stop_gracefully(['Unknown mesh type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end select

        ! Number of q-points
        nq = nextra + qpd%n_irr_point
        call mem%allocate(qvecs, [3, nq], persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        qvecs = 0.0_r8

        ! Get the list of q-points. First the odd ones:
        call lo_points_on_sphere(extraq)
        ! Massage them a little to make them a spiral close to Gamma
        do i = 1, nextra
            f0 = real(i, r8)/real(nextra, r8)
            f0 = f0*uc%bz%rmin*maxrad
            extraq(:, i) = extraq(:, i)*f0
        end do
        ! Add them together
        do i = 1, qpd%n_irr_point
            qvecs(:, i) = qpd%ip(i)%r
        end do
        do i = 1, nextra
            qvecs(:, i + qpd%n_irr_point) = extraq(:, i)
        end do

        ! and dispersions on this grid
        !call drd%generate(qpd,fc,uc,mw=mw,mem=mem,verbosity=opts%verbosity)

        call mem%allocate(sigmaRe, [opts%nf, dr%n_mode, nq], persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(sigmaIm, [opts%nf, dr%n_mode, nq], persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)

        ! Get self-energies on the grid
        !call get_selfenergy_on_grid(pd,qpd,drd,uc,fc,fct,fcf,ise,qp,dr,opts,sigmaRe,sigmaIm,mw,mem,opts%verbosity+1)
        allocate (wp(nq))
        call get_selfenergy_on_points(qvecs, wp, qp, dr, uc, fc, fct, fcf, ise, opts, sigmaRe, sigmaIm, tmr_calc, mw, mem, opts%verbosity + 1)

        ! Build interpolated model?
        call ise%generate(uc, dr, qvecs, wp, opts%maxf, sigmaRe, sigmaIm, opts%cutoff, opts%temperature, mw, mem, opts%verbosity + 1)
        !call ise%generate(uc,qpd,drd,opts%maxf,sigmaRe,sigmaIm,opts%cutoff,opts%temperature,mw,mem,opts%verbosity+1)
        call mem%deallocate(qvecs, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(sigmaRe, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(sigmaIm, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        deallocate (wp)

        call tmr_calc%stop()
        call tmr_calc%dump(mw, 'Timings:')

        ! Dump to file
        if (mw%talk) then
            call ise%write_to_hdf5('outfile.interpolated_selfenergy.hdf5')
        end if
    end block buildinterp
end if

! All done, print timings
if (mw%talk) then
    timer_total = walltime() - timer_total
    write (*, *) ' '
    write (*, *) ' Timings:'
    write (*, *) '              initialization:', timer_init
    write (*, *) '                       total:', timer_total
end if
! Kill MPI
call mw%destroy()

end program
