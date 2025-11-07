program project_happy_mole
use konstanter, only: r8, lo_exitcode_param, lo_pi, lo_freqtol

use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder

use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh, lo_read_qmesh_from_file, lo_get_small_group_of_qpoint
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use gottochblandat, only: open_file, walltime, lo_chop, lo_points_on_sphere, lo_does_file_exist, tochar, lo_trapezoid_integration, lo_lorentz
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer

!use dielscatter, only: lo_dielectric_response
use options, only: lo_opts
use lo_thermal_transport, only: lo_thermal_conductivity
use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lo_evaluate_phonon_self_energy, only: lo_phonon_selfenergy
use create_selfenergy_interpolation, only: generate_interpolated_selfenergy

implicit none
type(lo_opts) :: opts
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_timer) :: tmr_init !, tmr_calc, tmr_diel

type(lo_crystalstructure) :: uc !,ss
type(lo_forceconstant_secondorder) :: fc2
type(lo_forceconstant_thirdorder) :: fc3
type(lo_forceconstant_fourthorder) :: fc4

type(lo_phonon_dispersions) :: dr,ddr,kdr

class(lo_qpoint_mesh), allocatable :: qp,dqp,kqp

! Read information from file and work out the heuristics.
init: block
    !real(r8) :: t0

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

    call tmr_init%tock('read structures')

    call fc2%readfromfile(uc, 'infile.forceconstant', mem, opts%verbosity)
    if (mw%talk) write (*, *) '... read second order forceconstant'
    if (opts%thirdorder) then
        call fc3%readfromfile(uc, 'infile.forceconstant_thirdorder')
        if (mw%talk) write (*, *) '... read third order forceconstant'
    end if
    if (opts%fourthorder) then
        call fc4%readfromfile(uc, 'infile.forceconstant_fourthorder')
        if (mw%talk) write (*, *) '... read fourth order forceconstant'
    end if

    call tmr_init%tock('read forceconstants')

    ! Get a q-point meshes. We only want FFT meshes.
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
    call lo_generate_qmesh(dqp, uc, opts%qgrid_sigma, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
    call lo_generate_qmesh(kqp, uc, opts%qgrid_kappa, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)

    if ( mw%talk ) then
        write(*,*) '               q-mesh for phase space integrals: ',tochar(opts%qgrid)
        write(*,*) '         q-mesh the self-energy is evaluated on: ',tochar(opts%qgrid_sigma)
        write(*,*) 'q-mesh the thermal conductivity is evaluated on: ',tochar(opts%qgrid_kappa)
    endif

    ! And the initial harmonic dispersions
    call dr%generate(qp, fc2, uc, mw=mw, mem=mem, verbosity=opts%verbosity)
    call ddr%generate(dqp, fc2, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

    ! Now I can decide what the maximum frequency on the self-energy
    ! axis will be. This is quite a generous margin.
    opts%maxf = 3*(dr%omega_max*1.1_r8 + maxval(dr%default_smearing)*3)

    call tmr_init%tock('initial harmonic properties')

    call tmr_init%stop()
    call tmr_init%dump(mw, 'Initialization timings:')
end block init

! Think the first iteration will be special in many ways, so let's make
! that one its own thing.
firstiteration: block
    type(lo_interpolated_selfenergy_grid) :: ise
    type(lo_phonon_bandstructure) :: bs
    type(lo_phonon_dos) :: pd
    type(lo_thermal_conductivity) :: tc
    integer :: iter

    ! This is the zeroth iteration, or whatever I should call it. Here we always
    ! use adaptive Gaussian integration, because why not.
    call generate_interpolated_selfenergy('outfile.interpolated_selfenergy.hdf5',uc,fc2,fc3,fc4,ise,qp,dqp,dr,ddr, &
        opts%temperature, opts%maxf, opts%nf, 2, opts%sigma,&
        opts%isotopescattering, opts%thirdorder, opts%fourthorder, .false.,&
        mw, mem, opts%verbosity)

    ! For diagnostics I guess dumping the self-energy on a path makes sense?
    ! For that we first need the perfectly normal path for reference.
    call bs%generate(uc, fc2, timereversal=.true., mw=mw, mem=mem, verbosity=opts%verbosity, npts=100, readpathfromfile=.false.)

    ! Make sure the interpolated self-energy is nothing
    call ise%destroy()
    ! Read it from file?
    call ise%read_from_hdf5(uc,'outfile.interpolated_selfenergy.hdf5',mw,mem,opts%verbosity+1)
    call mw%barrier()

    ! Get spectral function on a path?
    if (mw%talk) then
        write(*,*) '... generating spectral function on path'
    endif
    call ise%spectral_function_along_path(bs,uc,mw,mem)

    ! Generate spectral function on a grid?
    if (mw%talk) then
        write(*,*) '... generating spectral function on a grid'
    endif
    call ise%spectral_function_on_grid(uc,fc2,kqp,opts%sigma,opts%temperature,tc,pd,kdr,mw,mem)

    if (mw%talk) then
        write (*, *) '... writing output'
        call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations_0.hdf5', mem)
        call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function_0.hdf5')
        call pd%write_to_hdf5(uc,opts%enhet,'outfile.spectral_function_dos_0.hdf5',mem)
    end if
    call tc%write_to_hdf5(kqp,kdr,uc,'outfile.thermal_conductivity_0.hdf5',opts%enhet,mw,mem)

    ! Then I guess we start to iterate, self-consistently?
    do iter=1,1
        ! Then I guess the next step is to get the self-energy again, but this time using
        ! a convolution integration instead?
        call generate_interpolated_selfenergy('outfile.interpolated_selfenergy.hdf5',uc,fc2,fc3,fc4,ise,qp,dqp,dr,ddr, &
            opts%temperature, opts%maxf, opts%nf, 4, opts%sigma,&
            opts%isotopescattering, opts%thirdorder, opts%fourthorder, .false.,&
            mw, mem, opts%verbosity)

        ! Make sure the intermediate things are cleaned:
        call ise%destroy()
        call kdr%destroy()
        call tc%destroy()
        call pd%destroy()
        ! Then we read the newly created interpolation and get a spectral function on a path?

        ! Read it from file?
        call ise%read_from_hdf5(uc,'outfile.interpolated_selfenergy.hdf5',mw,mem,opts%verbosity+1)
        ! Get spectral function on a path?
        call ise%spectral_function_along_path(bs,uc,mw,mem)
        ! Spectral function on a grid
        call ise%spectral_function_on_grid(uc,fc2,kqp,opts%sigma,opts%temperature,tc,pd,kdr,mw,mem)
        if (mw%talk) then
            write (*, *) '... writing output'
            call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations_'//tochar(iter)//'.hdf5', mem)
            call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function_'//tochar(iter)//'.hdf5')
            call pd%write_to_hdf5(uc,opts%enhet,'outfile.spectral_function_dos_'//tochar(iter)//'.hdf5',mem)
        end if
        call tc%write_to_hdf5(kqp,kdr,uc,'outfile.thermal_conductivity_'//tochar(iter)//'.hdf5',opts%enhet,mw,mem)
    enddo

end block firstiteration


! All done, print timings
if (mw%talk) then
    write (*, *) ' '
    write (*, *) 'All done! '
end if
! Kill MPI
call mw%destroy()

end program
