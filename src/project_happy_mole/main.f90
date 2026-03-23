program project_happy_mole
use konstanter, only: r8, lo_exitcode_param, lo_pi, lo_freqtol, lo_kappa_au_to_SI

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
use options, only: lo_opts

use lo_distributed_phonon_dispersion_relations, only: lo_distributed_phonon_dispersions
use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lo_evaluate_phonon_self_energy, only: lo_phonon_selfenergy
use create_selfenergy_interpolation, only: generate_interpolated_selfenergy
use lo_collision_matrix, only: lo_scattering_matrix
use bubble, only: bubble_only_transport

implicit none
type(lo_opts) :: opts
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_timer) :: tmr_init

type(lo_crystalstructure) :: uc
type(lo_forceconstant_secondorder) :: fc2
type(lo_forceconstant_thirdorder) :: fc3
type(lo_forceconstant_fourthorder) :: fc4

type(lo_phonon_dispersions) :: ddr
type(lo_distributed_phonon_dispersions) :: pdr

class(lo_qpoint_mesh), allocatable :: qp,dqp,kqp
type(lo_interpolated_selfenergy_grid) :: ise

! Read information from file and work out the heuristics.
init: block

    ! Init MPI!
    call mw%init()

    ! Init memory tracker
    call mem%init()

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

    ! Get q-point meshes. I will allow for the sigma mesh to not be FFT.
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
    if ( opts%readqmesh ) then
        call lo_read_qmesh_from_file(dqp, uc, 'infile.qgrid.hdf5', mem, verbosity=opts%verbosity)
    else
        call lo_generate_qmesh(dqp, uc, opts%qgrid_sigma, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
    endif
    call lo_generate_qmesh(kqp, uc, opts%qgrid_kappa, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)

    if ( mw%talk ) then
        write(*,*) '               q-mesh for phase space integrals: ',tochar(opts%qgrid)
        write(*,*) '         q-mesh the self-energy is evaluated on: ',tochar(opts%qgrid_sigma)
        write(*,*) 'q-mesh the thermal conductivity is evaluated on: ',tochar(opts%qgrid_kappa)
    endif

    ! And the initial harmonic dispersions
    call ddr%generate(dqp, fc2, uc, mw=mw, mem=mem, verbosity=opts%verbosity)
    call pdr%generate(qp, uc, fc2, opts%sigma, mw=mw, mem=mem, verbosity=opts%verbosity)

    ! Now I can decide the maximum frequency on the self-energy
    ! axis. This is quite a generous margin.
    opts%maxf = 3*(pdr%omega_max*1.1_r8 + maxval(pdr%default_smearing)*3)

    call tmr_init%tock('initial harmonic properties')

    call tmr_init%stop()
    call tmr_init%dump(mw, 'Initialization timings:')
end block init

! Think the first iteration will be special in many ways, so let's make
! that one its own thing.
if ( opts%readselfenergy) then
    readselfenergy: block
        call ise%destroy()
        call ise%read_from_hdf5(uc,fc2,'infile.interpolated_selfenergy.hdf5',mw,mem,opts%verbosity+1)
        call mw%barrier()
    end block readselfenergy
else
    calculateselfenergy: block
        integer, parameter :: max_n_iter=10
        type(lo_timer) :: tmr_sigma
        type(lo_phonon_bandstructure) :: bs
        type(lo_phonon_dos) :: pd
        real(r8), dimension(:,:,:), allocatable :: conv_kappa
        real(r8), dimension(:,:), allocatable :: conv_dos
        real(r8) :: f0,f1
        integer :: iter,i

        call tmr_sigma%start()

        ! Some space for convergence monitoring
        allocate(conv_kappa(3,3,max_n_iter))
        allocate(conv_dos(opts%nf,max_n_iter))
        conv_kappa=0.0_r8
        conv_dos=0.0_r8

        ! This is the zeroth iteration, or whatever I should call it. Here we
        ! always use adaptive Gaussian integration, because we have to use something.
        call generate_interpolated_selfenergy('outfile.interpolated_selfenergy.hdf5',uc,fc2,fc3,fc4,ise,qp,dqp,ddr,pdr, &
            opts%temperature, opts%maxf, opts%nf, 2, opts%sigma,&
            opts%isotopescattering, opts%thirdorder, opts%fourthorder, &
            mw, mem, opts%verbosity)
        call mw%barrier()
        call tmr_sigma%tock('integrate self-energy')

        ! For diagnostics I guess dumping the self-energy on a path makes sense?
        ! For that we first need the perfectly normal path for reference.
        call bs%generate(uc, fc2, timereversal=.true., mw=mw, mem=mem, verbosity=opts%verbosity, npts=100, readpathfromfile=.false.)

        ! Make sure the interpolated self-energy is nothing
        call ise%destroy()
        call ise%read_from_hdf5(uc,fc2,'outfile.interpolated_selfenergy.hdf5',mw,mem,opts%verbosity+1)
        call mw%barrier()
        call tmr_sigma%tock('initialize interpolation')

        ! Get spectral function on a path, for diagnostics
        if (mw%talk) then
            write(*,*) '... generating spectral function on path'
        endif
        call ise%spectral_function_along_path(bs,uc,mw,mem)
        call mw%barrier()
        call tmr_sigma%tock('interpolate to path')

        ! Generate spectral function on a grid?
        if (mw%talk) then
           write(*,*) '... generating spectral function on a grid'
        endif
        call ise%spectral_function_on_grid_rough(uc,fc2,kqp,opts%sigma,opts%temperature,pd,conv_kappa(:,:,1),mw,mem)
        conv_dos(:,1)=pd%dos

        call mw%barrier()
        call tmr_sigma%tock('interpolate to grid')

        if (mw%talk) then
            write (*, *) '... writing output'
            call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations_0.hdf5', mem)
            call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function_0.hdf5')
            call pd%write_to_hdf5(uc,opts%enhet,'outfile.spectral_function_dos_0.hdf5',mem)
        end if
        call mw%barrier()
        call tmr_sigma%tock('io')

        ! Then I guess we start to iterate, self-consistently?
        iterloop: do iter=1,max_n_iter-1
            ! Then I guess the next step is to get the self-energy again, but this time using
            ! a convolution integration instead?
            call generate_interpolated_selfenergy('outfile.interpolated_selfenergy.hdf5',uc,fc2,fc3,fc4,ise,qp,dqp,ddr,pdr, &
                opts%temperature, opts%maxf, opts%nf, 4, opts%sigma,&
                opts%isotopescattering, opts%thirdorder, opts%fourthorder, &
                mw, mem, opts%verbosity)
            call mw%barrier()
            call tmr_sigma%tock('integrate self-energy')


            ! Make sure the intermediate things are cleaned:
            call ise%destroy()
            call pd%destroy()

            ! Read the newly created spectral function from file
            call ise%destroy()
            call ise%read_from_hdf5(uc,fc2,'outfile.interpolated_selfenergy.hdf5',mw,mem,opts%verbosity+1)
            call mw%barrier()
            call tmr_sigma%tock('initialize interpolation')

            ! Get spectral function on a path?
            call ise%spectral_function_along_path(bs,uc,mw,mem)
            call mw%barrier()
            call tmr_sigma%tock('interpolate to path')

            ! Spectral function on a grid
            call ise%spectral_function_on_grid_rough(uc,fc2,kqp,opts%sigma,opts%temperature,pd,conv_kappa(:,:,iter+1),mw,mem)
            conv_dos(:,iter+1)=pd%dos
            call mw%barrier()
            call tmr_sigma%tock('interpolate to grid')

            if (mw%talk) then
                write (*, *) '... writing output'
                call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations_'//tochar(iter)//'.hdf5', mem)
                call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function_'//tochar(iter)//'.hdf5')
                call pd%write_to_hdf5(uc,opts%enhet,'outfile.spectral_function_dos_'//tochar(iter)//'.hdf5',mem)
            end if
            call mw%barrier()
            call tmr_sigma%tock('io')

            ! Here we should perhaps check for convergence. I wonder what to check.
            if ( mw%talk ) then
                write(*,*) 'Montiring convergence:'
                do i=1,iter+1
                    if ( i .gt. 1 ) then
                        f0=sum(abs(conv_kappa(:,:,i)-conv_kappa(:,:,i-1)))/sum(abs(conv_kappa(:,:,i)))
                        f1=sum(abs(conv_dos(:,i)-conv_dos(:,i-1)))/sum(abs(conv_dos(:,i)))
                        write(*,*) i,conv_kappa(1,1,i)*lo_kappa_au_to_SI,f0,f1
                    else
                        write(*,*) i,conv_kappa(1,1,i)*lo_kappa_au_to_SI
                    endif
                enddo
            endif

            ! What is a sensible criteria?

            f0=sum(abs(conv_kappa(:,:,iter+1)-conv_kappa(:,:,iter)))/sum(abs(conv_kappa(:,:,iter+1)))
            f1=sum(abs(conv_dos(:,iter+1)-conv_dos(:,iter)))/sum(abs(conv_dos(:,iter+1)))
            if ( f0 + f1 .lt. 1E-6_r8 ) then
                if ( mw%talk ) write(*,*) 'This seems converged!'
                exit iterloop
            endif

        enddo iterloop

        call tmr_sigma%stop()
        call tmr_sigma%dump(mw,"Self-consistent Green's function timings")

    end block calculateselfenergy
endif

postselfenergy: block

        ! ! Generate a bubble-only thermal transport (for now)
        ! if (mw%talk) then
        !     write(*,*) '... generating bubble-only thermal transport'
        ! endif
        ! call bubble_only_transport(kqp,psdr,uc,ise,opts%sigma,opts%temperature,mw,mem,opts%verbosity)

        ! if ( mw%talk ) write(*,*) 'done here ',__FILE__,__LINE__
        ! call mw%destroy()
        ! stop


        ! Create scattering matrix?
        !call scm%generate(uc,fc2,fc3,ise,kqp,mw,mem,opts%verbosity+5)

end block postselfenergy

! All done, print timings
if (mw%talk) then
    write (*, *) ' '
    write (*, *) 'All done! '
end if
! Kill MPI
call mw%destroy()

end program
