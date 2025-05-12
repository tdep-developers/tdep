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
use lineshape_helper, only: lo_spectralfunction_helper, evaluate_spectral_function
use scf_helper, only: return_new_bare_phonons
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
!type(lo_spectralfunction_helper) :: sf

type(lo_phonon_dispersions) :: dr,ddr

class(lo_qpoint_mesh), allocatable :: qp,dqp

real(r8) :: timer_init, timer_total
!integer :: iter

timer_total = walltime()
timer_init = walltime()

! How do I go about doing this? I guess the idea would be
! to first get the force constants, calculate the self-energy
! in mode space, then to dynamical matrices, then fit to new
! second order and so on.

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

    ! Get a q-point mesh. We only want FFT meshes.
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
    call lo_generate_qmesh(dqp, uc, opts%qgrid_sigma, 'fft', timereversal=opts%timereversal, headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)

    if ( mw%talk ) then
        write(*,*) 'q-mesh the self-energy is evaluated on:',opts%qgrid_sigma
    endif

    ! And the initial harmonic dispersions
    call dr%generate(qp, fc2, uc, mw=mw, mem=mem, verbosity=opts%verbosity)
    call ddr%generate(dqp, fc2, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

    ! Now I can decide what the maximum frequency on the self-energy axis will be.
    ! This is quite a generous margin.
    opts%maxf = 2*(dr%omega_max*1.1_r8 + maxval(dr%default_smearing)*3)

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

    ! For diagnostics I guess dumping the self-energy on a path makes sense? For that we
    ! first need the perfectly normal path for reference.
    call bs%generate(uc, fc2, timereversal=.true., mw=mw, mem=mem, verbosity=opts%verbosity, npts=100, readpathfromfile=.false.)

    ! Make sure the interpolated self-energy is nothing
    call ise%destroy()
    ! Read it from file?
    call ise%read_from_hdf5(uc,'outfile.interpolated_selfenergy.hdf5',mw,mem,opts%verbosity+1)
    ! Get spectral function on a path?
    call ise%spectral_function_along_path(bs,uc,mw,mem)

    ! Maybe get the spectral function on a grid? A third q-grid you say? Why not.
    call ise%spectral_function_on_grid(uc,fc2,[12,12,12],opts%sigma,opts%temperature,tc,pd,mw,mem)
    if ( mw%talk ) then
        call pd%write_to_hdf5(uc,opts%enhet,'outfile.phonon_dos.hdf5',mem)
    endif

    if (mw%talk) then
        write (*, *) '... writing output'
        call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations.hdf5', mem)
        call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function_0.hdf5')
    end if

    ! do iter=1,3
    !     ! Then I guess the next step is to get the self-energy again, but this time using
    !     ! a convolution integration instead?
    !     call generate_interpolated_selfenergy('outfile.interpolated_selfenergy.hdf5',uc,fc2,fc3,fc4,ise,qp,dqp,dr,ddr, &
    !         opts%temperature, opts%maxf, opts%nf, 4, opts%sigma,&
    !         opts%isotopescattering, opts%thirdorder, opts%fourthorder, .false.,&
    !         mw, mem, opts%verbosity)

    !     ! Make sure the interpolated self-energy is nothing
    !     call ise%destroy()
    !     ! Read it from file?
    !     call ise%read_from_hdf5(uc,'outfile.interpolated_selfenergy.hdf5',mw,mem,opts%verbosity+1)
    !     ! Get spectral function on a path?
    !     call bs%destroy()
    !     call bs%generate(uc, fc2, timereversal=.true., mw=mw, mem=mem, verbosity=opts%verbosity, npts=100, readpathfromfile=.false.)
    !     call ise%spectral_function_along_path(bs,uc,mw,mem)

    !     if (mw%talk) then
    !         write (*, *) '... writing output'
    !         call bs%write_to_hdf5(uc, opts%enhet, 'outfile.dispersion_relations.hdf5', mem)
    !         call bs%write_spectral_function_to_hdf5(opts%enhet, 'outfile.phonon_spectral_function_'//tochar(iter)//'.hdf5')
    !     end if
    ! enddo


end block firstiteration

! outerloop: do iter=1,-1

!     newbaseline: block
!         real(r8), dimension(:), allocatable :: x2,x3,x4
!         integer :: i

!         if ( map%xuc%nx_fc_pair .gt. 0 ) then
!             allocate(x2(map%xuc%nx_fc_pair))
!             x2=map%xuc%x_fc_pair
!         endif
!         if ( map%xuc%nx_fc_triplet .gt. 0 ) then
!             allocate(x3(map%xuc%nx_fc_triplet))
!             x3=map%xuc%x_fc_triplet
!         endif
!         if ( map%xuc%nx_fc_quartet .gt. 0 ) then
!             allocate(x4(map%xuc%nx_fc_quartet))
!             x4=map%xuc%x_fc_quartet
!         endif

!         !call return_new_bare_phonons(qp,dr,sf,map,uc,fc2,mw,mem,opts%verbosity)
!         ! Here is a good place for some mixing, I think.
!         !map%xuc%x_fc_pair = map%xuc%x_fc_pair*0.5_r8 + 0.5_r8*x2
!         !map%xuc%x_fc_pair = x2

!         !call lo_solve_for_irreducible_ifc_fastugly(map,uc,ss,sim,mw,mem,opts%verbosity,fix_secondorder=.true.)
!         !call lo_solve_for_irreducible_ifc_fastugly(map,uc,ss,sim,mw,mem,opts%verbosity,fix_secondorder=.false.)

!         if ( mw%talk ) then
!             write(*,*) 'Iteration',iter
!             write(*,*) 'irreducible pairs'
!             if ( map%xuc%nx_fc_pair .gt. 0 ) then
!                 do i=1,map%xuc%nx_fc_pair
!                     write(*,*) i,x2(i),map%xuc%x_fc_pair(i),abs( 100*(x2(i)-map%xuc%x_fc_pair(i))/x2(i) )
!                 enddo
!             endif
!             if ( map%xuc%nx_fc_triplet .gt. 0 ) then
!                 write(*,*) 'irreducible triplet'
!                 do i=1,map%xuc%nx_fc_triplet
!                     write(*,*) i,x3(i),map%xuc%x_fc_triplet(i),abs( 100*(x3(i)-map%xuc%x_fc_triplet(i))/x3(i) )
!                 enddo
!             endif
!             if ( map%xuc%nx_fc_quartet .gt. 0 ) then
!                 write(*,*) 'irreducible quartet'
!                 do i=1,map%xuc%nx_fc_quartet
!                     write(*,*) i,x4(i),map%xuc%x_fc_quartet(i),abs( 100*(x4(i)-map%xuc%x_fc_quartet(i))/x4(i) )
!                 enddo
!             endif
!         endif

!         ! Neat and efficient I think. Now what. Hmmm. Suppose I should
!         ! fetch new forceconstants.
!         if (map%have_fc_pair) then
!             call map%get_secondorder_forceconstant(uc, fc2, mem, opts%verbosity)
!             !if (mw%talk) call fc2%writetofile(uc, 'outfile.forceconstant')
!         end if
!         if (map%have_fc_triplet) then
!             call map%get_thirdorder_forceconstant(uc, fc3)
!             !if (mw%talk) call fc3%writetofile(uc, 'outfile.forceconstant_thirdorder')
!         end if
!         if (map%have_fc_quartet) then
!             call map%get_fourthorder_forceconstant(uc, fc4)
!             !if (mw%talk) call fc4%writetofile(uc, 'outfile.forceconstant_fourthorder')
!         end if

!         ! And we update the bare phonons.
!         call dr%generate(qp, fc2, uc, mw=mw, mem=mem, verbosity=opts%verbosity)
!     end block newbaseline

!     ! And then we need a new self-energy.
!     newselfenergy: block
!         type(lo_thermal_conductivity) :: tc
!         ! Calculate the self-energy on a closed grid
!         opts%integrationtype = 5
!         call get_selfenergy_on_closed_grid(sf, tc, pd, qp, dr, uc, fc2, fc3, fc4, ise, opts, tmr_calc, mw, mem, opts%verbosity + 1)

!         if (mw%r .eq. mw%n - 1) then
!             call sf%write_to_hdf5('infile.grid_spectral_function.hdf5')
!         end if

!     end block newselfenergy

! enddo outerloop


! All done, print timings
if (mw%talk) then
    write (*, *) ' '
    write (*, *) 'All done! '
end if
! Kill MPI
call mw%destroy()

end program
