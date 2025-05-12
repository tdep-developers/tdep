module create_selfenergy_interpolation
!!
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_exitcode_symmetry, lo_twopi, &
    lo_bohr_to_A, lo_freqtol, lo_exitcode_param, lo_pi, lo_imag, lo_groupvel_ms_to_Hartreebohr
use gottochblandat, only: tochar, walltime, lo_progressbar_init, lo_progressbar, lo_chop, &
        lo_linear_least_squares, lo_rsquare, lo_linspace, lo_linear_interpolation, lo_trapezoid_integration,&
        lo_clean_fractional_coordinates, lo_cross
use geometryfunctions, only: lo_inscribed_sphere_in_box, lo_plane
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_wedge_mesh, lo_qpoint, lo_read_qmesh_from_file
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm, lo_dgels

use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lineshape_helper, only: evaluate_spectral_function,taperfn_im
use lo_evaluate_phonon_self_energy, only: lo_phonon_selfenergy

implicit none

private
public :: generate_interpolated_selfenergy

contains

!> create the interpolated self-energy object and store it to file
subroutine generate_interpolated_selfenergy(filename,uc,fc,fct,fcf,ise,qp,dqp,dr,ddr, &
    temperature, max_energy, n_energy, integrationtype, sigma,&
    use_isotope, use_thirdorder, use_fourthorder, offdiagonal, &
    mw, mem, verbosity)
    !> filename
    character(len=*), intent(in) :: filename
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> tabulated and smeared spectral functions
    type(lo_interpolated_selfenergy_grid), intent(in) :: ise
    !> q-point mesh for self-energy integrations
    class(lo_qpoint_mesh), intent(in) :: qp
    !> q-point mesh the self-energy is evaluated on
    class(lo_qpoint_mesh), intent(in) :: dqp
    !> harmonic properties on the integration mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> harmonic propertors on the sigma mesh
    type(lo_phonon_dispersions), intent(in) :: ddr
    !> temperature
    real(r8), intent(in) :: temperature
    !> max energy to consider
    real(r8), intent(in) :: max_energy
    !> number of energies
    integer, intent(in) :: n_energy
    !> how are we going to integrate
    integer, intent(in) :: integrationtype
    !> adjustment of smearing parameter
    real(r8), intent(in) :: sigma
    !> include isotope scattering
    logical, intent(in) :: use_isotope
    !> include thirdorder scattering
    logical, intent(in) :: use_thirdorder
    !> include fourthorder scattering
    logical, intent(in) :: use_fourthorder
    !> include off-diagonal terms in the self-energy
    logical, intent(in) :: offdiagonal
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot
    integer, intent(in) :: verbosity

    type(lo_hdf5_helper) :: h5
    integer :: writerank


    ! We start by dumping some metadata to file
    init: block
        writerank=mw%n-1

        if ( mw%r .eq. writerank ) then
            call h5%init(__FILE__, __LINE__)
            call h5%open_file('write', trim(filename))
        endif

        ! Dump the q-grid, we need to have it embedded.
        if ( mw%r .eq. writerank ) then
            call h5%open_group('write','qmesh')
            call dqp%write_to_file(uc,'null',mem,verbosity,input_id=h5%group_id)
            call h5%close_group()
        endif

        if ( mw%r .eq. writerank ) then
            call h5%store_attribute(temperature, h5%file_id, 'temperature')
        endif

        if ( mw%talk ) then
            write(*,*) 'Evaluating self-energy on a rough mesh'
            write(*,*) '    N q-points self-energy:',dqp%n_full_point
            write(*,*) '    N q-points integration:',qp%n_full_point
        endif
    end block init

    ! Calculate the actual self-energy on the mesh. Let's consider
    ! the option of separate meshes later. Maybe that's a good idea,
    ! maybe not.
    selfenergy: block
        real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]
        type(lo_timer) :: tmr
        type(lo_phonon_selfenergy) :: se
        complex(r8), dimension(:,:), allocatable :: eig,inveig,sigma_mode,halfproduct,sigma_xyz
        complex(r8), dimension(:,:,:), allocatable :: bufc
        real(r8), dimension(:,:,:), allocatable :: bufr
        real(r8), dimension(:,:), allocatable :: buf_spectral
        real(r8), dimension(:), allocatable :: buf_taper,buf_sigmaIm,buf_sigmaRe
        real(r8) :: f0
        integer :: iq,ie,imode

        call tmr%start()
        do iq=1,dqp%n_irr_point

            if ( mw%talk ) then
                write(*,*) '... evaluating qpoint '//tochar(iq)//' out of '//tochar(dqp%n_irr_point)
                write(*,*) '... qv:',dqp%ip(iq)%r
            endif

            call se%generate(&
                dqp%ip(iq), qdir, ddr%iq(iq), uc, fc, fct, fcf, ise, qp, dr, &
                temperature, max_energy, n_energy, integrationtype, sigma,&
                use_isotope, use_thirdorder, use_fourthorder, offdiagonal, &
                tmr, mw, mem, verbosity=0)
            ! Dump the energy axis once we have it
            if ( iq == 1 ) then
            if ( mw%r .eq. writerank ) then
                call h5%store_data(se%energy_axis,h5%file_id,'omega')
            endif
            endif

            ! Convert self-energies to xyz coordinates and dump to file
            if ( mw%r .eq. writerank ) then
                allocate(bufc(se%n_mode,se%n_mode,se%n_energy))
                allocate(bufr(se%n_mode,se%n_mode,se%n_energy))
                allocate(eig(se%n_mode,se%n_mode))
                allocate(inveig(se%n_mode,se%n_mode))
                allocate(sigma_mode(se%n_mode,se%n_mode))
                allocate(sigma_xyz(se%n_mode,se%n_mode))
                allocate(halfproduct(se%n_mode,se%n_mode))
                bufc=0.0_r8
                bufr=0.0_r8
                eig=0.0_r8
                inveig=0.0_r8
                sigma_mode=0.0_r8
                sigma_xyz=0.0_r8
                halfproduct=0.0_r8

                eig=ddr%iq(iq)%egv
                inveig=transpose(conjg(eig))

                call h5%open_group('write','selfenergy_qpoint_'//tochar(iq))

                bufc=0.0_r8
                do ie=1,se%n_energy
                    sigma_mode=0.0_r8
                    do imode=1,se%n_mode
                        sigma_mode(imode,imode)=lo_imag*(se%im_3ph(ie,imode) + se%im_iso(ie,imode)) + se%re_3ph(ie,imode) + se%re_4ph(ie,imode)
                    enddo
                    call lo_gemm(eig,sigma_mode,halfproduct)
                    call lo_gemm(halfproduct,inveig,sigma_xyz)
                    bufc(:,:,ie)=sigma_xyz
                enddo
                bufr=aimag(bufc)
                call h5%store_data(bufr,h5%group_id,'sigma_Im')
                bufr=real(bufc,r8)
                call h5%store_data(bufr,h5%group_id,'sigma_Re')

                deallocate(bufr)
                deallocate(bufc)
                deallocate(eig)
                deallocate(inveig)
                deallocate(sigma_mode)
                deallocate(sigma_xyz)
                deallocate(halfproduct)

                ! While we are at it, we might as well dump the spectral functions as well?
                allocate(buf_spectral(se%n_energy,se%n_mode))
                allocate(buf_sigmaIm(se%n_energy))
                allocate(buf_sigmaRe(se%n_energy))
                allocate(buf_taper(se%n_energy))
                buf_spectral = 0.0_r8
                buf_sigmaIm = 0.0_r8
                buf_sigmaRe = 0.0_r8
                buf_taper = 0.0_r8
                do imode = 1, ddr%n_mode
                    ! Get the spectral functions
                    if (ddr%iq(iq)%omega(imode) .gt. lo_freqtol) then
                        buf_sigmaIm = se%im_3ph(:, imode) + se%im_iso(:, imode)
                        buf_sigmaRe = se%re_3ph(:, imode) + se%re_4ph(:, imode)
                        call taperfn_im(se%energy_axis, se%energy_axis(se%n_energy), ddr%iq(iq)%omega(imode), buf_taper)
                        buf_sigmaIm = buf_sigmaIm*buf_taper
                        call evaluate_spectral_function(se%energy_axis, buf_sigmaIm, buf_sigmaRe, ddr%iq(iq)%omega(imode), buf_spectral(:, imode))

                        ! Normalize both.
                        f0 = lo_trapezoid_integration(se%energy_axis, buf_spectral(:, imode))
                        buf_spectral(:, imode) = buf_spectral(:, imode)/f0
                    else
                        buf_spectral(:, imode) = 0.0_r8
                    end if
                end do
                call h5%store_data(buf_spectral,h5%group_id,'spectralfunction')

                deallocate(buf_spectral)
                deallocate(buf_sigmaIm)
                deallocate(buf_sigmaRe)
                deallocate(buf_taper)

                call h5%close_group()

            endif
        enddo
        call tmr%tock('convert and write to file')

        ! And we are done with writing
        if ( mw%r .eq. writerank ) then
            call h5%close_file()
            call h5%destroy(__FILE__, __LINE__)
        endif

        call tmr%stop()
        call tmr%dump(mw,'Self-energy timings')
    end block selfenergy

end subroutine

end module
