module create_selfenergy_interpolation
!!
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_exitcode_symmetry, lo_twopi, &
    lo_bohr_to_A, lo_freqtol, lo_exitcode_param, lo_pi, lo_imag, lo_groupvel_ms_to_Hartreebohr, lo_frequency_Hartree_to_THz
use gottochblandat, only: tochar, walltime, lo_progressbar_init, lo_progressbar, lo_chop, &
        lo_linear_least_squares, lo_rsquare, lo_linspace, lo_linear_interpolation, lo_trapezoid_integration,&
        lo_clean_fractional_coordinates, lo_cross
use geometryfunctions, only: lo_inscribed_sphere_in_box, lo_plane, lo_increment_dimensions
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_wedge_mesh, lo_qpoint, lo_read_qmesh_from_file, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm, lo_dgels
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_forcemap, only: lo_forcemap, lo_secondorder_rot_herm_huang
use ifc_solvers, only: lo_irreducible_forceconstant_from_qmesh_dynmat

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
    type(lo_crystalstructure), intent(inout) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> tabulated and smeared spectral functions
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> q-point mesh for self-energy integrations
    class(lo_qpoint_mesh), intent(in) :: qp
    !> q-point mesh the self-energy is evaluated on
    class(lo_qpoint_mesh), intent(inout) :: dqp
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
    type(lo_forceconstant_secondorder) :: aux_fc
    type(lo_phonon_dispersions) :: aux_dr
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
            call h5%store_attribute(fc%polar, h5%file_id, 'polar')
        endif

        if ( mw%talk ) then
            write(*,*) 'Evaluating self-energy on a rough mesh'
            write(*,*) '    N q-points self-energy:',dqp%n_full_point
            write(*,*) '    N q-points integration:',qp%n_full_point
        endif
    end block init

    ! If we have polar interactions we might need an auxiliary non-polar phonon dispersion
    auxiliaryfc: block
        real(r8), parameter :: safety_margin_for_cutoff=0.3_r8
        type(lo_crystalstructure) :: ss
        type(lo_interaction_tensors) :: slt
        type(lo_forcemap) :: map

        class(lo_qpoint_mesh), allocatable :: dummy_qp
        type(lo_qpoint), dimension(:), allocatable :: list_qp
        complex(r8), dimension(:,:,:), allocatable :: dynmat,lrdynmat
        real(r8), dimension(3,3) :: m0
        real(r8) :: rc2,f0
        integer, dimension(3) :: supercelldim
        integer :: i,nq

        ! Use the same cutoff as the input IFCs? Suppose that makes sense.
        rc2=fc%cutoff

        ! Then I need a supercell large enough to safely fit the cutoff
        supercelldim=1
        do
            do i = 1, 3
                m0(:, i) = uc%latticevectors(:,i)*supercelldim(i)
            end do
            f0 = lo_inscribed_sphere_in_box(m0)
            if (f0 .gt. rc2+safety_margin_for_cutoff ) exit
            supercelldim = lo_increment_dimensions(supercelldim, uc%latticevectors)
        end do
        supercelldim=supercelldim

        ! Sort out the symmetry thingies
        call uc%build_supercell(ss,supercelldim)
        call ss%classify('supercell',uc)
        call slt%generate(uc, ss, rc2, -1.0_r8, -1.0_r8, .false., mw, mem, -1)
        call map%generate(uc, ss, polarcorrectiontype=0, st=slt, mw=mw, mem=mem, verbosity=-1)

        ! Create a q-mesh?
        call lo_generate_qmesh(dummy_qp, uc, supercelldim+3, 'fft', timereversal=.true., headrankonly=.false., mw=mw, mem=mem, verbosity=-1)

        ! Generate a bunch of dynamical matrices?
        nq=0
        do i=1,dummy_qp%n_irr_point
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            nq=nq+1
        enddo

        if ( nq .gt. 0 ) then
            allocate(dynmat(uc%na*3,uc%na*3,nq))
            allocate(lrdynmat(uc%na*3,uc%na*3,nq))
            allocate(list_qp(nq))
            dynmat=0.0_r8
            lrdynmat=0.0_r8
        endif

        nq=0
        do i=1,dummy_qp%n_irr_point
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            nq=nq+1
            call fc%dynamicalmatrix(uc,dummy_qp%ip(i),dynmat(:,:,nq),mem)
            list_qp(nq)%r = dummy_qp%ip(i)%r
        enddo

        ! Fit non-polar IFC to dynamical matrices
        call lo_irreducible_forceconstant_from_qmesh_dynmat(map, uc, list_qp, nq, dynmat, lrdynmat, .true., .true., mw, verbosity, .true.)

        ! And finally, return a new non-polar second order forceconstant.
        call map%get_secondorder_forceconstant(uc,aux_fc,mem,-1)

        ! Makes sense to write this to file. Plaintext for now, have to fix hdf5 eventually
        if ( mw%talk ) then
            call aux_fc%writetofile(uc,'outfile.aux_forceconstant')
        endif

        ! And calculate the dispersions for the auxiliary IFCs
        call aux_dr%generate(dqp,aux_fc,uc,mw,mem,verbosity)
    end block auxiliaryfc

    ! Calculate the actual self-energy on the mesh. Let's consider
    ! the option of separate meshes later. Maybe that's a good idea,
    ! maybe not.
    selfenergy: block
        real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]
        type(lo_timer) :: tmr
        type(lo_phonon_selfenergy) :: se
        complex(r8), dimension(:,:), allocatable :: eig,inveig,sigma_mode,halfproduct,sigma_xyz,cm0,cm1,cm2
        complex(r8), dimension(:,:), allocatable :: aux_eig,aux_inveig,aux_trafo
        complex(r8), dimension(:,:), allocatable :: left_trf,right_trf
        complex(r8), dimension(:,:,:), allocatable :: bufc,symmetryop
        real(r8), dimension(:,:,:), allocatable :: bufr
        real(r8), dimension(:,:), allocatable :: buf_spectral
        real(r8), dimension(:), allocatable :: buf_taper,buf_sigmaIm,buf_sigmaRe
        real(r8) :: f0,f1,fm0,fm1
        integer :: iq,ie,imode,iatom,ialpha,ii,i,iop

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
                allocate(cm0(uc%na*3,uc%na*3))
                allocate(cm1(uc%na*3,uc%na*3))
                allocate(cm2(uc%na*3,uc%na*3))
                allocate(aux_eig(uc%na*3,uc%na*3))
                allocate(aux_inveig(uc%na*3,uc%na*3))
                allocate(aux_trafo(uc%na*3,uc%na*3))
                allocate(left_trf(uc%na*3,uc%na*3))
                allocate(right_trf(uc%na*3,uc%na*3))

                bufc=0.0_r8
                bufr=0.0_r8
                eig=0.0_r8
                inveig=0.0_r8
                sigma_mode=0.0_r8
                sigma_xyz=0.0_r8
                halfproduct=0.0_r8
                cm0=0.0_r8
                cm1=0.0_r8
                cm2=0.0_r8
                aux_eig=0.0_r8
                aux_inveig=0.0_r8
                aux_trafo=0.0_r8
                left_trf=0.0_r8
                right_trf=0.0_r8

                ! Construct transformation from mode-mode to xyz for the actual dispersions
                if ( fc%polar ) then
                    eig=ddr%iq(iq)%egv
                    inveig=transpose(conjg(eig))
                    do imode=1,ddr%n_mode
                        if ( ddr%iq(iq)%omega(imode) .gt. lo_freqtol ) then
                            f0=sqrt(ddr%iq(iq)%omega(imode))
                            f1=1.0_r8/f0
                        else
                            f0=0.0_r8
                            f1=0.0_r8
                        endif
                        do iatom=1,uc%na
                            fm0=uc%invsqrtmass(iatom)
                            fm1=1.0_r8/fm0
                            do ialpha=1,3
                                ii=(iatom-1)*3 + ialpha
                                eig(ii,imode)=eig(ii,imode)*f0
                                inveig(imode,ii)=inveig(imode,ii)*f0
                            enddo
                        enddo
                    enddo

                    ! Then we construct an additional transformation to frequency-weight the interpolation.
                    aux_eig=aux_dr%iq(iq)%egv
                    aux_inveig=transpose(conjg(aux_eig))
                    do imode=1,aux_dr%n_mode
                        if ( aux_dr%iq(iq)%omega(imode) .gt. lo_freqtol ) then
                            f0=sqrt(aux_dr%iq(iq)%omega(imode))
                            f1=1.0_r8/f0
                        else
                            f0=0.0_r8
                            f1=0.0_r8
                        endif
                        do iatom=1,uc%na
                            fm0=uc%invsqrtmass(iatom)
                            fm1=1.0_r8/fm0
                            do ialpha=1,3
                                ii=(iatom-1)*3 + ialpha
                                aux_eig(ii,imode)=aux_eig(ii,imode)*f1
                                aux_inveig(imode,ii)=aux_inveig(imode,ii)*f1
                            enddo
                        enddo
                    enddo
                    call lo_gemm(aux_eig,aux_inveig,aux_trafo)

                    call lo_gemm(aux_trafo,eig,left_trf)
                    call lo_gemm(inveig,aux_trafo,right_trf,transb='C')
                else

                    eig=ddr%iq(iq)%egv
                    inveig=transpose(conjg(eig))
                    do imode=1,ddr%n_mode
                        if ( ddr%iq(iq)%omega(imode) .gt. lo_freqtol ) then
                            f0=sqrt(ddr%iq(iq)%omega(imode))
                            f1=1.0_r8/f0
                        else
                            f0=0.0_r8
                            f1=0.0_r8
                        endif
                        do iatom=1,uc%na
                            fm0=uc%invsqrtmass(iatom)
                            fm1=1.0_r8/fm0
                            do ialpha=1,3
                                ii=(iatom-1)*3 + ialpha
                                eig(ii,imode)=eig(ii,imode)*f1
                                inveig(imode,ii)=inveig(imode,ii)*f1
                            enddo
                        enddo
                    enddo
                    left_trf=eig
                    right_trf=inveig
                endif

                call h5%open_group('write','selfenergy_qpoint_'//tochar(iq))

                ! Need some symmetry operations
                allocate(symmetryop(uc%na*3,uc%na*3,dqp%ip(iq)%n_invariant_operation))
                symmetryop=0.0_r8
                do i=1,dqp%ip(iq)%n_invariant_operation
                    iop = dqp%ip(iq)%invariant_operation(i)
                    call lo_eigenvector_transformation_matrix(symmetryop(:,:,i),uc%rcart,dqp%ip(iq)%r,uc%sym%op(iop))
                enddo

                bufc=0.0_r8
                do ie=1,se%n_energy
                    sigma_mode=0.0_r8
                    do imode=1,se%n_mode
                        sigma_mode(imode,imode)=lo_imag*(se%im_3ph(ie,imode) + se%im_iso(ie,imode)) + se%re_3ph(ie,imode) + se%re_4ph(ie,imode)
                    enddo

                    ! call lo_gemm(eig,sigma_mode,halfproduct)
                    ! call lo_gemm(halfproduct,inveig,sigma_xyz)

                    ! call lo_gemm(aux_trafo,sigma_xyz,halfproduct)
                    ! call lo_gemm(halfproduct,aux_trafo,sigma_xyz,transb='C')

                    call lo_gemm(left_trf,sigma_mode,halfproduct)
                    call lo_gemm(halfproduct,right_trf,sigma_xyz)

                    ! Ok, decent start, now we are in xyz-space. Now the next step is to

                    ! Then it seems we have to enforce symmetry? Or not?
                    cm0=0.0_r8
                    do i=1,dqp%ip(iq)%n_invariant_operation
                        iop = dqp%ip(iq)%invariant_operation(i)
                        call lo_gemm(symmetryop(:,:,i),sigma_xyz,cm1)
                        call lo_gemm(cm1,symmetryop(:,:,i),cm2,transb='C')
                        cm0=cm0+cm2
                    enddo
                    sigma_xyz=cm0/real(dqp%ip(iq)%n_invariant_operation,r8)

                    bufc(:,:,ie)=sigma_xyz
                enddo
                bufr=aimag(bufc)
                call h5%store_data(bufr,h5%group_id,'sigma_Im')
                bufr=real(bufc,r8)
                call h5%store_data(bufr,h5%group_id,'sigma_Re')

                deallocate(symmetryop)
                deallocate(cm0)
                deallocate(cm1)
                deallocate(cm2)

                deallocate(bufr)
                deallocate(bufc)
                deallocate(eig)
                deallocate(inveig)
                deallocate(sigma_mode)
                deallocate(sigma_xyz)
                deallocate(halfproduct)
                deallocate(aux_eig)
                deallocate(aux_inveig)
                deallocate(aux_trafo)
                deallocate(left_trf)
                deallocate(right_trf)

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
