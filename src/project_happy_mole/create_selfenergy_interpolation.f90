module create_selfenergy_interpolation
!!
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_exitcode_symmetry, lo_twopi, &
    lo_bohr_to_A, lo_freqtol, lo_exitcode_param, lo_pi, lo_imag, lo_groupvel_ms_to_Hartreebohr, lo_frequency_Hartree_to_THz
use gottochblandat, only: tochar, walltime, lo_progressbar_init, lo_progressbar, lo_chop, &
        lo_linear_least_squares, lo_rsquare, lo_linspace, lo_linear_interpolation, lo_trapezoid_integration,&
        lo_clean_fractional_coordinates, lo_cross
use geometryfunctions, only: lo_inscribed_sphere_in_box, lo_plane, lo_increment_dimensions, lo_bounding_sphere_of_box
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
use type_distancetable, only: lo_distancetable
use lo_voronoi, only: lo_voronoi_cell

use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lo_evaluate_phonon_self_energy, only: lo_phonon_selfenergy
use lo_distributed_phonon_dispersion_relations, only: lo_distributed_phonon_dispersions
use lo_spectralfunction_helpers, only: lo_evaluate_spectral_function,lo_tapering_function,lo_make_eigenvector_parallel,lo_permute_eigenpairs

implicit none

private
public :: generate_interpolated_selfenergy

contains

!> create the interpolated self-energy object and store it to file
subroutine generate_interpolated_selfenergy(filename,uc,fc,fct,fcf,ise,qp,dqp,dr,ddr,pdr, &
    temperature, max_energy, n_energy, integrationtype, sigma,&
    use_isotope, use_thirdorder, use_fourthorder, &
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
    !> distributed harmonic properties on the integration mesh
    type(lo_distributed_phonon_dispersions), intent(in) :: pdr
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
        writerank=0 !mw%n-1

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
    if ( fc%polar ) then
    auxiliaryfc: block
        call create_auxiliary_ifc(uc,fc,aux_fc,mw,mem)

        ! Makes sense to write this to file. Plaintext for now, have to fix hdf5 eventually
        if ( mw%r .eq. writerank ) then
            call aux_fc%writetofile(uc,'outfile.aux_forceconstant')
        endif

        ! And calculate the dispersions for the auxiliary IFCs
        call aux_dr%generate(dqp,aux_fc,uc,mw,mem,verbosity)
    end block auxiliaryfc
    endif

    ! Calculate the actual self-energy on the mesh.
    selfenergy: block
        real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]
        type(lo_timer) :: tmr
        type(lo_phonon_selfenergy) :: se
        complex(r8), dimension(:,:), allocatable :: eig,inveig,sigma_mode,halfproduct,sigma_xyz,cm0,cm1,cm2
        !complex(r8), dimension(:,:), allocatable :: aux_eig,aux_inveig,aux_trafo,aux_dynmat
        complex(r8), dimension(:,:), allocatable :: left_trf,right_trf,aux_dynmat
        complex(r8), dimension(:,:,:), allocatable :: bufc,symmetryop
        real(r8), dimension(:,:,:), allocatable :: bufr
        real(r8), dimension(:), allocatable :: aux_omega
        real(r8) :: f0,f1,fm0,fm1
        integer :: iq,ie,imode,iatom,ialpha,ii,i,iop

        call tmr%start()
        do iq=1,dqp%n_irr_point

            if ( mw%talk ) then
                write(*,*) '... evaluating qpoint '//tochar(iq)//' out of '//tochar(dqp%n_irr_point)
            endif

            call se%generate(&
                dqp%ip(iq), qdir, uc, fc, fct, fcf, ise, qp, pdr,&
                temperature, max_energy, n_energy, integrationtype, sigma,&
                use_isotope, use_thirdorder, use_fourthorder,&
                tmr, mw, mem, verbosity=0)

            ! Dump the energy axis once we have it
            if ( iq == 1 ) then
                if ( mw%r .eq. writerank ) then
                    call h5%store_data(se%energy_axis,h5%file_id,'omega')
                endif
            endif

            if ( mw%r .eq. writerank ) then
                call h5%open_group('write','se_qpoint_'//tochar(iq))
                    ! store reference harmonic
                    call h5%store_data(ddr%iq(iq)%omega, h5%group_id, 'omega')
                    call h5%store_data(aimag(ddr%iq(iq)%egv), h5%group_id, 'im_egv')
                    call h5%store_data(real(ddr%iq(iq)%egv,r8), h5%group_id, 're_egv')
                    ! store actual self-energies
                    call h5%store_data(se%im, h5%group_id, 'sigma_im')
                    call h5%store_data(se%re, h5%group_id, 'sigma_re')
                call h5%close_group()
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
                allocate(left_trf(uc%na*3,uc%na*3))
                allocate(right_trf(uc%na*3,uc%na*3))
                allocate(aux_dynmat(uc%na*3,uc%na*3))
                allocate(aux_omega(uc%na*3))

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
                left_trf=0.0_r8
                right_trf=0.0_r8

                ! Construct transformation from mode-mode to xyz for the actual dispersions
                if ( fc%polar ) then

                    ! Get somewhat reasonable smooth frequencies
                    call aux_fc%dynamicalmatrix(uc,dqp%ip(iq),aux_dynmat,mem)
                    call lo_gemm(aux_dynmat,ddr%iq(iq)%egv,cm1)
                    call lo_gemm(ddr%iq(iq)%egv,cm1,cm0,transa='C')
                    do imode=1,uc%na*3
                        aux_omega(imode)=sqrt( abs(cm0(imode,imode)) )
                    enddo

                    eig=ddr%iq(iq)%egv
                    inveig=transpose(conjg(eig))
                    do imode=1,ddr%n_mode
                        if ( ddr%iq(iq)%omega(imode) .gt. lo_freqtol ) then
                            f0=sqrt(ddr%iq(iq)%omega(imode))
                            f0=f0/aux_omega(imode)
                            f1=1.0_r8/f0
                        else
                            f0=0.0_r8
                            f1=0.0_r8
                        endif
                        do iatom=1,uc%na
                            do ialpha=1,3
                                ii=(iatom-1)*3 + ialpha
                                eig(ii,imode)=eig(ii,imode)*f0
                                inveig(imode,ii)=inveig(imode,ii)*f0
                            enddo
                        enddo
                    enddo
                    do imode=1,ddr%n_mode
                        do iatom=1,uc%na
                            fm0=uc%invsqrtmass(iatom)
                            fm1=1.0_r8/fm0
                            do ialpha=1,3
                                ii=(iatom-1)*3 + ialpha
                                eig(ii,imode)=eig(ii,imode)*fm1
                                inveig(imode,ii)=inveig(imode,ii)*fm1
                            enddo
                        enddo
                    enddo

                    left_trf=eig
                    right_trf=inveig
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

                ! Start with imaginary part of self-energy
                bufc=0.0_r8
                do ie=1,se%n_energy
                    sigma_mode=0.0_r8
                    do imode=1,se%n_mode
                        sigma_mode(imode,imode)=lo_imag*se%im(ie,imode)
                    enddo

                    call lo_gemm(left_trf,sigma_mode,halfproduct)
                    call lo_gemm(halfproduct,right_trf,sigma_xyz)

                    ! Seems like a reasonable thing to enforce symmetry? Not sure if it really helps.
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
                call h5%store_data(bufr,h5%group_id,'sigma_Im_i')
                bufr=real(bufc,r8)
                call h5%store_data(bufr,h5%group_id,'sigma_Im_r')

                ! Then real part
                bufc=0.0_r8
                do ie=1,se%n_energy
                    sigma_mode=0.0_r8
                    do imode=1,se%n_mode
                        sigma_mode(imode,imode)=se%re(ie,imode)
                    enddo

                    call lo_gemm(left_trf,sigma_mode,halfproduct)
                    call lo_gemm(halfproduct,right_trf,sigma_xyz)

                    ! Seems like a reasonable thing to enforce symmetry?
                    ! Not sure if it really helps.
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
                call h5%store_data(bufr,h5%group_id,'sigma_Re_i')
                bufr=real(bufc,r8)
                call h5%store_data(bufr,h5%group_id,'sigma_Re_r')

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
                deallocate(left_trf)
                deallocate(right_trf)
                deallocate(aux_dynmat)
                deallocate(aux_omega)

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

! subroutine lo_set_up_fourier_interpolation(fc,p,qp,mem)
!     !> forceconstant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> structure
!     type(lo_crystalstructure), intent(inout) :: p
!     !> q-mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> memory tracker
!     type(lo_mem_helper), intent(inout) :: mem

!     type(lo_voronoi_cell) :: voro
!     type(lo_crystalstructure) :: ss
!     integer, dimension(3) :: griddensity
!     integer, dimension(:,:), allocatable :: vectormapping

!     init: block
!         select type(qp)
!         type is(lo_fft_mesh)
!             griddensity=qp%griddensity
!         end select

!         write(*,*) 'density:',griddensity

!         call p%build_supercell(ss,dimensions=griddensity)
!         call ss%classify('supercell',p)
!     end block init

!     voronoicell: block
!         type(lo_distancetable) :: dt
!         real(r8), dimension(3,1) :: dummypos
!         real(r8) :: f0

!         f0=lo_bounding_sphere_of_box(ss%latticevectors)*3
!         dummypos=0.0_r8
!         call dt%generate(dummypos,ss%latticevectors,f0,-1)
!         call voro%generate(dt%particle(1),f0*2,1E-6_r8,mem)
!     end block voronoicell

!     mapvectors: block
!         real(r8), dimension(3) :: v0,v1,v2
!         integer, dimension(:,:,:,:,:), allocatable :: dj
!         integer, dimension(3) :: gi
!         integer :: j,k,l,a1,a2,ii,jj,kk

!         allocate(dj(griddensity(1),griddensity(2),griddensity(3),p%na,p%na))
!         dj=0

!         k=0
!         do a1=1,p%na
!         do j=1,ss%na
!             v0=ss%rcart(:,j)-p%rcart(:,a1)
!             a2=ss%info%index_in_unitcell(j)

!             ! Sanity check
!             v1=ss%info%cellindex(:,j)
!             v1=matmul(p%latticevectors,(v1-1.0_r8)+p%r(:,a2))-p%rcart(:,a1)
!             if ( norm2(v0-v1) .gt. 1E-6_r8 ) then
!                 call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
!             endif

!             l=0
!             do ii=-3,3
!             do jj=-3,3
!             do kk=-3,3
!                 v1=real([ii,jj,kk],r8)
!                 v1=matmul(ss%latticevectors,v1)
!                 v1=v1+v0
!                 if ( voro%is_point_inside(v1,1E-5_r8) ) then
!                     k=k+1
!                     l=l+1
!                 endif
!             enddo
!             enddo
!             enddo
!             if ( l .eq. 0 ) then
!                 call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
!             endif

!             ! Make a note of the multiplicity
!             dj( ss%info%cellindex(1,j),ss%info%cellindex(2,j),ss%info%cellindex(3,j),a1,a2)=l
!         enddo
!         enddo

!         do a1=1,p%na
!         do a2=1,p%na
!             write(*,*) 'nrm',a1,a2,sum(dj(:,:,:,a1,a2))
!         enddo
!         enddo

!         ! Then we do it again and store some things
!         allocate(vectormapping(9,k))
!         vectormapping=0

!         k=0
!         do a1=1,p%na
!         do j=1,ss%na

!             v0=ss%rcart(:,j)-p%rcart(:,a1)
!             a2=ss%info%index_in_unitcell(j)

!             l=0
!             do ii=-1,1
!             do jj=-1,1
!             do kk=-1,1
!                 v1=real([ii,jj,kk],r8)
!                 v1=matmul(ss%latticevectors,v1)
!                 v1=v1+v0
!                 if ( voro%is_point_inside(v1,1E-5_r8) ) then
!                     k=k+1
!                     ! Actual lattice vector?
!                     ! v1 = v(j) - v(i)
!                     ! v1 = lv(j) + r(a2) - r(a1)
!                     ! lv(j) = v1 -r(a2)+r(a1)
!                     v2=v1 - p%rcart(:,a2) + p%rcart(:,a1)
!                     v2=matmul(p%inv_latticevectors,v2)
!                     if ( norm2(v2-anint(v2)) .gt. 1E-6_r8 ) then
!                         call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
!                     else
!                         gi=int(anint(v2))
!                     endif

!                     ! store actual lattice vector
!                     vectormapping(1:3,k)=gi
!                     ! store unitcell atom indices
!                     vectormapping(4,k)=a1
!                     vectormapping(5,k)=a2
!                     ! store weight
!                     vectormapping(6,k)=dj( ss%info%cellindex(1,j),ss%info%cellindex(2,j),ss%info%cellindex(3,j),a1,a2)
!                     ! store Fourier indexing thing
!                     vectormapping(7:9,k)=ss%info%cellindex(:,j)
!                 endif
!             enddo
!             enddo
!             enddo
!         enddo
!         enddo

!         deallocate(dj)
!     end block mapvectors

!     inversetransform: block
!         complex(r8), dimension(:,:), allocatable :: dm0,dm4,rotmat,cm0
!         complex(r8), dimension(:,:,:,:,:), allocatable :: dm1
!         complex(r8), dimension(:,:,:), allocatable :: dm2
!         complex(r8) :: expikr
!         real(r8), dimension(3) :: v0
!         real(r8) :: kdotr,weight
!         integer :: iq,ir,i,j,k,a1,a2
!         integer :: jq,iop


!         allocate(dm0(p%na*3,p%na*3))
!         allocate(dm1(p%na*3,p%na*3,griddensity(1),griddensity(2),griddensity(3)))
!         allocate(rotmat(p%na*3,p%na*3))
!         allocate(cm0(p%na*3,p%na*3))
!         dm0=0.0_r8
!         dm1=0.0_r8
!         rotmat=0.0_r8
!         cm0=0.0_r8

!         do iq=1,qp%n_full_point
!             ! get dynamical matrix at this q
!             jq=qp%ap(iq)%irreducible_index
!             iop=qp%ap(iq)%operation_from_irreducible
!             if ( iop .gt. 0 ) then
!                 call lo_eigenvector_transformation_matrix(rotmat,p%rcart,qp%ip( jq )%r,p%sym%op(iop),inverseoperation=.false.)
!             else
!                 call lo_eigenvector_transformation_matrix(rotmat,p%rcart,qp%ip( jq )%r,p%sym%op(-iop),inverseoperation=.true.)
!             endif
!             call fc%dynamicalmatrix(p,qp%ip(jq),dm0,mem)
!             call lo_gemm(rotmat,dm0,cm0)
!             call lo_gemm(cm0,rotmat,dm0,transb='C')

!             !call fc%dynamicalmatrix(p,qp%ap(iq),dm0,mem)
!             ! naive fourier transform
!             do i=1,griddensity(1)
!             do j=1,griddensity(2)
!             do k=1,griddensity(3)
!                 v0=[i,j,k]-1.0_r8
!                 v0=matmul(p%latticevectors,v0)
!                 kdotr=-dot_product(v0,qp%ap(iq)%r)*lo_twopi
!                 expikr=cmplx(cos(kdotr),sin(kdotr),r8)
!                 dm1(:,:,i,j,k)=dm1(:,:,i,j,k)+dm0*expikr
!             enddo
!             enddo
!             enddo
!         enddo
!         dm1=dm1/real(product(griddensity),r8)

!         ! Rearrange somewhat so that we can interpolate?
!         allocate(dm2(3,3,size(vectormapping,2)))
!         dm2=0.0_r8
!         do ir=1,size(vectormapping,2)
!             a1=vectormapping(4,ir)
!             a2=vectormapping(5,ir)
!             weight=1.0_r8/real(vectormapping(6,ir),r8)
!             dm2(:,:,ir)=weight*dm1((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,vectormapping(7,ir),vectormapping(8,ir),vectormapping(9,ir))
!         enddo
!         dm2=real(dm2,r8)

!         ! Re-interpolate?
!         allocate(dm4(p%na*3,p%na*3))
!         do iq=1,qp%n_full_point
!             ! get dynamical matrix at this q
!             call fc%dynamicalmatrix(p,qp%ap(iq),dm0,mem)
!             ! naive fourier transform
!             dm4=0.0_r8
!             do ir=1,size(vectormapping,2)
!                 v0=vectormapping(1:3,ir)
!                 v0=matmul(p%latticevectors,v0)
!                 kdotr=dot_product(v0,qp%ap(iq)%r)*lo_twopi
!                 expikr=cmplx(cos(kdotr),sin(kdotr),r8)
!                 a1=vectormapping(4,ir)
!                 a2=vectormapping(5,ir)
!                 dm4((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3)=dm4((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3) + dm2(:,:,ir)*expikr
!             enddo

!             !write(*,*) iq,sum(abs(dm0)),sum(abs(dm0-dm4))
!         enddo

!         !write(*,*) 'sum',sum(abs(aimag(dm2))),sum(abs(real(dm2)))

!         ! write(*,*) 'done here for now!'
!         ! stop

!     end block inversetransform

! end subroutine

subroutine create_auxiliary_ifc(uc,fc,aux_fc,mw,mem)
    !> structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> original forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> non-polar smooth forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: aux_fc
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), parameter :: safety_margin_for_cutoff=0.3_r8
    type(lo_crystalstructure) :: ss
    type(lo_interaction_tensors) :: slt
    type(lo_forcemap) :: map

    class(lo_qpoint_mesh), allocatable :: dummy_qp
    type(lo_qpoint), dimension(:), allocatable :: list_qp
    complex(r8), dimension(:,:,:), allocatable :: dynmat,lrdynmat
    real(r8), dimension(:,:,:,:), allocatable :: rottensor
    real(r8), dimension(:), allocatable :: hermitian_rhs, huang_rhs, rotational_rhs
    real(r8), dimension(3, 3, 3, 3) :: bracket,lrbracket
    real(r8), dimension(3,3) :: m0
    real(r8), dimension(3) :: r
    real(r8) :: rc2,f0
    integer, dimension(3) :: supercelldim
    integer :: i,j,k,l,ii,nq
    integer :: a1,pair,al,be,gm,la

    rc2=1.5_r8*fc%cutoff
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

    ! Sort out the symmetry thingies
    call uc%build_supercell(ss,supercelldim)
    call ss%classify('supercell',uc)
    call slt%generate(uc, ss, rc2, -1.0_r8, -1.0_r8, .false., mw, mem, -1)
    call map%generate(uc, ss, polarcorrectiontype=0, st=slt, mw=mw, mem=mem, verbosity=-1)

    ! Create a q-mesh?
    call lo_generate_qmesh(dummy_qp, uc, supercelldim, 'fft', timereversal=.true., headrankonly=.false., mw=mw, mem=mem, verbosity=-1)

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

    ! Calculate the elastic constants bracket
    bracket=0.0_r8
    do a1 = 1, fc%na
    do pair = 1, fc%atom(a1)%n
        m0 = fc%atom(a1)%pair(pair)%m
        r = fc%atom(a1)%pair(pair)%r
        do la = 1, 3
        do gm = 1, 3
        do be = 1, 3
        do al = 1, 3
            bracket(al, be, gm, la) = bracket(al, be, gm, la) + m0(al, be)*r(gm)*r(la)
        end do
        end do
        end do
        end do
    end do
    end do
    ! and long-range bracket, if applicable
    if ( fc%loto%correctiontype .gt. 0 ) then
        allocate(rottensor(3,3,3,uc%na))
        rottensor=0.0_r8
        call fc%ew%longrange_elastic_constant_bracket(uc,fc%loto%born_effective_charges,fc%loto%eps,lrbracket,rottensor,reconly=.true.,mw=mw, verbosity=-1)
    endif
    ! Let's see if we can flatten this
    allocate(hermitian_rhs(9*uc%na))
    allocate(rotational_rhs(27*uc%na))
    allocate(huang_rhs(81))
    hermitian_rhs=0.0_r8
    rotational_rhs=0.0_r8
    huang_rhs=0.0_r8

    ii=0
    do i=1,3
    do j=1,3
    do k=1,3
    do l=1,3
        ii=ii+1
        huang_rhs(ii) = bracket(i,j,k,l) + lrbracket(i,j,k,l)
    enddo
    enddo
    enddo
    enddo

    call map%forceconstant_constraints(uc,rotational=.true.,huanginvariances=.true.,hermitian=.true.,elastic=.true.,hermitian_rhs=hermitian_rhs,huang_rhs=huang_rhs,rotational_rhs=rotational_rhs,verbosity=-1)

    ! Fit non-polar IFC to dynamical matrices
    call lo_irreducible_forceconstant_from_qmesh_dynmat(map, uc, list_qp, nq, dynmat, lrdynmat, .true., .true., mw, -1, .true.) !, weights=wts)
    ! And finally, return a new non-polar second order forceconstant.
    call map%get_secondorder_forceconstant(uc,aux_fc,mem,-1)

    ! ! Then again, to match eigenvectors somehow:
    ! allocate(dm0(uc%na*3,uc%na*3))
    ! allocate(dm1(uc%na*3,uc%na*3))
    ! dm0=0.0_r8
    ! dm1=0.0_r8
    ! nq=0

    ! do k=1,1

    ! do i=1,dummy_qp%n_irr_point
    !     if ( mod(i,mw%n) .ne. mw%r ) cycle
    !     nq=nq+1
    !     !call fc%dynamicalmatrix(uc,dummy_qp%ip(i),dynmat(:,:,nq),mem)
    !     call aux_wp1%generate(fc,uc,mem,dummy_qp%ip(i))
    !     call aux_wp2%generate(aux_fc,uc,mem,dummy_qp%ip(i))

    !     dm0=0.0_r8
    !     do j=1,uc%na*3
    !         dm0(j,j)=aux_wp2%omega(j)**2
    !     enddo
    !     call lo_gemm(aux_wp1%egv,dm0,dm1)
    !     call lo_gemm(dm1,aux_wp1%egv,dm0,transb='C')
    !     dynmat(:,:,nq)=dm0
    !     list_qp(nq)%r = dummy_qp%ip(i)%r
    ! enddo
    ! call lo_irreducible_forceconstant_from_qmesh_dynmat(map, uc, list_qp, nq, dynmat, lrdynmat, .true., .true., mw, verbosity, .true.) !, weights=wts)
    ! ! And finally, return a new non-polar second order forceconstant.
    ! call map%get_secondorder_forceconstant(uc,aux_fc,mem,-1)

    ! enddo
end subroutine

end module
