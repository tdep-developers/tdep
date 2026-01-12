module lo_selfenergy_interpolation
!! Allows us to evaluate the phonon self-energy at arbitrary q-points
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_tol, lo_sqtol, lo_exitcode_symmetry, lo_twopi, &
                      lo_bohr_to_A, lo_freqtol, lo_exitcode_param, lo_pi, lo_imag, lo_groupvel_ms_to_Hartreebohr
use gottochblandat, only: tochar, walltime, lo_progressbar_init, lo_progressbar, lo_chop, &
                          lo_linear_least_squares, lo_rsquare, lo_linspace, lo_linear_interpolation, lo_trapezoid_integration,&
                          lo_clean_fractional_coordinates, lo_cross, lo_mean, lo_complex_singular_value_decomposition,&
                          lo_real_singular_value_decomposition
use geometryfunctions, only: lo_inscribed_sphere_in_box, lo_plane, lo_bounding_sphere_of_box
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_wedge_mesh, lo_qpoint, lo_read_qmesh_from_file, lo_generate_qmesh, lo_get_small_group_of_qpoint
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm, lo_dgels, lo_zheev
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_distancetable, only: lo_distancetable
use lo_voronoi, only: lo_voronoi_cell

use lo_tetrahedron_interpolation, only: lo_linear_tetrahedron_interpolation
use type_phonon_dos, only: lo_phonon_dos
use lo_thermal_transport, only: lo_thermal_conductivity
use lo_spectralfunction_helpers, only: lo_evaluate_spectral_function,lo_gaussian_smear_spectral_function,lo_find_spectral_function_max_and_fwhm,lo_integrate_spectral_function,lo_tapering_function,lo_make_eigenvector_parallel,lo_permute_eigenpairs
implicit none

private
public :: lo_interpolated_selfenergy_grid

type lo_interpolated_selfenergy_grid
    !> point locator thingy
    type(lo_linear_tetrahedron_interpolation) :: box
    !> q-point mesh used to generate the self-energy
    class(lo_qpoint_mesh), allocatable :: qp
    !> number of energies
    integer :: n_energy=-lo_hugeint
    !> energy axis
    real(r8), dimension(:), allocatable :: omega
    !> real part of self-energy (xyz,xyz,energy,q)
    complex(r8), dimension(:,:,:,:), allocatable :: sigma_Re
    !> imaginary part of self-energy (xyz,xyz,energy,q)
    complex(r8), dimension(:,:,:,:), allocatable :: sigma_Im
    !> Auxiliary IFCs for interpolation
    type(lo_forceconstant_secondorder) :: aux_fc
    !> Is this a polar material?
    logical :: polar=.false.

    ! harmonic omega, per q
    real(r8), dimension(:,:), allocatable :: harm_omega
    complex(r8), dimension(:,:,:), allocatable :: harm_egv
    real(r8), dimension(:,:,:), allocatable :: sIm,sRe

    ! Fourier interpolation thingies
    integer :: n_rvec=-lo_hugeint
    real(r8), dimension(:,:,:,:), allocatable :: bre,bim
    complex(r8), dimension(:,:,:,:), allocatable :: bc
    real(r8), dimension(:,:), allocatable :: rvec
    integer, dimension(:,:), allocatable :: atomind

    complex(r8), dimension(:,:,:,:), allocatable :: rfc
    complex(r8), dimension(:,:,:,:), allocatable :: ifc
    contains
        procedure :: read_from_hdf5=>read_interpolated_selfenergy_from_hdf5
        procedure :: evaluate=>evaluate_self_energy
        procedure :: evaluate_smeared_J=>evaluate_self_energy
        procedure :: destroy=>destroy_interpolated_selfenergy
        procedure :: spectral_function_along_path=>spectral_function_path_interp
        procedure :: spectral_function_on_grid=>spectral_function_grid_interp
end type

interface ! to evaluate
    module subroutine evaluate_self_energy(ise,p,qv,omega,egv,sigma_Re,sigma_Im,mem,mw)
        class(lo_interpolated_selfenergy_grid), intent(inout) :: ise
        type(lo_crystalstructure), intent(in) :: p
        real(r8), dimension(3), intent(in) :: qv
        real(r8), dimension(:), intent(in) :: omega
        complex(r8), dimension(:,:), intent(in) :: egv
        real(r8), dimension(:,:), intent(out) :: sigma_Re
        real(r8), dimension(:,:), intent(out) :: sigma_Im
        type(lo_mem_helper), intent(inout) :: mem
        type(lo_mpi_helper), intent(inout), optional :: mw
    end subroutine
end interface

interface ! to path
    module subroutine spectral_function_path_interp(ise, bs, uc, mw, mem)
        class(lo_interpolated_selfenergy_grid), intent(inout) :: ise
        type(lo_phonon_bandstructure), intent(inout) :: bs
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    module subroutine spectral_function_grid_interp(ise, uc, fc, qp, smearing_prefactor, temperature, tc, pd, dr, mw, mem)
        class(lo_interpolated_selfenergy_grid), intent(inout) :: ise
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        class(lo_qpoint_mesh), intent(inout) :: qp
        real(r8), intent(in) :: smearing_prefactor
        real(r8), intent(in) :: temperature
        type(lo_thermal_conductivity), intent(out) :: tc
        type(lo_phonon_dos), intent(out) :: pd
        type(lo_phonon_dispersions), intent(out) :: dr
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
end interface

contains

subroutine read_interpolated_selfenergy_from_hdf5(ise,p,fc,filename,mw,mem,verbosity)
    !> self-energy
    class(lo_interpolated_selfenergy_grid), intent(out) :: ise
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> filename
    character(len=*), intent(in) :: filename
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! First we grab the raw data from file.
    readfile: block
        real(r8), dimension(:,:,:), allocatable :: rbuf0,rbuf1
        real(r8), dimension(:,:), allocatable :: rb2,rb3
        real(r8), dimension(:), allocatable :: rb1
        type(lo_hdf5_helper) :: h5
        integer :: iq
        ! complex(r8), dimension(3,3) :: cm0,cm1
        ! integer :: iq,jq,ie,a1,a2

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'Reading interpolated self-energy from file'
        endif

        call h5%init(__FILE__,__LINE__)
        call h5%open_file('read',trim(filename))

        ! First we read the q-point mesh from file
        call h5%open_group('read','qmesh')
        call lo_read_qmesh_from_file(ise%qp,p,'null',mem,0,h5%group_id)
        call h5%close_group()

        ! Read some auxiliary data
        call h5%read_attribute(ise%polar,h5%file_id,'polar')

        ! Read the auxiliary IFC from separate file, will put into hdf5 if it works
        if ( ise%polar ) then
            call ise%aux_fc%readfromfile(p,'outfile.aux_forceconstant',mem,-1)
        endif

        if ( verbosity .gt. 0 ) then
            write(*,*) '... read q-mesh'
        endif

        ! Then we grab the energy-axis, seems reasonable
        call h5%read_data(ise%omega,h5%file_id,'omega')
        ise%n_energy = size(ise%omega)

        ! Create buffer space
        allocate(ise%sigma_Im(p%na*3,p%na*3,size(ise%omega),ise%qp%n_irr_point))
        allocate(ise%sigma_Re(p%na*3,p%na*3,size(ise%omega),ise%qp%n_irr_point))
        ise%sigma_Im=0.0_r8
        ise%sigma_Re=0.0_r8
        do iq=1,ise%qp%n_irr_point
            call h5%open_group('read','selfenergy_qpoint_'//tochar(iq))

            call h5%read_data(rbuf0,h5%group_id,'sigma_Re_i')
            call h5%read_data(rbuf1,h5%group_id,'sigma_Re_r')
            ise%sigma_Re(:,:,:,iq)=cmplx(rbuf1,rbuf0,r8)
            deallocate(rbuf0)
            deallocate(rbuf1)

            call h5%read_data(rbuf0,h5%group_id,'sigma_Im_i')
            call h5%read_data(rbuf1,h5%group_id,'sigma_Im_r')
            ise%sigma_Im(:,:,:,iq)=cmplx(rbuf1,rbuf0,r8)
            deallocate(rbuf0)
            deallocate(rbuf1)

            call h5%close_group()
        enddo


        allocate(ise%harm_egv(p%na*3,p%na*3,ise%qp%n_irr_point))
        allocate(ise%harm_omega(p%na*3,ise%qp%n_irr_point))
        allocate(ise%sIm(size(ise%omega),p%na*3,ise%qp%n_irr_point))
        allocate(ise%sRe(size(ise%omega),p%na*3,ise%qp%n_irr_point))

        do iq=1,ise%qp%n_irr_point
            call h5%open_group('read','se_qpoint_'//tochar(iq))

            call h5%read_data(rb1,h5%group_id,'omega')
            ise%harm_omega(:,iq)=rb1
            deallocate(rb1)

            call h5%read_data(rb2,h5%group_id,'re_egv')
            call h5%read_data(rb3,h5%group_id,'im_egv')
            ise%harm_egv(:,:,iq)=cmplx(rb2,rb3,r8)
            deallocate(rb2)
            deallocate(rb3)

            call h5%read_data(rb2,h5%group_id,'sigma_im')
            call h5%read_data(rb3,h5%group_id,'sigma_re')
            ise%sIm(:,:,iq)=rb2
            ise%sRe(:,:,iq)=rb3
            deallocate(rb2)
            deallocate(rb3)

            call h5%close_group()
        enddo


        if ( verbosity .gt. 0 ) then
            write(*,*) '... read self-energy'
        endif

        call h5%close_file()
        call h5%destroy()

        ! Generate triangulation thingy
        call ise%box%generate(ise%qp,p)

        if ( verbosity .gt. 0 ) then
            write(*,*) 'Done reading self-energy from file'
        endif
    end block readfile

    ! Store the grid into a Fourier interpolation thingy?
    fourierinterpolation: block
        type(lo_crystalstructure) :: ss
        type(lo_distancetable) :: dt
        type(lo_voronoi_cell) :: voro

        complex(r8), dimension(:,:,:,:,:), allocatable :: dm0,dm1
        complex(r8), dimension(:,:), allocatable :: rotmat,cm0,cm1,cm2
        complex(r8) :: expikr
        real(r8), dimension(3,1) :: dummypos
        real(r8), dimension(3) :: v0,v1,v2
        real(r8) :: kdotr,f0,f1,weight
        integer, dimension(:,:,:,:,:), allocatable :: dj
        integer, dimension(:,:), allocatable :: vectormapping
        integer, dimension(3) :: griddensity,gi
        integer :: i,j,k,l,ii,jj,kk,ie,iq,jq,iop,a1,a2,ir

        ! Build supercell and Voronoi cell?
        select type(qp=>ise%qp)
        type is(lo_fft_mesh)
            griddensity=qp%griddensity
        end select
        call p%build_supercell(ss,dimensions=griddensity)
        call ss%classify('supercell',p)
        f0=lo_bounding_sphere_of_box(ss%latticevectors)*3
        dummypos=0.0_r8
        call dt%generate(dummypos,ss%latticevectors,f0,-1)
        call voro%generate(dt%particle(1),f0*2,1E-6_r8,mem)

        allocate(dj(griddensity(1),griddensity(2),griddensity(3),p%na,p%na))
        dj=0

        k=0
        do a1=1,p%na
        do j=1,ss%na
            v0=ss%rcart(:,j)-p%rcart(:,a1)
            a2=ss%info%index_in_unitcell(j)

            ! Sanity check
            v1=ss%info%cellindex(:,j)
            v1=matmul(p%latticevectors,(v1-1.0_r8)+p%r(:,a2))-p%rcart(:,a1)
            if ( norm2(v0-v1) .gt. 1E-6_r8 ) then
                call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
            endif

            l=0
            do ii=-3,3
            do jj=-3,3
            do kk=-3,3
                v1=real([ii,jj,kk],r8)
                v1=matmul(ss%latticevectors,v1)
                v1=v1+v0
                if ( voro%is_point_inside(v1,1E-5_r8) ) then
                    k=k+1
                    l=l+1
                endif
            enddo
            enddo
            enddo
            if ( l .eq. 0 ) then
                call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
            endif

            ! Make a note of the multiplicity
            dj( ss%info%cellindex(1,j),ss%info%cellindex(2,j),ss%info%cellindex(3,j),a1,a2)=l
        enddo
        enddo

        ! Then we do it again and store some things
        allocate(vectormapping(9,k))
        vectormapping=0

        k=0
        do a1=1,p%na
        do j=1,ss%na

            v0=ss%rcart(:,j)-p%rcart(:,a1)
            a2=ss%info%index_in_unitcell(j)

            l=0
            do ii=-1,1
            do jj=-1,1
            do kk=-1,1
                v1=real([ii,jj,kk],r8)
                v1=matmul(ss%latticevectors,v1)
                v1=v1+v0
                if ( voro%is_point_inside(v1,1E-5_r8) ) then
                    k=k+1
                    ! Actual lattice vector?
                    ! v1 = v(j) - v(i)
                    ! v1 = lv(j) + r(a2) - r(a1)
                    ! lv(j) = v1 -r(a2)+r(a1)
                    v2=v1 - p%rcart(:,a2) + p%rcart(:,a1)
                    v2=matmul(p%inv_latticevectors,v2)
                    if ( norm2(v2-anint(v2)) .gt. 1E-6_r8 ) then
                        call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
                    else
                        gi=int(anint(v2))
                    endif

                    ! store actual lattice vector
                    vectormapping(1:3,k)=gi
                    ! store unitcell atom indices
                    vectormapping(4,k)=a1
                    vectormapping(5,k)=a2
                    ! store weight
                    vectormapping(6,k)=dj( ss%info%cellindex(1,j),ss%info%cellindex(2,j),ss%info%cellindex(3,j),a1,a2)
                    ! store Fourier indexing thing
                    vectormapping(7:9,k)=ss%info%cellindex(:,j)
                endif
            enddo
            enddo
            enddo
        enddo
        enddo
        deallocate(dj)

        ise%n_rvec=size(vectormapping,2)
        allocate(ise%rvec(3,ise%n_rvec))
        allocate(ise%atomind(2,ise%n_rvec))
        allocate(ise%rfc(3,3,ise%n_rvec,ise%n_energy))
        allocate(ise%ifc(3,3,ise%n_rvec,ise%n_energy))
        ise%rvec=0.0_r8
        ise%atomind=0
        ise%rfc=0.0_r8
        ise%ifc=0.0_r8
        do ir=1,ise%n_rvec
            ise%atomind(:,ir)=vectormapping(4:5,ir)
            ise%rvec(:,ir)=matmul(p%latticevectors,real(vectormapping(1:3,ir),r8))
        enddo

        allocate(dm0(p%na*3,p%na*3,griddensity(1),griddensity(2),griddensity(3)))
        allocate(dm1(p%na*3,p%na*3,griddensity(1),griddensity(2),griddensity(3)))
        allocate(rotmat(p%na*3,p%na*3))
        allocate(cm0(p%na*3,p%na*3))
        allocate(cm1(p%na*3,p%na*3))
        allocate(cm2(p%na*3,p%na*3))
        dm0=0.0_r8
        dm1=0.0_r8
        rotmat=0.0_r8
        cm0=0.0_r8
        cm1=0.0_r8
        cm2=0.0_r8
        do ie=1,ise%n_energy
            if ( mod(ie,mw%n) .ne. mw%r ) cycle
            ! Naive inverse Fourier transform
            dm0=0.0_r8
            dm1=0.0_r8
            do iq=1,ise%qp%n_full_point
                jq=ise%qp%ap(iq)%irreducible_index
                iop=ise%qp%ap(iq)%operation_from_irreducible
                if ( iop .gt. 0 ) then
                    call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(iop),inverseoperation=.false.)
                else
                    call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(-iop),inverseoperation=.true.)
                endif

                cm0=ise%sigma_Re(:,:,ie,jq)
                call lo_gemm(rotmat,cm0,cm2)
                call lo_gemm(cm2,rotmat,cm0,transb='C')

                cm1=ise%sigma_Im(:,:,ie,jq)
                call lo_gemm(rotmat,cm1,cm2)
                call lo_gemm(cm2,rotmat,cm1,transb='C')

                do i=1,griddensity(1)
                do j=1,griddensity(2)
                do k=1,griddensity(3)
                    v0=[i,j,k]-1.0_r8
                    v0=matmul(p%latticevectors,v0)
                    kdotr=-dot_product(v0,ise%qp%ap(iq)%r)*lo_twopi
                    expikr=cmplx(cos(kdotr),sin(kdotr),r8)
                    dm0(:,:,i,j,k)=dm0(:,:,i,j,k)+cm0*expikr
                    dm1(:,:,i,j,k)=dm1(:,:,i,j,k)+cm1*expikr
                enddo
                enddo
                enddo
            enddo

            ! Sort into something resonable
            do ir=1,ise%n_rvec
                a1=ise%atomind(1,ir)
                a2=ise%atomind(2,ir)
                weight=1.0_r8/real(vectormapping(6,ir),r8)
                ise%rfc(:,:,ir,ie)=weight*dm0((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,vectormapping(7,ir),vectormapping(8,ir),vectormapping(9,ir))
                ise%ifc(:,:,ir,ie)=weight*dm1((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,vectormapping(7,ir),vectormapping(8,ir),vectormapping(9,ir))
            enddo

            ! Why not test the interpolation?
            ! f0=0.0_r8
            ! f1=0.0_r8
            ! do iq=1,ise%qp%n_full_point
            !     jq=ise%qp%ap(iq)%irreducible_index
            !     iop=ise%qp%ap(iq)%operation_from_irreducible
            !     if ( iop .gt. 0 ) then
            !         call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(iop),inverseoperation=.false.)
            !     else
            !         call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(-iop),inverseoperation=.true.)
            !     endif
            !     cm0=ise%sigma_Im(:,:,ie,jq) !+ lo_imag*ise%sigma_Re(:,:,ie,jq)
            !     call lo_gemm(rotmat,cm0,cm1)
            !     call lo_gemm(cm1,rotmat,cm0,transb='C')

            !     cm1=0.0_r8
            !     do ir=1,ise%n_rvec
            !         a1=ise%atomind(1,ir)
            !         a2=ise%atomind(2,ir)
            !         kdotr=dot_product(ise%qp%ap(iq)%r,ise%rvec(:,ir))*lo_twopi
            !         expikr=cmplx(cos(kdotr),sin(kdotr),r8)
            !         cm1((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3)=&
            !         cm1((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3)+expikr*ise%cfc(:,:,ir,ie)
            !     enddo
            !     cm1=cm1/real(product(griddensity),r8)
            !     f0=f0+sum(abs(cm0-cm1))
            !     f1=f1+sum(abs(cm0))
            ! enddo
            !write(*,*) 'ie',ie,f0,f1,f0/f1
        enddo
        call mw%allreduce('sum',ise%rfc)
        call mw%allreduce('sum',ise%ifc)
        ise%rfc=ise%rfc/real(product(griddensity),r8)
        ise%ifc=ise%ifc/real(product(griddensity),r8)

        if ( verbosity .gt. 0 ) then
            write(*,*) 'Created naive Fourier interpolation'
        endif

    end block fourierinterpolation

    ! testinterpolation: block
    !     type(lo_phonon_dispersions_qpoint) :: ompoint
    !     real(r8), dimension(:,:), allocatable :: buf_re,buf_im
    !     complex(r8), dimension(:,:,:), allocatable :: bim
    !     integer :: iq

    !     if ( verbosity .gt. 0 ) then
    !         write(*,*) 'Testing interpolation:'
    !     endif

    !     allocate(buf_im(ise%n_energy,p%na*3))
    !     allocate(buf_re(ise%n_energy,p%na*3))
    !     allocate(bim(p%na*3,p%na*3,ise%n_energy))

    !     do iq=1,ise%qp%n_irr_point
    !         call ompoint%generate(fc,p,mem,ise%qp%ip(iq))
    !         call ise%evaluate(p,ise%qp%ip(iq)%r,ise%harm_omega(:,iq),ise%harm_egv(:,:,iq),buf_re,buf_im,mem)

    !         if ( verbosity .gt. 0 ) then
    !             write(*,*) ''
    !             write(*,*) 'iq',iq,matmul(p%inv_reciprocal_latticevectors,ise%qp%ip(iq)%r)
    !             write(*,*) 'nop',ise%qp%ip(iq)%n_invariant_operation
    !             write(*,*) 'im',sum(abs(buf_im)),sum(abs(buf_im-ise%sIm(:,:,iq)))
    !             write(*,*) 'Fourier',sum(abs(ise%ifc(:,:,1,:)-ise%sigma_im(:,:,:,iq)))
    !             ! do i=1,p%na*3
    !             !     write(*,*) ompoint%omega(i),dot_product(ise%harm_egv(:,i,iq),ompoint%egv(:,i))
    !             ! enddo
    !             !write(*,*) 'omega',ise%harm_omega(:,iq)
    !             !write(*,*) 'do',sum(abs(ise%harm_omega(:,iq)-ompoint%omega))
    !             ! do i=1,p%na*3
    !             ! do j=1,p%na*3
    !             !     f0=sum(abs(buf_im(:,i)-ise%sIm(:,j,iq)))
    !             !     if ( f0 .lt. 1E-10_r8 ) then
    !             !         write(*,*) i,j,sum(abs(buf_im(:,i)-ise%sIm(:,j,iq)))
    !             !     endif
    !             ! enddo
    !             ! enddo
    !         endif
    !     enddo
    ! end block testinterpolation

end subroutine

subroutine destroy_interpolated_selfenergy(ise)
    class(lo_interpolated_selfenergy_grid), intent(inout) :: ise

    call ise%box%destroy()
    if ( allocated(ise%qp) ) then
        call ise%qp%destroy(ise%qp)
    endif
    ise%n_energy=-lo_hugeint
    if ( allocated(ise%omega   ) ) deallocate(ise%omega   )
    if ( allocated(ise%sigma_Re) ) deallocate(ise%sigma_Re)
    if ( allocated(ise%sigma_Im) ) deallocate(ise%sigma_Im)
    call ise%aux_fc%destroy()
end subroutine

end module
