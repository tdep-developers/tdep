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
use type_blas_lapack_wrappers, only: lo_gemm, lo_dgels, lo_zheev, lo_gemv, lo_zgels
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_distancetable, only: lo_distancetable
use lo_voronoi, only: lo_voronoi_cell
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_forcemap, only: lo_forcemap

use lo_tetrahedron_interpolation, only: lo_linear_tetrahedron_interpolation
use type_phonon_dos, only: lo_phonon_dos
use lo_thermal_transport, only: lo_thermal_conductivity
use lo_spectralfunction_helpers, only: lo_evaluate_spectral_function,lo_gaussian_smear_spectral_function,lo_find_spectral_function_max_and_fwhm,lo_integrate_spectral_function,lo_tapering_function,lo_make_eigenvector_parallel,lo_optical_manifold
implicit none

private
public :: lo_interpolated_selfenergy_grid
!public :: lo_dynamical_matrix_coefficient_matrix_for_single_q

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

    ! ! harmonic omega, per q
    ! real(r8), dimension(:,:), allocatable :: harm_omega
    ! complex(r8), dimension(:,:,:), allocatable :: harm_egv
    ! real(r8), dimension(:,:,:), allocatable :: sIm,sRe

    ! ! Fourier interpolation thingies
    ! integer :: n_rvec=-lo_hugeint
    ! real(r8), dimension(:,:,:,:), allocatable :: bre,bim
    ! complex(r8), dimension(:,:,:,:), allocatable :: bc
    ! real(r8), dimension(:,:), allocatable :: rvec
    ! integer, dimension(:,:), allocatable :: atomind

    ! complex(r8), dimension(:,:,:,:), allocatable :: rfc
    ! complex(r8), dimension(:,:,:,:), allocatable :: ifc

    ! ! TDEP interpolation thingies
    ! type(lo_forcemap) :: map
    ! complex(r8), dimension(:,:), allocatable :: irr_re
    ! complex(r8), dimension(:,:), allocatable :: irr_im

    ! Weird optical manifold thing for interpolation
    real(r8), dimension(:,:), allocatable :: optical_manifold
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
        type(lo_hdf5_helper) :: h5
        integer :: iq

        if ( mw%talk  ) then
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

        if ( verbosity .gt. 0 ) then
            write(*,*) '... read self-energy'
        endif

        call h5%close_file()
        call h5%destroy()

        ! Generate triangulation thingy
        call ise%box%generate(ise%qp,p)

        ! Might need the optical manifold thing
        allocate(ise%optical_manifold(p%na*3,p%na*3))
        ise%optical_manifold=0.0_r8
        call lo_optical_manifold(fc,p,ise%optical_manifold)

        if ( verbosity .gt. 0 ) then
            write(*,*) 'Done reading self-energy from file'
        endif
    end block readfile

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

        ! ! Store the grid into a Fourier interpolation thingy?
        ! fourierinterpolation: block
        !     type(lo_crystalstructure) :: ss
        !     type(lo_distancetable) :: dt
        !     type(lo_voronoi_cell) :: voro

        !     complex(r8), dimension(:,:,:,:,:), allocatable :: dm0,dm1
        !     complex(r8), dimension(:,:), allocatable :: rotmat,cm0,cm1,cm2
        !     complex(r8) :: expikr
        !     real(r8), dimension(3,1) :: dummypos
        !     real(r8), dimension(3) :: v0,v1,v2
        !     real(r8) :: kdotr,f0,f1,weight
        !     integer, dimension(:,:,:,:,:), allocatable :: dj
        !     integer, dimension(:,:), allocatable :: vectormapping
        !     integer, dimension(3) :: griddensity,gi
        !     integer :: i,j,k,l,ii,jj,kk,ie,iq,jq,iop,a1,a2,ir

        !     ! Build supercell and Voronoi cell?
        !     select type(qp=>ise%qp)
        !     type is(lo_fft_mesh)
        !         griddensity=qp%griddensity
        !     end select
        !     call p%build_supercell(ss,dimensions=griddensity)
        !     call ss%classify('supercell',p)
        !     f0=lo_bounding_sphere_of_box(ss%latticevectors)*3
        !     dummypos=0.0_r8
        !     call dt%generate(dummypos,ss%latticevectors,f0,-1)
        !     call voro%generate(dt%particle(1),f0*2,1E-6_r8,mem)

        !     allocate(dj(griddensity(1),griddensity(2),griddensity(3),p%na,p%na))
        !     dj=0

        !     k=0
        !     do a1=1,p%na
        !     do j=1,ss%na
        !         v0=ss%rcart(:,j)-p%rcart(:,a1)
        !         a2=ss%info%index_in_unitcell(j)

        !         ! Sanity check
        !         v1=ss%info%cellindex(:,j)
        !         v1=matmul(p%latticevectors,(v1-1.0_r8)+p%r(:,a2))-p%rcart(:,a1)
        !         if ( norm2(v0-v1) .gt. 1E-6_r8 ) then
        !             call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
        !         endif

        !         l=0
        !         do ii=-3,3
        !         do jj=-3,3
        !         do kk=-3,3
        !             v1=real([ii,jj,kk],r8)
        !             v1=matmul(ss%latticevectors,v1)
        !             v1=v1+v0
        !             if ( voro%is_point_inside(v1,1E-5_r8) ) then
        !                 k=k+1
        !                 l=l+1
        !             endif
        !         enddo
        !         enddo
        !         enddo
        !         if ( l .eq. 0 ) then
        !             call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
        !         endif

        !         ! Make a note of the multiplicity
        !         dj( ss%info%cellindex(1,j),ss%info%cellindex(2,j),ss%info%cellindex(3,j),a1,a2)=l
        !     enddo
        !     enddo

        !     ! Then we do it again and store some things
        !     allocate(vectormapping(9,k))
        !     vectormapping=0

        !     k=0
        !     do a1=1,p%na
        !     do j=1,ss%na

        !         v0=ss%rcart(:,j)-p%rcart(:,a1)
        !         a2=ss%info%index_in_unitcell(j)

        !         l=0
        !         do ii=-1,1
        !         do jj=-1,1
        !         do kk=-1,1
        !             v1=real([ii,jj,kk],r8)
        !             v1=matmul(ss%latticevectors,v1)
        !             v1=v1+v0
        !             if ( voro%is_point_inside(v1,1E-5_r8) ) then
        !                 k=k+1
        !                 ! Actual lattice vector?
        !                 ! v1 = v(j) - v(i)
        !                 ! v1 = lv(j) + r(a2) - r(a1)
        !                 ! lv(j) = v1 -r(a2)+r(a1)
        !                 v2=v1 - p%rcart(:,a2) + p%rcart(:,a1)
        !                 v2=matmul(p%inv_latticevectors,v2)
        !                 if ( norm2(v2-anint(v2)) .gt. 1E-6_r8 ) then
        !                     call lo_stop_gracefully(['I do not understand vectors'],lo_exitcode_param,__FILE__,__LINE__)
        !                 else
        !                     gi=int(anint(v2))
        !                 endif

        !                 ! store actual lattice vector
        !                 vectormapping(1:3,k)=gi
        !                 ! store unitcell atom indices
        !                 vectormapping(4,k)=a1
        !                 vectormapping(5,k)=a2
        !                 ! store weight
        !                 vectormapping(6,k)=dj( ss%info%cellindex(1,j),ss%info%cellindex(2,j),ss%info%cellindex(3,j),a1,a2)
        !                 ! store Fourier indexing thing
        !                 vectormapping(7:9,k)=ss%info%cellindex(:,j)
        !             endif
        !         enddo
        !         enddo
        !         enddo
        !     enddo
        !     enddo
        !     deallocate(dj)

        !     ise%n_rvec=size(vectormapping,2)
        !     allocate(ise%rvec(3,ise%n_rvec))
        !     allocate(ise%atomind(2,ise%n_rvec))
        !     allocate(ise%rfc(3,3,ise%n_rvec,ise%n_energy))
        !     allocate(ise%ifc(3,3,ise%n_rvec,ise%n_energy))
        !     ise%rvec=0.0_r8
        !     ise%atomind=0
        !     ise%rfc=0.0_r8
        !     ise%ifc=0.0_r8
        !     do ir=1,ise%n_rvec
        !         ise%atomind(:,ir)=vectormapping(4:5,ir)
        !         ise%rvec(:,ir)=matmul(p%latticevectors,real(vectormapping(1:3,ir),r8))
        !     enddo

        !     allocate(dm0(p%na*3,p%na*3,griddensity(1),griddensity(2),griddensity(3)))
        !     allocate(dm1(p%na*3,p%na*3,griddensity(1),griddensity(2),griddensity(3)))
        !     allocate(rotmat(p%na*3,p%na*3))
        !     allocate(cm0(p%na*3,p%na*3))
        !     allocate(cm1(p%na*3,p%na*3))
        !     allocate(cm2(p%na*3,p%na*3))
        !     dm0=0.0_r8
        !     dm1=0.0_r8
        !     rotmat=0.0_r8
        !     cm0=0.0_r8
        !     cm1=0.0_r8
        !     cm2=0.0_r8
        !     do ie=1,ise%n_energy
        !         if ( mod(ie,mw%n) .ne. mw%r ) cycle
        !         ! Naive inverse Fourier transform
        !         dm0=0.0_r8
        !         dm1=0.0_r8
        !         do iq=1,ise%qp%n_full_point
        !             jq=ise%qp%ap(iq)%irreducible_index
        !             iop=ise%qp%ap(iq)%operation_from_irreducible
        !             if ( iop .gt. 0 ) then
        !                 call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(iop),inverseoperation=.false.)
        !             else
        !                 call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(-iop),inverseoperation=.true.)
        !             endif

        !             cm0=ise%sigma_Re(:,:,ie,jq)
        !             call lo_gemm(rotmat,cm0,cm2)
        !             call lo_gemm(cm2,rotmat,cm0,transb='C')

        !             cm1=ise%sigma_Im(:,:,ie,jq)
        !             call lo_gemm(rotmat,cm1,cm2)
        !             call lo_gemm(cm2,rotmat,cm1,transb='C')

        !             do i=1,griddensity(1)
        !             do j=1,griddensity(2)
        !             do k=1,griddensity(3)
        !                 v0=[i,j,k]-1.0_r8
        !                 v0=matmul(p%latticevectors,v0)
        !                 kdotr=-dot_product(v0,ise%qp%ap(iq)%r)*lo_twopi
        !                 expikr=cmplx(cos(kdotr),sin(kdotr),r8)
        !                 dm0(:,:,i,j,k)=dm0(:,:,i,j,k)+cm0*expikr
        !                 dm1(:,:,i,j,k)=dm1(:,:,i,j,k)+cm1*expikr
        !             enddo
        !             enddo
        !             enddo
        !         enddo

        !         ! Sort into something resonable
        !         do ir=1,ise%n_rvec
        !             a1=ise%atomind(1,ir)
        !             a2=ise%atomind(2,ir)
        !             weight=1.0_r8/real(vectormapping(6,ir),r8)
        !             ise%rfc(:,:,ir,ie)=weight*dm0((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,vectormapping(7,ir),vectormapping(8,ir),vectormapping(9,ir))
        !             ise%ifc(:,:,ir,ie)=weight*dm1((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,vectormapping(7,ir),vectormapping(8,ir),vectormapping(9,ir))
        !         enddo

        !         ! Why not test the interpolation?
        !         ! f0=0.0_r8
        !         ! f1=0.0_r8
        !         ! do iq=1,ise%qp%n_full_point
        !         !     jq=ise%qp%ap(iq)%irreducible_index
        !         !     iop=ise%qp%ap(iq)%operation_from_irreducible
        !         !     if ( iop .gt. 0 ) then
        !         !         call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(iop),inverseoperation=.false.)
        !         !     else
        !         !         call lo_eigenvector_transformation_matrix(rotmat,p%rcart,ise%qp%ip( jq )%r,p%sym%op(-iop),inverseoperation=.true.)
        !         !     endif
        !         !     cm0=ise%sigma_Im(:,:,ie,jq) !+ lo_imag*ise%sigma_Re(:,:,ie,jq)
        !         !     call lo_gemm(rotmat,cm0,cm1)
        !         !     call lo_gemm(cm1,rotmat,cm0,transb='C')

        !         !     cm1=0.0_r8
        !         !     do ir=1,ise%n_rvec
        !         !         a1=ise%atomind(1,ir)
        !         !         a2=ise%atomind(2,ir)
        !         !         kdotr=dot_product(ise%qp%ap(iq)%r,ise%rvec(:,ir))*lo_twopi
        !         !         expikr=cmplx(cos(kdotr),sin(kdotr),r8)
        !         !         cm1((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3)=&
        !         !         cm1((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3)+expikr*ise%cfc(:,:,ir,ie)
        !         !     enddo
        !         !     cm1=cm1/real(product(griddensity),r8)
        !         !     f0=f0+sum(abs(cm0-cm1))
        !         !     f1=f1+sum(abs(cm0))
        !         ! enddo
        !         !write(*,*) 'ie',ie,f0,f1,f0/f1
        !     enddo
        !     call mw%allreduce('sum',ise%rfc)
        !     call mw%allreduce('sum',ise%ifc)
        !     ise%rfc=ise%rfc/real(product(griddensity),r8)
        !     ise%ifc=ise%ifc/real(product(griddensity),r8)

        !     if ( verbosity .gt. 0 ) then
        !         write(*,*) 'Created naive Fourier interpolation'
        !     endif

        ! end block fourierinterpolation

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

    ! ! Slightly smarter Fourier interpolation, possibly
    ! tdepinterpolation: block
    !     type(lo_interaction_tensors) :: slt
    !     type(lo_crystalstructure) :: ss
    !     complex(r8), dimension(:,:,:), allocatable :: cfm
    !     complex(r8), dimension(:,:), allocatable :: ATA,cm0,cm1
    !     complex(r8), dimension(:), allocatable :: cv0,cv1,ATBr,ATBi,cw0,cw1
    !     real(r8) :: cutoff
    !     integer, dimension(3) :: griddensity
    !     integer :: iq,ie,a1,a2,i,j,k,ii,jj

    !     ! Build supercell
    !     select type(qp=>ise%qp)
    !     type is(lo_fft_mesh)
    !         griddensity=qp%griddensity
    !     end select
    !     call p%build_supercell(ss,dimensions=griddensity)
    !     call ss%classify('supercell',p)

    !     ! Sort out the symmetry things
    !     cutoff=min(60.0_r8,ss%maxcutoff()*0.7_r8)
    !     call slt%generate(p, ss, cutoff, -1.0_r8, -1.0_r8, .false., mw, mem, -1)
    !     call ise%map%generate(p, ss, polarcorrectiontype=0, st=slt, mw=mw, mem=mem, verbosity=-1)

    !     if ( mw%talk ) then
    !         write(*,*) 'post-symmetry:'
    !         write(*,*) 'npair:  ',ise%map%xuc%n_fc_pair
    !         write(*,*) 'nirr:   ',ise%map%xuc%nx_fc_pair
    !         write(*,*) 'cutoff: ',cutoff
    !     endif

    !     ! Now for the tricky part, I have to fit dispersions to the self-energy
    !     ! in dynamical matrix space. Should not be too bad I hope.
    !     allocate( cfm( (p%na*3)**2,ise%map%xuc%nx_fc_pair,ise%qp%n_irr_point ) )
    !     allocate( ATA(ise%map%xuc%nx_fc_pair,ise%map%xuc%nx_fc_pair) )
    !     cfm=0.0_r8
    !     ATA=0.0_r8
    !     do iq=1,ise%qp%n_irr_point
    !         if ( mod(iq,mw%n) .ne. mw%r ) cycle
    !         call lo_dynamical_matrix_coefficient_matrix_for_single_q(ise%map,ise%qp%ip(iq)%r,cfm(:,:,iq))
    !         call lo_gemm(cfm(:,:,iq),cfm(:,:,iq),ATA,transa='T',beta=cmplx(1.0_r8,0.0_r8,r8))
    !     enddo
    !     call mw%allreduce('sum',cfm)
    !     call mw%allreduce('sum',ATA)

    !     ! Space for solution
    !     allocate(ise%irr_re(ise%map%xuc%nx_fc_pair,ise%n_energy))
    !     allocate(ise%irr_im(ise%map%xuc%nx_fc_pair,ise%n_energy))
    !     ise%irr_re=0.0_r8
    !     ise%irr_im=0.0_r8

    !     allocate(cv0( (p%na*3)**2 ))
    !     allocate(cv1( (p%na*3)**2 ))
    !     allocate(cw0( (p%na*3)**2 ))
    !     allocate(cw1( (p%na*3)**2 ))
    !     allocate(ATBr( ise%map%xuc%nx_fc_pair ))
    !     allocate(ATBi( ise%map%xuc%nx_fc_pair ))
    !     allocate(cm0(ise%map%xuc%nx_fc_pair,ise%map%xuc%nx_fc_pair) )
    !     allocate(cm1(ise%map%xuc%nx_fc_pair,1) )
    !     cv0=0.0_r8
    !     cv1=0.0_r8
    !     cw0=0.0_r8
    !     cw1=0.0_r8

    !     ATBr=0.0_r8
    !     ATBi=0.0_r8
    !     cm0=0.0_r8
    !     cm1=0.0_r8

    !     ! Ok, decent start.
    !     do ie=1,ise%n_energy

    !         if ( mod(ie,mw%n) .ne. mw%r ) cycle

    !         ATBr=0.0_r8
    !         ATBi=0.0_r8
    !         do iq=1,ise%qp%n_irr_point
    !             k=0
    !             do a1 = 1, p%na
    !             do a2 = 1, p%na
    !                 do i = 1, 3
    !                 do j = 1, 3
    !                     ii=(a1-1)*3 + i
    !                     jj=(a2-1)*3 + j
    !                     k=k+1
    !                     cv0(k) = ise%sigma_Re(jj,ii,ie,iq)
    !                     cv1(k) = ise%sigma_Im(jj,ii,ie,iq)
    !                 end do
    !                 end do
    !             end do
    !             end do
    !             ! ATBr=ATBr + matmul(transpose(cfm(:,:,iq)),cv0 )
    !             ! ATBi=ATBi + matmul(transpose(cfm(:,:,iq)),cv1 )
    !             call lo_gemv(cfm(:,:,iq),cv0,ATBr,trans='T',beta=cmplx(1.0_r8,0.0_r8,r8))
    !             call lo_gemv(cfm(:,:,iq),cv1,ATBi,trans='T',beta=cmplx(1.0_r8,0.0_r8,r8))
    !         enddo
    !         ! Solve linear system - complex I guess?
    !         cm0=ATA
    !         cm1(:,1)=ATBr + ATBi
    !         call lo_zgels(cm0,cm1)
    !         ise%irr_re(:,ie)=cm1(:,1)

    !         cm0=ATA
    !         cm1(:,1)=ATBi + ATBr
    !         call lo_zgels(cm0,cm1)
    !         ise%irr_im(:,ie)=cm1(:,1)

    !         ! if ( mw%talk ) then
    !         !     write(*,*) 're'
    !         !     do i=1,ise%map%xuc%nx_fc_pair
    !         !         write(*,*) i,ise%irr_re(i,ie)
    !         !     enddo
    !         !     write(*,*) 'im'
    !         !     do i=1,ise%map%xuc%nx_fc_pair
    !         !         write(*,*) i,ise%irr_im(i,ie)
    !         !     enddo
    !         ! endif

    !         ! See if we can reconstruct?
    !         do iq=1,ise%qp%n_irr_point

    !             k=0
    !             do a1 = 1, p%na
    !             do a2 = 1, p%na
    !                 do i = 1, 3
    !                 do j = 1, 3
    !                     ii=(a1-1)*3 + i
    !                     jj=(a2-1)*3 + j
    !                     k=k+1
    !                     cv0(k) = ise%sigma_Re(jj,ii,ie,iq)
    !                     cv1(k) = ise%sigma_Im(jj,ii,ie,iq)
    !                     !coefficientmatrix(jj, :) = lo_chop(Ck(i, j, a1, a2, :), lo_sqtol)
    !                 end do
    !                 end do
    !             end do
    !             end do
    !             cw0=matmul(cfm(:,:,iq),ise%irr_re(:,ie))
    !             cw1=matmul(cfm(:,:,iq),ise%irr_im(:,ie))
    !             if ( mw%talk ) then
    !                 !write(*,*) 'iq',sum(abs(cv0-cw0)),sum(abs(cv1-cw1)),sum(abs(cv0)),sum(abs(cv1))
    !                 !write(*,*) 'iq',sum(abs(cv0-cw0))/sum(abs(cv0)),sum(abs(cv1-cw1))/sum(abs(cv1)),sum(abs(cv0)),sum(abs(cv1))
    !                 ! do i=1,size(cv0)
    !                 !     write(*,*) i,cv0(i),cw0(i)
    !                 ! enddo
    !             endif
    !         enddo
    !     enddo
    !     call mw%allreduce('sum',ise%irr_im)
    !     call mw%allreduce('sum',ise%irr_re)

    !     if ( mw%talk ) then
    !         write(*,*) '... created TDEP fit'
    !         write(*,*) 're',sum(abs(real(ise%irr_re,r8))),sum(abs(aimag(ise%irr_re)))
    !         write(*,*) 'im',sum(abs(real(ise%irr_im,r8))),sum(abs(aimag(ise%irr_im)))
    !     endif

    !     ! Now do a fit at all energies?
    ! end block tdepinterpolation

! !> Construct the dynamical matrix coefficient matrix for a specific q-point
! subroutine lo_dynamical_matrix_coefficient_matrix_for_single_q(map, qv, coefficientmatrix, uc)
!     !> forcemap
!     type(lo_forcemap), intent(in) :: map
!     !> q-vector
!     real(r8), dimension(3), intent(in) :: qv
!     !> coefficient matrix
!     complex(r8), dimension(:, :), intent(out) :: coefficientmatrix
!     !> unitcell, in case I want the masses in there
!     type(lo_crystalstructure), intent(in), optional :: uc

!     complex(r8), dimension(:, :, :, :, :), allocatable :: Ck
!     complex(r8), dimension(9, 9) :: C1
!     complex(r8) :: expiqr
!     real(r8) :: k_dot_r
!     integer :: i, j, k, l, ii, jj, sh, o, a1, a2, ipair
!     integer :: nx, na, nfc

!     ! Size of things
!     na = map%n_atom_uc         ! number of atoms
!     nx = map%xuc%nx_fc_pair    ! dimensions of irreducible IFC

!     if (size(coefficientmatrix, 1) .ne. 3*3*na*na ) then
!         write (*, *) 'bad dimensions in dynmatrixcoeffM'
!         stop
!     end if
!     if (size(coefficientmatrix, 2) .ne. nx) then
!         write (*, *) 'bad dimensions in dynmatrixcoeffM'
!         stop
!     end if

!     allocate (Ck(3, 3, na, na, nx))
!     Ck = 0.0_r8
!     do ipair = 1, map%xuc%n_fc_pair
!         !write(*,*) 'DANGER INDEX CHECK THIS',__LINE__,__FILE__
!         a1 = map%xuc%fc_pair(ipair)%i1
!         a2 = map%xuc%fc_pair(ipair)%i2
!         sh = map%xuc%fc_pair(ipair)%irreducible_shell
!         o = map%xuc%fc_pair(ipair)%operation_from_shell
!         nfc = map%fc_pair_shell(sh)%nx
!         if (nfc .eq. 0) cycle
!         ! Get the Fourier transform thingy
!         k_dot_r = dot_product(map%xuc%fc_pair(ipair)%lv, qv)*lo_twopi
!         expiqr = cmplx(cos(k_dot_r), sin(k_dot_r), r8)
!         C1 = 0.0_r8
!         C1(:, 1:nfc) = matmul(map%op_pair(o)%sotr, map%fc_pair_shell(sh)%coeff)
!         ! Add to the self-term
!         do k = 1, nfc
!             l = map%fc_pair_shell(sh)%ind_global(k)
!             do i = 1, 3
!             do j = 1, 3
!                 ii = (i - 1)*3 + j
!                 Ck(i, j, a1, a1, l) = Ck(i, j, a1, a1, l) - C1(ii, k)
!             end do
!             end do
!         end do
!         ! not the self-term
!         C1 = C1*expiqr
!         do k = 1, nfc
!             l = map%fc_pair_shell(sh)%ind_global(k)
!             do i = 1, 3
!             do j = 1, 3
!                 ii = (i - 1)*3 + j
!                 Ck(i, j, a1, a2, l) = Ck(i, j, a1, a2, l) + C1(ii, k)
!             end do
!             end do
!         end do
!     end do

!     ! Put this in the right place
!     coefficientmatrix = 0.0_r8
!     jj = 0
!     do a1 = 1, na
!     do a2 = 1, na
!         do i = 1, 3
!         do j = 1, 3
!             jj = jj + 1
!             if (present(uc)) then
!                 coefficientmatrix(jj, :) = lo_chop(Ck(i, j, a1, a2, :), lo_sqtol)*uc%invsqrtmass(a1)*uc%invsqrtmass(a2)
!             else
!                 coefficientmatrix(jj, :) = lo_chop(Ck(i, j, a1, a2, :), lo_sqtol)
!             end if
!         end do
!         end do
!     end do
!     end do

!     deallocate (Ck)
! end subroutine

end module
