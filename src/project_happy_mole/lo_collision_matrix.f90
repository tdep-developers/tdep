module lo_collision_matrix
!! Container for things related to thermal transport. Not how it is actually calculated,
!! that is generated elsewhere, but how it is stored and written to file.
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_freqtol, lo_pi, lo_twopi, lo_exitcode_param, &
                      lo_groupvel_Hartreebohr_to_ms, lo_kb_Hartree, lo_kappa_au_to_SI, lo_frequency_Hartree_to_THz, &
                      lo_frequency_Hartree_to_icm, lo_frequency_Hartree_to_meV, lo_time_au_to_s, lo_bohr_to_A
use gottochblandat, only: lo_trapezoid_integration, lo_clean_fractional_coordinates, lo_chop
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh,lo_fft_mesh,lo_wedge_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions,lo_phonon_dispersions_qpoint
use type_symmetryoperation, only: lo_expandoperation_pair, lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm
use quadratures_stencils, only: lo_centraldifference

use lo_selfenergy_interpolation, only: lo_interpolated_selfenergy_grid
use lo_spectralfunction_convolution, only: lo_convolution_handle
use lo_spectralfunction_helpers, only: lo_evaluate_spectral_function,lo_gaussian_smear_spectral_function,lo_find_spectral_function_max_and_fwhm,lo_integrate_spectral_function
use lo_distributed_phonon_dispersion_relations, only: lo_distributed_phonon_dispersions_qpoint
use lo_verletboxes, only: lo_verletbox
implicit none

private
public :: lo_scattering_matrix

type lo_scattering_matrix
    contains
        procedure :: generate=>create_scattering_matrix
end type

type lo_scattering_matrix_convolutions
end type

contains

subroutine create_scattering_matrix(scm,p,fc2,fc3,ise,temperature,qp,adaptive_prefactor,mw,mem,verbosity)
    !> scattering matrix
    class(lo_scattering_matrix), intent(inout) :: scm
    !> structure
    type(lo_crystalstructure), intent(inout) :: p
    !> second order forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc2
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fc3
    !> interpolated self-energy
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> temperature
    real(r8), intent(in) :: temperature
    !> q-point grid
    class(lo_qpoint_mesh), intent(inout) :: qp
    !> smearing prefactor
    real(r8), intent(in) :: adaptive_prefactor
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_phonon_dispersions) :: dr
    type(lo_convolution_handle) :: ch
    integer, dimension(:,:), allocatable :: qmesh_permutation
    integer, dimension(:), allocatable :: qmesh_permutation_sign

    init: block
        if ( mw%talk ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'Creating scattering matrix'
        endif

        ! Going to need harmonic dispersions on the mesh
        call dr%generate(qp,fc2,p,mw,mem,verbosity=-1)

        ! And a convolution helper
        call ch%generate(ise%omega,temperature,dr%n_mode)

        if ( mw%talk ) then
            write(lo_iou,*) '   n irreducible q-points:',qp%n_irr_point
            write(lo_iou,*) '               n q-points:',qp%n_full_point
        endif

        ! We will eventually need the symmetry permutation matrices
        call qmesh_permutation_arrays(qp,p,qmesh_permutation,qmesh_permutation_sign,mw,mem)
    end block init

    buildmatrix: block
        type(lo_hdf5_helper) :: h5

        real(r8), dimension(:,:,:,:), allocatable :: halfmatrix
        real(r8), dimension(:,:,:,:), allocatable :: fullmatrix,fullmatrix2
        real(r8), dimension(:,:), allocatable :: bar_bubble,bar_transform
        real(r8), dimension(:,:), allocatable :: submatrix
        real(r8), dimension(:), allocatable :: dr0,dr1
        integer :: i,j,iq,jq,kq,lq,iop

        allocate(submatrix(dr%n_mode,dr%n_mode))
        allocate(halfmatrix(dr%n_mode,dr%n_mode,qp%n_irr_point,qp%n_full_point))
        submatrix=0.0_r8
        halfmatrix=0.0_r8

        fullqploop: do jq=1,qp%n_full_point
            if ( mod(jq,mw%n) .ne. mw%r ) cycle
            irrqploop: do iq=1,qp%n_irr_point
                if ( mw%talk ) write(*,*) 'submatrix',iq,jq
                call cubic_scattering_matrix_entry(iq,jq,p,fc2,fc3,qp,dr,ise,ch,adaptive_prefactor,submatrix,mem)
                halfmatrix(:,:,iq,jq)=submatrix
            enddo irrqploop
        enddo fullqploop
        call mw%allreduce('sum',halfmatrix)

        ! So, let's see if we can expand it into a full matrix.
        allocate(fullmatrix(dr%n_mode,dr%n_mode,qp%n_full_point,qp%n_full_point))
        fullmatrix=0.0_r8

        do iop=1,size(qmesh_permutation,2)
            ! Uh. Not sure how to think here.
            do iq=1,qp%n_irr_point
                ! what full point is this
                jq=qp%ip(iq)%index_full_point(1)
                ! then we permute with operation
                kq=qmesh_permutation(jq,iop)
                ! if it's not the same we full out the row?
                if ( iq .eq. kq ) cycle

                do i=1,qp%n_full_point
                    lq=qmesh_permutation(i,iop)
                    fullmatrix(:,:,kq,lq) = halfmatrix(:,:,iq,i)
                enddo
            enddo
        enddo

        allocate(fullmatrix2(dr%n_mode,dr%n_mode,qp%n_full_point,qp%n_full_point))
        fullmatrix2=0.0_r8
        do iop=1,size(qmesh_permutation,2)
            do i=1,qp%n_full_point
            do j=1,qp%n_full_point
                iq=qmesh_permutation(i,iop)
                jq=qmesh_permutation(j,iop)
                fullmatrix2(:,:,iq,jq)=fullmatrix2(:,:,iq,jq)+fullmatrix(:,:,i,j)
            enddo
            enddo
            if ( mw%talk ) then
                write(*,*) 'iop',iop,sum(abs(fullmatrix-fullmatrix2/real(iop,r8)))/sum(abs(fullmatrix))
            endif
        enddo
        fullmatrix=fullmatrix2/size(qmesh_permutation,2)

        ! Then I guess we might need the static bubble as well.
        allocate(bar_bubble(dr%n_mode,qp%n_full_point))
        allocate(bar_transform(dr%n_mode,qp%n_full_point))
        allocate(dr0(dr%n_mode))
        allocate(dr1(dr%n_mode))
        bar_bubble=0.0_r8
        bar_transform=0.0_r8
        do iq=1,qp%n_irr_point
            if ( mod(iq,mw%n) .ne. mw%r ) cycle
            call averaged_bubble_entry(iq,p,qp,dr,ise,ch,temperature,adaptive_prefactor,dr0,dr1,mem)
            do i=1,qp%ip(iq)%n_full_point
                jq=qp%ip(iq)%index_full_point(i)
                bar_bubble(:,jq)=dr0
                bar_transform(:,jq)=dr1
            enddo
        enddo
        call mw%allreduce('sum',bar_bubble)
        call mw%allreduce('sum',bar_transform)

        ! Dump to file for examination
        if ( mw%talk ) then
            !call h5%init(__FILE__, __LINE__)
            call h5%open_file('write', 'outfile.collision_matrix.hdf5')

            call h5%store_data(fullmatrix,h5%file_id,'collision_matrix')
            call h5%store_data(bar_bubble,h5%file_id,'bar_bubble')
            call h5%store_data(bar_transform,h5%file_id,'bar_transform')

            call h5%close_file()
            !call h5%destroy(__FILE__, __LINE__)
        endif
    end block buildmatrix

    !call ise%evaluate()

end subroutine

!> calculate a specific entry in the scattering matrix
subroutine cubic_scattering_matrix_entry(iq,jq,p,fc,fct,qp,dr,ise,ch,adaptive_prefactor,submatrix,mem)
    !> index of first q-point, irreducible
    integer, intent(in) :: iq
    !> index of second q-point, full
    integer, intent(in) :: jq
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-point grid
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> interpolated self-energy
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> convolution handle
    type(lo_convolution_handle), intent(inout) :: ch
    !> adaptive smearing prefactor
    real(r8), intent(in) :: adaptive_prefactor
    !> actual scattering matrix
    real(r8), dimension(:,:), intent(out) :: submatrix
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_distributed_phonon_dispersions_qpoint) :: op3
    complex(r8), dimension(:,:), allocatable :: nuv1,nuv2
    real(r8), dimension(3) :: qv1,qv2,qv3
    integer :: kq

    init: block
        real(r8), dimension(3) :: v0
        real(r8) :: f0,f1
        integer, dimension(3) :: gi
        integer :: imode,iatom,ialpha,ii

        ! First step is to prep all harmonic properties.
        allocate(nuv1(p%na*3,p%na*3))
        allocate(nuv2(p%na*3,p%na*3))
        nuv1=0.0_r8
        nuv2=0.0_r8

        qv1=qp%ip(iq)%r
        qv2=qp%ap(jq)%r
        qv3=-qv1-qv2
        qv3=qv3 + p%bz%gshift(qv3)
        call op3%generate(fc,p,mem,qvec=qv3)
        do imode=1,p%na*3
            op3%sigma(imode) = qp%adaptive_sigma( op3%vel(:,imode),dr%default_smearing(imode),adaptive_prefactor)
        enddo

        select type(qp)
        type is(lo_fft_mesh)
            ! Fetch the third q-point from the mesh
            v0=-qv1-qv2
            v0=matmul(p%inv_reciprocal_latticevectors,v0)
            gi=qp%index_from_coordinate(v0)
            kq=qp%gridind2ind(gi(1),gi(2),gi(3))
            qv3=qp%ap(kq)%r
            ! Replace eigenvectors/frequencies with those from the grid to stay consistent in gauge
            op3%omega = dr%aq(kq)%omega
            op3%vel = dr%aq(kq)%vel
            op3%egv = dr%aq(kq)%egv
        type is(lo_wedge_mesh)
            v0=-qv1-qv2
            qv3=v0 + p%bz%gshift(v0)

            call lo_stop_gracefully(['FIXME WEDGE MESH'],0,__FILE__,__LINE__)
        end select

        ! Pre-calculate the nu-vectors
        do imode=1,p%na*3
            if ( dr%iq(iq)%omega(imode) .gt. lo_freqtol ) then
                f0=1.0_r8/sqrt(2.0_r8*dr%iq(iq)%omega(imode) )
            else
                f0=0.0_r8
            endif
            do iatom=1,p%na
                f1=p%invsqrtmass(iatom)
                do ialpha=1,3
                    ii=(iatom-1)*3+ialpha
                    nuv1(ii,imode)=dr%iq(iq)%egv(ii,imode)*f0*f1
                enddo
            enddo

            if ( dr%aq(jq)%omega(imode) .gt. lo_freqtol ) then
                f0=1.0_r8/sqrt(2.0_r8*dr%aq(jq)%omega(imode) )
            else
                f0=0.0_r8
            endif
            do iatom=1,p%na
                f1=p%invsqrtmass(iatom)
                do ialpha=1,3
                    ii=(iatom-1)*3+ialpha
                    nuv2(ii,imode)=dr%aq(jq)%egv(ii,imode)*f0*f1
                enddo
            enddo

            if ( op3%omega(imode) .gt. lo_freqtol ) then
                f0=1.0_r8/sqrt(2.0_r8*op3%omega(imode) )
            else
                f0=0.0_r8
            endif
            do iatom=1,p%na
                f1=p%invsqrtmass(iatom)
                do ialpha=1,3
                    ii=(iatom-1)*3+ialpha
                    op3%nuvec(ii,imode)=op3%egv(ii,imode)*f0*f1
                enddo
            enddo
        enddo
    end block init

    interpolate: block
        complex(r8), dimension(:), allocatable :: ptf_phi,evp1,evp2
        complex(r8) :: c0
        real(r8), dimension(:,:,:), allocatable :: psisq_3ph,buf_element
        real(r8), dimension(:,:,:), allocatable :: buf_gg,buf_ll
        real(r8), dimension(:,:), allocatable :: buf_re,buf_im,buf_j
        real(r8), dimension(:), allocatable :: buf_integral
        real(r8) :: sigma,f0,f1,pref
        integer :: imode,jmode,b1,b2,b3

        allocate(psisq_3ph(dr%n_mode,dr%n_mode,dr%n_mode))
        allocate(buf_element(dr%n_mode,dr%n_mode,dr%n_mode))
        allocate(evp1(dr%n_mode**2))
        allocate(evp2(dr%n_mode**3))
        allocate(ptf_phi(dr%n_mode**3))
        allocate(buf_re(ise%n_energy,dr%n_mode))
        allocate(buf_im(ise%n_energy,dr%n_mode))
        allocate(buf_j(ise%n_energy,dr%n_mode))
        allocate(buf_gg(ise%n_energy,dr%n_mode,dr%n_mode))
        allocate(buf_ll(ise%n_energy,dr%n_mode,dr%n_mode))
        allocate(buf_integral(ise%n_energy))
        buf_re=0.0_r8
        buf_im=0.0_r8
        buf_j=0.0_r8
        buf_gg=0.0_r8
        buf_ll=0.0_r8
        psisq_3ph=0.0_r8
        buf_element=0.0_r8
        evp1=0.0_r8
        evp2=0.0_r8
        ptf_phi=0.0_r8
        buf_integral=0.0_r8

        ! Before actually integrating we need the matrix elements.
        call pretransform_phi(fct, qv2, qv3, ptf_phi)
        psisq_3ph=0.0_r8
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), nuv2(:,b2), 1, nuv1(:, b1), 1, evp1, dr%n_mode)
            do b3 = 1, dr%n_mode
                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), op3%nuvec(:, b3), 1, evp1, 1, evp2, dr%n_mode)
                evp2 = conjg(evp2)
                c0=dot_product(evp2, ptf_phi)
                psisq_3ph(b1, b2, b3)=abs(conjg(c0)*c0)
            end do
        enddo
        enddo
        ! Degeneracies are extremely annoying, but has to be dealt with.
        call cubic_degeneracy_fold_in_fold_out(dr%iq(iq)%omega,dr%aq(jq)%omega,op3%omega,psisq_3ph,lo_freqtol)

        ! Start with q',s', evaluate, smear and normalize spectral function
        call ise%evaluate(p,qv2,dr%aq(jq)%omega,dr%aq(jq)%egv,buf_re,buf_im,mem)
        buf_j=0.0_r8
        do imode=1,dr%n_mode
            if ( dr%aq(jq)%omega(imode) .lt. lo_freqtol ) cycle
            call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),dr%aq(jq)%omega(imode),buf_j(:,imode))
            sigma=qp%adaptive_sigma( dr%aq(jq)%vel(:,imode), dr%default_smearing(imode), adaptive_prefactor)
            call lo_gaussian_smear_spectral_function(ise%omega,sigma,buf_j(:,imode))
            f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
            buf_j(:,imode)=buf_j(:,imode)/f0
        enddo
        ! Fix degeneracies to be on the safe side
        call spectrum_degeneracy_fold_in_fold_out(dr%aq(jq)%omega,buf_j,lo_freqtol)
        call ch%buffer_and_transform_J(buf_j,2)

        ! Evaluate q'',s'' just as above
        call ise%evaluate(p,qv3,op3%omega,op3%egv,buf_re,buf_im,mem)
        buf_j=0.0_r8
        do imode=1,dr%n_mode
            if ( op3%omega(imode) .lt. lo_freqtol ) cycle
            call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),op3%omega(imode),buf_j(:,imode))
            sigma=qp%adaptive_sigma( op3%vel(:,imode), dr%default_smearing(imode), adaptive_prefactor)
            call lo_gaussian_smear_spectral_function(ise%omega,sigma,buf_j(:,imode))
            f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
            buf_j(:,imode)=buf_j(:,imode)/f0
        enddo
        call spectrum_degeneracy_fold_in_fold_out(op3%omega,buf_j,lo_freqtol)
        call ch%buffer_and_transform_J_to_cubic_kernel(buf_j)

        ! So, at this point
        ! ch%lesser2 holds J(q') in the time domain
        ! ch%cubic holds J(q'')*thermal prefactor in the time domain
        ! time to convolute and inverse transform the combinations
        buf_gg=0.0_r8
        do imode=1,dr%n_mode
        do jmode=1,dr%n_mode
            ch%cbuf = ch%lesser2(:,imode)*ch%cubic(:,jmode)
            call ch%inverse_transform_to_real(ch%cbuf,buf_gg(:,imode,jmode))
        enddo
        enddo

        ! Now we can fetch J_{qs}
        call ise%evaluate(p,qp%ip(iq)%r,dr%iq(iq)%omega,dr%iq(iq)%egv,buf_re,buf_im,mem)
        buf_j=0.0_r8
        do imode=1,dr%n_mode
            if ( dr%iq(iq)%omega(imode) .lt. lo_freqtol ) cycle
            call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),dr%iq(iq)%omega(imode),buf_j(:,imode))
            sigma=qp%adaptive_sigma( dr%iq(iq)%vel(:,imode), dr%default_smearing(imode), adaptive_prefactor)
            call lo_gaussian_smear_spectral_function(ise%omega,sigma,buf_j(:,imode))
            f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
            buf_j(:,imode)=buf_j(:,imode)/f0
        enddo
        call spectrum_degeneracy_fold_in_fold_out(dr%iq(iq)%omega,buf_j,lo_freqtol)

        ! All the integrals. Can be flipped around so that it's only N^2 integrals, but
        ! I want to get things right first.
        buf_element=0.0_r8
        do b1=1,dr%n_mode
        do b2=1,dr%n_mode
        do b3=1,dr%n_mode
            buf_integral=buf_j(:,b1)*( buf_gg(:,b2,b3) )*psisq_3ph(b1,b2,b3)
            f0=lo_trapezoid_integration(ise%omega,buf_integral)
            buf_element(b1,b2,b3)=f0
        enddo
        enddo
        enddo

        !call cubic_degeneracy_fold_in_fold_out(dr%iq(iq)%omega,dr%aq(jq)%omega,op3%omega,buf_element,lo_freqtol)

        pref=1.0_r8
        submatrix=0.0_r8
        do b1=1,dr%n_mode
        do b2=1,dr%n_mode
            f0=0.0_r8
            do b3=1,dr%n_mode
                f0=f0+buf_element(b1,b2,b3)
            enddo
            submatrix(b1,b2) = f0*pref
        enddo
        enddo

        deallocate(psisq_3ph)
        deallocate(evp1)
        deallocate(evp2)
        deallocate(ptf_phi)
        deallocate(buf_re)
        deallocate(buf_im)
        deallocate(buf_j)
        deallocate(buf_gg)
        deallocate(buf_ll)
        deallocate(buf_integral)
    end block interpolate

end subroutine

!> pre-transform three-phonon matrix element to xyz
subroutine pretransform_phi(fct, q2, q3, ptf)
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2, q3
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3
    real(r8) :: iqr
    integer :: a1, a2, a3, ia, ib, ic, t, nb

    nb = fct%na*3
    ptf = 0.0_r8
    do a1 = 1, fct%na
    do t = 1, fct%atom(a1)%n
        a2 = fct%atom(a1)%triplet(t)%i2
        a3 = fct%atom(a1)%triplet(t)%i3

        rv2 = fct%atom(a1)%triplet(t)%lv2
        rv3 = fct%atom(a1)%triplet(t)%lv3

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            l = (ia - 1)*nb*nb + (ib - 1)*nb + ic
            ptf(l) = ptf(l) + fct%atom(a1)%triplet(t)%m(i, j, k)*expiqr
        end do
        end do
        end do
    end do
    end do
end subroutine

!> permutation arrays from q-point mesh
subroutine qmesh_permutation_arrays(qp,p,permutation,reversal,mw,mem)
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-mesh permutations
    integer, dimension(:,:), allocatable, intent(out) :: permutation
    !> does permutation include time reversal
    integer, dimension(:), allocatable, intent(out) :: reversal
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem



    real(r8), parameter :: tol=1E-10_r8
    type(lo_verletbox) :: vb
    real(r8), dimension(:,:), allocatable :: r0
    real(r8), dimension(3) :: v0,v1
    real(r8) :: f0
    integer, dimension(:,:), allocatable :: buf_perm
    integer, dimension(:), allocatable :: di,buf_sign
    integer, dimension(3) :: boxdim
    integer :: i,j,k,l,iop,ctr,bi,bj,bk,ii,jj,kk,pm

    if ( mw%talk ) then
        write(lo_iou,*) '... preparing q-mesh permutation arrays'
    endif

    allocate(r0(3,qp%n_full_point))
    r0=0.0_r8
    select type(qp)
    type is(lo_fft_mesh)
        do i=1,qp%n_full_point
            v0=qp%ap(i)%r
            v0=matmul(p%inv_reciprocal_latticevectors,v0)
            v0=lo_clean_fractional_coordinates(v0)
            r0(:,i)=v0
        enddo
    type is(lo_wedge_mesh)
        do i=1,qp%n_full_point
            r0(:,i)=qp%ap(i)%r
        enddo
    class default
        call lo_stop_gracefully(['Unknown mesh type'],0,__FILE__,__LINE__)
    end select

    ! Sort points into verlet boxes
    boxdim=vb%boxdim(r0,min(qp%n_full_point,6000),1E-8_r8)
    boxdim=max(boxdim,4)
    call vb%generate(r0,boxdim,mem)

    allocate(buf_perm(qp%n_full_point,p%sym%n*2))
    allocate(buf_sign(p%sym%n*2))
    buf_perm=-1
    buf_sign=0
    ctr=0
    do pm=1,2 ! loop over +-, i.e. add time reversal symmetry
    do iop=1,p%sym%n
        ctr=ctr+1

        ! Note wether this is time reversal or not
        buf_sign(ctr)=pm

        do i=1,qp%n_full_point
            ! make it parallel
            if ( mod(i,mw%n) .ne. mw%r ) cycle

            ! Rotate the array
            select type(qp)
            type is(lo_fft_mesh)
                v0=r0(:,i)
                select case(pm)
                case(1)
                    v0=matmul(p%sym%op(iop)%rfm,v0)
                case(2)
                    v0=-matmul(p%sym%op(iop)%rfm,v0)
                end select
                v0=lo_clean_fractional_coordinates(v0)
            type is(lo_wedge_mesh)
                v0=r0(:,i)
                select case(pm)
                case(1)
                    v0=matmul(p%sym%op(iop)%m,v0)
                case(2)
                    v0=-matmul(p%sym%op(iop)%m,v0)
                end select
            class default
                call lo_stop_gracefully(['Unknown mesh type'],0,__FILE__,__LINE__)
            end select

            ! Check the central box first
            call vb%boxind(v0,bi,bj,bk)
            l=-1
            do j=1,vb%box(bi,bj,bk)%n
                k=vb%box(bi,bj,bk)%ind(j)
                v1=v0 - r0(:,k)
                f0=v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
                if ( f0 .lt. tol ) then
                    l=k
                    exit
                endif
            enddo
            ! If not there, check adjacent. This is a corner case
            ! and quite unlikely but still needs to happen.
            if ( l .lt. 0 ) then
                b1l: do ii=max(bi-1,1),min(bi+1,vb%nx)
                do jj=max(bj-1,1),min(bj+1,vb%ny)
                do kk=max(bk-1,1),min(bk+1,vb%nz)
                    do j=1,vb%box(ii,jj,kk)%n
                        k=vb%box(ii,jj,kk)%ind(j)
                        v1=v0 - r0(:,k)
                        f0=v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
                        if ( f0 .lt. tol ) then
                            l=k
                            exit b1l
                        endif
                    enddo
                enddo
                enddo
                enddo b1l
            endif

            ! Now we really should have an assignment.
            if ( l .gt. 0 ) then
                buf_perm(i,ctr )=l
            else
                call lo_stop_gracefully(['Mesh not compliant with symmetry'],0,__FILE__,__LINE__)
            endif
        enddo
    enddo
    enddo

    ! Sync across ranks
    call mw%allreduce('max',buf_perm)

    ! It only makes sense to return the unique? Or?
    allocate(di(p%sym%n*2))
    di=1
    do i=1,p%sym%n*2
        if ( di(i) .eq. 0 ) cycle
        do j=i+1,p%sym%n*2
            k=sum(abs(buf_perm(:,j)-buf_perm(:,i)))
            if ( k .eq. 0 ) di(j)=0
        enddo
    enddo

    ! Return only the unique
    i=sum(di)
    allocate(permutation(qp%n_full_point,i))
    allocate(reversal(i))
    permutation=0
    reversal=0
    j=0
    do i=1,p%sym%n*2
        if ( di(i) .eq. 0 ) cycle
        j=j+1
        permutation(:,j)=buf_perm(:,i)
        reversal(j)=buf_sign(i)
    enddo

    call vb%destroy(mem)
    deallocate(r0)
    deallocate(buf_perm)
    deallocate(buf_sign)

end subroutine

!> calculate a specific entry in the scattering matrix
subroutine averaged_bubble_entry(iq,p,qp,dr,ise,ch,temperature,adaptive_prefactor,bar_bubble,bar_transform,mem)
    !> index of first q-point, irreducible
    integer, intent(in) :: iq
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> interpolated self-energy
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> convolution helper
    type(lo_convolution_handle), intent(in) :: ch
    !> adaptive smearing prefactor
    real(r8), intent(in) :: temperature
    !> smearing prefactor
    real(r8), intent(in) :: adaptive_prefactor
    !> spectrally averaged bubble
    real(r8), dimension(:), intent(out) :: bar_bubble
    real(r8), dimension(:), intent(out) :: bar_transform
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    interpolate: block
        real(r8), parameter :: integraltol=1E-12_r8
        real(r8), dimension(:,:), allocatable :: buf_re,buf_im,buf_j
        real(r8), dimension(:), allocatable :: xlo,xmid,xhi,buf_n
        real(r8) :: f0,f1,f2,scalefactor,pref,sigma
        integer :: imode

        allocate(buf_re(ise%n_energy,dr%n_mode))
        allocate(buf_im(ise%n_energy,dr%n_mode))
        allocate(buf_j(ise%n_energy,dr%n_mode))
        allocate(buf_n(ise%n_energy))
        allocate(xlo(dr%n_mode))
        allocate(xmid(dr%n_mode))
        allocate(xhi(dr%n_mode))
        buf_re=0.0_r8
        buf_im=0.0_r8
        buf_j=0.0_r8
        buf_n=0.0_r8


        ! Now we can fetch J_{qs}
        call ise%evaluate(p,qp%ip(iq)%r,dr%iq(iq)%omega,dr%iq(iq)%egv,buf_re,buf_im,mem)
        !buf_j=0.0_r8
        do imode=1,dr%n_mode

            if ( dr%iq(iq)%omega(imode) .lt. lo_freqtol ) cycle
            ! Figure out extrema for integration
            call lo_find_spectral_function_max_and_fwhm(&
                dr%iq(iq)%omega(imode), &
                ise%omega, &
                buf_im(:,imode), &
                buf_re(:,imode), &
                xmid(imode),&
                xlo(imode),&
                xhi(imode))
            ! Actual integration
            scalefactor=1.0_r8
            call lo_integrate_spectral_function(ise%omega, &
                                                dr%iq(iq)%omega(imode), &
                                                buf_im(:, imode),&
                                                buf_re(:, imode), &
                                                xmid(imode), xlo(imode), xhi(imode), &
                                                scalefactor, temperature, integraltol, f0, f1, f2)
            bar_bubble(imode)=f1*lo_pi !*0.5_r8
            bar_transform(imode)=f2*lo_pi !*0.5_r8

            ! Get something for the similarity transform.
            call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),dr%iq(iq)%omega(imode),buf_j(:,imode))
            sigma=qp%adaptive_sigma( dr%iq(iq)%vel(:,imode),dr%default_smearing(imode), adaptive_prefactor )
            call lo_gaussian_smear_spectral_function(ise%omega,sigma,buf_j(:,imode))
            f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
            buf_j(:,imode)=buf_j(:,imode)/f0

            buf_n = ch%thermal_n(0:ch%n)
            bar_transform = lo_trapezoid_integration(ise%omega,buf_n*(buf_n+1.0_r8)*buf_j(:,imode))

            bar_bubble = lo_trapezoid_integration(ise%omega,buf_j(:,imode)*buf_j(:,imode))
        enddo

        deallocate(buf_re)
        deallocate(buf_im)
        deallocate(buf_j)
        deallocate(xlo)
        deallocate(xmid)
        deallocate(xhi)
    end block interpolate

end subroutine

!> figure out degeneracy fixer thingy
subroutine spectrum_degeneracy_fold_in_fold_out(om,buf,tol)
    !> frequencies
    real(r8), dimension(:), intent(in) :: om
    !> buffer to fix
    real(r8), dimension(:,:), intent(inout) :: buf
    !> tolerance
    real(r8), intent(in) :: tol

    real(r8), dimension(:,:), allocatable :: buf0
    integer :: nb,i,j
    integer :: ctr

    nb=size(om)
    allocate(buf0(size(buf,1),nb))
    buf0=0.0_r8
    do i=1,nb
        ctr=0
        do j=1,nb
            if ( abs(om(i)-om(j)) .lt. tol ) then
                ctr=ctr+1
                buf0(:,i)=buf0(:,i)+buf(:,j)
            endif
        enddo
        buf0(:,i)=buf0(:,i)/real(ctr,r8)
    enddo
    buf=buf0
    deallocate(buf0)
end subroutine

!> figure out degeneracy fixer thingy
subroutine cubic_degeneracy_fold_in_fold_out(om1,om2,om3,buf,tol)
    !> frequencies
    real(r8), dimension(:), intent(in) :: om1,om2,om3
    !> buffer to fix
    real(r8), dimension(:,:,:), intent(inout) :: buf
    !> tolerance
    real(r8), intent(in) :: tol

    real(r8), dimension(:,:,:), allocatable :: buf0
    integer, dimension(:,:,:), allocatable :: dj
    integer, dimension(:,:), allocatable :: di
    integer :: nb,i,j,k,ii,jj,kk
    integer :: ctr1,ctr2,ctr3

    nb=size(om1)

    allocate(di(nb,3))
    di=-1
    ctr1=1
    ctr2=1
    ctr3=1
    di(1,:)=1
    do i=2,nb
        if ( abs(om1(i)-om1(i-1)) .gt. tol ) ctr1=ctr1+1
        if ( abs(om2(i)-om2(i-1)) .gt. tol ) ctr2=ctr2+1
        if ( abs(om3(i)-om3(i-1)) .gt. tol ) ctr3=ctr3+1
        di(i,1)=ctr1
        di(i,2)=ctr2
        di(i,3)=ctr3
    enddo

    allocate(buf0(nb,nb,nb))
    allocate(dj(nb,nb,nb))
    ! Fold in
    buf0=0.0_r8
    dj=0
    do i=1,nb
    do j=1,nb
    do k=1,nb
        ii=di(i,1)
        jj=di(j,2)
        kk=di(k,3)
        buf0(ii,jj,kk)=buf0(ii,jj,kk)+buf(i,j,k)
        dj(ii,jj,kk)=dj(ii,jj,kk)+1
    enddo
    enddo
    enddo
    ! Fold out
    do i=1,nb
    do j=1,nb
    do k=1,nb
        ii=di(i,1)
        jj=di(j,2)
        kk=di(k,3)
        buf(i,j,k)=buf0(ii,jj,kk)/real(dj(ii,jj,kk),r8)
    enddo
    enddo
    enddo
    deallocate(di)
    deallocate(dj)
    deallocate(buf0)
end subroutine

end module
