submodule(lo_selfenergy_interpolation) lo_selfenergy_interpolation_evaluate
implicit none

contains

!> evaluate phonon self-energy at arbitrary point
module subroutine evaluate_self_energy(ise,p,qv,omega,egv,sigma_Re,sigma_Im,mem,mw)
    !> self-energy interpolation
    class(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point
    real(r8), dimension(3), intent(in) :: qv
    !> harmonic phonon frequencies at this q
    real(r8), dimension(:), intent(in) :: omega
    !> harmonic phonon eigenvectors at this q
    complex(r8), dimension(:,:), intent(in) :: egv
    !> real part of self-energy
    real(r8), dimension(:,:), intent(out) :: sigma_Re
    !> imaginary part of self-energy
    real(r8), dimension(:,:), intent(out) :: sigma_Im
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> perhaps we try to make it parallel?
    type(lo_mpi_helper), intent(inout), optional :: mw

    type(lo_qpoint) :: qpoint
    complex(r8), dimension(:,:), allocatable :: left_trf,right_trf
    real(r8), dimension(:,:,:), allocatable :: buf_re,buf_im
    integer :: n_mode

    ! Some prep
    init: block
        n_mode=size(omega)

        qpoint%r=qv
        call lo_get_small_group_of_qpoint(qpoint, p)

        allocate(buf_re(n_mode,n_mode,ise%n_energy))
        allocate(buf_im(n_mode,n_mode,ise%n_energy))
        buf_re=0.0_r8
        buf_im=0.0_r8
    end block init

    ! Linearly interpolate the raw self-energy to this q
    linear_interpolation: block
        complex(r8), dimension(:,:,:), allocatable :: rotmat
        complex(r8), dimension(:,:), allocatable :: cm0,cm1
        real(r8), dimension(4) :: weight
        integer, dimension(4) :: irr_ind,full_ind !,operation
        integer :: i,iq,ie,ii,iop

        call ise%box%indices_and_weights(ise%qp,p,qv,weight,irr_ind,full_ind)

        ! Have to symmetry-rotate self-energy
        allocate(rotmat(n_mode,n_mode,4))
        rotmat=0.0_r8
        do i=1,4
            ii=full_ind(i)
            iop=ise%qp%ap(ii)%operation_from_irreducible
            if ( iop .gt. 0 ) then
                call lo_eigenvector_transformation_matrix(rotmat(:,:,i),p%rcart,ise%qp%ip( irr_ind(i) )%r,p%sym%op(iop),inverseoperation=.false.)
            else
                call lo_eigenvector_transformation_matrix(rotmat(:,:,i),p%rcart,ise%qp%ip( irr_ind(i) )%r,p%sym%op(-iop),inverseoperation=.true.)
            endif
        enddo

        allocate(cm0(n_mode,n_mode))
        allocate(cm1(n_mode,n_mode))
        cm0=0.0_r8
        cm1=0.0_r8
        do i=1,4
            iq=irr_ind(i)
            do ie=1,ise%n_energy
                cm0=ise%sigma_Re(:,:,ie,iq) + lo_imag*ise%sigma_Im(:,:,ie,iq)
                call lo_gemm(rotmat(:,:,i),cm0,cm1)
                call lo_gemm(cm1,rotmat(:,:,i),cm0,transb='C')
                buf_re(:,:,ie)=buf_re(:,:,ie) + weight(i)*real(cm0,r8)
                buf_im(:,:,ie)=buf_im(:,:,ie) + weight(i)*aimag(cm0)
            enddo
        enddo
        deallocate(cm0)
        deallocate(cm1)
        deallocate(rotmat)
    end block linear_interpolation

    build_trafo: block
        type(lo_phonon_dispersions_qpoint) :: aux_wp
        complex(r8), dimension(:,:), allocatable :: aux_eig,aux_inveig,aux_trafo,eig,inveig
        real(r8) :: f0,f1
        integer :: imode,iatom,ialpha,ii

        allocate(left_trf(n_mode,n_mode))
        allocate(right_trf(n_mode,n_mode))

        allocate(eig(n_mode,n_mode))
        allocate(inveig(n_mode,n_mode))
        allocate(aux_eig(n_mode,n_mode))
        allocate(aux_inveig(n_mode,n_mode))
        allocate(aux_trafo(n_mode,n_mode))

        left_trf=0.0_r8
        right_trf=0.0_r8
        eig=0.0_r8
        inveig=0.0_r8
        aux_eig=0.0_r8
        aux_inveig=0.0_r8
        aux_trafo=0.0_r8

        ! Transformation with real dispersions
        if ( ise%polar ) then
            eig=egv
            inveig=transpose(conjg(eig))
            do imode=1,n_mode
                if ( omega(imode) .gt. lo_freqtol ) then
                    f0=sqrt(omega(imode))
                    f1=1.0_r8/f0
                else
                    f0=0.0_r8
                    f1=0.0_r8
                endif

                do iatom=1,p%na
                    do ialpha=1,3
                        ii=(iatom-1)*3 + ialpha
                        eig(ii,imode)=eig(ii,imode)*f1
                        inveig(imode,ii)=inveig(imode,ii)*f1
                    enddo
                enddo
            enddo

            ! Transformation with auxiliary dispersions
            call aux_wp%generate(ise%aux_fc,p,mem,qpoint=qpoint)
            aux_eig=aux_wp%egv
            aux_inveig=transpose(conjg(aux_eig))
            do imode=1,n_mode
                if ( aux_wp%omega(imode) .gt. lo_freqtol ) then
                    f0=sqrt(aux_wp%omega(imode))
                    f1=1.0_r8/f0
                else
                    f0=0.0_r8
                    f1=0.0_r8
                endif

                do iatom=1,p%na
                    do ialpha=1,3
                        ii=(iatom-1)*3 + ialpha
                        aux_eig(ii,imode)=aux_eig(ii,imode)*f0
                        aux_inveig(imode,ii)=aux_inveig(imode,ii)*f0
                    enddo
                enddo
            enddo
            call lo_gemm(aux_eig,aux_inveig,aux_trafo)

            ! Build the actual transformations
            call lo_gemm(inveig,aux_trafo,left_trf,transb='C')
            call lo_gemm(aux_trafo,eig,right_trf)
        else
            ! Non-polar interpolation, much easier
            eig=egv
            inveig=transpose(conjg(eig))
            do imode=1,n_mode
                if ( omega(imode) .gt. lo_freqtol ) then
                    f0=sqrt(omega(imode))
                    f1=1.0_r8/f0
                else
                    f0=0.0_r8
                    f1=0.0_r8
                endif

                do iatom=1,p%na
                    do ialpha=1,3
                        ii=(iatom-1)*3 + ialpha
                        eig(ii,imode)=eig(ii,imode)*f0
                        inveig(imode,ii)=inveig(imode,ii)*f0
                    enddo
                enddo
            enddo
            left_trf=inveig
            right_trf=eig
        endif

        deallocate(eig)
        deallocate(inveig)
        deallocate(aux_eig)
        deallocate(aux_inveig)
        deallocate(aux_trafo)
    end block build_trafo

    convert_to_mode: block
        complex(r8), dimension(:,:), allocatable :: cm0,cm1
        integer :: imode,ie

        ! Actuall allocate
        allocate(cm0(n_mode,n_mode))
        allocate(cm1(n_mode,n_mode))
        cm0=0.0_r8
        cm1=0.0_r8

        if ( present(mw) ) then
            sigma_re=0.0_r8
            sigma_im=0.0_r8
            do ie=1,ise%n_energy
                if ( mod(ie,mw%n) .ne. mw%r ) cycle
                cm1 = buf_re(:,:,ie) + lo_imag*buf_im(:,:,ie)
                call lo_gemm(left_trf,cm1,cm0)
                call lo_gemm(cm0,right_trf,cm1)
                do imode=1,n_mode
                    sigma_re(ie,imode)=real( cm1(imode,imode),r8)
                    sigma_im(ie,imode)=aimag( cm1(imode,imode) )
                enddo
            enddo
            call mw%allreduce('sum',sigma_Im)
            call mw%allreduce('sum',sigma_Re)
        else
            sigma_re=0.0_r8
            sigma_im=0.0_r8
            do ie=1,ise%n_energy
                ! Self-energy in weighted xyz coordinates
                cm1 = buf_re(:,:,ie) + lo_imag*buf_im(:,:,ie)

                call lo_gemm(left_trf,cm1,cm0)
                call lo_gemm(cm0,right_trf,cm1)
                do imode=1,n_mode
                    sigma_re(ie,imode)=real( cm1(imode,imode),r8)
                    sigma_im(ie,imode)=aimag( cm1(imode,imode) )
                enddo
            enddo
        endif

        deallocate(cm0)
        deallocate(cm1)
    end block convert_to_mode

    ! ! Enforce degeneracy?
    ! if ( maxval(wp%degeneracy) .gt. 1 ) then
    !     allocate(modefixed(p%na*3))
    !     allocate(bf0(ise%n_energy))
    !     allocate(bf1(ise%n_energy))
    !     bf0=0.0_r8
    !     bf1=0.0_r8
    !     modefixed=.false.
    !     do imode=1,p%na*3
    !         if ( modefixed(imode) ) cycle
    !         bf0=0.0_r8
    !         bf1=0.0_r8
    !         do i=1,wp%degeneracy(imode)
    !             jmode=wp%degenmode(i,imode)
    !             bf0=bf0 + sigma_Re(:,jmode)
    !             bf1=bf1 + sigma_Im(:,jmode)
    !         enddo
    !         bf0=bf0/real(wp%degeneracy(imode),r8)
    !         bf1=bf1/real(wp%degeneracy(imode),r8)
    !         do i=1,wp%degeneracy(imode)
    !             jmode=wp%degenmode(i,imode)
    !             sigma_Re(:,jmode)=bf0
    !             sigma_Im(:,jmode)=bf1
    !             modefixed(jmode)=.true.
    !         enddo
    !     enddo
    !     deallocate(modefixed)
    !     deallocate(bf0)
    !     deallocate(bf1)
    ! endif

end subroutine

!module subroutine evaluate_smeared_J(ise,p,qp,wp)
!end subroutine

end submodule
