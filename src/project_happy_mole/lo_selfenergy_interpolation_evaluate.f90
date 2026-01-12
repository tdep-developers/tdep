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
    complex(r8), dimension(:,:,:), allocatable :: buf_cre,buf_cim
    integer :: n_mode

    ! Some prep
    init: block
        n_mode=size(omega)

        qpoint%r=qv

        allocate(buf_cre(n_mode,n_mode,ise%n_energy))
        allocate(buf_cim(n_mode,n_mode,ise%n_energy))
        buf_cre=0.0_r8
        buf_cim=0.0_r8
    end block init

    ! ! Linearly interpolate the raw self-energy to this q
    ! linear_interpolation: block
    !     complex(r8), dimension(:,:,:), allocatable :: rotmat
    !     complex(r8), dimension(:,:), allocatable :: cm0,cm1
    !     real(r8), dimension(4) :: weight
    !     integer, dimension(4) :: irr_ind,full_ind
    !     integer :: i,iq,ie,ii,iop

    !     call ise%box%indices_and_weights(ise%qp,p,qv,weight,irr_ind,full_ind)

    !     ! Have to symmetry-rotate self-energy
    !     allocate(rotmat(n_mode,n_mode,4))
    !     rotmat=0.0_r8
    !     do i=1,4
    !         ii=full_ind(i)
    !         iop=ise%qp%ap(ii)%operation_from_irreducible
    !         if ( iop .gt. 0 ) then
    !             call lo_eigenvector_transformation_matrix(rotmat(:,:,i),p%rcart,ise%qp%ip( irr_ind(i) )%r,p%sym%op(iop),inverseoperation=.false.)
    !         else
    !             call lo_eigenvector_transformation_matrix(rotmat(:,:,i),p%rcart,ise%qp%ip( irr_ind(i) )%r,p%sym%op(-iop),inverseoperation=.true.)
    !         endif
    !     enddo

    !     allocate(cm0(n_mode,n_mode))
    !     allocate(cm1(n_mode,n_mode))
    !     cm0=0.0_r8
    !     cm1=0.0_r8

    !     buf_cre=0.0_r8
    !     buf_cim=0.0_r8
    !     do i=1,4
    !         iq=irr_ind(i)
    !         do ie=1,ise%n_energy

    !             cm0=ise%sigma_Re(:,:,ie,iq)
    !             call lo_gemm(rotmat(:,:,i),cm0,cm1)
    !             call lo_gemm(cm1,rotmat(:,:,i),cm0,transb='C')
    !             buf_cre(:,:,ie)=buf_cre(:,:,ie) + weight(i)*cm0

    !             cm0=ise%sigma_Im(:,:,ie,iq)
    !             call lo_gemm(rotmat(:,:,i),cm0,cm1)
    !             call lo_gemm(cm1,rotmat(:,:,i),cm0,transb='C')
    !             buf_cim(:,:,ie)=buf_cim(:,:,ie) + weight(i)*cm0
    !         enddo
    !     enddo
    !     deallocate(cm0)
    !     deallocate(cm1)
    !     deallocate(rotmat)
    ! end block linear_interpolation

    fourier_interpolation: block
        complex(r8), dimension(:), allocatable :: phasefactor
        real(r8) :: kdotr
        integer :: ie,ir,a1,a2

        allocate(phasefactor(ise%n_rvec))
        phasefactor=0.0_r8
        do ir=1,ise%n_rvec
            kdotr=dot_product(qv,ise%rvec(:,ir))*lo_twopi
            phasefactor(ir)=cmplx(cos(kdotr),sin(kdotr),r8)
        enddo

        buf_cim=0.0_r8
        buf_cre=0.0_r8
        do ie=1,ise%n_energy
            do ir=1,ise%n_rvec
                a1=ise%atomind(1,ir)
                a2=ise%atomind(2,ir)

                buf_cim((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,ie)=&
                buf_cim((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,ie)+&
                phasefactor(ir)*(ise%ifc(:,:,ir,ie))

                buf_cre((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,ie)=&
                buf_cre((a2-1)*3+1:a2*3,(a1-1)*3+1:a1*3,ie)+&
                phasefactor(ir)*(ise%rfc(:,:,ir,ie))
            enddo
        enddo
    end block fourier_interpolation

    build_trafo: block
        complex(r8), dimension(:,:), allocatable :: aux_eig,aux_inveig,aux_trafo,eig,inveig
        complex(r8), dimension(:,:), allocatable :: cm0,cm1
        real(r8), dimension(:), allocatable :: aux_omega
        real(r8) :: f0,f1,fm0,fm1
        integer :: imode,iatom,ialpha,ii

        allocate(left_trf(n_mode,n_mode))
        allocate(right_trf(n_mode,n_mode))

        allocate(eig(n_mode,n_mode))
        allocate(inveig(n_mode,n_mode))
        allocate(aux_eig(n_mode,n_mode))
        allocate(aux_inveig(n_mode,n_mode))
        allocate(aux_trafo(n_mode,n_mode))
        allocate(cm0(n_mode,n_mode))
        allocate(cm1(n_mode,n_mode))
        allocate(aux_omega(n_mode))

        left_trf=0.0_r8
        right_trf=0.0_r8
        eig=0.0_r8
        inveig=0.0_r8
        aux_eig=0.0_r8
        aux_inveig=0.0_r8
        aux_trafo=0.0_r8

        ! Transformation with real dispersions
        if ( ise%polar ) then

            call ise%aux_fc%dynamicalmatrix(p,qpoint,cm0,mem)
            call lo_gemm(cm0,egv,cm1)
            call lo_gemm(egv,cm1,cm0,transa='C')
            do imode=1,n_mode
                aux_omega(imode)=sqrt( abs(cm0(imode,imode)) )
            enddo

            eig=egv
            inveig=transpose(conjg(eig))
            do imode=1,n_mode
                if ( omega(imode) .gt. lo_freqtol ) then
                    f0=sqrt(omega(imode))
                    f0=f0/aux_omega(imode)
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
            do imode=1,n_mode
                do iatom=1,p%na
                    fm0=p%invsqrtmass(iatom)
                    fm1=1.0_r8/fm0
                    do ialpha=1,3
                        ii=(iatom-1)*3 + ialpha
                        eig(ii,imode)=eig(ii,imode)*fm0
                        inveig(imode,ii)=inveig(imode,ii)*fm0
                    enddo
                enddo
            enddo
            right_trf=eig
            left_trf=inveig
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
            write(*,*) 'FIXME PARALLEL',__FILE__,__LINE__
            stop
            ! sigma_re=0.0_r8
            ! sigma_im=0.0_r8
            ! do ie=1,ise%n_energy
            !     if ( mod(ie,mw%n) .ne. mw%r ) cycle
            !     cm1 = buf_re(:,:,ie) + lo_imag*buf_im(:,:,ie)
            !     call lo_gemm(left_trf,cm1,cm0)
            !     call lo_gemm(cm0,right_trf,cm1)
            !     do imode=1,n_mode
            !         sigma_re(ie,imode)=real( cm1(imode,imode),r8)
            !         sigma_im(ie,imode)=aimag( cm1(imode,imode) )
            !     enddo
            ! enddo
            ! call mw%allreduce('sum',sigma_Im)
            ! call mw%allreduce('sum',sigma_Re)
        else
            sigma_re=0.0_r8
            sigma_im=0.0_r8
            do ie=1,ise%n_energy
                ! Self-energy in weighted xyz coordinates
                cm1=buf_cim(:,:,ie) + buf_cre(:,:,ie)
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

    ! Enforce degeneracy?
    fixdegen: block
        real(r8), dimension(:), allocatable :: bf0,bf1
        logical, dimension(:), allocatable :: modefixed
        integer :: imode,jmode,k

        allocate(modefixed(p%na*3))
        allocate(bf0(ise%n_energy))
        allocate(bf1(ise%n_energy))
        bf0=0.0_r8
        bf1=0.0_r8
        modefixed=.false.
        do imode=1,p%na*3
            if ( modefixed(imode) ) cycle
            bf0=0.0_r8
            bf1=0.0_r8
            k=0
            do jmode=1,p%na*3
                if ( abs(omega(imode)-omega(jmode)) .lt. lo_freqtol ) then
                    k=k+1
                    bf0=bf0 + sigma_Re(:,jmode)
                    bf1=bf1 + sigma_Im(:,jmode)
                endif
            enddo
            bf0=bf0/real(k,r8)
            bf1=bf1/real(k,r8)
            do jmode=1,p%na*3
                if ( abs(omega(imode)-omega(jmode)) .lt. lo_freqtol ) then
                    modefixed(jmode)=.true.
                    sigma_Re(:,jmode)=bf0
                    sigma_Im(:,jmode)=bf1
                endif
            enddo
        enddo
        deallocate(modefixed)
        deallocate(bf0)
        deallocate(bf1)
    end block fixdegen

end subroutine


!module subroutine evaluate_smeared_J(ise,p,qp,wp)
!end subroutine

end submodule
