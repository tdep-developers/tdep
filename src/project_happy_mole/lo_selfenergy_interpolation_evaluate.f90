submodule(lo_selfenergy_interpolation) lo_selfenergy_interpolation_evaluate
implicit none

contains

!> evaluate phonon self-energy at arbitrary point
module subroutine evaluate_self_energy(ise,p,qv,wp,sigma_Re,sigma_Im,mw)
    !> self-energy interpolation
    class(lo_interpolated_selfenergy_grid), intent(in) :: ise
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point
    real(r8), dimension(3), intent(in) :: qv
    !> harmonic phonon properties at this q
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> real part of self-energy
    real(r8), dimension(:,:), intent(out) :: sigma_Re
    !> imaginary part of self-energy
    real(r8), dimension(:,:), intent(out) :: sigma_Im
    !> perhaps we try to make it parallel?
    type(lo_mpi_helper), intent(inout), optional :: mw

    complex(r8), dimension(:,:), allocatable :: cm0,cm1,cm2
    real(r8), dimension(:,:,:), allocatable :: buf_re,buf_im
    real(r8), dimension(:), allocatable :: bf0,bf1
    real(r8), dimension(4) :: weight
    integer, dimension(4) :: ind
    integer :: n_mode,i,iq,ie,imode,jmode
    logical, dimension(:), allocatable :: modefixed

    !write(*,*) ''
    !write(*,*) 'EVALUATE',matmul(p%inv_reciprocal_latticevectors,qv)

    call ise%box%indices_and_weights(ise%qp,p,qv,weight,ind)
    n_mode=size(wp%omega)

    ! Actuall allocate
    allocate(buf_re(n_mode,n_mode,ise%n_energy))
    allocate(buf_im(n_mode,n_mode,ise%n_energy))
    allocate(cm0(n_mode,n_mode))
    allocate(cm1(n_mode,n_mode))
    allocate(cm2(n_mode,n_mode))
    buf_re=0.0_r8
    buf_im=0.0_r8
    cm0=0.0_r8
    cm1=0.0_r8
    cm2=0.0_r8

    do i=1,4
        iq=ind(i)
        buf_re(:,:,:)=buf_re(:,:,:) + weight(i)*ise%sigma_Re(:,:,:,iq)
        buf_im(:,:,:)=buf_im(:,:,:) + weight(i)*ise%sigma_Im(:,:,:,iq)
    enddo

    cm0=0.0_r8
    do imode=1,n_mode
        if ( wp%omega(imode) .gt. lo_freqtol ) then
            cm0(:,imode)=wp%egv(:,imode) !/sqrt(wp%omega(imode))
        endif
    enddo

    if ( present(mw) ) then
        sigma_re=0.0_r8
        sigma_im=0.0_r8
        do ie=1,ise%n_energy
            if ( mod(ie,mw%n) .ne. mw%r ) cycle
            cm1 = buf_re(:,:,ie)
            cm2 = matmul( conjg(transpose(cm0)),cm1 )
            cm1 = matmul(cm2,cm0)
            do imode=1,n_mode
                sigma_re(ie,imode)=real( cm1(imode,imode),r8)
            enddo

            cm1 = buf_im(:,:,ie)
            cm2 = matmul( conjg(transpose(cm0)),cm1 )
            cm1 = matmul(cm2,cm0)
            do imode=1,n_mode
                sigma_im(ie,imode)=real( cm1(imode,imode),r8)
            enddo
        enddo
        call mw%allreduce('sum',sigma_Im)
        call mw%allreduce('sum',sigma_Re)
    else
        do ie=1,ise%n_energy
            cm1 = buf_re(:,:,ie)
            cm2 = matmul( conjg(transpose(cm0)),cm1 )
            cm1 = matmul(cm2,cm0)
            do imode=1,n_mode
                sigma_re(ie,imode)=real( cm1(imode,imode),r8)
            enddo

            cm1 = buf_im(:,:,ie)
            cm2 = matmul( conjg(transpose(cm0)),cm1 )
            cm1 = matmul(cm2,cm0)
            do imode=1,n_mode
                sigma_im(ie,imode)=real( cm1(imode,imode),r8)
            enddo
        enddo
    endif

    deallocate(cm0)
    deallocate(cm1)
    deallocate(cm2)

    ! Enforce degeneracy?
    if ( maxval(wp%degeneracy) .gt. 1 ) then
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
            do i=1,wp%degeneracy(imode)
                jmode=wp%degenmode(i,imode)
                bf0=bf0 + sigma_Re(:,jmode)
                bf1=bf1 + sigma_Im(:,jmode)
            enddo
            bf0=bf0/real(wp%degeneracy(imode),r8)
            bf1=bf1/real(wp%degeneracy(imode),r8)
            do i=1,wp%degeneracy(imode)
                jmode=wp%degenmode(i,imode)
                sigma_Re(:,jmode)=bf0
                sigma_Im(:,jmode)=bf1
                modefixed(jmode)=.true.
            enddo
        enddo
        deallocate(modefixed)
        deallocate(bf0)
        deallocate(bf1)
    endif

end subroutine

end submodule
