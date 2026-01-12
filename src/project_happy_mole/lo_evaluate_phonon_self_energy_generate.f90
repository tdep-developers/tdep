submodule(lo_evaluate_phonon_self_energy) lo_evaluate_phonon_self_energy_generate
implicit none
contains

!> Generate the self-energy at a single q-point
module subroutine generate(se, qpoint, qdir, uc, fc, fct, fcf, ise, qp, ddr,&
    temperature, max_energy, n_energy, integrationtype, sigma,&
    use_isotope, use_thirdorder, use_fourthorder, &
    tmr, mw, mem, verbosity)
    !> self energy
    class(lo_phonon_selfenergy), intent(out) :: se
    !> qpoint of interest
    class(lo_qpoint), intent(in) :: qpoint
    !> direction of probe?
    real(r8), dimension(3), intent(in) :: qdir
    !> harmonic properties at this q-point
    !type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> pre-calculated self-energies
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> distributed phonon dispersions
    type(lo_distributed_phonon_dispersions), intent(in) :: ddr
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
    !> timer (start and stop outside this routine)
    type(lo_timer), intent(inout) :: tmr
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot
    integer, intent(in) :: verbosity

    ! harmonic properties at Gamma
    type(lo_distributed_phonon_dispersions_qpoint) :: ompoint

    call mem%tick()
    call tmr%tick()
    ! First thing to do is to initialize all arrays, and set all options
    setopts: block
        integer :: i

        se%integrationtype = integrationtype  ! How to integrate

        ! Number of points on the energy axis
        if (se%integrationtype .eq. 5) then
            se%n_energy = ise%n_energy
        else
            se%n_energy = n_energy
        end if
        se%n_mode = uc%na*3                          ! Number of modes
        se%smearing_prefactor = sigma                ! Adjustment of default smearing parameter
        se%isotope_scattering = use_isotope          ! include isotope scattering?
        se%thirdorder_scattering = use_thirdorder    ! include threephonon term?
        se%fourthorder_scattering = use_fourthorder  ! include fourphonon term?
        se%qdir = qdir                               ! Make note on the probe direction
        ! Make space for all the necessary arrays
        allocate(se%energy_axis(se%n_energy))
        allocate(se%im(se%n_energy, se%n_mode))
        allocate(se%re(se%n_energy, se%n_mode))
        ! Get the energy range
        if (se%integrationtype .eq. 5) then
            se%energy_axis = ise%omega
        else
            se%energy_axis = 0.0_r8
            call lo_linspace(0.0_r8, max_energy, se%energy_axis)
        end if
        ! Set self-energies to nothing, for now.
        se%im = 0.0_r8
        se%re = 0.0_r8

        ! Harmonic properties the new fancy way
        call ompoint%generate(fc,uc,mem,qpoint,qdirection=qdir)

        ! And some space for auxiliary information
        allocate (se%xmid(se%n_mode))
        allocate (se%xlo(se%n_mode))
        allocate (se%xhi(se%n_mode))
        allocate (se%normalization_residual(se%n_mode))
        se%xmid = 0.0_r8
        se%xlo = 0.0_r8
        se%xhi = 0.0_r8
        se%normalization_residual = 0.0_r8

        ! Now things should be clean and nice. Maybe say what we are about to do?
        if (verbosity .gt. 1) then
            write (*, *) '         liso:', se%isotope_scattering
            write (*, *) '  lthirdorder:', se%thirdorder_scattering
            write (*, *) ' lfourthorder:', se%fourthorder_scattering
            write (*, *) '        sigma:', se%smearing_prefactor*lo_mean(ddr%default_smearing)*lo_frequency_hartree_to_thz, ' Thz'
            write (*, *) '  temperature: ', tochar(temperature), ' K'
            write (*, *) '  frequencies:'
            do i = 1, ddr%n_mode
                write (*, *) '    mode ', tochar(i, -3), ', omega: ', tochar(ompoint%omega(i)*lo_frequency_Hartree_to_THz, 6), ' THz'
            end do
        end if
    end block setopts
    call tmr%tock('self-energy initialization')

    ! Cubic + isotope anharmonicity
    call threephonon_imaginary_selfenergy(se,qpoint,ompoint,qp,ddr,ise,uc,fc,fct,temperature,mw,mem,verbosity)
    call tmr%tock('three-phonon integrals')

    ! Maybe fourthorder things
    if (se%fourthorder_scattering) then
        ! fourthorder: block
        !     real(r8), dimension(:), allocatable :: delta
        !     integer :: j

        !     allocate (delta(dr%n_mode))
        !     call fourphonon_real_selfenergy(qpoint, wp, gp, qp, uc, temperature, dr, fcf, delta, se%skipsym, mw, mem, verbosity)
        !     do j = 1, se%n_mode
        !         se%re(:, j) = se%re(:, j) + delta(j)
        !     end do
        !     deallocate (delta)
        ! end block fourthorder
        ! call tmr%tock('four-phonon self-energy')
    end if

    ! ! Fix symmetry, if broken?
    ! fixsymmetry: block
    !     complex(r8), dimension(:,:,:), allocatable :: symmetryop
    !     complex(r8), dimension(:,:), allocatable :: sigma_mode,sigma_xyz,cm0,cm1,cm2
    !     integer :: ie,i,iop,imode


    !     allocate(symmetryop(uc%na*3,uc%na*3,qpoint%n_invariant_operation))
    !     allocate(sigma_mode(uc%na*3,uc%na*3))
    !     allocate(sigma_xyz(uc%na*3,uc%na*3))
    !     allocate(cm0(uc%na*3,uc%na*3))
    !     allocate(cm1(uc%na*3,uc%na*3))
    !     allocate(cm2(uc%na*3,uc%na*3))
    !     symmetryop=0.0_r8
    !     sigma_mode=0.0_r8
    !     sigma_xyz=0.0_r8
    !     cm0=0.0_r8
    !     cm1=0.0_r8
    !     cm2=0.0_r8

    !     do i=1,qpoint%n_invariant_operation
    !         iop = qpoint%invariant_operation(i)
    !         call lo_eigenvector_transformation_matrix(symmetryop(:,:,i),uc%rcart,qpoint%r,uc%sym%op(iop))
    !     enddo

    !     do ie=1,se%n_energy
    !         sigma_mode=0.0_r8
    !         do imode=1,se%n_mode
    !             if ( ompoint%omega(imode) .gt. lo_freqtol ) then
    !                 sigma_mode(imode,imode)=lo_imag*se%im(ie,imode) + se%re(ie,imode)
    !                 sigma_mode(imode,imode)=sigma_mode(imode,imode)*ompoint%omega(imode)
    !             endif
    !         enddo
    !         call lo_gemm(ompoint%egv,sigma_mode,cm0)
    !         call lo_gemm(cm0,ompoint%egv,sigma_xyz,transb='C')

    !         cm0=0.0_r8
    !         do i=1,qpoint%n_invariant_operation
    !             iop = qpoint%invariant_operation(i)
    !             call lo_gemm(symmetryop(:,:,i),sigma_xyz,cm1)
    !             call lo_gemm(cm1,symmetryop(:,:,i),cm2,transb='C')
    !             cm0=cm0+cm2
    !         enddo
    !         cm0=cm0/real(qpoint%n_invariant_operation,r8)

    !         ! then transform back again?
    !         call lo_gemm(ompoint%egv,cm0,cm1,transa='C')
    !         call lo_gemm(cm1,ompoint%egv,cm2)
    !         do imode=1,se%n_mode
    !             if ( ompoint%omega(imode) .gt. lo_freqtol ) then
    !                 se%im(ie,imode)=aimag(cm2(imode,imode))/ompoint%omega(imode)
    !                 se%re(ie,imode)=real(cm2(imode,imode),r8)/ompoint%omega(imode)
    !             endif
    !         enddo
    !     enddo
    ! end block fixsymmetry

    ! Massage the self-energy ever so slightly
    massageselfenergy: block
        real(r8), parameter :: KK_prefactor = 2.0_r8/lo_pi
        complex(r8), dimension(:), allocatable :: z0
        complex(r8) :: eta
        real(r8), dimension(:), allocatable :: taper,Omega,Omegasquare,y0
        real(r8) :: dOmega,Omegaprime
        integer :: imode,ie,ctr

        ! First step is to ensure the imaginary part of the self-energy goes to zero
        ! at large Omega, and that it is strictly positive.
        call mem%allocate(taper, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        taper=0.0_r8
        call lo_tapering_function(se%energy_axis,taper)
        do imode=1,se%n_mode
            se%im(:,imode)=max(se%im(:,imode),0.0_r8)
            se%im(:,imode)=se%im(:,imode)*taper
        enddo
        call mem%deallocate(taper, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Now we KK-transform to get the raw real part of the self-energy:
        call mem%allocate(Omega, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(Omegasquare, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(z0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(y0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        Omega = se%energy_axis
        Omegasquare = Omega**2
        dOmega = Omega(2) - Omega(1)                  ! prefactor for riemann integral
        eta = lo_imag*(Omega(2) - Omega(1))*1E-8_r8   ! small imaginary thing to avoid divergence

        ctr = 0
        se%re = 0.0_r8
        do imode = 1, se%n_mode
        do ie = 1, se%n_energy
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            Omegaprime = Omega(ie)
            y0 = se%im(:, imode)*Omega
            z0 = (Omegaprime + eta)**2 - Omegasquare
            y0 = real(y0/z0, r8)
            y0(1) = y0(1)*0.5_r8
            y0(se%n_energy) = y0(se%n_energy)*0.5_r8
            se%re(ie, imode) = sum(y0)*dOmega*KK_prefactor
        end do
        end do
        call mw%allreduce('sum', se%re)

        ! Intermediate cleanup
        call mem%deallocate(Omega, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(Omegasquare, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(z0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(y0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Now for the tricky part, we need to get the spectral function to normalize properly.
        ! To that end I start by forcing the self-energy at zero to zero. This is a rigid shift.
        do imode=1,se%n_mode
            se%re(:,imode)=se%re(:,imode) - se%re(1,imode)
        enddo

        ! Then we add a quadratic term to fix normalization
        call add_quadratic_to_fix_normalization(ompoint%omega,se%re,se%im,se%energy_axis,se%normalization_residual,mw)

        call tmr%tock('normalize spectral function')
    end block massageselfenergy

    ! Make sure I did not mess up with memory management.
    call mem%tock(__FILE__, __LINE__, mw%comm)
end subroutine

end submodule