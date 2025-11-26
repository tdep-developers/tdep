submodule(lo_evaluate_phonon_self_energy) lo_evaluate_phonon_self_energy_threephonon
implicit none



contains

!> evaluate the isotope self-energy
module subroutine threephonon_imaginary_selfenergy(se,qpoint,ompoint,qp,dr,ise,p,fc2,fc3,temperature,mw,mem,verbosity)
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> q-point we are investigating
    class(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties at this q-point
    type(lo_distributed_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> distributed harmonic dispersions
    type(lo_distributed_phonon_dispersions), intent(in) :: dr
    !> self-energy interpolation
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc2
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fc3
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    select case (se%integrationtype)
    case (1:2)
        call threephonon_imaginary_selfenergy_gaussian_v0(qpoint, ompoint, se, qp, dr, p, fc2, fc3, temperature, mw, mem, verbosity)
    case (4)
        call threephonon_imaginary_selfenergy_convolution_v0(qpoint, ompoint, se, qp, dr, p, fc2, fc3, ise, temperature, mw, mem, verbosity)
    case default
        call lo_stop_gracefully(['Unknown integration type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end select
end subroutine

!> threephonon self energy gaussian, with on-the-fly matrix elements
subroutine threephonon_imaginary_selfenergy_gaussian_v0(qpoint, ompoint, se, qp, dr, p, fc2, fc3, temperature, mw, mem, verbosity)
    !> q-point we are investigating
    class(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties at q-point we are investigating
    type(lo_distributed_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_distributed_phonon_dispersions), intent(in) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc2
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fc3
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_distributed_phonon_dispersions_qpoint) :: op3
    complex(r8), dimension(:), allocatable :: ptf_phi,evp1,evp2
    real(r8), dimension(:,:,:), allocatable :: psisq_3ph
    real(r8), dimension(:,:), allocatable :: psisq_isotope
    real(r8), dimension(:), allocatable :: buf
    real(r8) :: t0,invf
    integer :: ipt,iq

    init: block
        t0 = walltime()
        invf = se%n_energy/se%energy_axis(se%n_energy)
        se%im = 0.0_r8

        call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ptf_phi, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(psisq_3ph,[dr%n_mode,dr%n_mode,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(psisq_isotope,[dr%n_mode,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf = 0.0_r8
        ptf_phi = 0.0_r8
        evp1 = 0.0_r8
        evp2 = 0.0_r8
        psisq_3ph = 0.0_r8
        psisq_isotope = 0.0_r8
    end block init

    qploop: do ipt=1,dr%n_full_qpoint_local

        ! Need to know the global q-point index
        iq=dr%aq(ipt)%global_full_index

        matrixelement: block
            complex(r8) :: c0
            real(r8), dimension(3) :: qv1, qv2, qv3
            integer :: b1,b2,b3

            qv1=qpoint%r
            qv2=qp%ap(iq)%r
            qv3=-qv1-qv2

            if ( se%thirdorder_scattering ) then
                call op3%generate(fc2,p,mem,qvec=qv3)
                do b1=1,se%n_mode
                    op3%sigma(b1) = qp%adaptive_sigma( qp%ap(iq)%radius,op3%vel(:,b1),dr%default_smearing(b1),se%smearing_prefactor)
                enddo

                call pretransform_phi(fc3, qv2, qv3, ptf_phi)

                psisq_3ph=0.0_r8
                do b1 = 1, dr%n_mode
                do b2 = 1, dr%n_mode
                    evp1 = 0.0_r8
                    !call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv2(:, b2, q), 1, fh%ugv1(:, b1), 1, fh%evp1, dr%n_mode)
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), dr%aq(ipt)%nuvec(:,b2), 1, ompoint%nuvec(:, b1), 1, evp1, dr%n_mode)
                    do b3 = 1, dr%n_mode
                        evp2 = 0.0_r8
                        !call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv3(:, b3, q), 1, fh%evp1, 1, fh%evp2, dr%n_mode)
                        call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), op3%nuvec(:, b3), 1, evp1, 1, evp2, dr%n_mode)
                        evp2 = conjg(evp2)
                        c0=dot_product(evp2, ptf_phi)
                        psisq_3ph(b1, b2, b3)=abs(conjg(c0)*c0)
                    end do
                enddo
                enddo
            else
                psisq_3ph=0.0_r8
            endif

            if ( se%isotope_scattering ) then
                do b1=1,dr%n_mode
                do b2=1,dr%n_mode
                    psisq_isotope(b1,b2)=isotope_scattering_strength(p,ompoint%egv(:,b1),dr%aq(ipt)%egv(:,b2),ompoint%omega(b1),dr%aq(ipt)%omega(b2))
                enddo
                enddo
            else
                psisq_isotope=0.0_r8
            endif

        end block matrixelement

        frequencyintegral: block
            real(r8) :: s2,s3,sigma,om2,om3,plf1,plf2
            integer :: b1,b2,b3
            integer :: ilo,ihi,ii,jj,i

            psisq_3ph=psisq_3ph*threephonon_prefactor*qp%ap(iq)%integration_weight
            psisq_isotope=psisq_isotope*isotope_prefactor*qp%ap(iq)%integration_weight

            if ( se%thirdorder_scattering ) then
                ! Slightly smarter version that should be a little faster.
                do b2 = 1, dr%n_mode
                do b3 = 1, dr%n_mode
                    ! reset buffer
                    buf = 0.0_r8
                    ilo = lo_hugeint
                    ihi = -lo_hugeint
                    ! Get the smearing parameter
                    s2 = dr%aq(ipt)%sigma(b2)
                    s3 = op3%sigma(b3)
                    sigma = sqrt(s2**2 + s3**2)
                    ! fetch frequencies
                    !om1 = ompoint%omega(b1)
                    om2 = dr%aq(ipt)%omega(b2)
                    om3 = op3%omega(b3)
                    ! Get the S-function
                    plf1 = 1.0_r8 + lo_planck(temperature, om2) + lo_planck(temperature, om3)
                    plf2 = lo_planck(temperature, om3) - lo_planck(temperature, om2)

                    ! delta(bigOM-om2-om3)
                    ii = max(floor((om2 + om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((om2 + om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) + plf1*lo_gauss(se%energy_axis(i), om2 + om3, sigma)
                    end do
                    ! delta(bigOM+om2+om3)
                    ii = max(floor((-om2 - om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((-om2 - om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) - plf1*lo_gauss(se%energy_axis(i), -om2 - om3, sigma)
                    end do
                    ! delta(bigOM-om2+om3)
                    ii = max(floor((om2 - om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((om2 - om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) + plf2*lo_gauss(se%energy_axis(i), om2 - om3, sigma)
                    end do
                    ! delta(bigOM+om2-om3)
                    ii = max(floor((-om2 + om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((-om2 + om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) - plf2*lo_gauss(se%energy_axis(i), -om2 + om3, sigma)
                    end do

                    ! Increment the self-energy
                    if (ilo .lt. ihi) then
                        ! Three-phonon scattering
                        do b1 = 1, dr%n_mode
                            se%im(ilo:ihi, b1) = se%im(ilo:ihi, b1) + buf(ilo:ihi)*psisq_3ph(b1,b2,b3)
                        enddo
                    end if
                end do
                end do
            endif

            if ( se%isotope_scattering ) then
                do b1 = 1,dr%n_mode
                do b2 = 1,dr%n_mode
                    sigma=dr%aq(ipt)%sigma(b2)
                    ii = max(floor((dr%aq(ipt)%omega(b2) - 4*sigma)*invf), 1)
                    jj = min(ceiling((dr%aq(ipt)%omega(b2) + 4*sigma)*invf), se%n_energy)
                    do i = ii, jj
                        se%im(i, b1) = se%im(i, b1) + lo_gauss(se%energy_axis(i), dr%aq(ipt)%omega(b2), sigma)*psisq_isotope(b1,b2)
                    end do
                enddo
                enddo
            endif

        end block frequencyintegral

        if (verbosity .gt. 0) then
            if (lo_trueNtimes(ipt, 127, dr%n_full_qpoint_local)) call lo_progressbar(' ... threephonon imaginary selfenergy', ipt, dr%n_full_qpoint_local)
        end if
    enddo qploop

    ! Sum it up
    call mw%allreduce('sum', se%im)

    ! Fix degeneracies? Why not.
    fixdegen: block
        integer :: b1,b2,i
        do b1 = 1, dr%n_mode
            buf = 0.0_r8
            do i = 1, ompoint%degeneracy(b1)
                b2 = ompoint%degenmode(i, b1)
                buf = buf + se%im(:, b2)
            end do
            buf = buf/real(ompoint%degeneracy(b1), r8)
            do i = 1, ompoint%degeneracy(b1)
                b2 = ompoint%degenmode(i, b1)
                se%im(:, b2) = buf
            end do
        end do
    end block fixdegen

    ! And some cleanup
    call mem%deallocate(buf,     persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(ptf_phi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1,    persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2,    persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(psisq_3ph,     persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(psisq_isotope, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (verbosity .gt. 0) call lo_progressbar(' ... threephonon imaginary selfenergy', dr%n_mode*qp%n_full_point, dr%n_mode*qp%n_full_point, walltime() - t0)
end subroutine

!> threephonon self energy, convolution based, with on-the-fly matrix elements
subroutine threephonon_imaginary_selfenergy_convolution_v0(qpoint, ompoint, se, qp, dr, p, fc2, fc3, ise, temperature, mw, mem, verbosity)
    !> q-point we are investigating
    class(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties at q-point we are investigating
    type(lo_distributed_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_distributed_phonon_dispersions), intent(in) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc2
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fc3
    !> interpolated self-energy
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_convolution_handle) :: ch
    type(lo_distributed_phonon_dispersions_qpoint) :: op3
    complex(r8), dimension(:,:), allocatable :: cbuf1
    complex(r8), dimension(:), allocatable :: cbuf0
    complex(r8), dimension(:), allocatable :: ptf_phi,evp1,evp2
    real(r8), dimension(:,:,:), allocatable :: psisq_3ph
    real(r8), dimension(:,:), allocatable :: psisq_isotope
    real(r8) :: t0
    integer :: ipt,iq

    init: block
        t0 = walltime()
        se%im = 0.0_r8

        call mem%allocate(ptf_phi, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(psisq_3ph,[dr%n_mode,dr%n_mode,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(psisq_isotope,[dr%n_mode,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ptf_phi = 0.0_r8
        evp1 = 0.0_r8
        evp2 = 0.0_r8
        psisq_3ph = 0.0_r8
        psisq_isotope = 0.0_r8

        call ch%generate(ise%omega,temperature,dr%n_mode)
        ! Need fft-compatible buffers to accumulate
        allocate(cbuf0(-ch%n:ch%n))
        allocate(cbuf1(-ch%n:ch%n,dr%n_mode))
        cbuf0 = 0.0_r8
        cbuf1 = 0.0_r8
    end block init

    qploop: do ipt=1,dr%n_full_qpoint_local

        ! Need to know the global q-point index
        iq=dr%aq(ipt)%global_full_index

        ! Evaluate the relevant matrix element(s)
        matrixelement: block
            complex(r8) :: c0
            real(r8), dimension(3) :: qv1, qv2, qv3
            integer :: b1,b2,b3

            qv1=qpoint%r
            qv2=qp%ap(iq)%r
            qv3=-qv1-qv2

            if ( se%thirdorder_scattering ) then
                call op3%generate(fc2,p,mem,qvec=qv3)
                do b1=1,se%n_mode
                    op3%sigma(b1) = qp%adaptive_sigma( qp%ap(iq)%radius,op3%vel(:,b1),dr%default_smearing(b1),se%smearing_prefactor)
                enddo

                call pretransform_phi(fc3, qv2, qv3, ptf_phi)

                psisq_3ph=0.0_r8
                do b1 = 1, dr%n_mode
                do b2 = 1, dr%n_mode
                    evp1 = 0.0_r8
                    !call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv2(:, b2, q), 1, fh%ugv1(:, b1), 1, fh%evp1, dr%n_mode)
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), dr%aq(ipt)%nuvec(:,b2), 1, ompoint%nuvec(:, b1), 1, evp1, dr%n_mode)
                    do b3 = 1, dr%n_mode
                        evp2 = 0.0_r8
                        !call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), fh%ugv3(:, b3, q), 1, fh%evp1, 1, fh%evp2, dr%n_mode)
                        call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), op3%nuvec(:, b3), 1, evp1, 1, evp2, dr%n_mode)
                        evp2 = conjg(evp2)
                        c0=dot_product(evp2, ptf_phi)
                        psisq_3ph(b1, b2, b3)=abs(conjg(c0)*c0)
                    end do
                enddo
                enddo
            else
                psisq_3ph=0.0_r8
            endif

            if ( se%isotope_scattering ) then
                do b1=1,dr%n_mode
                do b2=1,dr%n_mode
                    psisq_isotope(b1,b2)=isotope_scattering_strength(p,ompoint%egv(:,b1),dr%aq(ipt)%egv(:,b2),ompoint%omega(b1),dr%aq(ipt)%omega(b2))
                enddo
                enddo
            else
                psisq_isotope=0.0_r8
            endif
        end block matrixelement

        ! Evaluate spectral functions
        spectralfunction: block
            real(r8), dimension(:,:), allocatable :: buf_j,buf_re,buf_im
            real(r8), dimension(3) :: qv1, qv2, qv3
            real(r8) :: f0
            integer :: imode

            qv1=qpoint%r
            qv2=qp%ap(iq)%r
            qv3=-qv1-qv2
            call mem%allocate(buf_re,[ise%n_energy,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_im,[ise%n_energy,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_j,[ise%n_energy,dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf_re=0.0_r8
            buf_im=0.0_r8
            buf_j=0.0_r8

            ! Spectral function for q'
            call ise%evaluate(p,qv2,dr%aq(ipt)%omega,dr%aq(ipt)%egv,buf_re,buf_im,mem)
            ! Turn self-energy into smeared normalized spectral functions
            buf_j=0.0_r8
            do imode=1,dr%n_mode
                if ( dr%aq(ipt)%omega(imode) .lt. lo_freqtol ) cycle
                call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),dr%aq(ipt)%omega(imode),buf_j(:,imode))
                call lo_gaussian_smear_spectral_function(ise%omega,dr%aq(ipt)%sigma(imode),buf_j(:,imode))
                f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
                buf_j(:,imode)=buf_j(:,imode)/f0
            enddo

            ! Pre-buffer in convolution helper. Means I store greater and lesser for q' in buffer 1.
            call ch%buffer_and_transform_J_to_greater_lesser(buf_j,1)

            ! Also pre-buffer the bare spectral without any thermal prefactors for isotope scattering.
            ! Stored in ch%lesser3(:,:)
            call ch%buffer_and_transform_J(buf_j,3)

            ! Spectral function for q''
            call ise%evaluate(p,qv3,op3%omega,op3%egv,buf_re,buf_im,mem)
            ! Turn self-energy into smeared normalized spectral functions in the time domain
            buf_j=0.0_r8
            do imode=1,dr%n_mode
                if ( op3%omega(imode) .lt. lo_freqtol ) cycle
                call lo_evaluate_spectral_function(ise%omega,buf_im(:,imode),buf_re(:,imode),op3%omega(imode),buf_j(:,imode))
                call lo_gaussian_smear_spectral_function(ise%omega,op3%sigma(imode),buf_j(:,imode))
                f0=lo_trapezoid_integration(ise%omega,buf_j(:,imode))
                buf_j(:,imode)=buf_j(:,imode)/f0
            enddo

            ! Pre-buffer in convolution helper. Means I store greater and lesser for q'' in buffer 2.
            call ch%buffer_and_transform_J_to_greater_lesser(buf_j,2)

            call mem%deallocate(buf_re,persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_im,persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_j,persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block spectralfunction

        frequencyintegral: block
            integer :: b1,b2,b3

            ! Add the prefactor to the matrix elements
            psisq_3ph=psisq_3ph*threephonon_prefactor*qp%ap(iq)%integration_weight
            psisq_isotope=psisq_isotope*isotope_prefactor*qp%ap(iq)%integration_weight

            ! Make this faster and loop only over half of b2
            if ( se%thirdorder_scattering ) then
                do b3 = 1, dr%n_mode
                do b2 = 1, dr%n_mode
                    cbuf0 = ch%greater1(:,b2)*ch%greater2(:,b3) - ch%lesser1(:,b2)*ch%lesser2(:,b3)
                    do b1=1, dr%n_mode
                        cbuf1(:,b1)=cbuf1(:,b1) + psisq_3ph(b1,b2,b3)*cbuf0
                    enddo
                enddo
                enddo
            endif

            if ( se%isotope_scattering ) then
                do b2=1, dr%n_mode
                do b1=1, dr%n_mode
                    cbuf1(:,b1)=cbuf1(:,b1) + psisq_isotope(b1,b2)*ch%lesser3(:,b2)
                enddo
                enddo
            endif
        end block frequencyintegral

        if (verbosity .gt. 0) then
            if (lo_trueNtimes(ipt, 127, dr%n_full_qpoint_local)) call lo_progressbar(' ... threephonon imaginary selfenergy', ipt, dr%n_full_qpoint_local)
        end if
    enddo qploop

    ! Sum it up and go back to frequency domain (cbuf1 is in the time domain right now)
    inversetransform: block
        integer :: imode
        ! Sync across ranks
        call mw%allreduce('sum', cbuf1)
        ! Inverse transform
        se%im=0.0_r8
        do imode=1,dr%n_mode
            if ( mod(imode,mw%n) .ne. mw%r ) cycle
            call ch%inverse_transform_to_real(cbuf1(:,imode),se%im(:,imode))
        enddo
        call mw%allreduce('sum', se%im)
    end block inversetransform

    ! And some cleanup
    call ch%destroy()
    deallocate(cbuf0)
    deallocate(cbuf1)
    call mem%deallocate(ptf_phi,       persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1,          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2,          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(psisq_3ph,     persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(psisq_isotope, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Fix degeneracies? Why not.
    fixdegen: block
        real(r8), dimension(:), allocatable :: buf
        integer :: b1,b2,i

        call mem%allocate(buf,ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf=0.0_r8

        do b1 = 1, dr%n_mode
            buf = 0.0_r8
            do i = 1, ompoint%degeneracy(b1)
                b2 = ompoint%degenmode(i, b1)
                buf = buf + se%im(:, b2)
            end do
            buf = buf/real(ompoint%degeneracy(b1), r8)
            do i = 1, ompoint%degeneracy(b1)
                b2 = ompoint%degenmode(i, b1)
                se%im(:, b2) = buf
            end do
        end do

        call mem%deallocate(buf,persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block fixdegen


    if (verbosity .gt. 0) call lo_progressbar(' ... threephonon imaginary selfenergy', dr%n_mode*qp%n_full_point, dr%n_mode*qp%n_full_point, walltime() - t0)

    ! if ( mw%talk ) write(*,*) 'Done here for now ',__FILE__,__LINE__
    ! call mw%destroy()
    ! stop
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

end submodule