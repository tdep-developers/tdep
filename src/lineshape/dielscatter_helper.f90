submodule(dielscatter) dielscatter_helper
!! Helper routines for Raman/IR stuff
implicit none
contains

!> build the helper array with all matrix elements and stuff
module subroutine buildhelper(dh, wp, di, p, qp, dr, fc, fct, temperature, deltaeps, mw, mem, verbosity)
    !> helper container
    class(diel_helper), intent(out) :: dh
    !> harmonic properties at this q-point (which should be Gamma)
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> interaction tensors
    type(lo_dielectric_tensor), intent(in) :: di
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> temperature
    real(r8), intent(in) :: temperature
    !> shift in dielectric constant from second order polarizability terms
    real(r8), dimension(:, :), intent(out) :: deltaeps
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8) :: t0, t1
    complex(r8), dimension(:, :), allocatable :: buf_ir_1
    complex(r8), dimension(:, :, :), allocatable :: buf_raman_1
    complex(r8), dimension(:, :, :, :), allocatable :: buf_ir_2
    complex(r8), dimension(:, :, :, :), allocatable :: buf_raman_2
    complex(r8), dimension(:, :, :, :), allocatable :: buf_phi
    real(r8), dimension(:, :, :), allocatable :: sym9, sym81

    t0 = walltime()
    t1 = t0

    init: block
        integer :: i

        ! Space for matrix elements
        allocate (dh%cmp_irI(9, dr%n_mode))
        allocate (dh%cmp_irI_sph(9, dr%n_mode))
        allocate (dh%cmp_irII(9, dr%n_mode, dr%n_mode, qp%n_irr_point))
        allocate (dh%cmp_irIII(9, dr%n_mode, dr%n_mode, dr%n_mode, qp%n_irr_point))
        allocate (dh%cmp_rmI(81, dr%n_mode))
        allocate (dh%cmp_rmII(81, dr%n_mode, dr%n_mode, qp%n_irr_point))
        dh%cmp_irI = 0.0_r8
        dh%cmp_irI_sph = 0.0_r8
        dh%cmp_irII = 0.0_r8
        dh%cmp_irIII = 0.0_r8
        dh%cmp_rmI = 0.0_r8
        dh%cmp_rmII = 0.0_r8

        ! Space for intermediate buffers
        call mem%allocate(buf_ir_1, [3, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_ir_2, [3, dr%n_mode, dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_raman_1, [3, 3, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_raman_2, [9, dr%n_mode, dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_phi, [dr%n_mode, dr%n_mode, dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_ir_1 = 0.0_r8
        buf_ir_2 = 0.0_r8
        buf_raman_1 = 0.0_r8
        buf_raman_2 = 0.0_r8
        buf_phi = 0.0_r8

        ! Rotation matrices
        call mem%allocate(sym9, [9, 9, p%sym%n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(sym81, [81, 81, p%sym%n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        sym9 = 0.0_r8
        sym81 = 0.0_r8
        do i = 1, p%sym%n
            sym9(:, :, i) = lo_expandoperation_pair(p%sym%op(i)%m)
            sym81(:, :, i) = lo_expandoperation_quartet(p%sym%op(i)%m)
        end do
    end block init

    ! Fiddle with eigenvectors a little bit
    fixeigenvectors: block
        real(r8), parameter :: isqrt2 = 1.0_r8/sqrt(2.0_r8)
        real(r8) :: f0, f1
        integer :: imode, iatom, ix, iq, ialpha, ctr

        allocate (dh%ugv_gamma(p%na*3, p%na*3))
        allocate (dh%ugv(p%na*3, p%na*3, qp%n_irr_point))
        dh%ugv_gamma = 0.0_r8
        dh%ugv = 0.0_r8

        ! eigenvectors at gamma, epsilon to nu transformation
        do imode = 1, dr%n_mode
        do iatom = 1, p%na
        do ix = 1, 3
            if (wp%omega(imode) .gt. lo_freqtol) then
                f0 = 1.0_r8/sqrt(wp%omega(imode))
            else
                f0 = 0.0_r8
            end if
            ialpha = (iatom - 1)*3 + ix
            f1 = p%invsqrtmass(iatom)
            dh%ugv_gamma(ialpha, imode) = wp%egv(ialpha, imode)*f0*f1*isqrt2
        end do
        end do
        end do

        ! eigenvectors at general q-point
        ctr = 0
        do iq = 1, qp%n_irr_point
        do imode = 1, dr%n_mode
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle

            if (dr%iq(iq)%omega(imode) .gt. lo_freqtol) then
                f0 = 1.0_r8/sqrt(dr%iq(iq)%omega(imode))
            else
                f0 = 0.0_r8
            end if

            do iatom = 1, p%na
                f1 = p%invsqrtmass(iatom)
                do ix = 1, 3
                    ialpha = (iatom - 1)*3 + ix
                    dh%ugv(ialpha, imode, iq) = dr%iq(iq)%egv(ialpha, imode)*f0*f1*isqrt2
                end do
            end do
        end do
        end do
        call mw%allreduce('sum', dh%ugv)

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... transformed eigenvectors (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block fixeigenvectors

    ! Get the lowest order matrix elements, these will be used in several places.
    order1: block
        integer :: imode

        buf_ir_1 = 0.0_r8
        buf_raman_1 = 0.0_r8
        do imode = 1, dr%n_mode
            buf_ir_1(:, imode) = di%polarization_matrixelement_firstorder(dh%ugv_gamma(:, imode))
            buf_raman_1(:, :, imode) = di%epsilon_matrixelement_firstorder(dh%ugv_gamma(:, imode))
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... got lowest order matrix elements (', tochar(t1 - t0), 's)'
            t0 = t1
        end if

        ! Perhaps get the spherically averaged guys? Could make sense.
        call spherically_averaged_firstorder_IR_matrixelements(fc, p, di, sym9, dh%cmp_irI_sph, mw, mem)

    end block order1

    ! second order matrix elements, a lot of them
    order2: block
        complex(r8), dimension(:, :), allocatable :: pretransform_rm2, pretransform_ir2
        complex(r8), dimension(:), allocatable :: pretransform_phi3
        complex(r8), dimension(:), allocatable :: outerprod1, outerprod2
        complex(r8) cv9(9), cv3(3)
        real(r8), dimension(3) :: qv1, qv2
        integer :: ctr, ib1, ib2, ib3, iq

        call mem%allocate(pretransform_rm2, [9, dr%n_mode**2], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(pretransform_ir2, [3, dr%n_mode**2], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(pretransform_phi3, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(outerprod1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(outerprod2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        pretransform_rm2 = 0.0_r8
        pretransform_ir2 = 0.0_r8
        pretransform_phi3 = 0.0_r8
        outerprod1 = 0.0_r8
        outerprod2 = 0.0_r8

        if (verbosity .gt. 0) then
            call lo_progressbar_init()
        end if

        ctr = 0
        do iq = 1, qp%n_irr_point
            ! Start pre-transforming.
            qv1 = qp%ip(iq)%r
            qv2 = -qv1

            ! Makes the most sense to have parallelism here, before the pretransforms to not do them many times
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle

            if (di%n_eps_pair .gt. 0) then
                call di%pretransform_epsilon_secondorder(qv1, pretransform_rm2)
            end if

            if (di%n_Z_pair .gt. 0) then
                call di%pretransform_Z_secondorder(qv1, pretransform_ir2)
            end if

            ! If either of the pair terms are nonzero we need the three-phonon guys
            if (di%n_Z_pair .gt. 0 .or. di%n_eps_pair .gt. 0) then
                call fct%pretransform(qv1, qv2, pretransform_phi3)
            end if

            ! Now actually calculate the matrix elements. Do the IR and Raman ones first.
            do ib1 = 1, dr%n_mode
            do ib2 = 1, dr%n_mode
                ! Now ... fancy outer product thing for Raman 2 and IR 2, basically
                ! nu1 \otimes nu2, with conjugation.
                outerprod1 = 0.0_r8
                call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), dh%ugv(:, ib2, iq), 1, dh%ugv(:, ib1, iq), 1, outerprod1, dr%n_mode)
                outerprod1 = conjg(outerprod1)
                ! And get the actual matrix elements
                call zgemv('N', 9, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), pretransform_rm2, 9, outerprod1, 1, (0.0_r8, 0.0_r8), cv9, 1)
                call zgemv('N', 3, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), pretransform_ir2, 3, outerprod1, 1, (0.0_r8, 0.0_r8), cv3, 1)
                buf_ir_2(:, ib1, ib2, iq) = cv3
                buf_raman_2(:, ib1, ib2, iq) = cv9
            end do
            end do

            ! And now the three-phonon guys. There is a chance that this is correct.
            do ib3 = 1, dr%n_mode ! mode at Gamma!
            do ib1 = 1, dr%n_mode ! mode at q
                outerprod1 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), dh%ugv(:, ib1, iq), 1, dh%ugv_gamma(:, ib3), 1, outerprod1, dr%n_mode)
                do ib2 = 1, dr%n_mode ! mode at -q
                    if (wp%omega(ib2) .lt. lo_freqtol) cycle
                    outerprod2 = 0.0_r8
                    call zgerc(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), dh%ugv(:, ib2, iq), 1, outerprod1, 1, outerprod2, dr%n_mode)
                    buf_phi(ib1, ib2, ib3, iq) = dot_product(outerprod2, pretransform_phi3)
                    ! ! Calculate stupid matrix elements to see if I got this somewhat right? Seems like it.
                    ! omega(1)=wp%omega(ib3)
                    ! omega(2)=dr%iq(iq)%omega(ib1)
                    ! omega(3)=dr%iq(iq)%omega(ib2)
                    ! if ( minval(omega) .lt. lo_freqtol ) cycle
                    !
                    ! egv(:,1)=wp%egv(:,ib3)
                    ! egv(:,2)=dr%iq(iq)%egv(:,ib1)
                    ! egv(:,3)=conjg(dr%iq(iq)%egv(:,ib2))
                    ! c0=fct%scatteringamplitude(omega,egv,-qv1*lo_twopi,-qv2*lo_twopi)
                    ! c1=buf_phi(ib1,ib2,ib3,iq)*6.0_r8*2.0_r8**(3.0_r8/2.0_r8)
                    !
                    ! if ( abs(c0) .lt. 1E-15_r8 ) cycle
                    !
                    ! write(*,*) 'bbb',ib1,ib2,ib3
                    ! write(*,*) c0
                    ! write(*,*) c1
                    ! write(*,*) c1/c0
                    ! write(*,*) c0/c1
                end do
            end do
            end do

            if (verbosity .gt. 0 .and. iq .lt. qp%n_irr_point) then
                call lo_progressbar(' ... second order matrix elements', iq, qp%n_irr_point, walltime() - t0)
            end if
        end do

        call mw%allreduce('sum', buf_ir_2)
        call mw%allreduce('sum', buf_raman_2)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... second order matrix elements', qp%n_irr_point, qp%n_irr_point, t1 - t0)
            t0 = t1
        end if

        call mem%deallocate(pretransform_rm2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(pretransform_ir2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(pretransform_phi3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(outerprod1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(outerprod2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block order2

    ! build the compound matrix elements
    compound: block
        real(r8), dimension(81) :: m81a, m81b, m81c
        real(r8), dimension(9) :: m9a, m9b, m9c
        integer :: iq, ctr
        integer :: ix, iy, iz, iu, ia, ib, ii, iop, i, ib1, ib2, ib3

        if (verbosity .gt. 0) then
            call lo_progressbar_init()
        end if

        ! First the easy lowest order guys
        dh%cmp_irI = 0.0_r8
        dh%cmp_rmI = 0.0_r8
        do ib1 = 1, dr%n_mode
            if (mod(ib1, mw%n) .ne. mw%r) cycle
            do ix = 1, 3
            do iy = 1, 3
                ii = (ix - 1)*3 + iy
                dh%cmp_irI(ii, ib1) = real(buf_ir_1(ix, ib1)*conjg(buf_ir_1(iy, ib1)), r8)
            end do
            end do

            do ix = 1, 3
            do iy = 1, 3
            do iz = 1, 3
            do iu = 1, 3
                ii = iu + (iz - 1)*3 + (iy - 1)*9 + (ix - 1)*27
                dh%cmp_rmI(ii, ib1) = real(buf_raman_1(ix, iy, ib1)*conjg(buf_raman_1(iz, iu, ib1)), r8)
            end do
            end do
            end do
            end do
        end do
        call mw%allreduce('sum', dh%cmp_irI)
        call mw%allreduce('sum', dh%cmp_rmI)

        ! Build all the annoying guys
        dh%cmp_rmII = 0.0_r8
        dh%cmp_irII = 0.0_r8
        ctr = 0
        do iq = 1, qp%n_irr_point

            do ib1 = 1, dr%n_mode
            do ib2 = 1, dr%n_mode
                ! Think about where to put MPI guy
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                if (di%n_eps_pair .gt. 0) then
                    do ix = 1, 3
                    do iy = 1, 3
                    do iz = 1, 3
                    do iu = 1, 3
                        ia = (ix - 1)*3 + iy
                        ib = (iz - 1)*3 + iu
                        ii = iu + (iz - 1)*3 + (iy - 1)*9 + (ix - 1)*27
                        m81a(ii) = real(buf_raman_2(ia, ib1, ib2, iq)*conjg(buf_raman_2(ib, ib1, ib2, iq)), r8)
                        m81b(ii) = real(conjg(buf_raman_2(ia, ib1, ib2, iq))*buf_raman_2(ib, ib1, ib2, iq), r8)
                    end do
                    end do
                    end do
                    end do

                    m81c = 0.0_r8
                    do i = 1, qp%ip(iq)%n_full_point
                        iop = qp%ip(iq)%operation_full_point(i)
                        ! No idea how this works. Uuuh.
                        if (iop .gt. 0) then
                            call dgemv('N', 81, 81, 1.0_r8, sym81(:, :, iop), 81, m81a, 1, 1.0_r8, m81c, 1)
                        else
                            call dgemv('N', 81, 81, 1.0_r8, sym81(:, :, -iop), 81, m81b, 1, 1.0_r8, m81c, 1)
                        end if
                    end do
                    dh%cmp_rmII(:, ib1, ib2, iq) = m81c/real(qp%ip(iq)%n_full_point, r8)
                end if

                if (di%n_Z_pair .gt. 0) then
                    do ix = 1, 3
                    do iy = 1, 3
                        ia = (ix - 1)*3 + iy
                        m9a(ia) = real(buf_ir_2(ix, ib1, ib2, iq)*conjg(buf_ir_2(iy, ib1, ib2, iq)), r8)
                        m9b(ia) = real(conjg(buf_ir_2(ix, ib1, ib2, iq))*buf_ir_2(iy, ib1, ib2, iq), r8)
                    end do
                    end do

                    m9c = 0.0_r8
                    do i = 1, qp%ip(iq)%n_full_point
                        iop = qp%ip(iq)%operation_full_point(i)
                        ! No idea how this works. Uuuh.
                        if (iop .gt. 0) then
                            call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, iop), 9, m9a, 1, 1.0_r8, m9c, 1)
                        else
                            call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, -iop), 9, m9b, 1, 1.0_r8, m9c, 1)
                        end if
                    end do
                    dh%cmp_irII(:, ib1, ib2, iq) = m9c/real(qp%ip(iq)%n_full_point, r8)
                end if

                ! Annoying guys with three indices
                if (di%n_Z_singlet .gt. 0 .and. di%n_Z_pair .gt. 0) then
                    do ib3 = 1, dr%n_mode ! mode at Gamma
                        if (wp%omega(ib3) .lt. lo_freqtol) cycle
                        do ix = 1, 3
                        do iy = 1, 3
                            ia = (ix - 1)*3 + iy
                            m9a(ia) = real(buf_ir_2(ix, ib1, ib2, iq)*buf_ir_1(iy, ib3), r8)
                            m9b(ia) = aimag(buf_ir_2(iy, ib1, ib2, iq)*buf_ir_1(ix, ib3)) !,r8)
                        end do
                        end do
                        m9c = 0.0_r8
                        do i = 1, qp%ip(iq)%n_full_point
                            iop = qp%ip(iq)%operation_full_point(i)
                            ! No idea how this works. Uuuh.
                            if (iop .gt. 0) then
                                call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, iop), 9, m9a, 1, 1.0_r8, m9c, 1)
                            else
                                call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, -iop), 9, m9b, 1, 1.0_r8, m9c, 1)
                            end if
                        end do
                        m9c = m9c*real(buf_phi(ib1, ib2, ib3, iq), r8)/real(qp%ip(iq)%n_full_point, r8)
                        dh%cmp_irIII(:, ib1, ib2, ib3, iq) = m9c

                    end do
                end if
            end do
            end do

            if (verbosity .gt. 0 .and. iq .lt. qp%n_irr_point) then
                call lo_progressbar(' ... compound matrix elements', iq, qp%n_irr_point, walltime() - t0)
            end if
        end do

        ! Sync across ranks
        call mw%allreduce('sum', dh%cmp_rmII)
        call mw%allreduce('sum', dh%cmp_irII)
        call mw%allreduce('sum', dh%cmp_irIII)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... compound matrix elements', qp%n_irr_point, qp%n_irr_point, t1 - t0)
            t0 = t1
        end if
    end block compound

    ! While I have the matrix elements, calculate the thermal contribution
    thermalavg: block
        real(r8), dimension(9) :: v0, v1, w0
        real(r8) :: n
        integer :: iq, imode, ctr, i, iop

        if (verbosity .gt. 0) call lo_progressbar_init()

        ctr = 0
        w0 = 0.0_r8
        do iq = 1, qp%n_irr_point
        do imode = 1, dr%n_mode
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            n = lo_planck(temperature, dr%iq(iq)%omega(imode))
            v0 = 0.0_r8
            v1 = real(buf_raman_2(:, imode, imode, iq), r8)
            do i = 1, qp%ip(iq)%n_full_point
                iop = qp%ip(iq)%operation_full_point(i)
                ! No idea how this works. Uuuh.
                if (iop .gt. 0) then
                    call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, iop), 9, v1, 1, 1.0_r8, v0, 1)
                else
                    call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, -iop), 9, v1, 1, 1.0_r8, v0, 1)
                end if
            end do
            v0 = v0/real(qp%ip(iq)%n_full_point, r8)
            w0 = w0 + (2*n + 1)*v0*qp%ip(iq)%integration_weight
            if (verbosity .gt. 0 .and. ctr .lt. qp%n_irr_point*dr%n_mode) then
                call lo_progressbar(' ... thermal dielectric constant', ctr, qp%n_irr_point*dr%n_mode, walltime() - t0)
            end if
        end do
        end do

        call mw%allreduce('sum', w0)
        ! Average over symmetry operations to be on the safe side.
        v0 = 0.0_r8
        do iop = 1, p%sym%n
            call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, iop), 9, w0, 1, 1.0_r8, v0, 1)
        end do
        v0 = lo_chop(v0/real(p%sym%n, r8), 1E-12_r8)

        deltaeps = lo_unflatten_2tensor(v0)
        !m0=I3-m0*4*lo_pi/volume
        ! Convert to raw derivatives to epsilon?
        ! I do not add the +1 here, since otherwise it will be double-counted.
        ! it is already included in the baseline.
        deltaeps = -deltaeps*4*lo_pi/p%volume
        !do i=1,3
        !    deltaeps(i,i)=deltaeps(i,i)+1.0_r8
        !enddo

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... thermal dielectric constant', qp%n_irr_point, qp%n_irr_point, t1 - t0)
            t0 = t1
        end if

    end block thermalavg

    ! And some cleanup at the end
    call mem%deallocate(buf_ir_1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_ir_2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_raman_1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_raman_2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_phi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(sym9, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(sym81, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

subroutine spherically_averaged_firstorder_IR_matrixelements(fc, p, di, sym9, cmp, mw, mem)
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> dielectric stuff
    type(lo_dielectric_tensor), intent(in) :: di
    !> how many points for the spherical average
    !integer, intent(in) :: npts
    !> symmetry operations
    real(r8), dimension(:, :, :), intent(in) :: sym9
    !> first order compound matrix elements (xyz^2,n_mode)
    real(r8), dimension(:, :), intent(out) :: cmp
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer, parameter :: npts = 1000 ! Seems to be more than enough.
    type(lo_phonon_dispersions_qpoint) :: wp
    real(r8), dimension(:, :), allocatable :: buf_cmp1, buf_cmp2
    real(r8), dimension(:, :), allocatable :: qdir
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(9) :: w0
    real(r8), dimension(3) :: v0
    integer :: iq, imode, jmode, iatom, ix, ialpha, i, iop

    ! Get some points on a sphere
    call mem%allocate(qdir, [3, npts], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_cmp1, [9, p%na*3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_cmp2, [9, p%na*3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    qdir = 0.0_r8
    buf_cmp1 = 0.0_r8
    buf_cmp2 = 0.0_r8
    if (mw%r .eq. mw%n - 1) then
        call lo_points_on_sphere(qdir)
    end if
    call mw%bcast(qdir, from=mw%n - 1)

    ! Start averaging
    do iq = 1, npts
        if (mod(iq, mw%n) .ne. mw%r) cycle

        call wp%generate(fc, p, mem, qvec=[0.0_r8, 0.0_r8, 0.0_r8], qdirection=qdir(:, iq))

        buf_cmp1 = 0.0_r8
        do imode = 1, p%na*3
            if (wp%omega(imode) .lt. lo_freqtol) cycle
            ! Scale eigenvector
            do iatom = 1, p%na
            do ix = 1, 3
                ialpha = (iatom - 1)*3 + ix
                wp%egv(ialpha, imode) = wp%egv(ialpha, imode)*p%invsqrtmass(iatom)/sqrt(2.0_r8*wp%omega(imode))
            end do
            end do
            ! Get matrix element
            v0 = real(di%polarization_matrixelement_firstorder(wp%egv(:, imode)), r8)
            m0 = lo_outerproduct(v0, v0)
            buf_cmp1(:, imode) = lo_flattentensor(m0)
        end do
        ! Take care of degeneracies?
        do imode = 1, p%na*3
            if (wp%omega(imode) .lt. lo_freqtol) cycle
            w0 = 0.0_r8
            do i = 1, wp%degeneracy(imode)
                jmode = wp%degenmode(i, imode)
                w0 = w0 + buf_cmp1(:, jmode)/real(wp%degeneracy(imode), r8)
            end do
            do i = 1, wp%degeneracy(imode)
                jmode = wp%degenmode(i, imode)
                buf_cmp1(:, jmode) = w0
            end do
        end do
        ! Accumulate
        buf_cmp2 = buf_cmp2 + buf_cmp1
    end do

    ! Sync over ranks
    call mw%allreduce('sum', buf_cmp2)
    buf_cmp2 = buf_cmp2/real(npts, r8)

    ! Make sure we obey the proper symmetries?
    do imode = 1, p%na*3
        w0 = 0.0_r8
        do iop = 1, p%sym%n
            call dgemv('N', 9, 9, 1.0_r8, sym9(:, :, iop), 9, buf_cmp2(:, imode), 1, 1.0_r8, w0, 1)
        end do
        buf_cmp2(:, imode) = lo_chop(w0/real(p%sym%n, r8), 1E-10_r8)
    end do
    cmp = buf_cmp2

    ! And cleanup
    call mem%deallocate(qdir, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_cmp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_cmp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

end submodule
