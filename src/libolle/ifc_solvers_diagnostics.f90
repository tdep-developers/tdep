submodule(ifc_solvers) ifc_solvers_diagnostics
!!
!! Figure out some diagnostics related to the solution
!!
implicit none
contains

!> print a lot of diagnostics
module subroutine report_diagnostics(map, tp, ih, dh, mw, mem, verbosity)
    !> forcemap
    type(lo_forcemap), intent(in) :: map
    !> division of things
    type(lo_trainpred_helper), intent(in) :: tp
    !> force constant things
    type(lo_ifc_helper), intent(inout) :: ih
    !> dielectric things
    type(lo_dielectric_helper), intent(inout) :: dh
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk?
    integer, intent(in) :: verbosity

    ! check if the forceconstants make sense in the cross-validation sense.
    crossvalidation_fc: block
        real(r8), dimension(n_division, 4) :: force_rsq
        real(r8), dimension(n_division, 5) :: energy_dev
        real(r8), dimension(:, :), allocatable :: wA
        real(r8), dimension(:), allocatable :: wB, wC
        real(r8), dimension(:), allocatable :: f0, f1, f2, f3, f4, e0, e1, e2, e3, e4, u0, fp, ep
        real(r8) :: energypref
        integer :: div, nrow, i, nt, ii, jj
        character(len=1000) :: opf

        if (verbosity .gt. 0) then
            write (lo_iou, *) '... diagnostics from three-fold cross-validation:'
        end if

        ! Some temporary space
        call mem%allocate(wB, map%n_atom_ss*3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        call mem%allocate(wC, map%n_atom_ss*3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        wB = 0.0_r8
        wC = 0.0_r8

        force_rsq = 0.0_r8
        energy_dev = 0.0_r8
        do div = 1, n_division
            ! number of steps in this division
            nt = tp%nstep_pred(div)
            ! some space
            if (nt .gt. 0) then
                call mem%allocate(f1, nt*map%n_atom_ss*3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(f2, nt*map%n_atom_ss*3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(f3, nt*map%n_atom_ss*3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(f4, nt*map%n_atom_ss*3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e1, nt, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e2, nt, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e3, nt, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e4, nt, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                f1 = 0.0_r8
                f2 = 0.0_r8
                f3 = 0.0_r8
                f4 = 0.0_r8
                e1 = 0.0_r8
                e2 = 0.0_r8
                e3 = 0.0_r8
                e4 = 0.0_r8
            else
                call mem%allocate(f1, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(f2, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(f3, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(f4, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e1, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e2, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e3, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%allocate(e4, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                f1 = 0.0_r8
                f2 = 0.0_r8
                f3 = 0.0_r8
                f4 = 0.0_r8
                e1 = 0.0_r8
                e2 = 0.0_r8
                e3 = 0.0_r8
                e4 = 0.0_r8
            end if

            ! fetch the baseline, first the raw data
            call tp%collect(map, 0, div, .false., nrow, mem, dh, ih, values=f0, energies=e0, displacements=u0)

            ! then the polar data
            call tp%collect(map, -1, div, .false., nrow, mem, dh, ih, values=fp, energies=ep)

            ! collect forces
            if (map%have_fc_pair) then
                call tp%collect(map, 2, div, .false., nrow, mem, dh, ih, coeff=wA)
                if (nt .gt. 0) then
                    call lo_gemv(wA, ih%dx2(:, div), f2)
                    do i = 1, nt
                        ii = (i - 1)*3*map%n_atom_ss + 1
                        jj = i*3*map%n_atom_ss
                        wB = u0(ii:jj)
                        wC = f2(ii:jj)
                        e2(i) = -dot_product(wB, wC)*0.5_r8
                    end do
                end if
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                ! force predictive R^2
                force_rsq(div, 1) = distributed_rsquare(f2, f0, nrow, mw)
                force_rsq(div, 2) = distributed_rsquare(f2 + fp, f0, nrow, mw)
                ! energy standard deviation
                energy_dev(div, 1) = distributed_stddev(e0, nt, mw)
                energy_dev(div, 2) = distributed_stddev(e0 - e2, nt, mw)
                energy_dev(div, 3) = distributed_stddev(e0 - e2 - ep, nt, mw)
            end if

            if (map%have_fc_triplet) then
                call tp%collect(map, 3, div, .false., nrow, mem, dh, ih, coeff=wA)
                if (nt .gt. 0) then
                    call lo_gemv(wA, ih%dx3(:, div), f3)
                    do i = 1, nt
                        ii = (i - 1)*3*map%n_atom_ss + 1
                        jj = i*3*map%n_atom_ss
                        wB = u0(ii:jj)
                        wC = f3(ii:jj)
                        e3(i) = -dot_product(wB, wC)/3.0_r8
                    end do
                end if
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                ! force predictive R^2
                force_rsq(div, 3) = distributed_rsquare(f2 + fp + f3, f0, nrow, mw)
                ! energy deviation
                energy_dev(div, 4) = distributed_stddev(e0 - e2 - ep - e3, nt, mw)
            end if

            if (map%have_fc_quartet) then
                call tp%collect(map, 4, div, .false., nrow, mem, dh, ih, coeff=wA)
                if (nt .gt. 0) then
                    call lo_gemv(wA, ih%dx4(:, div), f4)
                    do i = 1, nt
                        ii = (i - 1)*3*map%n_atom_ss + 1
                        jj = i*3*map%n_atom_ss
                        wB = u0(ii:jj)
                        wC = f4(ii:jj)
                        e4(i) = -dot_product(wB, wC)*0.25_r8
                    end do
                end if
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                ! force predictive R^2
                force_rsq(div, 4) = distributed_rsquare(f2 + fp + f3 + f4, f0, nrow, mw)
                ! energy deviation
                energy_dev(div, 5) = distributed_stddev(e0 - e2 - ep - e3 - e4, nt, mw)
            end if

            ! cleanup
            call mem%deallocate(f1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(f2, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(f3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(f4, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(e1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(e2, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(e3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(e4, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)

            call mem%deallocate(f0, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(e0, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(u0, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(fp, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(ep, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
        end do

        call mem%deallocate(wB, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        call mem%deallocate(wC, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)

        ! report what we get
        if (verbosity .gt. 0) then
            energypref = lo_Hartree_to_eV*1000.0_r8/real(map%n_atom_ss, r8)
            opf = "(1X,A37,1X,A28)"
            write (lo_iou, opf) 'R^2', 'stddev energy (meV/atom)'

            opf = "(1X,A30,1X,F8.5,1X,F16.7)"
            if (map%have_fc_pair) then
                write (lo_iou, '(1X,A30,10X,F16.7)') '       input:', lo_mean(energy_dev(:, 1))*energypref
                write (lo_iou, opf) 'second order:', lo_mean(force_rsq(:, 1)), lo_mean(energy_dev(:, 2))*energypref
            end if
            if (map%have_fc_pair .and. map%polar .gt. 0 .and. map%xuc%nx_Z_singlet .gt. 0) then
                write (lo_iou, opf) 'second order + polar:', lo_mean(force_rsq(:, 2)), lo_mean(energy_dev(:, 3))*energypref
            end if
            if (map%have_fc_triplet) then
                write (lo_iou, opf) 'third order:', lo_mean(force_rsq(:, 3)), lo_mean(energy_dev(:, 4))*energypref
            end if
            if (map%have_fc_quartet) then
                write (lo_iou, opf) 'fourth order:', lo_mean(force_rsq(:, 4)), lo_mean(energy_dev(:, 5))*energypref
            end if
        end if
    end block crossvalidation_fc

    ! check how the born charges and dielectric constant thing went
    if (map%polar .gt. 1) then
        crossvalidation_diel: block
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: epsB, eps0, eps1, eps2, z0, z1, z2, z3
            real(r8), dimension(n_division, 4) :: Z_rsq, Z_dev, eps_rsq, eps_dev
            integer :: nrow, div, nt
            character(len=1000) :: opf

            Z_rsq = 0.0_r8
            Z_dev = 0.0_r8
            eps_rsq = 0.0_r8
            eps_dev = 0.0_r8
            do div = 1, n_division
                ! number of steps in this division
                nt = tp%nstep_pred(div)
                ! some space
                if (nt .gt. 0) then
                    call mem%allocate(z1, nt*map%n_atom_ss*9, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(z2, nt*map%n_atom_ss*9, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(z3, nt*map%n_atom_ss*9, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(eps0, nt*9, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(eps1, nt*9, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(eps2, nt*9, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                else
                    call mem%allocate(z1, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(z2, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(z3, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(eps0, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(eps1, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                    call mem%allocate(eps2, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                end if
                z1 = 0.0_r8
                z2 = 0.0_r8
                z3 = 0.0_r8
                eps0 = 0.0_r8
                eps1 = 0.0_r8
                eps2 = 0.0_r8
                ! fetch raw data
                call tp%collect(map, -2, div, .false., nrow, mem, dh, ih, values=epsB)
                call tp%collect(map, -3, div, .false., nrow, mem, dh, ih, values=z0)

                ! check Z
                if (map%xuc%nx_Z_singlet .gt. 0) then
                    call tp%collect(map, 21, div, .false., nrow, mem, dh, ih, coeff=wA)
                    if (nt .gt. 0) then
                        call lo_gemv(wA, dh%dx_Z1(:, div), z1)
                    end if
                    call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                    ! Z predictive R^2
                    Z_rsq(div, 2) = distributed_rsquare(z1, z0, nrow, mw)
                    ! energy standard deviation
                    Z_dev(div, 2) = distributed_stddev(z0 - z1, nrow, mw)
                end if
                if (map%have_Z_pair) then
                    call tp%collect(map, 22, div, .false., nrow, mem, dh, ih, coeff=wA)
                    if (nt .gt. 0) then
                        call lo_gemv(wA, dh%dx_Z2(:, div), z2)
                    end if
                    call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                    ! Z predictive R^2
                    Z_rsq(div, 3) = distributed_rsquare(z1 + z2, z0, nrow, mw)
                    ! energy standard deviation
                    Z_dev(div, 3) = distributed_stddev(z0 - z1 - z2, nrow, mw)
                end if
                if (map%have_Z_triplet) then
                    call tp%collect(map, 23, div, .false., nrow, mem, dh, ih, coeff=wA)
                    if (nt .gt. 0) then
                        call lo_gemv(wA, dh%dx_Z3(:, div), z3)
                    end if
                    call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                    ! Z predictive R^2
                    Z_rsq(div, 4) = distributed_rsquare(z1 + z2 + z3, z0, nrow, mw)
                    ! energy standard deviation
                    Z_dev(div, 4) = distributed_stddev(z0 - z1 - z2 - z3, nrow, mw)
                end if

                ! check eps
                if (map%polar .gt. 1) then
                    call tp%collect(map, 10, div, .false., nrow, mem, dh, ih, coeff=wA)
                    if (nt .gt. 0) then
                        call lo_gemv(wA, dh%dx_eps0(:, div), eps0)
                    end if
                    call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                    ! Z predictive R^2
                    eps_rsq(div, 2) = distributed_rsquare(eps0, epsB, nrow, mw)
                    ! energy standard deviation
                    eps_dev(div, 2) = distributed_stddev(epsB - eps0, nrow, mw)
                end if
                if (map%have_eps_singlet) then
                    call tp%collect(map, 11, div, .false., nrow, mem, dh, ih, coeff=wA)
                    if (nt .gt. 0) then
                        call lo_gemv(wA, dh%dx_eps1(:, div), eps1)
                    end if
                    call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                    ! Z predictive R^2
                    !eps_rsq(div,3)=distributed_rsquare(eps0+eps1,epsB,nrow,mw)
                    eps_rsq(div, 3) = distributed_rsquare(eps1, epsB - eps0, nrow, mw)
                    ! energy standard deviation
                    eps_dev(div, 3) = distributed_stddev(epsB - eps0 - eps1, nrow, mw)
                end if
                if (map%have_eps_pair) then
                    call tp%collect(map, 12, div, .false., nrow, mem, dh, ih, coeff=wA)
                    if (nt .gt. 0) then
                        call lo_gemv(wA, dh%dx_eps2(:, div), eps2)
                    end if
                    call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                    ! Z predictive R^2
                    eps_rsq(div, 4) = distributed_rsquare(eps1 + eps2, epsB - eps0, nrow, mw)
                    ! energy standard deviation
                    eps_dev(div, 4) = distributed_stddev(epsB - eps0 - eps1 - eps2, nrow, mw)
                end if

                call mem%deallocate(z1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(z2, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(z3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(eps0, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(eps1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(eps2, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)

                Z_dev(div, 1) = distributed_stddev(z0, nrow, mw)
                eps_dev(div, 1) = distributed_stddev(epsB, nrow, mw)
                call mem%deallocate(z0, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(epsB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)

            end do

            ! report what we get
            if (verbosity .gt. 0) then
                ! opf="(1X,A20,5X,A5,5X,A15)"
                ! write(lo_iou,opf) 'Born charge','R^2','stddev'
                opf = "(1X,A30,1X,F8.5,1X,F16.7)"
                write (lo_iou, '(1X,A30,10X,F16.7)') '       input:', lo_mean(Z_dev(:, 1))
                if (map%xuc%nx_Z_singlet .gt. 0) then
                    write (lo_iou, opf) 'Z singlet:', lo_mean(Z_rsq(:, 2)), lo_mean(Z_dev(:, 2))
                end if
                if (map%have_Z_pair) then
                    write (lo_iou, opf) 'Z pair:', lo_mean(Z_rsq(:, 3)), lo_mean(Z_dev(:, 3))
                end if
                if (map%have_Z_triplet) then
                    write (lo_iou, opf) 'Z triplet:', lo_mean(Z_rsq(:, 4)), lo_mean(Z_dev(:, 4))
                end if

                ! opf="(1X,A23,3X,A5,5X,A15)"
                ! write(lo_iou,opf) 'Dielectric constant','R^2','stddev'
                !opf="(1X,A25,1X,F8.5,1X,F16.7)"
                write (lo_iou, '(1X,A30,10X,F16.7)') '       input:', lo_mean(eps_dev(:, 1))
                write (lo_iou, opf) 'eps global:', lo_mean(eps_rsq(:, 2)), lo_mean(eps_dev(:, 2))
                if (map%have_eps_singlet) then
                    write (lo_iou, opf) 'eps singlet:', lo_mean(eps_rsq(:, 3)), lo_mean(eps_dev(:, 3))
                end if
                if (map%have_eps_pair) then
                    write (lo_iou, opf) 'eps pair:', lo_mean(eps_rsq(:, 4)), lo_mean(eps_dev(:, 4))
                end if

            end if
        end block crossvalidation_diel
    end if

    ! ! Get the statistics of the fit?
    ! fitst: block
    !     real(r8) :: rsq
    !
    !     rsq=distributed_rsquare(dh%eps0,dh%epso,tp%nt*9,mw)
    !     if ( verbosity .gt. 0 ) then
    !         write(*,*) 'R^2 eps 0',rsq
    !     endif
    !     rsq=distributed_rsquare(dh%eps1,dh%epso-dh%eps0,tp%nt*9,mw)
    !     if ( verbosity .gt. 0 ) then
    !         write(*,*) 'R^2 eps 1',rsq
    !     endif
    ! end block fitst

end subroutine

end submodule
