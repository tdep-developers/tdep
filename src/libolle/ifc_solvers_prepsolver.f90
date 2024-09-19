#include "precompilerdefinitions"
submodule (ifc_solvers) ifc_solvers_prepsolver
    use konstanter, only: lo_pi
    use type_forcemap, only: lo_coeffmatrix_singlet,lo_coeffmatrix_pair,lo_coeffmatrix_triplet,lo_coeffmatrix_quartet,&
                             lo_coeffmatrix_eps_pair, lo_coeffmatrix_eps_singlet, lo_coeffmatrix_supercell_z_singlet, &
                             lo_coeffmatrix_z_pair, lo_coeffmatrix_z_triplet
    use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
    use type_blas_lapack_wrappers, only: lo_gemm
    use gottochblandat, only: lo_linear_least_squares
! Collection of routines that rearrange data for the force constant solvers.
implicit none
contains

!> Create all the dielectric coefficient matrices, in parallel
module subroutine coefficients_dielectric(dh, tp, map, sim, ss, mw, mem, verbosity)
    !> coefficient matrices and other useful things
    type(lo_dielectric_helper), intent(out) :: dh
    !> division helper
    type(lo_trainpred_helper), intent(in) :: tp
    !> forcemap
    type(lo_forcemap), intent(in) :: map
    !> MD simulation
    type(lo_mdsim), intent(in) :: sim
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk?
    integer, intent(in) :: verbosity

    real(r8) :: timer, t0, t1

    ! start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! set up some things
    init: block
        ! Start making space for coefficient matrices and auxiliary things.
        ! The global epsilon and normal Born charges should always be there.
        allocate (dh%dx_eps0(map%xuc%nx_eps_global, n_division))
        dh%dx_eps0 = 0.0_r8
        if (map%xuc%nx_eps_singlet .gt. 0) then
            allocate (dh%dx_eps1(map%xuc%nx_eps_singlet, n_division))
            dh%dx_eps1 = 0.0_r8
        end if
        if (map%xuc%nx_eps_pair .gt. 0) then
            allocate (dh%dx_eps2(map%xuc%nx_eps_pair, n_division))
            dh%dx_eps2 = 0.0_r8
        end if
        if (map%xuc%nx_Z_singlet .gt. 0) then
            allocate (dh%dx_Z1(map%xuc%nx_Z_singlet, n_division))
            dh%dx_Z1 = 0.0_r8
        end if
        if (map%xuc%nx_Z_pair .gt. 0) then
            allocate (dh%dx_Z2(map%xuc%nx_Z_pair, n_division))
            dh%dx_Z2 = 0.0_r8
        end if
        if (map%xuc%nx_Z_triplet .gt. 0) then
            allocate (dh%dx_Z3(map%xuc%nx_Z_triplet, n_division))
            dh%dx_Z3 = 0.0_r8
        end if

        if (tp%nt .gt. 0) then
            allocate (dh%cm_eps0(tp%nt*9, map%xuc%nx_eps_global))
            dh%cm_eps0 = 0.0_r8
            if (map%have_eps_singlet) then
                allocate (dh%cm_eps1(tp%nt*9, map%xuc%nx_eps_singlet))
                dh%cm_eps1 = 0.0_r8
            end if
            if (map%have_eps_pair) then
                allocate (dh%cm_eps2(tp%nt*9, map%xuc%nx_eps_pair))
                dh%cm_eps2 = 0.0_r8
            end if

            if (map%xuc%nx_Z_singlet .gt. 0) then
                allocate (dh%cm_Z1(tp%nt*ss%na*9, map%xuc%nx_Z_singlet))
                dh%cm_Z1 = 0.0_r8
            end if
            if (map%have_Z_pair) then
                allocate (dh%cm_Z2(tp%nt*ss%na*9, map%xuc%nx_Z_pair))
                dh%cm_Z2 = 0.0_r8
            end if
            if (map%have_Z_triplet) then
                allocate (dh%cm_Z3(tp%nt*ss%na*9, map%xuc%nx_Z_triplet))
                dh%cm_Z3 = 0.0_r8
            end if
        else
            allocate (dh%cm_eps0(1, 1))
            allocate (dh%cm_eps1(1, 1))
            allocate (dh%cm_eps2(1, 1))
            allocate (dh%cm_Z1(1, 1))
            allocate (dh%cm_Z2(1, 1))
            allocate (dh%cm_Z3(1, 1))
            dh%cm_eps0 = -lo_huge
            dh%cm_eps1 = -lo_huge
            dh%cm_eps2 = -lo_huge
            dh%cm_Z1 = -lo_huge
            dh%cm_Z2 = -lo_huge
            dh%cm_Z3 = -lo_huge
        end if

        ! Make space for values
        if (tp%nt .gt. 0) then
            allocate (dh%Z0(tp%nt*ss%na*9))
            allocate (dh%Z1(tp%nt*ss%na*9))
            allocate (dh%Z2(tp%nt*ss%na*9))
            allocate (dh%Z3(tp%nt*ss%na*9))
            allocate (dh%epso(tp%nt*9))
            allocate (dh%eps0(tp%nt*9))
            allocate (dh%eps1(tp%nt*9))
            allocate (dh%eps2(tp%nt*9))
            dh%Z0 = 0.0_r8
            dh%Z1 = 0.0_r8
            dh%Z2 = 0.0_r8
            dh%Z3 = 0.0_r8
            dh%epso = 0.0_r8
            dh%eps0 = 0.0_r8
            dh%eps1 = 0.0_r8
            dh%eps2 = 0.0_r8
        else
            allocate (dh%Z0(1))
            allocate (dh%Z1(1))
            allocate (dh%Z2(1))
            allocate (dh%Z3(1))
            allocate (dh%epso(1))
            allocate (dh%eps0(1))
            allocate (dh%eps1(1))
            allocate (dh%eps2(1))
            dh%Z0 = -lo_huge
            dh%Z1 = -lo_huge
            dh%Z2 = -lo_huge
            dh%Z3 = -lo_huge
            dh%epso = -lo_huge
            dh%eps0 = -lo_huge
            dh%eps1 = -lo_huge
            dh%eps2 = -lo_huge
        end if
    end block init

    coeff: block
        real(r8), dimension(:, :), allocatable :: CMZ1, CMZ2, CMZ3
        real(r8), dimension(:, :), allocatable :: CMeps0, CMeps1, CMeps2
        real(r8), dimension(3, 3) :: m0
        integer :: i, j, k, l, ii, jj, nf, t, tt

        ! These two always exist
        if (map%xuc%nx_Z_singlet .gt. 0) then
            call mem%allocate(CMZ1, [9*ss%na, map%xuc%nx_Z_singlet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            CMZ1 = 0.0_r8
        end if

        if (map%have_Z_pair) then
            call mem%allocate(CMZ2, [9*ss%na, map%xuc%nx_Z_pair], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            CMZ2 = 0.0_r8
        end if
        if (map%have_Z_triplet) then
            call mem%allocate(CMZ3, [9*ss%na, map%xuc%nx_Z_triplet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            CMZ3 = 0.0_r8
        end if

        call mem%allocate(CMeps0, [9, map%xuc%nx_eps_global], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        CMeps0 = 0.0_r8

        if (map%have_eps_singlet) then
            call mem%allocate(CMeps1, [9, map%xuc%nx_eps_singlet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            CMeps1 = 0.0_r8
        end if
        if (map%have_eps_pair) then
            call mem%allocate(CMeps2, [9, map%xuc%nx_eps_pair], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            CMeps2 = 0.0_r8
        end if

        ! Start creating the coefficients
        if (verbosity .gt. 0) call lo_progressbar_init()
        tsloop: do tt = 1, tp%nt
            ! index to timestep (tt is a local index)
            t = tp%timestep_ind(tt)
            ! dielectric tensor
            CMeps0 = map%eps_global_shell%coeff
            ! born charges
            if (map%xuc%nx_Z_singlet .gt. 0) then
                call lo_coeffmatrix_supercell_Z_singlet(map, ss, CMZ1)
            end if
            if (map%have_Z_pair) then
                call lo_coeffmatrix_Z_pair(map, sim%u(:, :, t), CMZ2)
            end if
            if (map%have_Z_triplet) then
                call lo_coeffmatrix_Z_triplet(map, sim%u(:, :, t), CMZ3)
            end if
            if (map%have_eps_singlet) then
                call lo_coeffmatrix_eps_singlet(map, sim%u(:, :, t), CMeps1)
            end if
            if (map%have_eps_pair) then
                call lo_coeffmatrix_eps_pair(map, sim%u(:, :, t), CMeps2)
            end if

            ! Store, first epsilon things.
            ! Odd observation, in the simulation I have eps stored, but I think
            ! the reasonable thing to work with is derivatives with respect to
            ! electric field. It makes sense to convert it.
            m0 = sim%eps(:, :, t)
            do i = 1, 3
                m0(i, i) = 1.0_r8 - m0(i, i)
            end do
            m0 = m0*sim%crystalstructure%volume*0.25_r8/lo_pi

            nf = 9
            ii = (tt - 1)*nf + 1
            jj = tt*nf
            do j = 1, 3
            do i = 1, 3
                l = (tt - 1)*nf + (j - 1)*3 + i
                !dh%epso(l)=sim%eps(i,j,t)*tp%weight(tt)
                dh%epso(l) = m0(i, j)*tp%weight(tt)
            end do
            end do
            dh%cm_eps0(ii:jj, :) = CMeps0*tp%weight(tt)
            if (map%have_eps_singlet) then
                dh%cm_eps1(ii:jj, :) = CMeps1*tp%weight(tt)
            end if
            if (map%have_eps_pair) then
                dh%cm_eps2(ii:jj, :) = CMeps2*tp%weight(tt)
            end if

            ! Then store borncharge things
            nf = 9*ss%na
            ii = (tt - 1)*nf + 1
            jj = tt*nf
            if (map%xuc%nx_Z_singlet .gt. 0) then
                dh%cm_Z1(ii:jj, :) = CMZ1*tp%weight(tt)
            end if
            if (map%have_Z_pair) then
                dh%cm_Z2(ii:jj, :) = CMZ2*tp%weight(tt)
            end if
            if (map%have_Z_triplet) then
                dh%cm_Z3(ii:jj, :) = CMZ3*tp%weight(tt)
            end if

            ! store measured values of Z
            do k = 1, ss%na
            do j = 1, 3
            do i = 1, 3
                ii = ii + 1
                l = (tt - 1)*nf + (k - 1)*9 + (j - 1)*3 + i
                dh%Z0(l) = sim%Z(i, j, k, t)*tp%weight(tt)
            end do
            end do
            end do

            if (verbosity .gt. 0 .and. tt .lt. tp%nt) then
                call lo_progressbar(' ... creating polar coefficients', tt, tp%nt, walltime() - timer)
            end if
        end do tsloop

        if (map%xuc%nx_Z_singlet .gt. 0) then
            call mem%deallocate(CMZ1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (map%have_Z_pair) then
            call mem%deallocate(CMZ2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (map%have_Z_triplet) then
            call mem%deallocate(CMZ3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        call mem%deallocate(CMeps0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (map%have_eps_singlet) then
            call mem%deallocate(CMeps1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (map%have_eps_pair) then
            call mem%deallocate(CMeps2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        ! for sensible timings
        call mw%barrier()

        if (verbosity .gt. 0) then
            call lo_progressbar(' ... creating polar coefficients', tp%nt, tp%nt, walltime() - timer)
        end if
    end block coeff
end subroutine

!> Create all the coefficient matrices, in parallel
module subroutine coefficients_ifc(ih, tp, map, sim, uc, ss, mw, mem, verbosity)
    !> Coefficient matrices
    type(lo_ifc_helper), intent(out) :: ih
    !> division helper
    type(lo_trainpred_helper), intent(in) :: tp
    !> Forcemap
    type(lo_forcemap), intent(in) :: map
    !> MD simulation
    type(lo_mdsim), intent(in) :: sim
    !> Unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> Supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk?
    integer, intent(in) :: verbosity

    logical, dimension(81) :: rel_quartet_ntheta
    real(r8) :: timer, t0, t1

    ! start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! set up some thing
    init: block
        type(lo_forceconstant_secondorder) :: fc
        integer :: i, j, nf

        nf = tp%nt*ss%na*3
        ! Make space for coefficient matrices
        if (map%have_fc_singlet) then
            allocate (ih%dx1(map%xuc%nx_fc_singlet, n_division))
            ih%dx1 = 0.0_r8
            if (nf .gt. 0) then
                allocate (ih%cm1(nf, map%xuc%nx_fc_singlet))
                ih%cm1 = 0.0_r8
            else
                allocate (ih%cm1(1, 1))
                ih%cm1 = -lo_huge
            end if
        end if
        if (map%have_fc_pair) then
            allocate (ih%dx2(map%xuc%nx_fc_pair, n_division))
            ih%dx2 = 0.0_r8
            if (nf .gt. 0) then
                allocate (ih%cm2(nf, map%xuc%nx_fc_pair))
                ih%cm2 = 0.0_r8
            else
                allocate (ih%cm2(1, 1))
                ih%cm2 = -lo_huge
            end if
        end if
        if (map%have_fc_triplet) then
            allocate (ih%dx3(map%xuc%nx_fc_triplet, n_division))
            ih%dx3 = 0.0_r8
            if (nf .gt. 0) then
                allocate (ih%cm3(nf, map%xuc%nx_fc_triplet))
                ih%cm3 = 0.0_r8
            else
                allocate (ih%cm3(1, 1))
                ih%cm3 = -lo_huge
            end if
        end if
        if (map%have_fc_quartet) then
            allocate (ih%dx4(map%xuc%nx_fc_quartet, n_division))
            ih%dx4 = 0.0_r8
            if (nf .gt. 0) then
                allocate (ih%cm4(nf, map%xuc%nx_fc_quartet))
                ih%cm4 = 0.0_r8
            else
                allocate (ih%cm4(1, 1))
                ih%cm4 = -lo_huge
            end if
        end if
        ! Make space for forces and energies
        if (nf .gt. 0) then
            allocate (ih%f0(nf))
            allocate (ih%f1(nf))
            allocate (ih%f2(nf))
            allocate (ih%f3(nf))
            allocate (ih%f4(nf))
            allocate (ih%fp(nf))
            allocate (ih%u0(nf))
            allocate (ih%epot(tp%nt))
            allocate (ih%e2(tp%nt))
            allocate (ih%e3(tp%nt))
            allocate (ih%e4(tp%nt))
            allocate (ih%epol(tp%nt))
            allocate (ih%emag(tp%nt))
            ih%f0 = 0.0_r8
            ih%f1 = 0.0_r8
            ih%f2 = 0.0_r8
            ih%f3 = 0.0_r8
            ih%f4 = 0.0_r8
            ih%fp = 0.0_r8
            ih%u0 = 0.0_r8
            ih%epot = 0.0_r8
            ih%e2 = 0.0_r8
            ih%e3 = 0.0_r8
            ih%e4 = 0.0_r8
            ih%epol = 0.0_r8
            ih%emag = 0.0_r8
        else
            ! To avoid strange crashes
            allocate (ih%f0(1))
            allocate (ih%f1(1))
            allocate (ih%f2(1))
            allocate (ih%f3(1))
            allocate (ih%f4(1))
            allocate (ih%fp(1))
            allocate (ih%u0(1))
            allocate (ih%epot(1))
            allocate (ih%e2(1))
            allocate (ih%e3(1))
            allocate (ih%e4(1))
            allocate (ih%epol(1))
            allocate (ih%emag(1))
            ih%f0 = -lo_huge
            ih%f1 = -lo_huge
            ih%f2 = -lo_huge
            ih%f3 = -lo_huge
            ih%f4 = -lo_huge
            ih%fp = -lo_huge
            ih%u0 = -lo_huge
            ih%epot = -lo_huge
            ih%e2 = -lo_huge
            ih%e3 = -lo_huge
            ih%e4 = -lo_huge
            ih%epol = -lo_huge
            ih%emag = -lo_huge
        end if

        if (map%have_fc_quartet) then
            ! Check which ntheta are possible?
            rel_quartet_ntheta = .false.
            do i = 1, map%n_fc_quartet_shell
                j = map%fc_quartet_shell(i)%nx
                if (j .gt. 0) rel_quartet_ntheta(j) = .true.
            end do
        end if

        ! If it's polar, it might be a good idea to be able to remove those forces.
        if (map%polar .eq. 1) then
            ! Get some forceconstant thingy. This is a little confusing since
            ! things can be defined in different ways.
            call map%get_secondorder_forceconstant(uc, fc, mem, verbosity=-1)

            if (fc%polar .and. map%polarcorrectiontype .eq. 3) then
                allocate (ih%polar_fc(3, 3, ss%na, ss%na))
                ih%polar_fc = 0.0_r8
                call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, ih%polar_fc, 1E-15_r8)
            else if (fc%polar .and. map%polarcorrectiontype .eq. 4) then
                ! AA: polar forces for 2D materials
                write(*,*) ('Getting the LR IFCS in 2D to subtract the forces: ')
                allocate (ih%polar_fc(3, 3, ss%na, ss%na))
                ih%polar_fc = 0.0_r8
                call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, ih%polar_fc, 1E-15_r8)
            else
                allocate (ih%polar_fc(1, 1, 1, 1))
                ih%polar_fc = -lo_huge
            end if
        else
            allocate (ih%polar_fc(1, 1, 1, 1))
            ih%polar_fc = -lo_huge
       
        end if

       

    end block init

    coeff: block
        real(r8), dimension(:, :), allocatable :: CMD1, CMD2, CMD3, CMD4, pf
        real(r8) :: epol
        integer :: j, k, l, ii, jj, nf, t, tt

        nf = sim%na*3
        ! Temporary coefficient matrices
        if (map%have_fc_singlet) then
            call mem%allocate(CMD1, [nf, map%xuc%nx_fc_singlet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (map%have_fc_pair) then
            call mem%allocate(CMD2, [nf, map%xuc%nx_fc_pair], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (map%have_fc_triplet) then
            call mem%allocate(CMD3, [nf, map%xuc%nx_fc_triplet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (map%have_fc_quartet) then
            call mem%allocate(CMD4, [nf, map%xuc%nx_fc_quartet], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        call mem%allocate(pf, [3, sim%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        pf = 0.0_r8

        ! Start creating the coefficients
        if (verbosity .gt. 0) call lo_progressbar_init()
        do tt = 1, tp%nt
            t = tp%timestep_ind(tt)
            ! singlets
            if (map%have_fc_singlet) then
                call lo_coeffmatrix_singlet(CMD1, map)
            end if
            ! pairs
            if (map%have_fc_pair) then
                call lo_coeffmatrix_pair(sim%u(:, :, t), CMD2, map)
            end if
            ! triplets
            if (map%have_fc_triplet) then
                call lo_coeffmatrix_triplet(sim%u(:, :, t), CMD3, map)
            end if
            ! quartets
            if (map%have_fc_quartet) then
                call lo_coeffmatrix_quartet(sim%u(:, :, t), CMD4, map, rel_quartet_ntheta)
            end if

            ! And perhaps the polar force thingy
            ! AA: do not force the correction to be type 3!
            if (map%polar .gt. 0  .and. map%xuc%nx_Z_singlet .gt. 0) then
                pf = 0.0_r8
                epol = 0.0_r8
                do j = 1, ss%na
                    do k = 1, ss%na
                        pf(:, j) = pf(:, j) - matmul(ih%polar_fc(:, :, j, k), sim%u(:, k, t))
                    end do
                    epol = epol - dot_product(sim%u(:, j, t), pf(:, j))*0.5_r8
                end do
            else
                pf = 0.0_r8
                epol = 0.0_r8
            end if
            
            !write(*,*) ('I want to write the polar forces out:')

            !if (tt == 1) then
            !    open(unit=10, file='pf', status='replace', action='write')
            !else
            !    open(unit=10, file='pf', status='old', position='append', action='write')
            !endif
    
            !do j = 1, ss%na
            !    write(10, '(3F10.5)') pf(:, j)
            !end do

            !close(10)




            ! Store
            ii = (tt - 1)*nf + 1
            jj = tt*nf
            if (map%have_fc_singlet) ih%cm1(ii:jj, :) = -CMD1*tp%weight(tt)
            if (map%have_fc_pair) ih%cm2(ii:jj, :) = -CMD2*tp%weight(tt)
            if (map%have_fc_triplet) ih%cm3(ii:jj, :) = -CMD3*tp%weight(tt)
            if (map%have_fc_quartet) ih%cm4(ii:jj, :) = -CMD4*tp%weight(tt)
            ! And forces (subtracting the polar if applicable) + displacements
            do j = 1, sim%na
            do k = 1, 3
                l = (tt - 1)*nf + (j - 1)*3 + k
                ih%f0(l) = (sim%f(k, j, t) - pf(k, j))*tp%weight(tt)
                ih%fp(l) = pf(k, j)*tp%weight(tt)
                ih%u0(l) = sim%u(k, j, t)*tp%weight(tt)
            end do
            end do
            ! And energies

            ih%epot(tt) = sim%stat%potential_energy(t)*tp%weight(tt)
            ih%epol(tt) = epol*tp%weight(tt)
            if (verbosity .gt. 0 .and. tt .lt. tp%nt) then
                call lo_progressbar('... creating fc coefficients', tt, tp%nt, walltime() - timer)
            end if
        end do


        if (map%have_fc_singlet) call mem%deallocate(CMD1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (map%have_fc_pair) call mem%deallocate(CMD2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (map%have_fc_triplet) call mem%deallocate(CMD3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (map%have_fc_quartet) call mem%deallocate(CMD4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(pf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! make sure we are synced, for timings and things
        call mw%barrier()

        if (verbosity .gt. 0) then
            call lo_progressbar('... creating fc coefficients', tp%nt, tp%nt, walltime() - timer)
        end if
    end block coeff
end subroutine

!> Straight merge of coefficient matrices, works when not too large.
module subroutine mergematrices(sc, A, B, ne, nx, mw)
    !> container for scaling information
    class(lo_scaling_and_product), intent(out) :: sc
    !> coefficient matrix A, mpi-distributed
    real(r8), dimension(:, :), intent(in) :: A
    !> response vector B, mpi-distributed
    real(r8), dimension(:), intent(in) :: B
    !> number of rows on this rank
    integer, intent(in) :: ne
    !> number of variables
    integer, intent(in) :: nx
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    integer :: n

    call mw%size_and_offset(ne, n)

    allocate (sc%ATA(n, nx))
    allocate (sc%ATB(n))
    sc%ATA = 0.0_r8
    sc%ATB = 0.0_r8
    if (ne .gt. 0) then
        sc%ATA(mw%ctr_offset_per_rank(mw%r + 1) + 1:mw%ctr_offset_per_rank(mw%r + 1) + ne, :) = A
        sc%ATB(mw%ctr_offset_per_rank(mw%r + 1) + 1:mw%ctr_offset_per_rank(mw%r + 1) + ne) = B
    end if
    call mw%allreduce('sum', sc%ATA)
    call mw%allreduce('sum', sc%ATB)
end subroutine

!> Scale and shift coefficient matrices
module subroutine scale_and_shift(sc, A, B, ne, nx, mw, noscale, usereference, n_atom)
    !> container for scaling information
    class(lo_scaling_and_product), intent(out) :: sc
    !> coefficient matrix A, mpi-distributed
    real(r8), dimension(:, :), intent(in) :: A
    !> response vector B, mpi-distributed
    real(r8), dimension(:), intent(in) :: B
    !> number of rows on this rank
    integer, intent(in) :: ne
    !> number of variables
    integer, intent(in) :: nx
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Skip the scaling part
    logical, intent(in), optional :: noscale
    !> Skip the usage of reference positions
    logical, intent(in), optional :: usereference
    !> If I am to skip the use of reference positions, I need the number of atoms
    integer, intent(in), optional :: n_atom

    real(r8), dimension(:, :), allocatable :: wA
    real(r8), dimension(:), allocatable :: wB, d1, d2
    real(r8) :: f0, wtfact, sqrtwtfact
    integer :: i, j, l, nt, n
    logical :: scalethings
    logical :: skipreference

    if (present(noscale)) then
        scalethings = .not. noscale
    else
        scalethings = .true.
    end if

    if (present(usereference)) then
        skipreference = .not. usereference
    else
        skipreference = .false.
    end if

    ! small sanity check
    if (skipreference) then
        if (.not. present(n_atom)) then
            call lo_stop_gracefully(['Need number of atoms'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
        if (scalethings) then
            call lo_stop_gracefully(['Can not scale and skip reference at the same time'], &
                                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
    end if

    ! Make some space
    allocate (sc%xs(nx))
    allocate (sc%xm(nx))
    allocate (sc%ATA(nx, nx))
    allocate (sc%ATB(nx))
    sc%xs = 0.0_r8
    sc%xm = 0.0_r8
    sc%ys = 0.0_r8
    sc%ym = 0.0_r8
    sc%ATA = 0.0_r8
    sc%ATB = 0.0_r8

    if (scalethings) then
        if (ne .gt. 0) then
            allocate (wA(ne, nx))
            allocate (wB(ne))
            allocate (d1(ne))
            wA = A
            wB = B
            d1 = 0.0_r8
        end if
        allocate (d2(nx))
        d2 = 0.0_r8

        i = ne
        call mw%allreduce('sum', i, j)
        wtfact = 1.0_r8
        sqrtwtfact = 1.0_r8
        ! Start scaling, get the center
        if (ne .gt. 0) then
            d1 = 1.0_r8*sqrtwtfact
            call lo_gemv(wA, d1, sc%xm, trans='T')
        end if
        call mw%allreduce('sum', sc%xm)

        ! Shift them to the mean, I think
        d2 = 0.0_r8
        if (ne .gt. 0) then
            do i = 1, nx
                wA(:, i) = sqrtwtfact*(wA(:, i) - sc%xm(i))
                d2(i) = dot_product(wA(:, i), wA(:, i))
            end do
        end if
        call mw%allreduce('sum', d2)
        sc%xs = sqrt(d2)
        ! Then scale
        if (ne .gt. 0) then
            do i = 1, nx
                wA(:, i) = wA(:, i)/sc%xs(i)
            end do
        end if

        ! Now scale the other things
        if (ne .gt. 0) then
            sc%ym = sum(wB)*wtfact
        else
            sc%ym = 0.0_r8
        end if
        call mw%allreduce('sum', sc%ym)
        if (ne .gt. 0) then
            wB = sqrtwtfact*(wB - sc%ym)
            f0 = dot_product(wB, wB)
        else
            f0 = 0.0_r8
        end if
        call mw%allreduce('sum', f0)
        sc%ys = sqrt(f0)
        if (ne .gt. 0) then
            wB = wB/sc%ys
        end if
        ! And now the large matrix product
        if (ne .gt. 0) then
            ! And do the large matrix multiplications
            call lo_gemm(wA, wA, sc%ATA, transa='T')
            call lo_gemv(wA, wB, sc%ATB, trans='T')
        end if
        if (ne .gt. 0) then
            deallocate (wA)
            deallocate (wB)
            deallocate (d1)
        end if
        deallocate (d2)
    else
        if (skipreference .and. ne .gt. 0) then  !& fprettify  messes up the `ne`
            ! Ujuj, do fancy reference skipping thing?
            nt = ne/n_atom/3  ! number of timesteps
            n = nt*(nt - 1)/2  ! new number of configurations?
            allocate (wA(n*n_atom*3, nx))
            allocate (wB(n*n_atom*3))
            wA = 0.0_r8
            wB = 0.0_r8
            l = 0
            do i = 1, nt
            do j = i + 1, nt
                l = l + 1
                wA((l - 1)*n_atom*3 + 1:l*n_atom*3, :) &
                    = A((i - 1)*n_atom*3 + 1:i*n_atom*3, :) - A((j - 1)*n_atom*3 + 1:j*n_atom*3, :)
                wB((l - 1)*n_atom*3 + 1:l*n_atom*3) &
                    = B((i - 1)*n_atom*3 + 1:i*n_atom*3) - B((j - 1)*n_atom*3 + 1:j*n_atom*3)
            end do
            end do
            call lo_gemm(wA, wA, sc%ATA, transa='T')
            call lo_gemv(wA, wB, sc%ATB, trans='T')
        else
            ! This is the normal, boring case. Just multiply things.
            if (ne .gt. 0) then
                call lo_gemm(A, A, sc%ATA, transa='T')
                call lo_gemv(A, B, sc%ATB, trans='T')
            end if
        end if
    end if

    ! And accumulate over ranks
    call mw%allreduce('sum', sc%ATA)
    call mw%allreduce('sum', sc%ATB)
end subroutine

!> destroy
module subroutine destroy_scaleshift(sc)
    !> container for scaling information
    class(lo_scaling_and_product), intent(inout) :: sc

    if (allocated(sc%xs)) deallocate (sc%xs)
    if (allocated(sc%xm)) deallocate (sc%xm)
    if (allocated(sc%ATA)) deallocate (sc%ATA)
    if (allocated(sc%ATB)) deallocate (sc%ATB)
    sc%ys = -lo_huge
    sc%ym = -lo_huge
end subroutine

!> Statistical R^2, but distriuted
module function distributed_rsquare(predictions, values, n, mw) result(rsquare)
    !> predictions (per rank)
    real(r8), dimension(:), intent(in) :: predictions
    !> correct values (per rank)
    real(r8), dimension(:), intent(in) :: values
    !> how many measurements on this rank (can be zero)
    integer, intent(in) :: n
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> rsquare
    real(r8) :: rsquare

    real(r8) :: mean
    real(r8), dimension(2) :: v

    if (n .gt. 0) then
        v(1) = sum(values)
        v(2) = size(values)
    else
        v = 0.0_r8
    end if
    call mw%allreduce('sum', v)
    mean = v(1)/v(2)
    if (n .gt. 0) then
        v(1) = sum((values - predictions)**2)
        v(2) = sum((values - mean)**2)
    else
        v = 0.0_r8
    end if
    call mw%allreduce('sum', v)
    if (v(2) .gt. lo_sqtol) then
        rsquare = 1.0_r8 - v(1)/v(2)
    else
        rsquare = -1.0_r8
    end if
end function

!> Standard deviation, but distriuted
module function distributed_stddev(values, n, mw) result(dev)
    !> values (per rank)
    real(r8), dimension(:), intent(in) :: values
    !> how many measurements on this rank (can be zero)
    integer, intent(in) :: n
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> standard deviation
    real(r8) :: dev

    real(r8) :: mean
    real(r8), dimension(2) :: v

    if (n .gt. 0) then
        v(1) = sum(values)
        v(2) = size(values)
    else
        v = 0.0_r8
    end if
    call mw%allreduce('sum', v)
    mean = v(1)/v(2)
    if (n .gt. 0) then
        v(1) = sum((values - mean)**2)
        v(2) = size(values)
    else
        v = 0.0_r8
    end if
    call mw%allreduce('sum', v)
    if (v(2) .gt. lo_sqtol) then
        dev = sqrt(v(1)/v(2))
    else
        dev = -1.0_r8
    end if
end function

!> Free energy perturbation, but distributed
module function distributed_freepert(e0, e1, temperature, n, mw) result(deltaF)
    !> energies of H0
    real(r8), dimension(:), intent(in) :: e0
    !> energies of H1
    real(r8), dimension(:), intent(in) :: e1
    !> simulation temperature
    real(r8), intent(in) :: temperature
    !> how many measurements on this rank (can be zero)
    integer, intent(in) :: n
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> free energy difference
    real(r8) :: deltaF

    real(r8) :: mean0, mean1
    real(r8), dimension(4) :: u
    real(r8), dimension(2) :: v
    ! First we say the model must have the same mean as
    ! the real values. What it means in practice is that
    ! I subtract the mean from both.
    if (n .gt. 0) then
        u(1) = sum(e0)
        u(2) = sum(e1)
        u(3) = size(e0)
        u(4) = size(e1)
    else
        u = 0.0_r8
    end if
    call mw%allreduce('sum', u)
    mean0 = u(1)/u(3)
    mean1 = u(2)/u(4)
    ! I don't know how much cheating this is, but my gut feeling says it's ok.
    ! Anyway, now do the magical explog thing
    if (n .gt. 0) then
        v(1) = sum(exp(-((e1 - mean1) - (e0 - mean0))/temperature/lo_kb_Hartree))
        v(2) = size(e0)
    else
        v = 0.0_r8
    end if
    call mw%allreduce('sum', v)
    deltaF = -lo_kb_Hartree*temperature*log(v(1)/v(2))
end function

end submodule
