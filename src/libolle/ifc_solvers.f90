module ifc_solvers
!!
!! Routines that take various kinds of data and transforms it to
!! interaction tensor, in a large variety of ways.
!!
use konstanter, only: r8, lo_iou, lo_huge, lo_hugeint, lo_tiny, lo_status, lo_exitcode_param, &
                      lo_forceconstant_2nd_HartreeBohr_to_eVA, lo_forceconstant_3rd_HartreeBohr_to_eVA, &
                      lo_forceconstant_4th_HartreeBohr_to_eVA, lo_sqtol, lo_twopi, lo_exitcode_blaslapack, &
                      lo_frequency_THz_to_Hartree, lo_frequency_Hartree_to_THz, lo_Hartree_to_eV, &
                      lo_kb_Hartree, lo_temperaturetol
use gottochblandat, only: walltime, tochar, lo_stddev, lo_mean, lo_chop, lo_sqnorm, lo_linear_least_squares, &
                          lo_unflatten_2tensor, lo_flattentensor, open_file, lo_negsqrt, lo_outerproduct, lo_invert3x3matrix, &
                          lo_real_symmetric_eigenvalues_eigenvectors, lo_planck, lo_progressbar, lo_progressbar_init
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint
use type_forcemap, only: lo_forcemap,lo_coeffmatrix_unitcell_Z_singlet
use type_mdsim, only: lo_mdsim
use type_blas_lapack_wrappers, only: lo_gemv
use lo_longrange_electrostatics, only: lo_ewald_parameters

implicit none
private
public :: lo_solve_for_borncharges
public :: lo_solve_for_irreducible_ifc
public :: lo_irreducible_forceconstant_from_supercell_dynmat
public :: lo_supercell_dynmat_from_irreducible_forceconstant
public :: lo_irreducible_forceconstant_from_qmesh_dynmat
public :: lo_solve_for_irreducible_ifc_fastugly

! Fixed number of divisions for prediction and training. See no need to change, but maybe.
integer, parameter :: n_division = 3

!> scaling information for an AX=B system. Copied from glmnet, have to think what each thing means.
type lo_scaling_and_product
    real(r8) :: ys = -lo_huge
    real(r8) :: ym = -lo_huge
    real(r8), dimension(:), allocatable :: xs
    real(r8), dimension(:), allocatable :: xm
    real(r8), dimension(:, :), allocatable :: ATA
    real(r8), dimension(:), allocatable :: ATB
contains
    procedure :: generate => scale_and_shift
    procedure :: merge => mergematrices
    procedure :: destroy => destroy_scaleshift
end type

! Helper to chop the simulation into pieces for training/prediction
type lo_trainpred_helper
    !> number of timesteps on this rank
    integer :: nt = -lo_hugeint
    !> indices to timesteps on this rank
    integer, dimension(:), allocatable :: timestep_ind
    !> weight per timestep on this rank
    real(r8), dimension(:), allocatable :: weight
    !> timestep map, so that I can handle the divisions. Refers to local indices.
    logical, dimension(:, :), allocatable :: train_ts
    logical, dimension(:, :), allocatable :: pred_ts
    !> Stepcounter per division
    integer, dimension(n_division) :: nstep_pred = -lo_hugeint
    integer, dimension(n_division) :: nstep_train = -lo_hugeint
contains
    procedure :: collect => collect_training_prediction
end type

! Helper to solve for dielectric things
type lo_dielectric_helper
    !> number of independent
    real(r8), dimension(:, :), allocatable :: dx_eps0
    real(r8), dimension(:, :), allocatable :: dx_eps1
    real(r8), dimension(:, :), allocatable :: dx_eps2
    real(r8), dimension(:, :), allocatable :: dx_Z1
    real(r8), dimension(:, :), allocatable :: dx_Z2
    real(r8), dimension(:, :), allocatable :: dx_Z3
    !> measured values
    real(r8), dimension(:), allocatable :: Z0, Z1, Z2, Z3
    real(r8), dimension(:), allocatable :: epso, eps0, eps1, eps2
    !> coefficient matrices
    real(r8), dimension(:, :), allocatable :: cm_eps0, cm_eps1, cm_eps2
    real(r8), dimension(:, :), allocatable :: cm_Z1, cm_Z2, cm_Z3
    !> ewald parameters
    type(lo_ewald_parameters) :: ew
end type

! Helper type
type :: lo_ifc_helper
    !> coefficient matrices
    real(r8), dimension(:, :), allocatable :: cm1, cm2, cm3, cm4
    !> flattened forces
    real(r8), dimension(:), allocatable :: f0, f1, f2, f3, f4, fp
    !> flattened displacements
    real(r8), dimension(:), allocatable :: u0
    !> prediction Rsquare?
    real(r8) :: rsq1 = -lo_huge, rsq2 = -lo_huge, rsq3 = -lo_huge, rsq4 = -lo_huge
    !> Solution per division
    real(r8), dimension(:, :), allocatable :: dx1, dx2, dx3, dx4
    !> Energies
    real(r8), dimension(:), allocatable :: epot, e2, e3, e4, epol, emag
    !> polar force constant
    real(r8), dimension(:, :, :, :), allocatable :: polar_fc
end type

! Interfaces to ifc_solvers_prepsolver
interface
    module subroutine coefficients_ifc(ih, tp, map, sim, uc, ss, mw, mem, verbosity)
        type(lo_ifc_helper), intent(out) :: ih
        type(lo_trainpred_helper), intent(in) :: tp
        type(lo_forcemap), intent(in) :: map
        type(lo_mdsim), intent(in) :: sim
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine coefficients_dielectric(dh, tp, map, sim, ss, mw, mem, verbosity)
        type(lo_dielectric_helper), intent(out) :: dh
        type(lo_trainpred_helper), intent(in) :: tp
        type(lo_forcemap), intent(in) :: map
        type(lo_mdsim), intent(in) :: sim
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine grab_prediction_and_response(sp, div, order, training, nrow, A, B, C, rE, rP, fP)
        class(lo_ifc_helper), intent(in) :: sp
        integer, intent(in) :: div
        integer, intent(in) :: order
        logical, intent(in) :: training
        integer, intent(out) :: nrow
        real(r8), dimension(:, :), allocatable, intent(out), optional :: A
        real(r8), dimension(:), allocatable, intent(out), optional :: B
        real(r8), dimension(:), allocatable, intent(out), optional :: C
        real(r8), dimension(:), allocatable, intent(out), optional :: rE
        real(r8), dimension(:), allocatable, intent(out), optional :: rP
        real(r8), dimension(:), allocatable, intent(out), optional :: fP
    end subroutine
    module subroutine scale_and_shift(sc, A, B, ne, nx, mw, noscale, usereference, n_atom)
        class(lo_scaling_and_product), intent(out) :: sc
        real(r8), dimension(:, :), intent(in) :: A
        real(r8), dimension(:), intent(in) :: B
        integer, intent(in) :: ne
        integer, intent(in) :: nx
        type(lo_mpi_helper), intent(inout) :: mw
        logical, intent(in), optional :: noscale
        logical, intent(in), optional :: usereference
        integer, intent(in), optional :: n_atom
    end subroutine
    module subroutine mergematrices(sc, A, B, ne, nx, mw)
        class(lo_scaling_and_product), intent(out) :: sc
        real(r8), dimension(:, :), intent(in) :: A
        real(r8), dimension(:), intent(in) :: B
        integer, intent(in) :: ne
        integer, intent(in) :: nx
        type(lo_mpi_helper), intent(inout) :: mw
    end subroutine
    module subroutine destroy_scaleshift(sc)
        class(lo_scaling_and_product), intent(inout) :: sc
    end subroutine
    module function distributed_rsquare(predictions, values, n, mw) result(rsquare)
        real(r8), dimension(:), intent(in) :: predictions
        real(r8), dimension(:), intent(in) :: values
        integer, intent(in) :: n
        type(lo_mpi_helper), intent(inout) :: mw
        real(r8) :: rsquare
    end function
    module function distributed_stddev(values, n, mw) result(dev)
        real(r8), dimension(:), intent(in) :: values
        integer, intent(in) :: n
        type(lo_mpi_helper), intent(inout) :: mw
        real(r8) :: dev
    end function
    module function distributed_freepert(e0, e1, temperature, n, mw) result(deltaF)
        real(r8), dimension(:), intent(in) :: e0
        real(r8), dimension(:), intent(in) :: e1
        real(r8), intent(in) :: temperature
        integer, intent(in) :: n
        type(lo_mpi_helper), intent(inout) :: mw
        real(r8) :: deltaF
    end function
end interface
! Interfaces to ifc_solvers_diel
interface
    module subroutine lo_solve_for_borncharges(map, p, Z, eps, filename, mw, mem, verbosity)
        class(lo_forcemap), intent(inout) :: map
        type(lo_crystalstructure), intent(in) :: p
        real(r8), dimension(:, :, :), intent(in), optional :: Z
        real(r8), dimension(3, 3), intent(in), optional :: eps
        character(len=*), intent(in), optional :: filename
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine lsq_solve_diel(dh, tp, map, uc, ss, mw, mem, verbosity)
        type(lo_dielectric_helper), intent(inout) :: dh
        type(lo_trainpred_helper), intent(in) :: tp
        type(lo_forcemap), intent(inout) :: map
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface
! Interfaces to ifc_solvers_putget
interface
    module subroutine lo_irreducible_forceconstant_from_supercell_dynmat( &
        map, ss, dynmat_M, lrdynmat_M, dynmat_T, lrdynmat_T, enforce_constraints, &
        fullhavemass, longrangehavemass, subset, mw, verbosity)
        type(lo_forcemap), intent(inout) :: map
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), dimension(:, :), intent(in), optional :: dynmat_M
        real(r8), dimension(:, :, :, :), intent(in), optional :: dynmat_T
        real(r8), dimension(:, :), intent(in), optional :: lrdynmat_M
        real(r8), dimension(:, :, :, :), intent(in), optional :: lrdynmat_T
        logical, intent(in), optional :: enforce_constraints
        logical, intent(in) :: fullhavemass
        logical, intent(in) :: longrangehavemass
        integer, dimension(:), intent(in), optional :: subset
        type(lo_mpi_helper), intent(inout) :: mw
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine lo_supercell_dynmat_from_irreducible_forceconstant(map, uc, ss, x, dynmat, mw, mem, dmlr, shortrange)
        type(lo_forcemap), intent(in) :: map
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), dimension(:), intent(in) :: x
        real(r8), dimension(:, :), intent(out) :: dynmat
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        real(r8), dimension(:, :, :, :), intent(in), optional :: dmlr
        logical, intent(in), optional :: shortrange
    end subroutine
    module subroutine lo_irreducible_forceconstant_from_qmesh_dynmat( &
        map, uc, qp, nq, dynmat, lrdynmat, fullhavemasses, longrangehavemasses, mw, verbosity, enforce, weights)
        type(lo_forcemap), intent(inout) :: map
        type(lo_crystalstructure), intent(in) :: uc
        class(lo_qpoint), dimension(:), intent(in) :: qp
        integer, intent(in) :: nq
        complex(r8), dimension(:, :, :), intent(in) :: dynmat
        complex(r8), dimension(:, :, :), intent(in) :: lrdynmat
        logical, intent(in) :: fullhavemasses
        logical, intent(in) :: longrangehavemasses
        type(lo_mpi_helper), intent(inout) :: mw
        integer, intent(in) :: verbosity
        logical, intent(in), optional :: enforce
        real(r8), dimension(:), intent(in), optional :: weights
    end subroutine
end interface
! Interface to ifc_solvers_lsq
interface
    module subroutine lsq_solve(map, uc, ss, ih, tp, force_stable, usereference, maxdisp, temperature, mw, mem, verbosity)
        type(lo_forcemap), intent(inout) :: map
        type(lo_crystalstructure), intent(inout) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_ifc_helper), intent(inout) :: ih
        type(lo_trainpred_helper), intent(in) :: tp
        logical, intent(in) :: force_stable
        logical, intent(in) :: usereference
        real(r8), intent(in) :: maxdisp
        real(r8), intent(in) :: temperature
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface
! ifc_solvers_trainpred
interface
    module subroutine divide_simulation(tp, sim, wts, mw)
        type(lo_trainpred_helper), intent(out) :: tp
        type(lo_mdsim), intent(in) :: sim
        real(r8), dimension(:), intent(in) :: wts
        type(lo_mpi_helper), intent(inout) :: mw
    end subroutine
    module subroutine collect_training_prediction( &
        tp, map, interaction, div, training, nrow, mem, dh, ih, coeff, values, displacements, energies)
        class(lo_trainpred_helper), intent(in) :: tp
        type(lo_forcemap), intent(in) :: map
        integer, intent(in) :: interaction
        integer, intent(in) :: div
        logical, intent(in) :: training
        integer, intent(out) :: nrow
        type(lo_mem_helper), intent(inout) :: mem
        type(lo_dielectric_helper), intent(in), optional :: dh
        type(lo_ifc_helper), intent(in), optional :: ih
        real(r8), dimension(:, :), allocatable, intent(out), optional :: coeff
        real(r8), dimension(:), allocatable, intent(out), optional :: values
        real(r8), dimension(:), allocatable, intent(out), optional :: displacements
        real(r8), dimension(:), allocatable, intent(out), optional :: energies
    end subroutine
end interface
interface
    module subroutine report_diagnostics(map, tp, ih, dh, mw, mem, verbosity)
        type(lo_forcemap), intent(in) :: map
        type(lo_trainpred_helper), intent(in) :: tp
        type(lo_ifc_helper), intent(inout) :: ih
        type(lo_dielectric_helper), intent(inout) :: dh
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
end interface
! ifc_solvers_fastugly
interface
    module subroutine lo_solve_for_irreducible_ifc_fastugly(map, uc, ss, sim, mw, mem, verbosity, fix_secondorder)
        type(lo_forcemap), intent(inout) :: map
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_mdsim), intent(in) :: sim
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
        logical, intent(in) :: fix_secondorder
    end subroutine
end interface

contains

!> Solve for the irreducible representations in a variety of ways
subroutine lo_solve_for_irreducible_ifc(map, sim, uc, ss, solver, mw, mem, verbosity, &
                                        born_effective_charges, dielectric_tensor, &
                                        lotofilename, dev1, dev2, dev3, dev4, &
                                        generating_sigma, temperature, weights, stat, &
                                        prev_max_omega, max_displacement, noreference)
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> displacement-force data
    type(lo_mdsim), intent(inout) :: sim
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss
    !> Solver
    integer, intent(in) :: solver
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> Talk a lot?
    integer, intent(in) :: verbosity
    !> Born effective charges
    real(r8), dimension(:, :, :), intent(in), optional :: born_effective_charges
    !> Dielectric tensor
    real(r8), dimension(3, 3), intent(in), optional :: dielectric_tensor
    !> Filename that holds dielectric tensor and Born charges
    character(len=*), intent(in), optional :: lotofilename
    !> Standard deviation of irreducible
    real(r8), dimension(:), intent(out), optional :: dev1, dev2, dev3, dev4
    !> In case of reweighted least squares, I need the generating thermal displacement matrices
    real(r8), dimension(:, :, :, :), intent(in), optional :: generating_sigma
    !> For reweighted least squares, I also need the temperature
    real(r8), intent(in), optional :: temperature
    !> weight per configuration
    real(r8), dimension(:), intent(inout), optional :: weights
    !> Return some status things from the solver
    real(r8), dimension(6), intent(out), optional :: stat
    !> Provide the previous max omega
    real(r8), intent(in), optional :: prev_max_omega
    !> Largest allowed mean square displacement
    real(r8), intent(in), optional :: max_displacement
    !> Do a reference free fit for the second order
    logical, intent(in), optional :: noreference

    type(lo_crystalstructure) :: wuc, wss
    type(lo_trainpred_helper) :: tp
    type(lo_ifc_helper) :: ih
    type(lo_dielectric_helper) :: dh

    real(r8), dimension(6) :: dstat
    real(r8) :: timer, maxdisp
    logical :: usereference

    if (verbosity .gt. 0) then
        timer = walltime()
        write (lo_iou, *) ''
        write (lo_iou, *) 'SOLVING FOR INTERACTION TENSORS'
    end if

    ! A few general things that need to be sorted out first.
    preinit: block
        real(r8), dimension(:), allocatable :: wts
        ! Set the status to everything bad
        dstat = -lo_hugeint
        ! Set the output to nothing
        if (map%have_fc_singlet) map%xuc%x_fc_singlet = 0.0_r8
        if (map%have_fc_pair) map%xuc%x_fc_pair = 0.0_r8
        if (map%have_fc_triplet) map%xuc%x_fc_triplet = 0.0_r8
        if (map%have_fc_quartet) map%xuc%x_fc_quartet = 0.0_r8
        if (present(dev1)) dev1 = 0.0_r8
        if (present(dev2)) dev2 = 0.0_r8
        if (present(dev3)) dev3 = 0.0_r8
        if (present(dev4)) dev4 = 0.0_r8

        ! Sort out weights
        call mem%allocate(wts, sim%nt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        wts = 0.0_r8
        if (present(weights)) then
            wts = weights
        else
            wts = 1.0_r8
        end if

        ! Sort out the usage of reference positions for the second order
        if (present(noreference)) then
            usereference = .not. noreference
        else
            usereference = .true.
        end if

        ! Divide the simulation across ranks
        call divide_simulation(tp, sim, wts, mw)
        call mem%deallocate(wts, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        if (verbosity .gt. 0) then
            write (lo_iou, *) '... distributed data'
        end if
    end block preinit

    ! Before we tackle the force constants, we deal with the born charges and dielectric tensors to all orders.
    call mem%tick()
    dielectric: block
        ! Decide on what to do
        select case (map%polar)
        case (1)
            if (verbosity .gt. 0) then
                write (lo_iou, *) '... detected polar interactions'
            end if
            ! Some sanity checks, then just a simple solution.
            if (present(born_effective_charges) .and. present(dielectric_tensor)) then
                call lo_solve_for_borncharges(map, uc, Z=born_effective_charges, &
                                              eps=dielectric_tensor, mw=mw, mem=mem, verbosity=verbosity)
            elseif (present(lotofilename)) then
                call lo_solve_for_borncharges(map, uc, filename=trim(lotofilename), mw=mw, mem=mem, verbosity=verbosity)
            else
                call lo_stop_gracefully(['You have to provide Born effective charges and dielectric tensors somehow.'], &
                                        lo_exitcode_param, __FILE__, __LINE__)
            end if
        case (2)
            if (verbosity .gt. 0) then
                write (lo_iou, *) '... detected electric response as a function of position'
            end if
            ! Make sure the pertinent information actually exist.
            if (.not. allocated(sim%eps)) then
                call lo_stop_gracefully(['Simulation does not contain dielectric tensor for every step.'], &
                                        lo_exitcode_param, __FILE__, __LINE__)
            end if
            if (.not. allocated(sim%Z)) then
                call lo_stop_gracefully(['Simulation does not contain Born charges for every step.'], &
                                        lo_exitcode_param, __FILE__, __LINE__)
            end if
            ! pre-build the coefficient matrices
            call coefficients_dielectric(dh, tp, map, sim, ss, mw, mem, verbosity)
            ! solve
            call lsq_solve_diel(dh, tp, map, uc, ss, mw, mem, verbosity)
        case default
            ! nothing to worry about
            if (verbosity .gt. 0) then
                write (lo_iou, *) '... no polar interactions to worry about'
            end if
        end select
    end block dielectric
    call mem%tock(__FILE__, __LINE__, mw%comm)

    ! As soon as we know the Born charges we can construct the IFC constraints
    buildconstraints: block
        type(lo_ewald_parameters) :: ew
        complex(r8), dimension(:,:,:,:), allocatable :: D0
        real(r8), dimension(:,:,:,:), allocatable :: rottensor
        real(r8), dimension(:,:,:), allocatable :: wZ,wD
        real(r8), dimension(:,:), allocatable :: coeff_Z
        real(r8), dimension(:), allocatable :: v0
        real(r8), dimension(:), allocatable :: hermitian_rhs, huang_rhs, rotational_rhs
        real(r8), dimension(3,3,3,3) :: bracket
        real(r8), dimension(3,3) :: eps,m0,m1
        integer :: i,j,k,l,ii
        integer :: a1,a2

        allocate(hermitian_rhs(9*uc%na))
        allocate(rotational_rhs(27*uc%na))
        allocate(huang_rhs(81))
        hermitian_rhs=0.0_r8
        rotational_rhs=0.0_r8
        huang_rhs=0.0_r8

        if ( map%polar .gt. 0 ) then
            ! First thing we need is the polar dynamical matrix to know
            ! how we should build the Hermiticity constraints.
            eps=lo_unflatten_2tensor(matmul(map%eps_global_shell%coeff, map%xuc%x_eps_global))

            allocate(wZ(3,3,uc%na))
            allocate(wD(3,3,uc%na))
            allocate(v0(9*uc%na))
            allocate(coeff_Z(uc%na*9,map%xuc%nx_Z_singlet))
            allocate(D0(3,3,uc%na,uc%na))
            wZ=0.0_r8
            wD=0.0_r8
            v0=0.0_r8
            coeff_Z=0.0_r8
            D0=0.0_r8

            call lo_coeffmatrix_unitcell_Z_singlet(map, coeff_Z)
            v0=matmul(coeff_Z,map%xuc%x_Z_singlet)
            do a1 = 1, uc%na
                wZ(:, :, a1) = lo_unflatten_2tensor(v0((a1 - 1)*9 + 1:a1*9))
            end do

            call ew%set(uc, eps, 2, 1E-20_r8, verbosity=-1)
            call ew%longrange_dynamical_matrix(uc, [0.0_r8, 0.0_r8, 0.0_r8], wZ, wD, eps, D0, reconly=.true., chgmult=.true.)

            do a1=1,uc%na
                m0=0.0_r8
                m1=0.0_r8
                do a2=1,uc%na
                    m0=m0+real(d0(:,:,a1,a2),r8)
                    m1=m1+real(d0(:,:,a2,a1),r8)
                enddo
                m0=m0-m1
                hermitian_rhs( (a1-1)*9+1:a1*9 ) = -lo_flattentensor(m0)
            enddo
            hermitian_rhs = lo_chop(hermitian_rhs,1E-12_r8)

            ! Next up we do the huang + rotational invariances
            allocate(rottensor(3,3,3,uc%na))
            rottensor=0.0_r8
            call ew%longrange_elastic_constant_bracket( uc, wZ, eps, bracket, rottensor, reconly=.true., mw=mw, verbosity=verbosity)

            huang_rhs=0.0_r8
            ii=0
            do i=1,3
            do j=1,3
            do k=1,3
            do l=1,3
                ii=ii+1
                huang_rhs(ii) = bracket(k,l,i,j) - bracket(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            huang_rhs = lo_chop(huang_rhs,1E-12_r8)

            ii=0
            do a1=1,uc%na
            do i=1,3
            do j=1,3
            do k=1,3
                ii=ii+1
                rotational_rhs(ii) = rottensor(i,j,k,a1)-rottensor(i,k,j,a1)
            enddo
            enddo
            enddo
            enddo
            rotational_rhs = lo_chop(rotational_rhs,1E-12_r8)

        else
            ! We just have zeros
            huang_rhs=0.0_r8
            hermitian_rhs=0.0_r8
            rotational_rhs=0.0_r8
        endif

        call map%forceconstant_constraints(uc,rotational=.true.,huanginvariances=.true.,hermitian=.true.,elastic=.false.,hermitian_rhs=hermitian_rhs,huang_rhs=huang_rhs,rotational_rhs=rotational_rhs,verbosity=verbosity)
        !call map%forceconstant_constraints(uc,rotational=.false.,huanginvariances=.true.,hermitian=.true.,hermitian_rhs=hermitian_rhs,huang_rhs=huang_rhs,verbosity=verbosity)
    end block buildconstraints

    ! Get some things that are always needed
    !call mem%tick()
    ifc: block
        real(r8) :: temp

        ! Is there a minimum frequency specified?
        if (present(max_displacement)) then
            maxdisp = max_displacement
        else
            maxdisp = -lo_huge ! Or some other magic number
        end if

        ! Now, should I relax? Better to deal with it right away, I think.
        if (map%have_fc_singlet) then
            ! Well I wouldn't have used the first order if I did not need to relax?
            wuc = uc
            wss = ss
            !call relax_positions_and_maybe_cell(map,uc,ss,sim,mw,verbosity,wuc,wss,wts)
        else
            ! Make a copy of the existing structures instead? Think so, to be on the safe side.
            wuc = uc
            wss = ss
        end if

        ! Get the coefficient matrices and things like that
        call coefficients_ifc(ih, tp, map, sim, wuc, wss, mw, mem, verbosity + 2)

        ! Get the actual solution
        select case (solver)
        case (1)
            ! Normal least squares solver. Generally the best.
            call lsq_solve(map, uc, ss, ih, tp, force_stable=.false., &
                           usereference=usereference, maxdisp=maxdisp, temperature=-1.0_r8, mw=mw, mem=mem, &
                           verbosity=verbosity)
        case (2)
            ! Least squares solver, but always stable!
            ! Check that a max displacement is provided
            if (maxdisp .lt. 0.0_r8) then
                call lo_stop_gracefully(['For a stable solution you have to provide a largest allowed', &
                                         'displacement, in units of the nearest neighbour distance.  '], &
                                        lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end if
            if (present(temperature)) then
                if (temperature .gt. -lo_tiny) then
                    temp = temperature
                else
                    call lo_stop_gracefully(['For a stable solution you have to provide a temperature.'], &
                                            lo_exitcode_param, __FILE__, __LINE__, mw%comm)
                end if
            else
                call lo_stop_gracefully(['For a stable solution you have to provide a temperature.'], &
                                        lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end if
            call lsq_solve(map, uc, ss, ih, tp, force_stable=.true., &
                           usereference=usereference, maxdisp=maxdisp, temperature=temp, mw=mw, mem=mem, &
                           verbosity=verbosity + 1)
            ! case(3)
            !     ! Iteratively reweighted least squares. A little fiddly. Needs much
            !     ! more input, so make sure everything is here.
            !     lo_status=0
            !     if ( .not. present(generating_sigma) ) lo_status=lo_status+1
            !     if ( .not. present(prev_max_omega)   ) lo_status=lo_status+1
            !     if ( .not. present(stat)             ) lo_status=lo_status+1
            !     if ( .not. present(temperature)      ) lo_status=lo_status+1
            !     ! Iteratively reweighted stable least squares. This solver puts slightly rougher
            !     ! requirements on the input, so check that right away.
            !     ! if ( map%thirdorder ) then
            !     !     call lo_stop_gracefully(['Cannot use reweighted least squares with third order forceconstants'],&
            !     !                             lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            !     ! endif
            !     ! if ( map%have_fc_quartet ) then
            !     !     call lo_stop_gracefully(['Cannot use reweighted least squares with fourth order forceconstants'],&
            !     !                             lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            !     ! endif
            !     ! if ( .not.present(generating_sigma) ) then
            !     !     call lo_stop_gracefully(['Need to provide the generating thermal displacement matrices for IRWLSQ'],&
            !     !                             lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            !     ! endif
            !     ! if ( .not.present(temperature) ) then
            !     !     call lo_stop_gracefully(['Need to provide temperature for IRWLSQ'],&
            !     !                             lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            !     ! endif
            !     ! if ( .not.present(weights) ) then
            !     !     call lo_stop_gracefully(['Need to provide weights for IRWLSQ'],&
            !     !                             lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            !     ! endif
            !     ! call reweight_solve(map,wuc,wss,sp,mw,verbosity+2,&
            !     !                     temperature=temperature,stat=dstat,prev_max_omega=prev_max_omega)
            !     ! if ( present(stat) ) stat=dstat
        case default
            call lo_stop_gracefully(['Unknown solver type:'//tochar(solver)], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end select

        ! And we are done!
        if (verbosity .gt. 0) then
            write (lo_iou, *) '... finished solving for forceconstants'
        end if
    end block ifc
    !call mem%tock(__FILE__,__LINE__,mw%comm)

    ! Report some statistics
    finalize: block
        integer :: i

        ! check if we have to do something
        i = verbosity
        call mw%allreduce('max', i)
        if (i .gt. 0) then
            call report_diagnostics(map, tp, ih, dh, mw, mem, verbosity)
        end if
    end block finalize
    !     real(r8), dimension(n_division) :: prsq2,prsq3,prsq4
    !     real(r8), dimension(n_division) :: ersq2,ersq3,ersq4,ersq0
    !     real(r8), dimension(n_division) :: fpert2,fpert3,fpert4
    !     ! Some temporary workspace
    !     real(r8), dimension(:,:), allocatable :: wA
    !     real(r8), dimension(:), allocatable :: wf0,wu0,wB,wC,wD,wE,wepot,wepol,wfp
    !     !real(r8) :: f0,f1,f2,e0,e1,e2,beta,tomev
    !     real(r8) :: edev2,edev3,edev4,fdiff2,fdiff3,fdiff4,f0
    !     integer :: div,nf,nt,i,j,ii,jj
    !
    !     prsq2=0.0_r8
    !     prsq3=0.0_r8
    !     prsq4=0.0_r8
    !     ersq0=0.0_r8
    !     ersq2=0.0_r8
    !     ersq3=0.0_r8
    !     ersq4=0.0_r8
    !     fpert2=0.0_r8
    !     fpert3=0.0_r8
    !     fpert4=0.0_r8
    !     ! ! First get all the statistics from the different cross-validation sets
    !     ! do div=1,n_division
    !     !     ! First fetch all the forces and energies
    !     !     call grab_prediction_and_response(sp,div,0,.true.,nf,B=wf0,C=wu0,rE=wepot,rP=wepol,fP=wfp)
    !     !     ! Some work arrays
    !     !     if ( nf .gt. 0 ) then
    !     !         wf0=wf0+wfp
    !     !         nt=nf/sim%na/3
    !     !         allocate(wB(nf))
    !     !         allocate(wC(nf))
    !     !         allocate(wD(nt))
    !     !         allocate(wE(nt))
    !     !         wB=0.0_r8
    !     !         wC=0.0_r8
    !     !         wD=0.0_r8
    !     !         wE=0.0_r8
    !     !     else
    !     !         allocate(wB(1))
    !     !         allocate(wC(1))
    !     !         allocate(wD(1))
    !     !         allocate(wE(1))
    !     !         nt=0
    !     !         wB=-lo_huge
    !     !         wC=-lo_huge
    !     !         wD=-lo_huge
    !     !         wE=-lo_huge
    !     !     endif
    !     !     ! Baseline standard deviation
    !     !     ersq0(div)=distributed_stddev(wepot,nt,mw)
    !     !
    !     !     if ( map%polar .gt. 0 ) then
    !     !         ! If we have polar contributions, don't forget about them
    !     !         if ( nf .gt. 0) then
    !     !             wC=wC+wfp
    !     !             wE=wE+wepol
    !     !         endif
    !     !     endif
    !     !    if ( map%have_fc_pair ) then
    !     !         call grab_prediction_and_response(sp,div,2,.true.,nf,A=wA)
    !     !         if ( nf .gt. 0 ) then
    !     !             call lo_gemv(wA,ih%dx2(:,div),wB)
    !     !             do i=1,nt
    !     !                 ii=(i-1)*sim%na*3+1
    !     !                 jj=i*sim%na*3
    !     !                 wD(i)=-dot_product(wu0(ii:jj),wB(ii:jj))*0.5_r8
    !     !             enddo
    !     !             wC=wC+wB
    !     !             wE=wE+wD
    !     !         endif
    !     !         prsq2(div)=distributed_rsquare(wC,wf0,nf,mw)
    !     !         ersq2(div)=distributed_stddev(wepot-wE,nt,mw)
    !     !         fpert2(div)=distributed_freepert(wE,wepot,sim%temperature_thermostat,nt,mw)
    !     !         lo_deallocate(wA)
    !     !    endif
    !     !     if ( map%have_fc_triplet ) then
    !     !         call grab_prediction_and_response(sp,div,3,.true.,nf,A=wA)
    !     !         if ( nf .gt. 0 ) then
    !     !             call lo_gemv(wA,ih%dx3(:,div),wB)
    !     !             do i=1,nt
    !     !                 ii=(i-1)*sim%na*3+1
    !     !                 jj=i*sim%na*3
    !     !                 wD(i)=-dot_product(wu0(ii:jj),wB(ii:jj))/3.0_r8
    !     !             enddo
    !     !             wC=wC+wB
    !     !             wE=wE+wD
    !     !         endif
    !     !         prsq3(div)=distributed_rsquare(wC,wf0,nf,mw)
    !     !         ersq3(div)=distributed_stddev(wepot-wE,nt,mw)
    !     !         fpert3(div)=distributed_freepert(wE,wepot,sim%temperature_thermostat,nt,mw)
    !     !         lo_deallocate(wA)
    !     !     endif
    !     !     if ( map%have_fc_quartet ) then
    !     !         call grab_prediction_and_response(sp,div,4,.true.,nf,A=wA)
    !     !         if ( nf .gt. 0 ) then
    !     !             call lo_gemv(wA,ih%dx4(:,div),wB)
    !     !             do i=1,nt
    !     !                 ii=(i-1)*sim%na*3+1
    !     !                 jj=i*sim%na*3
    !     !                 wD(i)=-dot_product(wu0(ii:jj),wB(ii:jj))/4.0_r8
    !     !             enddo
    !     !             wC=wC+wB
    !     !             wE=wE+wD
    !     !         endif
    !     !         prsq4(div)=distributed_rsquare(wC,wf0,nf,mw)
    !     !         ersq4(div)=distributed_stddev(wepot-wE,nt,mw)
    !     !         fpert4(div)=distributed_freepert(wE,wepot,sim%temperature_thermostat,nt,mw)
    !     !         lo_deallocate(wA)
    !     !     endif
    !     !     ! And cleanup before the next division
    !     !     if ( allocated(wB) ) lo_deallocate(wB)
    !     !     if ( allocated(wC) ) lo_deallocate(wC)
    !     !     if ( allocated(wD) ) lo_deallocate(wD)
    !     !     if ( allocated(wE) ) lo_deallocate(wE)
    !     !     if ( allocated(wu0) ) lo_deallocate(wu0)
    !     !     if ( allocated(wf0) ) lo_deallocate(wf0)
    !     !     if ( allocated(wf0) ) lo_deallocate(wu0)
    !     !     if ( allocated(wfp) ) lo_deallocate(wfp)
    !     !     if ( allocated(wepot) ) lo_deallocate(wepot)
    !     !     if ( allocated(wepol) ) lo_deallocate(wepol)
    !     ! enddo
    !
    !     ! ! Also, I might as well use the full solution to get good estimates
    !     ! if ( sp%nt .gt. 0 ) then
    !     !     lo_allocate(wf0(ih%nf))
    !     !     lo_allocate(wu0(ih%nt))
    !     ! endif
    !     ! if ( map%have_fc_pair ) then
    !     !     if ( ih%nt .gt. 0 ) then
    !     !         call lo_gemv(ih%cm2,ih%x2,wf0)
    !     !         do i=1,ih%nt
    !     !             ii=(i-1)*sim%na*3+1
    !     !             jj=i*sim%na*3
    !     !             ih%e2(i)=-dot_product(ih%u0(ii:jj),wf0(ii:jj))*0.5_r8
    !     !         enddo
    !     !     endif
    !     !     ih%rsq2=distributed_rsquare(ih%f2,ih%f0,ih%nf,mw)
    !     !     wu0=ih%epot-ih%e2
    !     !     edev2=distributed_stddev(wu0,nt,mw)
    !     !     fdiff2=distributed_freepert(ih%e2,ih%epot,sim%temperature_thermostat,nt,mw)
    !     ! endif
    !     ! if ( map%have_fc_triplet ) then
    !     !     if ( ih%nt .gt. 0 ) then
    !     !         call lo_gemv(ih%cm3,ih%x3,wf0)
    !     !         do i=1,ih%nt
    !     !             ii=(i-1)*sim%na*3+1
    !     !             jj=i*sim%na*3
    !     !             ih%e3(i)=-dot_product(ih%u0(ii:jj),wf0(ii:jj))*0.5_r8
    !     !         enddo
    !     !     endif
    !     !     ih%rsq3=distributed_rsquare(ih%f3,ih%f0,ih%nf,mw)
    !     !     wu0=ih%epot-ih%e2-ih%e3
    !     !     edev3=distributed_stddev(wu0,nt,mw)
    !     !     wu0=ih%e2+ih%e3
    !     !     fdiff3=distributed_freepert(wu0,ih%epot,sim%temperature_thermostat,nt,mw)
    !     ! endif
    !     ! if ( map%have_fc_quartet ) then
    !     !     if ( ih%nt .gt. 0 ) then
    !     !         call lo_gemv(ih%cm4,ih%x4,wf0)
    !     !         do i=1,ih%nt
    !     !             ii=(i-1)*sim%na*3+1
    !     !             jj=i*sim%na*3
    !     !             ih%e4(i)=-dot_product(ih%u0(ii:jj),wf0(ii:jj))*0.5_r8
    !     !         enddo
    !     !     endif
    !     !
    !     !     ih%rsq4=distributed_rsquare(ih%f4,ih%f0,ih%nf,mw)
    !     !     wu0=ih%epot-ih%e2-ih%e3-ih%e4
    !     !     edev4=distributed_stddev(wu0,nt,mw)
    !     !     wu0=ih%e2+ih%e3+ih%e4
    !     !     fdiff4=distributed_freepert(wu0,ih%epot,sim%temperature_thermostat,nt,mw)
    !     ! endif
    !     ! if ( allocated(wf0) ) deallocate(wf0)
    !     ! if ( allocated(wu0) ) deallocate(wu0)
    !
    !     ! Report to stdout
    !     if ( verbosity .gt. 0 ) then
    !         write(*,*) ''
    !         write(*,*) 'Measures on forceconstants'
    !         if ( map%have_fc_pair ) then
    !             prsq2=max(prsq2,-1.0_r8)
    !             f0=100*lo_stddev(prsq2)/ih%rsq2
    !             write(*,"(1X,A,1X,F9.7,' (',F9.7,'%)')") '    R^2 secondorder:',ih%rsq2,f0
    !         endif
    !         if ( map%have_fc_triplet ) then
    !             prsq3=max(prsq3,-1.0_r8)
    !             write(*,"(1X,A,1X,F9.7,A,F9.7,A)") '    R^2 thirdorder :',lo_mean(prsq3),'(',lo_stddev(prsq3),')'
    !         endif
    !         if ( map%have_fc_quartet ) then
    !             prsq4=max(prsq4,-1.0_r8)
    !             write(*,"(1X,A,1X,F9.7,A,F9.7,A)") '    R^2 fourthorder:',lo_mean(prsq4),'(',lo_stddev(prsq4),')'
    !         endif
    !         if ( map%have_fc_pair ) then
    !             f0=lo_forceconstant_2nd_HartreeBohr_to_eVA*100
    !             write(*,"(1X,A)") '             theta(i)         stddev(theta(i))'
    !             do j=1,min(10,map%xuc%nx_fc_pair)
    !                 write(*,"(3X,I3,1X,F18.8,1X,F18.8)") j,map%xuc%x_fc_pair(j)*f0,lo_stddev(ih%dx2(j,:))*f0
    !             enddo
    !         endif
    !         if ( map%have_fc_triplet ) then
    !             f0=lo_forceconstant_3rd_HartreeBohr_to_eVA*100
    !             write(*,"(1X,A)") '             theta(i)         stddev(theta(i))'
    !             do j=1,min(10,map%xuc%nx_fc_triplet)
    !                 write(*,"(3X,I3,1X,F18.8,1X,F18.8)") j,map%xuc%x_fc_triplet(j)*f0,lo_stddev(ih%dx3(j,:))*f0
    !             enddo
    !         endif
    !         if ( map%have_fc_quartet ) then
    !             f0=lo_forceconstant_4th_HartreeBohr_to_eVA*100
    !             write(*,"(1X,A)") '             theta(i)         stddev(theta(i))'
    !             do j=1,min(10,map%xuc%nx_fc_quartet)
    !                 write(*,"(3X,I3,1X,F18.8,1X,F18.8)") j,map%xuc%x_fc_quartet(j)*f0,lo_stddev(ih%dx4(j,:))*f0
    !             enddo
    !         endif
    !         if ( map%polar .gt. 0 ) then
    !             write(*,"(1X,A)") '             theta(i)         stddev(theta(i))'
    !             do j=1,min(10,map%xuc%nx_eps_global)
    !                 write(*,*) j,map%xuc%x_eps_global(j)
    !             enddo
    !             write(*,"(1X,A)") '             theta(i)         stddev(theta(i))'
    !             do j=1,min(10,map%xuc%nx_Z_singlet)
    !                 write(*,*) j,map%xuc%x_Z_singlet(j)
    !             enddo
    !         endif
    !
    !         ! Say something about how close we are in free energy.
    !         ! Maybe later. But it's not a bad idea.
    !
    !         ! tomev=lo_Hartree_to_eV*1000/map%n_atom_ss
    !         ! if ( sim%temperature_thermostat .gt. lo_temperaturetol ) then
    !         !     beta=1.0_r8/sim%temperature_thermostat/lo_kb_Hartree
    !         ! else
    !         !     beta=0.0_r8
    !         ! endif
    !         !
    !         ! write(*,*) ''
    !         ! write(*,*) 'Measures of free energy corrections (all in meV/atom)'
    !         ! write(*,"(1X,5X,A,5X,A,5X,A)") 'stddev','deltaF1','deltaF2'
    !         ! if ( map%have_fc_pair ) then
    !         !     f0=edev2*tomev
    !         !     ersq2=ersq2/sim%na
    !         !     e0=lo_stddev(ersq2)*tomev*sim%na
    !         !     if ( sim%stochastic ) then
    !         !         f1=-beta*edev2**2*0.5_r8*tomev
    !         !         e1=-lo_stddev(beta*ersq2**2)*tomev*0.5_r8*sim%na
    !         !     else
    !         !         f1=beta*edev2**2*0.5_r8*tomev*sim%na
    !         !         e1=lo_stddev(beta*ersq2**2)*tomev*0.5_r8*sim%na
    !         !     endif
    !         !     f2=fdiff2*tomev
    !         !     e2=lo_stddev(fpert2)*tomev
    !         !     write(*,"(3X,6(1X,F12.5))") f0,e0,f1,e1,f2,e2
    !         ! endif
    !         ! if ( map%thirdorder ) then
    !         !     f0=edev3*tomev
    !         !     ersq3=ersq3/sim%na
    !         !     e0=lo_stddev(ersq3)*tomev*sim%na
    !         !     if ( sim%stochastic ) then
    !         !         f1=-beta*edev3**2*0.5_r8*tomev
    !         !         e1=-lo_stddev(beta*ersq3**2)*tomev*0.5_r8*sim%na
    !         !     else
    !         !         f1=beta*edev3**2*0.5_r8*tomev*sim%na
    !         !         e1=lo_stddev(beta*ersq3**2)*tomev*0.5_r8*sim%na
    !         !     endif
    !         !     f2=fdiff3*tomev
    !         !     e2=lo_stddev(fpert3)*tomev
    !         !     write(*,"(3X,6(1X,F12.5))") f0,e0,f1,e1,f2,e2
    !         ! endif
    !         ! if ( map%have_fc_quartet ) then
    !         !     f0=edev4*tomev
    !         !     ersq4=ersq4/sim%na
    !         !     e0=lo_stddev(ersq4)*tomev*sim%na
    !         !     if ( sim%stochastic ) then
    !         !         f1=-beta*edev4**2*0.5_r8*tomev
    !         !         e1=-lo_stddev(beta*ersq4**2)*tomev*0.5_r8*sim%na
    !         !     else
    !         !         f1=beta*edev4**2*0.5_r8*tomev*sim%na
    !         !         e1=lo_stddev(beta*ersq4**2)*tomev*0.5_r8*sim%na
    !         !     endif
    !         !     f2=fdiff4*tomev
    !         !     e2=lo_stddev(fpert4)*tomev
    !         !     write(*,"(3X,6(1X,F12.5))") f0,e0,f1,e1,f2,e2
    !         ! endif
    !     endif
    ! end block reportstat

end subroutine

end module
