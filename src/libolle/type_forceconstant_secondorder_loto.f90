submodule(type_forceconstant_secondorder) type_forceconstant_secondorder_loto
!!
!! Deals with longrange electrostatics. Most of that is moved to it's own thing,
!! here are just things calling those routines.
!!
use gottochblandat, only: lo_unflatten_2tensor
implicit none
contains

!> ensure the Hermiticity of the Born effective charges
module subroutine set_ewald_and_enforce_borncharge_hermiticity(fc, p, mem, verbosity, fixlambda)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(inout) :: fc
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> fix lambda parameter
    real(r8), intent(in), optional :: fixlambda

    ! Tolerance for the Ewald summation. This is pretty safe.
    real(r8), parameter :: ewaldtol = 1E-20_r8

    ! First a sanity test to make sure noone uses the deprecated things.
    select case (fc%loto%correctiontype)
    case (0)
        ! fine
    case (1:2)
        ! not fine
        call lo_stop_gracefully(['This polar correction type is deprecated'], lo_exitcode_param, __FILE__, __LINE__)
    case (3)
        ! also fine
    case default
        ! This should never happen
        call lo_stop_gracefully(['Undefined polar correction type'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    setewald: block
        ! set the Ewald parameters
        if (present(fixlambda)) then
            call fc%ew%set(p, fc%loto%eps, 3, ewaldtol, verbosity, forced_lambda=fixlambda)
        else
            call fc%ew%set(p, fc%loto%eps, 2, ewaldtol, verbosity)
        end if

        ! If we are to optimize, then do that.
        !call fc%ew%force_borncharges_Hermitian(fc%loto%x_Z, fc%loto%coeff_Z, fc%loto%eps, p, verbosity)
    end block setewald

    ! Now set the charges and so on.
    setchg: block
        real(r8), dimension(:), allocatable :: Zflat
        complex(r8), dimension(:, :, :, :), allocatable :: Dc
        real(r8), dimension(3, 3) :: m
        integer :: a1, a2

        ! Convert from flat to sensible charges
        call mem%allocate(Zflat, 9*fc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        Zflat = 0.0_r8
        call dgemv('N', 9*p%na, fc%loto%nx_Z, 1.0_r8, fc%loto%coeff_Z, 9*p%na, fc%loto%x_z, 1, 0.0_r8, Zflat, 1)
        m = 0.0_r8
        do a1 = 1, p%na
            fc%loto%born_effective_charges(:, :, a1) = lo_unflatten_2tensor(Zflat((a1 - 1)*9 + 1:a1*9))
            m = m + fc%loto%born_effective_charges(:, :, a1)
        end do
        if (sum(abs(m)) .gt. lo_sqtol) then
            call lo_stop_gracefully(['Born charges do not add up to zero'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
        call mem%deallocate(Zflat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Get the dynamical matrix at gamma, and subsequently use this to fix the onsite corrections:
        call mem%allocate(Dc, [3, 3, p%na, p%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        Dc = 0.0_r8
        select case (fc%loto%correctiontype)
        case (3)
            !fc%loto%born_onsite_correction=0.0_r8
            call fc%ew%longrange_dynamical_matrix(p, [0.0_r8, 0.0_r8, 0.0_r8], fc%loto%born_effective_charges, fc%loto%born_onsite_correction, fc%loto%eps, Dc, reconly=.true.)
        case default
            call lo_stop_gracefully(['Undefined polar correction type'], lo_exitcode_param, __FILE__, __LINE__)
        end select
        ! Slight sanity check
        if (sum(abs(aimag(Dc))) .gt. lo_sqtol*fc%na**2) then
            call lo_stop_gracefully(['Dynamical matrix at gamma not real.'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if

        ! Now I can determine the onsite correction term!
        do a1 = 1, fc%na
            m = 0.0_r8
            do a2 = 1, fc%na
                m = m + real(Dc(:, :, a1, a2), r8)
            end do
            fc%loto%born_onsite_correction(:, :, a1) = -m
        end do
        fc%loto%born_onsite_correction = lo_chop(fc%loto%born_onsite_correction, 1E-13_r8)
        call mem%deallocate(Dc, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block setchg
end subroutine

!> longrange dipole-dipole dynamical matrix
module subroutine longrange_dynamical_matrix(fc, D, p, q, Dx, Dy, Dz)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> dynamical matrix
    complex(r8), dimension(:, :, :, :), intent(out) :: D
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> q-vector
    real(r8), dimension(3), intent(in) :: q
    !> gradient of dynamical matrix
    complex(r8), dimension(:, :, :, :), intent(out), optional :: Dx, Dy, Dz

    ! I just send this along!
    if (fc%loto%correctiontype .ne. 3) then
        call lo_stop_gracefully(['Need to use polar correction type 3.'], lo_exitcode_param, __FILE__, __LINE__)
    end if
    if (present(Dx)) then
        call fc%ew%longrange_dynamical_matrix(p, q, &
                                              fc%loto%born_effective_charges, fc%loto%born_onsite_correction, fc%loto%eps, &
                                              D, Dx, Dy, Dz, reconly=.true.)
    else
        call fc%ew%longrange_dynamical_matrix(p, q, &
                                              fc%loto%born_effective_charges, fc%loto%born_onsite_correction, fc%loto%eps, &
                                              D, reconly=.true.)
    end if
end subroutine

!> Non-analytical contribution at Gamma
module subroutine nonanalytical_dynamical_matrix(fc, p, qdir, D)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> direction of approach to Gamma
    real(r8), dimension(3), intent(in) :: qdir
    !> non-analytical dynamical matrix
    real(r8), dimension(:, :, :, :), intent(out) :: D

    real(r8), dimension(3) :: q
    real(r8) :: f0
    integer :: a1, a2, i, j

    q = qdir/norm2(qdir)
    f0 = 0.0_r8
    do j = 1, 3
    do i = 1, 3
        f0 = f0 + q(i)*fc%loto%eps(i, j)*q(j)
    end do
    end do

    f0 = (1.0_r8/f0)*4.0_r8*lo_pi/p%volume
    D = 0.0_r8
    do a1 = 1, p%na
    do a2 = 1, p%na
        do i = 1, 3
        do j = 1, 3
            D(i, j, a1, a2) = dot_product(q, fc%loto%born_effective_charges(i, :, a1))*dot_product(q, fc%loto%born_effective_charges(j, :, a2))*f0
        end do
        end do
    end do
    end do
end subroutine

!> A specialized copy of the generalized routine, really fast for supercells. Only works at Gamma though.
module subroutine supercell_longrange_dynamical_matrix_at_gamma(fc, ss, dynmat, thres)
    !> forceconstant with LO-TO things, for the supercell
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> forceconstant
    real(r8), dimension(:, :, :, :), intent(out) :: dynmat
    !> ewald tolerance
    real(r8), intent(in) :: thres

    ! Some cheap sanity tests
    if (ss%info%supercell .eqv. .false.) then
        call lo_stop_gracefully(['To get the supercell longrange dynamical matrix the cell needs to be a supercell'], lo_exitcode_param, __FILE__, __LINE__)
    end if
    if (maxval(ss%info%index_in_unitcell) .ne. fc%na) then
        call lo_stop_gracefully(['Mismatch between forceconstant and cell.'], lo_exitcode_param, __FILE__, __LINE__)
    end if
    if (fc%loto%correctiontype .ne. 3) then
        call lo_stop_gracefully(['Need to use polar correction type 3.'], lo_exitcode_param, __FILE__, __LINE__)
    end if

    ! Again, just send it along for now. Will fix at some point, but I'm too lazy.
    call fc%ew%supercell_longrange_forceconstant(fc%loto%born_effective_charges, fc%loto%eps, ss, dynmat, thres)
end subroutine

end submodule
