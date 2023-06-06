submodule(type_forcemap) type_forcemap_coefficient_diel
!!
!! Build coefficient matrices for born charges and dielectric constants.
!!
use type_symmetryoperation, only: lo_expandoperation_pair, lo_expandoperation_triplet
implicit none
contains

!> coefficient matrix for Born charges in the unitcell
module subroutine lo_coeffmatrix_unitcell_Z_singlet(map, CM)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM

    real(r8), dimension(9, map%xuc%nx_Z_singlet) :: m0
    real(r8), dimension(9, 9) :: rotm
    integer :: a1, sh, op, nx

    CM = 0.0_r8
    do a1 = 1, map%n_atom_uc
        sh = map%xuc%Z_singlet(a1)%irreducible_shell
        op = map%xuc%Z_singlet(a1)%operation_from_shell
        rotm = lo_expandoperation_pair(map%op_singlet(op)%m3)
        nx = map%Z_singlet_shell(sh)%nx
        if (nx .eq. 0) cycle
        m0 = 0.0_r8
        call lo_gemm(rotm, map%Z_singlet_shell(sh)%coeff, m0(:, 1:nx))
        CM((a1 - 1)*9 + 1:a1*9, map%Z_singlet_shell(sh)%ind_global) = m0(:, 1:nx)
    end do
    CM = lo_chop(CM, 1E-13_r8)
end subroutine

!> coefficient matrix for Born charges in the supercell
module subroutine lo_coeffmatrix_supercell_Z_singlet(map, ss, CM)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM

    real(r8), dimension(9, map%xuc%nx_Z_singlet) :: m0
    real(r8), dimension(9, 9) :: rotm
    integer :: a1, a2, sh, op, nx

    CM = 0.0_r8
    do a1 = 1, map%n_atom_ss
        a2 = ss%info%index_in_unitcell(a1)
        sh = map%xuc%Z_singlet(a2)%irreducible_shell
        op = map%xuc%Z_singlet(a2)%operation_from_shell
        rotm = lo_expandoperation_pair(map%op_singlet(op)%m3)
        nx = map%Z_singlet_shell(sh)%nx
        if (nx .eq. 0) cycle
        m0 = 0.0_r8
        call lo_gemm(rotm, map%Z_singlet_shell(sh)%coeff, m0(:, 1:nx))
        CM((a1 - 1)*9 + 1:a1*9, map%Z_singlet_shell(sh)%ind_global) = m0(:, 1:nx)
    end do
    CM = lo_chop(CM, 1E-13_r8)
end subroutine

!> coefficient matrix for born charge pairs
module subroutine lo_coeffmatrix_Z_pair(map, UM, CM)
    !> symmetry things in shells
    type(lo_forcemap), intent(in) :: map
    !> displacements
    real(r8), dimension(:, :), intent(in) :: UM
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM

    real(r8), dimension(27, 27) :: B
    real(r8), dimension(9, 27) :: A, C
    real(r8), dimension(3) :: u, v
    integer, dimension(27) :: xind
    integer :: ipair, ish, iop, nx, a1, a2, ii, jj

    CM = 0.0_r8
    do ipair = 1, map%xss%n_Z_pair
        ish = map%xss%ind_Z_pair(1, ipair)
        iop = map%xss%ind_Z_pair(2, ipair)
        nx = map%Z_pair_shell(ish)%nx
        if (nx .eq. 0) cycle
        a1 = map%xss%ind_Z_pair(3, ipair)
        a2 = map%xss%ind_Z_pair(4, ipair)
        v = UM(:, a1)
        u = UM(:, a2)
        xind(1:nx) = map%Z_pair_shell(ish)%ind_global
        A = zpair_kronexpand(u, v)
        ! rotate coefficient matrix
        B = 0.0_r8
        call lo_gemm(map%op_triplet(iop)%sotr, map%Z_pair_shell(ish)%coeff, B(:, 1:nx))
        ! apply kronecker guy
        call lo_gemm(A, B(:, 1:nx), C(:, 1:nx))
        ! now C hold the coefficient matrix, store it in the right place
        ii = (a1 - 1)*9 + 1
        jj = a1*9
        CM(ii:jj, xind(1:nx)) = CM(ii:jj, xind(1:nx)) + C(:, 1:nx)
    end do
end subroutine

!> coefficient matrix for born charge triplets
module subroutine lo_coeffmatrix_Z_triplet(map, UM, CM)
    !> symmetry things in shells
    type(lo_forcemap), intent(in) :: map
    !> displacements
    real(r8), dimension(:, :), intent(in) :: UM
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM

    real(r8), dimension(81, 81) :: B
    real(r8), dimension(9, 81) :: A, C
    real(r8), dimension(3) :: u2, u3, v2, v3
    integer, dimension(81) :: xind
    integer :: itriplet, ish, iop, nx, a1, a2, a3, ii, jj

    CM = 0.0_r8
    do itriplet = 1, map%xss%n_Z_triplet
        ish = map%xss%ind_Z_triplet(1, itriplet)
        iop = map%xss%ind_Z_triplet(2, itriplet)
        nx = map%Z_triplet_shell(ish)%nx
        if (nx .eq. 0) cycle
        a1 = map%xss%ind_Z_triplet(3, itriplet)
        a2 = map%xss%ind_Z_triplet(4, itriplet)
        a3 = map%xss%ind_Z_triplet(5, itriplet)
        v2 = UM(:, a1)
        v3 = UM(:, a1)
        u2 = UM(:, a2)
        u3 = UM(:, a3)
        xind(1:nx) = map%Z_triplet_shell(ish)%ind_global
        A = ztriplet_kronexpand(u2, u3, v2, v3)
        ! rotate coefficient matrix
        B = 0.0_r8
        call lo_gemm(map%op_quartet(iop)%sotr, map%Z_triplet_shell(ish)%coeff, B(:, 1:nx))
        ! apply kronecker guy
        call lo_gemm(A, B(:, 1:nx), C(:, 1:nx))
        ! now C hold the coefficient matrix, store it in the right place
        ii = (a1 - 1)*9 + 1
        jj = a1*9
        CM(ii:jj, xind(1:nx)) = CM(ii:jj, xind(1:nx)) + C(:, 1:nx)
    end do
end subroutine

!> coefficient matrix for dielectric singlets in the supercell
module subroutine lo_coeffmatrix_eps_singlet(map, UM, CM)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> displacements
    real(r8), dimension(:, :), intent(in) :: UM
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM

    real(r8), dimension(27, map%xuc%nx_eps_singlet) :: B
    real(r8), dimension(9, map%xuc%nx_eps_singlet) :: C
    real(r8), dimension(27, 27) :: rotm
    real(r8), dimension(9, 27) :: A
    integer, dimension(map%xuc%nx_eps_singlet) :: xind
    integer :: a1, ish, iop, nx

    CM = 0.0_r8
    do a1 = 1, map%n_atom_ss
        ish = map%xss%ind_eps_singlet(1, a1)
        iop = map%xss%ind_eps_singlet(2, a1)
        nx = map%eps_singlet_shell(ish)%nx
        if (nx .eq. 0) cycle
        xind(1:nx) = map%eps_singlet_shell(ish)%ind_global
        ! rotate coefficient matrix
        rotm = lo_expandoperation_triplet(map%op_singlet(iop)%m3)
        call lo_gemm(rotm, map%eps_singlet_shell(ish)%coeff, B(:, 1:nx))
        ! kronecker displacement guy
        A = epssing_kronexpand(UM(:, a1))
        call lo_gemm(A, B(:, 1:nx), C(:, 1:nx))
        ! store in coefficient matrix
        CM(:, xind(1:nx)) = CM(:, xind(1:nx)) + C(:, 1:nx)
    end do
end subroutine

!> coefficient matrix for dielectric pairs
module subroutine lo_coeffmatrix_eps_pair(map, UM, CM)
    !> symmetry things in shells
    type(lo_forcemap), intent(in) :: map
    !> displacements
    real(r8), dimension(:, :), intent(in) :: UM
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM

    real(r8), dimension(81, 81) :: B
    real(r8), dimension(9, 81) :: A, C
    real(r8), dimension(3) :: u1, u2, v1, v2
    integer, dimension(81) :: xind
    integer :: ipair, ish, iop, nx, a1, a2

    CM = 0.0_r8
    do ipair = 1, map%xss%n_eps_pair
        ish = map%xss%ind_eps_pair(1, ipair)
        iop = map%xss%ind_eps_pair(2, ipair)
        nx = map%eps_pair_shell(ish)%nx
        if (nx .eq. 0) cycle
        a1 = map%xss%ind_eps_pair(3, ipair)
        a2 = map%xss%ind_eps_pair(4, ipair)
        v1 = UM(:, a1)
        v2 = UM(:, a1)
        u1 = UM(:, a1)
        u2 = UM(:, a2)
        xind(1:nx) = map%eps_pair_shell(ish)%ind_global
        A = epspair_kronexpand(u1, u2, v1, v2)
        ! rotate coefficient matrix
        B = 0.0_r8
        call lo_gemm(map%op_quartet(iop)%sotr, map%eps_pair_shell(ish)%coeff, B(:, 1:nx))
        ! apply kronecker guy
        call lo_gemm(A, B(:, 1:nx), C(:, 1:nx))
        ! now C hold the coefficient matrix, store it in the right place
        CM(:, xind(1:nx)) = CM(:, xind(1:nx)) + C(:, 1:nx)
    end do
end subroutine

! below are not exposed

! !> does the kronecker product I9 otimes u
pure function epssing_kronexpand(u) result(A)
    !> displacements
    real(r8), dimension(3), intent(in) :: u
    !> kronecker product
    real(r8), dimension(9, 27) :: A

    integer :: i
    A = 0.0_r8
    do i = 1, 9
        A(i, (i - 1)*3 + 1:i*3) = u
    end do
end function

pure function epspair_kronexpand(u1, u2, v1, v2) result(A)
    !> displacements
    real(r8), dimension(3), intent(in) :: u1, u2, v1, v2
    !> kronecker product
    real(r8), dimension(9, 81) :: A

    integer :: i, j, k, l
    A = 0.0_r8
    do i = 1, 9
        do j = 1, 3
        do k = 1, 3
            l = (i - 1)*9 + (j - 1)*3 + k
            A(i, l) = u1(j)*u2(k) - v1(j)*v2(k)
        end do
        end do
    end do
end function

!> does the kronecker product I9 otimes u
pure function zpair_kronexpand(u, v) result(A)
    !> displacements
    real(r8), dimension(3), intent(in) :: u, v
    !> kronecker product
    real(r8), dimension(9, 27) :: A

    integer :: i
    A = 0.0_r8
    do i = 1, 9
        A(i, (i - 1)*3 + 1:i*3) = u - v
    end do
end function

pure function ztriplet_kronexpand(u2, u3, v2, v3) result(A)
    !> displacements
    real(r8), dimension(3), intent(in) :: u2, u3, v2, v3
    !> kronecker product
    real(r8), dimension(9, 81) :: A

    integer :: i, j, k, l
    A = 0.0_r8
    do i = 1, 9
        do j = 1, 3
        do k = 1, 3
            l = (i - 1)*9 + (j - 1)*3 + k
            A(i, l) = u2(j)*u3(k) - v2(j)*v3(k)
        end do
        end do
    end do
end function

end submodule
