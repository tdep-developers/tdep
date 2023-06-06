#include "precompilerdefinitions"
submodule(type_forcemap) type_forcemap_coefficient_triplet
implicit none
contains

!> get the triplet coefficient matrix
module subroutine lo_coeffmatrix_triplet(UM, CM, map)
    !> the displacements
    real(r8), dimension(:, :), intent(in) :: UM
    !> the coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM
    !> symmetry stuff rearranged into coordination shells
    type(lo_forcemap), intent(in) :: map

    integer :: a1, a2, a3, i, k, ii
    integer :: find, uind, unsh, unop, nx
    integer, dimension(27) :: xind
    real(r8), dimension(3) :: u2, u3, v2, v3
    real(r8) :: wm27x1(27, 1), wm27x2(27, 2), wm27x3(27, 3), wm27x4(27, 4), wm27x5(27, 5), wm27x6(27, 6), wm27x7(27, 7), &
                wm27x8(27, 8), wm27x9(27, 9), wm27x10(27, 10), wm27x11(27, 11), wm27x12(27, 12), wm27x13(27, 13), &
                wm27x14(27, 14), wm27x15(27, 15), wm27x16(27, 16), wm27x17(27, 17), wm27x18(27, 18), wm27x19(27, 19), &
                wm27x20(27, 20), wm27x21(27, 21), wm27x22(27, 22), wm27x23(27, 23), wm27x24(27, 24), wm27x25(27, 25), &
                wm27x26(27, 26), wm27x27(27, 27)
    real(r8) :: cm1x3(1, 3), cm2x3(2, 3), cm3x3(3, 3), cm4x3(4, 3), cm5x3(5, 3), cm6x3(6, 3), cm7x3(7, 3), &
                cm8x3(8, 3), cm9x3(9, 3), cm10x3(10, 3), cm11x3(11, 3), cm12x3(12, 3), cm13x3(13, 3), cm14x3(14, 3), &
                cm15x3(15, 3), cm16x3(16, 3), cm17x3(17, 3), cm18x3(18, 3), cm19x3(19, 3), cm20x3(20, 3), cm21x3(21, 3), &
                cm22x3(22, 3), cm23x3(23, 3), cm24x3(24, 3), cm25x3(25, 3), cm26x3(26, 3), cm27x3(27, 3)

    CM = 0.0_r8
    do i = 1, map%xss%n_fc_triplet
        unsh = map%xss%ind_fc_triplet(1, i)
        unop = map%xss%ind_fc_triplet(2, i)
        a1 = map%xss%ind_fc_triplet(3, i)
        a2 = map%xss%ind_fc_triplet(4, i)
        a3 = map%xss%ind_fc_triplet(5, i)
        v2 = UM(:, a1)
        v3 = UM(:, a1)
        u2 = UM(:, a2)
        u3 = UM(:, a3)
        nx = map%fc_triplet_shell(unsh)%nx
        if (nx .eq. 0) cycle
        xind(1:nx) = map%fc_triplet_shell(unsh)%ind_global
        select case (nx)
        case (1)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x1)
            cm1x3 = coeff27x1(wm27x1, u2, u3) - coeff27x1(wm27x1, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm1x3(k, ii)
            end do
            end do
        case (2)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x2)
            cm2x3 = coeff27x2(wm27x2, u2, u3) - coeff27x2(wm27x2, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm2x3(k, ii)
            end do
            end do
        case (3)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x3)
            cm3x3 = coeff27x3(wm27x3, u2, u3) - coeff27x3(wm27x3, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm3x3(k, ii)
            end do
            end do
        case (4)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x4)
            cm4x3 = coeff27x4(wm27x4, u2, u3) - coeff27x4(wm27x4, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm4x3(k, ii)
            end do
            end do
        case (5)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x5)
            cm5x3 = coeff27x5(wm27x5, u2, u3) - coeff27x5(wm27x5, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm5x3(k, ii)
            end do
            end do
        case (6)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x6)
            cm6x3 = coeff27x6(wm27x6, u2, u3) - coeff27x6(wm27x6, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm6x3(k, ii)
            end do
            end do
        case (7)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x7)
            cm7x3 = coeff27x7(wm27x7, u2, u3) - coeff27x7(wm27x7, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm7x3(k, ii)
            end do
            end do
        case (8)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x8)
            cm8x3 = coeff27x8(wm27x8, u2, u3) - coeff27x8(wm27x8, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm8x3(k, ii)
            end do
            end do
        case (9)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x9)
            cm9x3 = coeff27x9(wm27x9, u2, u3) - coeff27x9(wm27x9, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm9x3(k, ii)
            end do
            end do
        case (10)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x10)
            cm10x3 = coeff27x10(wm27x10, u2, u3) - coeff27x10(wm27x10, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm10x3(k, ii)
            end do
            end do
        case (11)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x11)
            cm11x3 = coeff27x11(wm27x11, u2, u3) - coeff27x11(wm27x11, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm11x3(k, ii)
            end do
            end do
        case (12)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x12)
            cm12x3 = coeff27x12(wm27x12, u2, u3) - coeff27x12(wm27x12, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm12x3(k, ii)
            end do
            end do
        case (13)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x13)
            cm13x3 = coeff27x13(wm27x13, u2, u3) - coeff27x13(wm27x13, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm13x3(k, ii)
            end do
            end do
        case (14)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x14)
            cm14x3 = coeff27x14(wm27x14, u2, u3) - coeff27x14(wm27x14, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm14x3(k, ii)
            end do
            end do
        case (15)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x15)
            cm15x3 = coeff27x15(wm27x15, u2, u3) - coeff27x15(wm27x15, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm15x3(k, ii)
            end do
            end do
        case (16)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x16)
            cm16x3 = coeff27x16(wm27x16, u2, u3) - coeff27x16(wm27x16, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm16x3(k, ii)
            end do
            end do
        case (17)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x17)
            cm17x3 = coeff27x17(wm27x17, u2, u3) - coeff27x17(wm27x17, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm17x3(k, ii)
            end do
            end do
        case (18)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x18)
            cm18x3 = coeff27x18(wm27x18, u2, u3) - coeff27x18(wm27x18, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm18x3(k, ii)
            end do
            end do
        case (19)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x19)
            cm19x3 = coeff27x19(wm27x19, u2, u3) - coeff27x19(wm27x19, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm19x3(k, ii)
            end do
            end do
        case (20)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x20)
            cm20x3 = coeff27x20(wm27x20, u2, u3) - coeff27x20(wm27x20, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm20x3(k, ii)
            end do
            end do
        case (21)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x21)
            cm21x3 = coeff27x21(wm27x21, u2, u3) - coeff27x21(wm27x21, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm21x3(k, ii)
            end do
            end do
        case (22)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x22)
            cm22x3 = coeff27x22(wm27x22, u2, u3) - coeff27x22(wm27x22, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm22x3(k, ii)
            end do
            end do
        case (23)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x23)
            cm23x3 = coeff27x23(wm27x23, u2, u3) - coeff27x23(wm27x23, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm23x3(k, ii)
            end do
            end do
        case (24)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x24)
            cm24x3 = coeff27x24(wm27x24, u2, u3) - coeff27x24(wm27x24, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm24x3(k, ii)
            end do
            end do
        case (25)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x25)
            cm25x3 = coeff27x25(wm27x25, u2, u3) - coeff27x25(wm27x25, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm25x3(k, ii)
            end do
            end do
        case (26)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x26)
            cm26x3 = coeff27x26(wm27x26, u2, u3) - coeff27x26(wm27x26, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm26x3(k, ii)
            end do
            end do
        case (27)
            call lo_gemm(map%op_triplet(unop)%sotr, map%fc_triplet_shell(unsh)%coeff, wm27x27)
            cm27x3 = coeff27x27(wm27x27, u2, u3) - coeff27x27(wm27x27, v2, v3)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm27x3(k, ii)
            end do
            end do
        end select
    end do

contains

    pure function coeff27x1(cm, u2, u3) result(m)
        real(r8), dimension(27, 1), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(1, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x2(cm, u2, u3) result(m)
        real(r8), dimension(27, 2), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(2, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x3(cm, u2, u3) result(m)
        real(r8), dimension(27, 3), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(3, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x4(cm, u2, u3) result(m)
        real(r8), dimension(27, 4), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(4, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x5(cm, u2, u3) result(m)
        real(r8), dimension(27, 5), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(5, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x6(cm, u2, u3) result(m)
        real(r8), dimension(27, 6), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(6, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x7(cm, u2, u3) result(m)
        real(r8), dimension(27, 7), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(7, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x8(cm, u2, u3) result(m)
        real(r8), dimension(27, 8), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(8, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x9(cm, u2, u3) result(m)
        real(r8), dimension(27, 9), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(9, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x10(cm, u2, u3) result(m)
        real(r8), dimension(27, 10), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(10, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x11(cm, u2, u3) result(m)
        real(r8), dimension(27, 11), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(11, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x12(cm, u2, u3) result(m)
        real(r8), dimension(27, 12), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(12, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x13(cm, u2, u3) result(m)
        real(r8), dimension(27, 13), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(13, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x14(cm, u2, u3) result(m)
        real(r8), dimension(27, 14), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(14, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x15(cm, u2, u3) result(m)
        real(r8), dimension(27, 15), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(15, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x16(cm, u2, u3) result(m)
        real(r8), dimension(27, 16), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(16, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x17(cm, u2, u3) result(m)
        real(r8), dimension(27, 17), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(17, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x18(cm, u2, u3) result(m)
        real(r8), dimension(27, 18), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(18, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x19(cm, u2, u3) result(m)
        real(r8), dimension(27, 19), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(19, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x20(cm, u2, u3) result(m)
        real(r8), dimension(27, 20), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(20, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x21(cm, u2, u3) result(m)
        real(r8), dimension(27, 21), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(21, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(21, 1) = cm(1, 21)*u2x*u3x + cm(4, 21)*u2y*u3x + cm(7, 21)*u2z*u3x + cm(2, 21)*u2x*u3y + cm(5, 21)*u2y*u3y + cm(8, 21)*u2z*u3y + cm(3, 21)*u2x*u3z + cm(6, 21)*u2y*u3z + cm(9, 21)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(21, 2) = cm(10, 21)*u2x*u3x + cm(13, 21)*u2y*u3x + cm(16, 21)*u2z*u3x + cm(11, 21)*u2x*u3y + cm(14, 21)*u2y*u3y + cm(17, 21)*u2z*u3y + cm(12, 21)*u2x*u3z + cm(15, 21)*u2y*u3z + cm(18, 21)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m(21, 3) = cm(19, 21)*u2x*u3x + cm(22, 21)*u2y*u3x + cm(25, 21)*u2z*u3x + cm(20, 21)*u2x*u3y + cm(23, 21)*u2y*u3y + cm(26, 21)*u2z*u3y + cm(21, 21)*u2x*u3z + cm(24, 21)*u2y*u3z + cm(27, 21)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x22(cm, u2, u3) result(m)
        real(r8), dimension(27, 22), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(22, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(21, 1) = cm(1, 21)*u2x*u3x + cm(4, 21)*u2y*u3x + cm(7, 21)*u2z*u3x + cm(2, 21)*u2x*u3y + cm(5, 21)*u2y*u3y + cm(8, 21)*u2z*u3y + cm(3, 21)*u2x*u3z + cm(6, 21)*u2y*u3z + cm(9, 21)*u2z*u3z
        m(22, 1) = cm(1, 22)*u2x*u3x + cm(4, 22)*u2y*u3x + cm(7, 22)*u2z*u3x + cm(2, 22)*u2x*u3y + cm(5, 22)*u2y*u3y + cm(8, 22)*u2z*u3y + cm(3, 22)*u2x*u3z + cm(6, 22)*u2y*u3z + cm(9, 22)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(21, 2) = cm(10, 21)*u2x*u3x + cm(13, 21)*u2y*u3x + cm(16, 21)*u2z*u3x + cm(11, 21)*u2x*u3y + cm(14, 21)*u2y*u3y + cm(17, 21)*u2z*u3y + cm(12, 21)*u2x*u3z + cm(15, 21)*u2y*u3z + cm(18, 21)*u2z*u3z
        m(22, 2) = cm(10, 22)*u2x*u3x + cm(13, 22)*u2y*u3x + cm(16, 22)*u2z*u3x + cm(11, 22)*u2x*u3y + cm(14, 22)*u2y*u3y + cm(17, 22)*u2z*u3y + cm(12, 22)*u2x*u3z + cm(15, 22)*u2y*u3z + cm(18, 22)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m(21, 3) = cm(19, 21)*u2x*u3x + cm(22, 21)*u2y*u3x + cm(25, 21)*u2z*u3x + cm(20, 21)*u2x*u3y + cm(23, 21)*u2y*u3y + cm(26, 21)*u2z*u3y + cm(21, 21)*u2x*u3z + cm(24, 21)*u2y*u3z + cm(27, 21)*u2z*u3z
        m(22, 3) = cm(19, 22)*u2x*u3x + cm(22, 22)*u2y*u3x + cm(25, 22)*u2z*u3x + cm(20, 22)*u2x*u3y + cm(23, 22)*u2y*u3y + cm(26, 22)*u2z*u3y + cm(21, 22)*u2x*u3z + cm(24, 22)*u2y*u3z + cm(27, 22)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x23(cm, u2, u3) result(m)
        real(r8), dimension(27, 23), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(23, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(21, 1) = cm(1, 21)*u2x*u3x + cm(4, 21)*u2y*u3x + cm(7, 21)*u2z*u3x + cm(2, 21)*u2x*u3y + cm(5, 21)*u2y*u3y + cm(8, 21)*u2z*u3y + cm(3, 21)*u2x*u3z + cm(6, 21)*u2y*u3z + cm(9, 21)*u2z*u3z
        m(22, 1) = cm(1, 22)*u2x*u3x + cm(4, 22)*u2y*u3x + cm(7, 22)*u2z*u3x + cm(2, 22)*u2x*u3y + cm(5, 22)*u2y*u3y + cm(8, 22)*u2z*u3y + cm(3, 22)*u2x*u3z + cm(6, 22)*u2y*u3z + cm(9, 22)*u2z*u3z
        m(23, 1) = cm(1, 23)*u2x*u3x + cm(4, 23)*u2y*u3x + cm(7, 23)*u2z*u3x + cm(2, 23)*u2x*u3y + cm(5, 23)*u2y*u3y + cm(8, 23)*u2z*u3y + cm(3, 23)*u2x*u3z + cm(6, 23)*u2y*u3z + cm(9, 23)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(21, 2) = cm(10, 21)*u2x*u3x + cm(13, 21)*u2y*u3x + cm(16, 21)*u2z*u3x + cm(11, 21)*u2x*u3y + cm(14, 21)*u2y*u3y + cm(17, 21)*u2z*u3y + cm(12, 21)*u2x*u3z + cm(15, 21)*u2y*u3z + cm(18, 21)*u2z*u3z
        m(22, 2) = cm(10, 22)*u2x*u3x + cm(13, 22)*u2y*u3x + cm(16, 22)*u2z*u3x + cm(11, 22)*u2x*u3y + cm(14, 22)*u2y*u3y + cm(17, 22)*u2z*u3y + cm(12, 22)*u2x*u3z + cm(15, 22)*u2y*u3z + cm(18, 22)*u2z*u3z
        m(23, 2) = cm(10, 23)*u2x*u3x + cm(13, 23)*u2y*u3x + cm(16, 23)*u2z*u3x + cm(11, 23)*u2x*u3y + cm(14, 23)*u2y*u3y + cm(17, 23)*u2z*u3y + cm(12, 23)*u2x*u3z + cm(15, 23)*u2y*u3z + cm(18, 23)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m(21, 3) = cm(19, 21)*u2x*u3x + cm(22, 21)*u2y*u3x + cm(25, 21)*u2z*u3x + cm(20, 21)*u2x*u3y + cm(23, 21)*u2y*u3y + cm(26, 21)*u2z*u3y + cm(21, 21)*u2x*u3z + cm(24, 21)*u2y*u3z + cm(27, 21)*u2z*u3z
        m(22, 3) = cm(19, 22)*u2x*u3x + cm(22, 22)*u2y*u3x + cm(25, 22)*u2z*u3x + cm(20, 22)*u2x*u3y + cm(23, 22)*u2y*u3y + cm(26, 22)*u2z*u3y + cm(21, 22)*u2x*u3z + cm(24, 22)*u2y*u3z + cm(27, 22)*u2z*u3z
        m(23, 3) = cm(19, 23)*u2x*u3x + cm(22, 23)*u2y*u3x + cm(25, 23)*u2z*u3x + cm(20, 23)*u2x*u3y + cm(23, 23)*u2y*u3y + cm(26, 23)*u2z*u3y + cm(21, 23)*u2x*u3z + cm(24, 23)*u2y*u3z + cm(27, 23)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x24(cm, u2, u3) result(m)
        real(r8), dimension(27, 24), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(24, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(21, 1) = cm(1, 21)*u2x*u3x + cm(4, 21)*u2y*u3x + cm(7, 21)*u2z*u3x + cm(2, 21)*u2x*u3y + cm(5, 21)*u2y*u3y + cm(8, 21)*u2z*u3y + cm(3, 21)*u2x*u3z + cm(6, 21)*u2y*u3z + cm(9, 21)*u2z*u3z
        m(22, 1) = cm(1, 22)*u2x*u3x + cm(4, 22)*u2y*u3x + cm(7, 22)*u2z*u3x + cm(2, 22)*u2x*u3y + cm(5, 22)*u2y*u3y + cm(8, 22)*u2z*u3y + cm(3, 22)*u2x*u3z + cm(6, 22)*u2y*u3z + cm(9, 22)*u2z*u3z
        m(23, 1) = cm(1, 23)*u2x*u3x + cm(4, 23)*u2y*u3x + cm(7, 23)*u2z*u3x + cm(2, 23)*u2x*u3y + cm(5, 23)*u2y*u3y + cm(8, 23)*u2z*u3y + cm(3, 23)*u2x*u3z + cm(6, 23)*u2y*u3z + cm(9, 23)*u2z*u3z
        m(24, 1) = cm(1, 24)*u2x*u3x + cm(4, 24)*u2y*u3x + cm(7, 24)*u2z*u3x + cm(2, 24)*u2x*u3y + cm(5, 24)*u2y*u3y + cm(8, 24)*u2z*u3y + cm(3, 24)*u2x*u3z + cm(6, 24)*u2y*u3z + cm(9, 24)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(21, 2) = cm(10, 21)*u2x*u3x + cm(13, 21)*u2y*u3x + cm(16, 21)*u2z*u3x + cm(11, 21)*u2x*u3y + cm(14, 21)*u2y*u3y + cm(17, 21)*u2z*u3y + cm(12, 21)*u2x*u3z + cm(15, 21)*u2y*u3z + cm(18, 21)*u2z*u3z
        m(22, 2) = cm(10, 22)*u2x*u3x + cm(13, 22)*u2y*u3x + cm(16, 22)*u2z*u3x + cm(11, 22)*u2x*u3y + cm(14, 22)*u2y*u3y + cm(17, 22)*u2z*u3y + cm(12, 22)*u2x*u3z + cm(15, 22)*u2y*u3z + cm(18, 22)*u2z*u3z
        m(23, 2) = cm(10, 23)*u2x*u3x + cm(13, 23)*u2y*u3x + cm(16, 23)*u2z*u3x + cm(11, 23)*u2x*u3y + cm(14, 23)*u2y*u3y + cm(17, 23)*u2z*u3y + cm(12, 23)*u2x*u3z + cm(15, 23)*u2y*u3z + cm(18, 23)*u2z*u3z
        m(24, 2) = cm(10, 24)*u2x*u3x + cm(13, 24)*u2y*u3x + cm(16, 24)*u2z*u3x + cm(11, 24)*u2x*u3y + cm(14, 24)*u2y*u3y + cm(17, 24)*u2z*u3y + cm(12, 24)*u2x*u3z + cm(15, 24)*u2y*u3z + cm(18, 24)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m(21, 3) = cm(19, 21)*u2x*u3x + cm(22, 21)*u2y*u3x + cm(25, 21)*u2z*u3x + cm(20, 21)*u2x*u3y + cm(23, 21)*u2y*u3y + cm(26, 21)*u2z*u3y + cm(21, 21)*u2x*u3z + cm(24, 21)*u2y*u3z + cm(27, 21)*u2z*u3z
        m(22, 3) = cm(19, 22)*u2x*u3x + cm(22, 22)*u2y*u3x + cm(25, 22)*u2z*u3x + cm(20, 22)*u2x*u3y + cm(23, 22)*u2y*u3y + cm(26, 22)*u2z*u3y + cm(21, 22)*u2x*u3z + cm(24, 22)*u2y*u3z + cm(27, 22)*u2z*u3z
        m(23, 3) = cm(19, 23)*u2x*u3x + cm(22, 23)*u2y*u3x + cm(25, 23)*u2z*u3x + cm(20, 23)*u2x*u3y + cm(23, 23)*u2y*u3y + cm(26, 23)*u2z*u3y + cm(21, 23)*u2x*u3z + cm(24, 23)*u2y*u3z + cm(27, 23)*u2z*u3z
        m(24, 3) = cm(19, 24)*u2x*u3x + cm(22, 24)*u2y*u3x + cm(25, 24)*u2z*u3x + cm(20, 24)*u2x*u3y + cm(23, 24)*u2y*u3y + cm(26, 24)*u2z*u3y + cm(21, 24)*u2x*u3z + cm(24, 24)*u2y*u3z + cm(27, 24)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x25(cm, u2, u3) result(m)
        real(r8), dimension(27, 25), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(25, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(21, 1) = cm(1, 21)*u2x*u3x + cm(4, 21)*u2y*u3x + cm(7, 21)*u2z*u3x + cm(2, 21)*u2x*u3y + cm(5, 21)*u2y*u3y + cm(8, 21)*u2z*u3y + cm(3, 21)*u2x*u3z + cm(6, 21)*u2y*u3z + cm(9, 21)*u2z*u3z
        m(22, 1) = cm(1, 22)*u2x*u3x + cm(4, 22)*u2y*u3x + cm(7, 22)*u2z*u3x + cm(2, 22)*u2x*u3y + cm(5, 22)*u2y*u3y + cm(8, 22)*u2z*u3y + cm(3, 22)*u2x*u3z + cm(6, 22)*u2y*u3z + cm(9, 22)*u2z*u3z
        m(23, 1) = cm(1, 23)*u2x*u3x + cm(4, 23)*u2y*u3x + cm(7, 23)*u2z*u3x + cm(2, 23)*u2x*u3y + cm(5, 23)*u2y*u3y + cm(8, 23)*u2z*u3y + cm(3, 23)*u2x*u3z + cm(6, 23)*u2y*u3z + cm(9, 23)*u2z*u3z
        m(24, 1) = cm(1, 24)*u2x*u3x + cm(4, 24)*u2y*u3x + cm(7, 24)*u2z*u3x + cm(2, 24)*u2x*u3y + cm(5, 24)*u2y*u3y + cm(8, 24)*u2z*u3y + cm(3, 24)*u2x*u3z + cm(6, 24)*u2y*u3z + cm(9, 24)*u2z*u3z
        m(25, 1) = cm(1, 25)*u2x*u3x + cm(4, 25)*u2y*u3x + cm(7, 25)*u2z*u3x + cm(2, 25)*u2x*u3y + cm(5, 25)*u2y*u3y + cm(8, 25)*u2z*u3y + cm(3, 25)*u2x*u3z + cm(6, 25)*u2y*u3z + cm(9, 25)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(21, 2) = cm(10, 21)*u2x*u3x + cm(13, 21)*u2y*u3x + cm(16, 21)*u2z*u3x + cm(11, 21)*u2x*u3y + cm(14, 21)*u2y*u3y + cm(17, 21)*u2z*u3y + cm(12, 21)*u2x*u3z + cm(15, 21)*u2y*u3z + cm(18, 21)*u2z*u3z
        m(22, 2) = cm(10, 22)*u2x*u3x + cm(13, 22)*u2y*u3x + cm(16, 22)*u2z*u3x + cm(11, 22)*u2x*u3y + cm(14, 22)*u2y*u3y + cm(17, 22)*u2z*u3y + cm(12, 22)*u2x*u3z + cm(15, 22)*u2y*u3z + cm(18, 22)*u2z*u3z
        m(23, 2) = cm(10, 23)*u2x*u3x + cm(13, 23)*u2y*u3x + cm(16, 23)*u2z*u3x + cm(11, 23)*u2x*u3y + cm(14, 23)*u2y*u3y + cm(17, 23)*u2z*u3y + cm(12, 23)*u2x*u3z + cm(15, 23)*u2y*u3z + cm(18, 23)*u2z*u3z
        m(24, 2) = cm(10, 24)*u2x*u3x + cm(13, 24)*u2y*u3x + cm(16, 24)*u2z*u3x + cm(11, 24)*u2x*u3y + cm(14, 24)*u2y*u3y + cm(17, 24)*u2z*u3y + cm(12, 24)*u2x*u3z + cm(15, 24)*u2y*u3z + cm(18, 24)*u2z*u3z
        m(25, 2) = cm(10, 25)*u2x*u3x + cm(13, 25)*u2y*u3x + cm(16, 25)*u2z*u3x + cm(11, 25)*u2x*u3y + cm(14, 25)*u2y*u3y + cm(17, 25)*u2z*u3y + cm(12, 25)*u2x*u3z + cm(15, 25)*u2y*u3z + cm(18, 25)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m(21, 3) = cm(19, 21)*u2x*u3x + cm(22, 21)*u2y*u3x + cm(25, 21)*u2z*u3x + cm(20, 21)*u2x*u3y + cm(23, 21)*u2y*u3y + cm(26, 21)*u2z*u3y + cm(21, 21)*u2x*u3z + cm(24, 21)*u2y*u3z + cm(27, 21)*u2z*u3z
        m(22, 3) = cm(19, 22)*u2x*u3x + cm(22, 22)*u2y*u3x + cm(25, 22)*u2z*u3x + cm(20, 22)*u2x*u3y + cm(23, 22)*u2y*u3y + cm(26, 22)*u2z*u3y + cm(21, 22)*u2x*u3z + cm(24, 22)*u2y*u3z + cm(27, 22)*u2z*u3z
        m(23, 3) = cm(19, 23)*u2x*u3x + cm(22, 23)*u2y*u3x + cm(25, 23)*u2z*u3x + cm(20, 23)*u2x*u3y + cm(23, 23)*u2y*u3y + cm(26, 23)*u2z*u3y + cm(21, 23)*u2x*u3z + cm(24, 23)*u2y*u3z + cm(27, 23)*u2z*u3z
        m(24, 3) = cm(19, 24)*u2x*u3x + cm(22, 24)*u2y*u3x + cm(25, 24)*u2z*u3x + cm(20, 24)*u2x*u3y + cm(23, 24)*u2y*u3y + cm(26, 24)*u2z*u3y + cm(21, 24)*u2x*u3z + cm(24, 24)*u2y*u3z + cm(27, 24)*u2z*u3z
        m(25, 3) = cm(19, 25)*u2x*u3x + cm(22, 25)*u2y*u3x + cm(25, 25)*u2z*u3x + cm(20, 25)*u2x*u3y + cm(23, 25)*u2y*u3y + cm(26, 25)*u2z*u3y + cm(21, 25)*u2x*u3z + cm(24, 25)*u2y*u3z + cm(27, 25)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x26(cm, u2, u3) result(m)
        real(r8), dimension(27, 26), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(26, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(21, 1) = cm(1, 21)*u2x*u3x + cm(4, 21)*u2y*u3x + cm(7, 21)*u2z*u3x + cm(2, 21)*u2x*u3y + cm(5, 21)*u2y*u3y + cm(8, 21)*u2z*u3y + cm(3, 21)*u2x*u3z + cm(6, 21)*u2y*u3z + cm(9, 21)*u2z*u3z
        m(22, 1) = cm(1, 22)*u2x*u3x + cm(4, 22)*u2y*u3x + cm(7, 22)*u2z*u3x + cm(2, 22)*u2x*u3y + cm(5, 22)*u2y*u3y + cm(8, 22)*u2z*u3y + cm(3, 22)*u2x*u3z + cm(6, 22)*u2y*u3z + cm(9, 22)*u2z*u3z
        m(23, 1) = cm(1, 23)*u2x*u3x + cm(4, 23)*u2y*u3x + cm(7, 23)*u2z*u3x + cm(2, 23)*u2x*u3y + cm(5, 23)*u2y*u3y + cm(8, 23)*u2z*u3y + cm(3, 23)*u2x*u3z + cm(6, 23)*u2y*u3z + cm(9, 23)*u2z*u3z
        m(24, 1) = cm(1, 24)*u2x*u3x + cm(4, 24)*u2y*u3x + cm(7, 24)*u2z*u3x + cm(2, 24)*u2x*u3y + cm(5, 24)*u2y*u3y + cm(8, 24)*u2z*u3y + cm(3, 24)*u2x*u3z + cm(6, 24)*u2y*u3z + cm(9, 24)*u2z*u3z
        m(25, 1) = cm(1, 25)*u2x*u3x + cm(4, 25)*u2y*u3x + cm(7, 25)*u2z*u3x + cm(2, 25)*u2x*u3y + cm(5, 25)*u2y*u3y + cm(8, 25)*u2z*u3y + cm(3, 25)*u2x*u3z + cm(6, 25)*u2y*u3z + cm(9, 25)*u2z*u3z
        m(26, 1) = cm(1, 26)*u2x*u3x + cm(4, 26)*u2y*u3x + cm(7, 26)*u2z*u3x + cm(2, 26)*u2x*u3y + cm(5, 26)*u2y*u3y + cm(8, 26)*u2z*u3y + cm(3, 26)*u2x*u3z + cm(6, 26)*u2y*u3z + cm(9, 26)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(21, 2) = cm(10, 21)*u2x*u3x + cm(13, 21)*u2y*u3x + cm(16, 21)*u2z*u3x + cm(11, 21)*u2x*u3y + cm(14, 21)*u2y*u3y + cm(17, 21)*u2z*u3y + cm(12, 21)*u2x*u3z + cm(15, 21)*u2y*u3z + cm(18, 21)*u2z*u3z
        m(22, 2) = cm(10, 22)*u2x*u3x + cm(13, 22)*u2y*u3x + cm(16, 22)*u2z*u3x + cm(11, 22)*u2x*u3y + cm(14, 22)*u2y*u3y + cm(17, 22)*u2z*u3y + cm(12, 22)*u2x*u3z + cm(15, 22)*u2y*u3z + cm(18, 22)*u2z*u3z
        m(23, 2) = cm(10, 23)*u2x*u3x + cm(13, 23)*u2y*u3x + cm(16, 23)*u2z*u3x + cm(11, 23)*u2x*u3y + cm(14, 23)*u2y*u3y + cm(17, 23)*u2z*u3y + cm(12, 23)*u2x*u3z + cm(15, 23)*u2y*u3z + cm(18, 23)*u2z*u3z
        m(24, 2) = cm(10, 24)*u2x*u3x + cm(13, 24)*u2y*u3x + cm(16, 24)*u2z*u3x + cm(11, 24)*u2x*u3y + cm(14, 24)*u2y*u3y + cm(17, 24)*u2z*u3y + cm(12, 24)*u2x*u3z + cm(15, 24)*u2y*u3z + cm(18, 24)*u2z*u3z
        m(25, 2) = cm(10, 25)*u2x*u3x + cm(13, 25)*u2y*u3x + cm(16, 25)*u2z*u3x + cm(11, 25)*u2x*u3y + cm(14, 25)*u2y*u3y + cm(17, 25)*u2z*u3y + cm(12, 25)*u2x*u3z + cm(15, 25)*u2y*u3z + cm(18, 25)*u2z*u3z
        m(26, 2) = cm(10, 26)*u2x*u3x + cm(13, 26)*u2y*u3x + cm(16, 26)*u2z*u3x + cm(11, 26)*u2x*u3y + cm(14, 26)*u2y*u3y + cm(17, 26)*u2z*u3y + cm(12, 26)*u2x*u3z + cm(15, 26)*u2y*u3z + cm(18, 26)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m(21, 3) = cm(19, 21)*u2x*u3x + cm(22, 21)*u2y*u3x + cm(25, 21)*u2z*u3x + cm(20, 21)*u2x*u3y + cm(23, 21)*u2y*u3y + cm(26, 21)*u2z*u3y + cm(21, 21)*u2x*u3z + cm(24, 21)*u2y*u3z + cm(27, 21)*u2z*u3z
        m(22, 3) = cm(19, 22)*u2x*u3x + cm(22, 22)*u2y*u3x + cm(25, 22)*u2z*u3x + cm(20, 22)*u2x*u3y + cm(23, 22)*u2y*u3y + cm(26, 22)*u2z*u3y + cm(21, 22)*u2x*u3z + cm(24, 22)*u2y*u3z + cm(27, 22)*u2z*u3z
        m(23, 3) = cm(19, 23)*u2x*u3x + cm(22, 23)*u2y*u3x + cm(25, 23)*u2z*u3x + cm(20, 23)*u2x*u3y + cm(23, 23)*u2y*u3y + cm(26, 23)*u2z*u3y + cm(21, 23)*u2x*u3z + cm(24, 23)*u2y*u3z + cm(27, 23)*u2z*u3z
        m(24, 3) = cm(19, 24)*u2x*u3x + cm(22, 24)*u2y*u3x + cm(25, 24)*u2z*u3x + cm(20, 24)*u2x*u3y + cm(23, 24)*u2y*u3y + cm(26, 24)*u2z*u3y + cm(21, 24)*u2x*u3z + cm(24, 24)*u2y*u3z + cm(27, 24)*u2z*u3z
        m(25, 3) = cm(19, 25)*u2x*u3x + cm(22, 25)*u2y*u3x + cm(25, 25)*u2z*u3x + cm(20, 25)*u2x*u3y + cm(23, 25)*u2y*u3y + cm(26, 25)*u2z*u3y + cm(21, 25)*u2x*u3z + cm(24, 25)*u2y*u3z + cm(27, 25)*u2z*u3z
        m(26, 3) = cm(19, 26)*u2x*u3x + cm(22, 26)*u2y*u3x + cm(25, 26)*u2z*u3x + cm(20, 26)*u2x*u3y + cm(23, 26)*u2y*u3y + cm(26, 26)*u2z*u3y + cm(21, 26)*u2x*u3z + cm(24, 26)*u2y*u3z + cm(27, 26)*u2z*u3z
        m = m*0.5_r8
    end function

    pure function coeff27x27(cm, u2, u3) result(m)
        real(r8), dimension(27, 27), intent(in) :: cm
        real(r8), dimension(3), intent(in) :: u2, u3
        real(r8), dimension(27, 3) :: m
        real(r8) :: u2x, u2y, u2z
        real(r8) :: u3x, u3y, u3z
        u2x = u2(1); u2y = u2(2); u2z = u2(3)
        u3x = u3(1); u3y = u3(2); u3z = u3(3)
        m(1, 1) = cm(1, 1)*u2x*u3x + cm(4, 1)*u2y*u3x + cm(7, 1)*u2z*u3x + cm(2, 1)*u2x*u3y + cm(5, 1)*u2y*u3y + cm(8, 1)*u2z*u3y + cm(3, 1)*u2x*u3z + cm(6, 1)*u2y*u3z + cm(9, 1)*u2z*u3z
        m(2, 1) = cm(1, 2)*u2x*u3x + cm(4, 2)*u2y*u3x + cm(7, 2)*u2z*u3x + cm(2, 2)*u2x*u3y + cm(5, 2)*u2y*u3y + cm(8, 2)*u2z*u3y + cm(3, 2)*u2x*u3z + cm(6, 2)*u2y*u3z + cm(9, 2)*u2z*u3z
        m(3, 1) = cm(1, 3)*u2x*u3x + cm(4, 3)*u2y*u3x + cm(7, 3)*u2z*u3x + cm(2, 3)*u2x*u3y + cm(5, 3)*u2y*u3y + cm(8, 3)*u2z*u3y + cm(3, 3)*u2x*u3z + cm(6, 3)*u2y*u3z + cm(9, 3)*u2z*u3z
        m(4, 1) = cm(1, 4)*u2x*u3x + cm(4, 4)*u2y*u3x + cm(7, 4)*u2z*u3x + cm(2, 4)*u2x*u3y + cm(5, 4)*u2y*u3y + cm(8, 4)*u2z*u3y + cm(3, 4)*u2x*u3z + cm(6, 4)*u2y*u3z + cm(9, 4)*u2z*u3z
        m(5, 1) = cm(1, 5)*u2x*u3x + cm(4, 5)*u2y*u3x + cm(7, 5)*u2z*u3x + cm(2, 5)*u2x*u3y + cm(5, 5)*u2y*u3y + cm(8, 5)*u2z*u3y + cm(3, 5)*u2x*u3z + cm(6, 5)*u2y*u3z + cm(9, 5)*u2z*u3z
        m(6, 1) = cm(1, 6)*u2x*u3x + cm(4, 6)*u2y*u3x + cm(7, 6)*u2z*u3x + cm(2, 6)*u2x*u3y + cm(5, 6)*u2y*u3y + cm(8, 6)*u2z*u3y + cm(3, 6)*u2x*u3z + cm(6, 6)*u2y*u3z + cm(9, 6)*u2z*u3z
        m(7, 1) = cm(1, 7)*u2x*u3x + cm(4, 7)*u2y*u3x + cm(7, 7)*u2z*u3x + cm(2, 7)*u2x*u3y + cm(5, 7)*u2y*u3y + cm(8, 7)*u2z*u3y + cm(3, 7)*u2x*u3z + cm(6, 7)*u2y*u3z + cm(9, 7)*u2z*u3z
        m(8, 1) = cm(1, 8)*u2x*u3x + cm(4, 8)*u2y*u3x + cm(7, 8)*u2z*u3x + cm(2, 8)*u2x*u3y + cm(5, 8)*u2y*u3y + cm(8, 8)*u2z*u3y + cm(3, 8)*u2x*u3z + cm(6, 8)*u2y*u3z + cm(9, 8)*u2z*u3z
        m(9, 1) = cm(1, 9)*u2x*u3x + cm(4, 9)*u2y*u3x + cm(7, 9)*u2z*u3x + cm(2, 9)*u2x*u3y + cm(5, 9)*u2y*u3y + cm(8, 9)*u2z*u3y + cm(3, 9)*u2x*u3z + cm(6, 9)*u2y*u3z + cm(9, 9)*u2z*u3z
        m(10, 1) = cm(1, 10)*u2x*u3x + cm(4, 10)*u2y*u3x + cm(7, 10)*u2z*u3x + cm(2, 10)*u2x*u3y + cm(5, 10)*u2y*u3y + cm(8, 10)*u2z*u3y + cm(3, 10)*u2x*u3z + cm(6, 10)*u2y*u3z + cm(9, 10)*u2z*u3z
        m(11, 1) = cm(1, 11)*u2x*u3x + cm(4, 11)*u2y*u3x + cm(7, 11)*u2z*u3x + cm(2, 11)*u2x*u3y + cm(5, 11)*u2y*u3y + cm(8, 11)*u2z*u3y + cm(3, 11)*u2x*u3z + cm(6, 11)*u2y*u3z + cm(9, 11)*u2z*u3z
        m(12, 1) = cm(1, 12)*u2x*u3x + cm(4, 12)*u2y*u3x + cm(7, 12)*u2z*u3x + cm(2, 12)*u2x*u3y + cm(5, 12)*u2y*u3y + cm(8, 12)*u2z*u3y + cm(3, 12)*u2x*u3z + cm(6, 12)*u2y*u3z + cm(9, 12)*u2z*u3z
        m(13, 1) = cm(1, 13)*u2x*u3x + cm(4, 13)*u2y*u3x + cm(7, 13)*u2z*u3x + cm(2, 13)*u2x*u3y + cm(5, 13)*u2y*u3y + cm(8, 13)*u2z*u3y + cm(3, 13)*u2x*u3z + cm(6, 13)*u2y*u3z + cm(9, 13)*u2z*u3z
        m(14, 1) = cm(1, 14)*u2x*u3x + cm(4, 14)*u2y*u3x + cm(7, 14)*u2z*u3x + cm(2, 14)*u2x*u3y + cm(5, 14)*u2y*u3y + cm(8, 14)*u2z*u3y + cm(3, 14)*u2x*u3z + cm(6, 14)*u2y*u3z + cm(9, 14)*u2z*u3z
        m(15, 1) = cm(1, 15)*u2x*u3x + cm(4, 15)*u2y*u3x + cm(7, 15)*u2z*u3x + cm(2, 15)*u2x*u3y + cm(5, 15)*u2y*u3y + cm(8, 15)*u2z*u3y + cm(3, 15)*u2x*u3z + cm(6, 15)*u2y*u3z + cm(9, 15)*u2z*u3z
        m(16, 1) = cm(1, 16)*u2x*u3x + cm(4, 16)*u2y*u3x + cm(7, 16)*u2z*u3x + cm(2, 16)*u2x*u3y + cm(5, 16)*u2y*u3y + cm(8, 16)*u2z*u3y + cm(3, 16)*u2x*u3z + cm(6, 16)*u2y*u3z + cm(9, 16)*u2z*u3z
        m(17, 1) = cm(1, 17)*u2x*u3x + cm(4, 17)*u2y*u3x + cm(7, 17)*u2z*u3x + cm(2, 17)*u2x*u3y + cm(5, 17)*u2y*u3y + cm(8, 17)*u2z*u3y + cm(3, 17)*u2x*u3z + cm(6, 17)*u2y*u3z + cm(9, 17)*u2z*u3z
        m(18, 1) = cm(1, 18)*u2x*u3x + cm(4, 18)*u2y*u3x + cm(7, 18)*u2z*u3x + cm(2, 18)*u2x*u3y + cm(5, 18)*u2y*u3y + cm(8, 18)*u2z*u3y + cm(3, 18)*u2x*u3z + cm(6, 18)*u2y*u3z + cm(9, 18)*u2z*u3z
        m(19, 1) = cm(1, 19)*u2x*u3x + cm(4, 19)*u2y*u3x + cm(7, 19)*u2z*u3x + cm(2, 19)*u2x*u3y + cm(5, 19)*u2y*u3y + cm(8, 19)*u2z*u3y + cm(3, 19)*u2x*u3z + cm(6, 19)*u2y*u3z + cm(9, 19)*u2z*u3z
        m(20, 1) = cm(1, 20)*u2x*u3x + cm(4, 20)*u2y*u3x + cm(7, 20)*u2z*u3x + cm(2, 20)*u2x*u3y + cm(5, 20)*u2y*u3y + cm(8, 20)*u2z*u3y + cm(3, 20)*u2x*u3z + cm(6, 20)*u2y*u3z + cm(9, 20)*u2z*u3z
        m(21, 1) = cm(1, 21)*u2x*u3x + cm(4, 21)*u2y*u3x + cm(7, 21)*u2z*u3x + cm(2, 21)*u2x*u3y + cm(5, 21)*u2y*u3y + cm(8, 21)*u2z*u3y + cm(3, 21)*u2x*u3z + cm(6, 21)*u2y*u3z + cm(9, 21)*u2z*u3z
        m(22, 1) = cm(1, 22)*u2x*u3x + cm(4, 22)*u2y*u3x + cm(7, 22)*u2z*u3x + cm(2, 22)*u2x*u3y + cm(5, 22)*u2y*u3y + cm(8, 22)*u2z*u3y + cm(3, 22)*u2x*u3z + cm(6, 22)*u2y*u3z + cm(9, 22)*u2z*u3z
        m(23, 1) = cm(1, 23)*u2x*u3x + cm(4, 23)*u2y*u3x + cm(7, 23)*u2z*u3x + cm(2, 23)*u2x*u3y + cm(5, 23)*u2y*u3y + cm(8, 23)*u2z*u3y + cm(3, 23)*u2x*u3z + cm(6, 23)*u2y*u3z + cm(9, 23)*u2z*u3z
        m(24, 1) = cm(1, 24)*u2x*u3x + cm(4, 24)*u2y*u3x + cm(7, 24)*u2z*u3x + cm(2, 24)*u2x*u3y + cm(5, 24)*u2y*u3y + cm(8, 24)*u2z*u3y + cm(3, 24)*u2x*u3z + cm(6, 24)*u2y*u3z + cm(9, 24)*u2z*u3z
        m(25, 1) = cm(1, 25)*u2x*u3x + cm(4, 25)*u2y*u3x + cm(7, 25)*u2z*u3x + cm(2, 25)*u2x*u3y + cm(5, 25)*u2y*u3y + cm(8, 25)*u2z*u3y + cm(3, 25)*u2x*u3z + cm(6, 25)*u2y*u3z + cm(9, 25)*u2z*u3z
        m(26, 1) = cm(1, 26)*u2x*u3x + cm(4, 26)*u2y*u3x + cm(7, 26)*u2z*u3x + cm(2, 26)*u2x*u3y + cm(5, 26)*u2y*u3y + cm(8, 26)*u2z*u3y + cm(3, 26)*u2x*u3z + cm(6, 26)*u2y*u3z + cm(9, 26)*u2z*u3z
        m(27, 1) = cm(1, 27)*u2x*u3x + cm(4, 27)*u2y*u3x + cm(7, 27)*u2z*u3x + cm(2, 27)*u2x*u3y + cm(5, 27)*u2y*u3y + cm(8, 27)*u2z*u3y + cm(3, 27)*u2x*u3z + cm(6, 27)*u2y*u3z + cm(9, 27)*u2z*u3z
        m(1, 2) = cm(10, 1)*u2x*u3x + cm(13, 1)*u2y*u3x + cm(16, 1)*u2z*u3x + cm(11, 1)*u2x*u3y + cm(14, 1)*u2y*u3y + cm(17, 1)*u2z*u3y + cm(12, 1)*u2x*u3z + cm(15, 1)*u2y*u3z + cm(18, 1)*u2z*u3z
        m(2, 2) = cm(10, 2)*u2x*u3x + cm(13, 2)*u2y*u3x + cm(16, 2)*u2z*u3x + cm(11, 2)*u2x*u3y + cm(14, 2)*u2y*u3y + cm(17, 2)*u2z*u3y + cm(12, 2)*u2x*u3z + cm(15, 2)*u2y*u3z + cm(18, 2)*u2z*u3z
        m(3, 2) = cm(10, 3)*u2x*u3x + cm(13, 3)*u2y*u3x + cm(16, 3)*u2z*u3x + cm(11, 3)*u2x*u3y + cm(14, 3)*u2y*u3y + cm(17, 3)*u2z*u3y + cm(12, 3)*u2x*u3z + cm(15, 3)*u2y*u3z + cm(18, 3)*u2z*u3z
        m(4, 2) = cm(10, 4)*u2x*u3x + cm(13, 4)*u2y*u3x + cm(16, 4)*u2z*u3x + cm(11, 4)*u2x*u3y + cm(14, 4)*u2y*u3y + cm(17, 4)*u2z*u3y + cm(12, 4)*u2x*u3z + cm(15, 4)*u2y*u3z + cm(18, 4)*u2z*u3z
        m(5, 2) = cm(10, 5)*u2x*u3x + cm(13, 5)*u2y*u3x + cm(16, 5)*u2z*u3x + cm(11, 5)*u2x*u3y + cm(14, 5)*u2y*u3y + cm(17, 5)*u2z*u3y + cm(12, 5)*u2x*u3z + cm(15, 5)*u2y*u3z + cm(18, 5)*u2z*u3z
        m(6, 2) = cm(10, 6)*u2x*u3x + cm(13, 6)*u2y*u3x + cm(16, 6)*u2z*u3x + cm(11, 6)*u2x*u3y + cm(14, 6)*u2y*u3y + cm(17, 6)*u2z*u3y + cm(12, 6)*u2x*u3z + cm(15, 6)*u2y*u3z + cm(18, 6)*u2z*u3z
        m(7, 2) = cm(10, 7)*u2x*u3x + cm(13, 7)*u2y*u3x + cm(16, 7)*u2z*u3x + cm(11, 7)*u2x*u3y + cm(14, 7)*u2y*u3y + cm(17, 7)*u2z*u3y + cm(12, 7)*u2x*u3z + cm(15, 7)*u2y*u3z + cm(18, 7)*u2z*u3z
        m(8, 2) = cm(10, 8)*u2x*u3x + cm(13, 8)*u2y*u3x + cm(16, 8)*u2z*u3x + cm(11, 8)*u2x*u3y + cm(14, 8)*u2y*u3y + cm(17, 8)*u2z*u3y + cm(12, 8)*u2x*u3z + cm(15, 8)*u2y*u3z + cm(18, 8)*u2z*u3z
        m(9, 2) = cm(10, 9)*u2x*u3x + cm(13, 9)*u2y*u3x + cm(16, 9)*u2z*u3x + cm(11, 9)*u2x*u3y + cm(14, 9)*u2y*u3y + cm(17, 9)*u2z*u3y + cm(12, 9)*u2x*u3z + cm(15, 9)*u2y*u3z + cm(18, 9)*u2z*u3z
        m(10, 2) = cm(10, 10)*u2x*u3x + cm(13, 10)*u2y*u3x + cm(16, 10)*u2z*u3x + cm(11, 10)*u2x*u3y + cm(14, 10)*u2y*u3y + cm(17, 10)*u2z*u3y + cm(12, 10)*u2x*u3z + cm(15, 10)*u2y*u3z + cm(18, 10)*u2z*u3z
        m(11, 2) = cm(10, 11)*u2x*u3x + cm(13, 11)*u2y*u3x + cm(16, 11)*u2z*u3x + cm(11, 11)*u2x*u3y + cm(14, 11)*u2y*u3y + cm(17, 11)*u2z*u3y + cm(12, 11)*u2x*u3z + cm(15, 11)*u2y*u3z + cm(18, 11)*u2z*u3z
        m(12, 2) = cm(10, 12)*u2x*u3x + cm(13, 12)*u2y*u3x + cm(16, 12)*u2z*u3x + cm(11, 12)*u2x*u3y + cm(14, 12)*u2y*u3y + cm(17, 12)*u2z*u3y + cm(12, 12)*u2x*u3z + cm(15, 12)*u2y*u3z + cm(18, 12)*u2z*u3z
        m(13, 2) = cm(10, 13)*u2x*u3x + cm(13, 13)*u2y*u3x + cm(16, 13)*u2z*u3x + cm(11, 13)*u2x*u3y + cm(14, 13)*u2y*u3y + cm(17, 13)*u2z*u3y + cm(12, 13)*u2x*u3z + cm(15, 13)*u2y*u3z + cm(18, 13)*u2z*u3z
        m(14, 2) = cm(10, 14)*u2x*u3x + cm(13, 14)*u2y*u3x + cm(16, 14)*u2z*u3x + cm(11, 14)*u2x*u3y + cm(14, 14)*u2y*u3y + cm(17, 14)*u2z*u3y + cm(12, 14)*u2x*u3z + cm(15, 14)*u2y*u3z + cm(18, 14)*u2z*u3z
        m(15, 2) = cm(10, 15)*u2x*u3x + cm(13, 15)*u2y*u3x + cm(16, 15)*u2z*u3x + cm(11, 15)*u2x*u3y + cm(14, 15)*u2y*u3y + cm(17, 15)*u2z*u3y + cm(12, 15)*u2x*u3z + cm(15, 15)*u2y*u3z + cm(18, 15)*u2z*u3z
        m(16, 2) = cm(10, 16)*u2x*u3x + cm(13, 16)*u2y*u3x + cm(16, 16)*u2z*u3x + cm(11, 16)*u2x*u3y + cm(14, 16)*u2y*u3y + cm(17, 16)*u2z*u3y + cm(12, 16)*u2x*u3z + cm(15, 16)*u2y*u3z + cm(18, 16)*u2z*u3z
        m(17, 2) = cm(10, 17)*u2x*u3x + cm(13, 17)*u2y*u3x + cm(16, 17)*u2z*u3x + cm(11, 17)*u2x*u3y + cm(14, 17)*u2y*u3y + cm(17, 17)*u2z*u3y + cm(12, 17)*u2x*u3z + cm(15, 17)*u2y*u3z + cm(18, 17)*u2z*u3z
        m(18, 2) = cm(10, 18)*u2x*u3x + cm(13, 18)*u2y*u3x + cm(16, 18)*u2z*u3x + cm(11, 18)*u2x*u3y + cm(14, 18)*u2y*u3y + cm(17, 18)*u2z*u3y + cm(12, 18)*u2x*u3z + cm(15, 18)*u2y*u3z + cm(18, 18)*u2z*u3z
        m(19, 2) = cm(10, 19)*u2x*u3x + cm(13, 19)*u2y*u3x + cm(16, 19)*u2z*u3x + cm(11, 19)*u2x*u3y + cm(14, 19)*u2y*u3y + cm(17, 19)*u2z*u3y + cm(12, 19)*u2x*u3z + cm(15, 19)*u2y*u3z + cm(18, 19)*u2z*u3z
        m(20, 2) = cm(10, 20)*u2x*u3x + cm(13, 20)*u2y*u3x + cm(16, 20)*u2z*u3x + cm(11, 20)*u2x*u3y + cm(14, 20)*u2y*u3y + cm(17, 20)*u2z*u3y + cm(12, 20)*u2x*u3z + cm(15, 20)*u2y*u3z + cm(18, 20)*u2z*u3z
        m(21, 2) = cm(10, 21)*u2x*u3x + cm(13, 21)*u2y*u3x + cm(16, 21)*u2z*u3x + cm(11, 21)*u2x*u3y + cm(14, 21)*u2y*u3y + cm(17, 21)*u2z*u3y + cm(12, 21)*u2x*u3z + cm(15, 21)*u2y*u3z + cm(18, 21)*u2z*u3z
        m(22, 2) = cm(10, 22)*u2x*u3x + cm(13, 22)*u2y*u3x + cm(16, 22)*u2z*u3x + cm(11, 22)*u2x*u3y + cm(14, 22)*u2y*u3y + cm(17, 22)*u2z*u3y + cm(12, 22)*u2x*u3z + cm(15, 22)*u2y*u3z + cm(18, 22)*u2z*u3z
        m(23, 2) = cm(10, 23)*u2x*u3x + cm(13, 23)*u2y*u3x + cm(16, 23)*u2z*u3x + cm(11, 23)*u2x*u3y + cm(14, 23)*u2y*u3y + cm(17, 23)*u2z*u3y + cm(12, 23)*u2x*u3z + cm(15, 23)*u2y*u3z + cm(18, 23)*u2z*u3z
        m(24, 2) = cm(10, 24)*u2x*u3x + cm(13, 24)*u2y*u3x + cm(16, 24)*u2z*u3x + cm(11, 24)*u2x*u3y + cm(14, 24)*u2y*u3y + cm(17, 24)*u2z*u3y + cm(12, 24)*u2x*u3z + cm(15, 24)*u2y*u3z + cm(18, 24)*u2z*u3z
        m(25, 2) = cm(10, 25)*u2x*u3x + cm(13, 25)*u2y*u3x + cm(16, 25)*u2z*u3x + cm(11, 25)*u2x*u3y + cm(14, 25)*u2y*u3y + cm(17, 25)*u2z*u3y + cm(12, 25)*u2x*u3z + cm(15, 25)*u2y*u3z + cm(18, 25)*u2z*u3z
        m(26, 2) = cm(10, 26)*u2x*u3x + cm(13, 26)*u2y*u3x + cm(16, 26)*u2z*u3x + cm(11, 26)*u2x*u3y + cm(14, 26)*u2y*u3y + cm(17, 26)*u2z*u3y + cm(12, 26)*u2x*u3z + cm(15, 26)*u2y*u3z + cm(18, 26)*u2z*u3z
        m(27, 2) = cm(10, 27)*u2x*u3x + cm(13, 27)*u2y*u3x + cm(16, 27)*u2z*u3x + cm(11, 27)*u2x*u3y + cm(14, 27)*u2y*u3y + cm(17, 27)*u2z*u3y + cm(12, 27)*u2x*u3z + cm(15, 27)*u2y*u3z + cm(18, 27)*u2z*u3z
        m(1, 3) = cm(19, 1)*u2x*u3x + cm(22, 1)*u2y*u3x + cm(25, 1)*u2z*u3x + cm(20, 1)*u2x*u3y + cm(23, 1)*u2y*u3y + cm(26, 1)*u2z*u3y + cm(21, 1)*u2x*u3z + cm(24, 1)*u2y*u3z + cm(27, 1)*u2z*u3z
        m(2, 3) = cm(19, 2)*u2x*u3x + cm(22, 2)*u2y*u3x + cm(25, 2)*u2z*u3x + cm(20, 2)*u2x*u3y + cm(23, 2)*u2y*u3y + cm(26, 2)*u2z*u3y + cm(21, 2)*u2x*u3z + cm(24, 2)*u2y*u3z + cm(27, 2)*u2z*u3z
        m(3, 3) = cm(19, 3)*u2x*u3x + cm(22, 3)*u2y*u3x + cm(25, 3)*u2z*u3x + cm(20, 3)*u2x*u3y + cm(23, 3)*u2y*u3y + cm(26, 3)*u2z*u3y + cm(21, 3)*u2x*u3z + cm(24, 3)*u2y*u3z + cm(27, 3)*u2z*u3z
        m(4, 3) = cm(19, 4)*u2x*u3x + cm(22, 4)*u2y*u3x + cm(25, 4)*u2z*u3x + cm(20, 4)*u2x*u3y + cm(23, 4)*u2y*u3y + cm(26, 4)*u2z*u3y + cm(21, 4)*u2x*u3z + cm(24, 4)*u2y*u3z + cm(27, 4)*u2z*u3z
        m(5, 3) = cm(19, 5)*u2x*u3x + cm(22, 5)*u2y*u3x + cm(25, 5)*u2z*u3x + cm(20, 5)*u2x*u3y + cm(23, 5)*u2y*u3y + cm(26, 5)*u2z*u3y + cm(21, 5)*u2x*u3z + cm(24, 5)*u2y*u3z + cm(27, 5)*u2z*u3z
        m(6, 3) = cm(19, 6)*u2x*u3x + cm(22, 6)*u2y*u3x + cm(25, 6)*u2z*u3x + cm(20, 6)*u2x*u3y + cm(23, 6)*u2y*u3y + cm(26, 6)*u2z*u3y + cm(21, 6)*u2x*u3z + cm(24, 6)*u2y*u3z + cm(27, 6)*u2z*u3z
        m(7, 3) = cm(19, 7)*u2x*u3x + cm(22, 7)*u2y*u3x + cm(25, 7)*u2z*u3x + cm(20, 7)*u2x*u3y + cm(23, 7)*u2y*u3y + cm(26, 7)*u2z*u3y + cm(21, 7)*u2x*u3z + cm(24, 7)*u2y*u3z + cm(27, 7)*u2z*u3z
        m(8, 3) = cm(19, 8)*u2x*u3x + cm(22, 8)*u2y*u3x + cm(25, 8)*u2z*u3x + cm(20, 8)*u2x*u3y + cm(23, 8)*u2y*u3y + cm(26, 8)*u2z*u3y + cm(21, 8)*u2x*u3z + cm(24, 8)*u2y*u3z + cm(27, 8)*u2z*u3z
        m(9, 3) = cm(19, 9)*u2x*u3x + cm(22, 9)*u2y*u3x + cm(25, 9)*u2z*u3x + cm(20, 9)*u2x*u3y + cm(23, 9)*u2y*u3y + cm(26, 9)*u2z*u3y + cm(21, 9)*u2x*u3z + cm(24, 9)*u2y*u3z + cm(27, 9)*u2z*u3z
        m(10, 3) = cm(19, 10)*u2x*u3x + cm(22, 10)*u2y*u3x + cm(25, 10)*u2z*u3x + cm(20, 10)*u2x*u3y + cm(23, 10)*u2y*u3y + cm(26, 10)*u2z*u3y + cm(21, 10)*u2x*u3z + cm(24, 10)*u2y*u3z + cm(27, 10)*u2z*u3z
        m(11, 3) = cm(19, 11)*u2x*u3x + cm(22, 11)*u2y*u3x + cm(25, 11)*u2z*u3x + cm(20, 11)*u2x*u3y + cm(23, 11)*u2y*u3y + cm(26, 11)*u2z*u3y + cm(21, 11)*u2x*u3z + cm(24, 11)*u2y*u3z + cm(27, 11)*u2z*u3z
        m(12, 3) = cm(19, 12)*u2x*u3x + cm(22, 12)*u2y*u3x + cm(25, 12)*u2z*u3x + cm(20, 12)*u2x*u3y + cm(23, 12)*u2y*u3y + cm(26, 12)*u2z*u3y + cm(21, 12)*u2x*u3z + cm(24, 12)*u2y*u3z + cm(27, 12)*u2z*u3z
        m(13, 3) = cm(19, 13)*u2x*u3x + cm(22, 13)*u2y*u3x + cm(25, 13)*u2z*u3x + cm(20, 13)*u2x*u3y + cm(23, 13)*u2y*u3y + cm(26, 13)*u2z*u3y + cm(21, 13)*u2x*u3z + cm(24, 13)*u2y*u3z + cm(27, 13)*u2z*u3z
        m(14, 3) = cm(19, 14)*u2x*u3x + cm(22, 14)*u2y*u3x + cm(25, 14)*u2z*u3x + cm(20, 14)*u2x*u3y + cm(23, 14)*u2y*u3y + cm(26, 14)*u2z*u3y + cm(21, 14)*u2x*u3z + cm(24, 14)*u2y*u3z + cm(27, 14)*u2z*u3z
        m(15, 3) = cm(19, 15)*u2x*u3x + cm(22, 15)*u2y*u3x + cm(25, 15)*u2z*u3x + cm(20, 15)*u2x*u3y + cm(23, 15)*u2y*u3y + cm(26, 15)*u2z*u3y + cm(21, 15)*u2x*u3z + cm(24, 15)*u2y*u3z + cm(27, 15)*u2z*u3z
        m(16, 3) = cm(19, 16)*u2x*u3x + cm(22, 16)*u2y*u3x + cm(25, 16)*u2z*u3x + cm(20, 16)*u2x*u3y + cm(23, 16)*u2y*u3y + cm(26, 16)*u2z*u3y + cm(21, 16)*u2x*u3z + cm(24, 16)*u2y*u3z + cm(27, 16)*u2z*u3z
        m(17, 3) = cm(19, 17)*u2x*u3x + cm(22, 17)*u2y*u3x + cm(25, 17)*u2z*u3x + cm(20, 17)*u2x*u3y + cm(23, 17)*u2y*u3y + cm(26, 17)*u2z*u3y + cm(21, 17)*u2x*u3z + cm(24, 17)*u2y*u3z + cm(27, 17)*u2z*u3z
        m(18, 3) = cm(19, 18)*u2x*u3x + cm(22, 18)*u2y*u3x + cm(25, 18)*u2z*u3x + cm(20, 18)*u2x*u3y + cm(23, 18)*u2y*u3y + cm(26, 18)*u2z*u3y + cm(21, 18)*u2x*u3z + cm(24, 18)*u2y*u3z + cm(27, 18)*u2z*u3z
        m(19, 3) = cm(19, 19)*u2x*u3x + cm(22, 19)*u2y*u3x + cm(25, 19)*u2z*u3x + cm(20, 19)*u2x*u3y + cm(23, 19)*u2y*u3y + cm(26, 19)*u2z*u3y + cm(21, 19)*u2x*u3z + cm(24, 19)*u2y*u3z + cm(27, 19)*u2z*u3z
        m(20, 3) = cm(19, 20)*u2x*u3x + cm(22, 20)*u2y*u3x + cm(25, 20)*u2z*u3x + cm(20, 20)*u2x*u3y + cm(23, 20)*u2y*u3y + cm(26, 20)*u2z*u3y + cm(21, 20)*u2x*u3z + cm(24, 20)*u2y*u3z + cm(27, 20)*u2z*u3z
        m(21, 3) = cm(19, 21)*u2x*u3x + cm(22, 21)*u2y*u3x + cm(25, 21)*u2z*u3x + cm(20, 21)*u2x*u3y + cm(23, 21)*u2y*u3y + cm(26, 21)*u2z*u3y + cm(21, 21)*u2x*u3z + cm(24, 21)*u2y*u3z + cm(27, 21)*u2z*u3z
        m(22, 3) = cm(19, 22)*u2x*u3x + cm(22, 22)*u2y*u3x + cm(25, 22)*u2z*u3x + cm(20, 22)*u2x*u3y + cm(23, 22)*u2y*u3y + cm(26, 22)*u2z*u3y + cm(21, 22)*u2x*u3z + cm(24, 22)*u2y*u3z + cm(27, 22)*u2z*u3z
        m(23, 3) = cm(19, 23)*u2x*u3x + cm(22, 23)*u2y*u3x + cm(25, 23)*u2z*u3x + cm(20, 23)*u2x*u3y + cm(23, 23)*u2y*u3y + cm(26, 23)*u2z*u3y + cm(21, 23)*u2x*u3z + cm(24, 23)*u2y*u3z + cm(27, 23)*u2z*u3z
        m(24, 3) = cm(19, 24)*u2x*u3x + cm(22, 24)*u2y*u3x + cm(25, 24)*u2z*u3x + cm(20, 24)*u2x*u3y + cm(23, 24)*u2y*u3y + cm(26, 24)*u2z*u3y + cm(21, 24)*u2x*u3z + cm(24, 24)*u2y*u3z + cm(27, 24)*u2z*u3z
        m(25, 3) = cm(19, 25)*u2x*u3x + cm(22, 25)*u2y*u3x + cm(25, 25)*u2z*u3x + cm(20, 25)*u2x*u3y + cm(23, 25)*u2y*u3y + cm(26, 25)*u2z*u3y + cm(21, 25)*u2x*u3z + cm(24, 25)*u2y*u3z + cm(27, 25)*u2z*u3z
        m(26, 3) = cm(19, 26)*u2x*u3x + cm(22, 26)*u2y*u3x + cm(25, 26)*u2z*u3x + cm(20, 26)*u2x*u3y + cm(23, 26)*u2y*u3y + cm(26, 26)*u2z*u3y + cm(21, 26)*u2x*u3z + cm(24, 26)*u2y*u3z + cm(27, 26)*u2z*u3z
        m(27, 3) = cm(19, 27)*u2x*u3x + cm(22, 27)*u2y*u3x + cm(25, 27)*u2z*u3x + cm(20, 27)*u2x*u3y + cm(23, 27)*u2y*u3y + cm(26, 27)*u2z*u3y + cm(21, 27)*u2x*u3z + cm(24, 27)*u2y*u3z + cm(27, 27)*u2z*u3z
        m = m*0.5_r8
    end function

end subroutine

end submodule
