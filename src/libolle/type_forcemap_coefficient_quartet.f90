#include "precompilerdefinitions"
submodule(type_forcemap) type_forcemap_coefficient_quartet
implicit none
contains

!> quartet coefficient matrix
module subroutine lo_coeffmatrix_quartet(UM, CM, map, relntheta)
    !> the displacements
    real(r8), dimension(:, :), intent(in) :: UM
    !> the coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM
    !> symmetry stuff rearranged into coordination shells
    type(lo_forcemap), intent(in) :: map
    !> which arrays will be relevant
    logical, dimension(81), intent(in) :: relntheta

    integer :: a1, a2, a3, a4, i, k, ii
    integer :: find, uind, ntheta, unop, unsh
    integer, dimension(81) :: thetaind
    real(r8), dimension(3) :: u2, u3, u4, v2, v3, v4
    ! had to make these allocatable since I ran into lots of stack issues
    real(r8), dimension(:, :), allocatable :: &
        wm81x1, wm81x2, wm81x3, wm81x4, wm81x5, wm81x6, wm81x7, wm81x8, &
        wm81x9, wm81x10, wm81x11, wm81x12, wm81x13, wm81x14, wm81x15, &
        wm81x16, wm81x17, wm81x18, wm81x19, wm81x20, wm81x21, wm81x22, &
        wm81x23, wm81x24, wm81x25, wm81x26, wm81x27, wm81x28, wm81x29, &
        wm81x30, wm81x31, wm81x32, wm81x33, wm81x34, wm81x35, wm81x36, &
        wm81x37, wm81x38, wm81x39, wm81x40, wm81x41, wm81x42, wm81x43, &
        wm81x44, wm81x45, wm81x46, wm81x47, wm81x48, wm81x49, wm81x50, &
        wm81x51, wm81x52, wm81x53, wm81x54, wm81x55, wm81x56, wm81x57, &
        wm81x58, wm81x59, wm81x60, wm81x61, wm81x62, wm81x63, wm81x64, &
        wm81x65, wm81x66, wm81x67, wm81x68, wm81x69, wm81x70, wm81x71, &
        wm81x72, wm81x73, wm81x74, wm81x75, wm81x76, wm81x77, wm81x78, &
        wm81x79, wm81x80, wm81x81

    real(r8), dimension(:, :), allocatable ::  cm1x3, cm2x3, cm3x3, cm4x3, &
                                              cm5x3, cm6x3, cm7x3, cm8x3, cm9x3, &
                                              cm10x3, cm11x3, cm12x3, cm13x3, cm14x3, &
                                              cm15x3, cm16x3, cm17x3, cm18x3, cm19x3, &
                                              cm20x3, cm21x3, cm22x3, cm23x3, cm24x3, &
                                              cm25x3, cm26x3, cm27x3, cm28x3, cm29x3, &
                                              cm30x3, cm31x3, cm32x3, cm33x3, cm34x3, &
                                              cm35x3, cm36x3, cm37x3, cm38x3, cm39x3, &
                                              cm40x3, cm41x3, cm42x3, cm43x3, cm44x3, &
                                              cm45x3, cm46x3, cm47x3, cm48x3, cm49x3, &
                                              cm50x3, cm51x3, cm52x3, cm53x3, cm54x3, &
                                              cm55x3, cm56x3, cm57x3, cm58x3, cm59x3, &
                                              cm60x3, cm61x3, cm62x3, cm63x3, cm64x3, &
                                              cm65x3, cm66x3, cm67x3, cm68x3, cm69x3, &
                                              cm70x3, cm71x3, cm72x3, cm73x3, cm74x3, &
                                              cm75x3, cm76x3, cm77x3, cm78x3, cm79x3, &
                                              cm80x3, cm81x3

    ! Only allocate those that will actually be used, usually quite few
    if (relntheta(1)) allocate (wm81x1(81, 1))
    if (relntheta(2)) allocate (wm81x2(81, 2))
    if (relntheta(3)) allocate (wm81x3(81, 3))
    if (relntheta(4)) allocate (wm81x4(81, 4))
    if (relntheta(5)) allocate (wm81x5(81, 5))
    if (relntheta(6)) allocate (wm81x6(81, 6))
    if (relntheta(7)) allocate (wm81x7(81, 7))
    if (relntheta(8)) allocate (wm81x8(81, 8))
    if (relntheta(9)) allocate (wm81x9(81, 9))
    if (relntheta(10)) allocate (wm81x10(81, 10))
    if (relntheta(11)) allocate (wm81x11(81, 11))
    if (relntheta(12)) allocate (wm81x12(81, 12))
    if (relntheta(13)) allocate (wm81x13(81, 13))
    if (relntheta(14)) allocate (wm81x14(81, 14))
    if (relntheta(15)) allocate (wm81x15(81, 15))
    if (relntheta(16)) allocate (wm81x16(81, 16))
    if (relntheta(17)) allocate (wm81x17(81, 17))
    if (relntheta(18)) allocate (wm81x18(81, 18))
    if (relntheta(19)) allocate (wm81x19(81, 19))
    if (relntheta(20)) allocate (wm81x20(81, 20))
    if (relntheta(21)) allocate (wm81x21(81, 21))
    if (relntheta(22)) allocate (wm81x22(81, 22))
    if (relntheta(23)) allocate (wm81x23(81, 23))
    if (relntheta(24)) allocate (wm81x24(81, 24))
    if (relntheta(25)) allocate (wm81x25(81, 25))
    if (relntheta(26)) allocate (wm81x26(81, 26))
    if (relntheta(27)) allocate (wm81x27(81, 27))
    if (relntheta(28)) allocate (wm81x28(81, 28))
    if (relntheta(29)) allocate (wm81x29(81, 29))
    if (relntheta(30)) allocate (wm81x30(81, 30))
    if (relntheta(31)) allocate (wm81x31(81, 31))
    if (relntheta(32)) allocate (wm81x32(81, 32))
    if (relntheta(33)) allocate (wm81x33(81, 33))
    if (relntheta(34)) allocate (wm81x34(81, 34))
    if (relntheta(35)) allocate (wm81x35(81, 35))
    if (relntheta(36)) allocate (wm81x36(81, 36))
    if (relntheta(37)) allocate (wm81x37(81, 37))
    if (relntheta(38)) allocate (wm81x38(81, 38))
    if (relntheta(39)) allocate (wm81x39(81, 39))
    if (relntheta(40)) allocate (wm81x40(81, 40))
    if (relntheta(41)) allocate (wm81x41(81, 41))
    if (relntheta(42)) allocate (wm81x42(81, 42))
    if (relntheta(43)) allocate (wm81x43(81, 43))
    if (relntheta(44)) allocate (wm81x44(81, 44))
    if (relntheta(45)) allocate (wm81x45(81, 45))
    if (relntheta(46)) allocate (wm81x46(81, 46))
    if (relntheta(47)) allocate (wm81x47(81, 47))
    if (relntheta(48)) allocate (wm81x48(81, 48))
    if (relntheta(49)) allocate (wm81x49(81, 49))
    if (relntheta(50)) allocate (wm81x50(81, 50))
    if (relntheta(51)) allocate (wm81x51(81, 51))
    if (relntheta(52)) allocate (wm81x52(81, 52))
    if (relntheta(53)) allocate (wm81x53(81, 53))
    if (relntheta(54)) allocate (wm81x54(81, 54))
    if (relntheta(55)) allocate (wm81x55(81, 55))
    if (relntheta(56)) allocate (wm81x56(81, 56))
    if (relntheta(57)) allocate (wm81x57(81, 57))
    if (relntheta(58)) allocate (wm81x58(81, 58))
    if (relntheta(59)) allocate (wm81x59(81, 59))
    if (relntheta(60)) allocate (wm81x60(81, 60))
    if (relntheta(61)) allocate (wm81x61(81, 61))
    if (relntheta(62)) allocate (wm81x62(81, 62))
    if (relntheta(63)) allocate (wm81x63(81, 63))
    if (relntheta(64)) allocate (wm81x64(81, 64))
    if (relntheta(65)) allocate (wm81x65(81, 65))
    if (relntheta(66)) allocate (wm81x66(81, 66))
    if (relntheta(67)) allocate (wm81x67(81, 67))
    if (relntheta(68)) allocate (wm81x68(81, 68))
    if (relntheta(69)) allocate (wm81x69(81, 69))
    if (relntheta(70)) allocate (wm81x70(81, 70))
    if (relntheta(71)) allocate (wm81x71(81, 71))
    if (relntheta(72)) allocate (wm81x72(81, 72))
    if (relntheta(73)) allocate (wm81x73(81, 73))
    if (relntheta(74)) allocate (wm81x74(81, 74))
    if (relntheta(75)) allocate (wm81x75(81, 75))
    if (relntheta(76)) allocate (wm81x76(81, 76))
    if (relntheta(77)) allocate (wm81x77(81, 77))
    if (relntheta(78)) allocate (wm81x78(81, 78))
    if (relntheta(79)) allocate (wm81x79(81, 79))
    if (relntheta(80)) allocate (wm81x80(81, 80))
    if (relntheta(81)) allocate (wm81x81(81, 81))

    if (relntheta(1)) allocate (cm1x3(1, 3))
    if (relntheta(2)) allocate (cm2x3(2, 3))
    if (relntheta(3)) allocate (cm3x3(3, 3))
    if (relntheta(4)) allocate (cm4x3(4, 3))
    if (relntheta(5)) allocate (cm5x3(5, 3))
    if (relntheta(6)) allocate (cm6x3(6, 3))
    if (relntheta(7)) allocate (cm7x3(7, 3))
    if (relntheta(8)) allocate (cm8x3(8, 3))
    if (relntheta(9)) allocate (cm9x3(9, 3))
    if (relntheta(10)) allocate (cm10x3(10, 3))
    if (relntheta(11)) allocate (cm11x3(11, 3))
    if (relntheta(12)) allocate (cm12x3(12, 3))
    if (relntheta(13)) allocate (cm13x3(13, 3))
    if (relntheta(14)) allocate (cm14x3(14, 3))
    if (relntheta(15)) allocate (cm15x3(15, 3))
    if (relntheta(16)) allocate (cm16x3(16, 3))
    if (relntheta(17)) allocate (cm17x3(17, 3))
    if (relntheta(18)) allocate (cm18x3(18, 3))
    if (relntheta(19)) allocate (cm19x3(19, 3))
    if (relntheta(20)) allocate (cm20x3(20, 3))
    if (relntheta(21)) allocate (cm21x3(21, 3))
    if (relntheta(22)) allocate (cm22x3(22, 3))
    if (relntheta(23)) allocate (cm23x3(23, 3))
    if (relntheta(24)) allocate (cm24x3(24, 3))
    if (relntheta(25)) allocate (cm25x3(25, 3))
    if (relntheta(26)) allocate (cm26x3(26, 3))
    if (relntheta(27)) allocate (cm27x3(27, 3))
    if (relntheta(28)) allocate (cm28x3(28, 3))
    if (relntheta(29)) allocate (cm29x3(29, 3))
    if (relntheta(30)) allocate (cm30x3(30, 3))
    if (relntheta(31)) allocate (cm31x3(31, 3))
    if (relntheta(32)) allocate (cm32x3(32, 3))
    if (relntheta(33)) allocate (cm33x3(33, 3))
    if (relntheta(34)) allocate (cm34x3(34, 3))
    if (relntheta(35)) allocate (cm35x3(35, 3))
    if (relntheta(36)) allocate (cm36x3(36, 3))
    if (relntheta(37)) allocate (cm37x3(37, 3))
    if (relntheta(38)) allocate (cm38x3(38, 3))
    if (relntheta(39)) allocate (cm39x3(39, 3))
    if (relntheta(40)) allocate (cm40x3(40, 3))
    if (relntheta(41)) allocate (cm41x3(41, 3))
    if (relntheta(42)) allocate (cm42x3(42, 3))
    if (relntheta(43)) allocate (cm43x3(43, 3))
    if (relntheta(44)) allocate (cm44x3(44, 3))
    if (relntheta(45)) allocate (cm45x3(45, 3))
    if (relntheta(46)) allocate (cm46x3(46, 3))
    if (relntheta(47)) allocate (cm47x3(47, 3))
    if (relntheta(48)) allocate (cm48x3(48, 3))
    if (relntheta(49)) allocate (cm49x3(49, 3))
    if (relntheta(50)) allocate (cm50x3(50, 3))
    if (relntheta(51)) allocate (cm51x3(51, 3))
    if (relntheta(52)) allocate (cm52x3(52, 3))
    if (relntheta(53)) allocate (cm53x3(53, 3))
    if (relntheta(54)) allocate (cm54x3(54, 3))
    if (relntheta(55)) allocate (cm55x3(55, 3))
    if (relntheta(56)) allocate (cm56x3(56, 3))
    if (relntheta(57)) allocate (cm57x3(57, 3))
    if (relntheta(58)) allocate (cm58x3(58, 3))
    if (relntheta(59)) allocate (cm59x3(59, 3))
    if (relntheta(60)) allocate (cm60x3(60, 3))
    if (relntheta(61)) allocate (cm61x3(61, 3))
    if (relntheta(62)) allocate (cm62x3(62, 3))
    if (relntheta(63)) allocate (cm63x3(63, 3))
    if (relntheta(64)) allocate (cm64x3(64, 3))
    if (relntheta(65)) allocate (cm65x3(65, 3))
    if (relntheta(66)) allocate (cm66x3(66, 3))
    if (relntheta(67)) allocate (cm67x3(67, 3))
    if (relntheta(68)) allocate (cm68x3(68, 3))
    if (relntheta(69)) allocate (cm69x3(69, 3))
    if (relntheta(70)) allocate (cm70x3(70, 3))
    if (relntheta(71)) allocate (cm71x3(71, 3))
    if (relntheta(72)) allocate (cm72x3(72, 3))
    if (relntheta(73)) allocate (cm73x3(73, 3))
    if (relntheta(74)) allocate (cm74x3(74, 3))
    if (relntheta(75)) allocate (cm75x3(75, 3))
    if (relntheta(76)) allocate (cm76x3(76, 3))
    if (relntheta(77)) allocate (cm77x3(77, 3))
    if (relntheta(78)) allocate (cm78x3(78, 3))
    if (relntheta(79)) allocate (cm79x3(79, 3))
    if (relntheta(80)) allocate (cm80x3(80, 3))
    if (relntheta(81)) allocate (cm81x3(81, 3))

    CM = 0.0_r8
    do i = 1, map%xss%n_fc_quartet
        unsh = map%xss%ind_fc_quartet(1, i)
        unop = map%xss%ind_fc_quartet(2, i)
        a1 = map%xss%ind_fc_quartet(3, i)
        a2 = map%xss%ind_fc_quartet(4, i)
        a3 = map%xss%ind_fc_quartet(5, i)
        a4 = map%xss%ind_fc_quartet(6, i)
        ntheta = map%fc_quartet_shell(unsh)%nx
        if (ntheta .eq. 0) cycle
        v2 = UM(:, a1)
        v3 = UM(:, a1)
        v4 = UM(:, a1)
        u2 = UM(:, a2)
        u3 = UM(:, a3)
        u4 = UM(:, a4)
        ntheta = map%fc_quartet_shell(unsh)%nx
        thetaind(1:ntheta) = map%fc_quartet_shell(unsh)%ind_global
        select case (ntheta)
        case (1)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x1)
            cm1x3 = coeff81x1(wm81x1, u2, u3, u4) - coeff81x1(wm81x1, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm1x3(k, ii)
                end do; end do
        case (2)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x2)
            cm2x3 = coeff81x2(wm81x2, u2, u3, u4) - coeff81x2(wm81x2, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm2x3(k, ii)
                end do; end do
        case (3)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x3)
            cm3x3 = coeff81x3(wm81x3, u2, u3, u4) - coeff81x3(wm81x3, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm3x3(k, ii)
                end do; end do
        case (4)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x4)
            cm4x3 = coeff81x4(wm81x4, u2, u3, u4) - coeff81x4(wm81x4, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm4x3(k, ii)
                end do; end do
        case (5)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x5)
            cm5x3 = coeff81x5(wm81x5, u2, u3, u4) - coeff81x5(wm81x5, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm5x3(k, ii)
                end do; end do
        case (6)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x6)
            cm6x3 = coeff81x6(wm81x6, u2, u3, u4) - coeff81x6(wm81x6, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm6x3(k, ii)
                end do; end do
        case (7)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x7)
            cm7x3 = coeff81x7(wm81x7, u2, u3, u4) - coeff81x7(wm81x7, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm7x3(k, ii)
                end do; end do
        case (8)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x8)
            cm8x3 = coeff81x8(wm81x8, u2, u3, u4) - coeff81x8(wm81x8, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm8x3(k, ii)
                end do; end do
        case (9)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x9)
            cm9x3 = coeff81x9(wm81x9, u2, u3, u4) - coeff81x9(wm81x9, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm9x3(k, ii)
                end do; end do
        case (10)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x10)
            cm10x3 = coeff81x10(wm81x10, u2, u3, u4) - coeff81x10(wm81x10, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm10x3(k, ii)
                end do; end do
        case (11)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x11)
            cm11x3 = coeff81x11(wm81x11, u2, u3, u4) - coeff81x11(wm81x11, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm11x3(k, ii)
                end do; end do
        case (12)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x12)
            cm12x3 = coeff81x12(wm81x12, u2, u3, u4) - coeff81x12(wm81x12, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm12x3(k, ii)
                end do; end do
        case (13)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x13)
            cm13x3 = coeff81x13(wm81x13, u2, u3, u4) - coeff81x13(wm81x13, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm13x3(k, ii)
                end do; end do
        case (14)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x14)
            cm14x3 = coeff81x14(wm81x14, u2, u3, u4) - coeff81x14(wm81x14, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm14x3(k, ii)
                end do; end do
        case (15)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x15)
            cm15x3 = coeff81x15(wm81x15, u2, u3, u4) - coeff81x15(wm81x15, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm15x3(k, ii)
                end do; end do
        case (16)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x16)
            cm16x3 = coeff81x16(wm81x16, u2, u3, u4) - coeff81x16(wm81x16, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm16x3(k, ii)
                end do; end do
        case (17)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x17)
            cm17x3 = coeff81x17(wm81x17, u2, u3, u4) - coeff81x17(wm81x17, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm17x3(k, ii)
                end do; end do
        case (18)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x18)
            cm18x3 = coeff81x18(wm81x18, u2, u3, u4) - coeff81x18(wm81x18, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm18x3(k, ii)
                end do; end do
        case (19)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x19)
            cm19x3 = coeff81x19(wm81x19, u2, u3, u4) - coeff81x19(wm81x19, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm19x3(k, ii)
                end do; end do
        case (20)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x20)
            cm20x3 = coeff81x20(wm81x20, u2, u3, u4) - coeff81x20(wm81x20, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm20x3(k, ii)
                end do; end do
        case (21)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x21)
            cm21x3 = coeff81x21(wm81x21, u2, u3, u4) - coeff81x21(wm81x21, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm21x3(k, ii)
                end do; end do
        case (22)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x22)
            cm22x3 = coeff81x22(wm81x22, u2, u3, u4) - coeff81x22(wm81x22, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm22x3(k, ii)
                end do; end do
        case (23)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x23)
            cm23x3 = coeff81x23(wm81x23, u2, u3, u4) - coeff81x23(wm81x23, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm23x3(k, ii)
                end do; end do
        case (24)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x24)
            cm24x3 = coeff81x24(wm81x24, u2, u3, u4) - coeff81x24(wm81x24, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm24x3(k, ii)
                end do; end do
        case (25)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x25)
            cm25x3 = coeff81x25(wm81x25, u2, u3, u4) - coeff81x25(wm81x25, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm25x3(k, ii)
                end do; end do
        case (26)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x26)
            cm26x3 = coeff81x26(wm81x26, u2, u3, u4) - coeff81x26(wm81x26, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm26x3(k, ii)
                end do; end do
        case (27)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x27)
            cm27x3 = coeff81x27(wm81x27, u2, u3, u4) - coeff81x27(wm81x27, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm27x3(k, ii)
                end do; end do
        case (28)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x28)
            cm28x3 = coeff81x28(wm81x28, u2, u3, u4) - coeff81x28(wm81x28, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm28x3(k, ii)
                end do; end do
        case (29)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x29)
            cm29x3 = coeff81x29(wm81x29, u2, u3, u4) - coeff81x29(wm81x29, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm29x3(k, ii)
                end do; end do
        case (30)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x30)
            cm30x3 = coeff81x30(wm81x30, u2, u3, u4) - coeff81x30(wm81x30, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm30x3(k, ii)
                end do; end do
        case (31)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x31)
            cm31x3 = coeff81x31(wm81x31, u2, u3, u4) - coeff81x31(wm81x31, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm31x3(k, ii)
                end do; end do
        case (32)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x32)
            cm32x3 = coeff81x32(wm81x32, u2, u3, u4) - coeff81x32(wm81x32, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm32x3(k, ii)
                end do; end do
        case (33)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x33)
            cm33x3 = coeff81x33(wm81x33, u2, u3, u4) - coeff81x33(wm81x33, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm33x3(k, ii)
                end do; end do
        case (34)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x34)
            cm34x3 = coeff81x34(wm81x34, u2, u3, u4) - coeff81x34(wm81x34, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm34x3(k, ii)
                end do; end do
        case (35)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x35)
            cm35x3 = coeff81x35(wm81x35, u2, u3, u4) - coeff81x35(wm81x35, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm35x3(k, ii)
                end do; end do
        case (36)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x36)
            cm36x3 = coeff81x36(wm81x36, u2, u3, u4) - coeff81x36(wm81x36, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm36x3(k, ii)
                end do; end do
        case (37)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x37)
            cm37x3 = coeff81x37(wm81x37, u2, u3, u4) - coeff81x37(wm81x37, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm37x3(k, ii)
                end do; end do
        case (38)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x38)
            cm38x3 = coeff81x38(wm81x38, u2, u3, u4) - coeff81x38(wm81x38, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm38x3(k, ii)
                end do; end do
        case (39)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x39)
            cm39x3 = coeff81x39(wm81x39, u2, u3, u4) - coeff81x39(wm81x39, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm39x3(k, ii)
                end do; end do
        case (40)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x40)
            cm40x3 = coeff81x40(wm81x40, u2, u3, u4) - coeff81x40(wm81x40, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm40x3(k, ii)
                end do; end do
        case (41)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x41)
            cm41x3 = coeff81x41(wm81x41, u2, u3, u4) - coeff81x41(wm81x41, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm41x3(k, ii)
                end do; end do
        case (42)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x42)
            cm42x3 = coeff81x42(wm81x42, u2, u3, u4) - coeff81x42(wm81x42, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm42x3(k, ii)
                end do; end do
        case (43)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x43)
            cm43x3 = coeff81x43(wm81x43, u2, u3, u4) - coeff81x43(wm81x43, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm43x3(k, ii)
                end do; end do
        case (44)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x44)
            cm44x3 = coeff81x44(wm81x44, u2, u3, u4) - coeff81x44(wm81x44, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm44x3(k, ii)
                end do; end do
        case (45)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x45)
            cm45x3 = coeff81x45(wm81x45, u2, u3, u4) - coeff81x45(wm81x45, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm45x3(k, ii)
                end do; end do
        case (46)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x46)
            cm46x3 = coeff81x46(wm81x46, u2, u3, u4) - coeff81x46(wm81x46, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm46x3(k, ii)
                end do; end do
        case (47)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x47)
            cm47x3 = coeff81x47(wm81x47, u2, u3, u4) - coeff81x47(wm81x47, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm47x3(k, ii)
                end do; end do
        case (48)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x48)
            cm48x3 = coeff81x48(wm81x48, u2, u3, u4) - coeff81x48(wm81x48, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm48x3(k, ii)
                end do; end do
        case (49)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x49)
            cm49x3 = coeff81x49(wm81x49, u2, u3, u4) - coeff81x49(wm81x49, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm49x3(k, ii)
                end do; end do
        case (50)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x50)
            cm50x3 = coeff81x50(wm81x50, u2, u3, u4) - coeff81x50(wm81x50, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm50x3(k, ii)
                end do; end do
        case (51)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x51)
            cm51x3 = coeff81x51(wm81x51, u2, u3, u4) - coeff81x51(wm81x51, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm51x3(k, ii)
                end do; end do
        case (52)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x52)
            cm52x3 = coeff81x52(wm81x52, u2, u3, u4) - coeff81x52(wm81x52, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm52x3(k, ii)
                end do; end do
        case (53)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x53)
            cm53x3 = coeff81x53(wm81x53, u2, u3, u4) - coeff81x53(wm81x53, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm53x3(k, ii)
                end do; end do
        case (54)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x54)
            cm54x3 = coeff81x54(wm81x54, u2, u3, u4) - coeff81x54(wm81x54, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm54x3(k, ii)
                end do; end do
        case (55)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x55)
            cm55x3 = coeff81x55(wm81x55, u2, u3, u4) - coeff81x55(wm81x55, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm55x3(k, ii)
                end do; end do
        case (56)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x56)
            cm56x3 = coeff81x56(wm81x56, u2, u3, u4) - coeff81x56(wm81x56, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm56x3(k, ii)
                end do; end do
        case (57)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x57)
            cm57x3 = coeff81x57(wm81x57, u2, u3, u4) - coeff81x57(wm81x57, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm57x3(k, ii)
                end do; end do
        case (58)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x58)
            cm58x3 = coeff81x58(wm81x58, u2, u3, u4) - coeff81x58(wm81x58, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm58x3(k, ii)
                end do; end do
        case (59)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x59)
            cm59x3 = coeff81x59(wm81x59, u2, u3, u4) - coeff81x59(wm81x59, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm59x3(k, ii)
                end do; end do
        case (60)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x60)
            cm60x3 = coeff81x60(wm81x60, u2, u3, u4) - coeff81x60(wm81x60, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm60x3(k, ii)
                end do; end do
        case (61)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x61)
            cm61x3 = coeff81x61(wm81x61, u2, u3, u4) - coeff81x61(wm81x61, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm61x3(k, ii)
                end do; end do
        case (62)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x62)
            cm62x3 = coeff81x62(wm81x62, u2, u3, u4) - coeff81x62(wm81x62, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm62x3(k, ii)
                end do; end do
        case (63)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x63)
            cm63x3 = coeff81x63(wm81x63, u2, u3, u4) - coeff81x63(wm81x63, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm63x3(k, ii)
                end do; end do
        case (64)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x64)
            cm64x3 = coeff81x64(wm81x64, u2, u3, u4) - coeff81x64(wm81x64, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm64x3(k, ii)
                end do; end do
        case (65)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x65)
            cm65x3 = coeff81x65(wm81x65, u2, u3, u4) - coeff81x65(wm81x65, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm65x3(k, ii)
                end do; end do
        case (66)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x66)
            cm66x3 = coeff81x66(wm81x66, u2, u3, u4) - coeff81x66(wm81x66, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm66x3(k, ii)
                end do; end do
        case (67)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x67)
            cm67x3 = coeff81x67(wm81x67, u2, u3, u4) - coeff81x67(wm81x67, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm67x3(k, ii)
                end do; end do
        case (68)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x68)
            cm68x3 = coeff81x68(wm81x68, u2, u3, u4) - coeff81x68(wm81x68, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm68x3(k, ii)
                end do; end do
        case (69)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x69)
            cm69x3 = coeff81x69(wm81x69, u2, u3, u4) - coeff81x69(wm81x69, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm69x3(k, ii)
                end do; end do
        case (70)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x70)
            cm70x3 = coeff81x70(wm81x70, u2, u3, u4) - coeff81x70(wm81x70, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm70x3(k, ii)
                end do; end do
        case (71)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x71)
            cm71x3 = coeff81x71(wm81x71, u2, u3, u4) - coeff81x71(wm81x71, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm71x3(k, ii)
                end do; end do
        case (72)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x72)
            cm72x3 = coeff81x72(wm81x72, u2, u3, u4) - coeff81x72(wm81x72, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm72x3(k, ii)
                end do; end do
        case (73)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x73)
            cm73x3 = coeff81x73(wm81x73, u2, u3, u4) - coeff81x73(wm81x73, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm73x3(k, ii)
                end do; end do
        case (74)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x74)
            cm74x3 = coeff81x74(wm81x74, u2, u3, u4) - coeff81x74(wm81x74, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm74x3(k, ii)
                end do; end do
        case (75)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x75)
            cm75x3 = coeff81x75(wm81x75, u2, u3, u4) - coeff81x75(wm81x75, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm75x3(k, ii)
                end do; end do
        case (76)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x76)
            cm76x3 = coeff81x76(wm81x76, u2, u3, u4) - coeff81x76(wm81x76, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm76x3(k, ii)
                end do; end do
        case (77)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x77)
            cm77x3 = coeff81x77(wm81x77, u2, u3, u4) - coeff81x77(wm81x77, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm77x3(k, ii)
                end do; end do
        case (78)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x78)
            cm78x3 = coeff81x78(wm81x78, u2, u3, u4) - coeff81x78(wm81x78, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm78x3(k, ii)
                end do; end do
        case (79)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x79)
            cm79x3 = coeff81x79(wm81x79, u2, u3, u4) - coeff81x79(wm81x79, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm79x3(k, ii)
                end do; end do
        case (80)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x80)
            cm80x3 = coeff81x80(wm81x80, u2, u3, u4) - coeff81x80(wm81x80, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm80x3(k, ii)
                end do; end do
        case (81)
            call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x81)
            cm81x3 = coeff81x81(wm81x81, u2, u3, u4) - coeff81x81(wm81x81, v2, v3, v4)
            do k = 1, ntheta; do ii = 1, 3
                    uind = thetaind(k); find = (a1 - 1)*3 + ii
                    CM(find, uind) = CM(find, uind) + cm81x3(k, ii)
                end do; end do
        end select
    end do

end subroutine

end submodule
