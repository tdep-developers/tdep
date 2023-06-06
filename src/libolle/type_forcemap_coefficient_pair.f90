#include "precompilerdefinitions"
submodule(type_forcemap) type_forcemap_coefficient_pair
implicit none
contains

!> get the singlet coefficient matrix
module subroutine lo_coeffmatrix_singlet(CM, map)
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM
    !> symmetry stuff rearranged into coordination shells
    type(lo_forcemap), intent(in) :: map

    integer :: a1, k, ii
    integer :: find, uind, unsh, unop, nx
    real(r8), dimension(:, :), allocatable :: m0

    allocate (m0(3, map%xuc%nx_fc_singlet))
    m0 = 0.0_r8

    CM = 0.0_r8
    do a1 = 1, map%n_atom_ss
        ! coefficient matrix
        unsh = map%xss%ind_fc_singlet(1, a1)
        unop = map%xss%ind_fc_singlet(2, a1)
        nx = map%fc_singlet_shell(unsh)%nx
        if (nx .eq. 0) cycle
        ! transform to this atom
        m0(:, 1:nx) = matmul(map%op_singlet(unop)%m3, map%fc_singlet_shell(unsh)%coeff)
        do k = 1, nx
        do ii = 1, 3
            uind = map%fc_singlet_shell(unsh)%ind_global(k)
            find = (a1 - 1)*3 + ii
            CM(find, uind) = CM(find, uind) + m0(ii, k)
        end do
        end do
    end do

    deallocate (m0)
end subroutine

!> get the pair coefficient matrix
module subroutine lo_coeffmatrix_pair(UM, CM, map)
    !> displacements
    real(r8), dimension(:, :), intent(in) :: UM
    !> coefficient matrix
    real(r8), dimension(:, :), intent(out) :: CM
    !> symmetry stuff rearranged into coordination shells
    type(lo_forcemap), intent(in) :: map

    integer :: a1, a2, k, ii
    integer :: find, uind, unop, unsh, nx, ipair
    integer, dimension(9) :: xind
    real(r8), dimension(3) :: u, v
    real(r8) :: wm9x1(9, 1), wm9x2(9, 2), wm9x3(9, 3), wm9x4(9, 4), wm9x5(9, 5), wm9x6(9, 6), wm9x7(9, 7), wm9x8(9, 8), wm9x9(9, 9)
    real(r8) :: cm1x3(1, 3), cm2x3(2, 3), cm3x3(3, 3), cm4x3(4, 3), cm5x3(5, 3), cm6x3(6, 3), cm7x3(7, 3), cm8x3(8, 3), cm9x3(9, 3)

    CM = 0.0_r8
    do ipair = 1, map%xss%n_fc_pair
        unsh = map%xss%ind_fc_pair(1, ipair)
        unop = map%xss%ind_fc_pair(2, ipair)
        nx = map%fc_pair_shell(unsh)%nx
        if (nx .eq. 0) cycle
        a1 = map%xss%ind_fc_pair(3, ipair)
        a2 = map%xss%ind_fc_pair(4, ipair)
        v = UM(:, a1)
        u = UM(:, a2)
        xind(1:nx) = map%fc_pair_shell(unsh)%ind_global
        select case (nx)
        case (1)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x1)
            cm1x3 = coeff9x1(wm9x1, u) - coeff9x1(wm9x1, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm1x3(k, ii)
            end do
            end do
        case (2)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x2)
            cm2x3 = coeff9x2(wm9x2, u) - coeff9x2(wm9x2, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm2x3(k, ii)
            end do
            end do
        case (3)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x3)
            cm3x3 = coeff9x3(wm9x3, u) - coeff9x3(wm9x3, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm3x3(k, ii)
            end do
            end do
        case (4)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x4)
            cm4x3 = coeff9x4(wm9x4, u) - coeff9x4(wm9x4, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm4x3(k, ii)
            end do
            end do
        case (5)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x5)
            cm5x3 = coeff9x5(wm9x5, u) - coeff9x5(wm9x5, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm5x3(k, ii)
            end do
            end do
        case (6)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x6)
            cm6x3 = coeff9x6(wm9x6, u) - coeff9x6(wm9x6, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm6x3(k, ii)
            end do
            end do
        case (7)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x7)
            cm7x3 = coeff9x7(wm9x7, u) - coeff9x7(wm9x7, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm7x3(k, ii)
            end do
            end do
        case (8)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x8)
            cm8x3 = coeff9x8(wm9x8, u) - coeff9x8(wm9x8, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm8x3(k, ii)
            end do
            end do
        case (9)
            call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x9)
            cm9x3 = coeff9x9(wm9x9, u) - coeff9x9(wm9x9, v)
            do k = 1, nx
            do ii = 1, 3
                uind = xind(k)
                find = (a1 - 1)*3 + ii
                CM(find, uind) = CM(find, uind) + cm9x3(k, ii)
            end do
            end do
        end select
    end do
end subroutine

! !> does the kronecker product I3 otimes u
! pure function fc_pair_kronexpand(u,v) result(A)
!     !> displacements
!     real(r8), dimension(3), intent(in) :: u,v
!     !> kronecker product
!     real(r8), dimension(3,9) :: A
!
!     integer :: i
!     A=0.0_r8
!     do i=1,3
!         A(i,(i-1)*3+1:i*3)=u-v
!     enddo
! end function

! !> coefficient matrix for magnetic longitudinal fluctuations
! module subroutine lo_coeffmatrix_magnetic_longitudinal(UM,CM,map)
!     !> displacements
!     real(r8), dimension(:,:), intent(in) :: UM
!     !> coefficient matrix
!     real(r8), dimension(:,:), intent(out) :: CM
!     !> symmetry stuff rearranged into coordination shells
!     type(lo_forcemap), intent(in) :: map
!
!     integer :: a1,i
!     integer :: unsh,unop,ntheta
!     integer, dimension(3) :: thetaind
!     real(r8) :: wm3x1(3,1),wm3x2(3,2),wm3x3(3,3)
!     real(r8), dimension(3) :: v
!
!     write(*,*) 'FIXME MAGNETIC LONG'
!     stop
!     ! CM=0.0_r8
!     ! do a1=1,map%nss
!     ! do i=1,map%ss(a1)%nmagpair
!     !     ! second atom, displacement
!     !     v=UM(:,map%ss(a1)%magpairind(3,i))
!     !     unsh=map%ss(a1)%magpairind(1,i)
!     !     unop=map%ss(a1)%magpairind(2,i)
!     !     ntheta=map%magpairshell(unsh)%ntheta_qij
!     !     thetaind(1:ntheta)=map%magpairshell(unsh)%ind_global_qij
!     !     select case(ntheta)
!     !     case(1)
!     !         wm3x1=matmul(map%op_pair(unop)%m3,map%magpairshell(unsh)%coeff_qij)
!     !         CM(a1,thetaind(1:ntheta))=CM(a1,thetaind(1:ntheta))+lcoeff3x1(wm3x1,v)
!     !     case(2)
!     !         wm3x2=matmul(map%op_pair(unop)%m3,map%magpairshell(unsh)%coeff_qij)
!     !         CM(a1,thetaind(1:ntheta))=CM(a1,thetaind(1:ntheta))+lcoeff3x2(wm3x2,v)
!     !     case(3)
!     !         wm3x3=matmul(map%op_pair(unop)%m3,map%magpairshell(unsh)%coeff_qij)
!     !         CM(a1,thetaind(1:ntheta))=CM(a1,thetaind(1:ntheta))+lcoeff3x3(wm3x3,v)
!     !     end select
!     ! enddo
!     ! enddo
! end subroutine

! !> coefficient matrix for magnetic crossterm interactions
! module subroutine lo_coeffmatrix_magnetic_crossterm_energy(UM,MM,CM,map)
!     !> displacements
!     real(r8), dimension(:,:), intent(in) :: UM
!     !> magnetic moments
!     real(r8), dimension(:,:), intent(in) :: MM
!     !> coefficient matrix
!     real(r8), dimension(:), intent(out) :: CM
!     !> symmetry stuff rearranged into coordination shells
!     type(lo_forcemap), intent(in) :: map
!
!     integer :: a1,i
!     integer :: unsh,unop,ntheta
!     integer, dimension(9) :: thetaind
!     real(r8), dimension(3) :: u,v
!     real(r8) :: wm9x1(9,1),wm9x2(9,2),wm9x3(9,3),wm9x4(9,4),wm9x5(9,5),wm9x6(9,6),wm9x7(9,7),wm9x8(9,8),wm9x9(9,9)
!
!     write(*,*) 'FIXME MAGNETIC CROSSTERM ENERGY'
!     stop
!
!     ! CM=0.0_r8
!     ! do a1=1,map%nss
!     !     ! first atom, displacement
!     !     u=UM(:,a1)
!     !     do i=1,map%ss(a1)%nmagpair
!     !         ! second atom, magnetic moment
!     !         v=MM(:,map%ss(a1)%magpairind(3,i))
!     !         unsh=map%ss(a1)%magpairind(1,i)
!     !         unop=map%ss(a1)%magpairind(2,i)
!     !         ntheta=map%magpairshell(unsh)%ntheta_tij
!     !         thetaind(1:ntheta)=map%magpairshell(unsh)%ind_global_tij
!     !         select case(ntheta)
!     !         case(1)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x1 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x1(wm9x1,u,v)
!     !         case(2)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x2 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x2(wm9x2,u,v)
!     !         case(3)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x3 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x3(wm9x3,u,v)
!     !         case(4)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x4 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x4(wm9x4,u,v)
!     !         case(5)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x5 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x5(wm9x5,u,v)
!     !         case(6)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x6 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x6(wm9x6,u,v)
!     !         case(7)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x7 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x7(wm9x7,u,v)
!     !         case(8)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x8 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x8(wm9x8,u,v)
!     !         case(9)
!     !             call lo_gemm(map%op_pair(unop)%msotr,map%magpairshell(unsh)%coeff_tij, wm9x9 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x9(wm9x9,u,v)
!     !         end select
!     !     enddo
!     ! enddo
!     ! ! And factor 2 to be consistent with the phonon things
!     ! CM=CM*0.5_r8
! end subroutine

! !> get the pair coefficient matrix
! module subroutine lo_coeffmatrix_magnetic_crossterm_forces(UM,CM,map)
!     !> magnetic moments
!     real(r8), dimension(:,:), intent(in) :: UM
!     !> coefficient matrix
!     real(r8), dimension(:,:), intent(out) :: CM
!     !> symmetry stuff rearranged into coordination shells
!     type(lo_forcemap), intent(in) :: map
!
!     integer :: a1,i,k,ii
!     integer :: find,uind,unop,unsh,ntheta
!     integer, dimension(9) :: thetaind
!     real(r8), dimension(9,9) :: opm
!     real(r8), dimension(3) :: u,v
!     real(r8) :: wm9x1(9,1),wm9x2(9,2),wm9x3(9,3),wm9x4(9,4),wm9x5(9,5),wm9x6(9,6),wm9x7(9,7),wm9x8(9,8),wm9x9(9,9)
!     real(r8) :: cm1x3(1,3),cm2x3(2,3),cm3x3(3,3),cm4x3(4,3),cm5x3(5,3),cm6x3(6,3),cm7x3(7,3),cm8x3(8,3),cm9x3(9,3)
!
!     write(*,*) 'FIXME MAGNETIC CROSSTERM FORCES'
!     stop
!
!     ! CM=0.0_r8
!     ! do a1=1,map%nss
!     !     ! displacement of this atom, for the self-term
!     !     v=0.0_r8 !UM(:,a1)
!     !     do i=1,map%ss(a1)%nmagpair
!     !         ! the displacement
!     !         u=UM(:,map%ss(a1)%magpairind(3,i))
!     !         ! operation and shell
!     !         unsh=map%ss(a1)%magpairind(1,i)
!     !         unop=map%ss(a1)%magpairind(2,i)
!     !         ! get the operation with the correct sign
!     !         opm=map%op_pair(unop)%msotr
!     !         ntheta=map%magpairshell(unsh)%ntheta_tij
!     !         ! Skip if no thetas here
!     !         if ( ntheta .eq. 0 ) cycle
!     !         thetaind(1:ntheta)=map%magpairshell(unsh)%ind_global_tij
!     !         select case(ntheta)
!     !         case(1)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x1 )
!     !             cm1x3=coeff9x1(wm9x1,u)-coeff9x1(wm9x1,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm1x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(2)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x2 )
!     !             cm2x3=coeff9x2(wm9x2,u)-coeff9x2(wm9x2,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm2x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(3)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x3 )
!     !             cm3x3=coeff9x3(wm9x3,u)-coeff9x3(wm9x3,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm3x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(4)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x4 )
!     !             cm4x3=coeff9x4(wm9x4,u)-coeff9x4(wm9x4,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm4x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(5)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x5 )
!     !             cm5x3=coeff9x5(wm9x5,u)-coeff9x5(wm9x5,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm5x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(6)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x6 )
!     !             cm6x3=coeff9x6(wm9x6,u)-coeff9x6(wm9x6,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm6x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(7)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x7 )
!     !             cm7x3=coeff9x7(wm9x7,u)-coeff9x7(wm9x7,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm7x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(8)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x8 )
!     !             cm8x3=coeff9x8(wm9x8,u)-coeff9x8(wm9x8,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm8x3(k,ii)
!     !             enddo
!     !             enddo
!     !         case(9)
!     !             call lo_gemm(opm,map%magpairshell(unsh)%coeff_tij, wm9x9 )
!     !             cm9x3=coeff9x9(wm9x9,u)-coeff9x9(wm9x9,v)
!     !             do k=1,ntheta
!     !             do ii=1,3
!     !                 uind=thetaind(k)
!     !                 find=(a1-1)*3+ii
!     !                 CM(find,uind)=CM(find,uind)+cm9x3(k,ii)
!     !             enddo
!     !             enddo
!     !         end select
!     !     enddo
!     ! enddo
!     !
!     ! contains
!     ! ! Procedureally generated stuff
!     ! pure function coeff9x1(cl,u) result(m)
!     !     real(r8), dimension(9,1), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(1,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     ! end function
!     ! pure function coeff9x2(cl,u) result(m)
!     !     real(r8), dimension(9,2), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(2,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     ! end function
!     ! pure function coeff9x3(cl,u) result(m)
!     !     real(r8), dimension(9,3), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(3,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(3,1)=cl(1,3)*u2x + cl(2,3)*u2y + cl(3,3)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(3,2)=cl(4,3)*u2x + cl(5,3)*u2y + cl(6,3)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     !     m(3,3)=cl(7,3)*u2x + cl(8,3)*u2y + cl(9,3)*u2z
!     ! end function
!     ! pure function coeff9x4(cl,u) result(m)
!     !     real(r8), dimension(9,4), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(4,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(3,1)=cl(1,3)*u2x + cl(2,3)*u2y + cl(3,3)*u2z
!     !     m(4,1)=cl(1,4)*u2x + cl(2,4)*u2y + cl(3,4)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(3,2)=cl(4,3)*u2x + cl(5,3)*u2y + cl(6,3)*u2z
!     !     m(4,2)=cl(4,4)*u2x + cl(5,4)*u2y + cl(6,4)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     !     m(3,3)=cl(7,3)*u2x + cl(8,3)*u2y + cl(9,3)*u2z
!     !     m(4,3)=cl(7,4)*u2x + cl(8,4)*u2y + cl(9,4)*u2z
!     ! end function
!     ! pure function coeff9x5(cl,u) result(m)
!     !     real(r8), dimension(9,5), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(5,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(3,1)=cl(1,3)*u2x + cl(2,3)*u2y + cl(3,3)*u2z
!     !     m(4,1)=cl(1,4)*u2x + cl(2,4)*u2y + cl(3,4)*u2z
!     !     m(5,1)=cl(1,5)*u2x + cl(2,5)*u2y + cl(3,5)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(3,2)=cl(4,3)*u2x + cl(5,3)*u2y + cl(6,3)*u2z
!     !     m(4,2)=cl(4,4)*u2x + cl(5,4)*u2y + cl(6,4)*u2z
!     !     m(5,2)=cl(4,5)*u2x + cl(5,5)*u2y + cl(6,5)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     !     m(3,3)=cl(7,3)*u2x + cl(8,3)*u2y + cl(9,3)*u2z
!     !     m(4,3)=cl(7,4)*u2x + cl(8,4)*u2y + cl(9,4)*u2z
!     !     m(5,3)=cl(7,5)*u2x + cl(8,5)*u2y + cl(9,5)*u2z
!     ! end function
!     ! pure function coeff9x6(cl,u) result(m)
!     !     real(r8), dimension(9,6), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(6,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(3,1)=cl(1,3)*u2x + cl(2,3)*u2y + cl(3,3)*u2z
!     !     m(4,1)=cl(1,4)*u2x + cl(2,4)*u2y + cl(3,4)*u2z
!     !     m(5,1)=cl(1,5)*u2x + cl(2,5)*u2y + cl(3,5)*u2z
!     !     m(6,1)=cl(1,6)*u2x + cl(2,6)*u2y + cl(3,6)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(3,2)=cl(4,3)*u2x + cl(5,3)*u2y + cl(6,3)*u2z
!     !     m(4,2)=cl(4,4)*u2x + cl(5,4)*u2y + cl(6,4)*u2z
!     !     m(5,2)=cl(4,5)*u2x + cl(5,5)*u2y + cl(6,5)*u2z
!     !     m(6,2)=cl(4,6)*u2x + cl(5,6)*u2y + cl(6,6)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     !     m(3,3)=cl(7,3)*u2x + cl(8,3)*u2y + cl(9,3)*u2z
!     !     m(4,3)=cl(7,4)*u2x + cl(8,4)*u2y + cl(9,4)*u2z
!     !     m(5,3)=cl(7,5)*u2x + cl(8,5)*u2y + cl(9,5)*u2z
!     !     m(6,3)=cl(7,6)*u2x + cl(8,6)*u2y + cl(9,6)*u2z
!     ! end function
!     ! pure function coeff9x7(cl,u) result(m)
!     !     real(r8), dimension(9,7), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(7,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(3,1)=cl(1,3)*u2x + cl(2,3)*u2y + cl(3,3)*u2z
!     !     m(4,1)=cl(1,4)*u2x + cl(2,4)*u2y + cl(3,4)*u2z
!     !     m(5,1)=cl(1,5)*u2x + cl(2,5)*u2y + cl(3,5)*u2z
!     !     m(6,1)=cl(1,6)*u2x + cl(2,6)*u2y + cl(3,6)*u2z
!     !     m(7,1)=cl(1,7)*u2x + cl(2,7)*u2y + cl(3,7)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(3,2)=cl(4,3)*u2x + cl(5,3)*u2y + cl(6,3)*u2z
!     !     m(4,2)=cl(4,4)*u2x + cl(5,4)*u2y + cl(6,4)*u2z
!     !     m(5,2)=cl(4,5)*u2x + cl(5,5)*u2y + cl(6,5)*u2z
!     !     m(6,2)=cl(4,6)*u2x + cl(5,6)*u2y + cl(6,6)*u2z
!     !     m(7,2)=cl(4,7)*u2x + cl(5,7)*u2y + cl(6,7)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     !     m(3,3)=cl(7,3)*u2x + cl(8,3)*u2y + cl(9,3)*u2z
!     !     m(4,3)=cl(7,4)*u2x + cl(8,4)*u2y + cl(9,4)*u2z
!     !     m(5,3)=cl(7,5)*u2x + cl(8,5)*u2y + cl(9,5)*u2z
!     !     m(6,3)=cl(7,6)*u2x + cl(8,6)*u2y + cl(9,6)*u2z
!     !     m(7,3)=cl(7,7)*u2x + cl(8,7)*u2y + cl(9,7)*u2z
!     ! end function
!     ! pure function coeff9x8(cl,u) result(m)
!     !     real(r8), dimension(9,8), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(8,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(3,1)=cl(1,3)*u2x + cl(2,3)*u2y + cl(3,3)*u2z
!     !     m(4,1)=cl(1,4)*u2x + cl(2,4)*u2y + cl(3,4)*u2z
!     !     m(5,1)=cl(1,5)*u2x + cl(2,5)*u2y + cl(3,5)*u2z
!     !     m(6,1)=cl(1,6)*u2x + cl(2,6)*u2y + cl(3,6)*u2z
!     !     m(7,1)=cl(1,7)*u2x + cl(2,7)*u2y + cl(3,7)*u2z
!     !     m(8,1)=cl(1,8)*u2x + cl(2,8)*u2y + cl(3,8)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(3,2)=cl(4,3)*u2x + cl(5,3)*u2y + cl(6,3)*u2z
!     !     m(4,2)=cl(4,4)*u2x + cl(5,4)*u2y + cl(6,4)*u2z
!     !     m(5,2)=cl(4,5)*u2x + cl(5,5)*u2y + cl(6,5)*u2z
!     !     m(6,2)=cl(4,6)*u2x + cl(5,6)*u2y + cl(6,6)*u2z
!     !     m(7,2)=cl(4,7)*u2x + cl(5,7)*u2y + cl(6,7)*u2z
!     !     m(8,2)=cl(4,8)*u2x + cl(5,8)*u2y + cl(6,8)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     !     m(3,3)=cl(7,3)*u2x + cl(8,3)*u2y + cl(9,3)*u2z
!     !     m(4,3)=cl(7,4)*u2x + cl(8,4)*u2y + cl(9,4)*u2z
!     !     m(5,3)=cl(7,5)*u2x + cl(8,5)*u2y + cl(9,5)*u2z
!     !     m(6,3)=cl(7,6)*u2x + cl(8,6)*u2y + cl(9,6)*u2z
!     !     m(7,3)=cl(7,7)*u2x + cl(8,7)*u2y + cl(9,7)*u2z
!     !     m(8,3)=cl(7,8)*u2x + cl(8,8)*u2y + cl(9,8)*u2z
!     ! end function
!     ! pure function coeff9x9(cl,u) result(m)
!     !     real(r8), dimension(9,9), intent(in) :: cl
!     !     real(r8), dimension(3), intent(in) :: u
!     !     real(r8), dimension(9,3) :: m
!     !     real(r8) :: u2x,u2y,u2z
!     !     u2x=u(1); u2y=u(2); u2z=u(3)
!     !     m(1,1)=cl(1,1)*u2x + cl(2,1)*u2y + cl(3,1)*u2z
!     !     m(2,1)=cl(1,2)*u2x + cl(2,2)*u2y + cl(3,2)*u2z
!     !     m(3,1)=cl(1,3)*u2x + cl(2,3)*u2y + cl(3,3)*u2z
!     !     m(4,1)=cl(1,4)*u2x + cl(2,4)*u2y + cl(3,4)*u2z
!     !     m(5,1)=cl(1,5)*u2x + cl(2,5)*u2y + cl(3,5)*u2z
!     !     m(6,1)=cl(1,6)*u2x + cl(2,6)*u2y + cl(3,6)*u2z
!     !     m(7,1)=cl(1,7)*u2x + cl(2,7)*u2y + cl(3,7)*u2z
!     !     m(8,1)=cl(1,8)*u2x + cl(2,8)*u2y + cl(3,8)*u2z
!     !     m(9,1)=cl(1,9)*u2x + cl(2,9)*u2y + cl(3,9)*u2z
!     !     m(1,2)=cl(4,1)*u2x + cl(5,1)*u2y + cl(6,1)*u2z
!     !     m(2,2)=cl(4,2)*u2x + cl(5,2)*u2y + cl(6,2)*u2z
!     !     m(3,2)=cl(4,3)*u2x + cl(5,3)*u2y + cl(6,3)*u2z
!     !     m(4,2)=cl(4,4)*u2x + cl(5,4)*u2y + cl(6,4)*u2z
!     !     m(5,2)=cl(4,5)*u2x + cl(5,5)*u2y + cl(6,5)*u2z
!     !     m(6,2)=cl(4,6)*u2x + cl(5,6)*u2y + cl(6,6)*u2z
!     !     m(7,2)=cl(4,7)*u2x + cl(5,7)*u2y + cl(6,7)*u2z
!     !     m(8,2)=cl(4,8)*u2x + cl(5,8)*u2y + cl(6,8)*u2z
!     !     m(9,2)=cl(4,9)*u2x + cl(5,9)*u2y + cl(6,9)*u2z
!     !     m(1,3)=cl(7,1)*u2x + cl(8,1)*u2y + cl(9,1)*u2z
!     !     m(2,3)=cl(7,2)*u2x + cl(8,2)*u2y + cl(9,2)*u2z
!     !     m(3,3)=cl(7,3)*u2x + cl(8,3)*u2y + cl(9,3)*u2z
!     !     m(4,3)=cl(7,4)*u2x + cl(8,4)*u2y + cl(9,4)*u2z
!     !     m(5,3)=cl(7,5)*u2x + cl(8,5)*u2y + cl(9,5)*u2z
!     !     m(6,3)=cl(7,6)*u2x + cl(8,6)*u2y + cl(9,6)*u2z
!     !     m(7,3)=cl(7,7)*u2x + cl(8,7)*u2y + cl(9,7)*u2z
!     !     m(8,3)=cl(7,8)*u2x + cl(8,8)*u2y + cl(9,8)*u2z
!     !     m(9,3)=cl(7,9)*u2x + cl(8,9)*u2y + cl(9,9)*u2z
!     ! end function
! end subroutine

! !> coefficient matrix for magnetic pair interactions
! module subroutine lo_coeffmatrix_magnetic_pair(UM,CM,map)
!     !> magnetic moments
!     real(r8), dimension(:,:), intent(in) :: UM
!     !> coefficient matrix
!     real(r8), dimension(:), intent(out) :: CM
!     !> symmetry stuff rearranged into coordination shells
!     type(lo_forcemap), intent(in) :: map
!     !
!     integer :: a1,i
!     integer :: unsh,unop,ntheta
!     integer, dimension(9) :: thetaind
!     real(r8), dimension(3) :: u,v
!     real(r8) :: wm9x1(9,1),wm9x2(9,2),wm9x3(9,3),wm9x4(9,4),wm9x5(9,5),wm9x6(9,6),wm9x7(9,7),wm9x8(9,8),wm9x9(9,9)
!
!     write(*,*) 'FIXME MAGNETIC PAIR'
!     stop
!
!     ! CM=0.0_r8
!     ! do a1=1,map%nss
!     !     ! first magnetic moment
!     !     u=UM(:,a1)
!     !     do i=1,map%ss(a1)%nmagpair
!     !         ! second magnetic moment
!     !         v=UM(:,map%ss(a1)%magpairind(3,i))
!     !         unsh=map%ss(a1)%magpairind(1,i)
!     !         unop=map%ss(a1)%magpairind(2,i)
!     !         ntheta=map%magpairshell(unsh)%ntheta_jij
!     !         thetaind(1:ntheta)=map%magpairshell(unsh)%ind_global_jij
!     !         select case(ntheta)
!     !         case(1)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x1 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x1(wm9x1,u,v)
!     !         case(2)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x2 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x2(wm9x2,u,v)
!     !         case(3)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x3 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x3(wm9x3,u,v)
!     !         case(4)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x4 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x4(wm9x4,u,v)
!     !         case(5)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x5 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x5(wm9x5,u,v)
!     !         case(6)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x6 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x6(wm9x6,u,v)
!     !         case(7)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x7 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x7(wm9x7,u,v)
!     !         case(8)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x8 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x8(wm9x8,u,v)
!     !         case(9)
!     !             call lo_gemm(map%op_pair(unop)%sotr,map%magpairshell(unsh)%coeff_jij, wm9x9 )
!     !             CM(thetaind(1:ntheta))=CM(thetaind(1:ntheta))+magcoeff9x9(wm9x9,u,v)
!     !         end select
!     !     enddo
!     ! enddo
!     ! ! And factor 2 to be consistent with the phonon things
!     ! CM=CM*0.5_r8
! end subroutine

! Below is not exposed, procedurally generated things

! ! Procedurally generated coefficients below
! pure function magcoeff9x1(cm,u,v) result(m)
!     real(r8), dimension(9,1), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(1) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
! end function
!
! pure function magcoeff9x2(cm,u,v) result(m)
!     real(r8), dimension(9,2), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(2) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
! end function
!
! pure function magcoeff9x3(cm,u,v) result(m)
!     real(r8), dimension(9,3), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(3) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
!     m(3)=cm(1,3)*ux*vx + cm(4,3)*uy*vx + cm(7,3)*uz*vx + cm(2,3)*ux*vy + cm(5,3)*uy*vy + cm(8,3)*uz*vy + cm(3,3)*ux*vz + cm(6,3)*uy*vz + cm(9,3)*uz*vz
! end function
!
! pure function magcoeff9x4(cm,u,v) result(m)
!     real(r8), dimension(9,4), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(4) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
!     m(3)=cm(1,3)*ux*vx + cm(4,3)*uy*vx + cm(7,3)*uz*vx + cm(2,3)*ux*vy + cm(5,3)*uy*vy + cm(8,3)*uz*vy + cm(3,3)*ux*vz + cm(6,3)*uy*vz + cm(9,3)*uz*vz
!     m(4)=cm(1,4)*ux*vx + cm(4,4)*uy*vx + cm(7,4)*uz*vx + cm(2,4)*ux*vy + cm(5,4)*uy*vy + cm(8,4)*uz*vy + cm(3,4)*ux*vz + cm(6,4)*uy*vz + cm(9,4)*uz*vz
! end function
!
! pure function magcoeff9x5(cm,u,v) result(m)
!     real(r8), dimension(9,5), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(5) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
!     m(3)=cm(1,3)*ux*vx + cm(4,3)*uy*vx + cm(7,3)*uz*vx + cm(2,3)*ux*vy + cm(5,3)*uy*vy + cm(8,3)*uz*vy + cm(3,3)*ux*vz + cm(6,3)*uy*vz + cm(9,3)*uz*vz
!     m(4)=cm(1,4)*ux*vx + cm(4,4)*uy*vx + cm(7,4)*uz*vx + cm(2,4)*ux*vy + cm(5,4)*uy*vy + cm(8,4)*uz*vy + cm(3,4)*ux*vz + cm(6,4)*uy*vz + cm(9,4)*uz*vz
!     m(5)=cm(1,5)*ux*vx + cm(4,5)*uy*vx + cm(7,5)*uz*vx + cm(2,5)*ux*vy + cm(5,5)*uy*vy + cm(8,5)*uz*vy + cm(3,5)*ux*vz + cm(6,5)*uy*vz + cm(9,5)*uz*vz
! end function
!
! pure function magcoeff9x6(cm,u,v) result(m)
!     real(r8), dimension(9,6), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(6) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
!     m(3)=cm(1,3)*ux*vx + cm(4,3)*uy*vx + cm(7,3)*uz*vx + cm(2,3)*ux*vy + cm(5,3)*uy*vy + cm(8,3)*uz*vy + cm(3,3)*ux*vz + cm(6,3)*uy*vz + cm(9,3)*uz*vz
!     m(4)=cm(1,4)*ux*vx + cm(4,4)*uy*vx + cm(7,4)*uz*vx + cm(2,4)*ux*vy + cm(5,4)*uy*vy + cm(8,4)*uz*vy + cm(3,4)*ux*vz + cm(6,4)*uy*vz + cm(9,4)*uz*vz
!     m(5)=cm(1,5)*ux*vx + cm(4,5)*uy*vx + cm(7,5)*uz*vx + cm(2,5)*ux*vy + cm(5,5)*uy*vy + cm(8,5)*uz*vy + cm(3,5)*ux*vz + cm(6,5)*uy*vz + cm(9,5)*uz*vz
!     m(6)=cm(1,6)*ux*vx + cm(4,6)*uy*vx + cm(7,6)*uz*vx + cm(2,6)*ux*vy + cm(5,6)*uy*vy + cm(8,6)*uz*vy + cm(3,6)*ux*vz + cm(6,6)*uy*vz + cm(9,6)*uz*vz
! end function
!
! pure function magcoeff9x7(cm,u,v) result(m)
!     real(r8), dimension(9,7), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(7) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
!     m(3)=cm(1,3)*ux*vx + cm(4,3)*uy*vx + cm(7,3)*uz*vx + cm(2,3)*ux*vy + cm(5,3)*uy*vy + cm(8,3)*uz*vy + cm(3,3)*ux*vz + cm(6,3)*uy*vz + cm(9,3)*uz*vz
!     m(4)=cm(1,4)*ux*vx + cm(4,4)*uy*vx + cm(7,4)*uz*vx + cm(2,4)*ux*vy + cm(5,4)*uy*vy + cm(8,4)*uz*vy + cm(3,4)*ux*vz + cm(6,4)*uy*vz + cm(9,4)*uz*vz
!     m(5)=cm(1,5)*ux*vx + cm(4,5)*uy*vx + cm(7,5)*uz*vx + cm(2,5)*ux*vy + cm(5,5)*uy*vy + cm(8,5)*uz*vy + cm(3,5)*ux*vz + cm(6,5)*uy*vz + cm(9,5)*uz*vz
!     m(6)=cm(1,6)*ux*vx + cm(4,6)*uy*vx + cm(7,6)*uz*vx + cm(2,6)*ux*vy + cm(5,6)*uy*vy + cm(8,6)*uz*vy + cm(3,6)*ux*vz + cm(6,6)*uy*vz + cm(9,6)*uz*vz
!     m(7)=cm(1,7)*ux*vx + cm(4,7)*uy*vx + cm(7,7)*uz*vx + cm(2,7)*ux*vy + cm(5,7)*uy*vy + cm(8,7)*uz*vy + cm(3,7)*ux*vz + cm(6,7)*uy*vz + cm(9,7)*uz*vz
! end function
!
! pure function magcoeff9x8(cm,u,v) result(m)
!     real(r8), dimension(9,8), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(8) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
!     m(3)=cm(1,3)*ux*vx + cm(4,3)*uy*vx + cm(7,3)*uz*vx + cm(2,3)*ux*vy + cm(5,3)*uy*vy + cm(8,3)*uz*vy + cm(3,3)*ux*vz + cm(6,3)*uy*vz + cm(9,3)*uz*vz
!     m(4)=cm(1,4)*ux*vx + cm(4,4)*uy*vx + cm(7,4)*uz*vx + cm(2,4)*ux*vy + cm(5,4)*uy*vy + cm(8,4)*uz*vy + cm(3,4)*ux*vz + cm(6,4)*uy*vz + cm(9,4)*uz*vz
!     m(5)=cm(1,5)*ux*vx + cm(4,5)*uy*vx + cm(7,5)*uz*vx + cm(2,5)*ux*vy + cm(5,5)*uy*vy + cm(8,5)*uz*vy + cm(3,5)*ux*vz + cm(6,5)*uy*vz + cm(9,5)*uz*vz
!     m(6)=cm(1,6)*ux*vx + cm(4,6)*uy*vx + cm(7,6)*uz*vx + cm(2,6)*ux*vy + cm(5,6)*uy*vy + cm(8,6)*uz*vy + cm(3,6)*ux*vz + cm(6,6)*uy*vz + cm(9,6)*uz*vz
!     m(7)=cm(1,7)*ux*vx + cm(4,7)*uy*vx + cm(7,7)*uz*vx + cm(2,7)*ux*vy + cm(5,7)*uy*vy + cm(8,7)*uz*vy + cm(3,7)*ux*vz + cm(6,7)*uy*vz + cm(9,7)*uz*vz
!     m(8)=cm(1,8)*ux*vx + cm(4,8)*uy*vx + cm(7,8)*uz*vx + cm(2,8)*ux*vy + cm(5,8)*uy*vy + cm(8,8)*uz*vy + cm(3,8)*ux*vz + cm(6,8)*uy*vz + cm(9,8)*uz*vz
! end function
!
! pure function magcoeff9x9(cm,u,v) result(m)
!     real(r8), dimension(9,9), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u,v
!     real(r8), dimension(9) :: m
!     real(r8) :: ux,uy,uz,vx,vy,vz
!     ux=u(1); uy=u(2); uz=u(3); vx=v(1); vy=v(2); vz=v(3)
!     m(1)=cm(1,1)*ux*vx + cm(4,1)*uy*vx + cm(7,1)*uz*vx + cm(2,1)*ux*vy + cm(5,1)*uy*vy + cm(8,1)*uz*vy + cm(3,1)*ux*vz + cm(6,1)*uy*vz + cm(9,1)*uz*vz
!     m(2)=cm(1,2)*ux*vx + cm(4,2)*uy*vx + cm(7,2)*uz*vx + cm(2,2)*ux*vy + cm(5,2)*uy*vy + cm(8,2)*uz*vy + cm(3,2)*ux*vz + cm(6,2)*uy*vz + cm(9,2)*uz*vz
!     m(3)=cm(1,3)*ux*vx + cm(4,3)*uy*vx + cm(7,3)*uz*vx + cm(2,3)*ux*vy + cm(5,3)*uy*vy + cm(8,3)*uz*vy + cm(3,3)*ux*vz + cm(6,3)*uy*vz + cm(9,3)*uz*vz
!     m(4)=cm(1,4)*ux*vx + cm(4,4)*uy*vx + cm(7,4)*uz*vx + cm(2,4)*ux*vy + cm(5,4)*uy*vy + cm(8,4)*uz*vy + cm(3,4)*ux*vz + cm(6,4)*uy*vz + cm(9,4)*uz*vz
!     m(5)=cm(1,5)*ux*vx + cm(4,5)*uy*vx + cm(7,5)*uz*vx + cm(2,5)*ux*vy + cm(5,5)*uy*vy + cm(8,5)*uz*vy + cm(3,5)*ux*vz + cm(6,5)*uy*vz + cm(9,5)*uz*vz
!     m(6)=cm(1,6)*ux*vx + cm(4,6)*uy*vx + cm(7,6)*uz*vx + cm(2,6)*ux*vy + cm(5,6)*uy*vy + cm(8,6)*uz*vy + cm(3,6)*ux*vz + cm(6,6)*uy*vz + cm(9,6)*uz*vz
!     m(7)=cm(1,7)*ux*vx + cm(4,7)*uy*vx + cm(7,7)*uz*vx + cm(2,7)*ux*vy + cm(5,7)*uy*vy + cm(8,7)*uz*vy + cm(3,7)*ux*vz + cm(6,7)*uy*vz + cm(9,7)*uz*vz
!     m(8)=cm(1,8)*ux*vx + cm(4,8)*uy*vx + cm(7,8)*uz*vx + cm(2,8)*ux*vy + cm(5,8)*uy*vy + cm(8,8)*uz*vy + cm(3,8)*ux*vz + cm(6,8)*uy*vz + cm(9,8)*uz*vz
!     m(9)=cm(1,9)*ux*vx + cm(4,9)*uy*vx + cm(7,9)*uz*vx + cm(2,9)*ux*vy + cm(5,9)*uy*vy + cm(8,9)*uz*vy + cm(3,9)*ux*vz + cm(6,9)*uy*vz + cm(9,9)*uz*vz
! end function

! Procedureally generated stuff
pure function coeff9x1(cl, u) result(m)
    real(r8), dimension(9, 1), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(1, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
end function
pure function coeff9x2(cl, u) result(m)
    real(r8), dimension(9, 2), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(2, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
end function
pure function coeff9x3(cl, u) result(m)
    real(r8), dimension(9, 3), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(3, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(3, 1) = cl(1, 3)*u2x + cl(2, 3)*u2y + cl(3, 3)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(3, 2) = cl(4, 3)*u2x + cl(5, 3)*u2y + cl(6, 3)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
    m(3, 3) = cl(7, 3)*u2x + cl(8, 3)*u2y + cl(9, 3)*u2z
end function
pure function coeff9x4(cl, u) result(m)
    real(r8), dimension(9, 4), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(4, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(3, 1) = cl(1, 3)*u2x + cl(2, 3)*u2y + cl(3, 3)*u2z
    m(4, 1) = cl(1, 4)*u2x + cl(2, 4)*u2y + cl(3, 4)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(3, 2) = cl(4, 3)*u2x + cl(5, 3)*u2y + cl(6, 3)*u2z
    m(4, 2) = cl(4, 4)*u2x + cl(5, 4)*u2y + cl(6, 4)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
    m(3, 3) = cl(7, 3)*u2x + cl(8, 3)*u2y + cl(9, 3)*u2z
    m(4, 3) = cl(7, 4)*u2x + cl(8, 4)*u2y + cl(9, 4)*u2z
end function
pure function coeff9x5(cl, u) result(m)
    real(r8), dimension(9, 5), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(5, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(3, 1) = cl(1, 3)*u2x + cl(2, 3)*u2y + cl(3, 3)*u2z
    m(4, 1) = cl(1, 4)*u2x + cl(2, 4)*u2y + cl(3, 4)*u2z
    m(5, 1) = cl(1, 5)*u2x + cl(2, 5)*u2y + cl(3, 5)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(3, 2) = cl(4, 3)*u2x + cl(5, 3)*u2y + cl(6, 3)*u2z
    m(4, 2) = cl(4, 4)*u2x + cl(5, 4)*u2y + cl(6, 4)*u2z
    m(5, 2) = cl(4, 5)*u2x + cl(5, 5)*u2y + cl(6, 5)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
    m(3, 3) = cl(7, 3)*u2x + cl(8, 3)*u2y + cl(9, 3)*u2z
    m(4, 3) = cl(7, 4)*u2x + cl(8, 4)*u2y + cl(9, 4)*u2z
    m(5, 3) = cl(7, 5)*u2x + cl(8, 5)*u2y + cl(9, 5)*u2z
end function
pure function coeff9x6(cl, u) result(m)
    real(r8), dimension(9, 6), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(6, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(3, 1) = cl(1, 3)*u2x + cl(2, 3)*u2y + cl(3, 3)*u2z
    m(4, 1) = cl(1, 4)*u2x + cl(2, 4)*u2y + cl(3, 4)*u2z
    m(5, 1) = cl(1, 5)*u2x + cl(2, 5)*u2y + cl(3, 5)*u2z
    m(6, 1) = cl(1, 6)*u2x + cl(2, 6)*u2y + cl(3, 6)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(3, 2) = cl(4, 3)*u2x + cl(5, 3)*u2y + cl(6, 3)*u2z
    m(4, 2) = cl(4, 4)*u2x + cl(5, 4)*u2y + cl(6, 4)*u2z
    m(5, 2) = cl(4, 5)*u2x + cl(5, 5)*u2y + cl(6, 5)*u2z
    m(6, 2) = cl(4, 6)*u2x + cl(5, 6)*u2y + cl(6, 6)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
    m(3, 3) = cl(7, 3)*u2x + cl(8, 3)*u2y + cl(9, 3)*u2z
    m(4, 3) = cl(7, 4)*u2x + cl(8, 4)*u2y + cl(9, 4)*u2z
    m(5, 3) = cl(7, 5)*u2x + cl(8, 5)*u2y + cl(9, 5)*u2z
    m(6, 3) = cl(7, 6)*u2x + cl(8, 6)*u2y + cl(9, 6)*u2z
end function
pure function coeff9x7(cl, u) result(m)
    real(r8), dimension(9, 7), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(7, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(3, 1) = cl(1, 3)*u2x + cl(2, 3)*u2y + cl(3, 3)*u2z
    m(4, 1) = cl(1, 4)*u2x + cl(2, 4)*u2y + cl(3, 4)*u2z
    m(5, 1) = cl(1, 5)*u2x + cl(2, 5)*u2y + cl(3, 5)*u2z
    m(6, 1) = cl(1, 6)*u2x + cl(2, 6)*u2y + cl(3, 6)*u2z
    m(7, 1) = cl(1, 7)*u2x + cl(2, 7)*u2y + cl(3, 7)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(3, 2) = cl(4, 3)*u2x + cl(5, 3)*u2y + cl(6, 3)*u2z
    m(4, 2) = cl(4, 4)*u2x + cl(5, 4)*u2y + cl(6, 4)*u2z
    m(5, 2) = cl(4, 5)*u2x + cl(5, 5)*u2y + cl(6, 5)*u2z
    m(6, 2) = cl(4, 6)*u2x + cl(5, 6)*u2y + cl(6, 6)*u2z
    m(7, 2) = cl(4, 7)*u2x + cl(5, 7)*u2y + cl(6, 7)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
    m(3, 3) = cl(7, 3)*u2x + cl(8, 3)*u2y + cl(9, 3)*u2z
    m(4, 3) = cl(7, 4)*u2x + cl(8, 4)*u2y + cl(9, 4)*u2z
    m(5, 3) = cl(7, 5)*u2x + cl(8, 5)*u2y + cl(9, 5)*u2z
    m(6, 3) = cl(7, 6)*u2x + cl(8, 6)*u2y + cl(9, 6)*u2z
    m(7, 3) = cl(7, 7)*u2x + cl(8, 7)*u2y + cl(9, 7)*u2z
end function
pure function coeff9x8(cl, u) result(m)
    real(r8), dimension(9, 8), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(8, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(3, 1) = cl(1, 3)*u2x + cl(2, 3)*u2y + cl(3, 3)*u2z
    m(4, 1) = cl(1, 4)*u2x + cl(2, 4)*u2y + cl(3, 4)*u2z
    m(5, 1) = cl(1, 5)*u2x + cl(2, 5)*u2y + cl(3, 5)*u2z
    m(6, 1) = cl(1, 6)*u2x + cl(2, 6)*u2y + cl(3, 6)*u2z
    m(7, 1) = cl(1, 7)*u2x + cl(2, 7)*u2y + cl(3, 7)*u2z
    m(8, 1) = cl(1, 8)*u2x + cl(2, 8)*u2y + cl(3, 8)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(3, 2) = cl(4, 3)*u2x + cl(5, 3)*u2y + cl(6, 3)*u2z
    m(4, 2) = cl(4, 4)*u2x + cl(5, 4)*u2y + cl(6, 4)*u2z
    m(5, 2) = cl(4, 5)*u2x + cl(5, 5)*u2y + cl(6, 5)*u2z
    m(6, 2) = cl(4, 6)*u2x + cl(5, 6)*u2y + cl(6, 6)*u2z
    m(7, 2) = cl(4, 7)*u2x + cl(5, 7)*u2y + cl(6, 7)*u2z
    m(8, 2) = cl(4, 8)*u2x + cl(5, 8)*u2y + cl(6, 8)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
    m(3, 3) = cl(7, 3)*u2x + cl(8, 3)*u2y + cl(9, 3)*u2z
    m(4, 3) = cl(7, 4)*u2x + cl(8, 4)*u2y + cl(9, 4)*u2z
    m(5, 3) = cl(7, 5)*u2x + cl(8, 5)*u2y + cl(9, 5)*u2z
    m(6, 3) = cl(7, 6)*u2x + cl(8, 6)*u2y + cl(9, 6)*u2z
    m(7, 3) = cl(7, 7)*u2x + cl(8, 7)*u2y + cl(9, 7)*u2z
    m(8, 3) = cl(7, 8)*u2x + cl(8, 8)*u2y + cl(9, 8)*u2z
end function
pure function coeff9x9(cl, u) result(m)
    real(r8), dimension(9, 9), intent(in) :: cl
    real(r8), dimension(3), intent(in) :: u
    real(r8), dimension(9, 3) :: m
    real(r8) :: u2x, u2y, u2z
    u2x = u(1); u2y = u(2); u2z = u(3)
    m(1, 1) = cl(1, 1)*u2x + cl(2, 1)*u2y + cl(3, 1)*u2z
    m(2, 1) = cl(1, 2)*u2x + cl(2, 2)*u2y + cl(3, 2)*u2z
    m(3, 1) = cl(1, 3)*u2x + cl(2, 3)*u2y + cl(3, 3)*u2z
    m(4, 1) = cl(1, 4)*u2x + cl(2, 4)*u2y + cl(3, 4)*u2z
    m(5, 1) = cl(1, 5)*u2x + cl(2, 5)*u2y + cl(3, 5)*u2z
    m(6, 1) = cl(1, 6)*u2x + cl(2, 6)*u2y + cl(3, 6)*u2z
    m(7, 1) = cl(1, 7)*u2x + cl(2, 7)*u2y + cl(3, 7)*u2z
    m(8, 1) = cl(1, 8)*u2x + cl(2, 8)*u2y + cl(3, 8)*u2z
    m(9, 1) = cl(1, 9)*u2x + cl(2, 9)*u2y + cl(3, 9)*u2z
    m(1, 2) = cl(4, 1)*u2x + cl(5, 1)*u2y + cl(6, 1)*u2z
    m(2, 2) = cl(4, 2)*u2x + cl(5, 2)*u2y + cl(6, 2)*u2z
    m(3, 2) = cl(4, 3)*u2x + cl(5, 3)*u2y + cl(6, 3)*u2z
    m(4, 2) = cl(4, 4)*u2x + cl(5, 4)*u2y + cl(6, 4)*u2z
    m(5, 2) = cl(4, 5)*u2x + cl(5, 5)*u2y + cl(6, 5)*u2z
    m(6, 2) = cl(4, 6)*u2x + cl(5, 6)*u2y + cl(6, 6)*u2z
    m(7, 2) = cl(4, 7)*u2x + cl(5, 7)*u2y + cl(6, 7)*u2z
    m(8, 2) = cl(4, 8)*u2x + cl(5, 8)*u2y + cl(6, 8)*u2z
    m(9, 2) = cl(4, 9)*u2x + cl(5, 9)*u2y + cl(6, 9)*u2z
    m(1, 3) = cl(7, 1)*u2x + cl(8, 1)*u2y + cl(9, 1)*u2z
    m(2, 3) = cl(7, 2)*u2x + cl(8, 2)*u2y + cl(9, 2)*u2z
    m(3, 3) = cl(7, 3)*u2x + cl(8, 3)*u2y + cl(9, 3)*u2z
    m(4, 3) = cl(7, 4)*u2x + cl(8, 4)*u2y + cl(9, 4)*u2z
    m(5, 3) = cl(7, 5)*u2x + cl(8, 5)*u2y + cl(9, 5)*u2z
    m(6, 3) = cl(7, 6)*u2x + cl(8, 6)*u2y + cl(9, 6)*u2z
    m(7, 3) = cl(7, 7)*u2x + cl(8, 7)*u2y + cl(9, 7)*u2z
    m(8, 3) = cl(7, 8)*u2x + cl(8, 8)*u2y + cl(9, 8)*u2z
    m(9, 3) = cl(7, 9)*u2x + cl(8, 9)*u2y + cl(9, 9)*u2z
end function

! ! Rewrite the equations into a coefficient matrix times a vector.
! pure function coeff3x1(m,u) result(C)
!     real(r8), dimension(3,1), intent(in) :: m
!     real(r8), dimension(3), intent(in) :: u
!     real(r8), dimension(1,3) :: C
!     C(1,1)=m(1,1)*u(1)
!     C(1,2)=m(2,1)*u(2)
!     C(1,3)=m(3,1)*u(3)
! end function
! pure function coeff3x2(m,u) result(C)
!     real(r8), dimension(3,2), intent(in) :: m
!     real(r8), dimension(3), intent(in) :: u
!     real(r8), dimension(2,3) :: C
!     C(1,1)=m(1,1)*u(1)
!     C(1,2)=m(2,1)*u(2)
!     C(1,3)=m(3,1)*u(3)
!     C(2,1)=m(1,2)*u(1)
!     C(2,2)=m(2,2)*u(2)
!     C(2,3)=m(3,2)*u(3)
! end function
! pure function coeff3x3(m,u) result(C)
!     real(r8), dimension(3,3), intent(in) :: m
!     real(r8), dimension(3), intent(in) :: u
!     real(r8), dimension(3,3) :: C
!     C(1,1)=m(1,1)*u(1)
!     C(1,2)=m(2,1)*u(2)
!     C(1,3)=m(3,1)*u(3)
!     C(2,1)=m(1,2)*u(1)
!     C(2,2)=m(2,2)*u(2)
!     C(2,3)=m(3,2)*u(3)
!     C(3,1)=m(1,2)*u(1)
!     C(3,2)=m(2,2)*u(2)
!     C(3,3)=m(3,2)*u(3)
! end function

! pure function lcoeff3x1(cm,u) result(m)
!     real(r8), dimension(3,1), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u
!     real(r8), dimension(1) :: m
!     real(r8) :: ux,uy,uz
!     ux=u(1); uy=u(2); uz=u(3);
!     m(1)=cm(1,1)*ux + cm(2,1)*uy + cm(3,1)*uz
! end function
!
! pure function lcoeff3x2(cm,u) result(m)
!     real(r8), dimension(3,2), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u
!     real(r8), dimension(2) :: m
!     real(r8) :: ux,uy,uz
!     ux=u(1); uy=u(2); uz=u(3);
!     m(1)=cm(1,1)*ux + cm(2,1)*uy + cm(3,1)*uz
!     m(2)=cm(1,2)*ux + cm(2,2)*uy + cm(3,2)*uz
! end function
!
! pure function lcoeff3x3(cm,u) result(m)
!     real(r8), dimension(3,3), intent(in) :: cm
!     real(r8), dimension(3), intent(in) :: u
!     real(r8), dimension(3) :: m
!     real(r8) :: ux,uy,uz
!     ux=u(1); uy=u(2); uz=u(3);
!     m(1)=cm(1,1)*ux + cm(2,1)*uy + cm(3,1)*uz
!     m(2)=cm(1,2)*ux + cm(2,2)*uy + cm(3,2)*uz
!     m(3)=cm(1,3)*ux + cm(2,3)*uy + cm(3,3)*uz
! end function

end submodule
