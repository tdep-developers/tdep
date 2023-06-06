submodule(type_forcemap) type_forcemap_returntensors
!!
!! Converts from the irreducible representation to a meaningful representation
!!
use konstanter, only: lo_pi
use type_symmetryoperation, only: lo_expandoperation_triplet
use gottochblandat, only: lo_unflatten_3tensor, lo_unflatten_4tensor
implicit none

contains

!> get the first order force constant from the irreducible representation
module subroutine get_firstorder_forceconstant(map, p, fc)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> firstorder forceconstant
    type(lo_forceconstant_firstorder), intent(out) :: fc

    integer :: a1, unsh, unop
    real(r8), dimension(3) :: lfc

    fc%na = p%na
    allocate (fc%atom(fc%na))
    do a1 = 1, fc%na
        fc%atom(a1)%i1 = a1
        fc%atom(a1)%v1 = p%rcart(:, a1)
        fc%atom(a1)%lv1 = 0.0_r8
        unsh = map%xuc%fc_singlet(a1)%irreducible_shell
        unop = map%xuc%fc_singlet(a1)%operation_from_shell
        if (map%fc_singlet_shell(unsh)%nx .gt. 0) then
            lfc = matmul(map%fc_singlet_shell(unsh)%coeff, map%xuc%x_fc_singlet(map%fc_singlet_shell(unsh)%ind_global))
            fc%atom(a1)%m = lo_chop(matmul(map%op_singlet(unop)%m3, lfc), lo_sqtol)
        else
            fc%atom(a1)%m = 0.0_r8
        end if
    end do
end subroutine

!> get the secondorder forceconstants from the irreducible representation
module subroutine get_secondorder_forceconstant(map, p, fc, mem, verbosity)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> secondorder forceconstant
    type(lo_forceconstant_secondorder), intent(out) :: fc
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! Set up the structure, make sure everything is nothing.
    init: block
        integer :: a1, i
        ! Just init the basic stuff?
        fc%na = 0
        fc%cutoff = -1.0_r8
        fc%polar = .false.
        fc%npairshells = 0
        fc%npairop = 0
        fc%nifc = 0
        fc%nconstraints = 0
        fc%loto%correctiontype = 0
        fc%loto%eps = 0.0_r8
        fc%loto%nx_Z = 0
        ! Allocate space for the atoms and pairs, and make sure everything is nothing.
        fc%na = p%na
        allocate (fc%atom(fc%na))
        do a1 = 1, fc%na
            fc%atom(a1)%n = 0
        end do
        do i = 1, map%xuc%n_fc_pair
            a1 = map%xuc%fc_pair(i)%i1
            fc%atom(a1)%n = fc%atom(a1)%n + 1
        end do
        do a1 = 1, fc%na
            allocate (fc%atom(a1)%pair(fc%atom(a1)%n))
            do i = 1, fc%atom(a1)%n
                fc%atom(a1)%pair(i)%i1 = 0
                fc%atom(a1)%pair(i)%i2 = 0
                fc%atom(a1)%pair(i)%lv1 = 0.0_r8
                fc%atom(a1)%pair(i)%lv2 = 0.0_r8
                fc%atom(a1)%pair(i)%r = 0.0_r8
                fc%atom(a1)%pair(i)%m = 0.0_r8
                fc%atom(a1)%pair(i)%weight = 0.0_r8
                fc%atom(a1)%pair(i)%irreducible_shell = 0
                fc%atom(a1)%pair(i)%irreducible_operation = 0
            end do
            fc%atom(a1)%n = 0
        end do
    end block init

    ! Now store things
    storebasic: block
        real(r8), dimension(9) :: lfc, dv
        real(r8), dimension(3) :: r
        integer :: a1, a2, i, l, ii, jj, ipair
        integer :: unsh, unop

        do ipair = 1, map%xuc%n_fc_pair
            a1 = map%xuc%fc_pair(ipair)%i1
            a2 = map%xuc%fc_pair(ipair)%i2
            fc%atom(a1)%n = fc%atom(a1)%n + 1
            i = fc%atom(a1)%n
            ! What shell is this?
            unsh = map%xuc%fc_pair(ipair)%irreducible_shell
            unop = map%xuc%fc_pair(ipair)%operation_from_shell
            ! Fix the actual forceconstant
            if (map%xuc%fc_pair(ipair)%selfterm) then
                ! set self-terms to zero
                fc%atom(a1)%pair(i)%m = 0.0_r8
            else
                call lo_gemv(map%fc_pair_shell(unsh)%coeff, map%xuc%x_fc_pair(map%fc_pair_shell(unsh)%ind_global), dv)
                call lo_gemv(map%op_pair(unop)%sotr, dv, lfc)
                l = 0
                do ii = 1, 3
                do jj = 1, 3
                    l = l + 1
                    fc%atom(a1)%pair(i)%m(ii, jj) = lfc(l)
                end do
                end do
            end if
            ! A bunch of extra stuff
            fc%atom(a1)%pair(i)%i1 = a1
            fc%atom(a1)%pair(i)%i2 = a2
            r = matmul(p%latticevectors, map%xuc%fc_pair(ipair)%flv)
            fc%atom(a1)%pair(i)%lv1 = 0.0_r8
            fc%atom(a1)%pair(i)%lv2 = lo_chop(r, 1E-13_r8)
            fc%atom(a1)%pair(i)%r = lo_chop(r + p%rcart(:, a2) - p%rcart(:, a1), 1E-13_r8)
            ! and the cutoff
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%pair(i)%r) + 10*lo_tol)
        end do
        ! fix the acoustic sum rules
        call fc%setsumtozero()
        ! set tiny things to zero
        do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            fc%atom(a1)%pair(i)%m = lo_chop(fc%atom(a1)%pair(i)%m, lo_sqtol)
        end do
        end do
    end block storebasic

    ! ! get the coordination shell information thing
    ! getshells: block
    !     real(r8), dimension(:), allocatable :: shellrad
    !     integer, dimension(:), allocatable :: shellind,dum
    !     integer :: a1,i,j
    !
    !     ! get indices so that I can sort the shells by radius
    !     lo_allocate(shellrad(map%n_fc_pair_shell))
    !     lo_allocate(shellind(map%n_fc_pair_shell))
    !     lo_allocate(dum(map%n_fc_pair_shell))
    !     do a1=1,map%n_atom_uc
    !     do i=1,map%uc(a1)%n_fc_pair
    !         j=map%uc(a1)%fc_pair(i)%irreducible_shell
    !         shellrad(j)=lo_chop( norm2(fc%atom(a1)%pair(i)%r) , lo_sqtol )
    !     enddo
    !     enddo
    !     ! Invert the index list
    !     call qsort(shellrad,dum)
    !     do i=1,map%n_fc_pair_shell
    !     do j=1,map%n_fc_pair_shell
    !         if ( dum(i) .eq. j ) shellind(j)=i
    !     enddo
    !     enddo
    !
    !     ! create some space for the shells
    !     fc%npairshells=map%n_fc_pair_shell
    !     allocate(fc%pairshell( fc%npairshells ))
    !     do i=1,fc%npairshells
    !         fc%pairshell(i)%n=0
    !     enddo
    !
    !     ! count number of pairs per shell
    !     do a1=1,map%n_atom_uc
    !     do i=1,map%uc(a1)%n_fc_pair
    !         j=shellind( map%uc(a1)%fc_pair(i)%irreducible_shell )
    !         fc%pairshell(j)%n=fc%pairshell(j)%n+1
    !     enddo
    !     enddo
    !
    !     ! Make space for indices, vectors and stuff
    !     do i=1,fc%npairshells
    !         lo_allocate(fc%pairshell(i)%vec( 3, fc%pairshell(i)%n ))
    !         lo_allocate(fc%pairshell(i)%atind( fc%pairshell(i)%n ))
    !         lo_allocate(fc%pairshell(i)%pairind( fc%pairshell(i)%n ))
    !         fc%pairshell(i)%n=0
    !     enddo
    !
    !     ! Now fill the shells with atoms and stuff
    !     do a1=1,map%n_atom_uc
    !     do i=1,map%uc(a1)%n_fc_pair
    !         j=shellind( map%uc(a1)%fc_pair(i)%irreducible_shell )
    !         fc%pairshell(j)%n=fc%pairshell(j)%n+1
    !         fc%pairshell(j)%vec( :,fc%pairshell(j)%n )=lo_chop(fc%atom(a1)%pair(i)%r,lo_sqtol)
    !         fc%pairshell(j)%atind( fc%pairshell(j)%n )=a1
    !         fc%pairshell(j)%pairind( fc%pairshell(j)%n )=i
    !         fc%pairshell(j)%rad=lo_chop( norm2(fc%atom(a1)%pair(i)%r) , lo_sqtol)
    !         fc%pairshell(j)%norm=lo_frobnorm( fc%atom(a1)%pair(i)%m )
    !     enddo
    !     enddo
    !     lo_deallocate(shellrad)
    !     lo_deallocate(shellind)
    !     lo_deallocate(dum)
    ! endblock getshells

    ! store the polar information here
    if (map%polar .gt. 0) then
        polarstuff: block
            ! Yup, this material has polar corrections
            fc%polar = .true.
            ! First set the dielectric constant
            fc%loto%eps = lo_unflatten_2tensor(matmul(map%eps_global_shell%coeff, map%xuc%x_eps_global))
            fc%loto%eps = lo_chop(fc%loto%eps, lo_sqtol)
            ! Now there could be a rare case where we think it is polar but it is actually not.
            if (map%xuc%nx_Z_singlet .le. 0) then
                fc%loto%nx_Z = 0
                fc%loto%correctiontype = 0
                fc%polar = .false.
                allocate (fc%loto%born_effective_charges(3, 3, fc%na))
                allocate (fc%loto%born_onsite_correction(3, 3, fc%na))
                fc%loto%born_effective_charges = 0.0_r8
                fc%loto%born_onsite_correction = 0.0_r8
            else
                ! Whis is the perfectly normal case, set the polar things
                fc%loto%nx_Z = map%xuc%nx_Z_singlet
                allocate (fc%loto%x_Z(fc%loto%nx_Z))
                allocate (fc%loto%coeff_Z(fc%na*9, fc%loto%nx_Z))
                fc%loto%x_Z = map%xuc%x_Z_singlet
                fc%loto%coeff_Z = 0.0_r8
                ! confusing, I acknowledge, but it should be here.
                call lo_coeffmatrix_unitcell_Z_singlet(map, fc%loto%coeff_Z)
                ! What kind of correction are we talking about
                fc%loto%correctiontype = map%polarcorrectiontype
                ! space for charges and corrections
                allocate (fc%loto%born_effective_charges(3, 3, fc%na))
                allocate (fc%loto%born_onsite_correction(3, 3, fc%na))
                fc%loto%born_effective_charges = 0.0_r8
                fc%loto%born_onsite_correction = 0.0_r8
                ! Set the charges, ewald parameters, and minimize Hermitian error and get the on-site corrections.
                call fc%set_ewald_and_enforce_borncharge_hermiticity(p, mem, verbosity)
            end if
        end block polarstuff
    end if
end subroutine

!> get the thirdorder forceconstants from the irreducible representation
module subroutine get_thirdorder_forceconstant(map, p, fc)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> secondorder forceconstant
    type(lo_forceconstant_thirdorder), intent(out) :: fc

    integer :: a1, a2, a3, i, j, l
    integer :: ii, jj, kk, unsh, unop
    real(r8), dimension(3, 3, 3) :: m0
    real(r8), dimension(27) :: lfc, dv
    real(r8), dimension(3) :: v2, v3

    ! Make some space
    fc%na = p%na
    fc%cutoff = 0.0_r8
    allocate (fc%atom(fc%na))
    do a1 = 1, fc%na
        fc%atom(a1)%n = 0
    end do
    do i = 1, map%xuc%n_fc_triplet
        a1 = map%xuc%fc_triplet(i)%i1
        fc%atom(a1)%n = fc%atom(a1)%n + 1
    end do
    do a1 = 1, fc%na
        allocate (fc%atom(a1)%triplet(fc%atom(a1)%n))
        do i = 1, fc%atom(a1)%n
            fc%atom(a1)%triplet(i)%i1 = 0
            fc%atom(a1)%triplet(i)%i2 = 0
            fc%atom(a1)%triplet(i)%i3 = 0
            fc%atom(a1)%triplet(i)%lv1 = 0.0_r8
            fc%atom(a1)%triplet(i)%lv2 = 0.0_r8
            fc%atom(a1)%triplet(i)%lv3 = 0.0_r8
            fc%atom(a1)%triplet(i)%rv1 = 0.0_r8
            fc%atom(a1)%triplet(i)%rv2 = 0.0_r8
            fc%atom(a1)%triplet(i)%rv3 = 0.0_r8
            fc%atom(a1)%triplet(i)%m = 0.0_r8
            fc%atom(a1)%triplet(i)%mwm = 0.0_r8
        end do
        fc%atom(a1)%n = 0
    end do

    ! Store things
    do i = 1, map%xuc%n_fc_triplet
        a1 = map%xuc%fc_triplet(i)%i1
        a2 = map%xuc%fc_triplet(i)%i2
        a3 = map%xuc%fc_triplet(i)%i3
        unsh = map%xuc%fc_triplet(i)%irreducible_shell
        unop = map%xuc%fc_triplet(i)%operation_from_shell
        fc%atom(a1)%n = fc%atom(a1)%n + 1
        j = fc%atom(a1)%n
        fc%atom(a1)%triplet(j)%i1 = map%xuc%fc_triplet(i)%i1
        fc%atom(a1)%triplet(j)%i2 = map%xuc%fc_triplet(i)%i2
        fc%atom(a1)%triplet(j)%i3 = map%xuc%fc_triplet(i)%i3
        v2 = matmul(p%latticevectors, map%xuc%fc_triplet(i)%flv2)
        v3 = matmul(p%latticevectors, map%xuc%fc_triplet(i)%flv3)
        fc%atom(a1)%triplet(j)%lv1 = 0.0_r8
        fc%atom(a1)%triplet(j)%lv2 = v2
        fc%atom(a1)%triplet(j)%lv3 = v3
        v2 = v2 - p%rcart(:, a1) + p%rcart(:, a2)
        v3 = v3 - p%rcart(:, a1) + p%rcart(:, a3)
        fc%atom(a1)%triplet(j)%rv1 = 0.0_r8
        fc%atom(a1)%triplet(j)%rv2 = v2
        fc%atom(a1)%triplet(j)%rv3 = v3
        if (map%xuc%fc_triplet(i)%selfterm) then
            ! set self-terms to zero
            fc%atom(a1)%triplet(j)%m = 0.0_r8
        else
            ! get the normal ones
            call lo_gemv(map%fc_triplet_shell(unsh)%coeff, map%xuc%x_fc_triplet(map%fc_triplet_shell(unsh)%ind_global), dv)
            call lo_gemv(map%op_triplet(unop)%sotr, dv, lfc)
            l = 0
            do ii = 1, 3
            do jj = 1, 3
            do kk = 1, 3
                l = l + 1
                fc%atom(a1)%triplet(j)%m(ii, jj, kk) = lfc(l)
            end do
            end do
            end do
        end if
        fc%cutoff = max(fc%cutoff, norm2(v2) + 10*lo_tol)
        fc%cutoff = max(fc%cutoff, norm2(v3) + 10*lo_tol)
        fc%cutoff = max(fc%cutoff, norm2(v3 - v2) + 10*lo_tol)
    end do
    ! fix the acoustic sum rules
    do a1 = 1, fc%na
        m0 = 0.0_r8
        j = 0
        do i = 1, fc%atom(a1)%n
            if (lo_sqnorm(fc%atom(a1)%triplet(i)%rv2) + lo_sqnorm(fc%atom(a1)%triplet(i)%rv3) .gt. lo_sqtol) then
                m0 = m0 - fc%atom(a1)%triplet(i)%m
            else
                j = i
            end if
        end do
        if (j .gt. 0) then
            fc%atom(a1)%triplet(j)%m = m0
        else
            call lo_stop_gracefully(['Could not find self-term'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
    end do
    ! set tiny things to zero
    do a1 = 1, fc%na
    do i = 1, fc%atom(a1)%n
        fc%atom(a1)%triplet(i)%m = lo_chop(fc%atom(a1)%triplet(i)%m, lo_sqtol)
        a2 = fc%atom(a1)%triplet(i)%i2
        a3 = fc%atom(a1)%triplet(i)%i3
        fc%atom(a1)%triplet(i)%mwm = fc%atom(a1)%triplet(i)%m*p%invsqrtmass(a1)*p%invsqrtmass(a2)*p%invsqrtmass(a3)
    end do
    end do
end subroutine

!> get the fourth order forceconstants from the irreducible representation
module subroutine get_fourthorder_forceconstant(map, p, fc)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> secondorder forceconstant
    type(lo_forceconstant_fourthorder), intent(out) :: fc

    integer :: a1, a2, a3, a4, i, j, l
    integer :: ii, jj, kk, ll, unsh, unop
    real(r8), dimension(3, 3, 3, 3) :: m0
    real(r8), dimension(81) :: lfc, dv
    real(r8), dimension(3) :: v2, v3, v4

    ! Make some space
    fc%na = p%na
    fc%cutoff = 0.0_r8
    allocate (fc%atom(fc%na))
    do a1 = 1, fc%na
        fc%atom(a1)%n = 0
    end do
    do i = 1, map%xuc%n_fc_quartet
        a1 = map%xuc%fc_quartet(i)%i1
        fc%atom(a1)%n = fc%atom(a1)%n + 1
    end do
    do a1 = 1, fc%na
        allocate (fc%atom(a1)%quartet(fc%atom(a1)%n))
        do i = 1, fc%atom(a1)%n
            fc%atom(a1)%quartet(i)%i1 = 0
            fc%atom(a1)%quartet(i)%i2 = 0
            fc%atom(a1)%quartet(i)%i3 = 0
            fc%atom(a1)%quartet(i)%i4 = 0
            fc%atom(a1)%quartet(i)%lv1 = 0.0_r8
            fc%atom(a1)%quartet(i)%lv2 = 0.0_r8
            fc%atom(a1)%quartet(i)%lv3 = 0.0_r8
            fc%atom(a1)%quartet(i)%lv4 = 0.0_r8
            fc%atom(a1)%quartet(i)%rv1 = 0.0_r8
            fc%atom(a1)%quartet(i)%rv2 = 0.0_r8
            fc%atom(a1)%quartet(i)%rv3 = 0.0_r8
            fc%atom(a1)%quartet(i)%rv4 = 0.0_r8
            fc%atom(a1)%quartet(i)%m = 0.0_r8
            fc%atom(a1)%quartet(i)%mwm = 0.0_r8
        end do
        fc%atom(a1)%n = 0
    end do

    ! Store things
    do i = 1, map%xuc%n_fc_quartet
        a1 = map%xuc%fc_quartet(i)%i1
        a2 = map%xuc%fc_quartet(i)%i2
        a3 = map%xuc%fc_quartet(i)%i3
        a4 = map%xuc%fc_quartet(i)%i4
        unsh = map%xuc%fc_quartet(i)%irreducible_shell
        unop = map%xuc%fc_quartet(i)%operation_from_shell
        fc%atom(a1)%n = fc%atom(a1)%n + 1
        j = fc%atom(a1)%n
        fc%atom(a1)%quartet(j)%i1 = map%xuc%fc_quartet(i)%i1
        fc%atom(a1)%quartet(j)%i2 = map%xuc%fc_quartet(i)%i2
        fc%atom(a1)%quartet(j)%i3 = map%xuc%fc_quartet(i)%i3
        fc%atom(a1)%quartet(j)%i4 = map%xuc%fc_quartet(i)%i4
        v2 = matmul(p%latticevectors, map%xuc%fc_quartet(i)%flv2)
        v3 = matmul(p%latticevectors, map%xuc%fc_quartet(i)%flv3)
        v4 = matmul(p%latticevectors, map%xuc%fc_quartet(i)%flv4)
        fc%atom(a1)%quartet(j)%lv1 = 0.0_r8
        fc%atom(a1)%quartet(j)%lv2 = v2
        fc%atom(a1)%quartet(j)%lv3 = v3
        fc%atom(a1)%quartet(j)%lv4 = v4
        v2 = v2 - p%rcart(:, a1) + p%rcart(:, a2)
        v3 = v3 - p%rcart(:, a1) + p%rcart(:, a3)
        v4 = v4 - p%rcart(:, a1) + p%rcart(:, a4)
        fc%atom(a1)%quartet(j)%rv1 = 0.0_r8
        fc%atom(a1)%quartet(j)%rv2 = v2
        fc%atom(a1)%quartet(j)%rv3 = v3
        fc%atom(a1)%quartet(j)%rv4 = v4
        if (map%xuc%fc_quartet(i)%selfterm) then
            ! set self-terms to zero
            fc%atom(a1)%quartet(j)%m = 0.0_r8
        else
            ! get the normal ones
            call lo_gemv(map%fc_quartet_shell(unsh)%coeff, map%xuc%x_fc_quartet(map%fc_quartet_shell(unsh)%ind_global), dv)
            call lo_gemv(map%op_quartet(unop)%sotr, dv, lfc)
            l = 0
            do ii = 1, 3
            do jj = 1, 3
            do kk = 1, 3
            do ll = 1, 3
                l = l + 1
                fc%atom(a1)%quartet(j)%m(ii, jj, kk, ll) = lfc(l)
            end do
            end do
            end do
            end do
        end if
        fc%cutoff = max(fc%cutoff, norm2(v2) + 10*lo_tol)
        fc%cutoff = max(fc%cutoff, norm2(v3) + 10*lo_tol)
        fc%cutoff = max(fc%cutoff, norm2(v4) + 10*lo_tol)
        fc%cutoff = max(fc%cutoff, norm2(v3 - v2) + 10*lo_tol)
        fc%cutoff = max(fc%cutoff, norm2(v4 - v2) + 10*lo_tol)
        fc%cutoff = max(fc%cutoff, norm2(v4 - v3) + 10*lo_tol)
    end do
    ! fix the acoustic sum rules
    do a1 = 1, fc%na
        m0 = 0.0_r8
        j = 0
        do i = 1, fc%atom(a1)%n
            if (lo_sqnorm(fc%atom(a1)%quartet(i)%rv2) + lo_sqnorm(fc%atom(a1)%quartet(i)%rv3) + lo_sqnorm(fc%atom(a1)%quartet(i)%rv4) .gt. lo_sqtol) then
                m0 = m0 - fc%atom(a1)%quartet(i)%m
            else
                j = i
            end if
        end do
        if (j .gt. 0) then
            fc%atom(a1)%quartet(j)%m = m0
        else
            call lo_stop_gracefully(['Could not find self-term'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
    end do
    ! set tiny things to zero
    do a1 = 1, fc%na
    do i = 1, fc%atom(a1)%n
        fc%atom(a1)%quartet(i)%m = lo_chop(fc%atom(a1)%quartet(i)%m, lo_sqtol)
        a2 = fc%atom(a1)%quartet(i)%i2
        a3 = fc%atom(a1)%quartet(i)%i3
        a4 = fc%atom(a1)%quartet(i)%i4
        fc%atom(a1)%quartet(i)%mwm = fc%atom(a1)%quartet(i)%m*p%invsqrtmass(a1)*p%invsqrtmass(a2)*p%invsqrtmass(a3)*p%invsqrtmass(a4)
    end do
    end do
end subroutine

!> get the dielectric interactions to higher order
module subroutine get_dielectric_tensors(map, p, di)
    !> forcemap
    class(lo_forcemap), intent(in) :: map
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> dielectric interactions
    type(lo_dielectric_tensor), intent(out) :: di

    ! First set general things
    di%n_atom = map%n_atom_uc

    ! Set the global polarizability?
    di%eps_inf = lo_unflatten_2tensor(matmul(map%eps_global_shell%coeff, map%xuc%x_eps_global))
    di%eps_inf = lo_chop(di%eps_inf, lo_sqtol)

    if (map%xuc%nx_Z_singlet .gt. 0) then
        zsing: block
            real(r8), dimension(:, :), allocatable :: coeff
            real(r8), dimension(:), allocatable :: zf
            integer :: a1

            di%n_Z_singlet = p%na
            allocate (coeff(p%na*9, map%xuc%nx_Z_singlet))
            allocate (zf(p%na*9))
            coeff = 0.0_r8
            zf = 0.0_r8
            call lo_coeffmatrix_unitcell_Z_singlet(map, coeff)
            zf = matmul(coeff, map%xuc%x_Z_singlet)

            ! Set the Born charges
            allocate (di%Z_singlet(3, 3, p%na))
            di%Z_singlet = 0.0_r8
            do a1 = 1, p%na
                di%Z_singlet(:, :, a1) = lo_unflatten_2tensor(zf((a1 - 1)*9 + 1:a1*9))
            end do

            deallocate (coeff)
            deallocate (zf)
        end block zsing
    else
        ! No Born charges
        di%n_Z_singlet = 0
    end if

    if (map%have_Z_pair) then
        zpair: block
            real(r8), dimension(3, 3, 3) :: m0
            real(r8), dimension(27) :: vA, vB
            real(r8), dimension(3) :: v0
            real(r8), dimension(map%xuc%nx_Z_pair) :: x
            integer, dimension(map%xuc%nx_Z_pair) :: xind
            integer :: ipair, ish, iop, nx, a1, jpair

            ! populate and unfold tensors to the relevant pairs
            di%n_Z_pair = map%xuc%n_Z_pair
            allocate (di%Z_pair(di%n_Z_pair))
            do ipair = 1, di%n_Z_pair
                ish = map%xuc%Z_pair(ipair)%irreducible_shell
                iop = map%xuc%Z_pair(ipair)%operation_from_shell
                ! store simple things
                di%Z_pair(ipair)%a1 = map%xuc%Z_pair(ipair)%i1
                di%Z_pair(ipair)%a2 = map%xuc%Z_pair(ipair)%i2
                v0 = map%xuc%Z_pair(ipair)%flv
                v0 = matmul(p%latticevectors, v0)
                di%Z_pair(ipair)%lv = v0
                v0 = v0 + p%rcart(:, di%Z_pair(ipair)%a2) - p%rcart(:, di%Z_pair(ipair)%a1)
                di%Z_pair(ipair)%v = v0
                nx = map%Z_pair_shell(ish)%nx
                if (nx .eq. 0) then
                    di%Z_pair(ipair)%m = 0.0_r8
                else
                    xind(1:nx) = map%Z_pair_shell(ish)%ind_global
                    x(1:nx) = map%xuc%x_Z_pair(xind(1:nx))
                    ! convert to flattened tensor
                    call lo_gemv(map%Z_pair_shell(ish)%coeff, x(1:nx), vA)
                    ! rotate tensor
                    call lo_gemv(map%op_triplet(iop)%sotr, vA, vB)
                    ! unflatten
                    di%Z_pair(ipair)%m = lo_unflatten_3tensor(vB)
                    ! remove tiny things
                    di%Z_pair(ipair)%m = lo_chop(di%Z_pair(ipair)%m, 1E-12_r8)
                end if
            end do

            ! fix the self-terms
            do a1 = 1, p%na
                ! get a list of relevant pairs
                jpair = 0
                m0 = 0.0_r8
                do ipair = 1, di%n_Z_pair
                    if (di%Z_pair(ipair)%a1 .ne. a1) cycle

                    if (lo_sqnorm(di%Z_pair(ipair)%v) .lt. lo_sqtol) then
                        ! sanity test, there is only supposed to be one self-term
                        if (jpair .ne. 0) then
                            call lo_stop_gracefully(['Found more than one self-term, should be impossible'], lo_exitcode_symmetry, __FILE__, __LINE__)
                        else
                            jpair = ipair
                        end if
                    else
                        m0 = m0 + di%Z_pair(ipair)%m
                    end if
                end do
                ! add things together
                if (jpair .eq. 0) then
                    call lo_stop_gracefully(['Could not find self-term, should be impossible'], lo_exitcode_symmetry, __FILE__, __LINE__)
                else
                    di%Z_pair(jpair)%m = -lo_chop(m0, 1E-12_r8)
                end if
            end do
        end block zpair
    else
        ! nope, no pairs
        di%n_Z_pair = 0
    end if

    if (map%have_Z_triplet) then
        ztriplet: block
            real(r8), dimension(3, 3, 3, 3) :: m0
            real(r8), dimension(81) :: vA, vB
            real(r8), dimension(3) :: v2, v3
            real(r8), dimension(map%xuc%nx_Z_triplet) :: x
            integer, dimension(map%xuc%nx_Z_triplet) :: xind
            integer :: itriplet, ish, iop, nx, a1, jtriplet

            ! populate and unfold tensors to the relevant triplets
            di%n_Z_triplet = map%xuc%n_Z_triplet
            allocate (di%Z_triplet(di%n_Z_triplet))
            do itriplet = 1, di%n_Z_triplet
                ish = map%xuc%Z_triplet(itriplet)%irreducible_shell
                iop = map%xuc%Z_triplet(itriplet)%operation_from_shell
                ! store simple things
                di%Z_triplet(itriplet)%a1 = map%xuc%Z_triplet(itriplet)%i1
                di%Z_triplet(itriplet)%a2 = map%xuc%Z_triplet(itriplet)%i2
                di%Z_triplet(itriplet)%a3 = map%xuc%Z_triplet(itriplet)%i3
                v2 = map%xuc%Z_triplet(itriplet)%flv2
                v3 = map%xuc%Z_triplet(itriplet)%flv3
                v2 = matmul(p%latticevectors, v2)
                v3 = matmul(p%latticevectors, v3)
                di%Z_triplet(itriplet)%lv2 = v2
                di%Z_triplet(itriplet)%lv3 = v3
                v2 = v2 + p%rcart(:, di%Z_triplet(itriplet)%a2) - p%rcart(:, di%Z_triplet(itriplet)%a1)
                v3 = v3 + p%rcart(:, di%Z_triplet(itriplet)%a3) - p%rcart(:, di%Z_triplet(itriplet)%a1)
                di%Z_triplet(itriplet)%v2 = v2
                di%Z_triplet(itriplet)%v3 = v3
                nx = map%Z_triplet_shell(ish)%nx
                if (nx .eq. 0) then
                    di%Z_triplet(itriplet)%m = 0.0_r8
                else
                    xind(1:nx) = map%Z_triplet_shell(ish)%ind_global
                    x(1:nx) = map%xuc%x_Z_triplet(xind(1:nx))
                    ! convert to flattened tensor
                    call lo_gemv(map%Z_triplet_shell(ish)%coeff, x(1:nx), vA)
                    ! rotate tensor
                    call lo_gemv(map%op_quartet(iop)%sotr, vA, vB)
                    ! unflatten
                    di%Z_triplet(itriplet)%m = lo_unflatten_4tensor(vB)
                    ! remove tiny things
                    di%Z_triplet(itriplet)%m = lo_chop(di%Z_triplet(itriplet)%m, 1E-12_r8)
                end if
            end do

            ! fix the self-terms
            do a1 = 1, p%na
                ! get a list of relevant triplets
                jtriplet = 0
                m0 = 0.0_r8
                do itriplet = 1, di%n_Z_triplet
                    if (di%Z_triplet(itriplet)%a1 .ne. a1) cycle

                    if (lo_sqnorm(di%Z_triplet(itriplet)%v2) + lo_sqnorm(di%Z_triplet(itriplet)%v3) .lt. lo_sqtol) then
                        ! sanity test, there is only supposed to be one self-term
                        if (jtriplet .ne. 0) then
                            call lo_stop_gracefully(['Found more than one self-term, should be impossible'], lo_exitcode_symmetry, __FILE__, __LINE__)
                        else
                            jtriplet = itriplet
                        end if
                    else
                        m0 = m0 + di%Z_triplet(itriplet)%m
                    end if
                end do
                ! add things together
                if (jtriplet .eq. 0) then
                    call lo_stop_gracefully(['Could not find self-term, should be impossible'], lo_exitcode_symmetry, __FILE__, __LINE__)
                else
                    di%Z_triplet(jtriplet)%m = -lo_chop(m0, 1E-12_r8)
                end if
            end do
        end block ztriplet
    else
        ! nope, no pairs
        di%n_Z_triplet = 0
    end if

    if (map%have_eps_singlet) then
        esing: block
            real(r8), dimension(27, 27) :: rotm
            real(r8), dimension(27) :: vA, vB
            real(r8), dimension(map%xuc%nx_eps_singlet) :: x
            integer, dimension(map%xuc%nx_eps_singlet) :: xind
            integer :: a1, ish, iop, nx

            di%n_eps_singlet = map%xuc%n_eps_singlet
            allocate (di%eps_singlet(di%n_eps_singlet))

            do a1 = 1, p%na
                ish = map%xuc%eps_singlet(a1)%irreducible_shell
                iop = map%xuc%eps_singlet(a1)%operation_from_shell
                nx = map%eps_singlet_shell(ish)%nx
                if (nx .eq. 0) then
                    di%eps_singlet(a1)%m = 0.0_r8
                else
                    xind(1:nx) = map%eps_singlet_shell(ish)%ind_global
                    x(1:nx) = map%xuc%x_eps_singlet(xind(1:nx))
                    rotm = lo_expandoperation_triplet(map%op_singlet(iop)%m3)
                    call lo_gemv(map%eps_singlet_shell(ish)%coeff, x(1:nx), vA)
                    call lo_gemv(rotm, vA, vB)
                    di%eps_singlet(a1)%m = lo_unflatten_3tensor(vB)
                    di%eps_singlet(a1)%m = lo_chop(di%eps_singlet(a1)%m, 1E-13_r8)
                end if
            end do
        end block esing
    else
        ! none of these
        di%n_eps_singlet = 0
    end if

    if (map%have_eps_pair) then
        epspair: block
            real(r8), dimension(3, 3, 3, 3) :: m0
            real(r8), dimension(81) :: vA, vB
            real(r8), dimension(3) :: v0
            real(r8), dimension(map%xuc%nx_eps_pair) :: x
            integer, dimension(map%xuc%nx_eps_pair) :: xind
            integer :: ipair, ish, iop, nx, a1, jpair

            ! populate and unfold tensors to the relevant pairs
            di%n_eps_pair = map%xuc%n_eps_pair
            allocate (di%eps_pair(di%n_eps_pair))
            do ipair = 1, di%n_eps_pair
                ish = map%xuc%eps_pair(ipair)%irreducible_shell
                iop = map%xuc%eps_pair(ipair)%operation_from_shell
                ! store simple things
                di%eps_pair(ipair)%a1 = map%xuc%eps_pair(ipair)%i1
                di%eps_pair(ipair)%a2 = map%xuc%eps_pair(ipair)%i2
                v0 = map%xuc%eps_pair(ipair)%flv
                v0 = matmul(p%latticevectors, v0)
                di%eps_pair(ipair)%lv = v0
                v0 = v0 + p%rcart(:, di%eps_pair(ipair)%a2) - p%rcart(:, di%eps_pair(ipair)%a1)
                di%eps_pair(ipair)%v = v0
                nx = map%eps_pair_shell(ish)%nx
                if (nx .eq. 0) then
                    di%eps_pair(ipair)%m = 0.0_r8
                else
                    xind(1:nx) = map%eps_pair_shell(ish)%ind_global
                    x(1:nx) = map%xuc%x_eps_pair(xind(1:nx))
                    ! convert to flattened tensor
                    call lo_gemv(map%eps_pair_shell(ish)%coeff, x(1:nx), vA)
                    ! rotate tensor
                    call lo_gemv(map%op_quartet(iop)%sotr, vA, vB)
                    ! unflatten
                    di%eps_pair(ipair)%m = lo_unflatten_4tensor(vB)
                    ! remove tiny things
                    di%eps_pair(ipair)%m = lo_chop(di%eps_pair(ipair)%m, 1E-12_r8)
                end if
            end do

            ! fix the self-terms
            do a1 = 1, p%na
                ! get a list of relevant pairs
                jpair = 0
                m0 = 0.0_r8
                do ipair = 1, di%n_eps_pair
                    if (di%eps_pair(ipair)%a1 .ne. a1) cycle

                    if (lo_sqnorm(di%eps_pair(ipair)%v) .lt. lo_sqtol) then
                        ! sanity test, there is only supposed to be one self-term
                        if (jpair .ne. 0) then
                            call lo_stop_gracefully(['Found more than one self-term, should be impossible'], lo_exitcode_symmetry, __FILE__, __LINE__)
                        else
                            jpair = ipair
                        end if
                    else
                        m0 = m0 + di%eps_pair(ipair)%m
                    end if
                end do
                ! add things together
                if (jpair .eq. 0) then
                    call lo_stop_gracefully(['Could not find self-term, should be impossible'], lo_exitcode_symmetry, __FILE__, __LINE__)
                else
                    di%eps_pair(jpair)%m = -lo_chop(m0, 1E-12_r8)
                end if
            end do
        end block epspair
    else
        ! nope, no pairs
        di%n_eps_pair = 0
    end if
end subroutine

! !> get the secondorder magnetic interactions from the irreducible representation
! subroutine get_secondorder_jij( map,uc,jij )
!     !> the forcemap
!     class(lo_forcemap), intent(in) :: map
!     !> the crystalstructure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> the secondorder forceconstant
!     type(lo_jij_secondorder), intent(out) :: jij
!
!     write(*,*) 'FIXME MAGNETIC JIJ'
!     stop
! !
! !     ! Set up the structure, make sure everything is nothing.
! !     init: block
! !         integer :: a1,i
! !         ! Just init the basic stuff
! !         jij%na=0
! !         jij%verbosity=-1
! !         jij%cutoff=-1.0_r8
! !         jij%na=uc%na
! !         lo_allocate(jij%atom( jij%na ))
! !         do a1=1,jij%na
! !             jij%atom(a1)%m0=map%uc(a1)%J%m0
! !             jij%atom(a1)%m0dev=map%uc(a1)%J%m0dev
! !             jij%atom(a1)%n=map%uc(a1)%nmagpair
! !             lo_allocate(jij%atom(a1)%pair( jij%atom(a1)%n ))
! !             do i=1,jij%atom(a1)%n
! !                 jij%atom(a1)%pair(i)%i1=0
! !                 jij%atom(a1)%pair(i)%i2=0
! !                 jij%atom(a1)%pair(i)%lv1=0.0_r8
! !                 jij%atom(a1)%pair(i)%lv2=0.0_r8
! !                 jij%atom(a1)%pair(i)%r=0.0_r8
! !                 jij%atom(a1)%pair(i)%J=0.0_r8
! !                 jij%atom(a1)%pair(i)%irreducible_shell=0
! !                 jij%atom(a1)%pair(i)%irreducible_operation=0
! !             enddo
! !         enddo
! !     end block init
! !
! !     ! Now store things
! !     storebasic: block
! !         real(r8), dimension(9) :: ljij,ltij,dv !,lqij
! !         real(r8), dimension(3) :: r,lqij
! !         integer :: a1,a2,i,j,l,ii,jj
! !         integer :: unsh,unop
! !
! !         do a1=1,jij%na
! !             ! Insert dumping of self-term here
! !             do i=1,map%uc(a1)%nmagpair
! !                 ! What shell is this
! !                 unsh=map%uc(a1)%magpair(i)%irreducible_shell
! !                 unop=map%uc(a1)%magpair(i)%operation_from_shell
! !                 ! Fix the actual tensors
! !                 if ( map%magpairshell(unsh)%ntheta_jij .gt. 0 ) then
! !                     call lo_gemv( map%magpairshell(unsh)%coeff_jij,map%theta_magpair_jij(map%magpairshell(unsh)%ind_global_jij),dv )
! !                     call lo_gemv( map%pairop(unop)%sotr,dv,ljij )
! !                 else
! !                     ljij=0.0_r8
! !                 endif
! !                 if ( map%magpairshell(unsh)%ntheta_tij .gt. 0 ) then
! !                     call lo_gemv( map%magpairshell(unsh)%coeff_tij,map%theta_magpair_tij(map%magpairshell(unsh)%ind_global_tij),dv )
! !                     call lo_gemv( map%pairop(unop)%msotr,dv,ltij )
! !                 else
! !                     ltij=0.0_r8
! !                 endif
! !                 if ( map%magpairshell(unsh)%ntheta_qij .gt. 0 ) then
! !                     lqij=matmul(map%magpairshell(unsh)%coeff_qij,map%theta_magpair_qij(map%magpairshell(unsh)%ind_global_qij))
! !                     lqij=matmul(map%pairop(unop)%m3,lqij)
! !                     !call lo_gemv( map%magpairshell(unsh)%coeff_qij,map%theta_magpair_qij(map%magpairshell(unsh)%ind_global_qij),dv )
! !                     !call lo_gemv( map%pairop(unop)%sotr,dv,lqij )
! !                 else
! !                     lqij=0.0_r8
! !                 endif
! !
! !                 l=0
! !                 do ii=1,3
! !                 do jj=1,3
! !                     l=l+1
! !                     jij%atom(a1)%pair(i)%J(ii,jj)=ljij(l)
! !                     jij%atom(a1)%pair(i)%bird(ii,jj)=ltij(l)
! !                     jij%atom(a1)%pair(i)%Q(ii)=lqij(ii)
! !                     !jij%atom(a1)%pair(i)%Q(ii,jj)=lqij(l)
! !                 enddo
! !                 enddo
! !                 jij%atom(a1)%pair(i)%dude=-lo_huge
! !                 !
! !                 jij%atom(a1)%pair(i)%irreducible_shell=unsh
! !                 jij%atom(a1)%pair(i)%irreducible_operation=unop
! !                 ! A bunch of extra stuff
! !                 a2=map%uc(a1)%magpair(i)%i2
! !                 jij%atom(a1)%pair(i)%i1=a1
! !                 jij%atom(a1)%pair(i)%i2=a2
! !
! !                 r=matmul(uc%latticevectors,map%uc(a1)%magpair(i)%flv)
! !                 jij%atom(a1)%pair(i)%lv1=0.0_r8
! !                 jij%atom(a1)%pair(i)%lv2=r
! !                 jij%atom(a1)%pair(i)%r=r+uc%rcart(:,a2)-uc%rcart(:,a1)
! !                 ! and the cutoff
! !                 jij%cutoff=max(jij%cutoff,norm2(jij%atom(a1)%pair(i)%r)+10*lo_tol)
! !             enddo
! !         enddo
! !
! !         ! Now find the dude-terms from the transpose
! !         do a1=1,jij%na
! !         do i=1,jij%atom(a1)%n
! !             a2=jij%atom(a1)%pair(i)%i2
! !             do j=1,jij%atom(a2)%n
! !                 if ( jij%atom(a2)%pair(j)%i2 .ne. a1 ) cycle
! !                 if ( lo_sqnorm(jij%atom(a1)%pair(i)%r+jij%atom(a2)%pair(j)%r) .gt. lo_sqtol ) cycle
! !                 jij%atom(a1)%pair(i)%dude=transpose(jij%atom(a2)%pair(j)%bird)
! !             enddo
! !         enddo
! !         enddo
! !
! !         ! And set the bird & dude self-terms
! !         do a1=1,jij%na
! !             jij%atom(a1)%selfterm_bird=0.0_r8
! !             jij%atom(a1)%selfterm_dude=0.0_r8
! !             do i=1,map%uc(a1)%nmagpair
! !                 jij%atom(a1)%selfterm_bird=jij%atom(a1)%selfterm_bird-jij%atom(a1)%pair(i)%bird
! !                 jij%atom(a1)%selfterm_dude=jij%atom(a1)%selfterm_dude-jij%atom(a1)%pair(i)%dude
! !             enddo
! !             jij%atom(a1)%selfterm_bird=lo_chop( jij%atom(a1)%selfterm_bird, 1E-15_r8 )
! !             jij%atom(a1)%selfterm_dude=lo_chop( jij%atom(a1)%selfterm_dude, 1E-15_r8 )
! !         enddo
! !
! !         ! set tiny things to zero
! !         do a1=1,jij%na
! !         do i=1,jij%atom(a1)%n
! !             jij%atom(a1)%pair(i)%J=lo_chop(jij%atom(a1)%pair(i)%J,lo_sqtol)
! !             jij%atom(a1)%pair(i)%bird=lo_chop(jij%atom(a1)%pair(i)%bird,lo_sqtol)
! !             jij%atom(a1)%pair(i)%dude=lo_chop(jij%atom(a1)%pair(i)%dude,lo_sqtol)
! !             jij%atom(a1)%pair(i)%Q=lo_chop(jij%atom(a1)%pair(i)%Q,lo_sqtol)
! !         enddo
! !         enddo
! !     end block storebasic
! end subroutine

end submodule
