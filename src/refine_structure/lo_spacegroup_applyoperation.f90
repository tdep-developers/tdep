submodule(lo_spacegroup) lo_spacegroup_applyoperation
!! Generates symmetry operations
implicit none
contains

!> Apply an operation to a vector
module pure function lo_operate_on_vector(op, v, reciprocal, inverse, fractional) result(w)
    !> operation
    type(lo_spacegroup_operation), intent(in) :: op
    !> vector
    real(r8), dimension(3), intent(in) :: v
    !> should the inverse operation be applied?
    logical, intent(in), optional :: inverse
    !> is it an operation in reciprocal space
    logical, intent(in), optional :: reciprocal
    !> should it be done in fractional coordinates?
    logical, intent(in), optional :: fractional
    !> rotated vector
    real(r8), dimension(3) :: w

    logical :: rs, inv, fr
    real(r8), dimension(3, 3) :: m
    real(r8), dimension(3) :: t

    ! Ok, eight possibilites, forward or inverse, realspace or reciprocal, cartesian or fractional
    if (present(inverse)) then
        inv = inverse
    else
        inv = .false.
    end if
    if (present(reciprocal)) then
        rs = reciprocal
    else
        rs = .false.
    end if
    if (present(fractional)) then
        fr = fractional
    else
        fr = .false.
    end if

    ! Get the rotation part.
    if (fr) then
        if (inv) then
            if (rs) then
                m = op%irfm ! inverse, reciprocal,fractional
            else
                m = op%ifm  ! inverse, realspace, fractional
            end if
        else
            if (rs) then
                m = op%rfm  ! forward, reciprocal, fractional
            else
                m = op%fm   ! forward, realspace, fractional
            end if
        end if
    else
        if (inv) then
            m = op%im       ! Cartesian, real/reciprocal, inverse
        else
            m = op%m        ! Cartesian, real/reciprocal, forward
        end if
    end if
    ! get the translation part
    if (rs) then
        t = 0.0_r8
    else
        if (fr) then
            t = op%ftr
        else
            t = op%tr
        end if
    end if
    ! Begin operation "actual operation"
    if (inv) then
        w = v - t
        w = matmul(m, w)
    else
        w = matmul(m, v)
        w = w + t
    end if
end function

!> Apply operation to a 3x3 tensor
module pure function lo_operate_on_secondorder_tensor(op, m) result(n)
    !> The operation
    type(lo_spacegroup_operation), intent(in) :: op
    !> The original tensor
    real(r8), dimension(3, 3), intent(in) :: m
    !> The rotated tensor
    real(r8), dimension(3, 3) :: n

    integer :: i, j, ii, jj
    n = 0.0_r8
    do j = 1, 3
    do i = 1, 3
        do jj = 1, 3
        do ii = 1, 3
            n(i, j) = n(i, j) + m(ii, jj)*op%m(i, ii)*op%m(j, jj)
        end do
        end do
    end do
    end do
end function

!> Transformation matrix for phonon eigenvectors corresponding to a symmetry operation
module subroutine lo_eigenvector_transformation_matrix(gm, rcart, qv, op, inverseoperation)
    !> transformation matrix
    complex(r8), dimension(:, :), intent(out) :: gm
    !> position vectors in Cartesian coordinates, generally crystalstructure%rcart
    real(r8), dimension(:, :), intent(in) :: rcart
    !> the q-vector
    real(r8), dimension(3), intent(in) :: qv
    !> the point operation
    type(lo_spacegroup_operation), intent(in) :: op
    !> should the inverse transformation be used?
    logical, intent(in), optional :: inverseoperation

    complex(r8) :: expiqr
    real(r8), dimension(3) :: v0
    real(r8) :: iqvv
    integer :: a1, a2, na
    logical :: inv

    ! Forward or backwards operation
    if (present(inverseoperation)) then
        inv = inverseoperation
    else
        inv = .false.
    end if

    na = size(rcart, 2)
    !     ! Small sanity check
    ! #ifdef AGRESSIVE_SANITY
    !     if ( na*3 .ne. size(gm,1) ) then
    !         call lo_stop_gracefully(&
    !             ['inconsistent number of atoms, rotation matrix has size '//tochar(size(gm,1))//'x'//tochar(size(gm,2))//' and Na='//tochar(na)],&
    !             lo_exitcode_baddim,__FILE__,__LINE__)
    !     endif
    ! #endif

    gm = 0.0_r8
    do a2 = 1, na
        a1 = op%fmap(a2)
        v0 = lo_operate_on_vector(op, rcart(:, a1), inverse=.true.) - rcart(:, a2)
        iqvv = -dot_product(v0, qv)*lo_twopi
        expiqr = cmplx(cos(iqvv), sin(iqvv), r8)
        if (inv) then
            gm((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) = op%im*expiqr
        else
            gm((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) = op%m*expiqr
        end if
    end do
end subroutine

! !> expand a 3x3 symmetry operation to a 9x9 operation, to operate on flattened tensors
! pure function lo_expandoperation_pair(o) result(bigo)
!     !> the original 3x3 operation
!     real(r8), dimension(3,3), intent(in) :: o
!     !> the resulting 9x9 operation
!     real(r8), dimension(9,9) :: bigo
!     !
!     !integer :: i,j,ii,jj,k,l
!     !
!     !bigo=0.0_r8
!     !k=0
!     !do i=1,3
!     !do ii=1,3
!     !    k=k+1
!     !    l=0
!     !    do j=1,3
!     !    do jj=1,3
!     !        l=l+1
!     !        bigo(k,l)=o(ii,jj)*o(i,j)
!     !    enddo
!     !    enddo
!     !enddo
!     !enddo
!     bigo(1,1)=o(1,1)*o(1,1)
!     bigo(1,2)=o(1,2)*o(1,1)
!     bigo(1,3)=o(1,3)*o(1,1)
!     bigo(1,4)=o(1,1)*o(1,2)
!     bigo(1,5)=o(1,2)*o(1,2)
!     bigo(1,6)=o(1,3)*o(1,2)
!     bigo(1,7)=o(1,1)*o(1,3)
!     bigo(1,8)=o(1,2)*o(1,3)
!     bigo(1,9)=o(1,3)*o(1,3)
!     bigo(2,1)=o(2,1)*o(1,1)
!     bigo(2,2)=o(2,2)*o(1,1)
!     bigo(2,3)=o(2,3)*o(1,1)
!     bigo(2,4)=o(2,1)*o(1,2)
!     bigo(2,5)=o(2,2)*o(1,2)
!     bigo(2,6)=o(2,3)*o(1,2)
!     bigo(2,7)=o(2,1)*o(1,3)
!     bigo(2,8)=o(2,2)*o(1,3)
!     bigo(2,9)=o(2,3)*o(1,3)
!     bigo(3,1)=o(3,1)*o(1,1)
!     bigo(3,2)=o(3,2)*o(1,1)
!     bigo(3,3)=o(3,3)*o(1,1)
!     bigo(3,4)=o(3,1)*o(1,2)
!     bigo(3,5)=o(3,2)*o(1,2)
!     bigo(3,6)=o(3,3)*o(1,2)
!     bigo(3,7)=o(3,1)*o(1,3)
!     bigo(3,8)=o(3,2)*o(1,3)
!     bigo(3,9)=o(3,3)*o(1,3)
!     bigo(4,1)=o(1,1)*o(2,1)
!     bigo(4,2)=o(1,2)*o(2,1)
!     bigo(4,3)=o(1,3)*o(2,1)
!     bigo(4,4)=o(1,1)*o(2,2)
!     bigo(4,5)=o(1,2)*o(2,2)
!     bigo(4,6)=o(1,3)*o(2,2)
!     bigo(4,7)=o(1,1)*o(2,3)
!     bigo(4,8)=o(1,2)*o(2,3)
!     bigo(4,9)=o(1,3)*o(2,3)
!     bigo(5,1)=o(2,1)*o(2,1)
!     bigo(5,2)=o(2,2)*o(2,1)
!     bigo(5,3)=o(2,3)*o(2,1)
!     bigo(5,4)=o(2,1)*o(2,2)
!     bigo(5,5)=o(2,2)*o(2,2)
!     bigo(5,6)=o(2,3)*o(2,2)
!     bigo(5,7)=o(2,1)*o(2,3)
!     bigo(5,8)=o(2,2)*o(2,3)
!     bigo(5,9)=o(2,3)*o(2,3)
!     bigo(6,1)=o(3,1)*o(2,1)
!     bigo(6,2)=o(3,2)*o(2,1)
!     bigo(6,3)=o(3,3)*o(2,1)
!     bigo(6,4)=o(3,1)*o(2,2)
!     bigo(6,5)=o(3,2)*o(2,2)
!     bigo(6,6)=o(3,3)*o(2,2)
!     bigo(6,7)=o(3,1)*o(2,3)
!     bigo(6,8)=o(3,2)*o(2,3)
!     bigo(6,9)=o(3,3)*o(2,3)
!     bigo(7,1)=o(1,1)*o(3,1)
!     bigo(7,2)=o(1,2)*o(3,1)
!     bigo(7,3)=o(1,3)*o(3,1)
!     bigo(7,4)=o(1,1)*o(3,2)
!     bigo(7,5)=o(1,2)*o(3,2)
!     bigo(7,6)=o(1,3)*o(3,2)
!     bigo(7,7)=o(1,1)*o(3,3)
!     bigo(7,8)=o(1,2)*o(3,3)
!     bigo(7,9)=o(1,3)*o(3,3)
!     bigo(8,1)=o(2,1)*o(3,1)
!     bigo(8,2)=o(2,2)*o(3,1)
!     bigo(8,3)=o(2,3)*o(3,1)
!     bigo(8,4)=o(2,1)*o(3,2)
!     bigo(8,5)=o(2,2)*o(3,2)
!     bigo(8,6)=o(2,3)*o(3,2)
!     bigo(8,7)=o(2,1)*o(3,3)
!     bigo(8,8)=o(2,2)*o(3,3)
!     bigo(8,9)=o(2,3)*o(3,3)
!     bigo(9,1)=o(3,1)*o(3,1)
!     bigo(9,2)=o(3,2)*o(3,1)
!     bigo(9,3)=o(3,3)*o(3,1)
!     bigo(9,4)=o(3,1)*o(3,2)
!     bigo(9,5)=o(3,2)*o(3,2)
!     bigo(9,6)=o(3,3)*o(3,2)
!     bigo(9,7)=o(3,1)*o(3,3)
!     bigo(9,8)=o(3,2)*o(3,3)
!     bigo(9,9)=o(3,3)*o(3,3)
!     bigo=lo_chop(bigo,1E-11_r8)
! end function
!
! !> expand a 3x3 symmetry operation to a 27x27 operation, to operate on vector format of tensors
! pure function lo_expandoperation_triplet(o) result(bigo)
!     !> the original 3x3 operation
!     real(r8), dimension(3,3), intent(in) :: o
!     !> the resulting 27x27 operation
!     real(r8), dimension(27,27) :: bigo
!
!     integer :: i,ii,iii,j,jj,jjj,l,m
!
!     bigo=0.0_r8
!     l=0
!     do i=1,3
!     do ii=1,3
!     do iii=1,3
!         l=l+1
!         m=0
!         do j=1,3
!         do jj=1,3
!         do jjj=1,3
!             m=m+1
!             bigo(l,m)=o(iii,jjj)*o(ii,jj)*o(i,j)
!         enddo
!         enddo
!         enddo
!     enddo
!     enddo
!     enddo
!     bigo=lo_chop(bigo,1E-11_r8)
! end function
!
! !> expand a 3x3 symmetry operation to a 81x81 operation, to operate on vector format of tensors
! pure function lo_expandoperation_quartet(o) result(bigo)
!     !> the original 3x3 operation
!     real(r8), dimension(3,3), intent(in) :: o
!     !> the resulting 27x27 operation
!     real(r8), dimension(81,81) :: bigo
!
!     integer :: i,ii,iii,iiii,j,jj,jjj,jjjj,l,m
!     bigo=0.0_r8
!     l=0
!     do i=1,3
!     do ii=1,3
!     do iii=1,3
!     do iiii=1,3
!         l=l+1
!         m=0
!         do j=1,3
!         do jj=1,3
!         do jjj=1,3
!         do jjjj=1,3
!             m=m+1
!             bigo(l,m)=o(iiii,jjjj)*o(iii,jjj)*o(ii,jj)*o(i,j)
!         enddo
!         enddo
!         enddo
!         enddo
!     enddo
!     enddo
!     enddo
!     enddo
!     bigo=lo_chop(bigo,1E-11_r8)
! end function

end submodule
