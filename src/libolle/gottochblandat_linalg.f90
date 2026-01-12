submodule (gottochblandat) gottochblandat_linalg
implicit none
contains

!> Square root that produces negative results instead of imaginary for negative numbers
module elemental function lo_negsqrt(x) result(y)
    real(r8), intent(in) :: x
    real(r8) :: y
    y=sqrt(abs(x))*int(anint(sign(1.0_r8,x)))
end function

!> Fortran seems to not have a simple sign function so I just wrapped it neatly since I never remember the syntax.
module elemental function lo_sgn(x) result(sgn)
    real(r8), intent(in) :: x
    integer :: sgn
    sgn=int(anint(sign(1.0_r8,x)))
end function

!> unsigned volume of a tetrahedron
module pure function lo_unsigned_tetrahedron_volume(nodes) result(volume)
    !> the nodes of the tetrahedron
    real(r8), dimension(3,4), intent(in) :: nodes
    !> the volume
    real(r8) :: volume

    real(r8), dimension(3) :: a,b,c,d,v
    a=nodes(:,1)
    b=nodes(:,2)
    c=nodes(:,3)
    d=nodes(:,4)
    a=a-d
    b=b-d
    c=c-d

    v=lo_cross(b,c)
    volume=abs(dot_product(a,v)/6.0_r8)
end function

!> signed volume of a tetrahedron
module pure function lo_signed_tetrahedron_volume(nodes) result(volume)
    !> the nodes of the tetrahedron
    real(r8), dimension(3,4), intent(in) :: nodes
    !> the signed volume
    real(r8) :: volume

    real(r8), dimension(3) :: a,b,c,d,v
    a=nodes(:,1)
    b=nodes(:,2)
    c=nodes(:,3)
    d=nodes(:,4)
    a=a-d
    b=b-d
    c=c-d

    v=lo_cross(b,c)
    volume=dot_product(a,v)/6.0_r8
end function

!> The Frobenius norm of a matrix @todo retire in favor of intrinsic
module pure function lo_frobnorm(m) result(nrm)
    !> matrix
    real(r8), dimension(:,:), intent(in) :: m
    !> norm
    real(r8) :: nrm

    integer :: i,j
    nrm=0.0_r8
    do j=1,size(m,2)
    do i=1,size(m,1)
        nrm=nrm+m(i,j)*m(i,j)
    enddo
    enddo
    nrm=sqrt(nrm)
end function

!> Trace of a matrix
module pure function lo_trace(m) result(tr)
    !> matrix
    real(r8), dimension(:,:), intent(in) :: m
    !> trace
    real(r8) :: tr

    integer :: i
    tr=0.0_r8
    do i=1,size(m,1)
        tr=tr+m(i,i)
    enddo
end function

!> Squared of norm
pure module function lo_sqnorm(a) result(nrm)
    !> vector
    real(r8), dimension(3), intent(in) :: a
    !> squared norm
    real(r8) :: nrm

    nrm=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
end function

!> Cross product
module pure function lo_cross(b,c) result(a)
    !> first vector
    real(r8), dimension(3), intent(in) :: b
    !> second vector
    real(r8), dimension(3), intent(in) :: c
    !> cross product
    real(r8), dimension(3) :: a

    a(1)=b(2)*c(3)-b(3)*c(2)
    a(2)=b(3)*c(1)-b(1)*c(3)
    a(3)=b(1)*c(2)-b(2)*c(1)
end function

!> vector outer (dyadic) product, three-dimensional
module pure function lo_real_outerproduct(a,b) result(m)
    !> first vector
    real(r8), dimension(3), intent(in) :: a
    !> second vector
    real(r8), dimension(3), intent(in) :: b
    !> cross product
    real(r8), dimension(3,3) :: m

    m(1,1)=a(1)*b(1)
    m(2,1)=a(1)*b(2)
    m(3,1)=a(1)*b(3)
    m(1,2)=a(2)*b(1)
    m(2,2)=a(2)*b(2)
    m(3,2)=a(2)*b(3)
    m(1,3)=a(3)*b(1)
    m(2,3)=a(3)*b(2)
    m(3,3)=a(3)*b(3)
end function

!> vector outer product, three-dimensional
module pure function lo_complex_outerproduct(a,b) result(m)
    !> first vector
    complex(r8), dimension(3), intent(in) :: a
    !> second vector
    complex(r8), dimension(3), intent(in) :: b
    !> cross product
    complex(r8), dimension(3,3) :: m

    m(1,1)=conjg(a(1))*b(1)
    m(2,1)=conjg(a(1))*b(2)
    m(3,1)=conjg(a(1))*b(3)
    m(1,2)=conjg(a(2))*b(1)
    m(2,2)=conjg(a(2))*b(2)
    m(3,2)=conjg(a(2))*b(3)
    m(1,3)=conjg(a(3))*b(1)
    m(2,3)=conjg(a(3))*b(2)
    m(3,3)=conjg(a(3))*b(3)
end function

!> Determinant of a 3x3 matrix, doubles
module pure function lo_determ_real(a) result(det)
    !> 3x3 matrix
    real(r8), dimension(3,3), intent(in) :: a
    !> the determinant
    real(r8) :: det
    !
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
        + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
end function

!> Determinant of a 3x3 matrix, integers
module pure function lo_determ_int(a) result(det)
    !> 3x3 matrix
    integer, dimension(3,3), intent(in) :: a
    !> the determinant
    integer :: det
    !
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
        + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
end function

!> Invert 3x3 matrix, without lapack
pure module function lo_invert3x3matrix(m) result(n)
    !> 3x3 matrix to invert
    real(r8), dimension(3,3), intent(in) :: m
    !> the inverse
    real(r8), dimension(3,3) :: n

    real(r8) :: det

    det =  m(1,1)*m(2,2)*m(3,3) - m(1,1)*m(2,3)*m(3,2)&
         - m(1,2)*m(2,1)*m(3,3) + m(1,2)*m(2,3)*m(3,1)&
         + m(1,3)*m(2,1)*m(3,2) - m(1,3)*m(2,2)*m(3,1)
    det=1.0_r8/det
    ! Calculate the inverse of the matrix
    n(1,1)=+det*(m(2,2)*m(3,3)-m(2,3)*m(3,2))
    n(2,1)=-det*(m(2,1)*m(3,3)-m(2,3)*m(3,1))
    n(3,1)=+det*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
    n(1,2)=-det*(m(1,2)*m(3,3)-m(1,3)*m(3,2))
    n(2,2)=+det*(m(1,1)*m(3,3)-m(1,3)*m(3,1))
    n(3,2)=-det*(m(1,1)*m(3,2)-m(1,2)*m(3,1))
    n(1,3)=+det*(m(1,2)*m(2,3)-m(1,3)*m(2,2))
    n(2,3)=-det*(m(1,1)*m(2,3)-m(1,3)*m(2,1))
    n(3,3)=+det*(m(1,1)*m(2,2)-m(1,2)*m(2,1))
end function

!> Sqrt of a 3x3 matrix
module function lo_sqrt3x3matrix(matrix) result(matrix_sqrt)
    !> 3x3 matrix to sqrt
    real(r8), dimension(3,3), intent(in) :: matrix
    !> the sqrt
    real(r8), dimension(3,3) :: tmp, matrix_sqrt

    !> eigenvalues
    real(r8), dimension(3) :: eigenvalues
    !> eigenvectors
    real(r8), dimension(3,3) :: eigenvectors
    !> other things
    integer :: ii, jj, kk

    ! init
    tmp = 0.0_r8
    matrix_sqrt = 0.0_r8

    ! we want to compute sqrt(m)
    ! say m = \sum_i m_i |i><i|
    ! where |i> are the eigenvectors and m_i the eigenvalues
    ! then sqrt(m) = \sum_i sqrt(m_i) |i><i|

    ! let's compute m_i and |i>

    call lo_symmetric_eigensystem_3x3matrix(matrix, eigenvalues, eigenvectors)

    ! do the sqrt and sanity check if we got eigenvalue decomposition right
    do ii = 1, 3
        do jj = 1, 3
            do kk = 1, 3
                associate (&
                    mij => tmp(ii, jj), &
                    sij => matrix_sqrt(ii, jj), &
                    ek => eigenvalues(kk), &
                    vki => eigenvectors(kk, ii), &
                    vkj => eigenvectors(kk, jj) &
                    )
                    mij = mij + ek * vki * vkj
                    sij = sij + sqrt(ek) * vki * vkj
                end associate
            end do
            if (abs(tmp(ii, jj) - matrix(ii, jj)) > lo_tol) then
                write(*,*) tmp(ii, jj)
                write(*,*) matrix(ii, jj)
                call lo_stop_gracefully(['Got the eigenbasis of matrix wrong. HELP.'], lo_exitcode_blaslapack)
            end if
        end do
    end do
    ! made it here, lucky day
end function

!> Gram-Schmidt orthoginalization of a real matrix
module pure subroutine lo_real_gram_schmidt(X)
    !> matrix to orthoginalize. Result should be X^TX=I
    real(r8), dimension(:,:), intent(inout) :: X

    real(r8), dimension(:,:), allocatable :: Q,R
    integer :: i,k,nr,nc

    nr=size(X,1)
    nc=size(X,2)
    allocate(Q(nr,nc))
    allocate(R(nc,nc))
    Q=0.0_r8
    R=0.0_r8
    do k = 1,nc
        Q(:,k) = X(:,k)
        do i = 1,k-1
            R(i,k)=dot_product(Q(:,i),Q(:,k))
            Q(:,k)=Q(:,k)-R(i,k)*Q(:,i)
        enddo
        R(k,k) = norm2(Q(:,k))
        if ( abs(R(k,k)) .gt. lo_sqtol ) then
            Q(:,k)=Q(:,k)/R(k,k)
        else
            Q(:,k)=0.0_r8
        endif
    enddo
    X=lo_chop(Q,1E-13_r8)
    deallocate(Q,R)
end subroutine

!> Gram-Shmidt orthoginalization of a complex matrix
module pure subroutine lo_complex_gram_schmidt(X)
    !> matrix to orthoginalize. Result should be X^TX=I
    complex(r8), dimension(:,:), intent(inout) :: X

    complex(r8), dimension(:,:), allocatable :: Q,R
    integer :: i,k,nr,nc

    nr=size(X,1)
    nc=size(X,2)
    allocate(Q(nr,nc))
    allocate(R(nc,nc))
    Q=0.0_r8
    R=0.0_r8
    do k = 1,nc
        Q(:,k) = X(:,k)
        do i = 1,k-1
            R(i,k)=dot_product(Q(:,i),Q(:,k))
            Q(:,k)=Q(:,k)-R(i,k)*Q(:,i)
        enddo
        R(k,k) = sqrt( dot_product(Q(:,k),Q(:,k)) )
        if ( abs(R(k,k)) .gt. lo_sqtol ) then
            Q(:,k)=Q(:,k)/R(k,k)
        else
            Q(:,k)=0.0_r8
        endif
    enddo
    X=lo_chop(Q,1E-13_r8)
    deallocate(Q,R)
end subroutine

!> create identity matrices
module pure subroutine lo_identitymatrix(m)
    real(r8), dimension(:,:), intent(out) :: m

    integer :: i
    m=0.0_r8
    do i=1,size(m,1)
        m(i,i)=1.0_r8
    enddo
end subroutine

!> wrapper around gesvd. Will add things to return the singular vectors and that stuff when I need it.
module subroutine lo_real_singular_value_decomposition(A,S,U,V)
    !> matrix A
    real(r8), dimension(:,:), intent(in) :: A
    !> singular values
    real(r8), dimension(:), allocatable, intent(out) :: S
    !> left singular vectors
    real(r8), dimension(:,:), allocatable, intent(out), optional :: U
    !> right singular vectors
    real(r8), dimension(:,:), allocatable, intent(out), optional :: V

    real(r8), dimension(:,:), allocatable :: Sigma
    real(r8), dimension(:), allocatable :: wrk
    real(r8), dimension(1,1) :: dV,dU
    real(r8), dimension(1) :: dwrk
    integer :: nc,nr,ns,lwork

    nc=size(A,1)
    nr=size(A,2)
    ns=min(nc,nr)
    allocate(Sigma(nc,nr))
    allocate(S(ns))
    Sigma=A
    S=0.0_r8
    lwork=-1

    if ( present(U) .and. present(V) ) then
        ! Return left and right vectors
        allocate(U(nc,nc))
        allocate(V(nr,nr))
        U=0.0_r8
        V=0.0_r8
        call dgesvd('A','A',nc,nr,Sigma,nc,S,U,nc,V,nr,dwrk,lwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['dgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        lwork=int(anint(dwrk(1)))
        allocate(wrk(lwork))
        call dgesvd('A','A',nc,nr,Sigma,nc,S,U,nc,V,nr,wrk,lwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['dgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    else
        ! No vectors returned
        call dgesvd('N','N',nc,nr,Sigma,nc,S,dU,nc,dV,1,dwrk,lwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['dgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        lwork=int(anint(dwrk(1)))
        allocate(wrk(lwork))
        call dgesvd('N','N',nc,nr,Sigma,nc,S,dU,nc,dV,1,wrk,lwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['dgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    endif
    deallocate(wrk)
    deallocate(Sigma)
end subroutine

!> wrapper around gesvd. Will add things to return the singular vectors and that stuff when I need it.
module subroutine lo_complex_singular_value_decomposition(A,S,U,V)
    !> matrix A
    complex(r8), dimension(:,:), intent(in) :: A
    !> singular values
    real(r8), dimension(:), allocatable, intent(out) :: S
    !> left singular vectors
    complex(r8), dimension(:,:), allocatable, intent(out), optional :: U
    !> right singular vectors
    complex(r8), dimension(:,:), allocatable, intent(out), optional :: V

    complex(r8), dimension(:,:), allocatable :: Sigma
    real(r8), dimension(:), allocatable :: wrk,rwork
    complex(r8), dimension(1,1) :: dV,dU
    real(r8), dimension(1) :: dwrk
    integer :: nc,nr,ns,lwork

    nc=size(A,1)
    nr=size(A,2)
    ns=min(nc,nr)
    allocate(Sigma(nc,nr))
    allocate(S(ns))
    allocate(rwork(5*max(nr,nc)))
    Sigma=A
    S=0.0_r8
    rwork=0.0_r8
    lwork=-1
    if ( present(U) .and. present(V) ) then
        ! Return left and right vectors
        allocate(U(nc,nc))
        allocate(V(nr,nr))
        U=0.0_r8
        V=0.0_r8
        call zgesvd('A','A',nc,nr,Sigma,nc,S,U,nc,V,nr,dwrk,lwork,rwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['zgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        lwork=int(anint(dwrk(1)))
        allocate(wrk(lwork))
        call zgesvd('A','A',nc,nr,Sigma,nc,S,U,nc,V,nr,wrk,lwork,rwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['zgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    else
        ! Don't return vectors
        call zgesvd('N','N',nc,nr,Sigma,nc,S,dU,nc,dV,1,dwrk,lwork,rwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['zgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        lwork=int(anint(dwrk(1)))
        allocate(wrk(lwork))
        call zgesvd('N','N',nc,nr,Sigma,nc,S,dU,nc,dV,1,wrk,lwork,rwork,lo_status)
        if ( lo_status .ne. 0 ) call lo_stop_gracefully(['zgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    endif
    deallocate(wrk)
    deallocate(rwork)
    deallocate(Sigma)
end subroutine

!> Use lagrange multipliers to enforce linear constraints: will modify X such that A*X=0 by adding a correction proportional to X^2
module subroutine lo_enforce_linear_constraints(A,X,tolerance,nosquare)
    !> constraint matrix A
    real(r8), dimension(:,:), intent(in) :: A
    !> vector X
    real(r8), dimension(:), intent(inout) :: X
    !> tolerance
    real(r8), intent(in), optional :: tolerance
    !> don't make the correction proportional
    logical, intent(in), optional :: nosquare

    real(r8), parameter :: largefactor=1E20_r8
    integer, parameter :: blascrossover=200

    real(r8), dimension(:,:), allocatable :: bigM,bigV
    real(r8), dimension(:), allocatable :: work
    real(r8) :: f0,s_work(1),tol,largenum
    integer :: nx,neq,i,j
    integer :: m,n,nrhs,lda,ldb,lwork
    logical :: addsquare

    ! Set the tolerance
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_sqtol
    endif

    ! Add the inverse square thing
    if ( present(nosquare) ) then
        addsquare=.not.nosquare
    else
        addsquare=.true.
    endif

    nx=size(X,1)  ! number of x-values
    neq=size(A,1) ! number of constraints
    ! Make sure it makes sense to try to enforce the constraints
    if ( neq .gt. nx ) then
        call lo_stop_gracefully(['You have more constraints than variables you are trying to constrain, that makes no sense.'],lo_exitcode_param)
    endif
    if ( size(A,2) .ne. nx ) then
        call lo_stop_gracefully(['Incompatible dimensions for constraints, makes no sense.'],lo_exitcode_param)
    endif

    ! Maybe I don't even need to enforce them?
    if ( nx .gt. blascrossover ) then
        allocate(work(nx))
        work=0.0_r8
        call dgemv('N',neq,nx,1.0_r8,A,neq,X,1,0.0_r8,work,1)
        f0=sum(abs(work))
        deallocate(work)
    else
        f0=sum(abs(matmul(A,X)))
    endif
    if ( f0 .lt. tol ) then
        ! The constraints are already satisfied.
        return
    endif

    ! Ok if we made it here, I actually have to do something. First, find a large number.
    ! Don't want to use huge since that causes overflows and errors and stuff in odd cases.
    f0=lo_huge
    do i=1,nx
        if ( abs(X(i)) .gt. tol ) f0=min(f0,abs(X(i)))
    enddo
    largenum=largefactor/f0**2

    ! Now build the enlarged system of equations to determine the Lagrange multipliers.
    n=nx+neq
    allocate(bigM(n,n))
    allocate(bigV(n,1))
    bigM=0.0_r8
    do i=1,nx
    do j=1,neq
        bigM(nx+j,i)=A(j,i)
        bigM(i,nx+j)=A(j,i)
    enddo
    enddo
    do i=1,nx
        if ( addsquare ) then
            f0=X(i)
            if ( abs(f0) .gt. tol ) then
                bigM(i,i)=1.0_r8/f0**2
            else
                bigM(i,i)=largenum
            endif
        else
            bigM(i,i)=1.0_r8
        endif
    enddo
    bigV=0.0_r8
    bigV(nx+1:n,1)=matmul(A,X)

    lda=size(bigM,1)
    ldb=size(bigV,1)
    m=size(bigM,1)
    n=size(bigM,2)
    nrhs=1
    lwork = -1
    call dgels('N',m,n,nrhs,bigM,lda,bigV,ldb,s_work,lwork,lo_status)
    if ( lo_status .ne. 0 ) then
        call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_param,__FILE__,__LINE__)
    endif
    lwork = int(s_work(1))
    allocate(work(lwork))
    call dgels('N',m,n,nrhs,bigM,lda,bigV,ldb,work,lwork,lo_status)
    if ( lo_status .ne. 0 ) then
        call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! And here we actually enforce the constraints
    X=X-bigV(1:nx,1)
    ! And some cleanup.
    deallocate(work)
    deallocate(bigV)
    deallocate(bigM)
end subroutine

!> Make the coefficient matrix a bit tidier, essentially reduced row echelon form.
pure module subroutine lo_make_coeffmatrix_tidy(m,tolerance)
    !> the input coefficient matrix
    real(r8), dimension(:,:), intent(inout) :: m
    !> tolerance?
    real(r8), intent(in), optional :: tolerance

    !> the clean and nice matrix
    real(r8), dimension(size(m,1),size(m,2)) :: n

    integer :: i,j,k,col,nc,nr
    real(r8) :: f0,f1,tol

    ! Start with a Gram-Schmidt, always a reasonable thing to do
    call lo_real_gram_schmidt(m)

    ! dimensions of matrix
    nc=size(m,2)
    nr=size(m,1)
    ! tolerance?
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_sqtol
    endif

    ! Remove tiny numbers
    n=lo_chop(m,tol)

    ! Forward reduce
    do j=1,nc
        ! skip 0-vectors
        if ( sum(abs(n(:,j))) .lt. tol ) cycle
        ! find the first non-zero number in this column
        col=0
        f1=lo_huge
        do k=1,nr
            if ( abs(n(k,j)) .gt. tol ) then
                f1=n(k,j)
                col=k
                exit
            endif
        enddo
        ! go to the other columns. If they have a number in this column, reduce them
        do i=j+1,nc
            ! is there a number?
            if ( abs(n(col,i)) .gt. tol ) then
                f0=n(col,i)/f1
                n(:,i)=n(:,i)-n(:,j)*f0
            endif
        enddo
    enddo

    ! Backwards reduce
    do j=nc,1,-1
        ! skip 0-vectors
        if ( sum(abs(n(:,j))) .lt. tol ) cycle
        ! find the first non-zero number in this column
        col=0
        f1=lo_huge
        do k=1,nc
            if ( abs(n(k,j)) .gt. tol ) then
                f1=n(k,j)
                col=k
                exit
            endif
        enddo
        ! go to the other columns. If they have a number in this column, reduce them
        do i=j-1,1,-1
            ! is there a number?
            if ( abs(n(col,i)) .gt. tol ) then
                f0=n(col,i)/f1
                n(:,i)=n(:,i)-n(:,j)*f0
            endif
        enddo
    enddo
    ! Remove tiny numbers again, and normalize
    n=lo_chop(n,tol)
    do i=1,nc
        f0=minval(n(:,i))
        f1=maxval(n(:,i))
        if ( abs(f0)+abs(f1) .lt. tol ) cycle
        if ( abs(f0) .gt. abs(f1) ) then
            n(:,i)=n(:,i)/f0
        else
            n(:,i)=n(:,i)/f1
        endif
    enddo
    ! And slice out tiny things again, to make it neat.
    m=lo_chop(n,tol)
end subroutine

!> Make the coefficient matrix a bit tidier, essentially reduced row echelon form.
module pure subroutine lo_make_complex_coeffmatrix_tidy(m,tolerance)
    !> input coefficient matrix
    complex(r8), dimension(:,:), intent(inout) :: m
    !> tolerance?
    real(r8), intent(in), optional :: tolerance

    !> the clean and nice matrix
    complex(r8), dimension(size(m,1),size(m,2)) :: n

    integer :: i,j,k,col,nc,nr
    complex(r8) :: f0,f1
    real(r8) :: tol

    ! Start with a Gram-Schmidt, always a reasonable thing to do
    call lo_complex_gram_schmidt(m)

    ! dimensions of matrix
    nc=size(m,2)
    nr=size(m,1)
    ! tolerance?
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_sqtol
    endif

    ! Remove tiny numbers
    n=lo_chop(m,tol)

    ! Forward reduce
    do j=1,nc
        ! skip 0-vectors
        if ( sum(abs(n(:,j))) .lt. tol ) cycle
        ! find the first non-zero number in this column
        col=0
        f1=lo_huge
        do k=1,nr
            if ( abs(n(k,j)) .gt. tol ) then
                f1=n(k,j)
                col=k
                exit
            endif
        enddo
        if ( col .eq. 0 ) cycle
        ! go to the other columns. If they have a number in this column, reduce them
        do i=j+1,nc
            ! is there a number?
            if ( abs(n(col,i)) .gt. tol ) then
                f0=n(col,i)/f1
                n(:,i)=n(:,i)-n(:,j)*f0
            endif
        enddo
    enddo

    ! Backwards reduce
    do j=nc,1,-1
        ! skip 0-vectors
        if ( sum(abs(n(:,j))) .lt. tol ) cycle
        ! find the first non-zero number in this column
        col=0
        f1=lo_huge
        do k=1,nc
            if ( abs(n(k,j)) .gt. tol ) then
                f1=n(k,j)
                col=k
                exit
            endif
        enddo
        if ( col .eq. 0 ) cycle
        ! go to the other columns. If they have a number in this column, reduce them
        do i=j-1,1,-1
            ! is there a number?
            if ( abs(n(col,i)) .gt. tol ) then
                f0=n(col,i)/f1
                n(:,i)=n(:,i)-n(:,j)*f0
            endif
        enddo
    enddo

    ! Remove tiny numbers again, and normalize
    do i=1,nc
        f0=dot_product(n(:,i),n(:,i))
        if ( abs(f0) .gt. tol ) then
            n(:,i)=n(:,i)/sqrt(f0)
        endif
    enddo

    ! And slice out tiny things again, to make it neat.
    m=lo_chop(n,tol)
end subroutine

!> compress a set of equations, as in the "A" in A*x=B, to minimal linearly independent form
module subroutine lo_compress_equations(equations,n_compressed_equations,compressed_equations,trans,tolerance)
    !> original equations
    real(r8), dimension(:,:), intent(in) :: equations
    !> how many equations did it become? Note that it could be zero.
    integer, intent(out) :: n_compressed_equations
    !> compressed equations
    real(r8), dimension(:,:), allocatable, intent(out) :: compressed_equations
    !> are the equations transposed on entry?
    logical, intent(in) :: trans
    !> tolerance
    real(r8), intent(in), optional :: tolerance

    real(r8), dimension(:,:), allocatable :: Sigma,U
    real(r8), dimension(:), allocatable :: S,wrk
    real(r8), dimension(1,1) :: V
    real(r8), dimension(1) :: dwrk
    real(r8) :: tol
    integer :: i,l,nx,ne,ns,lwork

    ! set the tolerance
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_sqtol
    endif

    ! get the dimensions, and figure out the transpose thingy
    if ( trans ) then
        ne=size(equations,2)
        nx=size(equations,1)
        allocate(Sigma(nx,ne))
        Sigma=lo_chop(equations,tol)
    else
        ne=size(equations,1)
        nx=size(equations,2)
        allocate(Sigma(nx,ne))
        Sigma=lo_chop(transpose(equations),tol)
    endif

    ns=min(nx,ne)
    allocate(S(ns))
    allocate(U(nx,nx))
    S=0.0_r8
    U=0.0_r8
    V=0.0_r8
    ! I want a singular value decomposition, Sigma=U*S*V, first query workspace,
    ! then do the actual solution
    lwork=-1
    call dgesvd('A','N',nx,ne,Sigma,nx,S,U,nx,V,1,dwrk,lwork,lo_status)
    lwork=int(anint(dwrk(1)))
    allocate(wrk(lwork))
    call dgesvd('A','N',nx,ne,Sigma,nx,S,U,nx,V,1,wrk,lwork,lo_status)

    ! count number of compressed equations
    l=0
    do i=1,ns
        if ( abs(s(i)) .gt. tol ) l=l+1
    enddo
    n_compressed_equations=l
    if ( n_compressed_equations .gt. 0 ) then
        allocate(compressed_equations(n_compressed_equations,nx))
        compressed_equations=U(1:n_compressed_equations,1:nx)
        compressed_equations=lo_chop( compressed_equations, tol )
    endif

    deallocate(Sigma)
    deallocate(U)
    deallocate(S)
    deallocate(wrk)
end subroutine

!> Solve for the nullspace of a matrix, used to get irreducible representations. Not fast, but that should not matter.
module subroutine lo_nullspace_coefficient_matrix(invarM,coeffM,nvar,varind,tolerance)
    !> matrix that leaves the desired quantity invariant
    real(r8), dimension(:,:), intent(in) :: invarM
    !> coefficient matrix
    real(r8), dimension(:,:), allocatable, intent(out) :: coeffM
    !> number of independent variable
    integer, intent(out) :: nvar
    !> where are the irreducible in terms of the full
    integer, dimension(:), allocatable, intent(out), optional :: varind
    !> tolerance
    real(r8), intent(in), optional :: tolerance

    real(r8), dimension(:,:), allocatable :: Sigma,U,V
    real(r8), dimension(:), allocatable :: S,wrk
    real(r8), dimension(1) :: dwrk
    real(r8) :: tol
    integer, dimension(:), allocatable :: zeroSV
    integer :: n,lwork,nc,i,j

    ! set the tolerance
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_sqtol
    endif

    ! make some space
    n=size(invarM,1)
    allocate(Sigma(n,n))
    allocate(U(n,n))
    allocate(V(n,n))
    allocate(S(n))
    allocate(zeroSV(n))
    Sigma=lo_chop( invarM, tol )
    U=0.0_r8
    V=0.0_r8
    S=0.0_r8
    zeroSV=0
    ! I want a singular value decomposition, Sigma=U*S*V, first query workspace, then do the actual decomposition
    lwork=-1
    call dgesvd('A','A',n,n,Sigma,n,S,U,n,V,n,dwrk,lwork,lo_status)
    if ( lo_status .ne. 0 ) call lo_stop_gracefully(['dgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    lwork=int(anint(dwrk(1)))
    allocate(wrk(lwork))
    call dgesvd('A','A',n,n,Sigma,n,S,U,n,V,n,wrk,lwork,lo_status)
    if ( lo_status .ne. 0 ) call lo_stop_gracefully(['dgesvd exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    deallocate(wrk)

    ! Count the number of zero singular values
    nc=0
    do i=1,n
        if ( abs(S(i)) .lt. tol ) then
            nc=nc+1
            zeroSV(nc)=i
        endif
    enddo

    ! now return the irreducible represantation.
    if ( nc .eq. 0 ) then
        ! There could be zero independent things. Or? Yes it could.
        nvar=0
    elseif ( nc .eq. n ) then
        ! It could also be that the invarM did not help at all.
        allocate(coeffM(n,n))
        call lo_identitymatrix(coeffM)
        nvar=n
        if ( present(varind) ) then
            allocate(varind(n))
            do i=1,n
                varind(i)=i
            enddo
        endif
    else
        ! If I made it here, it means it's a non-trivial solution.
        ! I reuse the Sigma array, now I it's the full-rank nullspace
        Sigma=matmul(U(:,zeroSV(1:nc)),V(zeroSV(1:nc),:))
        ! get the coefficient matrix to reduced row echelon form
        call lo_make_coeffmatrix_tidy(Sigma,tol)

        ! Now return the non-zero parts of this
        nvar=nc
        allocate(coeffM(n,nvar))
        coeffM=0.0_r8
        if ( present(varind) ) then
            allocate(varind(nvar))
            varind=0
        endif
        j=0
        do i=1,n
            if ( sum(abs(Sigma(:,i))) .gt. tol ) then
                j=j+1
                if ( present(varind) ) varind(j)=i
                coeffM(:,j)=Sigma(:,i)
            endif
        enddo
    endif

    ! And some cleanup
    deallocate(Sigma)
    deallocate(U)
    deallocate(V)
    deallocate(S)
    deallocate(zeroSV)
end subroutine

!> Solve for the nullspace of a matrix, used to get irreducible representations. Not fast, but that should not matter.
module subroutine lo_complex_nullspace_coefficient_matrix(invarM,invariant_operations,coeff,nvar,varind,tolerance)
    !> matrix that leaves the desired quantity invariant
    complex(r8), dimension(:,:), intent(in), optional :: invarM
    !> list of operations
    complex(r8), dimension(:,:,:), intent(in), optional :: invariant_operations
    !> coefficient matrix
    complex(r8), dimension(:,:), allocatable, intent(out) :: coeff
    !> number of independent variable
    integer, intent(out) :: nvar
    !> where are the irreducible in terms of the full
    integer, dimension(:), allocatable, intent(out), optional :: varind
    !> tolerance
    real(r8), intent(in), optional :: tolerance

    complex(r8), dimension(:,:), allocatable :: IM
    real(r8) :: tol
    integer :: n

    ! Sort out the heuristics
    init: block
        complex(r8), dimension(:,:,:), allocatable :: ops
        integer :: i,j,ctr

        ! set the tolerance
        if ( present(tolerance) ) then
            tol=tolerance
        else
            tol=1E-9_r8
        endif

        ! check input
        if ( present(invariant_operations) .and. present(invarM) ) then
            call lo_stop_gracefully(['Choose either invarM or invariant_operations as input, not both'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! Get the invariant matrix, either as input or from the operations
        if ( present(invarM) ) then
            ! First check that we don't have conflicting input
            n=size(invarM,1)
            allocate(IM(n,n))
            IM=lo_chop(invarM,tol)
        else
            ! Be extra careful and don't add redundant operations.
            n=size(invariant_operations,1)
            i=size(invariant_operations,3)
            allocate(IM(n,n))
            allocate(ops(n,n,i))
            ops=-1.0_r8
            ctr=0
            opl: do i=1,size(invariant_operations,3)
                do j=1,ctr
                    if ( sum(abs(invariant_operations(:,:,i)-ops(:,:,j))) .lt. tol ) then
                        cycle opl
                    endif
                enddo
                ctr=ctr+1
                ops(:,:,ctr)=invariant_operations(:,:,i)
            enddo opl
            ! Now build the invariant matrix
            IM=0.0_r8
            do i=1,ctr
                IM=IM+ops(:,:,i)
            enddo
            do i=1,n
                IM(i,i)=IM(i,i)-ctr
            enddo
            deallocate(ops)
        endif

        ! Check for trivial solution
        if ( sum(abs(IM)) .lt. tol ) then
            allocate(coeff(n,n))
            coeff=0.0_r8
            do i=1,n
                coeff(i,i)=1.0_r8
            enddo
            nvar=n
            if ( present(varind) ) then
                allocate(varind(n))
                do i=1,n
                    varind(i)=i
                enddo
            endif
            ! and return early
            return
        endif
    end block init

    slv: block
        complex(r8), dimension(:,:), allocatable :: U,V,dU,dV
        real(r8), dimension(:), allocatable :: S
        integer, dimension(:), allocatable :: zeroSV
        integer :: nc,i,j

        call lo_complex_singular_value_decomposition(IM,S,U,V)

        allocate(zeroSV(n))
        zeroSV=0
        ! Count the number of zero singular values
        nc=0
        do i=1,n
            if ( abs(S(i)) .lt. tol ) then
                nc=nc+1
                zeroSV(nc)=i
            endif
        enddo

        ! now return the irreducible represantation.
        if ( nc .eq. 0 ) then
            ! There could be zero independent things. Or? Yes it could.
            nvar=0
        elseif ( nc .eq. n ) then
            ! It could also be that the invarM did not help at all.
            allocate(coeff(n,n))
            coeff=0.0_r8
            do i=1,n
                coeff(i,i)=1.0_r8
            enddo
            nvar=n
            if ( present(varind) ) then
                allocate(varind(n))
                do i=1,n
                    varind(i)=i
                enddo
            endif
        else
            ! If I made it here, it means it's a non-trivial solution.
            ! I reuse the Sigma array, now I it's the full-rank nullspace
            allocate(dU(n,nc))
            allocate(dV(nc,n))
            dU=U(:,zeroSV(1:nc))
            dV=V(zeroSV(1:nc),:)
            IM=0.0_r8
            call zgemm('N','N',n,n,nc,cmplx(1.0_r8,0.0_r8,r8),dU,n,dV,nc,cmplx(0.0_r8,0.0_r8,r8),IM,n)
            deallocate(dU)
            deallocate(dV)
            ! get the coefficient matrix to reduced row echelon form
            call lo_make_complex_coeffmatrix_tidy(IM,tol)

            ! Now return the non-zero parts of this
            nvar=nc
            allocate(coeff(n,nvar))
            coeff=0.0_r8
            if ( present(varind) ) then
                allocate(varind(nvar))
                varind=0
            endif
            j=0
            do i=1,n
                if ( sum(abs(IM(:,i))) .gt. tol ) then
                    j=j+1
                    if ( present(varind) ) varind(j)=i
                    coeff(:,j)=IM(:,i)
                endif
            enddo
        endif
        deallocate(IM)
    end block slv
end subroutine

!> Solve for the nullspace of a matrix, used to get irreducible representations. Not fast, but that should not matter.
module subroutine lo_real_nullspace_coefficient_matrix(invarM,invariant_operations,coeff,nvar,varind,tolerance)
    !> matrix that leaves the desired quantity invariant
    real(r8), dimension(:,:), intent(in), optional :: invarM
    !> list of operations
    real(r8), dimension(:,:,:), intent(in), optional :: invariant_operations
    !> coefficient matrix
    real(r8), dimension(:,:), allocatable, intent(out) :: coeff
    !> number of independent variables
    integer, intent(out) :: nvar
    !> where are the irreducible in terms of the full
    integer, dimension(:), allocatable, intent(out), optional :: varind
    !> tolerance
    real(r8), intent(in), optional :: tolerance

    real(r8), dimension(:,:), allocatable :: IM
    real(r8) :: tol
    integer :: n

    ! Sort out the heuristics
    init: block
        real(r8), dimension(:,:,:), allocatable :: ops
        integer :: i,j,ctr

        ! set the tolerance
        if ( present(tolerance) ) then
            tol=tolerance
        else
            tol=lo_sqtol
        endif

        ! check input
        if ( present(invariant_operations) .and. present(invarM) ) then
            call lo_stop_gracefully(['Choose either invarM or invariant_operations as input, not both'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! Get the invariant matrix, either as input or from the operations
        if ( present(invarM) ) then
            ! First check that we don't have conflicting input
            n=size(invarM,1)
            allocate(IM(n,n))
            IM=lo_chop(invarM,tol)
        elseif ( present(invariant_operations) ) then
            ! Be extra careful and don't add redundant operations.
            n=size(invariant_operations,1)
            i=size(invariant_operations,3)
            allocate(IM(n,n))
            allocate(ops(n,n,i))
            ops=-1.0_r8
            ctr=0
            opl: do i=1,size(invariant_operations,3)
                do j=1,ctr
                    if ( sum(abs(invariant_operations(:,:,i)-ops(:,:,j))) .lt. tol ) then
                        cycle opl
                    endif
                enddo
                ctr=ctr+1
                ops(:,:,ctr)=invariant_operations(:,:,i)
            enddo opl
            ! Now build the invariant matrix
            IM=0.0_r8
            do i=1,ctr
                IM=IM+ops(:,:,i)
            enddo
            do i=1,n
                IM(i,i)=IM(i,i)-ctr
            enddo
            deallocate(ops)
        else
            call lo_stop_gracefully(['Provide either invarM or invariant_operations as input.'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! Check for trivial solution
        if ( sum(abs(IM)) .lt. tol ) then
            allocate(coeff(n,n))
            coeff=0.0_r8
            do i=1,n
                coeff(i,i)=1.0_r8
            enddo
            nvar=n
            if ( present(varind) ) then
                allocate(varind(n))
                do i=1,n
                    varind(i)=i
                enddo
            endif
            ! and return early
            return
        endif
    end block init

    slv: block
        real(r8), dimension(:,:), allocatable :: U,V,dU,dV
        real(r8), dimension(:), allocatable :: S
        integer, dimension(:), allocatable :: zeroSV
        integer :: nc,i,j

        call lo_real_singular_value_decomposition(IM,S,U,V)

        allocate(zeroSV(n))
        zeroSV=0
        ! Count the number of zero singular values
        nc=0
        do i=1,n
            if ( abs(S(i)) .lt. tol ) then
                nc=nc+1
                zeroSV(nc)=i
            endif
        enddo

        ! now return the irreducible represantation.
        if ( nc .eq. 0 ) then
            ! There could be zero independent things. Or? Yes it could.
            nvar=0
        elseif ( nc .eq. n ) then
            ! It could also be that the invarM did not help at all.
            allocate(coeff(n,n))
            coeff=0.0_r8
            do i=1,n
                coeff(i,i)=1.0_r8
            enddo
            nvar=n
            if ( present(varind) ) then
                allocate(varind(n))
                do i=1,n
                    varind(i)=i
                enddo
            endif
        else
            ! If I made it here, it means it's a non-trivial solution.
            ! I reuse the Sigma array, now it's the full-rank nullspace
            allocate(dU(n,nc))
            allocate(dV(nc,n))
            dU=U(:,zeroSV(1:nc))
            dV=V(zeroSV(1:nc),:)
            IM=0.0_r8
            call dgemm('N','N',n,n,nc,1.0_r8,dU,n,dV,nc,0.0_r8,IM,n)
            deallocate(dU)
            deallocate(dV)
            ! get the coefficient matrix to reduced row echelon form
            call lo_make_coeffmatrix_tidy(IM,tol)

            ! Now return the non-zero parts of this
            nvar=nc
            allocate(coeff(n,nvar))
            coeff=0.0_r8
            if ( present(varind) ) then
                allocate(varind(nvar))
                varind=0
            endif
            j=0
            do i=1,n
                if ( sum(abs(IM(:,i))) .gt. tol ) then
                    j=j+1
                    if ( present(varind) ) varind(j)=i
                    coeff(:,j)=IM(:,i)
                endif
            enddo
        endif
    end block slv
end subroutine

!> create matrix representation of a transposition
pure module subroutine lo_transpositionmatrix(tm)
    !> transposition matrix
    real(r8), dimension(:,:), intent(out) :: tm

    integer :: i,j,k,l,n

    tm=0.0_r8
    n=lo_intsqrt(size(tm,1))
    l=0
    do i=1,n
    do j=1,n
        l=l+1
        k=(j-1)*n+i
        tm(l,k)=1.0_r8
    enddo
    enddo
end subroutine

!> wrapper around LAPACK so that it is not destructive. Solves Ax=B, perhaps with constraints such that Cx=0, with a bunch of fancy options.
module subroutine lo_linear_least_squares(A,B,x,constraint_C,constraint_D,nconstraints,tolerance,subset,weights,gramified)
    !> matrix A
    real(r8), dimension(:,:), intent(in) :: A
    !> vector B
    real(r8), dimension(:), intent(in) :: B
    !> solution
    real(r8), dimension(:), intent(out) :: x
    !> linear constraints
    real(r8), dimension(:,:), intent(in), optional :: constraint_C
    !> D in Cx=D
    real(r8), dimension(:), intent(in), optional :: constraint_D
    !> number of constraints
    integer, intent(in), optional :: nconstraints
    !> specify a tolerance
    real(r8), intent(in), optional :: tolerance
    !> specify a subset of variables to consider
    integer, intent(in), dimension(:), optional :: subset
    !> optional weights, for weighted least squares
    real(r8), dimension(:), intent(in), optional :: weights
    !> am I solving A^T A X = A^T B instead of AX=B?
    logical, intent(in), optional :: gramified

    real(r8), dimension(:,:), allocatable :: wA,wC
    real(r8), dimension(:), allocatable :: sol,vB
    real(r8) :: tol
    integer, dimension(:), allocatable :: xind
    integer :: slvtype,nx,nc
    logical :: sub,gram

    init: block
        real(r8), dimension(:,:), allocatable :: vC
        integer :: rankA,i,j,l,n

        ! First thing to do is figure out a whole bunch of heuristics. Like how do I want to solve it, what variables to solve for,
        ! add a ridge-penalty, add constraints and how the constraints should be added. A little messy, but should be robust, at
        ! least more robust than just a plain application.

        ! Initial sanity test for the constraints
        if ( present(constraint_C) ) then
            if ( .not.present(nconstraints) ) then
                call lo_stop_gracefully(['Need to provide the number of constraints to the llsq'],lo_exitcode_param,__FILE__,__LINE__)
            elseif ( nconstraints .gt. 0 ) then
                if ( nconstraints .ne. size(constraint_C,1) ) then
                    call lo_stop_gracefully(['Inconsistent number of constraints'],lo_exitcode_param,__FILE__,__LINE__)
                endif
            endif
            nc=nconstraints
        else
            nc=0
        endif

        ! We can not have both a subset and gramian form, does not work.
        if ( present(gramified) .and. present(subset) ) then
            call lo_stop_gracefully(['Can not solve for a subset with the gramian form.'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Sanity test, I can not use weights and the Gramian form at the same time.
        if ( present(gramified) .and. present(weights) ) then
            call lo_stop_gracefully(['Can not use weights with the gramian form.'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Sanity test, double-check the dimensions?
        if ( size(A,2) .ne. size(x) ) then
            call lo_stop_gracefully(['Dimension mismatch for linear least squares.'],lo_exitcode_baddim,__FILE__,__LINE__)
        endif
        if ( size(A,1) .ne. size(B) ) then
            call lo_stop_gracefully(['Dimension mismatch for linear least squares.'],lo_exitcode_baddim,__FILE__,__LINE__)
        endif

        ! Figure out if I only want to solve for a subset of the variables.
        if ( present(subset) ) then
            ! I only count indices > 0 as something to care about.
            l=count(subset>0)
            if ( l .eq. 0 ) then
                call lo_stop_gracefully(['Makes no sense to solve for empty subset.'],lo_exitcode_param,__FILE__,__LINE__)
            endif
            if ( l .gt. size(A,2) ) then
                call lo_stop_gracefully(['Subset has to be smaller than the actual set. Duh.'],lo_exitcode_param,__FILE__,__LINE__)
            endif
            allocate(xind(l))
            xind=0
            l=0
            do i=1,size(subset)
                if ( subset(i) .le. 0 ) cycle
                l=l+1
                xind(l)=subset(i)
            enddo
            ! If the subset is the full set, I should switch the subset thing off.
            if ( l .lt. size(A,2) ) then
                sub=.true.
                nx=size(xind)
            else
                sub=.false.
                nx=size(A,2)
                deallocate(xind)
            endif
        else
            ! Not a subset
            sub=.false.
            ! Now I know the number of variables
            nx=size(A,2)
        endif

        ! Are we solving the problem transformed to Gramian matrix form?
        if ( present(gramified) ) then
            gram=gramified
            ! If there are no constraints, we don't have to care about
            ! wether it's on the gramian form or not.
            if ( nc .eq. 0 ) gram=.false.
        else
            gram=.false.
        endif

        ! Now make a copy of the input coefficient matrix to not destroy it.
        ! I might either make a straight copy, or just the relevant subset.
        if ( sub ) then
            ! Make a smaller matrix with just the rows I need
            allocate(wA(size(A,1),size(xind)))
            do i=1,size(xind)
                wA(:,i)=A(:,xind(i))
            enddo
            ! Constraints are tricky with a subset. Test a little. An alternative is to
            ! always turn them off, might do that in the future.
            if ( nc .gt. 0 ) then
                allocate(vC(nc,size(xind)))
                vC=0.0_r8
                do i=1,size(xind)
                    vC(:,i)=constraint_C(:,xind(i))
                enddo
                call lo_compress_equations(vC,l,wC,trans=.false.,tolerance=lo_sqtol)
                if ( l .ne. nc ) then
                    ! If I messed with the constraints, better to turn them off and fix afterwards.
                    nc=0
                    ! Kill the constraint matrix to make sure it is never used
                    deallocate(wC)
                endif
                deallocate(vC)
            endif
            ! And a copy of the prediction vector
            allocate(vB(size(B)))
            vB=B
        elseif ( gram ) then
            ! Make a copy, and pad it with the constraints
            n=nx+nc
            allocate(wA(n,n))
            allocate(vB(n))
            ! first copy the response vector
            vB=0.0_r8
            vB(1:nx)=B
            do i=1,nc
                vB(nx+i)=constraint_D(i)
            enddo
            ! then the full thingy
            wA=0.0_r8
            wA(1:nx,1:nx)=A
            do j=1,nx
            do i=1,nc
                l=nx+i
                wA(l,j)=constraint_C(i,j)
                wA(j,l)=constraint_C(i,j)
            enddo
            enddo
        else
            ! Just a straight copy
            allocate(wA(size(A,1),size(A,2)))
            wA=A
            if ( nc .gt. 0 ) then
                allocate(wC(nc,size(A,2)))
                wC=constraint_C
            endif
            ! And a copy of the prediction vector
            allocate(vB(size(B)))
            vB=B
        endif

        ! Some sanity checks
        if ( size(wA,1) .ne. size(vB) ) then
            call lo_stop_gracefully(['Incompatible dimensions for least squares'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! Add weights?
        if ( present(weights) ) then
            do i=1,size(wA,1)
                wA(i,:)=wA(i,:)*weights(i)
                vB(i)=vB(i)*weights(i)
            enddo
        endif

        ! Set the tolerance
        if ( present(tolerance) ) then
            tol=tolerance
        else
            tol=lo_sqtol
        endif

        ! Now there are a few possibilities how we can solve this. Getting kind of messy.
        if ( nc .gt. 0 ) then
            rankA=min(nx,size(wA,1))
            if ( rankA .ge. nx .and. gram .eqv. .false. ) then
                slvtype=2
            else
                slvtype=1
            endif
        else
            slvtype=1
        endif

        ! And some space for the solution
        allocate(sol(nx))
        sol=0.0_r8
        X=0.0_r8
    end block init

    ! Now actually solve it. A few variants.
    select case(slvtype)
    case(1)
    ! This is the normal, SVD-based least squares solution with no constraints.
    sl1: block
        real(r8), dimension(:,:), allocatable :: wB
        real(r8), dimension(:), allocatable :: work,o_s
        real(r8) :: s_work(1),o_rcond
        integer :: nrhs,lwork,o_rank,lda,ldb,m,n

        nrhs = 1
        lwork = -1
        m = size(wA,1)
        n = size(wA,2)
        o_rcond=100*epsilon(1.0_r8)
        allocate(wB(max(m,n),1))
        allocate(o_s(min(m,n)))
        o_s=0.0_r8
        wB=0.0_r8
        wB(1:m,1)=vB
        lda = max(1,size(wA,1))
        ldb = max(1,size(wB,1))
        call dgelss(m,n,nrhs,wA,lda,wB,ldb,o_s,o_rcond,o_rank,s_work,lwork,lo_status)
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['dgelss exit status'//tochar(lo_status)],lo_exitcode_param,__FILE__,__LINE__)
        endif
        lwork = int(anint( s_work(1)))
        allocate(work(lwork))
        work=0.0_r8
        call dgelss(m,n,nrhs,wA,lda,wB,ldb,o_s,o_rcond,o_rank,work,lwork,lo_status)
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['dgelss exit status'//tochar(lo_status)],lo_exitcode_param,__FILE__,__LINE__)
        endif
        sol=wB(1:nx,1)
        if ( sub ) then
            x(xind)=sol
        else
            x=sol
        endif
        deallocate(wA)
        deallocate(wB)
        deallocate(o_s)
        deallocate(work)
    end block sl1
    case(2)
    ! Overdetermined system, with constraints
    sl2: block
        real(r8), dimension(:), allocatable :: vA,vD,work
        real(r8) :: s_work(1)
        integer :: lwork,ne
        ne=size(wA,1)
        allocate(vA(ne))
        allocate(vD(nx))
        allocate(sol(nx))
        vA=vB
        vD=0.0_r8
        sol=0.0_r8
        lwork = -1
        call dgglse(ne,nx,nc,wA,ne,wC,nc,vA,vD,sol,s_work,lwork,lo_status)
        if ( lo_status .eq. 0 ) then
            lwork = int(anint(s_work(1)))
            allocate(work(lwork))
            call dgglse(ne,nx,nc,wA,ne,wC,nc,vA,vD,sol,work,lwork,lo_status)
            deallocate(work)
        endif
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['dgglse exit status '//tochar(lo_status)],lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( sub ) then
            x(xind)=sol
        else
            x=sol
        endif
        deallocate(vA)
        deallocate(vD)
        deallocate(sol)
    end block sl2
    end select

    ! Now for the final touch.
    finalize: block
        !real(r8), dimension(:), allocatable :: viol

        ! ! Check if the constraints are satisfied
        ! if ( present(zeroconstraints) ) then
        ! if ( nconstraints .gt. 0 ) then
        !     allocate(viol(size(zeroconstraints,2)))
        !     viol=0.0_r8
        !     ! call dgemv('N',size(zeroconstraints,1),size(zeroconstraints,2),1.0_r8,&
        !     !            zeroconstraints,size(zeroconstraints,1),X,1,0.0_r8,viol,1)
        !     ! Check if they are not satisfied. If so, fix it.
        !     if ( sum(abs(viol)) .gt. tol ) then
        !         write(*,*) 'Constraints not satisfied!'
        !         stop
        !         !call lo_enforce_linear_constraints(zeroconstraints,x,tol,nosquare=.true.)
        !     endif
        !     deallocate(viol)
        ! endif
        ! endif
        ! Chop away small numbers
        X=lo_chop(X,tol*1E-2_r8)
    end block finalize
end subroutine

!> Real symmetric eigenvalue problem, with careful checks and some neat options
module subroutine lo_real_symmetric_eigenvalues_eigenvectors(A,eigenvalues,eigenvectors,careful,tolerance,nzeros)
    !> matrix to be decomposed
    real(r8), dimension(:,:), intent(in) :: A
    !> eigenvalues
    real(r8), dimension(:), intent(out) :: eigenvalues
    !> eigenvectors
    real(r8), dimension(:,:), intent(out), optional :: eigenvectors
    !> check things carefully
    logical, intent(in), optional :: careful
    !> numerical tolerance
    real(r8), intent(in), optional :: tolerance
    !> make sure there are at least N eigenvalues that are exactly zero
    integer, intent(in), optional :: nzeros

    real(r8), dimension(:,:), allocatable :: wA
    real(r8) :: tol
    integer :: n
    logical :: check

    init: block
        real(r8) :: f0
        ! How careful should I check
        if ( present(careful) ) then
            check=careful
        else
            check=.false.
        endif
        ! Size of problem
        n=size(A,1)

        ! Set the tolerance
        if ( present(tolerance) ) then
            tol=tolerance
        else
            ! Get a reasonable tolerance. Gerschgorin maybe? Only necessary if we do careful tests, but I might as well.
            if ( check ) then
                f0=sum(abs(A))/size(A)
                tol=10*f0*lo_sqtol
            else
                tol=lo_sqtol
            endif
        endif

        ! Zero out the things to be returned
        eigenvalues=0.0_r8
        if ( present(eigenvectors) ) then
            eigenvectors=0.0_r8
        endif
        ! And a copy of the matrix
        allocate(wA(n,n))
        wA=A
    end block init

    ! Actually diagonalize
    diag: block
        ! Dummy things for lapack
        character(len=1) :: o_uplo,jobz
        real(r8), allocatable :: work(:),rwork(:)
        real(r8) :: s_work(1),o_z(1),s_rwork(1),o_vl,o_vu
        integer, allocatable :: iwork(:),o_isuppz(:)
        integer :: s_iwork(1),o_il,o_iu,o_m,o_info,lwork,lrwork,liwork,l_stat_alloc

        o_info=0
        o_uplo = 'U'
        o_vl = -lo_huge
        o_vu = lo_huge
        if( present(eigenvectors) ) then
            jobz = 'V'
        else
            jobz = 'N'
        endif
        o_iu = n
        l_stat_alloc = 0
        if( .not.present(eigenvectors) ) then
            allocate(o_isuppz(1))
            o_isuppz=0
        else
            allocate(o_isuppz(2*n), stat=l_stat_alloc)
        endif

        liwork = -1
        lrwork = -1
        lwork = -1
        if ( present(eigenvectors) ) then
            call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,n,0.0_r8,o_m,eigenvalues,&
                 eigenvectors,n,o_isuppz,s_work,lwork,s_iwork,liwork,o_info)
        else
            call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,n,0.0_r8,o_m,eigenvalues,&
                 o_z,1,o_isuppz,s_work,lwork,s_iwork,liwork,o_info)
        endif
        if (o_info == 0 ) then
            liwork = int(s_iwork(1))
            lrwork = int(s_rwork(1))
            lwork = int(s_work(1))
            allocate(iwork(liwork))
            allocate(rwork(lrwork))
            allocate(work(lwork))
            if ( present(eigenvectors) ) then
                call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,o_iu,0.0_r8,o_m,&
                     eigenvalues,eigenvectors,n,o_isuppz,work,lwork,iwork,liwork,o_info)
            else
                call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,o_iu,0.0_r8,o_m,&
                     eigenvalues,o_z,1,o_isuppz,work,lwork,iwork,liwork,o_info)
            endif
            deallocate(iwork)
            deallocate(rwork)
            deallocate(work)
        endif
        if ( allocated(o_isuppz) ) deallocate(o_isuppz)
        if( o_info <= -1000 ) then
            call xerbla('HEEVR',-o_info)
            call lo_stop_gracefully(['zheevr exit code '//tochar(-o_info)],lo_exitcode_blaslapack)
        endif
        ! Ok, now w should hold the eigenvalues, eigenvectors in eigenvectors.
    end block diag

    ! Put the final touches, and test some things to be on the safe side
    finalize: block
        real(r8), dimension(:), allocatable :: dr
        real(r8) :: f0
        integer, dimension(:), allocatable :: ind
        integer :: i,j
        logical :: sorted,orth

        ! First, I might have to force N of the eigenvalues to be zero.
        if ( present(nzeros) ) then
        if ( nzeros .gt. 0 ) then
            ! I will pick the N ones that are closest to zero in absolute values, I think.
            allocate(dr(n))
            allocate(ind(n))
            dr=abs(eigenvalues)
            call qsort(dr,ind)
            do i=1,nzeros
                eigenvalues(ind(i))=0.0_r8
            enddo
            deallocate(dr)
            deallocate(ind)
        endif
        endif

        ! Now test a bunch more things, to be on the safe side
        if ( check ) then
            ! Ok, now w should hold the eigenvalues, eigenvectors in eigenvectors. First test if it is sorted?
            sorted=.true.
            do i=1,n-1
                if ( eigenvalues(i) .gt. eigenvalues(i+1) ) then
                    sorted=.false.
                    exit
                endif
            enddo
            if ( sorted .eqv. .false. ) then
                allocate(ind(n))
                call qsort(eigenvalues,ind)
                if ( present(eigenvectors) ) then
                    wA=eigenvectors
                    do i=1,n
                        eigenvectors(:,i)=wA(:,ind(i))
                    enddo
                endif
                deallocate(ind)
            endif
            ! Ok, now eigenvalues and eigenvectors are sorted. Test something else?
            if ( present(eigenvectors) ) then
                ! Test if the eigenvectors are orthonormal?
                orth=.true.
                do i=1,n
                    f0=dot_product(eigenvectors(:,i),eigenvectors(:,i))
                    if ( abs(f0-1.0_r8) .gt. tol ) then
                        orth=.false.
                        exit
                    endif
                enddo
                if ( orth ) then
                    l1: do i=1,n
                    do j=i+1,n
                        f0=abs(dot_product(eigenvectors(:,i),eigenvectors(:,i)))
                        if ( f0 .gt. tol ) then
                            orth=.false.
                            exit l1
                        endif
                    enddo
                    enddo l1
                endif

                ! Orthogonalize
                if ( orth .eqv. .false. ) then
                    call lo_real_gram_schmidt(eigenvectors)
                endif
            endif
        endif
        ! And chop off annoying things
        eigenvalues=lo_chop(eigenvalues,tol)
        if ( present(eigenvectors) ) eigenvectors=lo_chop(eigenvectors,tol)
    end block finalize
end subroutine

!> Solve 3x3 symmetric eigenproblem quite fast.
module subroutine lo_symmetric_eigensystem_3x3matrix(matrix,eigenvalues,eigenvectors)
    !> matrix
    real(r8), dimension(3,3), intent(in) :: matrix
    !> eigenvalues
    real(r8), dimension(3), intent(out) :: eigenvalues
    !> eigenvectors
    real(r8), dimension(3,3), intent(out) :: eigenvectors

    real(r8), dimension(3) :: e

    tridiag: block
        real(r8), dimension(3) :: u,p
        real(r8) :: omega,f,k,h,g
        integer :: i,j
        ! initialize q to the identitity matrix
        eigenvectors=reshape([1.0_r8,0.0_r8,0.0_r8,0.0_r8,1.0_r8,0.0_r8,0.0_r8,0.0_r8,1.0_r8],[3,3])
        ! bring first row and column to the desired form
        h = matrix(1,2)**2 + matrix(1,3)**2
        if ( matrix(1,2) .gt. 0.0_r8 ) then
            g = -sqrt(h)
        else
            g = sqrt(h)
        endif
        e(1)  = g
        f     = g * matrix(1,2)
        u(2)  = matrix(1,2) - g
        u(3)  = matrix(1,3)
        omega = h - f
        if ( omega .gt. 0.0_r8 ) then
            omega = 1.0_r8/omega
            k     = 0.0d0
            do i=2,3
                f = matrix(2,i)*u(2) + matrix(i,3)*u(3)
                p(i) = omega * f
                k    = k + u(i) * f
            enddo
            k=0.5_r8*k*omega**2
            do i=2,3
                p(i)=p(i)-k*u(i)
            enddo
            eigenvalues(1) = matrix(1,1)
            eigenvalues(2) = matrix(2,2) - 2.0_r8 * p(2) * u(2)
            eigenvalues(3) = matrix(3,3) - 2.0_r8 * p(3) * u(3)
            ! store inverse householder transformation in q
            do j = 2, 3
                f=omega*u(j)
                do i = 2, 3
                    eigenvectors(i,j) = eigenvectors(i,j) - f * u(i)
                enddo
            enddo
            ! calculated updated a(2, 3) and store it in e(2)
            e(2) = matrix(2,3) - p(2)*u(3) - u(2)*p(3)
        else
            do i=1,3
                eigenvalues(i) = matrix(i,i)
                e(2) = matrix(2,3)
            enddo
        endif
    end block tridiag

    ! Calculate eigensystem of the remaining real symmetric tridiagonal
    ! matrix with the QL method
    ql: block
        real(r8) :: g, r, p, f, b, s, c, t
        integer :: l, m, i, j, k

        ! Loop over all off-diagonal elements
        outloop: do l=1,2
        inloop: do i=1,50
            ! Check for convergence and exit iteration loop if off-diagonal
            ! element E(L) is zero
            do m=l,2
                g = abs(eigenvalues(m)) + abs(eigenvalues(m+1))
                if ( abs(e(m)) + g .eq. g ) exit
            enddo
            if ( m .eq. l ) cycle outloop
            ! Calculate G = D(M) - K
            g = (eigenvalues(l+1) - eigenvalues(l)) / (2.0_r8 * e(l))
            r = sqrt(1.0_r8 + g**2)
            if ( g .ge. 0.0_r8 ) then
                g = eigenvalues(m) - eigenvalues(l) + e(l)/(g + r)
            else
                g = eigenvalues(m) - eigenvalues(l) + e(l)/(g - r)
            end if
            S = 1.0_r8
            C = 1.0_r8
            P = 0.0_r8
            do j=m-1,l,-1
                f = s*e(j)
                b = c*e(j)
                if ( abs(f) .gt. abs(g) ) then
                    c      = g / f
                    r      = sqrt(1.0_r8 + c**2)
                    e(j+1) = f * r
                    s      = 1.0_r8/r
                    c      = c * s
                else
                    s      = f / g
                    r      = sqrt(1.0_r8 + s**2)
                    e(j+1) = g * r
                    c      = 1.0_r8 / r
                    s      = s * c
                end if
                g      = eigenvalues(j+1) - p
                r      = (eigenvalues(j) - g) * s + 2.0_r8 * c * b
                p      = s * r
                eigenvalues(j+1) = g + p
                g      = c * r - b
                ! Form eigenvectors
                do k=1,3
                    t        = eigenvectors(k,j+1)
                    eigenvectors(k,j+1) = s*eigenvectors(k,j)+c*t
                    eigenvectors(k,j)   = c*eigenvectors(k,j)-s*t
                enddo
            enddo
            eigenvalues(l) = eigenvalues(l) - p
            e(l) = g
            e(m) = 0.0_r8
        enddo inloop
        enddo outloop
    end block ql
end subroutine

!> Complex eigenvalue problem, with careful checks and some neat options
module subroutine lo_complex_hermitian_eigenvalues_eigenvectors(A,eigenvalues,eigenvectors,careful,tolerance,nzeros)
    !> matrix to be decomposed
    complex(r8), dimension(:,:), intent(in) :: A
    !> eigenvalues
    real(r8), dimension(:), intent(out) :: eigenvalues
    !> eigenvectors
    complex(r8), dimension(:,:), intent(out), optional :: eigenvectors
    !> check things carefully
    logical, intent(in), optional :: careful
    !> numerical tolerance
    real(r8), intent(in), optional :: tolerance
    !> make sure there are at least N eigenvalues that are exactly zero
    integer, intent(in), optional :: nzeros

    complex(r8), dimension(:,:), allocatable :: wA
    real(r8) :: tol
    integer :: n
    logical :: check

    ! Some heuristics
    init: block
        complex(r8) :: c0
        real(r8) :: f0
        integer :: i,j
        ! How careful should I check
        if ( present(careful) ) then
            check=careful
        else
            check=.false.
        endif
        ! Size of problem
        n=size(A,1)

        ! Set the tolerance
        if ( present(tolerance) ) then
            tol=tolerance
        else
            tol=lo_sqtol
            ! Get a reasonable tolerance. Gerschgorin maybe? Only necessary if we do careful tests, but I might as well.
            if ( check ) then
                f0=sum(abs(A))/size(A)
                tol=10*f0*lo_sqtol
            else
                tol=lo_sqtol
            endif
        endif

        ! Zero out the things to be returned
        eigenvalues=0.0_r8
        if ( present(eigenvectors) ) then
            eigenvectors=0.0_r8
        endif
        ! And a copy of the matrix
        allocate(wA(n,n))
        if ( check ) then
            do i=1,n
            do j=i,n
                c0=(A(j,i)+conjg(A(i,j)))*0.5_r8
                wA(j,i)=c0
                wA(i,j)=conjg(c0)
            enddo
            enddo
        else
            wA=A
        endif

    end block init

    ! Actually diagonalize things
    hermdiag: block
        ! Below uses zheev
        !complex(r8), dimension(:), allocatable :: work
        !complex(r8) :: s_work(1)
        !real(r8), dimension(:), allocatable :: rwork
        !integer :: lwork
        !allocate(rwork(max(1,3*n-2)))
        !rwork=0.0_r8
        !lwork = -1
        !if ( present(eigenvectors) ) then
        !    call zheev('V','U',n,wA,n,eigenvalues,s_work,lwork,rwork,lo_status)
        !else
        !    call zheev('N','U',n,wA,n,eigenvalues,s_work,lwork,rwork,lo_status)
        !endif
        !if ( lo_status .ne. 0 ) then
        !    call lo_stop_gracefully(['zheev exit code '//tochar(-lo_status)],lo_exitcode_blaslapack)
        !endif
        !!
        !lwork = int(anint(real(s_work(1))))
        !allocate(work(lwork))
        !if ( present(eigenvectors) ) then
        !    call zheev('V','U',n,wA,n,eigenvalues,work,lwork,rwork,lo_status)
        !    eigenvectors=wA
        !else
        !    call zheev('N','U',n,wA,n,eigenvalues,work,lwork,rwork,lo_status)
        !endif
        !if ( lo_status .ne. 0 ) then
        !    call lo_stop_gracefully(['zheev exit code '//tochar(-lo_status)],lo_exitcode_blaslapack)
        !endif
        !deallocate(work)
        !deallocate(rwork)
        ! Below is zheevr, should be better than zheev they say
        character(len=1) :: o_uplo,jobz
        complex(r8), allocatable :: work(:)
        complex(r8) :: s_work(1),o_z(1)
        real(r8), allocatable :: rwork(:)
        real(r8) :: s_rwork(1),o_vl,o_vu
        integer, allocatable :: iwork(:),o_isuppz(:)
        integer :: s_iwork(1),o_il,o_iu,o_m,o_info,lwork,lrwork,liwork,l_stat_alloc

        ! LAPACK stuffs
        o_info=0
        o_uplo = 'U'
        o_vl = -lo_huge
        o_vu = lo_huge
        if( present(eigenvectors) ) then
            jobz = 'V'
        else
            jobz = 'N'
        endif
        o_iu = n
        l_stat_alloc = 0
        if( .not.present(eigenvectors) ) then
            allocate(o_isuppz(1))
            o_isuppz=0
        else
            allocate(o_isuppz(2*n), stat=l_stat_alloc)
        endif

        liwork = -1
        lrwork = -1
        lwork = -1
        if ( present(eigenvectors) ) then
            call zheevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,n,0.0_r8,o_m,eigenvalues,eigenvectors,&
                        n,o_isuppz,s_work,lwork,s_rwork,lrwork,s_iwork,liwork,o_info)
        else
            call zheevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,n,0.0_r8,o_m,eigenvalues,o_z,&
                        1,o_isuppz,s_work,lwork,s_rwork,lrwork,s_iwork,liwork,o_info)
        endif
        if (o_info == 0 ) then
            liwork = int(s_iwork(1))
            lrwork = int(s_rwork(1))
            lwork = int(s_work(1))
            allocate(iwork(liwork))
            allocate(rwork(lrwork))
            allocate(work(lwork))
            if ( present(eigenvectors) ) then
                call zheevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,o_iu,0.0_r8,o_m,eigenvalues,eigenvectors,&
                            n,o_isuppz,work,lwork,rwork,lrwork,iwork,liwork,o_info)
            else
                call zheevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,o_iu,0.0_r8,o_m,eigenvalues,o_z,&
                            1,o_isuppz,work,lwork,rwork,lrwork,iwork,liwork,o_info)
            endif
            deallocate(iwork)
            deallocate(rwork)
            deallocate(work)
        endif
        if ( allocated(o_isuppz) ) deallocate(o_isuppz)
        if( o_info <= -1000 ) then
            call xerbla('HEEVR',-o_info)
            call lo_stop_gracefully(['zheevr exit code '//tochar(-o_info)],lo_exitcode_blaslapack)
        endif
    end block hermdiag

    ! Put the final touches, and test some things to be on the safe side
    finalize: block
        real(r8), dimension(:), allocatable :: dr
        real(r8) :: f0
        integer, dimension(:), allocatable :: ind
        integer :: i,j
        logical :: sorted,orth

        ! First, I might have to force N of the eigenvalues to be zero.
        if ( present(nzeros) ) then
        if ( nzeros .gt. 0 ) then
            ! I will pick the N ones that are closest to zero in absolute values, I think.
            allocate(dr(n))
            allocate(ind(n))
            dr=abs(eigenvalues)
            call qsort(dr,ind)
            do i=1,nzeros
                eigenvalues(ind(i))=0.0_r8
            enddo
            deallocate(dr)
            deallocate(ind)
        endif
        endif

        ! Now test a bunch more things, to be on the safe side
        if ( check ) then
            ! Ok, now w should hold the eigenvalues, eigenvectors in eigenvectors. First test if it is sorted?
            ! Also how do you sort complex numbers? I'll sort by the real part, for some reason?
            sorted=.true.
            do i=1,n-1
                if ( eigenvalues(i) .gt. eigenvalues(i+1) ) then
                    sorted=.false.
                    exit
                endif
            enddo
            if ( sorted .eqv. .false. ) then
                allocate(ind(n))
                call qsort(eigenvalues,ind)
                if ( present(eigenvectors) ) then
                    wA=eigenvectors
                    do i=1,n
                        eigenvectors(:,i)=wA(:,ind(i))
                    enddo
                endif
                deallocate(ind)
            endif
            ! Ok, now eigenvalues and eigenvectors are sorted. Test something else?
            if ( present(eigenvectors) ) then
                ! Test if the eigenvectors are orthonormal?
                orth=.true.
                do i=1,n
                    f0=abs(dot_product(eigenvectors(:,i),eigenvectors(:,i)))
                    if ( abs(f0-1.0_r8) .gt. tol ) then
                        orth=.false.
                        exit
                    endif
                enddo
                if ( orth ) then
                    l1: do i=1,n
                    do j=i+1,n
                        f0=abs(dot_product(eigenvectors(:,i),eigenvectors(:,j)))
                        if ( abs(f0) .gt. tol ) then
                            orth=.false.
                            exit l1
                        endif
                    enddo
                    enddo l1
                endif

                ! Orthogonalize ( should never happen, I hope)
                if ( orth .eqv. .false. ) then
                    call lo_complex_gram_schmidt(eigenvectors)
                endif
            endif
        endif

        ! And chop off annoying things
        eigenvalues=lo_chop(eigenvalues,tol)
        if ( present(eigenvectors) ) eigenvectors=lo_chop(eigenvectors,tol)

    end block finalize
end subroutine

!> get left and right eigenvectors + eigenvalues of a general matrix
module subroutine lo_general_real_eigenvalues_eigenvectors(A,eigenvalues,vec_left,vec_right,orthogonal)
    !> matrix to check
    real(r8), dimension(:,:), intent(in) :: A
    !> eigenvalues
    complex(r8), dimension(:), intent(out) :: eigenvalues
    !> left and right eigenvectors. left are (:,i), right are also (:,i)
    real(r8), dimension(:,:), intent(out) :: vec_left,vec_right
    !> make sure it is orthogonal?
    logical, intent(in), optional :: orthogonal

    logical :: orth
    integer :: n

    ! Set some stuff
    init: block
        ! Set size of problem
        if ( size(A,1) .ne. size(A,2) ) then
            call lo_stop_gracefully(['Can only get eigenvectors of square matrices'],lo_exitcode_param,__FILE__,__LINE__)
        else
            n=size(A,1)
        endif
        vec_left=0.0_r8
        vec_right=0.0_r8
        eigenvalues=0.0_r8
        ! Should we make sure it is bi-orthogonal?
        if ( present(orthogonal) ) then
            orth=orthogonal
        else
            orth=.false.
        endif
    end block init

    ! First do the actual solution
    solve: block
        real(r8), dimension(:,:), allocatable :: rA
        real(r8), dimension(:), allocatable :: rwork,wr,wi
        real(r8), dimension(1) :: drwork
        integer :: lwork

        ! size of problem and space for stuff
        ! solve original eigenproblem, temporary space
        allocate(rA(n,n))
        allocate(wr(n))
        allocate(wi(n))
        rA=0.0_r8
        wr=0.0_r8
        wi=0.0_r8
        lwork=-1
        call dgeev('V','V', n, rA, n, wr, wi, vec_left, n, vec_right, n, drwork, lwork, lo_status)
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['dgeev exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        endif
        lwork=int(drwork(1))
        allocate(rwork(lwork))
        rwork=0.0_r8
        rA=A
        call dgeev('V','V', n, rA, n, wr, wi, vec_left, n, vec_right, n, rwork, lwork, lo_status)
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['dgeev exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        endif
        eigenvalues=cmplx(wr,wi,r8)
        deallocate(wr,wi,rA,rwork)
    end block solve

    ! then maybe orthogonalize
    if ( orth ) then
    ortho: block
        real(r8), dimension(:,:), allocatable :: wC,wU,wL,wA
        integer, dimension(:), allocatable :: ipiv
        integer :: i,j

        ! Attempt 1: Just do the whole thing and see what happens
        allocate(wC(n,n))
        allocate(wU(n,n))
        allocate(wL(n,n))
        allocate(wA(n,n))

        allocate(ipiv(n))
        ! Create VL*VR
        call dgemm('T','N',n,n,n,1.0_r8,vec_left,n,vec_right,n,0.0_r8,wA,n)
        ! LU factorize this
        wC=wA
        call dgetrf( n, n, wC, n, ipiv, lo_status )
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['zgetrf exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        endif
        ! check that it did not pivot. Don't think it ever should, in this particular case at least.
        do i=1,n
            if ( ipiv(i) .ne. i ) then
                call lo_stop_gracefully(['Think about pivoting'],lo_exitcode_blaslapack,__FILE__,__LINE__)
            endif
        enddo

        ! Extract upper and lower triangular guys. Not sure what to do if pivoting.
        ! Maybe pivot never happens. Can hope for that. Otherwise permute
        wU=0.0_r8
        wL=0.0_r8
        do i=1,n
            do j=i,n
                wU(i,j)=wC(i,j)
                wL(j,i)=wC(j,i)
            enddo
            wL(i,i)=1.0_r8
        enddo

        ! Check that I did this right
        call dgemm('N','N',n,n,n,1.0_r8,wU,n,wL,n,0.0_r8,wC,n)
        if ( sum(abs(wC-wA)) .gt. lo_sqtol ) then
            call lo_stop_gracefully(['Got the LU decomposition wrong'],lo_exitcode_blaslapack,__FILE__,__LINE__)
        endif

        ! Invert the triangular matrices
        call dtrtri( 'U', 'N', n, wU, n, lo_status )
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['dtrtri exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        endif
        call dtrtri( 'L', 'N', n, wL, n, lo_status )
        if ( lo_status .ne. 0 ) then
            call lo_stop_gracefully(['dtrtri exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        endif

        ! Correct the eigenvectors.
        call dgemm('N','T',n,n,n,1.0_r8,wU,n,vec_left,n,0.0_r8,wC,n)
        vec_left=transpose(wC)
        call dgemm('N','T',n,n,n,1.0_r8,vec_right,n,wL,n,0.0_r8,wC,n)
        vec_right=wC
    end block ortho
    endif
end subroutine

!> Matrix inversion. Returns the inverse in the matrix.
module subroutine lo_invert_real_matrix(A,iA)
    !> matrix to invert
    real(r8), dimension(:,:), intent(in) :: A
    !> inverse of matrix
    real(r8), dimension(:,:), intent(out) :: iA

    real(r8), dimension(:), allocatable :: work,ipiv
    integer :: n

    n=size(A,1)
    ! Some workspace
    allocate(work(n))
    allocate(ipiv(n))
    iA=A
    work=0.0_r8
    ipiv=0
    call dgetrf(n,n,iA,n,ipiv,lo_status)
    if ( lo_status .ne. 0 ) then
       call lo_stop_gracefully(['dgetrf exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    endif
    call dgetri(n,iA,n,ipiv,work,n,lo_status)
    if ( lo_status .ne. 0 ) then
       call lo_stop_gracefully(['dgetri exit code '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__)
    endif
    deallocate(work)
    deallocate(ipiv)
end subroutine

!> calculate the pseudoinverse via SVD
module subroutine lo_real_pseudoinverse(A,B,tolerance)
    !> matrix to pseudoinvert, n x m
    real(r8), dimension(:,:), intent(in) :: A
    !> pseudoinverse, m x n
    real(r8), dimension(:,:), intent(out) :: B
    !> optional tolerance
    real(r8), intent(in), optional :: tolerance

    real(r8), dimension(:,:), allocatable :: U,V,dm
    real(r8), dimension(:), allocatable :: S
    real(r8) :: tol
    integer :: i,n,m

    ! No idea how I chose the default here.
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=1E-12_r8
    endif
    n=size(A,1)
    m=size(A,2)
    allocate(dm(m,n))
    dm=0.0_r8
    call lo_real_singular_value_decomposition(A,S,U,V)
    ! Not sure about tolerance here. Think it's fine.
    do i=1,size(S)
        if ( abs(S(i)) .gt. tol ) then
            S(i)=1.0_r8/S(i)
        else
            S(i)=0.0_r8
        endif
    enddo
    do i=1,size(S)
        dm(i,:)=S(i)*U(:,i)
    enddo
    call dgemm('T','N',m,n,m,1.0_r8,V,m,DM,m,0.0_r8,B,m)
    B=lo_chop( B, 1E-13_r8 )
    ! Check that it actually worked? Nah. Cleanup instead.
    deallocate(dm)
    deallocate(S)
    deallocate(U)
    deallocate(V)
end subroutine

!> calculates D = A*B*C in one go
module subroutine lo_triplegemm(A,B,C,D,transa,transb,transc)
    !> matrices
    real(r8), dimension(:,:), intent(in) :: A,B,C
    !> output
    real(r8), dimension(:,:), intent(out) :: D
    !> transpose stuff?
    character(len=1), intent(in), optional :: transa,transb,transc

    real(r8), dimension(:,:), allocatable :: BC
    character(len=1) :: ta,tb,tc

    integer :: lda,ldb,ldc,ldd,ldbc
    integer :: ai,aj,bi,bj,ci,cj,di,dj

    ! Sort out all the ways things can be transposed
    if ( present(transa) ) then
        ta=transa
    else
        ta='N'
    endif
    if ( present(transb) ) then
        tb=transb
    else
        tb='N'
    endif
    if ( present(transc) ) then
        tc=transc
    else
        tc='N'
    endif

    ! Then all the funny dimensions, depending on how stuff got transposed
    if ( ta .eq. 'T' ) then
        ai=size(A,2)
        aj=size(A,1)
    else
        ai=size(A,1)
        aj=size(A,2)
    endif
    if ( tb .eq. 'T' ) then
        bi=size(B,2)
        bj=size(B,1)
    else
        bi=size(B,1)
        bj=size(B,2)
    endif
    if ( tc .eq. 'T' ) then
        ci=size(C,2)
        cj=size(C,1)
    else
        ci=size(C,1)
        cj=size(C,2)
    endif
    di=size(D,1)
    dj=size(D,2)

    ! Now I know enough to make temporary space
    allocate(BC(bi,cj))
    BC=0.0_r8
    lda=size(A,1)
    ldb=size(B,1)
    ldc=size(C,1)
    ldbc=size(BC,1)
    ldd=size(D,1)

    ! Begin operation actual operation
    call dgemm(tb,tc,bi,cj,bj,1.0_r8,B,ldb,C,ldc,0.0_r8,BC,ldbc)
    call dgemm(ta,'N',di,dj,aj,1.0_r8,A,lda,BC,ldbc,0.0_r8,D,ldd)
    deallocate(BC)
end subroutine

end submodule
