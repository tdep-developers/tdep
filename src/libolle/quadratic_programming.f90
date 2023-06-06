#include "precompilerdefinitions"
!> Solve standardized quadratic programming problems
module quadratic_programming
    use konstanter
    use gottochblandat
    !use mpi_wrappers ! Maybe remove
    use type_blas_lapack_wrappers
    implicit none

    private
    public :: lo_solve_quadratic_program

    interface
        module subroutine qpgen1(dmat, dvec, fddmat, n, sol, lagr, crval, amat, iamat,bvec, q, meq, iact, nact, iter, work, ierr)
            real(flyt), dimension(:,:), intent(inout) :: dmat
            real(flyt), dimension(:), intent(inout) :: dvec
            integer, intent(in) :: fddmat
            integer, intent(in) :: n
            real(flyt), dimension(:), intent(out) :: sol
            real(flyt), dimension(:), intent(out) :: lagr
            real(flyt), intent(inout) :: crval
            real(flyt), dimension(:,:), intent(inout) :: amat
            integer, dimension(:,:), intent(in) :: iamat
            real(flyt), dimension(:), intent(inout) :: bvec
            integer, intent(in) :: q
            integer, intent(in) :: meq
            integer, dimension(:), intent(out) :: iact
            integer, intent(inout) :: nact
            integer, dimension(2), intent(out) :: iter
            real(flyt), dimension(:), intent(inout) :: work
            integer, intent(inout) :: ierr
        end subroutine
    end interface

contains

!> Wrapper around the ancient solver.
subroutine lo_solve_quadratic_program(Q,c,x,A_equal,A_inequal,B_equal,B_inequal,n_equal,n_inequal,verbosity,tolerance)
    !> Hessian matrix Q. Must be positive definite.
    real(flyt), dimension(:,:), intent(in) :: Q
    !> Linear term c.
    real(flyt), dimension(:), intent(in) :: c
    !> Solution
    real(flyt), dimension(:), intent(out) :: x
    !> Equality and inequality constraints Ax=b, Ax >=b
    real(flyt), dimension(:,:), intent(in) :: A_equal,A_inequal
    real(flyt), dimension(:), intent(in) :: B_equal,B_inequal
    !> Number of constraints of each type. Set to zero to ignore.
    integer, intent(in) :: n_equal,n_inequal
    !> Talk a lot?
    integer, intent(in) :: verbosity
    !> Tolerance for when something is zero
    real(flyt), intent(in), optional :: tolerance

    integer :: nx,nc,nr
    real(flyt) :: timer,tol
    real(flyt), dimension(:,:), allocatable :: wA
    real(flyt), dimension(:), allocatable :: wb
    integer, dimension(:,:), allocatable :: iamat

    ! Start timer
    timer=walltime()
    ! Size of problem
    nx=size(Q,1)
    ! Tolerance
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=1E-11_flyt
    endif
    ! Total number of constraints
    nc=n_equal+n_inequal

    !@TODO Throw warnings/switch to linear solver if no inequality constraints, I suppose.
    !@TODO Check that C is positive definite
    !@TODO Check that constraints don't contradict each other
    !@TODO Check sizes of all arrays

    ! The solver likes the constraints in a sparse format, fix that.
    if ( nc .gt. 0 ) then
    packconstraints: block
        integer :: i,j,k,l

        ! First thing is to pack up the constraint A-matrices
        l=0
        do i=1,n_equal
            l=max(l,count( abs(A_equal(i,:))>tol) )
        enddo
        do i=1,n_inequal
            l=max(l,count( abs(A_inequal(i,:))>tol) )
        enddo
        ! Number of rows in the packed format, and total number of constraints
        nr=l
        nc=n_equal+n_inequal
        lo_allocate(iamat(nr+1,nc))
        lo_allocate(wA(nr,nc))
        iamat=0
        wA=0.0_flyt
        ! Actual packing
        k=0
        do i=1,n_equal
            k=k+1
            l=0
            do j=1,nx
                if ( abs(A_equal(i,j)) > tol ) then
                    l=l+1
                    iamat(l+1,k)=j
                    wA(l,k)=A_equal(i,j)
                endif
            enddo
            iamat(1,k)=l
        enddo
        do i=1,n_inequal
            k=k+1
            l=0
            do j=1,nx
                if ( abs(A_inequal(i,j)) > tol ) then
                    l=l+1
                    iamat(l+1,k)=j
                    wA(l,k)=A_inequal(i,j)
                endif
            enddo
            iamat(1,k)=l
        enddo
        ! Then stack the constraint B-vectors
        lo_allocate(wb(nc))
        k=0
        do i=1,n_equal
            k=k+1
            wb(k)=B_equal(i)
        enddo
        do i=1,n_inequal
            k=k+1
            wb(k)=B_inequal(i)
        enddo
    end block packconstraints
    endif

    solve: block
        ! Work arrays the solver needs
        real(flyt), dimension(:,:), allocatable :: wQ
        real(flyt), dimension(:), allocatable :: work,wc,wlagr
        real(flyt) :: crval
        integer, dimension(:), allocatable :: iact
        integer :: nact,ierr,iter(2),l
        ! Make space for work arrays
        lo_allocate(wQ(nx,nx))
        lo_allocate(wc(nx))
        wQ=Q
        wc=c
        if ( nc .gt. 0 ) then
            lo_allocate(wlagr(nc))
            lo_allocate(iact(nc))
            l=min(nx,nc)
            l=2*nx+l*(l+5)/2+2*nc+1+10 ! Additional 10 just to be sure.
            lo_allocate(work(l))
        else
            lo_allocate(wlagr(1))
            lo_allocate(iact(1))
            lo_allocate(work(2*nx+10))
        endif
        wlagr=0.0_flyt
        iact=0
        work=0.0_flyt
        ! Begin operation actual operation
        ierr=0
        call qpgen1(wQ,wc,nx,nx,X,wlagr,crval,wA,iamat,wb,nc,n_equal,iact,nact,iter,work,ierr)
        if ( ierr .ne. 0 ) then
            ! Did not go well
            call lo_stop_gracefully(['Quadratic programming failed, exit code '//tochar(ierr)],lo_exitcode_blaslapack,__FILE__,__LINE__)
        endif
    end block solve

    if ( verbosity .gt. 0 ) then
        write(*,*) 'Solved quadratic program (',tochar(walltime()-timer),'s)'
    endif

end subroutine

end module
