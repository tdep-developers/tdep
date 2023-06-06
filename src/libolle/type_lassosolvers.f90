#include "precompilerdefinitions"
module type_lassosolvers
use konstanter, only: flyt,lo_huge,lo_tiny,lo_status
use gottochblandat, only: walltime
use type_blas_lapack_wrappers, only: lo_dgemm,lo_dgemv,lo_dgels

implicit none
private
public :: lo_lasso_solve_GaussSiedel

contains

!> just the ridge regression, to start from or something
subroutine lo_ridgeregression(A,B,lambda,x)
    !> The coefficient matrix
    real(flyt), dimension(:,:), intent(in) :: A
    !> The observations
    real(flyt), dimension(:), intent(in) :: B
    !> The coupling parameter thing
    real(flyt), intent(in) :: lambda
    !> The solution
    real(flyt), dimension(:), intent(inout) :: x
    !
    real(flyt), dimension(:,:), allocatable :: AA,dum
    real(flyt), dimension(:), allocatable :: AB
    integer :: i,ni,nj
    !
    x=0.0_flyt
    ni=size(A,1) 
    nj=size(A,2)
    if ( size(x,1) .ne. nj ) then
        write(*,*) 'bad dimensions:'
        stop
    endif
    !
    lo_allocate(AA(nj,nj))
    AA=0.0_flyt
    do i=1,nj
        AA=1.0_flyt
    enddo
    ! AA = A^T * A + lambda*I
    call lo_dgemm(A,A,AA,transa='T',transb='N',alpha=1.0_flyt,beta=lambda)
    ! AB = A*B
    lo_allocate(AB(nj))
    AB=0.0_flyt
    call lo_dgemv(A,B,AB,trans='T')
    ! Solve stuff
    lo_allocate(dum(nj,1))
    dum(:,1)=AB
    call lo_dgels(AA,dum)
    x=dum(:,1)
    !
    lo_deallocate(AA)
    lo_deallocate(AB)
    lo_deallocate(dum)
    !
end subroutine

!> something for the GaussSiedel
subroutine compute_violation(alpha,gradient,viol,threshold,lambda)
    !> the guess?
    real(flyt), dimension(:), intent(in) :: alpha
    !> the gradient
    real(flyt), dimension(:), intent(in) :: gradient
    !> the violation
    real(flyt), dimension(:), intent(out) :: viol
    !> the tolerance
    real(flyt), intent(in) :: threshold
    !> the coupling parameter thingy
    real(flyt), intent(in) :: lambda
    !
    integer :: i
    !
    viol=0.0_flyt
    do i=1,size(alpha,1)
        if ( alpha(i) .gt. threshold ) then
            viol(i)=abs(lambda+gradient(i))
        elseif ( alpha(i) .lt. -threshold ) then
            viol(i)=abs(lambda-gradient(i))
        else
            viol(i)=max( max( -gradient(i)-lambda, -lambda+gradient(i) ), 0.0_flyt )
        endif
    enddo
end subroutine

!> LASSO L1 norm approximation thingy. Version 1.
subroutine lo_lasso_solve_GaussSiedel(X, y, lambda, solution)
    !> the coefficient matrix
    real(flyt), dimension(:,:), intent(in) :: X
    !> the values
    real(flyt), dimension(:), intent(in) :: y
    !> the coupling coefficient thing
    real(flyt), intent(in) :: lambda
    !> solution
    real(flyt), dimension(:), intent(out) :: solution
    !
    real(flyt), dimension(:), allocatable :: alpha
    real(flyt), dimension(:), allocatable :: Xy,gradient,viol
    real(flyt), dimension(:,:), allocatable :: XX 
    real(flyt) :: f0,tau,t0
    real(flyt) :: slope_maxviol,threshold,alpha_original
    integer :: iteration,line_mins
    integer :: i,j
    integer :: nf,nu,i_maxviol
    !
    ! I think I want to solve
    !
    ! || X*w -y ||^2 + lambda*|| w ||^1
    !
    ! I stole the implementation from some dude on the internet.
    ! There should probably be some kind of preconditioning step here,
    ! to make sure the tolerances and stuff works
    !
    threshold=1E-15_flyt
    nf=size(X,1)
    nu=size(X,2)
    ! Ridge regression starting guess for alpha, the solution
    lo_allocate(alpha(nu))
    call lo_ridgeregression(X,y,lambda,alpha)
    
    ! Xy = X^T * y
    lo_allocate(Xy(nu))
    Xy=0.0_flyt
    call lo_dgemv(X,y,Xy,trans='T')
    ! XX = X^T * X
    lo_allocate(XX(nu,nu))
    XX=0.0_flyt
    call lo_dgemm(X,X,XX,transa='T',transb='N',alpha=1.0_flyt,beta=0.0_flyt)
    ! gradient = XX*alpha - Xy
    lo_allocate(gradient(nu))
    gradient=Xy
    call lo_dgemv(XX,alpha,gradient,trans='T',alpha=1.0_flyt,beta=-1.0_flyt)    
    ! Compute violation?
    lo_allocate(viol(nu))
    call compute_violation(alpha,gradient,viol,threshold,lambda)
    ! Find worst violator
    f0=-lo_huge
    i_maxviol=0
    do i=1,nu
        if ( abs(viol(i)) .gt. f0 ) then
            f0=abs(viol(i))
            i_maxviol=i
        endif
    enddo
    
    ! Start minimizing
    iteration=0
    line_mins=0
    t0=walltime()
    tau=1E-14_flyt ! some kind of tolerance
    outerloop: do iteration=1,2*nu
        !
        innerloop: do i=1,10*nu
            ! get some kind of slope
            slope_maxviol=varSlope(alpha,i_maxviol,lambda*0.5_flyt,gradient,threshold)
            alpha_original=alpha(i_maxviol)
            ! optimize this alpha
            alpha(i_maxviol) = alpha(i_maxviol) - slope_maxviol/XX(i_maxviol,i_maxviol)
            if ( abs(alpha_original) .gt. lo_tiny ) then
            if ( sgn(alpha_original) .ne. sgn(alpha(i_maxviol)) ) then
                alpha(i_maxviol)=0.0_flyt
            endif
            endif
            ! get a new gradient and max violator
            gradient=Xy
            call lo_dgemv(XX,alpha,gradient,trans='T',alpha=1.0_flyt,beta=-1.0_flyt)
            call compute_violation(alpha,gradient,viol,threshold,lambda*0.5_flyt)
            f0=-lo_huge
            do j=1,nu
                if ( abs(alpha(j)) .gt. threshold ) then
                if ( viol(j) .gt. f0 ) then
                    f0=viol(j)
                    i_maxviol=j
                endif
                endif
            enddo
            ! break if no more violators
            if ( viol(i_maxviol) .lt. tau ) then
                write(*,*) 'converged inner loop'
                exit innerloop 
            endif
        enddo innerloop
        
        ! new new kind of violator
        f0=-lo_huge
        i_maxviol=0
        do j=1,nu
            if ( abs(alpha(j)) .lt. threshold ) then
            if ( viol(j) .gt. f0 ) then
                f0=viol(j)
                i_maxviol=j
            endif
            endif
        enddo
        !
        if ( abs(f0) .lt. tau ) then
            write(*,*) 'converged LASSO in',walltime()-t0,'seconds'
            exit outerloop
        endif
        !
        !opf="("//trim(int2char(nu))//"(2X,E17.9))"
        !write(*,opf) alpha
        write(*,*) 'iter',iteration,'tau',abs(f0),minval(abs(alpha)),maxval(abs(alpha))
        !
    enddo outerloop

    ! Return the solution
    solution=alpha
    
    contains
    !> signum function. Apparently it does not exist in fortran by default.
    integer function sgn(x)
        real(flyt), intent(in) :: x
        sgn=int(anint(sign(1.0_flyt,x)))
    end function
    !> have to figure out what this really does.
    function varSlope(alpha,i,lambda,gradient,threshold) result(slope)
        real(flyt), dimension(:), intent(in) :: alpha,gradient
        real(flyt), intent(in) :: lambda,threshold
        integer, intent(in) :: i
        real(flyt) :: slope
        !
        if ( alpha(i) .gt. lo_tiny .or. & 
        ( abs(alpha(i)) .gt. threshold .and. lambda+gradient(i) .lt. 0.0_flyt ) ) then
            slope=lambda+gradient(i)
        elseif ( alpha(i) .lt. -lo_tiny .or. ( abs(alpha(i)) .gt. threshold .and. -lambda+gradient(i) .gt. 0.0_flyt ) ) then
            slope=-lambda+gradient(i)
        else
            slope=0.0_flyt
        endif
    end function
end subroutine

end module
