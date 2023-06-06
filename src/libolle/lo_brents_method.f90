module lo_brents_method
!!
!! Reverse communication interface to Brent's method of finding the zero of a function.
!!
!! It's real easy to use:
!!
!! ! pick x0 and x1 so that they safely bracket 0
!! y0=function(x0)
!! y1=function(x1)
!! ! initialize the solver
!! call bt%init( x0,x1,y0,y1,tol )
!! ! iterate to find zero
!! do iter=1,maxiter
!!     ! check for convergence
!!     if ( abs(bt%current_y()) .lt. bt%tolerance ) then
!!         ! found a zero
!!         x=bt%current_x()
!!         exit
!!     endif
!!     ! pick next x-value
!!     x=bt%next_x()
!!     ! evaluate function
!!     y=function(x)
!!     ! tell the solver what happened
!!     call bt%update(x,y)
!! enddo
!!
use konstanter, only: r8,lo_hugeint,lo_huge

implicit none
private
public :: lo_brent_helper

!> helper for Brent's method of finding zeros
type lo_brent_helper
    !> x-value, lower bound
    real(r8), private :: xlo=-lo_huge
    !> x-value, upper bound
    real(r8), private :: xhi=-lo_huge
    !> y-value, lower bound
    real(r8), private :: ylo=-lo_huge
    !> y-value, upper bound
    real(r8), private :: yhi=-lo_huge
    !> tolerance
    real(r8) :: tolerance=-lo_huge
    !> iteration counter
    integer, private :: niter=-lo_hugeint
    !> helper variables
    logical, private :: mflag=.false.

    real(r8), private :: xc=-lo_huge
    real(r8), private :: fc=-lo_huge
    real(r8), private :: fd=-lo_huge
    real(r8), private :: s =-lo_huge
    real(r8), private :: xd=-lo_huge
    real(r8), private :: fs=-lo_huge
    real(r8), private :: sd=-lo_huge
    contains
        !> initialize the helper
        procedure :: init
        !> get the current function value/argument
        procedure :: current_y
        procedure :: current_x
        !> next argument to evaluate
        procedure :: next_x
        !> update solver with results
        procedure :: update
end type

contains

!> initialize Brent's method
subroutine init(bt,xlo,xhi,ylo,yhi,tolerance)
    !> Brent's method handle
    class(lo_brent_helper), intent(out) :: bt
    !> safe lower bound
    real(r8), intent(in) :: xlo
    !> safe upper bound
    real(r8), intent(in) :: xhi
    !> safe lower bound
    real(r8), intent(in) :: ylo
    !> safe upper bound
    real(r8), intent(in) :: yhi
    !> tolerance when to quit
    real(r8), intent(in) :: tolerance

    ! Get them in the right order
    if ( abs(ylo) .lt. abs(yhi) ) then
        bt%xlo=xlo
        bt%xhi=xhi
        bt%ylo=ylo
        bt%yhi=yhi
    else
        bt%xlo=xhi
        bt%xhi=xlo
        bt%ylo=yhi
        bt%yhi=ylo
    endif

    ! Set the helper variables
    bt%tolerance=tolerance
    bt%niter=0
    bt%xc=bt%xlo
    bt%fc=bt%ylo
    bt%mflag=.true.
    bt%s=0.0_r8
    bt%xd=0.0_r8
    bt%fs=-lo_huge
    bt%fd=-lo_huge
end subroutine

!> fetch the next argument to evaluate function at
function next_x(bt) result(x)
    !> Brent's method handle
    class(lo_brent_helper), intent(inout) :: bt
    !> next x-value to evaluate the function at
    real(r8) :: x

    real(r8) :: f0,f1
    logical :: sflag

    ! Figure out where to put the next point
    if ( abs(bt%ylo-bt%fc) .gt. bt%tolerance .and. abs(bt%yhi-bt%fc) .gt. bt%tolerance ) then
        x=bt%xlo*bt%yhi*bt%fc/( (bt%ylo-bt%yhi)*(bt%ylo-bt%fc) ) + &
        bt%xhi*bt%ylo*bt%fc/( (bt%yhi-bt%ylo)*(bt%yhi-bt%fc) )+&
        bt%xc*bt%ylo*bt%yhi/( (bt%fc-bt%ylo)*(bt%fc-bt%yhi) )
    else
        x=bt%xhi-bt%yhi*(bt%xhi-bt%xlo)/(bt%yhi-bt%ylo)
    endif

    ! Also try to figure out is that was a bad idea, and revert to just slicing the interval in half.
    sflag=.false.
    f0=min((3*bt%xlo+bt%xhi)*0.25_r8,bt%xhi)
    f1=max((3*bt%xlo+bt%xhi)*0.25_r8,bt%xhi)
    if ( x .lt. f0      .or.  x                 .gt. f1                       ) sflag=.true.
    if ( bt%mflag       .and. abs(x-bt%xhi)     .ge. abs(bt%xhi-bt%xc)*0.5_r8 ) sflag=.true.
    if ( .not.bt%mflag  .and. abs(x-bt%xhi)     .ge. abs(bt%xc-bt%xd)*0.5_r8  ) sflag=.true.
    if ( bt%mflag       .and. abs(bt%xhi-bt%xc) .lt. bt%tolerance             ) sflag=.true.
    if ( .not.bt%mflag  .and. abs(bt%xc-bt%xd)  .lt. bt%tolerance             ) sflag=.true.

    sflag=.true.
    if ( sflag ) then
        x=(bt%xlo+bt%xhi)*0.5_r8
        bt%mflag=.true.
    else
        bt%mflag=.false.
    endif
end function

!> update solver based on function value
subroutine update(bt,x,y)
    !> Brent's method handle
    class(lo_brent_helper), intent(inout) :: bt
    !> function argument
    real(r8), intent(in) :: x
    !> function value
    real(r8), intent(in) :: y

    real(r8) :: f0

    ! Update the solver settings
    bt%xd=bt%xc
    bt%fd=bt%fc
    bt%xc=bt%xhi
    bt%fc=bt%yhi
    if ( bt%ylo*y .lt. 0.0_r8 ) then
        bt%xhi=x
        bt%yhi=y
    else
        bt%xlo=x
        bt%ylo=y
    endif

    ! Maybe flip?
    if ( abs(bt%ylo) .lt. abs(bt%ylo) ) then
        f0=bt%ylo
        bt%ylo=bt%yhi
        bt%yhi=f0
        f0=bt%xlo
        bt%xlo=bt%xhi
        bt%xhi=f0
    endif

    ! Count iterations?
    bt%niter=bt%niter+1
end subroutine

!> fetch current argument
function current_x(bt) result(x)
    !> Brent's method handle
    class(lo_brent_helper), intent(in) :: bt
    !> function argument
    real(r8) :: x

    ! check for convergence
    if ( abs(bt%xhi) .lt. abs(bt%fs) ) then
        x=bt%xhi
    else
        x=bt%s
    endif
end function

!> fetch current function value
function current_y(bt) result(y)
    !> Brent's method handle
    class(lo_brent_helper), intent(in) :: bt
    !> function value
    real(r8) :: y

    ! check for convergence
    if ( abs(bt%xhi) .lt. abs(bt%fs) ) then
        y=bt%yhi
    else
        y=bt%fs
    endif
end function

end module
