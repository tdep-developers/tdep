#include "precompilerdefinitions"
module type_polynomial_interpolation
!! Several kinds of polynmial interpolations/fittings
use konstanter, only : flyt,lo_huge,lo_hugeint,lo_tol,lo_sqtol,lo_exitcode_param
use gottochblandat, only: tochar,lo_stop_gracefully
use type_blas_lapack_wrappers, only: lo_dgels
implicit none

private
public :: lo_polynomial
public :: lo_grid_interpolation

!> settings for a polynomial fit
type lo_polynomial
    !> how many dimensions
    integer :: ndim=-lo_hugeint
    !> order of polynomial
    integer :: order=-lo_hugeint
    !> how many coefficients
    integer :: ncoeff=-lo_hugeint
    !> polynomical exponents
    integer, dimension(:,:), allocatable :: exponents
    !> coefficient matrix
    real(flyt), dimension(:,:), allocatable :: coeffM
    !> informal names of the variables
    character(len=20), dimension(:), allocatable :: variablename
    contains
        !> set up the polynomial fitting thing
        procedure :: init=>initialize_polynomial
        !> evalute the polynomial at some point
        procedure :: eval=>evaluate_polynomial
#ifdef AGRESSIVE_SANITY
        !> storage size
        procedure :: size_in_mem=>memsize_polynomial
#endif
end type

!> helper type that creates grid-interpolations
type lo_grid_interpolation
    !> how many dimensions
    integer :: ndim
    !> how many points
    integer :: npts
    !> how many variables does it hold the interpolation for?
    integer :: nvar
    !> grid coordinates
    real(flyt), dimension(:,:), allocatable :: gridcoord
    !> grid values
    real(flyt), dimension(:,:), allocatable :: gridvals
    !> shift of coordinate
    real(flyt), dimension(:), allocatable :: coord_shift
    !> scaling of coordinate
    real(flyt), dimension(:), allocatable :: coord_scale
    !> characteristic distance
    real(flyt) :: coord_distance
    !> work matrices
    real(flyt), dimension(:,:), allocatable :: wM,wV
    real(flyt), dimension(:), allocatable :: weights
    !> polynomial to use
    type(lo_polynomial) :: pl
    contains
        !> create interpolation
        procedure :: generate=>create_grid_interpolation
        !> evaluate interpolation
        procedure :: eval=>evaluate_interpolation
#ifdef AGRESSIVE_SANITY
        !> storage size
        procedure :: size_in_mem=>memsize_interpolation
#endif
end type

contains

! Some sanity things
#ifdef AGRESSIVE_SANITY
#include "type_polynomial_interpolation_memory.f90"
#endif
!> set up a polynomial for fitting
subroutine initialize_polynomial( pl,order,ndim,xivals,xinames,order_per_dim )
    !> N-D polynomical
    class(lo_polynomial), intent(out) :: pl
    !> order of polynomial
    integer, intent(in) :: order
    !> number of dimensions
    integer, intent(in) :: ndim
    !> x1,x2,x3 ... values for the fit
    real(flyt), dimension(:,:), intent(in) :: xivals
    !> informal names of the variables. Helps debugging a lot.
    character(len=*), dimension(:), intent(in) :: xinames
    !> max order per dimensions
    integer, dimension(ndim), intent(in), optional :: order_per_dim

    ! Some default parameters I have to fix in the near future:
    integer, parameter :: maxdim=7      ! 6 structural parameters, 1 temperature. Easy to fix more.
    integer, parameter :: maxorder=4    ! reasonable hard limit for now
    integer, parameter :: maxcoeff=6435 ! enough for 7 dimensions to 8th order
    
    real(flyt) :: f0
    integer, dimension(:,:), allocatable :: di
    integer, dimension(maxdim) :: maxexp
    integer :: i1,i2,i3,i4,i5,i6,i7
    integer :: i,j,k,l,nxi

    ! check that the input is fine
    if ( ndim .gt. maxdim ) then
        call lo_stop_gracefully(['Current max number of dimensions is ',tochar(maxdim),'. Can fix more if you nag.'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( order .gt. maxorder ) then
        call lo_stop_gracefully(['Current max order of fit is ',tochar(maxorder),', because I say so.'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    pl%order=order
    pl%ndim=ndim
    maxexp=pl%order
    if ( present(order_per_dim) ) then
        maxexp(1:ndim)=order_per_dim
    endif

    lo_allocate(di(maxdim,maxcoeff))
    di=0
    select case(pl%ndim)
        case(1)
            l=0
            di=-1
            do i1=0,maxexp(1)
                if ( i1 .le. pl%order ) then
                    l=l+1
                    di(1,l)=i1
                endif
            enddo
        case(2)
            l=0
            di=-1
            do i1=0,maxexp(1) 
            do i2=0,maxexp(2) 
                if ( i1+i2 .le. pl%order ) then
                    l=l+1
                    di(1:2,l)=[i1,i2]
                endif
            enddo
            enddo
        case(3)
            l=0
            di=-1
            do i1=0,maxexp(1) 
            do i2=0,maxexp(2) 
            do i3=0,maxexp(3) 
                if ( i1+i2+i3 .le. pl%order ) then
                    l=l+1
                    di(1:3,l)=[i1,i2,i3]
                endif
            enddo
            enddo
            enddo
        case(4)
            l=0
            di=-1
            do i1=0,maxexp(1) 
            do i2=0,maxexp(2) 
            do i3=0,maxexp(3) 
            do i4=0,maxexp(4) 
                if ( i1+i2+i3+i4 .le. pl%order ) then
                    l=l+1
                    di(1:4,l)=[i1,i2,i3,i4]
                endif
            enddo
            enddo
            enddo
            enddo
        case(5)
            l=0
            di=-1
            do i1=0,maxexp(1) 
            do i2=0,maxexp(2) 
            do i3=0,maxexp(3) 
            do i4=0,maxexp(4) 
            do i5=0,maxexp(5) 
                if ( i1+i2+i3+i4+i5 .le. pl%order ) then
                    l=l+1
                    di(1:5,l)=[i1,i2,i3,i4,i5]
                endif
            enddo
            enddo
            enddo
            enddo
            enddo
        case(6)
            l=0
            di=-1
            do i1=0,maxexp(1) 
            do i2=0,maxexp(2) 
            do i3=0,maxexp(3) 
            do i4=0,maxexp(4) 
            do i5=0,maxexp(5) 
            do i6=0,maxexp(6) 
                if ( i1+i2+i3+i4+i5+i6 .le. pl%order ) then
                    l=l+1
                    di(1:6,l)=[i1,i2,i3,i4,i5,i6]
                endif
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
        case(7)
            l=0
            di=-1
            do i1=0,maxexp(1) 
            do i2=0,maxexp(2) 
            do i3=0,maxexp(3) 
            do i4=0,maxexp(4) 
            do i5=0,maxexp(5) 
            do i6=0,maxexp(6) 
            do i7=0,maxexp(7) 
                if ( i1+i2+i3+i4+i5+i6+i7 .le. pl%order ) then
                    l=l+1
                    di(1:7,l)=[i1,i2,i3,i4,i5,i6,i7]
                endif
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
        case default
            call lo_stop_gracefully(['Not done with ',tochar(pl%ndim),'-dimensional fits.'],lo_exitcode_param,__FILE__,__LINE__)
    end select

    ! Ok that might have been unnecessarily general, but whatever. Store the number of coefficients?
    pl%ncoeff=l
    lo_allocate(pl%exponents(pl%ndim,pl%ncoeff))
    pl%exponents=di(1:pl%ndim,1:pl%ncoeff)

    ! Store the x1,x2,x3,... datapoints
    nxi=size(xivals,2)
    ! Check right away that I have enough data-points for a fit
    if ( nxi .lt. pl%ncoeff+2 ) then
        call lo_stop_gracefully(['You have '//tochar(nxi)//' data points for '//tochar(pl%ncoeff)//' polynomical coefficients. I have decided that is too little.'],lo_exitcode_param)
    endif

    ! Store the names of the variables
    lo_allocate(pl%variablename(pl%ndim))
    do i=1,ndim
        pl%variablename(i)=trim(adjustl(xinames(i)))
    enddo

    ! Get the coefficient matrix
    lo_allocate(pl%coeffM(nxi,pl%ncoeff))
    pl%coeffM=0.0_flyt
    do i=1,nxi
        do j=1,pl%ncoeff
            f0=1.0_flyt
            do k=1,pl%ndim
                f0=f0*xivals(k,i)**(pl%exponents(k,j))
            enddo
            pl%coeffM(i,j)=f0
        enddo
    enddo

end subroutine

!> evaluate the polynomial. Should not have to be fast, only robust.
function evaluate_polynomial(pl,depvar,coefficients) result(val)
    !> grid simulation thingy
    class(lo_polynomial), intent(in) :: pl
    !> input variables
    real(flyt), dimension(:), intent(in) :: depvar
    !> polynomial coefficients
    real(flyt), dimension(:), intent(in) :: coefficients
    !> resulting value
    real(flyt) :: val

    real(flyt) :: f0
    integer :: i,j

    ! Some sanity checks are always neat.
    if ( size(coefficients,1) .ne. pl%ncoeff ) then
        call lo_stop_gracefully(['Trying to evaluare polynomial with ',tochar(pl%ncoeff),' coefficients, but ',tochar(size(coefficients,1)),' was provided.'],&
            lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( size(depvar,1) .ne. pl%ndim ) then
        call lo_stop_gracefully(['Trying to evaluare polynomial with ',tochar(pl%ndim),' dimensions, but ',tochar(size(depvar,1)),' was provided.'],&
            lo_exitcode_param,__FILE__,__LINE__)
    endif

    val=0.0_flyt
    do i=1,pl%ncoeff
        f0=coefficients(i)
        do j=1,pl%ndim
            f0=f0*depvar(j)**pl%exponents(j,i)
        enddo
        val=val+f0
    enddo
end function

!> create a grid interpolation
subroutine create_grid_interpolation(ip,gridcoord,gridval,order,coordscale,maxorder)
    !> interpolator
    class(lo_grid_interpolation), intent(out) :: ip
    !> grid coordinates
    real(flyt), dimension(:,:), intent(in) :: gridcoord
    !> values on the grid
    real(flyt), dimension(:,:), intent(in) :: gridval
    !> order of interpolation
    integer, intent(in) :: order
    !> characteristic length scale
    real(flyt), intent(in) :: coordscale
    !> different orders in different directions?
    integer, dimension(:), intent(in) :: maxorder

    character(len=2), dimension(size(gridcoord,1)) :: vname
    integer :: i

    ip%ndim=size(gridcoord,1)
    ip%npts=size(gridcoord,2)
    ip%nvar=size(gridval,2)
    lo_allocate(ip%gridcoord(ip%ndim,ip%npts))
    lo_allocate(ip%gridvals(ip%npts,ip%nvar))
    lo_allocate(ip%coord_shift(ip%ndim))
    lo_allocate(ip%coord_scale(ip%ndim))

    ! get the coordinate transformation
    do i=1,ip%ndim
        ip%coord_shift(i)=minval(gridcoord(i,:))
        ip%coord_scale(i)=maxval(gridcoord(i,:))-minval(gridcoord(i,:))
        if ( ip%coord_scale(i) .gt. lo_tol ) then
            ip%coord_scale(i)=1.0_flyt/ip%coord_scale(i)
        else
            ip%coord_scale(i)=0.0_flyt
        endif
    enddo

    ! Store grid values
    do i=1,ip%npts
        ip%gridcoord(:,i)=(gridcoord(:,i)-ip%coord_shift)*ip%coord_scale
        ip%gridvals(i,:)=gridval(i,:)
    enddo
    ! My characteristic distance is 0.1 in these scaled coordinates. Perhaps should not be hard-coded.
    ip%coord_distance=coordscale
    ! Create the polynomial
    do i=1,ip%ndim
        vname(i)='x'//tochar(i)
    enddo
    call ip%pl%init(order,ip%ndim,ip%gridcoord,vname,maxorder)
    ! Space for coefficient matrix and solution
    lo_allocate(ip%wM( ip%npts,ip%pl%ncoeff ))
    lo_allocate(ip%wV( ip%npts,1 ))
    lo_allocate(ip%weights( ip%npts ))
    ip%wM=0.0_flyt
    ip%wV=0.0_flyt
    ip%weights=0.0_flyt
end subroutine

!> evaluate the interpolation
function evaluate_interpolation(ip,var,gridcoord,rexp) result(val)
    !> the interpolation
    class(lo_grid_interpolation), intent(inout) :: ip
    !> which variable
    integer, intent(in) :: var
    !> coordinates to interpolate
    real(flyt), dimension(ip%ndim), intent(in) :: gridcoord
    !> what exponent to use when evaluating
    integer, intent(in), optional :: rexp
    !> value
    real(flyt) :: val

    real(flyt), dimension(ip%ndim) :: scaledcoord
    real(flyt) :: f0
    integer :: i,ee

    ! what exponent on the inverse distance?
    if ( present(rexp) ) then
        ee=rexp
    else
        ee=2
    endif

    scaledcoord=(gridcoord-ip%coord_shift)*ip%coord_scale
    ! Get the weights
    do i=1,ip%npts
        f0=1.0_flyt/( norm2(scaledcoord-ip%gridcoord(:,i))**ee+ip%coord_distance )
        ip%wM(i,:)=ip%pl%coeffM(i,:)*f0
        ip%wV(i,1)=ip%gridvals(i,var)*f0
    enddo
    ! Solve for the coefficients
    call lo_dgels(ip%wM,ip%wV)
    ! Evaluate
    val=ip%pl%eval(scaledcoord,ip%wV(1:ip%pl%ncoeff,1))
end function

end module
