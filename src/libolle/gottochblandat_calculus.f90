#include "precompilerdefinitions"
submodule (gottochblandat) gottochblandat_calculus
implicit none
contains

!> Given a function y=f(x), you want it defined as y'=f(x'), where x' is some other array of values than x (but should cover the whole thing). Assumes that both x and xi are linearly spaced.
module subroutine lo_put_function_on_new_axis(x,y,xi,yi,sigma,preservenorm)
    !> original x-axis
    real(r8), dimension(:), intent(in) :: x
    !> original y-axis
    real(r8), dimension(:), intent(in) :: y
    !> new x-axis
    real(r8), dimension(:), intent(in) :: xi
    !> new y-axis
    real(r8), dimension(:), intent(out) :: yi
    !> smearing parameter
    real(r8), intent(in), optional :: sigma
    !> make sure the norm is preserved
    logical, intent(in), optional :: preservenorm

    integer :: i,j,n,ni,ii,jj
    real(r8) :: f0,sig,xmin,xmax,inv_xrange,fivesigma

    ! First some sanity checks:
    if ( minval(xi) .gt. maxval(x) .or. maxval(xi) .lt. minval(x) ) then
        write(*,*) 'No overlap between x and xi, something is strange'
        stop
    endif

    ! Figure out a sigma:
    if ( present(sigma) ) then
        sig=sigma
    else
        ! the smearing is guessed from the largest spacing, in
        ! either the old or the new array.
        sig=max(x(2)-x(1),xi(2)-xi(1))*1.0_r8
    endif

    ! sizes of the arrays
    n=size(x,1)
    ni=size(xi,1)

    ! Some helper number to make the convolution O(N)
    xmin=minval(x)
    xmax=maxval(x)
    inv_xrange=1.0_r8/(xmax-xmin)
    fivesigma=5.0_r8*sig

    ! Do the actual convolution
    yi=0.0_r8
    do i=1,ni
        ! get the convolution loop indices
        ii=floor((xi(i)-xmin-fivesigma)*inv_xrange*n)
        ii=min(max(1,ii),n)
        jj=ceiling((xi(i)-xmin+fivesigma)*inv_xrange*n)
        jj=min(max(1,jj),n)
        ! do the convolution
        do j=ii,jj
            yi(i)=yi(i)+y(j)*lo_gauss(xi(i),x(j),sig)
        enddo
    enddo
    ! make sure it normalizes ok
    yi=yi*(x(2)-x(1))

    ! Double-check that it integrates to the same thing, sort of, maybe
    if ( present(preservenorm) ) then
        if ( preservenorm ) then
            f0=lo_trapezoid_integration(x,y)
            yi=yi*f0/lo_trapezoid_integration(xi,yi)
        endif
    endif
end subroutine

!> Linear interpolation of data. Interpolate y=f(x) to point xi. Assumes that x,y are ordered properly, increasing monotonically. If xi is outside the interval of x, I return 0. Could be made faster if I would assume that x is linearly spaced, but this is not used for anything performance sensitive anyway.
module pure function lo_linear_interpolation(x,y,xi,threshold) result(yi)
    !> function x-values
    real(r8), dimension(:), intent(in) :: x
    !> function y-values
    real(r8), dimension(:), intent(in) :: y
    !> point x to interpolate to
    real(r8), intent(in) :: xi
    !> interpolated value
    real(r8) :: yi
    !> threshold to check if xi is inside the interval
    real(r8), intent(in), optional :: threshold
    !
    integer :: i,n,ii,jj
    real(r8) :: thres,A
    !
    n=size(x,1)
    ii=0
    jj=0

    ! simple edge case
    if ( n .eq. 2 ) then
        A=(xi-x(1))/(x(2)-x(1))
        if ( A .le. 0.0_r8 ) then
            yi=0.0_r8
            return
        elseif ( A .ge. 1.0_r8 ) then
            yi=0.0_r8
            return
        else
            yi=y(1)+(1.0_r8-A)*y(2)
            return
        endif
    endif

    ! Alter the threshold, perhaps
    if ( present(threshold) ) then
        thres=threshold
    else
        thres=lo_sqtol
    endif

    ! Check some simple thresholds, are we exactly at the endpoints?
    if ( abs(xi-x(1)) .lt. thres ) then
        yi=y(1)
        return
    endif
    if ( abs(xi-x(n)) .lt. thres ) then
        yi=y(n)
        return
    endif

    ! If I am outside the interval, return 0
    if ( minval(x)-xi .gt. 2*thres .or. abs(maxval(x))-xi .lt. -2*thres ) then
        yi=0.0_r8
        return
    endif
    ! Find the relevant interval
    ii=n
    do i=2,n
        if ( x(i)-xi .gt. thres ) then
            ii=i
            jj=i-1
            exit
        endif
    enddo
    ! check the edge cases
    if ( ii .lt. 2 ) then
        ii=2
        jj=1
    endif

    if ( ii .ge. n ) then
        ii=n
        jj=n-1
    endif
    ! Do the interpolation
    A=(x(ii)-xi)/(x(ii)-x(jj))
    yi=A*y(jj)+(1.0_r8-A)*y(ii)
end function

!> distribute n points kind of uniform on the unit sphere
module pure subroutine lo_points_on_sphere(r,randomseed)
    !> the points
    real(r8), dimension(:,:), intent(out) :: r
    !> a random 0-1 seed to shift the points around
    real(r8), intent(in), optional :: randomseed

    real(r8) :: offset,increment,x,y,z,rad,phi,rnd
    integer :: i,n

    if ( present(randomseed) ) then
        rnd=randomseed
    else
        rnd=0.0_r8
    endif

    n=size(r,2)
    r=0.0_r8
    offset=2.0_r8/(1.0_r8*n)
    increment=lo_pi*(3.0_r8-sqrt(5.0_r8))
    do i=1,n
        ! Negative 2 here compared to what you might find, since
        ! Fortran is 1-indexed.
        y=((i-1)*offset-1)+offset*0.5_r8;
        rad=sqrt( max(1.0_r8-y**2,0.0_r8) );
        phi=mod(real(i+rnd,r8),real(n,r8))*increment;
        x=cos(phi)*rad;
        z=sin(phi)*rad;
        r(:,i)=[x,y,z]/norm2([x,y,z]);
    enddo
end subroutine

!> Return linearly spaced points in an array. Assumes the array is allocated to the correct length.
module pure subroutine lo_linspace(minv,maxv,x)
    !> linearly spaced values
    real(r8), dimension(:), intent(out) :: x
    !> minimum value
    real(r8), intent(in) :: minv
    !> maximum value
    real(r8), intent(in) :: maxv
    !
    integer :: i,n
    real(r8) :: f0,f1
    !
    n=size(x,1)
    ! sanity check
    if ( n .eq. 1 ) then
        ! I just want one number, I guess
        x(1)=(maxv+minv)*0.5_r8
    else
        ! More than one number
        f1=maxv-minv
        do i=1,n
            f0=(i*1.0_r8-1.0_r8)/(n*1.0_r8-1.0_r8)
            x(i)=minv+f0*f1
        enddo
        x(1)=minv
        x(n)=maxv
    endif
end subroutine

!> Return logarithmically spaced points in an array. Assumes the array is allocated to the correct length.
module pure subroutine lo_logspace(minv,maxv,x)
    !> logspaced values returned in x
    real(r8), dimension(:), intent(out) :: x
    !> minimum value
    real(r8), intent(in) :: minv
    !> maximum value
    real(r8), intent(in) :: maxv

    real(r8), dimension(:), allocatable :: dum
    real(r8) :: shift,f0,f1
    integer :: i,n

    n=size(x,1)
    ! if just one number, just get the mean
    if ( n .eq. 1 ) then
        x(1)=(maxv+minv)*0.5_r8
    else
        allocate(dum(n))
        ! Get a shift, so that I don't take log(0)
        if ( minv .lt. lo_tol ) then
            shift=1.0_r8-minv
        else
            shift=0.0_r8
        endif
        f0=log(minv+shift)
        f1=log(maxv+shift)
        call lo_linspace(f0,f1,dum)
        do i=1,n
            x(i)=exp(dum(i))-shift
        enddo
        deallocate(dum)
    endif
end subroutine

!> Trapezoid integration of y=f(x)
module pure function lo_trapezoid_integration(x,y) result(f)
    !> x-values
    real(r8), dimension(:), intent(in) :: x
    !> function values
    real(r8), dimension(:), intent(in) :: y
    !> integral
    real(r8) :: f

    integer :: i
    f=0.0_r8
    do i=1,size(x,1)-1
        f=f+( y(i)+y(i+1) )*0.5_r8*( x(i+1)-x(i) )
    enddo
end function

!> Gaussian distribution
module elemental function lo_gauss(x,mu,sigma) result(g)
    !> point to evaluate
    real(r8), intent(in) :: x
    !> mean
    real(r8), intent(in) :: mu
    !> standard deviation
    real(r8), intent(in) :: sigma
    !> value
    real(r8) :: g

    ! 1/sqrt(2*pi)
    real(r8), parameter :: sqrt2piinv=0.398942280401433_r8
    real(r8) :: invsig
    if ( abs(x-mu) .gt. 8*sigma ) then
        ! after 8 sigma everything is basically zero anyway.
        g=0.0_r8
    else
        invsig=1.0_r8/sigma
        g=exp(-(x-mu)*(x-mu)*0.5_r8*invsig*invsig)*sqrt2piinv*invsig
    endif
end function

!> Lorentzian distribution
module elemental function lo_lorentz(x,mu,sigma) result(l)
    !> point to evaluate
    real(r8), intent(in) :: x
    !> mean
    real(r8), intent(in) :: mu
    !> Full width half maximum
    real(r8), intent(in) :: sigma
    !> value
    real(r8) :: l
    !
    real(r8), parameter :: invpi=0.318309886183791_r8
    real(r8) :: s2
    s2=sigma*0.5_r8
    l=(s2*invpi)/( (x-mu)**2 + s2**2)
end function

! !> Damped harmonic oscillator distribution
! module pure function lo_dampedoscillator(x,mu,sigma) result(d)
!     !> point to evaluate
!     real(r8), intent(in) :: x
!     !> mean
!     real(r8), intent(in) :: mu
!     !> width
!     real(r8), intent(in) :: sigma
!     !> value
!     real(r8) :: d
!     !
!     real(r8) :: f0,f1,f2
!     ! Don't remember why I add a small number here.
!     f0=(sigma/x)**2 + (x/mu-mu/x)**2 + 1E-10_r8
!     f1=sigma/(mu*x*f0)
!     f2=(lo_pi+2.0_r8*atan(mu/sigma))*0.25_r8
!     d=f1/f2
! end function

!> Mean of an array
module pure function lo_mean(x) result(m)
    !> values
    real(r8), dimension(:), intent(in) :: x
    !> mean
    real(r8) :: m
    m=sum(x)/(size(x,1)*1.0_r8)
end function

!> Standard deviation
module pure function lo_stddev_1d(x) result(s)
    !> values
    real(r8), dimension(:), intent(in) :: x
    !> standard deviation
    real(r8) :: s

    real(r8) :: mean
    integer :: i

    mean=sum(x)/real(size(x),r8)
    s=0.0_r8
    do i=1,size(x)
        s=s+(mean-x(i))**2
    enddo
    s=sqrt(s/real(size(x),r8))
end function
module pure function lo_stddev_3d(x) result(s)
    !> values
    real(r8), dimension(:,:,:), intent(in) :: x
    !> standard deviation
    real(r8) :: s

    real(r8) :: mean
    integer :: i,j,k

    mean=sum(x)/real(size(x),r8)
    s=0.0_r8
    do k=1,size(x,3)
    do j=1,size(x,2)
    do i=1,size(x,1)
        s=s+(mean-x(i,j,k))**2
    enddo
    enddo
    enddo
    s=sqrt(s/real(size(x),r8))
end function

!> statistical Rsquare thing
module pure function lo_rsquare_1d(values,predictions) result(rsq)
    !> true values
    real(r8), dimension(:), intent(in) :: values
    !> predicted values
    real(r8), dimension(:), intent(in) :: predictions
    !> metric
    real(r8) :: rsq

    integer :: i
    real(r8) :: f0,f1,f2

    f0=sum(values)/real(size(values),r8)
    f1=0.0_r8
    f2=0.0_r8
    do i=1,size(values,1)
        f1=f1+(values(i)-f0)**2
        f2=f2+(values(i)-predictions(i))**2
    enddo
    rsq=1.0_r8-f2/f1
end function
module pure function lo_rsquare_3d(values,predictions) result(rsq)
    !> true values
    real(r8), dimension(:,:,:), intent(in) :: values
    !> predicted values
    real(r8), dimension(:,:,:), intent(in) :: predictions
    !> metric
    real(r8) :: rsq

    integer :: i,j,k
    real(r8) :: f0,f1,f2

    f0=sum(values)/real(size(values),r8)
    f1=0.0_r8
    f2=0.0_r8
    do k=1,size(values,3)
    do j=1,size(values,2)
    do i=1,size(values,1)
        f1=f1+(values(i,j,k)-f0)**2
        f2=f2+(values(i,j,k)-predictions(i,j,k))**2
    enddo
    enddo
    enddo
    rsq=1.0_r8-f2/f1
end function


!> Given four function values and four corners of a tetrahedron, return the gradient inside the tetrahedron
module pure function lo_linear_gradient_in_tetrahedron(corners,values,tolerance) result(gradient)
    !> four corners of tetrahedron
    real(r8), dimension(3,4), intent(in) :: corners
    !> function values at corners
    real(r8), dimension(4), intent(in) :: values
    !> tolerance for zero volume
    real(r8), intent(in) :: tolerance
    !> gradient
    real(r8), dimension(3) :: gradient

    real(r8) :: v1x,v1y,v1z
    real(r8) :: v2x,v2y,v2z
    real(r8) :: v3x,v3y,v3z
    real(r8) :: e1,e2,e3,f0

    ! Just set it up and solved in mathematica. Should be as fast as it gets.
    v1x=corners(1,2)-corners(1,1)
    v1y=corners(2,2)-corners(2,1)
    v1z=corners(3,2)-corners(3,1)
    v2x=corners(1,3)-corners(1,1)
    v2y=corners(2,3)-corners(2,1)
    v2z=corners(3,3)-corners(3,1)
    v3x=corners(1,4)-corners(1,1)
    v3y=corners(2,4)-corners(2,1)
    v3z=corners(3,4)-corners(3,1)
    e1=values(2)-values(1)
    e2=values(3)-values(1)
    e3=values(4)-values(1)

    ! Common prefactor
    f0=(v1z*v2y*v3x - v1y*v2z*v3x - v1z*v2x*v3y + v1x*v2z*v3y + v1y*v2x*v3z - v1x*v2y*v3z)
    if ( abs(f0) .gt. tolerance ) then
        f0=1.0_r8/f0
    else
        f0=0.0_r8
    endif

    ! Actual gradient
    gradient(1)= e3*v1z*v2y - e3*v1y*v2z - e2*v1z*v3y + e1*v2z*v3y + e2*v1y*v3z - e1*v2y*v3z
    gradient(2)=-e3*v1z*v2x + e3*v1x*v2z + e2*v1z*v3x - e1*v2z*v3x - e2*v1x*v3z + e1*v2x*v3z
    gradient(3)= e3*v1y*v2x - e3*v1x*v2y - e2*v1y*v3x + e1*v2y*v3x + e2*v1x*v3y - e1*v2x*v3y
    gradient=gradient*f0
end function

!!> Use some clever way to fit a distribution to my data. I do not remember where I stole it.
!subroutine lo_fit_distribution_to_data(xin,yin,distribution,mu,sig,tolerance)
!    !> x-values
!    real(r8), dimension(:), intent(in) :: xin
!    !> y-values
!    real(r8), dimension(:), intent(in) :: yin
!    !> what kind of distribution
!    character(len=*), intent(in) :: distribution
!    !> the mean
!    real(r8), intent(out) :: mu
!    !> standard deviation
!    real(r8), intent(out) :: sig
!    !> the tolerance when to stop
!    real(r8), intent(in), optional :: tolerance
!    !
!    integer :: i,n
!    integer :: iter
!    real(r8), dimension(:), allocatable :: x,y,y0
!    real(r8), dimension(:), allocatable :: dyds,dydm
!    real(r8), dimension(:,:), allocatable :: F,Ft
!    real(r8), dimension(2) :: a0,b,c
!    real(r8), dimension(2,2) :: a
!    real(r8) :: f0,f1
!    real(r8) :: xshift,xscale,yscale,tol
!
!    ! how many values?
!    n=size(xin,1)
!
!    ! Some sanity checks
!    if ( n .lt. 10 ) then
!        write(*,*) 'I need more than 10 values to fit the distribution'
!        stop
!    endif
!    !
!    if ( n .ne. size(yin,1) ) then
!        write(*,*) 'x and y needs to be of equal length when doing a distribution fit.'
!        stop
!    endif
!    !
!    if ( present(tolerance) ) then
!        tol=tolerance
!    else
!        tol=1E-15_r8
!    endif
!
!    ! Scale and normalize the values
!    lo_allocate(x(n))
!    lo_allocate(y(n))
!    ! get all x-values between 0 and 1
!    xshift=minval(xin)
!    xscale=maxval(xin)-xshift
!    x=(xin-xshift)/xscale
!    ! and it should normalize to 1
!    yscale=lo_trapezoid_integration(x,yin)
!    y=yin/yscale
!
!    ! Get some reasonable guesses for the standard deviation and the mean
!    f0=0.0_r8
!    f1=0.0_r8
!    do i=1,n
!        f0=f0+y(i)*x(i)
!        f1=f1+y(i)
!    enddo
!    mu=f0/f1
!
!    f0=0.0_r8
!    f1=0.0_r8
!    do i=1,n
!        f0=f0+y(i)*( x(i)-mu )**2
!        f1=f1+y(i)
!    enddo
!    sig=sqrt(f0/f1)
!
!    ! Seems like an ok estimate, now iteratively find the values.
!    lo_allocate(y0(n))
!    lo_allocate(dyds(n))
!    lo_allocate(dydm(n))
!    lo_allocate(F(n,2))
!    lo_allocate(Ft(2,n))
!    f0=mu
!    f1=sig
!    do iter=1,50
!        ! the current distribution, and the derivative with respect to sigma and mu
!        select case(distribution)
!        case('gaussian')
!            do i=1,n
!                y0(i)=lo_gauss(x(i),mu,sig)
!                dyds(i)=gauss_derivative_sig(x(i),mu,sig)*y0(i)
!                dydm(i)=gauss_derivative_mu(x(i),mu,sig)*y0(i)
!            enddo
!        case('lorentzian')
!            do i=1,n
!                y0(i)=lo_lorentz(x(i),mu,sig)
!                dyds(i)=lorentz_derivative_sig(x(i),mu,sig)*y0(i)*y0(i)
!                dydm(i)=lorentz_derivative_mu(x(i),mu,sig)*y0(i)*y0(i)
!            enddo
!        end select
!        ! ok, don't remember any of this. But get a matrix
!        F(:,1)=dyds
!        F(:,2)=dydm
!        ! and it's transpose, just to be safe
!        Ft(1,:)=dyds
!        Ft(2,:)=dydm
!        ! The starter values
!        a0=[sig,mu]
!        ! Fancy matrix, to be inverted
!        a=matmul(Ft,F)
!        call lo_invertmatrix_lapack(a)
!        ! get b
!        b=matmul(Ft,y-y0)
!        c=matmul(a,b)+a0
!        !
!        sig=c(1)
!        mu=c(2)
!        !
!        if ( abs(mu-f0) .lt. 1E-15_r8 ) then
!            exit
!        else
!            f0=mu
!            f1=sig
!        endif
!    enddo
!
!    ! Scale it back to the input units
!    mu=mu*xscale+xshift
!    sig=sig*xscale
!
!    ! Cleanup
!    lo_deallocate(x)
!    lo_deallocate(y)
!    lo_deallocate(y0)
!    lo_deallocate(F)
!    lo_deallocate(Ft)
!    lo_deallocate(dyds)
!    lo_deallocate(dydm)
!    contains
!    ! derivative of a gaussian with respect to sigma, divided by a gaussian
!    pure function gauss_derivative_sig(x,mu,sig) result(y)
!        real(r8), intent(in) :: x,mu,sig
!        real(r8) :: y
!        y=(x**2-2*x*mu+mu**2-sig**2)/(sig**3)
!    end function
!    ! derivative of a gaussian with respect to mu, divided by a gaussian
!    pure function gauss_derivative_mu(x,mu,sig) result(y)
!        real(r8), intent(in) :: x,mu,sig
!        real(r8) :: y
!        y=(x-mu)/(sig**2)
!    end function
!    ! derivative of a lorentzian with respect to sigma, divided by squared lorentzian
!    pure function lorentz_derivative_sig(x,mu,sig) result(y)
!        real(r8), intent(in) :: x,mu,sig
!        real(r8) :: y
!        y=-(lo_pi*0.5_r8)+lo_twopi*( (x-mu)**2 )/( sig**2 )
!    end function
!    ! derivative of a lorentzian with respect to mu, divided by squared lorentzian
!    pure function lorentz_derivative_mu(x,mu,sig) result(y)
!        real(r8), intent(in) :: x,mu,sig
!        real(r8) :: y
!        y=2*lo_twopi*(x-mu)/sig
!    end function
!end subroutine

end submodule
