
!> use the 2D convex hull to get the limits for volume at each temperature. Best would be to use a triangulation, but I don't have the energy to write a triangulation routine so I have to rely on CGAL, and people don't want to compile that.
subroutine volume_limits_from_convex_hull(grid_coordinates,dim_volume,dim_temperature,tmin,tmax,ntpts,tempspacing,minvol,maxvol,verbosity)
    !> grid coordinates
    real(flyt), dimension(:,:), intent(in) :: grid_coordinates
    !> which dimension is volume
    integer, intent(in) :: dim_volume
    !> which dimension is temperature
    integer, intent(in) :: dim_temperature
    !> min temperature
    real(flyt), intent(in) :: tmin
    !> min temperature
    real(flyt), intent(in) :: tmax
    !> number of temperature points
    integer, intent(in) :: ntpts
    !> what kind of spacing in temperature
    character(len=3), intent(in) :: tempspacing
    !> min volume per temperature
    real(flyt), dimension(ntpts), intent(out) :: minvol
    !> max volume per temperature
    real(flyt), dimension(ntpts), intent(out) :: maxvol
    !> talk?
    integer, intent(in) :: verbosity

    real(flyt), parameter :: tol=1E-3_flyt
    real(flyt), dimension(:), allocatable :: temperatures
    real(flyt), dimension(:,:), allocatable :: r,rhull,dr
    integer :: nptot,nhull

    init: block
        integer, dimension(:), allocatable :: hullind
        integer :: i
        ! rearrange it so that volume is always the first dimensions, gets really confusing otherwise.
        ! perhaps this can also be used to project down from 3D-grids? Who knows.

        if ( verbosity .gt. 0 ) then
            write(*,*) '... grabbing volumes from the convex hull'
        endif

        nptot=size(grid_coordinates,2)
        allocate(r(2,nptot))
        do i=1,nptot
            r(1,i)=grid_coordinates(dim_volume,i)
            r(2,i)=grid_coordinates(dim_temperature,i)
        enddo

        ! get the temperature range
        lo_allocate(temperatures(ntpts))
        select case( tempspacing )
        case('lin') ! linearly spaced in temperature
            call lo_linspace(tmin,tmax,temperatures)
        case('log') ! log=spaced
            call logspace(tmin,tmax,temperatures)
        end select

        ! get the convex hull
        call convex_hull_2d(r,1E-3_r8,hullind)
        nhull=size(hullind)
        if ( nhull .lt. 4 ) call lo_stop_gracefully(['Convex hull of gridpoints has less than 4 points, makes no sense'],lo_exitcode_symmetry,__FILE__,__LINE__)
        ! get the points on the hull, and the vector that points to the next point on the hull.
        lo_allocate(rhull(2,nhull))
        lo_allocate(dr(2,nhull))
        rhull=0.0_flyt
        dr=0.0_flyt
        do i=1,nhull
            rhull(:,i)=r(:,hullind(i))
        enddo
        do i=1,nhull-1
            dr(:,i)=rhull(:,i+1)-rhull(:,i)
        enddo
        dr(:,nhull)=rhull(:,1)-rhull(:,nhull)
    end block init

    ! now get the volume limits from this
    findlim: block
        integer :: i,j,t
        real(flyt), dimension(:), allocatable :: dmin,dmax
        real(flyt) :: y,dy,s,x

        minvol=maxval(r(1,:))+1.0_flyt ! -lo_huge
        maxvol=-1.0_flyt !lo_huge
        do t=1,ntpts
            y=temperatures(t)
            ! compare with each line segment
            do i=1,nhull
                ! some things can happen. First is that the line segment is flat with respect to temperature.
                if ( abs( rhull(2,i)-y ) .lt. lo_tol ) then
                    ! catch the case of endpoint of segment on temperature slice
                    minvol(t)=min(minvol(t),rhull(1,i))
                    maxvol(t)=max(maxvol(t),rhull(1,i))
                else
                    ! figure out the intersection between the line
                    ! y1 = rhull(:,i)+s*dr(:,i)
                    ! and constant temperature.
                    dy=dr(2,i)
                    if ( abs(dy) .gt. lo_tol ) then
                        s=(y-rhull(2,i))/dy
                    else
                        s=-1.0_flyt
                    endif
                    if ( s .gt. -lo_tol .and. s .lt. 1.0_flyt+lo_tol ) then
                        x=rhull(1,i)+s*dr(1,i)
                        minvol(t)=min(minvol(t),x)
                        maxvol(t)=max(maxvol(t),x)
                    endif
                endif
            enddo
        enddo

        ! chop this to unique values
        call lo_return_unique(minvol,dmin,lo_tol)
        call lo_return_unique(maxvol,dmax,lo_tol)
        do i=1,size(minvol)
        do j=1,size(dmin)
            if ( abs(minvol(i)-dmin(j)) .lt. lo_tol ) minvol(i)=dmin(j)
        enddo
        enddo
        do i=1,size(maxvol)
        do j=1,size(dmax)
            if ( abs(maxvol(i)-dmax(j)) .lt. lo_tol ) maxvol(i)=dmax(j)
        enddo
        enddo
    end block findlim
end subroutine

!> logarithmically spaced points, slightly smarter in how the shifts/spacings are decided
subroutine logspace(minv,maxv,x)
    !> lower bound
    real(flyt), intent(in) :: minv
    !> upper bound
    real(flyt), intent(in) :: maxv
    !> logarithmically spaced points
    real(flyt), dimension(:), intent(out) :: x

    integer :: n,i
    real(flyt) :: dl
    real(flyt) :: f0,f1,f2,f3
    real(flyt) :: shift,rng,target_dl,inc

    n=size(x,1)
    if ( n .eq. 1 ) then
        x(1)=(minv+maxv)*0.5_flyt
        return
    endif

    rng=maxv-minv
    dl=1.0_flyt/(n-1.0_flyt)
    target_dl=rng/(n*10)
    ! decide on a reasonable shift
    shift=1.0_flyt
    inc=2.0_flyt

    shift=lo_sqtol
    ! what separation do I get between the first two points from this?
    do i=1,10000
        f0=log(shift)
        f1=log(rng+shift)
        f2=(f1-f0)*dl
        f3=exp(f0+f2)-exp(f0)

        if ( abs(f3-target_dl) .lt. 1E-14_flyt ) then
            exit
        endif

        if ( f3 .lt. target_dl ) then
            shift=shift*inc
        else
            shift=shift/inc
            inc=(inc+1)*0.5_flyt
        endif
    enddo

    f0=log(shift)
    f1=log(rng+shift)
    do i=1,n
        x(i)=(f1-f0)*dl*(i-1)+f0
        x(i)=exp(x(i))-shift+minv
    enddo
    x(1)=minv
    x(n)=maxv
end subroutine

!> Calculates the convex hull for a set of 2d points. The algorithm is a Graham scan @cite Graham1972a and checks for counterclockwise turns within a tolerance. The output is an index array to the points on the hull, and optionally an additional array with the points inside the hull.
subroutine convex_hull_2d(points,tol,points_on_hull)
    !> the points
    real(flyt), dimension(:,:), intent(in) :: points
    !> the tolerance
    real(flyt), intent(in) :: tol
    !> the points on the hull
    integer, dimension(:), allocatable, intent(out) :: points_on_hull

    real(r8), dimension(:,:), allocatable :: sorted_points
    integer :: i,np

    np=size(points,2)

    simplestop: block
        integer :: i
        ! maybe we can stop early
        if ( np .le. 3 ) then
            ! The hull is trivial if there are three or less points.
            allocate(points_on_hull(np))
            do i=1,np
                points_on_hull(i)=i
            enddo
            return
        endif
    end block simplestop

    ! Sort points sensibly
    presort: block
        real(r8), dimension(:,:), allocatable :: dr1
        real(r8), dimension(:), allocatable :: dr0,uny,una
        real(r8), dimension(2) :: p0
        real(r8) :: f0
        integer, dimension(:), allocatable :: di
        !integer, dimension(:), allocatable ::
        integer :: i,j,l,ctr
        ! First sort according to y-coordinate?

        allocate(sorted_points(2,np))
        sorted_points=-lo_huge

        ! Grab the unique y-coordinates, sorted by size
        allocate(dr0(np))
        dr0=points(2,:)
        call lo_return_unique(dr0,uny)
        call qsort(uny)

        ctr=0
        do i=1,size(uny)
            l=0
            do j=1,np
                if ( abs(uny(i)-points(2,j)) .lt. lo_tol ) then
                    l=l+1
                    dr0(l)=points(1,j)
                endif
            enddo
            call qsort(dr0(1:l))
            do j=1,l
                ctr=ctr+1
                sorted_points(:,ctr)=[dr0(j),uny(i)]
            enddo
        enddo

        ! Check I got it right
        !do i=1,np
        !    ctr=0
        !    do j=1,np
        !        if ( norm2(points(:,i)-sorted_points(:,j)) .lt. lo_tol ) ctr=ctr+1
        !    enddo
        !enddo

        ! calculate angle to first point
        p0=sorted_points(:,1)
        do i=1,np
            dr0(i)=atan2( sorted_points(2,i)-p0(2),sorted_points(1,i)-p0(1) )
        enddo
        allocate(dr1(2,np))
        allocate(di(np))
        ! Sort according to this angle
        call qsort(dr0,di)
        dr1=sorted_points

        ! make sure we are strong-sorted
        deallocate(uny)
        call lo_return_unique(dr0,uny,1E-10_r8)
        call qsort(uny)

        ctr=0
        do i=1,size(uny)
            l=0
            do j=1,np
                f0=atan2( dr1(2,j)-p0(2),dr1(1,j)-p0(1) )
                if ( abs(f0-uny(i)) .lt. 1E-10_r8 ) then
                    l=l+1
                    di(l)=j
                endif

            enddo
            call qsort(di(1:l))
            do j=1,l
                ctr=ctr+1
                sorted_points(:,ctr)=dr1(:,di(j))
            enddo
        enddo

        !do i=1,np
        !    ctr=0
        !    do j=1,np
        !        if ( norm2(points(:,i)-sorted_points(:,j)) .lt. lo_tol ) ctr=ctr+1
        !    enddo
        !    f0=atan2( sorted_points(2,i)-p0(2),sorted_points(1,i)-p0(1) )
        !    write(*,*) i,ctr,f0
        !enddo
    end block presort

    hull: block
        real(r8), dimension(:,:), allocatable :: dum
        real(r8), dimension(2) :: p0
        integer, dimension(:), allocatable :: ind_on

        integer :: i,j,l

        allocate(dum(2,np))
        allocate(ind_on(np))
        ind_on=-1
        ind_on(1)=1
        dum=-lo_huge
        dum(:,1)=sorted_points(:,1)

        l=1
        do i=2,np
            p0=sorted_points(:,i)
            do
                if ( l .lt. 2 ) exit
                if ( cross2d(dum(:,l-1),dum(:,l),p0,tol) .lt. 0 ) then
                    l=l-1
                    cycle
                else
                    exit
                endif
            enddo
            l=l+1
            dum(:,l)=p0
        enddo

        ! Return the points that are on the hull.
        allocate(points_on_hull(l))
        do i=1,l
            do j=1,np
                if ( norm2(dum(:,i)-points(:,j)) .lt. lo_tol ) then
                    points_on_hull(i)=j
                endif
            enddo
        enddo
    end block hull

    contains
    !> Check if three points make a clockwise or counterclockwise turn
    function cross2d(p1,p2,p3,tolerance) result(r)
        !> first point
        real(flyt), dimension(2), intent(in) :: p1
        !> second point
        real(flyt), dimension(2), intent(in) :: p2
        !> third point
        real(flyt), dimension(2), intent(in) :: p3
        !> tolerance for what is 0
        real(flyt), intent(in) :: tolerance
        !> A value < 0 means it's counterclockwise, > 0 clockwise and 0 that the points are on the same line.
        integer :: r

        real(flyt), dimension(2) :: v1,v2
        real(flyt) :: alpha

        v1=p2-p1
        v2=p2-p3
        ! signed angle between vectors
        alpha=( atan2(v2(2),v2(1))-atan2(v1(2),v1(1)) )*180.0_flyt/lo_pi
        alpha=mod(alpha+360.0_flyt,360.0_flyt)
        if ( abs(alpha-180.0_flyt) .lt. tolerance .or. abs(alpha) .lt. tolerance ) then
            r=0
        elseif ( alpha .gt. 180_flyt ) then
            r=1
        else
            r=-1
        endif
    end function
end subroutine
