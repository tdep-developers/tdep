
!> group calculations in a magnetic meaningful way
subroutine group_distributed_calculations(gs,groupind,groupctr,confctr,ngroup)
    !> grid simulation stuff
    type(lo_gridsim), intent(in) :: gs
    !> group information things
    integer, dimension(:,:,:), allocatable, intent(out) :: groupind
    integer, dimension(:,:), allocatable, intent(out) :: groupctr,confctr
    integer, dimension(:), allocatable, intent(out) :: ngroup

    integer, dimension(:,:), allocatable :: ddi
    integer, dimension(:), allocatable :: ctr,di,dj
    integer, parameter :: grthres=2
    integer :: i,j,k,l,ii,jj,ngmax

    ! max possible number of group members
    lo_allocate(di(gs%nsim))
    di=0
    do i=1,gs%raw%nconf
        di( gs%raw%gridind(i) )=di( gs%raw%gridind(i) )+1
    enddo
    ngmax=maxval(di)
    lo_deallocate(di)

    ! Some initial counting
    lo_allocate(ctr(gs%raw%nconf))
    lo_allocate(di(gs%raw%nconf))
    lo_allocate(ddi(ngmax,gs%raw%nconf))
    ctr=0
    ddi=0
    di=1
    do i=1,gs%raw%nconf
        if ( di(i) .eq. 0 ) cycle
        do j=i+1,gs%raw%nconf
            ii=gs%raw%gridind(i)
            jj=gs%raw%gridind(j)
            if ( gs%raw%gridind(i) .ne. gs%raw%gridind(j) ) cycle
            if ( abs(gs%raw%u(1,1,i)-gs%raw%u(1,1,j)) .gt. lo_sqtol ) cycle
            if ( sum(abs(gs%raw%u(:,:,i)-gs%raw%u(:,:,j)))/gs%na_ss/3 .gt. lo_sqtol ) cycle
            ! ok, these are equal
            di(j)=0
            ctr(i)=ctr(i)+1
            ddi(ctr(i),i)=j
        enddo
    enddo
    lo_deallocate(di)

    ! Build the actual groups, first count how many groups there are per gridpoint.
    lo_allocate(ngroup(gs%nsim))
    ngroup=0
    do i=1,gs%raw%nconf
        ii=gs%raw%gridind(i)
        if ( ctr(i) .gt. grthres ) ngroup(ii)=ngroup(ii)+1
    enddo

    ! Then sort everything out into groups
    lo_allocate(groupctr(maxval(ngroup),gs%nsim))
    lo_allocate(groupind(maxval(ctr),maxval(ngroup),gs%nsim))
    lo_allocate(di(gs%nsim))
    lo_allocate(confctr(maxval(ngroup),gs%nsim))
    groupctr=0
    groupind=0
    confctr=0
    di=0
    do i=1,gs%raw%nconf
        if ( ctr(i) .le. grthres ) cycle
        ii=gs%raw%gridind(i)
        di(ii)=di(ii)+1
        groupctr(di(ii),ii)=ctr(i)
        groupind( 1:groupctr(di(ii),ii),di(ii),ii )=ddi( 1:ctr(i),i )
    enddo

    ! Finally count the number of configurations per group
    do ii=1,gs%nsim
    do jj=1,ngroup(ii)
        k=0
        do i=1,groupctr(jj,ii)
        do j=i+1,groupctr(jj,ii)
            k=k+1
        enddo
        enddo
        confctr(jj,ii)=k
    enddo
    enddo

end subroutine



!> optimize a reference temperature for the magnetic things. Don't try to understand.
subroutine optimize_tref_for_sigma(gs,sigma,tref)
    !> grid simulation stuff
    type(lo_gridsim), intent(in) :: gs
    !> sigma's to optimize
    real(flyt), dimension(:,:), intent(in) :: sigma
    !> reference temperature
    real(flyt), intent(inout) :: tref

    integer, parameter :: nsweep=50
    integer, parameter :: polyorder=2
    real(flyt), dimension(nsweep) :: tsweep
    real(flyt) :: grad,hess,e0,e1
    
    integer :: i,ctr

    ! Starting guess, sweep across temperatures
    call lo_linspace(10.0_flyt,2000.0_flyt,tsweep)
    e0=lo_huge
    do i=1,nsweep
        e1=deviation(gs%grid_coordinates,sigma,tsweep(i))
        if ( e1 .lt. e0 ) then
            tref=tsweep(i)
            e0=e1
        endif
    enddo

    ! Newton minimize from there
    ctr=0
    do i=1,100
        e0=deviation(gs%grid_coordinates,sigma,tref)
        call gradhess(gs%grid_coordinates,sigma,tref,grad,hess)
        if ( abs(grad/hess) .lt. 1E-3_flyt ) then
            ctr=ctr+1
        else
            ctr=0
        endif
        if ( ctr .eq. 4 ) exit
        tref=tref+grad/hess
    enddo
    contains

    subroutine gradhess(coord,vals,tref,grad,hess)
        real(flyt), dimension(:,:), intent(in) :: coord
        real(flyt), dimension(:,:), intent(in) :: vals
        real(flyt), intent(in) :: tref
        real(flyt), intent(out) :: grad
        real(flyt), intent(out) :: hess

        integer, parameter :: n=6,m=2*n+1
        real(flyt), dimension(m,4) :: sc
        real(flyt), dimension(m) :: fv
        integer :: i

        call lo_centraldifference(n,tref,1E-3_flyt,sc)
        do i=1,m
            fv(i)=deviation(coord,vals,sc(i,1))
        enddo
        grad=-sum(sc(:,2)*fv)
        hess=sum(sc(:,3)*fv)
    end subroutine

    function deviation(coord,vals,tref) result(error)
        real(flyt), dimension(:,:), intent(in) :: coord
        real(flyt), dimension(:,:), intent(in) :: vals
        real(flyt), intent(in) :: tref
        real(flyt) :: error

        type(lo_polynomial) :: poly
        real(flyt), dimension(:,:), allocatable :: transformed_coordinates,coeff
        real(flyt), dimension(:,:), allocatable :: transformed_values
        real(flyt) :: f0,f1
        integer :: i,j,k,l

        ! transform the coordinates
        lo_allocate(transformed_coordinates(size(coord,1),size(coord,2)))
        lo_allocate(transformed_values(size(vals,1),size(vals,2)))
        transformed_coordinates=coord
        if ( gs%info%dim_temperature .gt. 0 ) then
            do i=1,size(coord,2)
                transformed_values(i,:)=vals(i,:)*magtempscaler(coord(gs%info%dim_temperature,i),tref)
                if ( gs%info%temperature_scale .gt. 0.0_flyt ) then
                    transformed_coordinates(gs%info%dim_temperature,i)=tempscaler(coord(gs%info%dim_temperature,i),gs%info%temperature_scale)
                endif
            enddo
        endif 

        ! Init polynomial
        call poly%init( polyorder,gs%ndim,transformed_coordinates,gs%info%dimension_names,gs%info%order_per_dim )
        lo_allocate(coeff(poly%ncoeff,size(vals,2)))
        ! Fit polynomial
        do i=1,size(vals,2)
            call linear_least_squares(poly%coeffM,transformed_values(:,i),coeff(:,i))
        enddo
        ! Evaluate
        error=0.0_flyt
        f0=0.0_flyt
        f1=0.0_flyt
        do j=1,size(vals,2)
        do i=1,size(coord,2)
            f0=f0+(transformed_values(i,j)-poly%eval(transformed_coordinates(:,i),coeff(:,j)))**2
            f1=f1+transformed_values(i,j)**2
        enddo
        enddo
        error=f0/f1
    end function

end subroutine

!> Take a list of magnetic moments and fit the onsite model to that.
function fit_onsite_model_to_histogram(magnetic_moments) result(par)
    !> list of magnetic moments
    real(flyt), dimension(:), intent(in) :: magnetic_moments
    !> model parameters
    real(flyt), dimension(2) :: par

    integer, parameter :: nhist=100
    real(flyt), dimension(nhist) :: hx,hy,cdfy,cdfx

    ! First build the histogram. Not sure how to do it smoothly, but I'll figure it out.
    hist: block
        real(flyt) :: xmin,xmax,sigma,foursigma,stddev,invf
        real(flyt) :: invh
        integer :: i,j,k,l,ii,jj,u

        ! Build the cumulative distribution function
        xmin=minval(magnetic_moments)*0.98_flyt
        xmax=maxval(magnetic_moments)*1.02_flyt
        call lo_linspace(xmin,xmax,cdfx)
        do i=1,nhist
            cdfy(i)=sum(magnetic_moments,magnetic_moments<cdfx(i))
        enddo
        cdfy=cdfy/cdfy(nhist)

        ! Get the pdf from this.
        hx=cdfx
        hy=0.0_flyt
        invh=1.0_flyt/(hx(2)-hx(1))
        ! get the edges first
        hy(1)=(cdfy(2)-cdfy(1))*invh
        hy(nhist)=( cdfy(nhist)-cdfy(nhist-1) )*invh
        ! then all the middle guys
        do i=2,nhist-1
            hy(i)= (cdfy(i+1)-cdfy(i-1))*0.5_flyt*invh
        enddo
    end block hist

    ! Fit the model.
    ! Now try to fit it to something
    fithist: block
        integer, parameter :: nalg=4
        real(flyt) :: f0,beta1,beta2
        real(flyt) :: mu0,sig0,e1,e2
        real(flyt), dimension(2,2) :: hess
        real(flyt), dimension(2) :: grad
        real(flyt), dimension(2,nalg) :: steps
        real(flyt), dimension(3) :: y
        integer :: i,j,iter,ctr

        ! Starting guess for the mean
        mu0=lo_trapezoid_integration(hx,hx*hy)
        ! Starting guess for sigma. Hmm.
        e1=sqrt(lo_trapezoid_integration(hx,hy*(hx-mu0)**2))
        e1=e1*2*Sqrt(2*log(2.0_flyt))
        sig0=-4*Log(2.0_flyt)/( e1**4-4*e1**2*mu0**2)

        ! Get the starting point
        par=[mu0,sig0]
        ctr=0
        do iter=1,100
            grad=lsqgrad(hx,hy,par)
            hess=lsqhess(hx,hy,par)
            e1=lsqerr(hx,hy,par)
            ! Check eigenvalues of Hessian. If negative, replace with identity.
            if ( minval(eigval2x2matrix(hess)) .lt. 1E-10_flyt ) then
                f0=norm2(hess)
                call lo_identitymatrix(hess)
                hess=hess*f0
            endif

            ! Design a few different steps to try. First one is just steepest descent
            steps(:,1)=grad/norm2(grad)
            steps(:,2)=matmul(inv2x2matrix(hess),grad)
            steps(:,2)=steps(:,2)/norm2(steps(:,2))
            steps(:,3)=steps(:,1)*[1.0_flyt,0.0_flyt]
            steps(:,4)=steps(:,1)*[0.0_flyt,1.0_flyt]

            ! Decide on step length
            do i=1,nalg
                beta1=dot_product(grad,steps(:,i))
                beta2=dot_product(steps(:,i),matmul(hess,steps(:,i)))
                if ( abs(beta2) .lt. 1E-10_flyt ) then
                    !write(*,*) iter,norm2(grad),norm2(hess),beta2,norm2(steps(:,i))
                    !write(*,*) 'THINKING IS HARD, CAN NOT HAPPEN, OPTIMIZING HISTOGRAM'
                    beta2=1.0_flyt
                    !stop
                    !else
                endif
                steps(:,i)=-steps(:,i)*beta1/beta2
                ! Just make sure no step is really big
                if ( norm2(steps(:,i)) .gt. 1E-2*norm2(par) ) then
                    steps(:,i)=steps(:,i)*1E-2_flyt*norm2(par)/norm2(steps)
                endif
            enddo

            ! Check which one is the best!
            j=0
            e2=lo_huge
            do i=1,nalg
                f0=lsqerr(hx,hy,par+steps(:,i))
                if ( f0 .lt. e2 ) then
                    j=i
                    e2=f0
                endif
            enddo
            ! Update position
            par=par+steps(:,j)

            if ( norm2(steps(:,j))/norm2(par) .lt. 1E-5_flyt ) then
                ctr=ctr+1
            else
                ctr=0
            endif
!write(*,*) 'iter',iter,norm2(grad),e1

            if ( ctr .eq. 4 ) exit

        enddo
    end block fithist
contains

!> invert a 2x2-matrix
function inv2x2matrix(m) result(n)
    real(flyt), dimension(2,2), intent(in) :: m
    real(flyt), dimension(2,2) :: n

    real(flyt) :: a,b,c,d,dt

    a=m(1,1)
    b=m(2,1)
    c=m(1,2)
    d=m(2,2)
    dt=a*d-c*b
    if ( abs(dt) .gt. lo_sqtol ) then
        dt=1.0_flyt/dt
        n(1,1)=d*dt
        n(2,1)=-c*dt
        n(1,2)=-b*dt
        n(2,2)=a*dt
    else
        n(1,1)=1.0_flyt
        n(2,1)=0.0_flyt
        n(1,2)=0.0_flyt
        n(2,2)=1.0_flyt
    endif
end function

!> eigenvalues of a 2x2-matrix
function eigval2x2matrix(m) result(eig)
    real(flyt), dimension(2,2), intent(in) :: m
    real(flyt), dimension(2) :: eig

    real(flyt) :: a,b,c,d,dt
    a=m(1,1)
    b=m(2,1)
    c=m(1,2)
    d=m(2,2)
    dt=sqrt(4*b*c+(a-d)**2)
    eig(1)=(a+d-dt)*0.5_flyt
    eig(2)=(a+d+dt)*0.5_flyt
end function

!> error function
function lsqerr(hx,hy,par) result(errorsq)
    real(flyt), dimension(:), intent(in) :: hx,hy
    real(flyt), dimension(2), intent(in) :: par
    real(flyt) :: errorsq

    real(flyt) :: mu,sig,gmu,gsig,nrm
    real(flyt) :: f0,f1,f2
    integer :: i

    nrm=normfactor(par)

    ! Squared error and gradient
    mu=par(1)
    sig=par(2)
    errorsq=0.0_flyt
    do i=1,nhist
        f0=distr(hx(i),mu,sig)*nrm
        errorsq=errorsq+( hy(i) - f0 )**2
    enddo
end function

!> gradient of the error function
function lsqgrad(hx,hy,par) result(grad)
    real(flyt), dimension(:), intent(in) :: hx,hy
    real(flyt), dimension(2), intent(in) :: par
    real(flyt), dimension(2) :: grad

    integer, parameter :: n=6,m=2*n+1
    real(flyt), dimension(m,4) :: sc
    real(flyt), dimension(m) :: fv
    integer :: i

    call lo_centraldifference(n,par(1),1E-5_flyt,sc)
    do i=1,m
        fv(i)=lsqerr(hx,hy,[sc(i,1),par(2)])
    enddo
    grad(1)=sum(sc(:,2)*fv)

    call lo_centraldifference(n,par(2),1E-5_flyt,sc)
    do i=1,m
        fv(i)=lsqerr(hx,hy,[par(1),sc(i,1)])
    enddo
    grad(2)=sum(sc(:,2)*fv)
end function

!> hessian of the error function
function lsqhess(hx,hy,par) result(hess)
    real(flyt), dimension(:), intent(in) :: hx,hy
    real(flyt), dimension(2), intent(in) :: par
    real(flyt), dimension(2,2) :: hess

    integer, parameter :: n=4,m=2*n+1
    real(flyt), parameter :: delta=1E-6_flyt
    real(flyt), dimension(m*m,6) :: dm
    real(flyt), dimension(m*m,1) :: fv
    real(flyt) :: x,y,dx,dy
    integer :: i,j,l

    l=0
    do i=-n,n
    do j=-n,n
        l=l+1
        dx=i*delta
        dy=j*delta
        x=par(1)+dx
        y=par(2)+dy
        fv(l,1)=lsqerr(hx,hy,[x,y])
        dm(l,:)=[1.0_flyt,dx,dy,dx*dy,dx*dx,dy*dy]
    enddo
    enddo
    call lo_dgels(dm,fv)

    hess(1,1)=2*fv(5,1)
    hess(2,2)=2*fv(6,1)
    hess(1,2)=fv(4,1)
    hess(2,1)=fv(4,1)
end function

!> normalization factor for my distribution
function normfactor(par) result(nrm)
    real(flyt), dimension(2), intent(in) :: par
    real(flyt) :: nrm

    integer, parameter :: nquad=30
    real(flyt), dimension(2,nquad) :: gq
    real(flyt) :: x0,x1,x2,x3,z1,z2,z3,w1,w2,w3
    integer :: i

    ! Get the quadrature weights
    call lo_gaussianquadrature(nquad,0.0_flyt,1.0_flyt,gq)

    x0=0.0_flyt
    z1=sqrt(log(2.0_flyt)/par(2))
    if ( par(1)**2 .gt. z1 ) then
        x1=sqrt(par(1)**2-z1)
    else
        x1=par(1)*0.5_flyt
    endif
    x2=sqrt(par(1)**2+z1)
    x3=sqrt(par(1)**2+sqrt(log(1E14_flyt)/par(2)))

    nrm=0.0_flyt
    do i=1,nquad
        z1=gq(1,i)*(x1-x0)+x0
        z2=gq(1,i)*(x2-x1)+x1
        z3=gq(1,i)*(x3-x2)+x2
        w1=gq(2,i)*(x1-x0)
        w2=gq(2,i)*(x2-x1)
        w3=gq(2,i)*(x3-x2)
        nrm=nrm+distr( z1,par(1),par(2) )*w1
        nrm=nrm+distr( z2,par(1),par(2) )*w2
        nrm=nrm+distr( z3,par(1),par(2) )*w3
    enddo
    nrm=1.0_flyt/nrm
end function

function distr(x,mu,sig) result(y)
    real(flyt) :: x,mu,sig
    real(flyt) :: y

    real(flyt) :: arg
    arg=-sig*(x**2-mu**2)**2
    if ( arg .lt. -40.0_flyt ) then
        y=0.0_flyt
    else
        y=exp(arg)
    endif
end function

end function

