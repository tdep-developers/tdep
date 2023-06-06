#include "precompilerdefinitions"
module glmnet_wrappers
!! Wrapper around the ancient fortran in GLMNET
use konstanter, only: flyt,lo_sqtol,lo_status,lo_exitcode_param,lo_tiny,lo_huge,lo_hugeint
use gottochblandat, only: walltime,lo_mean,lo_stddev,lo_progressbar_init,lo_progressbar,lo_linspace,&
                          lo_linear_least_squares,lo_return_unique,qsort,lo_enforce_linear_constraints
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use type_blas_lapack_wrappers, only: lo_gemv,lo_gemm
implicit none
private

public :: lo_elastic_net_helper
!public :: glmnet_regression_prefixed
public :: glmnet_regression
public :: glmnet_varselection
public :: glmnet_linear_system
public :: elast1

! Default number of lambda-values in the regression
integer, parameter :: nlambda=500
! These are parameters I have chosen to make fixed for the elastic net
! regression. As soon as I need to loosen them and use them, I will fix.
real(flyt), parameter :: thr=1E-15_flyt
real(flyt), parameter :: flmindef=1E-13_flyt
integer, parameter :: isd=1
integer, parameter :: intr=1
integer, parameter :: ka=1
integer, parameter :: maxit=1000000

!> Wrapper to contain all the settings of an elastic net regression, I think.
type lo_elastic_net_helper
    !> number of observations (no)
    integer :: nobs=-lo_hugeint
    !> number of variables (ni)
    integer :: nx=-lo_hugeint
    !> number of lambda-values
    integer :: nlambda=-lo_hugeint
    !> user-supplied lambda or automatic?
    logical :: autolambda=.false.
    !> alpha-parameter for the switching
    real(flyt) :: alpha=-lo_huge
    !> smallest allowed lambda-value
    real(flyt) :: smallest_lambda=-lo_huge
    !> threshold
    real(flyt) :: threshold=-lo_huge
    !> max number of iterations
    integer :: maxit
    !> List of solutions
    real(flyt), dimension(:,:), allocatable :: solutions
    !> Resulting lambdavalues
    real(flyt), dimension(:), allocatable :: lambdavals
    !> Rsquare per lambda
    real(flyt), dimension(:), allocatable :: rsq
    !> Input lambda
    real(flyt), dimension(:), allocatable :: input_lambda
    !> order of variables
    integer, dimension(:), allocatable :: varorder
    !> number of variables per lambda
    integer, dimension(:), allocatable :: nvar
    contains
        procedure :: init=>init_elastic_net
end type

contains

!> Linear system thingy with glmnet
subroutine glmnet_regression(A,B,alpha,varind,solutions,lambdavals,inlambda)
    !> Coefficient matrix
    real(flyt), dimension(:,:), intent(in) :: A
    !> Response vector
    real(flyt), dimension(:), intent(in) :: B
    !> Variable that goes between ridge and lasso. 0 is ridge, 1 is lasso.
    real(flyt), intent(in) :: alpha
    !> Possible variable selection schemes
    integer, dimension(:,:), allocatable, intent(out), optional :: varind
    !> Solution vectors
    real(flyt), dimension(:,:), allocatable, intent(out), optional :: solutions
    !> Lambda-values used for these solutions
    real(flyt), dimension(:), allocatable, intent(out), optional :: lambdavals
    !> Input lambda values!
    real(flyt), dimension(:), intent(in), optional :: inlambda

    type(lo_elastic_net_helper) :: eh
    integer, dimension(:), allocatable :: order_of_importance

    slv: block
        ! Things GMLnet wants
        real(flyt), dimension(:,:), allocatable :: wX,ca,cl
        real(flyt), dimension(:), allocatable :: wW,wY
        real(flyt), dimension(:), allocatable :: vp,ulam,a0,rsq,alm
        real(flyt) :: parm,flmin
        integer, dimension(:), allocatable :: ia,nin
        integer :: no,ni,nx,nlam
        integer :: ne,lmu,nlp,jerr
        ! My things
        real(flyt) :: f0,f1
        integer :: i,j,k

        ! Try to set all the input variables. I have tried to sort them into two categories:
        ! 1) Those that are no choice
        ! 2) The ones that matter
        ! 3) The ones that seem to do nothing

        ! So, first the ones that are no choice
        no=size(A,1)       ! Number of observations
        ni=size(A,2)       ! Number of independent
        if ( present(inlambda) ) then
            nlam=size(inlambda)
        else
            nlam=nlambda       ! Number of lambda-values for the solution path, max
        endif
        ne=ni+10            ! not sure what this is, but someone set this as the default
        nx=ni              ! max number of variables to return
        lo_allocate(wX(no,ni))
        lo_allocate(wW(no))
        lo_allocate(wY(no))
        wW=1.0_flyt ! Weights
        wX=A
        wY=B
        ! alpha between LASSO and Thikonov
        parm=alpha

        ! These don't seem to matter too much.
        lo_allocate(cl(2,ni)) ! Some kind of range. Will not use.
        cl(1,:)=-lo_huge
        cl(2,:)=lo_huge
        lo_allocate(vp(ni))   ! Strange penalty thing. Will not use.
        vp=1.0_flyt
        lmu=-1
        lo_allocate(ulam(nlam))
        if ( present(inlambda) ) then
            flmin=2.0_flyt   ! Means user-supplied lambda
            ulam=inlambda
        else
            flmin=flmindef   ! Smallest lambda, I think.
            ulam=-1.0_flyt
        endif

        ! I think these things are output.
        !===========
        lo_allocate(ca(nx,nlam))
        lo_allocate(a0(nlam))
        lo_allocate(rsq(nlam))
        lo_allocate(alm(nlam))
        ca=1.0_flyt
        a0=-1.0_flyt
        rsq=-1.0_flyt
        alm=-1.0_flyt
        ! Integer things
        lo_allocate(ia(nx))
        lo_allocate(nin(nlam))
        ia=-1
        nin=-1
        jerr=0

        !call eh%init(ka,parm,no,ni,wX,wY,wW,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
        call eh%init(parm,no,ni,wX,wY,wW,nlam,flmin,ulam,thr,intr,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
        ! Run the magical GLMNET thingy
        ! call elnet(ka,parm,no,ni,wX,wY,wW,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
        if ( jerr .ne. 0 ) then
            write(*,*) 'GLMNET error code:',jerr
        endif

        ! Collect useful information.
        lo_allocate(order_of_importance(ni))
        order_of_importance=ia

        if ( present(solutions) ) then
            lo_allocate(solutions(nx,lmu))
            solutions=0.0_flyt
            do i=1,lmu
                do j=1,nin(i)
                    k=ia(j)
                    ! I'm not 100% sure about the intercept thingy. Maybe it should be ignored.
                    solutions(k,i)=ca(j,i)+a0(i)
                enddo
            enddo
        endif

        if ( present(lambdavals) ) then
            lo_allocate(lambdavals(lmu))
            lambdavals=0.0_flyt
            ! Fix the first lambda
            f0=log(alm(2))
            f1=log(alm(3))
            lambdavals(1)=exp(2*f0-f1)
            do i=2,lmu
                lambdavals(i)=alm(i)
            enddo
        endif

    end block slv

    ! Return useful things
    if ( present(varind) ) then
    retstuff: block
        integer :: i,j,l
        integer :: nimp,nx,nperm

        ! Number of possible variables
        nx=size(A,2)
        ! Number I found reasonably important
        nimp=count(order_of_importance>0)
        ! How many permutations of variables I will return
        if ( nimp .lt. nx ) then
            nperm=nimp+1
        else
            nperm=nimp
        endif
        ! Start with nothing
        lo_allocate(varind(nx,nperm))
        varind=-1
        l=0
        ! First one will be all the variables
        if ( nperm .gt. nimp ) then
            l=l+1
            do i=1,nx
                varind(i,l)=i
            enddo
        endif
        do i=nimp,1,-1
            l=l+1
            do j=1,i
                varind(j,l)=order_of_importance(j)
            enddo
            call qsort(varind(1:i,l))
        enddo
    end block retstuff
    endif
end subroutine

!> Linear system with specific lambda and alpha-values
subroutine glmnet_linear_system(A,B,alpha,lambda,X,maxlambda,constr,nconstr)
    !> Coefficient matrix
    real(flyt), dimension(:,:), intent(in) :: A
    !> Response vector
    real(flyt), dimension(:), intent(in) :: B
    !> Variable that goes between ridge and lasso. 0 is ridge, 1 is lasso.
    real(flyt), intent(in) :: alpha
    !> Regularization constant
    real(flyt), intent(in) :: lambda
    !> Solution
    real(flyt), dimension(:), intent(out) :: X
    !> Optionally, max lambda
    real(flyt), intent(in), optional :: maxlambda
    !> Optionally, linear zeroconstraints
    real(flyt), dimension(:,:), intent(in), optional :: constr
    !> Optionally, number of linear zeroconstraints
    integer, intent(in), optional :: nconstr

    integer, parameter :: nlm=50 ! Number of Lambda in the path. Chosen at random.
    real(flyt), parameter :: largelambda=1E10_flyt ! What is a large value of lambda? Have to think about this one
    real(flyt), dimension(:,:), allocatable :: solutions
    real(flyt), dimension(nlm) :: inlambda
    integer :: nc

    ! Figure out constraints
    if ( present(constr) ) then
        if ( .not.present(nconstr) ) then
            call lo_stop_gracefully(['Need to provide the number of constraints'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        nc=nconstr
        if ( nc .gt. 0 ) then
        if ( nc .ne. size(constr,1) ) then
            call lo_stop_gracefully(['Inconsistent number of constraints'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        endif
    else
        nc=0
    endif

    ! Set the solution path.
    if ( lambda .lt. lo_tiny ) then
        call lo_stop_gracefully(['Lambda needs to be larger than zero.'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    if ( present(maxlambda) ) then
        if ( maxlambda .le. lambda ) then
            call lo_stop_gracefully(['Max lambda needs to be larger lambda. Much larger.'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        call lo_linspace(log(maxlambda),log(lambda),inlambda)
    else
        call lo_linspace(log(lambda*largelambda),log(lambda),inlambda)
    endif
    inlambda=exp(inlambda)
    inlambda(nlm)=lambda

    ! Solve with this path
    call glmnet_regression(A,B,alpha,inlambda=inlambda,solutions=solutions)
    ! Grab the solution we want
    X=solutions(:,nlm)
    ! Fix the constraints
    if ( nc .gt. 0 ) then
        call lo_enforce_linear_constraints(constr,X,nosquare=.true.)
    endif
    ! And we are done!
    deallocate(solutions)
end subroutine

!> Get a resonable list for variable selection.
subroutine glmnet_varselection(A,B,nalpha,varselection,mw,verbosity,tolerance)
    !> Coefficient matrix A
    real(flyt), dimension(:,:), intent(in) :: A
    !> Response vector B
    real(flyt), dimension(:), intent(in) :: B
    !> How many alpha-values to check
    integer, intent(in) :: nalpha
    !> Variable selection thingy
    integer, dimension(:,:), allocatable, intent(out) :: varselection
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Talk a lot?
    integer, intent(in) :: verbosity
    !> Add a tolerance
    real(flyt), intent(in), optional :: tolerance

    real(flyt), parameter :: maxerr=1.03_flyt ! Seems reasonable? No?
    real(flyt), dimension(:), allocatable :: alphavals,viol
    real(flyt) :: lsqerr,lsqdev,tol,timer
    integer :: ne,nx,nprog

    init: block
        real(flyt), dimension(:), allocatable :: X

        ! Shorthand for the sizes
        ne=size(A,1)
        nx=size(A,2)

        ! Start reporting progress
        if ( verbosity .gt. 0 ) then
            nprog=nalpha+2
            timer=walltime()
            call lo_progressbar_init()
        endif

        ! Set the alpha values
        lo_allocate(alphavals(nalpha))
        call lo_linspace(0.0_flyt,1.0_flyt,alphavals)

        ! Make space for a response vector?
        lo_allocate(viol(ne))
        lo_allocate(X(nx))
        X=0.0_flyt
        viol=0.0_flyt

        ! Get the least squares baseline
        call lo_linear_least_squares(A,B,X)
        ! Get the error
        call lo_gemv(A,X,viol)
        ! Least squares error
        viol=(B-viol)**2
        lsqerr=lo_mean(viol)
        lsqdev=lo_stddev(viol)

        ! And set the tolerance
        if ( present(tolerance) ) then
            tol=tolerance
        else
            tol=lo_sqtol
        endif

        ! Report some progress
        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... variable selection',1,nprog,walltime()-timer)
        endif

    end block init

    sel: block
        real(flyt), dimension(:,:), allocatable :: XX
        real(flyt) :: f0
        integer, dimension(:,:), allocatable :: di,dj
        integer, dimension(:), allocatable :: ind1,ind2
        integer :: ial,i,j,k,l

        ! Space for the all possible variable selections
        lo_allocate(di(nx,nalpha*nlambda+1))
        di=-1

        do ial=1,nalpha
            ! Make it MPI parallel
            if ( mod(ial,mw%n) .ne. mw%r ) cycle
            ! Do the elastic net regression thingy
            call glmnet_regression(A,B,alphavals(ial),solutions=XX)
            ! Check which ones are reasonable
            do i=1,size(XX,2),1
                ! Calculate violation?
                call lo_gemv(A,XX(:,i),viol)
                viol=(B-viol)**2
                f0=lo_mean(viol)/lsqerr

                ! Add resonable ones
                if ( f0 .lt. maxerr ) then
                    ! Store this variable selection
                    l=(ial-1)*(nlambda)+i
                    k=0
                    do j=1,nx
                        if ( abs(XX(j,i)) .gt. lo_tiny ) then
                            k=k+1
                            di(k,l)=j
                        endif
                    enddo
                    ! And sort the selection in order so that I can pick out the unique later.
                    if ( k .gt. 0 ) call qsort(di(1:k,l))
                endif
            enddo
            ! Report some progress
            if ( verbosity .gt. 0 ) then
                call lo_progressbar(' ... variable selection',ial+1,nprog,walltime()-timer)
            endif
        enddo
        ! Now collect all of them
        call mw%allreduce('max',di)
        ! Make sure that the one with all variables is in there
        do i=1,nx
            di(i,nalpha*nlambda+1)=i
        enddo
        ! Return the unique
        call lo_return_unique(di,dj)
        lo_allocate(ind1(size(dj,2)))
        lo_allocate(ind2(size(dj,2)))
        ! Pick out the meaningful
        do i=1,size(dj,2)
            ind1(i)=count(dj(:,i)>0)
        enddo
        ind1=-ind1
        call qsort(ind1,ind2)
        ind1=-ind1
        ! Prepare output
        l=count(ind1>0)
        lo_allocate(varselection(nx,l))
        varselection=0
        l=0
        do i=1,size(ind1)
            if ( ind1(i) .le. 0 ) cycle
            l=l+1
            varselection(:,l)=dj(:,ind2(i))
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... variable selection',nprog,nprog,walltime()-timer)
        endif

    end block sel
end subroutine

!> Initialize the elastic net solver
subroutine init_elastic_net(eh,alpha,nobs,nxtot,matrixA,vectorB,weights,nlam,flmin,ulam,thr,intr,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
    !> Solver helper
    class(lo_elastic_net_helper), intent(out) :: eh
    ! Input
    !> Number of observations, in total
    integer, intent(in) :: nobs
    !> Number of variables, in total
    integer, intent(in) :: nxtot
    !> Switching parameter between ridge and lasse
    real(flyt), intent(in) :: alpha
    integer, intent(in) :: nlam,intr
    real(flyt), intent(in) :: thr,flmin
    real(flyt), dimension(:), intent(in) :: ulam
    ! Output
    real(flyt), dimension(:,:), intent(out) :: ca
    real(flyt), dimension(:), intent(out) :: a0,rsq,alm
    integer, dimension(:), intent(out) :: ia,nin
    integer, intent(out) :: jerr,lmu,nlp
    ! Inout
    !> Coefficient matrix
    real(flyt), dimension(:,:), intent(inout) :: matrixA
    !> Measured values
    real(flyt), dimension(:), intent(inout) :: vectorB
    !> Weights per observation
    real(flyt), dimension(:), intent(inout) :: weights

    ! Set the parameters
    init: block
        !integer, dimension(:), allocatable :: ju
        !integer :: nlambda_out
        ! Total number of observations
        eh%nobs=nobs
        ! Number of equations/observations
        eh%nx=nxtot
        ! Number of lambda-values
        eh%nlambda=nlam
        ! Alpha-parameter for the switching
        eh%alpha=alpha
        ! Smallest lambda
        eh%smallest_lambda=flmin !1E-13_flyt !flmin
        ! Tolerance
        eh%threshold=thr ! 1E-15_flyt
        ! Max number of iterations
        eh%maxit=eh%nx*100000
        ! Space for dummy arrays
        lo_allocate(eh%input_lambda(eh%nlambda))
        lo_allocate(eh%solutions(eh%nx,eh%nlambda))
        lo_allocate(eh%rsq(eh%nlambda))
        lo_allocate(eh%lambdavals(eh%nlambda))
        lo_allocate(eh%varorder(eh%nx))
        lo_allocate(eh%nvar(eh%nlambda))
        eh%input_lambda=0.0_flyt
        eh%solutions=0.0_flyt
        eh%rsq=0.0_flyt
        eh%lambdavals=0.0_flyt
        eh%varorder=0
        eh%nvar=0
        if ( flmin .gt. 1.0_flyt ) then
            eh%autolambda=.false.
            eh%input_lambda=ulam
        else
            eh%autolambda=.true.
            eh%input_lambda=0.0_flyt
        endif

        !! So what I think chkvars does is check wether the correct number for a variable is zero.
        !! think this can be excluded and done before I do anything at all, like in removed from
        !! the coefficient matrix from the start.
    end block init

    enu: block
        real(flyt), dimension(:,:), allocatable :: ATA
        real(flyt), dimension(:), allocatable :: xm,xs,xv,vlam,g,v
        real(flyt) :: ym,ys
        integer :: j,k,l,nk

            lo_allocate(ATA(eh%nx,eh%nx))
            lo_allocate(g (eh%nx))
            lo_allocate(xm(eh%nx))
            lo_allocate(xs(eh%nx))
            lo_allocate(xv(eh%nx))
            lo_allocate(vlam(eh%nlambda))
            g =0.0_flyt
            xm=0.0_flyt
            xs=0.0_flyt
            xv=0.0_flyt
            vlam=0.0_flyt

            ! Start standard. This replaces the call below. This standardizes the data.
            !call standard(eh%nobs,eh%nx,matrixA,vectorB,weights,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)
            lo_allocate(v(eh%nobs))
            ! Standardize the weights
            ! weights=weights/sum(weights)
            ! v=sqrt(weights)
            weights=1.0_flyt
            v=1.0_flyt
            call lo_gemv(matrixA,v,xm,trans='T')
            do j=1,eh%nx
                ! Shift them to the mean, I think
                matrixA(:,j)=v*(matrixA(:,j)-xm(j))
                ! some kind of norm per variable
                xv(j)=dot_product(matrixA(:,j),matrixA(:,j))
                xs(j)=sqrt(xv(j))
            enddo

            do j=1,eh%nx
                matrixA(:,j)=matrixA(:,j)/xs(j)
            enddo

            xv=1.0_flyt
            ym=dot_product(weights,vectorB)
            vectorB=v*(vectorB-ym)
            ys=sqrt(dot_product(vectorB,vectorB))
            vectorB=vectorB/ys
            g=0.0_flyt
            call lo_gemv(matrixA,vectorB,g,trans='T')
            deallocate(v)
            ! done with standardization

            ! Some kind of sanity test
            if ( jerr .ne. 0 ) return

            ! Remember to scale lambda!
            if ( eh%autolambda .eqv. .false. ) eh%input_lambda=eh%input_lambda/ys

            ! Get the Gramian
            ATA=0.0_flyt
            call lo_gemm(matrixA,matrixA,ATA,transa='T')

            ! Now things should be scaled properly. Call the big routine thingy.
            call elast1(eh,ATA,g,lmu,nlp,jerr)
            if( jerr .gt. 0 ) then
                ! Throw error code
            endif
            ! Store solution in the normal arrays
            ia=eh%varorder
            nin=eh%nvar
            ca=eh%solutions
            rsq=eh%rsq
            alm=eh%lambdavals

            ! Unscale the solution?
            do k=1,lmu
                alm(k)=ys*alm(k)
                nk=nin(k)
                do l=1,nk
                    ca(l,k)=ys*ca(l,k)/xs(ia(l))
                enddo
                a0(k)=0.0_flyt
                if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))
            enddo
            deallocate(xm,xs,g,xv,vlam)
            return
    end block enu
end subroutine

!> routine yanked from glmnet and translated to readable fortran.
subroutine elast1(eh,ATA,ATB,lmu,nlp,jerr)
    ! elastic net helper
    class(lo_elastic_net_helper), intent(inout) :: eh
    !> How it went
    integer, intent(out) :: jerr,lmu,nlp
    ! Inout
    real(flyt), dimension(:,:), intent(in) :: ATA
    real(flyt), dimension(:), intent(inout) :: ATB

    real(flyt), dimension(:,:), allocatable :: c
    real(flyt), dimension(:), allocatable :: a,da
    integer, dimension(:), allocatable :: mm,ia
    real(flyt) :: ab,ak,alf,alm,dem,dlx,eqs,bta,del,omb,rsq,rsq0,u,v
    integer :: iz,jz,j,k,l,m,me,mnl,nin

    ! Magical tolerances. Will play with them to see what happens.
    real(flyt), parameter :: sml=1.0E-12_flyt ! default 1.0E-5_flyt
    real(flyt), parameter :: eps=1.0E-13_flyt ! default 1.0E-6_flyt
    real(flyt), parameter :: big=9.9E35_flyt
    integer, parameter :: mnlam=5
    real(flyt), parameter :: rsqmax=0.999999_flyt
    !real(flyt), parameter :: pmin=1E-12_flyt
    !real(flyt), parameter :: exmx=250.0_flyt
    ! Tolerances. Kinda tricky to figure out but will try.
    lmu=0
    nlp=0
    jerr=0

    allocate(c(eh%nx,eh%nx))
    allocate(a(eh%nx))
    allocate(mm(eh%nx))
    allocate(da(eh%nx))
    allocate(ia(eh%nx))
    ia=0
    c=0.0_flyt
    a=0.0_flyt
    mm=0
    da=0.0_flyt
    bta=eh%alpha
    omb=1.0_flyt-bta
    if ( eh%autolambda ) then
        eqs=max(eps,eh%smallest_lambda)
        alf=eqs**(1.0_flyt/real(eh%nlambda-1,flyt))
    endif
    rsq=0.0_flyt
    a=0.0_flyt
    mm=0
    nlp=0
    nin=nlp
    iz=0
    mnl=min(mnlam,eh%nlambda)
    alm=0.0_flyt

    ! Now major loop
    lambdaloop: do m=1,eh%nlambda

        ! Decide on lambda value
        if ( eh%autolambda ) then
            if ( m .gt. 2 ) then
                alm=alm*alf
            elseif ( m .eq. 2 ) then
                alm=0.0_flyt
                do j=1,eh%nx
                    alm=max(alm,abs( ATB(j) ))
                enddo
                alm=alf*alm/max(bta,1.0E-3_flyt)
            elseif ( m .eq. 1 ) then
                alm=big
            endif
        else
            alm=eh%input_lambda(m)
        endif

        dem=alm*omb
        ab=alm*bta
        rsq0=rsq
        jz=1

        inloop: do
            ! So something for some reason. Maybe the first iteration or something.
            if( iz*jz .eq. 0 ) then
                nlp=nlp+1
                dlx=0.0_flyt
                kloop: do k=1,eh%nx
                    ak=a(k)
                    u=ATB(k)+ak
                    v=abs(u)-ab
                    a(k)=0.0_flyt
                    if( v .gt. 0.0_flyt ) then
                        a(k)=sign(v,u)/(1.0_flyt+dem)
                    endif
                    if( a(k) .eq. ak ) cycle kloop
                    if ( mm(k) .eq. 0 ) then
                        nin=nin+1
                        if( nin .gt. eh%nx ) exit kloop
                        jl2: do j=1,eh%nx
                            if( mm(j) .ne. 0 ) then
                                c(j,nin)=c(k,mm(j))
                                cycle jl2
                            endif
                            if ( j .eq. k ) then
                                c(j,nin)=1.0_flyt
                                cycle jl2
                            endif
                            c(j,nin)=ATA(j,k)
                        enddo jl2
                        mm(k)=nin
                        ia(nin)=k
                    endif
                    del=a(k)-ak
                    rsq=rsq+del*(2.0_flyt*ATB(k)-del)
                    dlx=max(del**2,dlx)
                    do j=1,eh%nx
                        ATB(j)=ATB(j)-c(j,mm(k))*del
                    enddo
                enddo kloop
                if( dlx .lt. eh%threshold ) exit inloop
                if( nin .gt. eh%nx ) exit inloop
                if( nlp .gt. eh%maxit ) then
                    jerr=-m
                    return
                endif
            endif
            iz=1
            da(1:nin)=a(ia(1:nin))
            ! Then do more things, for possibly other reasons
            ininloop: do
                nlp=nlp+1
                dlx=0.0_flyt
                ll1: do l=1,nin
                    k=ia(l)
                    ak=a(k)
                    u=ATB(k)+ak
                    v=abs(u)-ab
                    a(k)=0.0_flyt
                    if( v .gt. 0.0_flyt ) then
                        a(k)=sign(v,u)/(1.0_flyt+dem)
                    endif
                    ! replace with tolerance?
                    if ( a(k) .eq. ak ) cycle ll1
                    !if ( abs(a(k)-ak) .lt. 0.0_flyt ) cycle ll1
                    del=a(k)-ak
                    rsq=rsq+del*(2.0_flyt*ATB(k)-del)
                    dlx=max(del**2,dlx)
                    do j=1,nin
                        ATB(ia(j))=ATB(ia(j))-c(ia(j),mm(k))*del
                    enddo
                enddo ll1

                if( dlx .lt. eh%threshold ) then
                    exit ininloop
                endif
                if( nlp .gt. eh%maxit ) then
                    jerr=-m
                    return
                endif
            enddo ininloop

            da(1:nin)=a(ia(1:nin))-da(1:nin)
            jl3: do j=1,eh%nx
                if( mm(j) .ne. 0 ) cycle jl3
                ATB(j)=ATB(j)-dot_product(da(1:nin),c(j,1:nin))
            enddo jl3
            jz=0
        enddo inloop

        if ( nin .gt. eh%nx ) then
            jerr=-10000-m
            exit lambdaloop
        endif
        ! This stores the output, I think
        if( nin .gt. 0 ) eh%solutions(1:nin,m)=a(ia(1:nin))
        eh%nvar(m)=nin
        eh%rsq(m)=rsq
        eh%lambdavals(m)=alm
        eh%varorder=ia
        lmu=m
        if( m .lt. mnl ) cycle lambdaloop
        if ( eh%autolambda .eqv. .false. ) cycle lambdaloop
        me=0
        do j=1,nin
            if( eh%solutions(j,m) .ne. 0.0_flyt ) me=me+1
        enddo
        if( rsq-rsq0 .lt. sml*rsq ) exit lambdaloop
        if( rsq .gt. rsqmax ) exit lambdaloop
    enddo lambdaloop
    deallocate(a,mm,c,da)
end subroutine

end module
