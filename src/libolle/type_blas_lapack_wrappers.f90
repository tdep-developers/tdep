#include "precompilerdefinitions"
module type_blas_lapack_wrappers
    ! Not complete by any means. I will add stuff when I need it.
    use konstanter, only: flyt
    implicit none

    !> Matrix-matrix multiplication
    !> C = alpha A*B + beta*C, where A or B can be transposed
    interface lo_gemm
        module procedure :: lo_dgemm,lo_zgemm
    end interface

    !> Matrix-vector multiplication
    !> Y = alpha*A*Y + beta*Y
    interface lo_gemv
        module procedure :: lo_dgemv,lo_zgemv
    end interface
contains

subroutine lo_dgemm(a,b,c,transa,transb,alpha,beta)
    real(flyt), intent(in) :: a(:,:)
    real(flyt), intent(in) :: b(:,:)
    real(flyt), intent(inout) :: c(:,:)
    character(len=1), intent(in), optional :: transa
    character(len=1), intent(in), optional :: transb
    real(flyt), intent(in), optional :: alpha
    real(flyt), intent(in), optional :: beta
    !
    character(len=1) :: o_transa,o_transb
    real(flyt) :: o_alpha,o_beta
    integer :: m,n,k,lda,ldb,ldc
    !
    if(present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1
    endif
    if(present(beta)) then
        o_beta = beta
    else
        o_beta = 0
    endif
    if(present(transa)) then
        o_transa = transa
    else
        o_transa = 'n'
    endif
    if(present(transb)) then
        o_transb = transb
    else
        o_transb = 'n'
    endif
    if ((o_transa .eq. 'n' .or. o_transa .eq. 'N')) then
        k = size(a,2)
    else
        k = size(a,1)
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call dgemm(o_transa,o_transb,m,n,k,o_alpha,a,lda,b,ldb,o_beta,c,ldc)
end subroutine

subroutine lo_zgemm(a,b,c,transa,transb,alpha,beta)
    complex(flyt), intent(in) :: a(:,:)
    complex(flyt), intent(in) :: b(:,:)
    complex(flyt), intent(inout) :: c(:,:)
    character(len=1), intent(in), optional :: transa
    character(len=1), intent(in), optional :: transb
    complex(flyt), intent(in), optional :: alpha
    complex(flyt), intent(in), optional :: beta
    !
    character(len=1) :: o_transa,o_transb
    complex(flyt) :: o_alpha,o_beta
    integer :: m,n,k,lda,ldb,ldc
    !
    if(present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1
    endif
    if(present(beta)) then
        o_beta = beta
    else
        o_beta = 0
    endif
    if(present(transa)) then
        o_transa = transa
    else
        o_transa = 'n'
    endif
    if(present(transb)) then
        o_transb = transb
    else
        o_transb = 'n'
    endif
    if((o_transa.eq.'n'.or.o_transa.eq.'N')) then
        k = size(a,2)
    else
        k = size(a,1)
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call zgemm(o_transa,o_transb,m,n,k,o_alpha,a,lda,b,ldb,o_beta,c,ldc)
end subroutine

subroutine lo_dgemv(a,x,y,alpha,beta,trans)
    real(flyt), intent(in) :: a(:,:)
    real(flyt), intent(in) :: x(:)
    real(flyt), intent(inout) :: y(:)
    real(flyt), intent(in), optional :: alpha
    real(flyt), intent(in), optional :: beta
    character(len=1), intent(in), optional :: trans
    !
    real(flyt) :: o_alpha,o_beta
    character(len=1) :: o_trans
    integer :: incx,incy,m,n,lda
    !
    if(present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1
    endif
    if(present(beta)) then
        o_beta = beta
    else
        o_beta = 0
    endif
    if(present(trans)) then
        o_trans = trans
    else
        o_trans = 'N'
    endif
    incx = 1
    incy = 1
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    call dgemv(o_trans,m,n,o_alpha,a,lda,x,incx,o_beta,y,incy)
end subroutine

subroutine lo_zgemv(a,x,y,alpha,beta,trans)
    complex(flyt), intent(in) :: a(:,:)
    complex(flyt), intent(in) :: x(:)
    complex(flyt), intent(inout) :: y(:)
    complex(flyt), intent(in), optional :: alpha
    complex(flyt), intent(in), optional :: beta
    character(len=1), intent(in), optional :: trans
    !
    complex(flyt) :: o_alpha,o_beta
    character(len=1) :: o_trans
    integer :: incx,incy,m,n,lda
    !
    if(present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1
    endif
    if(present(beta)) then
        o_beta = beta
    else
        o_beta = 0
    endif
    if(present(trans)) then
        o_trans = trans
    else
        o_trans = 'N'
    endif
    incx = 1
    incy = 1
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    call zgemv(o_trans,m,n,o_alpha,a,lda,x,incx,o_beta,y,incy)
end subroutine

subroutine lo_zgeev(a,w,vl,vr,info)
    complex(flyt), intent(inout) :: a(:,:)
    complex(flyt), intent(out) :: w(:)
    complex(flyt), intent(out), optional, target :: vl(:,:)
    complex(flyt), intent(out), optional, target :: vr(:,:)
    integer, intent(out), optional :: info
    !
    character(len=4), parameter :: srname = 'GEEV'
    character(len=1) :: jobvl,jobvr
    integer :: n,lda,ldvl,ldvr,lwork,l_stat_alloc,l_stat_dealloc,o_info
    complex(flyt), pointer :: o_vl(:,:),o_vr(:,:),work(:)
    real(flyt), pointer :: rwork(:)
    complex(flyt) :: s_work(1)
    complex(flyt), target :: l_a2_comp(1,1)
    !
    if(present(vl)) then
        jobvl = 'V'
    else
        jobvl = 'N'
    endif
    if(present(vr)) then
        jobvr = 'V'
    else
        jobvr = 'N'
    endif
    lda = max(1,size(a,1))
    if(present(vl)) then
        ldvl = max(1,size(vl,1))
    else
        ldvl = 1
    endif
    if(present(vr)) then
        ldvr = max(1,size(vr,1))
    else
        ldvr = 1
    endif
    n = size(a,2)
    l_stat_alloc = 0
    if(present(vl)) then
        o_vl => vl
    else
        o_vl => l_a2_comp
    endif
    if(present(vr)) then
        o_vr => vr
    else
        o_vr => l_a2_comp
    endif
    allocate(rwork(2*n), stat=l_stat_alloc)
    lwork = -1
    call zgeev(jobvl,jobvr,n,a,lda,w,o_vl,ldvl,o_vr,ldvr,s_work,lwork,rwork,o_info)

    if ( o_info == 0 ) then
        !lwork = s_work(1)
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call zgeev(jobvl,jobvr,n,a,lda,w,o_vl,ldvl,o_vr,ldvr,work,lwork,rwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    deallocate(rwork, stat=l_stat_dealloc)
    ! error
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zheev(a,w,jobz,uplo,info)
    complex(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: w(:)
    character(len=1), intent(in), optional :: jobz
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info
    !
    character(len=4), parameter :: srname = 'HEEV'
    character(len=1) :: o_jobz
    character(len=1) :: o_uplo
    integer :: o_info,n,lda,lwork,l_stat_alloc,l_stat_dealloc
    complex(flyt), pointer :: work(:)
    real(flyt), pointer :: rwork(:)
    complex(flyt) :: s_work(1)
    !
    if(present(jobz)) then
        o_jobz = jobz
    else
        o_jobz = 'N'
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    lda = max(1,size(a,1))
    n = size(a,2)
    l_stat_alloc = 0
    allocate(rwork(max(1,3*n-2)), stat=l_stat_alloc)
    lwork = -1
    call zheev(o_jobz,o_uplo,n,a,lda,w,s_work,lwork,rwork,o_info)
    !
    if ( o_info == 0 ) then
        !lwork = s_work(1)
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        ! call lapack routine
        if(l_stat_alloc==0) then
            call zheev(o_jobz,o_uplo,n,a,lda,w,work,lwork,rwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    deallocate(rwork, stat=l_stat_dealloc)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dgels(a,b,trans,info)
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(inout) :: b(:,:)
    character(len=1), intent(in), optional :: trans
    integer, intent(out), optional :: info
    !
    character(len=4), parameter :: srname = 'GELS'
    character(len=1) :: o_trans
    integer :: o_info,m,n,nrhs,lda,ldb,lwork,l_stat_alloc,l_stat_dealloc
    real(flyt), pointer :: work(:)
    real(flyt) :: s_work(1)
    !
    if(present(trans)) then
        o_trans = trans
    else
        o_trans = 'N'
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    m = size(a,1)
    n = size(a,2)
    nrhs = size(b,2)
    l_stat_alloc = 0
    lwork = -1
    call dgels(o_trans,m,n,nrhs,a,lda,b,ldb,s_work,lwork,o_info)
    if ( o_info == 0 ) then
        lwork = int(anint(s_work(1)))
        allocate(work(lwork), stat=l_stat_alloc)
        if(l_stat_alloc==0) then
            call dgels(o_trans,m,n,nrhs,a,lda,b,ldb,work,lwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dgelss(a,b,rank,s,rcond,info)
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(inout) :: b(:,:)
    real(flyt), intent(out), optional, target :: s(:)
    integer, intent(out), optional :: rank
    real(flyt), intent(in), optional :: rcond
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'GELSS'
    integer :: o_rank,o_info,m,n,nrhs,lda,ldb,lwork,l_stat_alloc,l_stat_dealloc
    real(flyt) :: o_rcond
    real(flyt), pointer :: o_s(:)
    real(flyt), pointer :: work(:)
    real(flyt) :: s_work(1)
    !
    if(present(rcond)) then
        o_rcond = rcond
    else
        o_rcond = 100*epsilon(1.0_flyt)
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    m = size(a,1)
    n = size(a,2)
    nrhs = size(b,2)
    l_stat_alloc = 0
    if(present(s)) then
        o_s => s
    else
        allocate(o_s(min(m,n)), stat=l_stat_alloc)
    endif
    lwork = -1
    call dgelss(m,n,nrhs,a,lda,b,ldb,o_s,o_rcond,o_rank,s_work,lwork,o_info)
    if ( o_info == 0 ) then
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call dgelss(m,n,nrhs,a,lda,b,ldb,o_s,o_rcond,o_rank,work,lwork,o_info)
        else; o_info = -1000
        endif
        if(present(rank)) then
            rank = o_rank
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(.not. present(s)) then
        deallocate(o_s, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zgelss(a,b,rank,s,rcond,info)
    complex(flyt), intent(inout) :: a(:,:)
    complex(flyt), intent(inout) :: b(:,:)
    real(flyt), intent(out), optional, target :: s(:)
    integer, intent(out), optional :: rank
    real(flyt), intent(in), optional :: rcond
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'GELSS'
    complex(flyt), pointer :: work(:)
    complex(flyt) :: s_work(1)
    real(flyt), pointer :: o_s(:)
    real(flyt) :: o_rcond
    integer :: o_rank
    integer :: o_info,m,n,nrhs,lda,ldb,lwork, l_stat_alloc, l_stat_dealloc
    real(flyt), pointer :: rwork(:)
    !
    if(present(rcond)) then
        o_rcond = rcond
    else
        o_rcond = 100*epsilon(1.0_flyt)
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    m = size(a,1)
    n = size(a,2)
    nrhs = size(b,2)
    l_stat_alloc = 0
    if(present(s)) then
        o_s => s
    else
        allocate(o_s(min(m,n)), stat=l_stat_alloc)
    endif
    if(l_stat_alloc==0) then
        allocate(rwork(5*min(m,n)), stat=l_stat_alloc)
    endif
    lwork = -1
    call zgelss(m,n,nrhs,a,lda,b,ldb,o_s,o_rcond,o_rank,s_work,lwork,rwork,o_info)
    if ( o_info == 0 ) then
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call zgelss(m,n,nrhs,a,lda,b,ldb,o_s,o_rcond,o_rank,work,lwork,rwork,o_info)
        else
            o_info = -1000
        endif
        if(present(rank)) then
            rank = o_rank
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(.not. present(s)) then
        deallocate(o_s, stat=l_stat_dealloc)
    endif
    deallocate(rwork, stat=l_stat_dealloc)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dsyev(a,w,jobz,uplo,info)
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: w(:)
    character(len=1), intent(in), optional :: jobz
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info
    !
    character(len=4), parameter :: srname = 'SYEV'
    character(len=1) :: o_jobz,o_uplo
    integer :: o_info,n,lda,lwork,l_stat_alloc,l_stat_dealloc
    real(flyt), pointer :: work(:)
    real(flyt) :: s_work(1)

    if(present(jobz)) then
        o_jobz = jobz
    else
        o_jobz = 'N'
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    lda = max(1,size(a,1))
    n = size(a,2)
    l_stat_alloc = 0
    lwork = -1
    call dsyev(o_jobz,o_uplo,n,a,lda,w,s_work,lwork,o_info)
    if ( o_info .eq. 0 ) then
        lwork = int(s_work(1))
        allocate(work(lwork), stat=l_stat_alloc)
        if(l_stat_alloc==0) then
            call dsyev(o_jobz,o_uplo,n,a,lda,w,work,lwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dgesvd(a,s,u,vt,ww,job,info)
    character(len=1), intent(in), optional :: job
    integer, intent(out), optional :: info
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: s(:)
    real(flyt), intent(out), optional, target :: u(:,:)
    real(flyt), intent(out), optional, target :: vt(:,:)
    real(flyt), intent(out), optional :: ww(:)
    !
    character(len=5), parameter :: srname = 'GESVD'
    character(len=1) :: o_job
    integer :: o_info
    character(len=1) :: jobu
    character(len=1) :: jobvt
    integer :: m,n,lda,ldu,ldvt,lwork,l_stat_alloc,l_stat_dealloc
    real(flyt), pointer :: o_u(:,:)
    real(flyt), pointer :: o_vt(:,:)
    real(flyt), pointer :: work(:)
    real(flyt) :: s_work(1)
    real(flyt), target :: l_a2_real(1,1)
    !
    if(present(job)) then
        o_job = job
    else
        o_job = 'N'
    endif
    lda = max(1,size(a,1))
    if(present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    endif
    if(present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    endif
    m = size(a,1)
    n = size(a,2)
    if(present(u)) then
        if(size(u,2)==m) then
            jobu = 'A'
        else
            jobu = 'S'
        endif
    else
        if((o_job.eq.'U'.or.o_job.eq.'u')) then
            jobu = 'O'
        else
            jobu = 'N'
        endif
    endif
    if(present(vt)) then
        if(size(vt,1)==n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        endif
    else
        if((o_job.eq.'V'.or.o_job.eq.'v')) then
            jobvt = 'O'
        else
            jobvt = 'N'
        endif
    endif
    l_stat_alloc = 0
    if(present(u)) then
        o_u => u
    else
        o_u => l_a2_real
    endif
    if(present(vt)) then
        o_vt => vt
    else
        o_vt => l_a2_real
    endif
    lwork = -1
    call dgesvd(jobu,jobvt,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,s_work,lwork,o_info)
    if ( o_info == 0 ) then
        !lwork = s_work(1)
        lwork = int(s_work(1))
        allocate(work(lwork), stat=l_stat_alloc)
        if(l_stat_alloc==0) then
            call dgesvd(jobu,jobvt,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,work,lwork,o_info)
        else
            o_info = -1000
        endif
        if(present(ww)) then
            ww = work(2:min(m,n))
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dgesdd(a,s,u,vt,jobz,info)
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: s(:)
    real(flyt), intent(out), optional, target :: u(:,:)
    real(flyt), intent(out), optional, target :: vt(:,:)
    character(len=1), intent(in), optional :: jobz
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'GESDD'
    character(len=1) :: o_jobz
    integer :: o_info,m,n,lda,ldu,ldvt,lwork,l_stat_alloc,l_stat_dealloc
    real(flyt), pointer :: o_u(:,:)
    real(flyt), pointer :: o_vt(:,:)
    real(flyt), pointer :: work(:)
    integer, pointer :: iwork(:)
    real(flyt) :: s_work(1)
    real(flyt), target :: l_a2_real(1,1)
    !
    if(present(jobz)) then
        o_jobz = jobz
    else
        o_jobz = 'N'
    endif
    lda = max(1,size(a,1))
    if(present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    endif
    if(present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    endif
    m = size(a,1)
    n = size(a,2)
    l_stat_alloc = 0
    if(present(u)) then
        o_u => u
    else
        o_u => l_a2_real
    endif
    if(present(vt)) then
        o_vt => vt
    else
        o_vt => l_a2_real
    endif
    allocate(iwork(8*min(m,n)), stat=l_stat_alloc)
    lwork = -1
    call dgesdd(o_jobz,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,s_work,lwork,iwork,o_info)
    if ( o_info == 0 ) then
        !lwork = s_work(1)
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call dgesdd(o_jobz,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,work,lwork,iwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    deallocate(iwork, stat=l_stat_dealloc)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dgglse(a,b,c,d,x,info)
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(inout) :: b(:,:)
    real(flyt), intent(inout) :: c(:)
    real(flyt), intent(inout) :: d(:)
    real(flyt), intent(out) :: x(:)
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'GGLSE'
    integer :: o_info,m,n,p,lda,ldb,lwork,l_stat_alloc,l_stat_dealloc
    real(flyt), pointer :: work(:)
    real(flyt) :: s_work(1)
    !
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    m = size(a,1)
    n = size(a,2)
    p = size(b,1)
    l_stat_alloc = 0
    lwork = -1
    call dgglse(m,n,p,a,lda,b,ldb,c,d,x,s_work,lwork,o_info)
    if ( o_info == 0 ) then
        !lwork = s_work(1)
        lwork = int(s_work(1))
        allocate(work(lwork), stat=l_stat_alloc)
        if(l_stat_alloc==0) then
            call dgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zgglse(a,b,c,d,x,info)
    complex(flyt), intent(inout) :: a(:,:)
    complex(flyt), intent(inout) :: b(:,:)
    complex(flyt), intent(inout) :: c(:)
    complex(flyt), intent(inout) :: d(:)
    complex(flyt), intent(out) :: x(:)
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'GGLSE'
    integer :: o_info
    integer :: m,n,p,lda,ldb,lwork,l_stat_alloc,l_stat_dealloc
    complex(flyt), pointer :: work(:)
    complex(flyt) :: s_work(1)
    !
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    m = size(a,1)
    n = size(a,2)
    p = size(b,1)
    l_stat_alloc = 0
    lwork = -1
    call zgglse(m,n,p,a,lda,b,ldb,c,d,x,s_work,lwork,o_info)
    if(o_info == 0) then
        lwork = int(s_work(1))
        allocate(work(lwork), stat=l_stat_alloc)
        if(l_stat_alloc==0) then
            call zgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dsyevd(a,w,jobz,uplo,info)
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: w(:)
    character(len=1), intent(in), optional :: jobz
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'SYEVD'
    character(len=1) :: o_jobz,o_uplo
    integer :: o_info,n,lda,lwork,liwork,l_stat_alloc,l_stat_dealloc
    real(flyt), pointer :: work(:)
    integer, pointer :: iwork(:)
    integer :: s_iwork(1)
    real(flyt) :: s_work(1)
    !
    if(present(jobz)) then
        o_jobz = jobz
    else
        o_jobz = 'N'
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    lda = max(1,size(a,1))
    n = size(a,2)
    l_stat_alloc = 0
    liwork = -1
    lwork = -1
    call dsyevd(o_jobz,o_uplo,n,a,lda,w,s_work,lwork,s_iwork,liwork,o_info)
    if ( o_info == 0 ) then
        !liwork = s_iwork(1)
        !lwork = s_work(1)
        liwork = int(s_iwork(1))
        lwork = int(s_work(1))
        allocate(iwork(liwork), stat=l_stat_alloc)
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call dsyevd(o_jobz,o_uplo,n,a,lda,w,work,lwork,iwork,liwork, o_info)
        else
            o_info = -1000
        endif
        deallocate(iwork, stat=l_stat_dealloc)
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zheevr(a,w,uplo,z,vl,vu,il,iu,m,isuppz,abstol,info)
    complex(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: w(:)
    complex(flyt), intent(out), optional, target :: z(:,:)
    integer, intent(out), optional, target :: isuppz(:)
    character(len=1), intent(in), optional :: uplo
    real(flyt), intent(in), optional :: vl
    real(flyt), intent(in), optional :: vu
    integer, intent(in), optional :: il
    integer, intent(in), optional :: iu
    integer, intent(out), optional :: m
    real(flyt), intent(in), optional :: abstol
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'HEEVR'
    character(len=1) :: o_uplo,jobz,range
    real(flyt) :: o_vl,o_vu,o_abstol
    integer :: o_il,o_iu,o_m,o_info,n,lda,ldz,lwork,lrwork,liwork,l_stat_alloc,l_stat_dealloc
    complex(flyt), pointer :: o_z(:,:)
    !integer, pointer :: o_isuppz(:)
    integer, allocatable :: o_isuppz(:)
    complex(flyt), pointer :: work(:)
    real(flyt), pointer :: rwork(:)
    integer, pointer :: iwork(:)
    integer :: s_iwork(1)
    real(flyt) :: s_rwork(1)
    complex(flyt) :: s_work(1)
    integer, target :: l_a1_inte(1)
    complex(flyt), target :: l_a2_comp(1,1)
    !
    o_info=0
    if(present(abstol)) then
        o_abstol = abstol
    else
        o_abstol = 0.0_flyt
    endif
    if(present(il)) then
        o_il = il
    else
        o_il = 1
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    if(present(vl)) then
        o_vl = vl
    else
        o_vl = -huge(vl)
    endif
    if(present(vu)) then
        o_vu = vu
    else
        o_vu = huge(vl)
    endif
    if(present(z)) then
        jobz = 'V'
    else
        jobz = 'N'
    endif
    lda = max(1,size(a,1))
    if(present(z)) then
        ldz = max(1,size(z,1))
    else
        ldz = 1
    endif
    n = size(a,2)
    if((present(vl).or.present(vu)).and.(present(il).or.present(iu))) then
        o_info=-1001
    elseif((present(vl).or.present(vu))) then
        range = 'V'
    elseif((present(il).or.present(iu))) then
        range = 'I'
    else
        range = 'A'
    endif
    if(present(iu)) then
        o_iu = iu
    else
        o_iu = n
    endif
    !
    l_stat_alloc = 0
    if(.not.present(z)) then
        if(present(isuppz)) then
            o_info=-1001
        else
            allocate(o_isuppz(1))
            o_isuppz=l_a1_inte
        endif
    else
        if(present(isuppz)) then
            allocate(o_isuppz(1))
            o_isuppz=l_a1_inte
        else
            allocate(o_isuppz(2*n), stat=l_stat_alloc)
        endif
    endif
    if(present(z)) then
        o_z => z
    else
        o_z => l_a2_comp
    endif
    !
    liwork = -1
    lrwork = -1
    lwork = -1
    if ( o_info == 0 ) then
        call zheevr(jobz,range,o_uplo,n,a,lda,o_vl,o_vu,o_il,o_iu,o_abstol,o_m,w,o_z,ldz,o_isuppz,s_work,lwork,s_rwork,lrwork,s_iwork,liwork,o_info)
    endif
    if (o_info == 0 ) then
        !liwork = s_iwork(1)
        !lrwork = s_rwork(1)
        !lwork = s_work(1)
        liwork = int(s_iwork(1))
        lrwork = int(s_rwork(1))
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(iwork(liwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            allocate(rwork(lrwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call zheevr(jobz,range,o_uplo,n,a,lda,o_vl,o_vu,o_il,o_iu,o_abstol,o_m,w,o_z,ldz,o_isuppz,work,lwork,rwork,lrwork,iwork,liwork,o_info)
        else
            o_info = -1000
        endif
        if(present(m)) then
            m = o_m
        endif
        deallocate(iwork, stat=l_stat_dealloc)
        deallocate(rwork, stat=l_stat_dealloc)
        deallocate(work, stat=l_stat_dealloc)
    endif
    if ( allocated(o_isuppz) ) deallocate(o_isuppz, stat=l_stat_dealloc)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zgesvd(a,s,u,vt,ww,job,info)
    complex(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: s(:)
    complex(flyt), intent(out), optional, target :: u(:,:)
    complex(flyt), intent(out), optional, target :: vt(:,:)
    real(flyt), intent(out), optional :: ww(:)
    character(len=1), intent(in), optional :: job
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'GESVD'
    character(len=1) :: o_job
    integer :: o_info
    character(len=1) :: jobu
    character(len=1) :: jobvt
    integer :: m,n,lda,ldu,ldvt,lwork,l_stat_alloc,l_stat_dealloc
    complex(flyt), pointer :: o_u(:,:)
    complex(flyt), pointer :: o_vt(:,:)
    complex(flyt), pointer :: work(:)
    real(flyt), pointer :: rwork(:)
    complex(flyt) :: s_work(1)
    complex(flyt), target :: l_a2_comp(1,1)
    !
    if(present(job)) then
        o_job = job
    else
        o_job = 'N'
    endif
    lda = max(1,size(a,1))
    if(present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    endif
    if(present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    endif
    m = size(a,1)
    n = size(a,2)
    if(present(u)) then
        if(size(u,2)==m) then
            jobu = 'A'
        else
            jobu = 'S'
        endif
    else
        if((o_job.eq.'U'.or.o_job.eq.'u')) then
            jobu = 'O'
        else
            jobu = 'N'
        endif
    endif
    if(present(vt)) then
        if(size(vt,1)==n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        endif
    else
        if((o_job.eq.'V'.or.o_job.eq.'v')) then
            jobvt = 'O'
        else
            jobvt = 'N'
        endif
    endif
    l_stat_alloc = 0
    if(present(u)) then
        o_u => u
    else
        o_u => l_a2_comp
    endif
    if(present(vt)) then
        o_vt => vt
    else
        o_vt => l_a2_comp
    endif
    allocate(rwork(5*min(m,n)), stat=l_stat_alloc)
    lwork = -1
    call zgesvd(jobu,jobvt,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,s_work,lwork,rwork,o_info)
    if ( o_info == 0 ) then
        !lwork = s_work(1)
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call zgesvd(jobu,jobvt,m,n,a,lda,s,o_u,ldu,o_vt,ldvt,work,lwork,rwork,o_info)
        else
            o_info = -1000
        endif
        if(present(ww)) then
            ww = rwork(1:min(m,n)-1)
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    deallocate(rwork, stat=l_stat_dealloc)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zgels(a,b,trans,info)
    complex(flyt), intent(inout) :: a(:,:)
    complex(flyt), intent(inout) :: b(:,:)
    character(len=1), intent(in), optional :: trans
    integer, intent(out), optional :: info
    !
    character(len=4), parameter :: srname = 'GELS'
    character(len=1) :: o_trans
    integer :: o_info
    integer :: m
    integer :: n
    integer :: nrhs
    integer :: lda
    integer :: ldb
    integer :: lwork
    integer :: l_stat_alloc, l_stat_dealloc
    complex(flyt), pointer :: work(:)
    complex(flyt) :: s_work(1)
    !
    if(present(trans)) then
        o_trans = trans
    else
        o_trans = 'N'
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    m = size(a,1)
    n = size(a,2)
    nrhs = size(b,2)
    l_stat_alloc = 0
    lwork = -1
    call zgels(o_trans,m,n,nrhs,a,lda,b,ldb,s_work,lwork,o_info)
    if ( o_info == 0 ) then
        !lwork = s_work(1)
        lwork = int(s_work(1))
        allocate(work(lwork), stat=l_stat_alloc)
        if(l_stat_alloc==0) then
            call zgels(o_trans,m,n,nrhs,a,lda,b,ldb,work,lwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dgesv(a,b,ipiv,info)
    integer, intent(out), optional :: info
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(inout) :: b(:,:)
    integer, intent(out), optional, target :: ipiv(:)
    !
    character(len=4), parameter :: srname = 'GESV'
    integer :: o_info, n, nrhs, lda, ldb, l_stat_alloc, l_stat_dealloc
    integer, pointer :: o_ipiv(:)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    l_stat_alloc = 0
    if(present(ipiv)) then
        o_ipiv => ipiv
    else
        allocate(o_ipiv(n), stat=l_stat_alloc)
    endif
    if(l_stat_alloc==0) then
        call dgesv(n,nrhs,a,lda,o_ipiv,b,ldb,o_info)
    else
        o_info = -1000
    endif
    if( .not.present(ipiv)) then
        deallocate(o_ipiv, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dsyevr(a,w,uplo,z,vl,vu,il,iu,m,isuppz,abstol,info)
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(out) :: w(:)
    character(len=1), intent(in), optional :: uplo
    real(flyt), intent(out), optional, target :: z(:,:)
    real(flyt), intent(in), optional :: vl
    real(flyt), intent(in), optional :: vu
    integer, intent(in), optional :: il
    integer, intent(in), optional :: iu
    integer, intent(out), optional :: m
    integer, intent(out), optional, target :: isuppz(:)
    real(flyt), intent(in), optional :: abstol
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'SYEVR'
    character(len=1) :: o_uplo,jobz,range
    real(flyt), pointer :: o_z(:,:),work(:)
    real(flyt), target :: l_a2_real(1,1)
    real(flyt) :: s_work(1),o_vl,o_vu,o_abstol
    integer, pointer :: o_isuppz(:),iwork(:)
    integer :: o_il,o_iu,o_m,o_info,n,lda,ldz,lwork,liwork,l_stat_alloc, l_stat_dealloc,s_iwork(1)
    !
    if(present(abstol)) then
        o_abstol = abstol
    else
        o_abstol = 0.0_flyt
    endif
    if(present(il)) then
        o_il = il
    else
        o_il = 1
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    if(present(vl)) then
        o_vl = vl
    else
        o_vl = -huge(vl)
    endif
    if(present(vu)) then
        o_vu = vu
    else
        o_vu = huge(vl)
    endif
    if(present(z)) then
        jobz = 'V'
    else
        jobz = 'V'
    endif
    lda = max(1,size(a,1))
    if(present(z)) then
        ldz = max(1,size(z,1))
    else
        ldz = 1
    endif
    n = size(a,2)
    if ( (present(vl).or.present(vu)) .and. (present(il).or.present(iu)) ) then
        o_info=-1001
        ! send to error
        if(present(info)) then
            info = o_info
        elseif(o_info <= -1000) then
            call xerbla(srname,-o_info)
        endif
        return
    elseif((present(vl).or.present(vu))) then
        range = 'V'
    elseif((present(il).or.present(iu))) then
        range = 'I'
    else
        range = 'A'
    endif
    if(present(iu)) then
        o_iu = iu
    else
        o_iu = n
    endif

    l_stat_alloc = 0
    if(present(isuppz)) then
        o_isuppz => isuppz
    else
        allocate(o_isuppz(2*n), stat=l_stat_alloc)
    endif
    if(present(z)) then
        o_z => z
    else
        o_z => l_a2_real
    endif
    liwork = -1
    lwork = -1
    call dsyevr(jobz,range,o_uplo,n,a,lda,o_vl,o_vu,o_il,o_iu,o_abstol,o_m,w,o_z,ldz,o_isuppz,s_work,lwork,s_iwork,liwork,o_info)
    if ( o_info .eq. 0 ) then
        liwork = int(s_iwork(1))
        lwork = int(s_work(1))
        if(l_stat_alloc==0) then
            allocate(iwork(liwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call dsyevr(jobz,range,o_uplo,n,a,lda,o_vl,o_vu,o_il,o_iu,o_abstol,o_m,w,o_z,ldz,o_isuppz,work,lwork,iwork,liwork,o_info)
        else
            o_info = -1000
        endif
        if(present(m)) then
            m = o_m
        endif
        deallocate(iwork, stat=l_stat_dealloc)
        deallocate(work, stat=l_stat_dealloc)
    endif

    if(.not. present(isuppz)) then
        deallocate(o_isuppz, stat=l_stat_dealloc)
    endif
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dsygst(a,b,itype,uplo,info)
    real(flyt), dimension(:,:), intent(inout) :: a
    real(flyt), dimension(:,:), intent(in) :: b
    integer, intent(in), optional :: itype
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'SYGST'
    character(len=1) :: o_uplo
    integer :: o_itype,o_info,n,lda,ldb

    if(present(itype)) then
        o_itype = itype
    else
        o_itype = 1
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    call dsygst(o_itype,o_uplo,n,a,lda,b,ldb,o_info)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zhegst(a,b,itype,uplo,info)
    complex(flyt), intent(inout) :: a(:,:)
    complex(flyt), intent(in) :: b(:,:)
    integer, intent(in), optional :: itype
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'HEGST'
    character(len=1) :: o_uplo
    integer :: o_itype,o_info,n,lda,ldb
    if(present(itype)) then
        o_itype = itype
    else
        o_itype = 1
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'u'
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    call zhegst(o_itype,o_uplo,n,a,lda,b,ldb,o_info)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zpotrf(a,uplo,info)
    complex(flyt), intent(inout) :: a(:,:)
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info
    !
    character(len=5), parameter :: srname = 'POTRF'
    character(len=1) :: o_uplo
    integer :: o_info,n,lda
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    lda = max(1,size(a,1))
    n = size(a,2)
    call zpotrf(o_uplo,n,a,lda,o_info)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zhegv(a,b,w,itype,jobz,uplo,info)
    complex(flyt), dimension(:,:), intent(inout) :: a
    complex(flyt), dimension(:,:), intent(inout) :: b
    real(flyt), dimension(:), intent(out) :: w
    integer, intent(in), optional :: itype
    character(len=1), intent(in), optional :: jobz
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info

    character(len=4), parameter :: srname = 'HEGV'
    character(len=1) :: o_jobz,o_uplo
    integer :: o_itype,o_info,n,lda,ldb,lwork,l_stat_alloc
    complex(flyt), pointer :: work(:)
    real(flyt), pointer :: rwork(:)
    complex(flyt) :: s_work(1)
    if(present(itype)) then
        o_itype = itype
    else
        o_itype = 1
    endif
    if(present(jobz)) then
        o_jobz = jobz
    else
        o_jobz = 'N'
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    l_stat_alloc = 0
    allocate(rwork(max(1,3*n-2)))
    lwork = -1
    call zhegv(o_itype,o_jobz,o_uplo,n,a,lda,b,ldb,w,s_work,lwork,rwork,o_info)
    if ( o_info .eq. 0 ) then
        lwork = int(s_work(1))
        ! <<< allocate work arrays with requested sizes >>>
        if(l_stat_alloc==0) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if(l_stat_alloc==0) then
            call zhegv(o_itype,o_jobz,o_uplo,n,a,lda,b,ldb,w,work,lwork,rwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(work)
    endif
    deallocate(rwork)
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_zhegvx(a,b,w,itype,uplo,z,vl,vu,il,iu,m,ifail,abstol,info)
    complex(flyt), intent(inout) :: a(:,:)
    complex(flyt), intent(inout) :: b(:,:)
    real(flyt), intent(out) :: w(:)
    complex(flyt), intent(out), optional, target :: z(:,:)
    integer, intent(out), optional, target :: ifail(:)
    integer, intent(in), optional :: itype
    character(len=1), intent(in), optional :: uplo
    real(flyt), intent(in), optional :: vl
    real(flyt), intent(in), optional :: vu
    integer, intent(in), optional :: il
    integer, intent(in), optional :: iu
    integer, intent(out), optional :: m
    real(flyt), intent(in), optional :: abstol
    integer, intent(out), optional :: info

    character(len=5), parameter :: srname = 'HEGVX'
    character(len=1) :: o_uplo,jobz,range
    real(flyt) :: o_vl, o_vu, o_abstol
    integer :: o_itype,o_il,o_iu,o_m,o_info,n,lda,ldb,ldz,lwork,l_stat_alloc,l_stat_dealloc
    complex(flyt), pointer :: o_z(:,:)
    integer, pointer :: o_ifail(:)
    complex(flyt), pointer :: work(:)
    real(flyt), pointer :: rwork(:)
    integer, pointer :: iwork(:)
    complex(flyt) :: s_work(1)
    integer, target :: l_a1_inte(1)
    complex(flyt), target :: l_a2_comp(1,1)
    !
    if(present(abstol)) then
        o_abstol = abstol
    else
        o_abstol = 0.0_flyt
    endif
    if(present(il)) then
        o_il = il
    else
        o_il = 1
    endif
    if(present(itype)) then
        o_itype = itype
    else
        o_itype = 1
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    if(present(vl)) then
        o_vl = vl
    else
        o_vl = -huge(vl)
    endif
    if(present(vu)) then
        o_vu = vu
    else
        o_vu = huge(vl)
    endif
    if(present(z)) then
        jobz = 'V'
    else
        jobz = 'N'
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    if(present(z)) then
        ldz = max(1,size(z,1))
    else
        ldz = 1
    endif
    n = size(a,2)
    if((present(vl).or.present(vu)).and.(present(il).or.present(iu))) then
        o_info=-1001; goto 1001
    elseif((present(vl).or.present(vu))) then
        range = 'V'
    elseif((present(il).or.present(iu))) then
        range = 'I'
    else
        range = 'A'
    endif
    if(present(iu)) then
        o_iu = iu
    else
        o_iu = n
    endif
    l_stat_alloc = 0
    if(.not.present(z)) then
        if(present(ifail)) then
            o_info=-1001; goto 1001
        else
            o_ifail => l_a1_inte
        endif
    else
        if(present(ifail)) then
            o_ifail => ifail
        else
            allocate(o_ifail(n), stat=l_stat_alloc)
        endif
    endif
    if(present(z)) then
        o_z => z
    else
        o_z => l_a2_comp
    endif
    if(l_stat_alloc==0) then
        allocate(iwork(5*n), stat=l_stat_alloc)
    endif
    if(l_stat_alloc==0) then
        allocate(rwork(7*n), stat=l_stat_alloc)
    endif
    lwork = -1
    call zhegvx(o_itype,jobz,range,o_uplo,n,a,lda,b,ldb,o_vl,o_vu,o_il,o_iu,&
        o_abstol,o_m,w,o_z,ldz,s_work,lwork,rwork,iwork,o_ifail,o_info)
    if(o_info /= 0) then
        goto 200
    endif
    lwork = int(s_work(1))
    if(l_stat_alloc==0) then
        allocate(work(lwork), stat=l_stat_alloc)
    endif
    ! <<< call lapack77 routine >>>
    if(l_stat_alloc==0) then
        call zhegvx(o_itype,jobz,range,o_uplo,n,a,lda,b,ldb,o_vl,o_vu,o_il,o_iu,&
            o_abstol,o_m,w,o_z,ldz,work,lwork,rwork,iwork,o_ifail,o_info)
    else
        o_info = -1000
    endif
    if(present(m)) then
        m = o_m
    endif
    deallocate(work, stat=l_stat_dealloc)
200    continue
    if(present(z) .and..not. present(ifail)) then
        deallocate(o_ifail, stat=l_stat_dealloc)
    endif
    deallocate(iwork, stat=l_stat_dealloc)
    deallocate(rwork, stat=l_stat_dealloc)
1001    continue
    if(present(info)) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla(srname,-o_info)
    endif
end subroutine

subroutine lo_dsygvd(a,b,w,itype,jobz,uplo,info)
    integer, intent(in), optional :: itype
    character(len=1), intent(in), optional :: jobz
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info
    real(flyt), intent(inout) :: a(:,:)
    real(flyt), intent(inout) :: b(:,:)
    real(flyt), intent(out) :: w(:)

    real(flyt), dimension(:), allocatable :: work
    real(flyt) :: s_work(1)
    character(len=1) :: o_jobz,o_uplo
    integer, dimension(:), allocatable :: iwork
    integer :: o_itype,o_info,n,lda,ldb,lwork,liwork,l_stat_alloc, l_stat_dealloc,s_iwork(1)
    if(present(itype)) then
        o_itype = itype
    else
        o_itype = 1
    endif
    if(present(jobz)) then
        o_jobz = jobz
    else
        o_jobz = 'N'
    endif
    if(present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    l_stat_alloc = 0
    liwork = -1
    lwork = -1
    call dsygvd(o_itype,o_jobz,o_uplo,n,a,lda,b,ldb,w,s_work,lwork,s_iwork,liwork,o_info)
    if ( o_info .eq. 0 ) then
        liwork = s_iwork(1)
        lwork = int(anint(s_work(1)))
        allocate(iwork(liwork), stat=l_stat_alloc)
        if ( l_stat_alloc==0 ) then
            allocate(work(lwork), stat=l_stat_alloc)
        endif
        if ( l_stat_alloc==0 ) then
            call dsygvd(o_itype,o_jobz,o_uplo,n,a,lda,b,ldb,w,work,lwork,iwork,liwork,o_info)
        else
            o_info = -1000
        endif
        deallocate(iwork,stat=l_stat_dealloc)
        deallocate(work,stat=l_stat_dealloc)
    endif
    if( present(info) ) then
        info = o_info
    elseif(o_info <= -1000) then
        call xerbla('sygvd',-o_info)
    endif
end subroutine

end module
