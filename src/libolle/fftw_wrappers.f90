#include "precompilerdefinitions"
module fftw_wrappers
!! I never remember the syntax of fftw so I wrapped it nicely. Also, it does only in-place transforms.
use konstanter, only: r8,lo_imag,lo_pi
use gottochblandat, only : lo_mean
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f'

private
public :: lo_fft
public :: lo_ifft
public :: lo_acf
public :: lo_convolution
public :: lo_abcd_convolution
!public :: lo_kramerkronig_odd

!> Generic interface for in-place forward fft
interface lo_fft
    module procedure lo_complex_inplace_forward_fft_1d
    module procedure lo_complex_inplace_forward_fft_2d
    module procedure lo_complex_inplace_forward_fft_3d
    module procedure lo_real_inplace_forward_fft_1d
    module procedure lo_real_inplace_forward_fft_3d
end interface

!> Generic interface for in-place backwards fft
interface lo_ifft
    module procedure lo_complex_inplace_backward_fft_1d
    module procedure lo_complex_inplace_backward_fft_2d
    module procedure lo_complex_inplace_backward_fft_3d
    module procedure lo_real_inplace_backward_fft_1d
end interface

!> Generic interface for in-place autocorrelation functions
interface lo_acf
    module procedure lo_real_inplace_autocorrelation_1d
end interface

!> Interface for convolutions
interface lo_convolution
    module procedure lo_real_convolution
end interface

contains

!> Kramers-Kronig transform of an odd function on a uniform x-axis.
! Actually, it was fast but not accurate enough, aliasing thing that breaks everything.
! If I have the time I can work with the convolutions a bit to see if I can make it work.
! subroutine lo_kramerkronig_odd(x,im,re)
!     !> x-axis
!     real(r8), dimension(:), intent(in) :: x
!     !> imaginary part
!     real(r8), dimension(:), intent(in) :: im
!     !> real part
!     real(r8), dimension(:), intent(out) :: re
!
!     type(c_ptr) :: plan
!     complex(c_double_complex), dimension(:), allocatable :: dy1,dy2,dy3,df1,df2,df3
!     integer :: n,i,j
!
!     n=size(x)
!     allocate(dy1(2*n-1))
!     allocate(df1(2*n-1))
!     dy1=0.0_r8
!     df1=0.0_r8
!     ! Make it clear the function is odd
!     do i=1,n
!         j=n-i+1
!         dy1(j)=-im(i)
!         j=n+i-1
!         dy1(j)=im(i)
!     enddo
!     ! Fourier transform
!     call dfftw_plan_dft_1d(plan,2*n-1,dy1,df1,FFTW_FORWARD,FFTW_ESTIMATE)
!     call dfftw_execute_dft(plan,dy1,df1)
!     call dfftw_destroy_plan(plan)
!     ! Multiply by signum
!     do i=1,n !(-1)? not sure
!         df1(i)=-df1(i)
!     enddo
!     df1(n)=0.0_r8
!     df1=df1*lo_imag
!
!     ! Inverse transform?
!     call dfftw_plan_dft_1d(plan,2*n-1,df1,dy1,FFTW_BACKWARD,FFTW_ESTIMATE)
!     call dfftw_execute_dft(plan,df1,dy1)
!     call dfftw_destroy_plan(plan)
!     ! return the proper interval
!     re=real(dy1(n:2*n-1),r8)/real((2*n-1),r8)
! end subroutine

!> Convolution of distributions
subroutine lo_real_convolution(y1,y2,dx,normalize,y3)
    !> functions to be convoluted. Must be on uniform x-axis, and the same x-axis
    real(r8), dimension(:), intent(in) :: y1,y2
    !> step on the x-axis
    real(r8), intent(in) :: dx
    !> should the result be normalized?
    logical, intent(in) :: normalize
    !> resulting convolution
    real(r8), dimension(:), intent(out) :: y3

    integer, parameter :: brute_force_crossover=700 ! Quick and dirty tests seems to have it this way.
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:), allocatable :: dy1,dy2,df1,df2
    real(r8) :: f0
    integer :: i,j,n,m

    n=size(y1)
    if ( n > brute_force_crossover ) then
        ! Large enough so that it is worth bothering with FFTs
        allocate(dy1(n))
        allocate(dy2(n))

        allocate(df1(n))
        allocate(df2(n))

        ! Start with a forward transform, and subtract the mean for some reason
        dy1=y1
        dy2=y2
        df1=0.0_r8
        df2=0.0_r8
        call dfftw_plan_dft_1d(plan,n,dy1,df1,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,dy1,df1)
        call dfftw_execute_dft(plan,dy2,df2)
        call dfftw_destroy_plan(plan)
        ! Multiply?
        dy1=df1*df2
        df1=0.0_r8
        ! Inverse transform it
        call dfftw_plan_dft_1d(plan,n,dy1,df1,FFTW_BACKWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,dy1,df1)
        call dfftw_destroy_plan(plan)
        ! Do the FFT shift thingy
        if ( mod(n,2) .eq. 0 ) then
            ! Even number of points
            m=n/2
            y3(1:m)=real(df1( m+1:n ),r8)
            y3(m+1:n)=real(df1( 1:m ),r8)
        else
            ! Odd number of points
            m=(n-1)/2
            y3(1:m) = real(df1( m+2:n ),r8)
            y3(m+1:n) = real(df1( 1:m+1 ),r8)
        endif
        ! Scale it to the first thing.
        y3=y3*dx/real(n,r8)
        ! cleanup
        deallocate(dy1)
        deallocate(dy2)
        deallocate(df1)
        deallocate(df2)
    else
        ! Small convolution, just brute-force it
        y3=0.0_r8
        do j=1,n
            do i=j+1,n
                y3(i)=y3(i)+y1(j)*y2(i-j)
            enddo
        enddo
        y3=y3*dx
    endif

    ! Perhaps normalize it?
    if ( normalize ) then
        f0=( sum(y3)-(y3(1)+y3(n))*0.5_r8 )*dx
        y3=y3/f0
    endif
end subroutine

!> Calculates y1 o y2 + z1 o z2.
subroutine lo_abcd_convolution(y1,y2,z1,z2,dx,s)
    !> functions to be convoluted
    real(r8), dimension(:), intent(in) :: y1,y2,z1,z2
    !> step on the x-axis
    real(r8), intent(in) :: dx
    !> resulting convolution
    real(r8), dimension(:), intent(out) :: s

    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:), allocatable :: dy1,dy2,dz1,dz2,df1,df2,df3,df4
    integer :: n

    n=size(y1)
    allocate(dy1(n))
    allocate(dy2(n))
    allocate(dz1(n))
    allocate(dz2(n))
    allocate(df1(n))
    allocate(df2(n))
    allocate(df3(n))
    allocate(df4(n))
    dy1=y1
    dy2=y2
    dz1=z1
    dz2=z2
    df1=0.0_r8
    df2=0.0_r8
    df3=0.0_r8
    df4=0.0_r8
    ! Forward transform
    call dfftw_plan_dft_1d(plan,n,dy1,df1,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dy1,df1)
    call dfftw_execute_dft(plan,dy2,df2)
    call dfftw_execute_dft(plan,dz1,df3)
    call dfftw_execute_dft(plan,dz2,df4)
    call dfftw_destroy_plan(plan)
    ! Convolution
    dy1=df1*df2 + df3*df4
    df1=0.0_r8
    ! Inverse transforms
    call dfftw_plan_dft_1d(plan,n,dy1,df1,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dy1,df1)
    call dfftw_destroy_plan(plan)
    ! Add scaling factors
    s=real(df1,r8)*dx/real(n,r8)
    ! cleanup
    deallocate(dy1)
    deallocate(dy2)
    deallocate(dz1)
    deallocate(dz2)
    deallocate(df1)
    deallocate(df2)
    deallocate(df3)
    deallocate(df4)
end subroutine

!> In-place autocorrelation function
subroutine lo_real_inplace_autocorrelation_1d(y)
    real(r8), dimension(:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:), allocatable :: dy,df
    integer :: i,n

    n=size(y,1)
    allocate(dy(n))
    allocate(df(n))
    ! Start with a forward transform, and subtract the mean for some reason
    dy=y-lo_mean(y)
    df=0.0_r8
    call dfftw_plan_dft_1d(plan,n,dy,df,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dy,df)
    call dfftw_destroy_plan(plan)
    ! Get the norm squared
    do i=1,n
        dy(i)=abs(conjg(df(i))*df(i))
    enddo
    df=0.0_r8
    ! Inverse transform it
    call dfftw_plan_dft_1d(plan,n,dy,df,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dy,df)
    call dfftw_destroy_plan(plan)
    ! Scale it to the first thing.
    y=real(df)/n
    ! cleanup
    deallocate(dy)
    deallocate(df)
end subroutine

subroutine lo_real_inplace_forward_fft_1d(y,win,wout,wplan)
    real(r8), dimension(:), intent(inout) :: y
    real(c_double), dimension(:), intent(inout), optional :: win,wout
    type(c_ptr), intent(in), optional :: wplan
    !
    type(c_ptr) :: plan
    real(c_double), dimension(:), allocatable :: dumin,dumut
    integer :: n1
    n1=size(y,1)
    if ( present(win) .and. present(wout) .and. present(wplan) ) then
        ! if I have provided workspace and a plan, use that.
        win=y
        call dfftw_execute_r2r(wplan,win,wout)
        y=wout/sqrt(n1*1.0_r8)
    else
        ! if not, create temporary space and a plan
        allocate(dumin(n1),dumut(n1))
        dumin=y
        call dfftw_plan_r2r_1d(plan,n1,dumin,dumut,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_r2r(plan,dumin,dumut)
        call dfftw_destroy_plan(plan)
        y=dumut/sqrt(n1*1.0_r8)
        deallocate(dumin,dumut)
    endif
end subroutine
subroutine lo_real_inplace_forward_fft_3d(y,win,wout,wplan)
    real(r8), dimension(:,:,:), intent(inout) :: y
    real(c_double), dimension(:,:,:), intent(inout), optional :: win,wout
    type(c_ptr), intent(in), optional :: wplan
    !
    type(c_ptr) :: plan
    real(c_double), dimension(:,:,:), allocatable :: dumin,dumut
    integer :: n1,n2,n3
    n1=size(y,1)
    n2=size(y,2)
    n3=size(y,3)
    if ( present(win) .and. present(wout) .and. present(wplan) ) then
        ! if I have provided workspace and a plan, use that.
        win=y
        call dfftw_execute_r2r(wplan,win,wout)
        y=wout/sqrt(n1*n2*n3*1.0_r8)
    else
        ! if not, create temporary space and a plan
        allocate(dumin(n1,n2,n3),dumut(n1,n2,n3))
        dumin=y
        dumut=0.0_r8
        call dfftw_plan_r2r_3d(plan,n1,n2,n3,dumin,dumut,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_r2r(plan,dumin,dumut)
        call dfftw_destroy_plan(plan)
        y=dumut/sqrt(n1*n2*n3*1.0_r8)
        deallocate(dumin,dumut)
    endif
end subroutine



subroutine lo_real_inplace_backward_fft_1d(y)
    real(r8), dimension(:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    real(c_double), dimension(:), allocatable :: dumin,dumut
    integer :: n1
    n1=size(y,1)
    allocate(dumin(n1),dumut(n1))
    dumin=y
    call dfftw_plan_r2r_1d(plan,n1,dumin,dumut,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_r2r(plan,dumin,dumut)
    call dfftw_destroy_plan(plan)
    y=dumut/sqrt(n1*1.0_r8)
    deallocate(dumin,dumut)
end subroutine

subroutine lo_complex_inplace_forward_fft_1d(y)
    complex(r8), dimension(:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:), allocatable :: dumin,dumut
    integer :: n1
    n1=size(y,1)
    allocate(dumin(n1),dumut(n1))
    dumin=y
    call dfftw_plan_dft_1d(plan,n1,dumin,dumut,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dumin,dumut)
    call dfftw_destroy_plan(plan)
    y=dumut/sqrt(n1*1.0_r8)
    deallocate(dumin,dumut)
end subroutine
subroutine lo_complex_inplace_backward_fft_1d(y)
    complex(r8), dimension(:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:), allocatable :: dumin,dumut
    integer :: n1
    n1=size(y,1)
    allocate(dumin(n1),dumut(n1))
    dumin=y
    call dfftw_plan_dft_1d(plan,n1,dumin,dumut,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dumin,dumut)
    call dfftw_destroy_plan(plan)
    y=dumut/sqrt(n1*1.0_r8)
    deallocate(dumin,dumut)
end subroutine

subroutine lo_complex_inplace_forward_fft_2d(y)
    !> data to be transformed
    complex(r8), dimension(:,:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:,:), allocatable :: dumin,dumut
    integer :: n1,n2
    n1=size(y,1)
    n2=size(y,2)
    allocate(dumin(n1,n2),dumut(n1,n2))
    dumin=y
    call dfftw_plan_dft_2d(plan,n1,n2,dumin,dumut,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dumin,dumut)
    call dfftw_destroy_plan(plan)
    y=dumut/sqrt(n1*n2*1.0_r8)
    deallocate(dumin,dumut)
end subroutine
subroutine lo_complex_inplace_backward_fft_2d(y)
    !> data to be transformed
    complex(r8), dimension(:,:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:,:), allocatable :: dumin,dumut
    !
    integer :: n1,n2
    !
    n1=size(y,1)
    n2=size(y,2)
    allocate(dumin(n1,n2),dumut(n1,n2))
    !
    dumin=y
    call dfftw_plan_dft_2d(plan,n1,n2,dumin,dumut,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dumin,dumut)
    call dfftw_destroy_plan(plan)
    y=dumut/sqrt(n1*n2*1.0_r8)
    deallocate(dumin,dumut)
end subroutine

subroutine lo_complex_inplace_forward_fft_3d(y)
    !> data to be transformed
    complex(r8), dimension(:,:,:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:,:,:), allocatable :: dumin,dumut
    !
    integer :: n1,n2,n3
    !
    n1=size(y,1)
    n2=size(y,2)
    n3=size(y,3)
    allocate(dumin(n1,n2,n3),dumut(n1,n2,n3))
    !
    dumin=y
    call dfftw_plan_dft_3d(plan,n1,n2,n3,dumin,dumut,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dumin,dumut)
    call dfftw_destroy_plan(plan)
    y=dumut/sqrt(n1*n2*n3*1.0_r8)
    deallocate(dumin,dumut)
end subroutine
subroutine lo_complex_inplace_backward_fft_3d(y)
    !> data to be transformed
    complex(r8), dimension(:,:,:), intent(inout) :: y
    !
    type(c_ptr) :: plan
    complex(c_double_complex), dimension(:,:,:), allocatable :: dumin,dumut
    !
    integer :: n1,n2,n3
    !
    n1=size(y,1)
    n2=size(y,2)
    n3=size(y,3)
    allocate(dumin(n1,n2,n3),dumut(n1,n2,n3))
    !
    dumin=y
    call dfftw_plan_dft_3d(plan,n1,n2,n3,dumin,dumut,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dumin,dumut)
    call dfftw_destroy_plan(plan)
    y=dumut/sqrt(n1*n2*n3*1.0_r8)
    deallocate(dumin,dumut)
end subroutine

end module
