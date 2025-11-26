module lo_spectralfunction_convolution
use konstanter, only: r8,lo_huge,lo_hugeint,lo_exitcode_param
use gottochblandat, only: lo_planck,lo_stop_gracefully
use, intrinsic :: iso_c_binding, only: c_ptr
implicit none
include 'fftw3.f'

type lo_convolution_handle
    !> number of points on frequency grid (in one direction, 0:max)
    integer :: n=-lo_hugeint
    !> step in frequency axis
    real(r8) :: deltaomega=-lo_huge
    !> frequency axis
    real(r8), dimension(:), allocatable :: x
    !> bose einstein factor
    real(r8), dimension(:), allocatable :: thermal_n
    real(r8), dimension(:), allocatable :: thermal_n_plus_one
    !> buffers for greater and lesser Green's functions
    complex(r8), dimension(:,:), allocatable :: greater1
    complex(r8), dimension(:,:), allocatable :: greater2
    complex(r8), dimension(:,:), allocatable :: greater3
    complex(r8), dimension(:,:), allocatable :: lesser1
    complex(r8), dimension(:,:), allocatable :: lesser2
    complex(r8), dimension(:,:), allocatable :: lesser3
    !> temporary buffers
    complex(r8), dimension(:), allocatable :: rbuf
    complex(r8), dimension(:), allocatable :: cbuf
    !> FFTW plans
    type(c_ptr) :: plan_fwd, plan_bwd
    contains
        procedure :: generate=>create_convolution_handle
        procedure :: buffer_and_transform_J_to_greater_lesser
        procedure :: buffer_and_transform_J
        procedure :: destroy=>destroy_convolution_helper
        procedure :: inverse_transform_to_real
end type

private
public :: lo_convolution_handle

contains

subroutine create_convolution_handle(ch,omega,temperature,n_mode)
    !> convolution helper
    class(lo_convolution_handle), intent(inout) :: ch
    !> frequency grid
    real(r8), dimension(:), intent(in) :: omega
    !> temperature
    real(r8), intent(in) :: temperature
    !> number of modes
    integer, intent(in) :: n_mode

    integer :: i

    ! Number of points we are interested in?
    ch%n = size(omega) - 1
    ! Some space for buffers
    allocate(ch%x(-ch%n:ch%n))
    allocate(ch%thermal_n(-ch%n:ch%n))
    allocate(ch%thermal_n_plus_one(-ch%n:ch%n))
    allocate(ch%greater1(-ch%n:ch%n, n_mode))
    allocate(ch%greater2(-ch%n:ch%n, n_mode))
    allocate(ch%greater3(-ch%n:ch%n, n_mode))
    allocate(ch%lesser1(-ch%n:ch%n, n_mode))
    allocate(ch%lesser2(-ch%n:ch%n, n_mode))
    allocate(ch%lesser3(-ch%n:ch%n, n_mode))
    allocate(ch%cbuf(-ch%n:ch%n))
    allocate(ch%rbuf(-ch%n:ch%n))

    ch%x=0.0_r8
    ch%thermal_n=0.0_r8
    ch%thermal_n_plus_one=0.0_r8
    ch%greater1=0.0_r8
    ch%greater2=0.0_r8
    ch%greater3=0.0_r8
    ch%lesser1=0.0_r8
    ch%lesser2=0.0_r8
    ch%lesser3=0.0_r8
    ch%cbuf=0.0_r8
    ch%rbuf=0.0_r8

    ! Get an x-axis, might come handy. x(0)=0
    ch%x(0) = 0.0_r8
    do i = 1, ch%n
        ch%x(i) = omega(i + 1)
        ch%x(-i) = -ch%x(i)
    end do
    ch%deltaomega = ch%x(1) - ch%x(0)

    ! Pre-calculate the planck factors?
    do i = 1, ch%n
        ch%thermal_n(i) = lo_planck(temperature, ch%x(i))
        ch%thermal_n(-i) = -lo_planck(temperature, ch%x(i))-1.0_r8
    end do
    ch%thermal_n_plus_one=ch%thermal_n+1.0_r8

    ! Pre-plan the FFTs
    call dfftw_plan_dft_1d(ch%plan_fwd, size(ch%cbuf), ch%cbuf, ch%cbuf, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(ch%plan_bwd, size(ch%cbuf), ch%cbuf, ch%cbuf, FFTW_BACKWARD, FFTW_ESTIMATE)
end subroutine

!> inverse transform with the shifts and appropriate prefactors
subroutine inverse_transform_to_real(ch,tbuf,fbuf)
    !> convolution helper
    class(lo_convolution_handle), intent(inout) :: ch
    !> complex buffer in time domain (-n:n)
    complex(r8), dimension(:), intent(in) :: tbuf
    !> real buffer in frequency domain (1:n+1)
    real(r8), dimension(:), intent(out) :: fbuf

    ! First we inverse transform
    call dfftw_execute_dft(ch%plan_bwd, tbuf, ch%cbuf)

    ! Do the FFT-shift thingy
    fbuf(1:ch%n + 1) = real(ch%cbuf(-ch%n:0), r8)
    ! Scale to sensible thing
    fbuf = fbuf*ch%deltaomega/real(size(ch%x), r8)
end subroutine

!> Just fourier transform the spectral function and store in the "lesser" buffer. Not great naming.
subroutine buffer_and_transform_J(ch,spectralfunction,whereto)
    !> convolution helper
    class(lo_convolution_handle), intent(inout) :: ch
    !> spectral function
    real(r8), dimension(:,:), intent(in) :: spectralfunction
    !> where should it be stored
    integer, intent(in) :: whereto

    integer :: imode,n_mode,ie

    n_mode=size(spectralfunction,2)
    do imode=1,n_mode
        ! buffer entire spectral function, forwards and backwards
        ch%rbuf(0:ch%n)=spectralfunction(:,imode)
        do ie=1,ch%n
            ch%rbuf(-ie)=-spectralfunction(ie+1,imode)
        enddo
        ! No prefactor here
        select case(whereto)
        case(1)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%lesser1(:, imode))
        case(2)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%lesser2(:, imode))
        case(3)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%lesser3(:, imode))
        case default
            call lo_stop_gracefully(['Unknown spectral function destination'],lo_exitcode_param,__FILE__,__LINE__)
        end select
    enddo
end subroutine

subroutine buffer_and_transform_J_to_greater_lesser(ch,spectralfunction,whereto)
    !> convolution helper
    class(lo_convolution_handle), intent(inout) :: ch
    !> spectral function
    real(r8), dimension(:,:), intent(in) :: spectralfunction
    !> where should it be stored
    integer, intent(in) :: whereto

    integer :: imode,n_mode,ie

    n_mode=size(spectralfunction,2)
    do imode=1,n_mode
        ! buffer entire spectral function, forwards and backwards
        ch%rbuf(0:ch%n)=spectralfunction(:,imode)
        do ie=1,ch%n
            ch%rbuf(-ie)=-spectralfunction(ie+1,imode)
        enddo
        ! Create lesser greens function
        ch%cbuf=ch%rbuf*ch%thermal_n
        select case(whereto)
        case(1)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%lesser1(:, imode))
        case(2)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%lesser2(:, imode))
        case(3)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%lesser3(:, imode))
        case default
            call lo_stop_gracefully(['Unknown spectral function destination'],lo_exitcode_param,__FILE__,__LINE__)
        end select
        ! Create greater greens function
        ch%cbuf=ch%rbuf*ch%thermal_n_plus_one
        select case(whereto)
        case(1)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%greater1(:, imode))
        case(2)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%greater2(:, imode))
        case(3)
            call dfftw_execute_dft(ch%plan_fwd, ch%cbuf, ch%greater3(:, imode))
        case default
            call lo_stop_gracefully(['Unknown spectral function destination'],lo_exitcode_param,__FILE__,__LINE__)
        end select
    enddo
end subroutine

subroutine destroy_convolution_helper(ch)
    !> convolution helper
    class(lo_convolution_handle), intent(inout) :: ch

    call dfftw_destroy_plan(ch%plan_fwd)
    call dfftw_destroy_plan(ch%plan_bwd)

    if ( allocated(ch%x                 ) ) deallocate(ch%x                 )
    if ( allocated(ch%thermal_n         ) ) deallocate(ch%thermal_n         )
    if ( allocated(ch%thermal_n_plus_one) ) deallocate(ch%thermal_n_plus_one)
    if ( allocated(ch%greater1          ) ) deallocate(ch%greater1          )
    if ( allocated(ch%greater2          ) ) deallocate(ch%greater2          )
    if ( allocated(ch%greater3          ) ) deallocate(ch%greater3          )
    if ( allocated(ch%lesser1           ) ) deallocate(ch%lesser1           )
    if ( allocated(ch%lesser2           ) ) deallocate(ch%lesser2           )
    if ( allocated(ch%lesser3           ) ) deallocate(ch%lesser3           )
    if ( allocated(ch%rbuf              ) ) deallocate(ch%rbuf              )
    if ( allocated(ch%cbuf              ) ) deallocate(ch%cbuf              )
end subroutine

end module