#include "precompilerdefinitions"
module type_lennardjones
    !! Lennard-Jones interatomic potential
    use konstanter, only: flyt
    implicit none

    private
    public :: lo_lj 

    !> lennard-jones force-field
    type lo_lj
        !> sigma
        real(flyt) :: sig        
        real(flyt) :: eps
        !> soft cutoff
        real(flyt) :: rc1
        !> hard cutoff
        real(flyt) :: rc2
        !> some shorthand
        real(flyt) :: A
        real(flyt) :: B
        real(flyt) :: twelveA
        real(flyt) :: sixB
        ! spline coefficients
        real(flyt) :: cA,cB,cC,cD
        real(flyt) :: c0,c1,c2,c3,c4,c5
        real(flyt) :: ljshift
        !
        contains
            !> generate spline coefficients and things
            procedure :: setup
    end type

contains

!> precalculate as much as possible
subroutine setup(lj,cutoff,sigma,eps)
    !> the Lennard-Jones potential
    class(lo_lj), intent(out) :: lj
    !> some parameters
    real(flyt), intent(in) :: cutoff,sigma,eps
    !
    real(flyt) :: x1,x2,V1,V2,V3,dx
    !
    lj%eps=eps   ! no idea, picked one at random
    lj%sig=sigma ! no idea, also random
    !
    ! it's neat to convert these to just two number, to the form
    ! A/r^12-B/r^6
    !
    lj%A=4.0_flyt*lj%eps*(lj%sig**12)
    lj%B=4.0_flyt*lj%eps*(lj%sig**6)
    lj%twelveA=12.0_flyt*lj%A ! just to not have to do the multiplication all the time
    lj%sixB=6.0_flyt*lj%B     ! same here
    !
    ! It's good to not have a hard cutoff, but smoothly join it to 0 with a cubic spline. Makes things a lot
    ! more robust with forces and energies that are continous across space.
    !
    lj%rc2=cutoff
    !
    ! the soft cutoff, where I switch to a spline. I make it a little smaller than the hard, by a random number.
    ! this value should be checked that it does not give strange results.
    !
    lj%rc1=cutoff-0.25_flyt
    ! the spline coefficients (thanks mathematica)
    !A -> -((2*V1 - V2*x1 + V2*x2)/(x1 - x2)^3)
    !B -> -((-3*V1*x1 + V2*x1^2 - 3*V1*x2 + V2*x1*x2 -  2*V2*x2^2)/(x1 - x2)^3)
    !C -> (x2*(-6*V1*x1 + 2*V2*x1^2 - V2*x1*x2 - V2*x2^2))/(x1 - x2)^3
    !D -> -((-3*V1*x1*x2^2 + V2*x1^2*x2^2 + V1*x2^3 - V2*x1*x2^3)/(x1 - x2)^3)}}
    x1=lj%rc1
    x2=lj%rc2
    lj%ljshift=lj%A/(x2**12)-lj%B/(x2**6) ! lift the energies
    dx=x1-x2
    V1=lj%A/(x1**12)-lj%B/(x1**6)-lj%ljshift
    V2=-12.0_flyt*lj%A/(x1**13)+6*lj%B/(x1**7)
    V3=156*lj%A/(x1**14)-42*lj%B/(x1**8)
    !
    lj%cA = -((2*V1 - V2*x1 + V2*x2)/(x1 - x2)**3)
    lj%cB = -((-3*V1*x1 + V2*x1**2 - 3*V1*x2 + V2*x1*x2 -  2*V2*x2**2)/(x1 - x2)**3)
    lj%cC = (x2*(-6*V1*x1 + 2*V2*x1**2 - V2*x1*x2 - V2*x2**2))/(x1 - x2)**3
    lj%cD = -((-3*V1*x1*x2**2 + V2*x1**2*x2**2 + V1*x2**3 - V2*x1*x2**3)/(x1 - x2)**3)
    !
    lj%c0=-((x2**3*(2*V1*(10*x1**2-5*x1*x2+x2**2)+dx*x1*(dx*V3*x1+2*V2*(-4*x1+x2))))/(2*dx**5))
    lj%c1=(x2**2*(60*V1*x1**2-dx*(2*V2*(6*x1-x2)*(2*x1+x2)+V3*x1*(-3*x1**2+x1*x2+2*x2**2))))/(2*dx**5)
    lj%c2=(x2*(-60*V1*x1*(x1+x2)+dx*(12*V2*x1*(2*x1+3*x2)-dx*V3*(3*x1**2+6*x1*x2+x2**2))))/(2*dx**5)
    lj%c3=(20*V1*(x1**2+4*x1*x2+x2**2)+dx*(dx*V3*(x1**2+6*x1*x2+3*x2**2)-4*V2*(2*x1**2+10*x1*x2+3*x2**2)))/(2*dx**5)
    lj%c4=(-30*V1*(x1+x2)+dx*(-dx*V3*(2*x1+3*x2)+2*V2*(7*x1+8*x2)))/(2*dx**5)
    lj%c5=(12*V1+dx*(-6*V2+dx*V3))/(2*dx**5)

    ! Square the cutoffs
    lj%rc1=lj%rc1**2
    lj%rc2=lj%rc2**2

    !
end subroutine


end module
