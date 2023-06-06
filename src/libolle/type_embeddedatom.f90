#include "precompilerdefinitions"
module type_embeddedatom
    !! Finnis-Sinclair type classical potentials
    use konstanter, only: flyt
    implicit none

    private
    public :: lo_eam

    type lo_eam
        real(flyt) :: a, eps, C, m, n
        real(flyt) :: r1,r2
        real(flyt) :: A1,A2,A3,A4,A5,A6,energyshift,splineenergyshift
        real(flyt) :: B1,B2,B3,B4,B5,B6,embedshift,splineembedshift
        !
        contains
            procedure :: setup
            procedure :: pairf
            procedure :: pairfd
            procedure :: embedf
            procedure :: embedfd
    end type

contains

!> sets the hard and soft cutoff, as well as all the spline coefficients and shifts
subroutine setup(eam,a,eps,C,m,n,cutoff)
    !> eam potential
    class(lo_eam), intent(out) :: eam
    !> eq lattice parameter
    real(flyt), intent(in) :: a
    !> energy scale
    real(flyt), intent(in) :: eps
    !> density scale
    real(flyt), intent(in) :: C
    !> exponent for repulsion
    real(flyt), intent(in) :: n
    !> exponent for embedding function
    real(flyt), intent(in) :: m
    !> hard cutoff
    real(flyt), intent(in) :: cutoff
    !
    real(flyt) :: r1,r2,f0,f1,f2
    real(flyt) :: A1,A2,A3,A4,A5,A6
    real(flyt) :: B1,B2,B3,B4,B5,B6
    !
    r2=cutoff              ! hard cutoff
    r1=cutoff-a*0.1_flyt  ! soft cutoff
    !
    eam%a=a
    eam%eps=eps
    eam%C=C
    eam%m=m
    eam%n=n
    eam%r1=r1
    eam%r2=r2

    ! Get spline coefficients for the pair potential
    A1 = (n*(a/r1)**n*r2**3*((30 + 11*n + n**2)*r1**2 - 2*(12 + 8*n + n**2)*r1*r2 + (6 + 5*n + n**2)*r2**2))/(2*r1*(r1 - r2)**5)
    A2 = -(n*(a/r1)**n*r2**2*(3*(30 + 11*n + n**2)*r1**3 - 4*(6 + 7*n + n**2)*r1**2*r2 - (12 + 13*n + n**2)*r1*r2**2 + 2*(3 + 4*n + n**2)*r2**3))/(4*r1**2*(r1 - r2)**5)
    A3 = (n*(a/r1)**n*r2*(3*(30 + 11*n + n**2)*r1**4 + 12*(6 + n)*r1**3*r2 - 4*(13 + 15*n + 2*n**2)*r1**2*r2**2 + 4*(2 + 3*n + n**2)*r1*r2**3 + (2 + 3*n + n**2)*r2**4))/(6*r1**3*(r1 - r2)**5)
    A4 = -(n*(a/r1)**n*((30 + 11*n + n**2)*r1**4 + 4*(30 + 11*n + n**2)*r1**3*r2 - 4*(6 + 13*n + 2*n**2)*r1**2*r2**2 - 12*(1 + n)*r1*r2**3 + 3*(2 + 3*n + n**2)*r2**4))/(8*r1**3*(r1 - r2)**5)
    A5 = (n*(a/r1)**n*(2*(24 + 10*n + n**2)*r1**3 - (-30 + n + n**2)*r1**2*r2 - 4*(6 + 7*n + n**2)*r1*r2**2 + 3*(2 + 3*n + n**2)*r2**3))/(10*r1**3*(r1 - r2)**5)
    A6 = -(n*(a/r1)**n*((20 + 9*n + n**2)*r1**2 - 2*(5 + 6*n + n**2)*r1*r2 + (2 + 3*n + n**2)*r2**2))/(12*r1**3*(r1 - r2)**5)
    ! Get the energy shift thing
    f0=A6*r1**6+A5*r1**5+A4*r1**4+A3*r1**3+A2*r1**2+A1*r1
    f1=A6*r2**6+A5*r2**5+A4*r2**4+A3*r2**3+A2*r2**2+A1*r2
    f2=(a/r1)**n
    eam%energyshift=-f2+f0-f1
    eam%splineenergyshift=-f1
    eam%A1=A1
    eam%A2=A2
    eam%A3=A3
    eam%A4=A4
    eam%A5=A5
    eam%A6=A6

    ! Get spline coefficients for the embedding function
    B1 = (m*(a/r1)**m*r2**3*((30 + 11*m + m**2)*r1**2 - 2*(12 + 8*m + m**2)*r1*r2 + (6 + 5*m + m**2)*r2**2))/(2*r1*(r1 - r2)**5)
    B2 = -(m*(a/r1)**m*r2**2*(3*(30 + 11*m + m**2)*r1**3 - 4*(6 + 7*m + m**2)*r1**2*r2 - (12 + 13*m + m**2)*r1*r2**2 + 2*(3 + 4*m + m**2)*r2**3))/(4*r1**2*(r1 - r2)**5)
    B3 = (m*(a/r1)**m*r2*(3*(30 + 11*m + m**2)*r1**4 + 12*(6 + m)*r1**3*r2 - 4*(13 + 15*m + 2*m**2)*r1**2*r2**2 + 4*(2 + 3*m + m**2)*r1*r2**3 + (2 + 3*m + m**2)*r2**4))/(6*r1**3*(r1 - r2)**5)
    B4 = -(m*(a/r1)**m*((30 + 11*m + m**2)*r1**4 + 4*(30 + 11*m + m**2)*r1**3*r2 - 4*(6 + 13*m + 2*m**2)*r1**2*r2**2 - 12*(1 + m)*r1*r2**3 + 3*(2 + 3*m + m**2)*r2**4))/(8*r1**3*(r1 - r2)**5)
    B5 = (m*(a/r1)**m*(2*(24 + 10*m + m**2)*r1**3 - (-30 + m + m**2)*r1**2*r2 - 4*(6 + 7*m + m**2)*r1*r2**2 + 3*(2 + 3*m + m**2)*r2**3))/(10*r1**3*(r1 - r2)**5)
    B6 = -(m*(a/r1)**m*((20 + 9*m + m**2)*r1**2 - 2*(5 + 6*m + m**2)*r1*r2 + (2 + 3*m + m**2)*r2**2))/(12*r1**3*(r1 - r2)**5)
    ! Shift thing
    f0=B6*r1**6+B5*r1**5+B4*r1**4+B3*r1**3+B2*r1**2+B1*r1
    f1=B6*r2**6+B5*r2**5+B4*r2**4+B3*r2**3+B2*r2**2+B1*r2
    f2=(a/r1)**m
    eam%embedshift=-f2+f0-f1
    eam%splineembedshift=-f1
    eam%B1=B1
    eam%B2=B2
    eam%B3=B3
    eam%B4=B4
    eam%B5=B5
    eam%B6=B6
    !
end subroutine

!> the eam pair potential
pure function pairf(eam,r) result(U)
    class(lo_eam), intent(in) :: eam
    real(flyt), intent(in) :: r
    real(flyt) :: U
    !
    if ( r .ge. eam%r2 ) then
        U=0.0_flyt
        return
    elseif( r .gt. eam%r1 ) then
        U=eam%A6*r**6+eam%A5*r**5+eam%A4*r**4+eam%A3*r**3+eam%A2*r**2+eam%A1*r+eam%splineenergyshift
        U=U*eam%eps
    else
        U=(eam%a/r)**eam%n+eam%energyshift
        U=U*eam%eps
    endif
end function

!> the eam embedding function
pure function embedf(eam,r) result(U)
    class(lo_eam), intent(in) :: eam
    real(flyt), intent(in) :: r
    real(flyt) :: U
    !
    if ( r .ge. eam%r2 ) then
        U=0.0_flyt
        return
    elseif( r .gt. eam%r1 ) then
        U=eam%B6*r**6+eam%B5*r**5+eam%B4*r**4+eam%B3*r**3+eam%B2*r**2+eam%B1*r+eam%splineembedshift
    else
        U=(eam%a/r)**eam%m+eam%embedshift
    endif
end function

!> the eam pair potential derivative
pure function pairfd(eam,r) result(U)
    class(lo_eam), intent(in) :: eam
    real(flyt), intent(in) :: r
    real(flyt) :: U
    !
    if ( r .ge. eam%r2 ) then
        U=0.0_flyt
        return
    elseif( r .gt. eam%r1 ) then
        U=6*eam%A6*(r**5)+5*eam%A5*(r**4)+4*eam%A4*(r**3)+3*eam%A3*(r**2)+2*eam%A2*r+eam%A1
        U=U*eam%eps
    else
        U=-eam%n*((eam%a/r)**eam%n)/r
        U=U*eam%eps
    endif
end function

!> the eam embedding derivative
pure function embedfd(eam,r) result(U)
    class(lo_eam), intent(in) :: eam
    real(flyt), intent(in) :: r
    real(flyt) :: U
    !
    if ( r .ge. eam%r2 ) then
        U=0.0_flyt
        return
    elseif( r .gt. eam%r1 ) then
        U=6*eam%B6*r**5+5*eam%B5*r**4+4*eam%B4*r**3+3*eam%B3*r**2+2*eam%B2*r+eam%B1
        !U=U*eam%C
    else
        U=-eam%m*((eam%a/r)**eam%m)/r
        !U=U*eam%C
    endif
end function

end module
