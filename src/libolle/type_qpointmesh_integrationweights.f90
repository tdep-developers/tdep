submodule (type_qpointmesh) type_qpointmesh_integrationweights
use konstanter, only: lo_kb_hartree,lo_tiny
use lo_voronoi, only: lo_voronoi_cell
use lo_verletboxes, only: lo_verletbox
use lo_sorting, only: lo_qsort
use type_distancetable, only: lo_distancetable
use gottochblandat, only: lo_return_unique
implicit none
contains

!> The Lambin & Vigneron tetrahedron method for the resolvent operator
module pure function lo_LV_complex_tetrahedron_weights(ein,z,tol) result(weights)
    !> dispersion values at the corners of the tetrahedron
    real(r8), dimension(4), intent(in) :: ein
    !> energy to evaulate the weights at
    complex(r8), intent(in) :: z
    !> tolerance for when something is zero
    real(r8), intent(in) :: tol
    !> integration weights
    complex(r8), dimension(4) :: weights

    real(r8), parameter :: half=0.5_r8,third=1.0_r8/3.0_r8
    complex(r8), dimension(4) :: cw
    complex(r8) :: EZ1,EZ2,EZ3,EZ4
    real(r8), dimension(4) :: ei
    real(r8) :: Emin,Emax,Zdist
    real(r8) :: E12,E13,E14,E23,E24,E34
    real(r8) :: a,b,c,d,e,f
    integer, dimension(4) :: ind
    integer :: whatkind

    ! Start by sorting the tetrahedron
    ind=sort_four_numbers(ein)
    ei=ein(ind)
    ! Min and max energy
    Emin=ei(1)
    Emax=ei(4)

    ! First the complex energy differences
    EZ1=z-ei(1)
    EZ2=z-ei(2)
    EZ3=z-ei(3)
    EZ4=z-ei(4)
    ! Smallest distance |z-Ei|, to determine wether I should switch to the
    ! asymptotic behavior for numerical stability.
    Zdist=lo_huge
    Zdist=min(Zdist,abs(EZ1))
    Zdist=min(Zdist,abs(EZ2))
    Zdist=min(Zdist,abs(EZ3))
    Zdist=min(Zdist,abs(EZ4))
    !@TODO add asymptotic thing with continued fractions

    ! Then the energy differences, for the coefficients. Must always be positive, I hope.
    E12=ei(2)-ei(1)
    E13=ei(3)-ei(1)
    E14=ei(4)-ei(1)
    E23=ei(3)-ei(2)
    E24=ei(4)-ei(2)
    E34=ei(4)-ei(3)
    a=0.0_r8
    b=0.0_r8
    c=0.0_r8
    d=0.0_r8
    e=0.0_r8
    f=0.0_r8
    if ( E12 .gt. tol ) a=1.0_r8/E12
    if ( E13 .gt. tol ) b=1.0_r8/E13
    if ( E14 .gt. tol ) c=1.0_r8/E14
    if ( E23 .gt. tol ) d=1.0_r8/E23
    if ( E24 .gt. tol ) e=1.0_r8/E24
    if ( E34 .gt. tol ) f=1.0_r8/E34
    ! Try to figure out all the special cases:
    ! case(1) e1<e2<e3<e4
    ! case(2) e1=e2
    ! case(3)    e2=e3
    ! case(4)       e3=e4
    ! case(5) e1=e2 e3=e4
    ! case(6) e1=e2=e3
    ! case(7)    e2=e3=e4
    ! case(8) e1=e2=e3=e4
    whatkind=0
    if ( E12+E23+E34 .lt. tol ) then
        whatkind=8
    elseif ( E23+E34 .lt. tol ) then
        whatkind=7
    elseif ( E12+E23 .lt. tol ) then
        whatkind=6
    elseif ( E12+E34 .lt. tol ) then
        whatkind=5
    elseif ( E34 .lt. tol ) then
        whatkind=4
    elseif ( E23 .lt. tol ) then
        whatkind=3
    elseif ( E12 .lt. tol ) then
        whatkind=2
    else
        whatkind=1
    endif
    ! Now get the actual weights
    cw=0.0_r8
    select case(whatkind)
    case(1) ! e1<e2<e3<e4
        cw(1)=a**2*d*e*EZ2**3*Log(EZ2/EZ1) - b**2*d*EZ3**3*f*Log(EZ3/EZ1) + c*(a*b*EZ1**2 + c*e*EZ4**3*f*Log(EZ4/EZ1))
        cw(2)=a**2*b*c*EZ1**3*Log(EZ1/EZ2) - b*d**2*EZ3**3*f*Log(EZ3/EZ2) + e*(-(a*d*EZ2**2) + c*e*EZ4**3*f*Log(EZ4/EZ2))
        cw(3)=a*b**2*c*EZ1**3*Log(EZ1/EZ3) - a*d**2*e*EZ2**3*Log(EZ2/EZ3) + f*(b*d*EZ3**2 + c*e*EZ4**3*f*Log(EZ4/EZ3))
        cw(4)=a*b*c**2*EZ1**3*Log(EZ1/EZ4) - a*d*e**2*EZ2**3*Log(EZ2/EZ4) + f*(-(c*e*EZ4**2) + b*d*EZ3**3*f*Log(EZ3/EZ4))
    case(2) ! e1=e2
        cw(1)=b*c*EZ1*(half - b*EZ3 - c*EZ4 + (b**2*EZ3**2 + b*c*EZ3*EZ4 + c**2*EZ4**2)*Log(EZ1)) - b**3*EZ3**3*f*Log(EZ3) + c**3*EZ4**3*f*Log(EZ4)
        cw(2)=b*c*EZ1*(half - b*EZ3 - c*EZ4 + (b**2*EZ3**2 + b*c*EZ3*EZ4 + c**2*EZ4**2)*Log(EZ1)) - b**3*EZ3**3*f*Log(EZ3) + c**3*EZ4**3*f*Log(EZ4)
        cw(3)=-(b**2*c*EZ1**2*(-1 + (2*b*EZ3 + c*EZ4)*Log(EZ1))) + b**2*EZ3**2*f*(1 + (2*b*EZ1 - EZ4*f)*Log(EZ3)) + c**2*EZ4**3*f**2*Log(EZ4)
        cw(4)=-(b*c**2*EZ1**2*(-1 + (b*EZ3 + 2*c*EZ4)*Log(EZ1))) + b**2*EZ3**3*f**2*Log(EZ3) - c**2*EZ4**2*f*(1 + (2*c*EZ1 + EZ3*f)*Log(EZ4))
    case(3) ! e2=e3
        cw(1)=-(a**2*c*EZ1**2*(-1 + (2*a*EZ2 + c*EZ4)*Log(EZ1))) + a**2*e*EZ2**2*(1 + (2*a*EZ1 - e*EZ4)*Log(EZ2)) + c**2*e**2*EZ4**3*Log(EZ4)
        cw(2)=a**3*c*EZ1**3*Log(EZ1) - a*e*EZ2*(half + a*EZ1 - e*EZ4 + (a**2*EZ1**2 - a*e*EZ1*EZ4 + e**2*EZ4**2)*Log(EZ2)) + c*e**3*EZ4**3*Log(EZ4)
        cw(3)=a**3*c*EZ1**3*Log(EZ1) - a*e*EZ2*(half + a*EZ1 - e*EZ4 + (a**2*EZ1**2 - a*e*EZ1*EZ4 + e**2*EZ4**2)*Log(EZ2)) + c*e**3*EZ4**3*Log(EZ4)
        cw(4)=a**2*c**2*EZ1**3*Log(EZ1) - a*e**2*EZ2**2*(1 + (a*EZ1 - 2*e*EZ4)*Log(EZ2)) - c*e**2*EZ4**2*(1 + (c*EZ1 + 2*e*EZ2)*Log(EZ4))
    case(4) ! e3=e4
        cw(1)=-(a*b**2*EZ1**2*(-1 + (a*EZ2 + 2*b*EZ3)*Log(EZ1))) + a**2*d**2*EZ2**3*Log(EZ2) - b**2*d*EZ3**2*(1 + (2*b*EZ1 + d*EZ2)*Log(EZ3))
        cw(2)=a**2*b**2*EZ1**3*Log(EZ1) - a*d**2*EZ2**2*(1 + (a*EZ1 - 2*d*EZ3)*Log(EZ2)) - b*d**2*EZ3**2*(1 + (b*EZ1 + 2*d*EZ2)*Log(EZ3))
        cw(3)=a*b**3*EZ1**3*Log(EZ1) - a*d**3*EZ2**3*Log(EZ2) + b*d*EZ3*(half + b*EZ1 + d*EZ2 + (b**2*EZ1**2 + b*d*EZ1*EZ2 + d**2*EZ2**2)*Log(EZ3))
        cw(4)=a*b**3*EZ1**3*Log(EZ1) - a*d**3*EZ2**3*Log(EZ2) + b*d*EZ3*(half + b*EZ1 + d*EZ2 + (b**2*EZ1**2 + b*d*EZ1*EZ2 + d**2*EZ2**2)*Log(EZ3))
    case(5) ! e1=e2 e3=e4
        cw(1)=-d - (3*d**2*EZ2)*half + 3*d**3*EZ2*EZ3 + 3*d**4*EZ2*EZ3**2*Log(EZ2/EZ3)
        cw(2)=-d - (3*d**2*EZ2)*half + 3*d**3*EZ2*EZ3 + 3*d**4*EZ2*EZ3**2*Log(EZ2/EZ3)
        cw(3)=d - (3*d**2*EZ3)*half - 3*d**3*EZ2*EZ3 + 3*d**4*EZ2**2*EZ3*Log(EZ3/EZ2)
        cw(4)=d - (3*d**2*EZ3)*half - 3*d**3*EZ2*EZ3 + 3*d**4*EZ2**2*EZ3*Log(EZ3/EZ2)
    case(6) ! e1=e2=e3
        cw(1)=f*third - (EZ4*f**2)*half + EZ4**2*f**3 + EZ4**3*f**4*Log(EZ4/EZ3)
        cw(2)=f*third - (EZ4*f**2)*half + EZ4**2*f**3 + EZ4**3*f**4*Log(EZ4/EZ3)
        cw(3)=f*third - (EZ4*f**2)*half + EZ4**2*f**3 + EZ4**3*f**4*Log(EZ4/EZ3)
        cw(4)=-f + (3*EZ3*f**2)*half - 3*EZ3*EZ4*f**3 + 3*EZ3*EZ4**2*f**4*Log(EZ3/EZ4)
    case(7) ! e2=e3=e4
        cw(1)=-a - (3*a**2*EZ2)*half + 3*a**3*EZ1*EZ2 + 3*a**4*EZ1**2*EZ2*Log(EZ2/EZ1)
        cw(2)=-a - (3*a**2*EZ2)*half + 3*a**3*EZ1*EZ2 + 3*a**4*EZ1**2*EZ2*Log(EZ2/EZ1)
        cw(3)=-a - (3*a**2*EZ2)*half + 3*a**3*EZ1*EZ2 + 3*a**4*EZ1**2*EZ2*Log(EZ2/EZ1)
        cw(4)=-a*third + (a**2*EZ1)*half - a**3*EZ1**2 + a**4*EZ1**3*Log(EZ2/EZ1)
    case(8) ! e1=e2=e3=e4
        cw(1)=0.25_r8/EZ1
        cw(2)=0.25_r8/EZ2
        cw(3)=0.25_r8/EZ3
        cw(4)=0.25_r8/EZ4
    end select
    ! Sort it back again
    weights(ind)=cw
end function

!> The Lambin & Vigneron tetrahedron method for the resolvent operator
module pure function lo_LV_tetrahedron_weights(ein,z,tol,sigma) result(weights)
    !> dispersion values at the corners of the tetrahedron
    real(r8), intent(in), dimension(4) :: ein
    !> energy to evaulate the weights at
    real(r8), intent(in) :: z
    !> tolerance for when something is zero
    real(r8), intent(in) :: tol
    !> sigma for when things are completely degenerate
    real(r8), intent(in) :: sigma
    !> integration weights
    real(r8), dimension(4) :: weights

    real(r8) :: EZ1,EZ2,EZ3,EZ4
    real(r8), dimension(4) :: ei,wt
    real(r8) :: Emin,Emax
    real(r8) :: E12,E13,E14,E23,E24,E34
    real(r8) :: a,b,c,d,e,f,ff0,ff1,ff2,ff3,gg0,gg1,gg2,gg3,hh0,hh1,hh2,hh3,ii0,ii1,ii2,ii3
    integer, dimension(4) :: ind
    integer :: whatkind

    ! Start by sorting the tetrahedron
    ind=sort_four_numbers(ein)
    ei=ein(ind)
    ! Min and max energy
    Emin=ei(1)
    Emax=ei(4)
    weights=0.0_r8
    wt=0.0_r8
    ! Possibility for an easy way out
    if ( z .le. Emin .or. z .ge. Emax ) then
        return
    endif
    whatkind=0
    ! Figure out what to do
    if ( Emax-Emin .lt. tol ) then
        ! completely degenerate, e1=e2=e3=e4=eps
        whatkind=1
    elseif ( z .lt. ei(2) ) then
        ! interval e1 < eps < e2
        whatkind=2
    elseif ( z .lt. ei(3) ) then
        ! interval e2 < eps < e3
        whatkind=3
    elseif ( z .lt. ei(4) ) then
        ! interval e3 < eps < e3
        whatkind=4
    endif

    ! Now the actual integration
    select case(whatkind)
    case(1) ! Gaussian or something
        EZ1=z-ei(1)
        EZ2=z-ei(2)
        EZ3=z-ei(3)
        EZ4=z-ei(4)
        wt(1)=lo_gauss(EZ1,0.0_r8,sigma)*0.25_r8
        wt(2)=lo_gauss(EZ2,0.0_r8,sigma)*0.25_r8
        wt(3)=lo_gauss(EZ3,0.0_r8,sigma)*0.25_r8
        wt(4)=lo_gauss(EZ4,0.0_r8,sigma)*0.25_r8
    case(2) ! e1<e<e2<e3<e4
        EZ1=z-ei(1)
        EZ2=z-ei(2)
        EZ3=z-ei(3)
        EZ4=z-ei(4)
        E12=ei(2)-ei(1)
        E13=ei(3)-ei(1)
        E14=ei(4)-ei(1)
        a=1.0_r8/E12
        b=1.0_r8/E13
        c=1.0_r8/E14
        wt(1)=a*b*c*EZ1**2*(-a*EZ2 - b*EZ3 - c*EZ4)
        wt(2)=a**2*b*c*EZ1**3
        wt(3)=a*b**2*c*EZ1**3
        wt(4)=a*b*c**2*EZ1**3
    case(3) ! e1<e2<e<e3<e4
        EZ1=z-ei(1)
        EZ2=z-ei(2)
        EZ3=z-ei(3)
        EZ4=z-ei(4)
        E13=ei(3)-ei(1)
        E14=ei(4)-ei(1)
        E23=ei(3)-ei(2)
        E24=ei(4)-ei(2)
        b=1.0_r8/E13
        c=1.0_r8/E14
        d=1.0_r8/E23
        e=1.0_r8/E24
        ff0=-b**2*EZ3
        ff2=-c**2*EZ4
        gg0=-d**2*EZ3
        gg2=-e**2*EZ4
        hh0=d**2*EZ2
        hh2=b**2*EZ1
        ii0=e**2*EZ2
        ii2=c**2*EZ1
        ff1=-c*d*EZ1*EZ3-d*e*EZ2*EZ3-c*e*EZ1*EZ4
        ff3=-b*d*EZ1*EZ3-b*e*EZ1*EZ4-d*e*EZ2*EZ4
        gg1=-b*c*EZ1*EZ3-b*e*EZ2*EZ3-c*e*EZ2*EZ4
        gg3=-b*d*EZ2*EZ3-b*c*EZ1*EZ4-c*d*EZ2*EZ4
        hh1=-b*c*EZ1*EZ3-b*e*EZ2*EZ3-c*e*EZ2*EZ4
        hh3=-c*d*EZ1*EZ3-d*e*EZ2*EZ3-c*e*EZ1*EZ4
        ii1=-b*d*EZ2*EZ3-b*c*EZ1*EZ4-c*d*EZ2*EZ4
        ii3=-b*d*EZ1*EZ3-b*e*EZ1*EZ4-d*e*EZ2*EZ4
        wt(1)=0.5_r8*(ff0*ff1+ff2*ff3)
        wt(2)=0.5_r8*(gg0*gg1+gg2*gg3)
        wt(3)=0.5_r8*(hh0*hh1+hh2*hh3)
        wt(4)=0.5_r8*(ii0*ii1+ii2*ii3)
    case(4) ! e1<e2<e3<e<e4
        EZ1=z-ei(1)
        EZ2=z-ei(2)
        EZ3=z-ei(3)
        EZ4=z-ei(4)
        E14=ei(4)-ei(1)
        E24=ei(4)-ei(2)
        E34=ei(4)-ei(3)
        c=1.0_r8/E14
        e=1.0_r8/E24
        f=1.0_r8/E34
        wt(1)=-(c**2*e*EZ4**3*f)
        wt(2)=-(c*e**2*EZ4**3*f)
        wt(3)=-(c*e*EZ4**3*f**2)
        wt(4)=c*e*EZ4**2*f*(c*EZ1 + e*EZ2 + EZ3*f)
    end select
    ! Sort it back again
    weights=wt(ind)
end function

!> The Lambin & Vigneron tetrahedron method for a Heaviside function
module pure function lo_LV_tetrahedron_heaviside(ein,z,tol) result(weight)
    !> dispersion values at the corners of the tetrahedron
    real(r8), intent(in), dimension(4) :: ein
    !> energy to evaulate the weights at
    real(r8), intent(in) :: z
    !> tolerance for when something is zero
    real(r8), intent(in) :: tol
    !> integration weights
    real(r8) :: weight

    real(r8) :: e1,e2,e3,e4
    integer :: whatkind

    ! First get the energies in order
    sorttet: block
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2

        if ( ein(1) <= ein(2) ) then
            low1 = 1
            high1 = 2
        else
            low1 = 2
            high1 = 1
        endif

        if ( ein(3) <= ein(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( ein(low1) <= ein(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( ein(high1) >= ein(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( ein(middle1) < ein(middle2) ) then
            e1=ein(lowest)
            e2=ein(middle1)
            e3=ein(middle2)
            e4=ein(highest)
        else
            e1=ein(lowest)
            e2=ein(middle2)
            e3=ein(middle1)
            e4=ein(highest)
        endif
    end block sorttet

    ! Figure out what to do
    weight=0.0_r8
    whatkind=0
    if ( z .le. e1 ) then
        ! Easy escape, weight is zero
        weight=0.0_r8
        return
    elseif ( e4-e1 .lt. tol ) then
        ! completely degenerate, e1=e2=e3=e4=z
        whatkind=1
    elseif ( z .lt. e2 ) then
        ! interval e1 < z < e2
        whatkind=2
    elseif ( z .lt. e3 ) then
        ! interval e2 < z < e3
        whatkind=3
    elseif ( z .lt. e4 ) then
        ! interval e3 < z < e4
        whatkind=4
    else
        ! interval z > e4
        whatkind=5
    endif

    ! Now the actual integration
    select case(whatkind)
    case(1) ! Gaussian or something? Not sure yet.
        weight=0.0_r8
        ! EZ1=z-ei(1)
        ! EZ2=z-ei(2)
        ! EZ3=z-ei(3)
        ! EZ4=z-ei(4)
        ! wt(1)=lo_gauss(EZ1,0.0_r8,sigma)*0.25_r8
        ! wt(2)=lo_gauss(EZ2,0.0_r8,sigma)*0.25_r8
        ! wt(3)=lo_gauss(EZ3,0.0_r8,sigma)*0.25_r8
        ! wt(4)=lo_gauss(EZ4,0.0_r8,sigma)*0.25_r8
    case(2) ! e1<z<e2<e3<e4
    interval2: block
        real(r8) :: E12,E13,E14
        real(r8) :: a,b,c
        real(r8) :: d1A,d1B,d1C

        E12=e2-e1
        E13=e3-e1
        E14=e4-e1
        a=1.0_r8/E12
        b=1.0_r8/E13
        c=1.0_r8/E14
        ! Coefficients
        d1A = (3*e1**2)*(a*b*c)
        d1B = -6*e1*(a*b*c)
        d1C = 3*(a*b*c)
        ! No weight from previous interval!
        ! Weight in this interval
        weight = d1A*(z-e1) + d1B*(z**2-e1**2)/2.0_r8 + d1C*(z**3-e1**3)/3.0_r8
    end block interval2
    case(3) ! e1<e2<z<e3<e4
    interval3: block
        real(r8) :: E13,E14,E23,E24
        real(r8) :: b,c,d,e
        real(r8) :: d2A,d2B,d2C

        E13=e3-e1
        E14=e4-e1
        E23=e3-e2
        E24=e4-e2
        b=1.0_r8/E13
        c=1.0_r8/E14
        d=1.0_r8/E23
        e=1.0_r8/E24

        d2A = -3*((e1 + e2)*(e3*e4) - e1*e2*(e3 + e4))*(b*c*d*e)
        d2B = 6*(-e1*e2 + e3*e4)*(b*c*d*e)
        d2C = 3*(e1 + e2 - e3 - e4)*(b*c*d*e)

        ! Weight from previous interval
        weight = b*c*(e2-e1)**2
        ! Weight from this interval
        weight = weight + d2A*(z-e2) + d2B*(z**2-e2**2)/2.0_r8 + d2C*(z**3-e2**3)/3.0_r8
    end block interval3
    case(4) ! e1<e2<e3<z<e4
    interval4: block
        real(r8) :: E14,E24,E34,c,e,f
        real(r8) :: d3A,d3B,d3C

        E14=e4-e1
        E24=e4-e2
        E34=e4-e3
        c=1.0_r8/E14
        e=1.0_r8/E24
        f=1.0_r8/E34
        d3A =  3*e4**2*(c*e*f)
        d3B = -6*e4*(c*e*f)
        d3C =  3*(c*e*f)
        ! Weight from previous intervals
        weight = (-e3**2 + e1*(e2 - e4) - e2*e4 + 2*e3*e4)*c*e
        weight = weight + d3A*(z-e3) + d3B*(z**2-e3**2)/2.0_r8 + d3C*(z**3-e3**3)/3.0_r8
    end block interval4
    case(5) ! e1<e2<e3<e4<z
        ! Simple enough!
        weight=1.0_r8
    end select
end function

!> The Lambin & Vigneron tetrahedron method for a Fermi function
module pure function lo_LV_tetrahedron_fermi(ein,z,temperature,tol) result(weight)
    !> dispersion values at the corners of the tetrahedron
    real(r8), intent(in), dimension(4) :: ein
    !> energy to evaulate the weights at
    real(r8), intent(in) :: z
    !> temperature
    real(r8), intent(in) :: temperature
    !> tolerance for when something is zero
    real(r8), intent(in) :: tol
    !> integration weights
    real(r8) :: weight

    real(r8) :: e1,e2,e3,e4,kbt,ikbt
    real(r8) :: d1A,d1B,d1C
    real(r8) :: d2A,d2B,d2C
    real(r8) :: d3A,d3B,d3C

    ! First get the energies in order
    sorttet: block
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2

        if ( ein(1) <= ein(2) ) then
            low1 = 1
            high1 = 2
        else
            low1 = 2
            high1 = 1
        endif

        if ( ein(3) <= ein(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( ein(low1) <= ein(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( ein(high1) >= ein(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( ein(middle1) < ein(middle2) ) then
            e1=ein(lowest)
            e2=ein(middle1)
            e3=ein(middle2)
            e4=ein(highest)
        else
            e1=ein(lowest)
            e2=ein(middle2)
            e3=ein(middle1)
            e4=ein(highest)
        endif
    end block sorttet

    ! The temperature factor
    kbt=temperature*lo_kb_Hartree
    ikbt=1.0_r8/kbt
    ! There could be cheap exits (36 sigma is where the Fermi function is 0 or 1 to 15 digits.)
    if ( z .lt. e1-36*kbt ) then
        weight=0.0_r8
        return
    endif
    if ( z .gt. e4+36*kbt ) then
        weight=1.0_r8
        return
    endif

    ! First get polynomial coefficients, for all three intervals
    if ( e2-e1 .gt. tol ) then
        interval2: block
            real(r8) :: E12,E13,E14
            real(r8) :: a,b,c
            E12=e2-e1
            E13=e3-e1
            E14=e4-e1
            a=1.0_r8/E12
            b=1.0_r8/E13
            c=1.0_r8/E14
            ! Coefficients
            d1A = (3*e1**2)*(a*b*c)
            d1B = -6*e1*(a*b*c)
            d1C = 3*(a*b*c)
        end block interval2
    else
        ! nothing in this interval
        d1A=0.0_r8
        d1B=0.0_r8
        d1C=0.0_r8
    endif
    if ( e3-e2 .gt. tol ) then
        interval3: block
            real(r8) :: E13,E14,E23,E24
            real(r8) :: b,c,d,e
            E13=e3-e1
            E14=e4-e1
            E23=e3-e2
            E24=e4-e2
            b=1.0_r8/E13
            c=1.0_r8/E14
            d=1.0_r8/E23
            e=1.0_r8/E24
            d2A = -3*((e1 + e2)*(e3*e4) - e1*e2*(e3 + e4))*(b*c*d*e)
            d2B = 6*(-e1*e2 + e3*e4)*(b*c*d*e)
            d2C = 3*(e1 + e2 - e3 - e4)*(b*c*d*e)
        end block interval3
    else
        ! nothing in this interval
        d2A=0.0_r8
        d2B=0.0_r8
        d2C=0.0_r8
    endif
    if ( e4-e3 .gt. tol ) then
        interval4: block
            real(r8) :: E14,E24,E34,c,e,f
            E14=e4-e1
            E24=e4-e2
            E34=e4-e3
            c=1.0_r8/E14
            e=1.0_r8/E24
            f=1.0_r8/E34
            d3A =  3*e4**2*(c*e*f)
            d3B = -6*e4*(c*e*f)
            d3C =  3*(c*e*f)
        end block interval4
    else
        ! nothing in this interval
        d3A=0.0_r8
        d3B=0.0_r8
        d3C=0.0_r8
    endif

    ! Neat. Now get try to evaluate this awful function.
    getwts: block
        real(r8) :: a11,a12,a13
        real(r8) :: a21,a22,a23
        real(r8) :: a31,a32,a33
        real(r8) :: xa,xb,xf

        xf=z*ikbt
        ! First interval
        xa=e1*ikbt
        xb=e2*ikbt
        a11=kbt*( (xb-xa)+log( (1+exp(xa-xf))/(1+exp(xb-xf)) ) )
        a12=kbt**2*(&
            xa*Log(1 + exp(-xa + xf)) - xb*Log(1 + exp(-xb + xf))&
            - lo_polylog_Li2(-exp(-xa + xf)) + lo_polylog_Li2(-exp(-xb + xf))&
        )
        a13=kbt**3*(&
            xa**2*Log(1 + exp(-xa + xf)) - xb**2*Log(1 + exp(-xb + xf)) - &
            2*xa*lo_polylog_Li2(-exp(-xa + xf)) + 2*xb*lo_polylog_Li2(-exp(-xb + xf)) - &
            2*lo_polylog_Li3(-exp(-xa + xf)) + 2*lo_polylog_Li3(-exp(-xb + xf))&
        )
        ! Second interval
        xa=e2*ikbt
        xb=e3*ikbt
        a21=kbt*( (xb-xa)+log( (1+exp(xa-xf))/(1+exp(xb-xf)) ) )
        a22=kbt**2*(&
            xa*Log(1 + exp(-xa + xf)) - xb*Log(1 + exp(-xb + xf))&
            - lo_polylog_Li2(-exp(-xa + xf)) + lo_polylog_Li2(-exp(-xb + xf))&
        )
        a23=kbt**3*(&
            xa**2*Log(1 + exp(-xa + xf)) - xb**2*Log(1 + exp(-xb + xf)) - &
            2*xa*lo_polylog_Li2(-exp(-xa + xf)) + 2*xb*lo_polylog_Li2(-exp(-xb + xf)) - &
            2*lo_polylog_Li3(-exp(-xa + xf)) + 2*lo_polylog_Li3(-exp(-xb + xf))&
        )
        ! Third interval
        xa=e3*ikbt
        xb=e4*ikbt
        a31=kbt*( (xb-xa)+log( (1+exp(xa-xf))/(1+exp(xb-xf)) ) )
        a32=kbt**2*(&
            xa*Log(1 + exp(-xa + xf)) - xb*Log(1 + exp(-xb + xf))&
            - lo_polylog_Li2(-exp(-xa + xf)) + lo_polylog_Li2(-exp(-xb + xf))&
        )
        a33=kbt**3*(&
            xa**2*Log(1 + exp(-xa + xf)) - xb**2*Log(1 + exp(-xb + xf)) - &
            2*xa*lo_polylog_Li2(-exp(-xa + xf)) + 2*xb*lo_polylog_Li2(-exp(-xb + xf)) - &
            2*lo_polylog_Li3(-exp(-xa + xf)) + 2*lo_polylog_Li3(-exp(-xb + xf))&
        )
        ! Add together to get weights
        weight= a11*d1A + a12*d1B + a13*d1C + a21*d2A + a22*d2B + a23*d2C + a31*d3A + a32*d3B + a33*d3C
    end block getwts
end function

!> The smearing parameter for adaptive gaussian smearing
module pure function smearingparameter(qp,gradient,sigma,adaptiveparameter) result(w)
    !
    ! I add a small safety feature, the sigma is my usual sensible default
    ! and I don't want the adaptive one to get too far from this since that
    ! can create some strange spikes
    !
    !> The q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> The gradient at this point
    real(r8), dimension(3), intent(in) :: gradient
    !> The base smearing
    real(r8), intent(in) :: sigma
    !> The scaling factor
    real(r8), intent(in) :: adaptiveparameter
    !> the resulting sigma
    real(r8) :: w
    !
    integer :: i
    !
    w=0.0_r8
    do i=1,3
        w=w+dot_product(abs(gradient),qp%scaledrecbasis(:,i))**2
    enddo
    w=sqrt(w*0.083333333333_r8)*adaptiveparameter
    ! I make sure the smearing does not get too large or small
    w=min(w,adaptiveparameter*sigma*4.0_r8)
    w=max(w,adaptiveparameter*sigma*0.25_r8)
end function

!> The smearing parameter for adaptive gaussian smearing
module pure function adaptive_sigma(radius,gradient,default_sigma,scale) result(sigma)
    !> radius of this q-point that give appropriate volume
    real(r8), intent(in) :: radius
    !> gradient at this point
    real(r8), dimension(3), intent(in) :: gradient
    !> baseline sensible smearing
    real(r8), intent(in) :: default_sigma
    !> scaling factor
    real(r8), intent(in) :: scale
    !> resulting smearing parameter
    real(r8) :: sigma

    real(r8), parameter :: largefactor=4.0_r8              ! largest multiple of baseline
    real(r8), parameter :: smallfactor=1.0_r8/4.0_r8       ! smallest multiple of baseline
    real(r8), parameter :: prefactor=lo_twopi/sqrt(2.0_r8) ! prefactor that takes care of 2pi and stuff.

    sigma=scale*prefactor*radius*norm2(gradient)
    sigma=max(default_sigma*smallfactor*scale,sigma)
    sigma=min(default_sigma*largefactor*scale,sigma)
end function

!> Integration weights for a tetrahedron Adapted from Matthieus code. And by adapted I mean copy-pasted. I add a small delta, this is to get rid of the pathological case of probing a completely degenerate tetrahedron. The original is from @cite Lehmann1972, optionally with the corrections from @cite Blochl1994
module subroutine lo_integration_weights_for_one_tetrahedron(tet,c,eps,sigma,tol,weights,blochlcorrections,hweights)
    !> The tetrahedron
    type(lo_qptetrahedron), intent(in) :: tet
    !> Energies at the corners
    real(r8), dimension(4), intent(in) :: c
    !> The probing energy
    real(r8), intent(in) :: eps
    !> Smearing parameter for the pathological cases
    real(r8), intent(in) :: sigma
    !> Tolerance for when a tetrahedron is degenerate
    real(r8), intent(in) :: tol
    !> The integration weights
    real(r8), intent(out), dimension(4) :: weights
    !> Should Blöchl corrections be added
    logical, intent(in), optional :: blochlcorrections
    !> Weights for Heaviside function
    real(r8), intent(out), dimension(4), optional :: hweights

    ! mine
    real(r8), dimension(4) :: eigen_1tetra
    real(r8), dimension(4) :: dtweightde,tweight
    real(r8) :: delta
    integer, dimension(4) :: tetorder
    integer :: bracket, i
    logical :: corrections
    ! stolen ones
    real(r8) :: cc,cc1,cc2,cc3,dcc1de,dcc2de,dcc3de,dccde
    real(r8) :: epsilon21,epsilon31,epsilon32,epsilon41,epsilon42,epsilon43
    real(r8) :: inv_epsilon21,inv_epsilon31
    real(r8) :: inv_epsilon32,inv_epsilon41,inv_epsilon42,inv_epsilon43
    real(r8) :: deleps1, deleps2, deleps3, deleps4
    real(r8) :: invepsum, cc_pre, dccde_pre
    real(r8) :: cc1_pre, cc2_pre, cc3_pre
    real(r8) :: cc_tmp, dccde_tmp
    real(r8) :: dcc1de_pre, dcc2de_pre, dcc3de_pre
    real(r8) :: volconst_mult

    ! Blöchl corrections?
    if ( present(blochlcorrections) ) then
        corrections=blochlcorrections
    else
        corrections=.false.
    endif
    ! then the weight. No idea what the factor 4 is here.
    volconst_mult=tet%integration_weight*0.25_r8

    ! Sort values according to "energy"
    tetorder=sort_four_numbers(c)
    eigen_1tetra=c(tetorder)

    ! The threshold of degeneracy for switching to gaussian integration
    delta=sigma/50.0_r8

    ! Start with zero weights
    weights=0.0_r8
    if ( present(hweights) ) hweights=0.0_r8

    ! Figure out what bracket I am in
    if ( sum(abs(eigen_1tetra-eps)) .lt. delta*4 ) then
        ! completely degenerate, e1=e2=e3=e4=eps
        bracket=1
    elseif ( eps .gt. eigen_1tetra(1) .and. eps .lt. eigen_1tetra(2) ) then
        ! interval e1 < eps < e2
        ! if this is too small, there will be not contributions
        if ( eps-eigen_1tetra(1) .lt. tol ) return
        bracket=2
        !
    elseif ( eps .gt. eigen_1tetra(2) .and. eps .lt. eigen_1tetra(3) ) then
        ! interval e2 < eps < e3
        ! if this is too small, there will be not contributions
        if ( eps-eigen_1tetra(1) .lt. tol ) return
        if ( eps-eigen_1tetra(2) .lt. tol ) return
        if ( eigen_1tetra(3)-eps .lt. tol ) return
        if ( eigen_1tetra(4)-eps .lt. tol ) return
        bracket=3
        !
    elseif ( eps .gt. eigen_1tetra(3) .and. eps .lt. eigen_1tetra(4) ) then
        ! interval e3 < eps < e3
        ! if this is too small, there will be not contributions
        if ( eigen_1tetra(4)-eps .lt. tol ) return
        bracket=4
        !
    elseif ( eps .lt. eigen_1tetra(1) ) then
        ! completely below
        bracket=0
        return
    else
        ! completely occupied
        bracket=0
        if ( present(hweights) ) then
            hweights=volconst_mult
        endif
        return
    endif
    !
    ! eigen_1tetra holds the sorted energies on the corners.
    ! eps is the probing energy
    !
    ! Some constants that are always needed
    epsilon21 = eigen_1tetra(2)-eigen_1tetra(1)
    epsilon31 = eigen_1tetra(3)-eigen_1tetra(1)
    epsilon41 = eigen_1tetra(4)-eigen_1tetra(1)
    epsilon32 = eigen_1tetra(3)-eigen_1tetra(2)
    epsilon42 = eigen_1tetra(4)-eigen_1tetra(2)
    epsilon43 = eigen_1tetra(4)-eigen_1tetra(3)
    inv_epsilon21 = 0.0_r8
    inv_epsilon31 = 0.0_r8
    inv_epsilon41 = 0.0_r8
    inv_epsilon32 = 0.0_r8
    inv_epsilon42 = 0.0_r8
    inv_epsilon43 = 0.0_r8

    ! Perhaps I should normalize everyting, these tols are
    ! quite meaningless in THz.
    if (epsilon21 > tol) inv_epsilon21=1.0_r8/epsilon21
    if (epsilon31 > tol) inv_epsilon31=1.0_r8/epsilon31
    if (epsilon41 > tol) inv_epsilon41=1.0_r8/epsilon41
    if (epsilon32 > tol) inv_epsilon32=1.0_r8/epsilon32
    if (epsilon42 > tol) inv_epsilon42=1.0_r8/epsilon42
    if (epsilon43 > tol) inv_epsilon43=1.0_r8/epsilon43

    ! Now the actual routine starts
    dtweightde=0.0_r8
    tweight=0.0_r8
    select case(bracket)
    case(1)
        ! e1=e2=e3=e4=eps
        ! This is the pathological case of a totally degenerate tetrahedron
        ! in the future I should perhaps do this as an adaptive gaussian
        ! instead of a flat, but it's extremely rare.
        do i=1,4
            dtweightde(i)=lo_gauss(eps,eigen_1tetra(i),sigma)*volconst_mult
        enddo
        ! In the case of a Heaviside, not sure if I should add all of them,
        ! or half, or what? Add all of them for now.
        tweight=0.0_r8 !volconst_mult ! 0.0_r8
    case(2)
        !  interval e1 < eps < e2
        deleps1=eps-eigen_1tetra(1)
        cc_pre=volconst_mult*inv_epsilon21*inv_epsilon31*inv_epsilon41
        invepsum=inv_epsilon21+inv_epsilon31+inv_epsilon41
        dccde_pre=3.0_r8*volconst_mult*inv_epsilon21*inv_epsilon31*inv_epsilon41
        cc=cc_pre*deleps1*deleps1*deleps1
        !
        ! Weights for Dirac delta function
        !
        dccde = dccde_pre * deleps1*deleps1
        dtweightde(1) = dccde*(4.0_r8-deleps1*invepsum)-cc*invepsum
        dtweightde(2) = (dccde*deleps1+cc) * inv_epsilon21
        dtweightde(3) = (dccde*deleps1+cc) * inv_epsilon31
        dtweightde(4) = (dccde*deleps1+cc) * inv_epsilon41
        !
        ! Weights for Heaviside function
        !
        tweight(1) = cc*(4.0_r8-deleps1*invepsum)
        tweight(2) = cc*deleps1*inv_epsilon21
        tweight(3) = cc*deleps1*inv_epsilon31
        tweight(4) = cc*deleps1*inv_epsilon41
        !
        ! Blöchl corrections: Bin Xu
        ! Olle: replaced 8/40 with 0.2
        if ( corrections ) then
            ! Dirac
            dtweightde(1)=dtweightde(1)+0.2_r8*dccde_pre*deleps1*(epsilon21+epsilon31+epsilon41)
            dtweightde(2)=dtweightde(2)+0.2_r8*dccde_pre*deleps1*(-epsilon21+epsilon32+epsilon42)
            dtweightde(3)=dtweightde(3)+0.2_r8*dccde_pre*deleps1*(-epsilon31-epsilon32+epsilon43)
            dtweightde(4)=dtweightde(4)+0.2_r8*dccde_pre*deleps1*(-epsilon41-epsilon42-epsilon43)
            ! Heaviside
            tweight(1) = tweight(1) + dccde_pre*deleps1*deleps1*(epsilon21+epsilon31+epsilon41)*0.1_r8
            tweight(2) = tweight(2) + dccde_pre*deleps1*deleps1*(-epsilon21+epsilon32+epsilon42)*0.1_r8
            tweight(3) = tweight(3) + dccde_pre*deleps1*deleps1*(-epsilon31-epsilon32+epsilon43)*0.1_r8
            tweight(4) = tweight(4) + dccde_pre*deleps1*deleps1*(-epsilon41-epsilon42-epsilon43)*0.1_r8
        endif
    case(3)
        ! e2 < eps < e3
        ! This is when the tetrahedron is sliced into a squarish thing
        deleps1 = eps-eigen_1tetra(1)
        deleps2 = eps-eigen_1tetra(2)
        deleps3 = eigen_1tetra(3)-eps
        deleps4 = eigen_1tetra(4)-eps

        cc1_pre = volconst_mult*inv_epsilon31*inv_epsilon41
        cc2_pre = volconst_mult*inv_epsilon41*inv_epsilon32*inv_epsilon31
        cc3_pre = volconst_mult*inv_epsilon42*inv_epsilon32*inv_epsilon41
        dcc1de_pre = 2.0_r8*cc1_pre
        dcc2de_pre = cc2_pre
        dcc3de_pre = cc3_pre
        cc1 = cc1_pre*deleps1*deleps1
        cc2 = cc2_pre*deleps1*deleps2*deleps3
        cc3 = cc3_pre*deleps2*deleps2*deleps4
        dcc1de = dcc1de_pre * deleps1
        dcc2de = dcc2de_pre * (-deleps1*deleps2 + deleps1*deleps3 + deleps2*deleps3)
        dcc3de = dcc3de_pre * (2.0_r8*deleps2*deleps4 - deleps2*deleps2)
        ! Dirac
        dtweightde(1) = + dcc1de + ((dcc1de+dcc2de)*deleps3 -(cc1+cc2)) * inv_epsilon31 + ((dcc1de+dcc2de+dcc3de)*deleps4 -(cc1+cc2+cc3)) * inv_epsilon41
        dtweightde(2) = + dcc1de+dcc2de+dcc3de + ((dcc2de+dcc3de)*deleps3 -(cc2+cc3) ) * inv_epsilon32 + (dcc3de*deleps4  -cc3 ) * inv_epsilon42
        dtweightde(3) = + ((dcc1de+dcc2de)*deleps1 + (cc1+cc2) ) * inv_epsilon31 + ((dcc2de+dcc3de)*deleps2 + (cc2+cc3) ) * inv_epsilon32
        dtweightde(4) = + ((dcc1de+dcc2de+dcc3de)*deleps1 + (cc1+cc2+cc3) ) * inv_epsilon41 + (dcc3de*deleps2 + cc3) * inv_epsilon42
        ! Heaviside
        tweight(1) = cc1 + (cc1+cc2)*deleps3*inv_epsilon31 + (cc1+cc2+cc3)*deleps4*inv_epsilon41
        tweight(2) = cc1+cc2+cc3+(cc2+cc3)*deleps3*inv_epsilon32 + cc3*deleps4*inv_epsilon42
        tweight(3) = (cc1+cc2)*deleps1*inv_epsilon31 + (cc2+cc3)*deleps2*inv_epsilon32
        tweight(4) = (cc1+cc2+cc3)*deleps1*inv_epsilon41 + cc3*deleps2*inv_epsilon42

        ! Blöchl corrections, bxu
        if ( corrections ) then
            ! Dirac
            dtweightde(1) = dtweightde(1) + 0.1_r8*cc1_pre*(6.0_r8-6.0_r8*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)*(epsilon21+epsilon31+epsilon41)
            dtweightde(2) = dtweightde(2) + 0.1_r8*cc1_pre*(6.0_r8-6.0_r8*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)*(-epsilon21+epsilon32+epsilon42)
            dtweightde(3) = dtweightde(3) + 0.1_r8*cc1_pre*(6.0_r8-6.0_r8*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)*(-epsilon31-epsilon32+epsilon43)
            dtweightde(4) = dtweightde(4) + 0.1_r8*cc1_pre*(6.0_r8-6.0_r8*(epsilon31+epsilon42)*deleps2*inv_epsilon32*inv_epsilon42)*(-epsilon41-epsilon42-epsilon43)
            ! Heaviside
            tweight(1) = tweight(1) + 0.1_r8*cc1_pre*(3.0_r8*epsilon21+6.0_r8*deleps2-3.0_r8*(epsilon31+epsilon42)*deleps2**2*inv_epsilon32*inv_epsilon42)*(epsilon21+epsilon31+epsilon41)
            tweight(2) = tweight(2) + 0.1_r8*cc1_pre*(3.0_r8*epsilon21+6.0_r8*deleps2-3.0_r8*(epsilon31+epsilon42)*deleps2**2*inv_epsilon32*inv_epsilon42)*(-epsilon21+epsilon32+epsilon42)
            tweight(3) = tweight(3) + 0.1_r8*cc1_pre*(3.0_r8*epsilon21+6.0_r8*deleps2-3.0_r8*(epsilon31+epsilon42)*deleps2**2*inv_epsilon32*inv_epsilon42)*(-epsilon31-epsilon32+epsilon43)
            tweight(4) = tweight(4) + 0.1_r8*cc1_pre*(3.0_r8*epsilon21+6.0_r8*deleps2-3.0_r8*(epsilon31+epsilon42)*deleps2**2*inv_epsilon32*inv_epsilon42)*(-epsilon41-epsilon42-epsilon43)
        endif
    case(4)
        ! interval e3 < eps < e4
        deleps4 = eigen_1tetra(4)-eps
        cc_pre = volconst_mult*inv_epsilon41*inv_epsilon42*inv_epsilon43
        invepsum = inv_epsilon41+inv_epsilon42+inv_epsilon43
        dccde_pre = -3.0_r8*cc_pre
        cc = cc_pre * deleps4*deleps4*deleps4
        cc_tmp = cc * deleps4
        dccde = dccde_pre * deleps4*deleps4
        dccde_tmp = -dccde*deleps4 + cc
        ! Dirac
        dtweightde(1) =  dccde_tmp * inv_epsilon41
        dtweightde(2) =  dccde_tmp * inv_epsilon42
        dtweightde(3) =  dccde_tmp * inv_epsilon43
        dtweightde(4) = -dccde*4.0_r8 - dccde_tmp*invepsum
        ! Heaviside
        tweight(1) =  volconst_mult - cc_tmp*inv_epsilon41
        tweight(2) =  volconst_mult - cc_tmp*inv_epsilon42
        tweight(3) =  volconst_mult - cc_tmp*inv_epsilon43
        tweight(4) =  volconst_mult - cc*4.0_r8 + cc_tmp*invepsum

        ! Blöchl corrections
        ! Olle: replaced 12/40 with 0.6
        if ( corrections ) then
            ! Dirac
            dtweightde(1) = dtweightde(1)-0.6_r8*cc_pre*deleps4*(epsilon21+epsilon31+epsilon41)
            dtweightde(2) = dtweightde(2)-0.6_r8*cc_pre*deleps4*(-epsilon21+epsilon32+epsilon42)
            dtweightde(3) = dtweightde(3)-0.6_r8*cc_pre*deleps4*(-epsilon31-epsilon32+epsilon43)
            dtweightde(4) = dtweightde(4)-0.6_r8*cc_pre*deleps4*(-epsilon41-epsilon42-epsilon43)
            ! Heaviside
            tweight(1) = tweight(1) + 0.6_r8*cc_pre*deleps4*deleps4*(epsilon21+epsilon31+epsilon41)
            tweight(2) = tweight(2) + 0.6_r8*cc_pre*deleps4*deleps4*(-epsilon21+epsilon32+epsilon42)
            tweight(3) = tweight(3) + 0.6_r8*cc_pre*deleps4*deleps4*(-epsilon31-epsilon32+epsilon43)
            tweight(4) = tweight(4) + 0.6_r8*cc_pre*deleps4*deleps4*(-epsilon41-epsilon42-epsilon43)
        endif
    case default
        ! do nothing
        dtweightde=0.0_r8
        tweight=0.0_r8
    end select

    weights(tetorder)=dtweightde
    if ( present(hweights) ) then
        hweights(tetorder)=tweight
    endif
end subroutine

!> Pretty fast way of sorting four numbers. Stole it from stack overflow.
pure function sort_four_numbers(n) result(ind)
    !> the four numbers
    real(r8), dimension(4), intent(in) :: n
    !> the output order so that n(ind) is sorted
    integer, dimension(4) :: ind
    !
    integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2

    if ( n(1) <= n(2) ) then
        low1 = 1
        high1 = 2
    else
        low1 = 2
        high1 = 1
    endif

    if ( n(3) <= n(4) ) then
        low2 = 3
        high2 = 4
    else
        low2 = 4
        high2 = 3
    endif

    if ( n(low1) <= n(low2) ) then
        lowest = low1
        middle1 = low2
    else
        lowest = low2
        middle1 = low1
    endif

    if ( n(high1) >= n(high2) ) then
        highest = high1
        middle2 = high2
    else
        highest = high2
        middle2 = high1
    endif

    if ( n(middle1) < n(middle2) ) then
        ind=[lowest,middle1,middle2,highest]
    else
        ind=[lowest,middle2,middle1,highest]
    endif
end function

!> Calculate integration weights based on a Voronoi tesselation
module subroutine voronoi_integration_weights(qp,p,mw,mem,integration_error,verbosity)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(inout) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> integration error
    real(r8), intent(out) :: integration_error
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_distancetable) :: dt
    integer, dimension(:), allocatable :: map_to_unique
    integer :: nfull
    real(r8), dimension(:,:), allocatable :: ir,ar
    real(r8) :: cutoff
    real(r8) :: timer,t0,t1

    ! Start timer
    timer=walltime()
    t0=timer
    t1=timer

    init: block
        real(r8), parameter :: degentol=1E-18_r8
        type(lo_verletbox) :: vb
        real(r8), dimension(:,:), allocatable :: dr
        real(r8), dimension(3,4) :: tetr
        real(r8), dimension(3) :: v0
        real(r8) :: f0
        integer, dimension(:), allocatable :: di
        integer, dimension(3) :: boxdim
        integer :: ikp,jkp,itet,icrn,jcrn,l
        integer :: bi,bj,bk,ci,cj,ck,ei,ej,ek,i

        ! How many points should this rank deal with?
        dt%np=0
        do ikp=1,qp%n_irr_point
            if ( mod(ikp,mw%n) .ne. mw%r ) cycle
            dt%np=dt%np+1
        enddo

        ! Here's an annoying thing: there is a possibility that some of the points
        ! on the full grid are redundant. That screws with everything, have to deal
        ! with that in some sensible way. Work out which the unique points are.
        ! In a pretty convoluted way, since there are many points.
        call mem%allocate(dr,[3,qp%n_full_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(di,qp%n_full_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(map_to_unique,qp%n_full_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        dr=0.0_r8
        di=0
        map_to_unique=0
        do ikp=1,qp%n_full_point
            dr(:,ikp)=lo_clean_fractional_coordinates( matmul(p%inv_reciprocal_latticevectors,qp%ap(ikp)%r) )
        enddo
        boxdim=vb%boxdim(dr,5000,1E-10_r8)
        call vb%generate(dr,boxdim,mem)

        ! Return the unique points
        map_to_unique=0
        di=0
        l=0
        do ikp=1,qp%n_full_point
            if ( map_to_unique(ikp) .ne. 0 ) then
                ! this is not a unique point
                cycle
            else
                ! Index unique points
                l=l+1
                di(ikp)=l
            endif
            call vb%boxind(dr(:,ikp),bi,bj,bk)
            do ci=bi-1,bi+1
            do cj=bj-1,bj+1
            do ck=bk-1,bk+1
                ei=lo_index_in_periodic_array(ci,vb%nx)
                ej=lo_index_in_periodic_array(cj,vb%ny)
                ek=lo_index_in_periodic_array(ck,vb%nz)
                do i=1,vb%box(ei,ej,ek)%n
                    jkp=vb%box(ei,ej,ek)%ind(i)
                    if ( jkp .lt. ikp ) cycle
                    v0=lo_clean_fractional_coordinates( dr(:,ikp)-dr(:,jkp)+0.5_r8 )-0.5_r8
                    f0=v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3)
                    if ( f0 .lt. degentol ) then
                        map_to_unique(jkp)=l
                    endif
                enddo
            enddo
            enddo
            enddo
        enddo
        nfull=l
        call mem%allocate(ar,[3,nfull],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ar=0.0_r8
        l=0
        do ikp=1,qp%n_full_point
            if ( di(ikp) .eq. 0 ) cycle
            l=l+1
            ar(:,l)=dr(:,ikp)
        enddo
        ! Small sanity check
        if ( l .ne. nfull ) then
            call lo_stop_gracefully(['Clearly I do not understand points'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        endif

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... preparing weights, reduced from ',tochar(qp%n_full_point),' to ',tochar(nfull),' points'
            write(lo_iou,*) '... pre-processed points (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Some workspace
        if ( dt%np .gt. 0 ) then
            call mem%allocate(ir,[3,dt%np],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            ir=0.0_r8
            allocate(dt%particle(dt%np))
        endif

        ! store coordinates of irreducible point
        l=0
        do ikp=1,qp%n_irr_point
            if ( mod(ikp,mw%n) .ne. mw%r ) cycle
            l=l+1
            ir(:,l)=lo_clean_fractional_coordinates( matmul(p%inv_reciprocal_latticevectors,qp%ip(ikp)%r) )
        enddo

        ! Stuff the full set of points into Verlet boxes. But to do that I need an estimate
        ! of the neighbour distance. This is confusing. A safe bet should be tetrahedron
        ! edges, and some multiple of that. In fractional coordinates.
        cutoff=0.0_r8
        do itet=1,qp%n_irr_tet
            if ( mod(itet,mw%n) .ne. mw%r ) cycle
            ! Collect the tetrahedron
            v0=dr(:,qp%it(itet)%full_index(1))
            do icrn=1,4
                ikp=qp%it(itet)%full_index(icrn)
                tetr(:,icrn)=lo_clean_fractional_coordinates( dr(:,ikp)-v0 +0.5_r8 )-0.5_r8
                tetr(:,icrn)=matmul(p%reciprocal_latticevectors,tetr(:,icrn))
            enddo

            ! Get the largest distance?
            do icrn=1,4
            do jcrn=icrn+1,4
                cutoff=max(cutoff,norm2(tetr(:,icrn)-tetr(:,jcrn)))
            enddo
            enddo
        enddo
        ! And sync across ranks
        call mw%allreduce('max',cutoff)

        ! Should be safe to slightly more than double this? Will update if it breaks.
        !cutoff=cutoff*1.800001_r8
        cutoff=cutoff*2.60001_r8

        ! Intermediate cleanup
        call vb%destroy(mem)
        call mem%deallocate(dr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... transformed points (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block init

    ! Specialized version of calculate distance table.
    distancetable: block
        type(lo_verletbox) :: vb
        real(r8), dimension(:,:), allocatable :: dr
        real(r8), dimension(:), allocatable :: ds
        real(r8), dimension(3,3) :: m0
        real(r8), dimension(3) :: v0
        real(r8) :: f0,f1,cutoffsq
        integer, dimension(:), allocatable :: di
        integer, dimension(3) :: boxdim,bdi
        integer :: iter,ia,ja,ka,ipt,jpt,i,ctr
        integer :: bi,bj,bk,ci,cj,ck,ei,ej,ek

        ! First try to stuff into boxes.
        boxdim=1
        do iter=1,1000
            ! which dimension should I increase?
            f0=0.0_r8
            ja=0
            do ia=1,3
                bdi=boxdim
                bdi(ia)=bdi(ia)+1
                do ka=1,3
                    m0(:,ka)=p%reciprocal_latticevectors(:,ka)/bdi(ka)
                enddo
                f1=lo_inscribed_sphere_in_box(m0)
                if ( f1 .gt. f0 ) then
                    ja=ia
                    f0=f1
                endif
            enddo
            ! Check if we have divided enough?
            if ( f0 .lt. cutoff ) then
                ! yup, what we have is enough
                exit
            else
                ! Increment the number of divisions
                boxdim(ja)=boxdim(ja)+1
            endif
        enddo
        if ( maxval(boxdim) .le. 3 ) then
            ! No point in using Verlet boxes
            boxdim=-1
        else
            ! Makes sense to use Verlet boxes!
            call vb%generate(ar,boxdim,mem)
        endif

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... decided on box division: ',tochar(boxdim),' (',tochar(t1-t0),'s)'
            t0=t1
        endif

        cutoffsq=cutoff**2
        if ( verbosity .gt. 0 ) call lo_progressbar_init()

        ! So, now proceed to calculate the distance table, with or without Verlet boxes.
        do ipt=1,dt%np
            ! Count neighbours for this point
            ctr=0
            if ( boxdim(1) .lt. 0 ) then
                do jpt=1,nfull
                    v0=lo_clean_fractional_coordinates( ar(:,jpt)-ir(:,ipt)+0.5_r8 )-0.5_r8
                    v0=matmul(p%reciprocal_latticevectors,v0)
                    f0=v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3)
                    if ( f0 .lt. 1E-14_r8 ) cycle
                    if ( f0 .lt. cutoffsq ) then
                        ctr=ctr+1
                    endif
                enddo
            else
                ! box index of this point
                call vb%boxind(ir(:,ipt),bi,bj,bk)
                do ci=bi-1,bi+1
                do cj=bj-1,bj+1
                do ck=bk-1,bk+1
                    ei=lo_index_in_periodic_array(ci,vb%nx)
                    ej=lo_index_in_periodic_array(cj,vb%ny)
                    ek=lo_index_in_periodic_array(ck,vb%nz)
                    do i=1,vb%box(ei,ej,ek)%n
                        jpt=vb%box(ei,ej,ek)%ind(i)
                        v0=lo_clean_fractional_coordinates( ar(:,jpt)-ir(:,ipt)+0.5_r8 )-0.5_r8
                        v0=matmul(p%reciprocal_latticevectors,v0)
                        f0=v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3)
                        if ( f0 .lt. 1E-14_r8 ) cycle
                        if ( f0 .lt. cutoffsq ) then
                            ctr=ctr+1
                        endif
                    enddo
                enddo
                enddo
                enddo
            endif

            ! Space for intermediate storage
            if ( ctr .gt. 0 ) then
                call mem%allocate(dr,[3,ctr],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                call mem%allocate(ds,ctr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                call mem%allocate(di,ctr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                dr=0.0_r8
                ds=0.0_r8
                di=0
            else
                call lo_stop_gracefully(['No neighbours for this point, must be wrong'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            endif

            ! Store intermediate
            ctr=0
            if ( boxdim(1) .lt. 0 ) then
                do jpt=1,nfull
                    v0=lo_clean_fractional_coordinates( ar(:,jpt)-ir(:,ipt)+0.5_r8 )-0.5_r8
                    v0=matmul(p%reciprocal_latticevectors,v0)
                    f0=v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3)
                    if ( f0 .lt. 1E-14_r8 ) cycle
                    if ( f0 .lt. cutoffsq ) then
                        ctr=ctr+1
                        dr(:,ctr)=v0
                        ds(ctr)=f0
                    endif
                enddo
            else
                ! box index of this point
                call vb%boxind(ir(:,ipt),bi,bj,bk)
                do ci=bi-1,bi+1
                do cj=bj-1,bj+1
                do ck=bk-1,bk+1
                    ei=lo_index_in_periodic_array(ci,vb%nx)
                    ej=lo_index_in_periodic_array(cj,vb%ny)
                    ek=lo_index_in_periodic_array(ck,vb%nz)
                    do i=1,vb%box(ei,ej,ek)%n
                        jpt=vb%box(ei,ej,ek)%ind(i)
                        v0=lo_clean_fractional_coordinates( ar(:,jpt)-ir(:,ipt)+0.5_r8 )-0.5_r8
                        v0=matmul(p%reciprocal_latticevectors,v0)
                        f0=v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3)
                        if ( f0 .lt. 1E-14_r8 ) cycle
                        if ( f0 .lt. cutoffsq ) then
                            ctr=ctr+1
                            dr(:,ctr)=v0
                            ds(ctr)=f0
                        endif
                    enddo
                enddo
                enddo
                enddo
            endif

            ! Sort neighbours by distance
            call lo_qsort(ds,di)

            ! Store clean and neat information (sorted by distance)
            dt%particle(ipt)%n=ctr+1
            allocate(dt%particle(ipt)%d(ctr+1) )
            allocate(dt%particle(ipt)%v(3,ctr+1) )
            dt%particle(ipt)%d=0.0_r8
            dt%particle(ipt)%v=0.0_r8
            do i=2,dt%particle(ipt)%n
                dt%particle(ipt)%d(i)=sqrt(ds(i-1))
                dt%particle(ipt)%v(:,i)=dr(:,di(i-1))
            enddo

            ! And cleanup
            call mem%deallocate(dr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(ds,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

            ! And report?
            if ( verbosity .gt. 0 ) then
                t1=walltime()
                call lo_progressbar(' ... distancetable',ipt,dt%np,t1-t0)
            endif
        enddo
        t0=t1

        ! Intermediate cleanup, won't need the Verlet boxes anymore.
        call vb%destroy(mem)
    end block distancetable

    ! Now, finally, calculate the Voronoi cell for each point, and use the volume
    ! of that cell as the integration weight.
    vorocells: block
        type(lo_voronoi_cell) :: vc
        real(r8), dimension(:), allocatable :: wts,wtf
        real(r8) :: f0
        integer, dimension(:), allocatable :: multiplicity,di
        integer :: ipt,jpt

        ! Figure out the multiplicity thing.
        call mem%allocate(multiplicity,qp%n_full_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(di,nfull,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        multiplicity=0
        di=0
        do ipt=1,qp%n_full_point
            jpt=map_to_unique(ipt)
            di(jpt)=di(jpt)+1
        enddo
        do ipt=1,qp%n_full_point
            jpt=map_to_unique(ipt)
            multiplicity(ipt)=di(jpt)
        enddo

        ! Dummy space for weights
        call mem%allocate(wts,qp%n_irr_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(wtf,qp%n_full_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        wts=0.0_r8
        wtf=0.0_r8

        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        ipt=0
        do jpt=1,qp%n_irr_point
            if ( mod(jpt,mw%n) .ne. mw%r ) cycle
            ipt=ipt+1
            ! Generate the voronoi diagram
            call vc%generate(dt%particle(ipt),cutoff,cutoff*1E-7_r8,mem)
            ! Store weights?
            wts(jpt)=vc%volume
            if ( verbosity .gt. 0 .and. ipt .lt. dt%np ) then
                call lo_progressbar(' ... voronoi tesselation',ipt,dt%np,walltime()-t0)
            endif
        enddo
        ! sync
        call mw%allreduce('sum',wts)

        ! Before we deal with multiplicities, make sure I store away the relevant radius per point
        do ipt=1,qp%n_irr_point
            f0=(3.0_r8*wts(ipt)/4.0_r8/lo_pi)**(1.0_r8/3.0_r8)
            qp%ip(ipt)%radius=f0
        enddo
        do ipt=1,qp%n_full_point
            jpt=qp%ap(ipt)%irreducible_index
            f0=(3.0_r8*wts(jpt)/4.0_r8/lo_pi)**(1.0_r8/3.0_r8)
            qp%ap(ipt)%radius=f0
        enddo

        ! Take care of multiplicity
        do ipt=1,qp%n_full_point
            jpt=qp%ap(ipt)%irreducible_index
            wtf(ipt)=wts(jpt)/multiplicity(ipt)
        enddo

        ! Check if things add neatly to 1
        integration_error=abs(sum(wtf)*p%volume-1.0_r8)
        if ( integration_error .gt. 1E-11_r8 ) then
            if ( mw%talk ) then
                write(*,*) 'WARNING, bad integration weights:',sum(wtf)*p%volume
            endif
        endif

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            call lo_progressbar(' ... voronoi tesselation',dt%np,dt%np,t1-t0)
            write(*,*) '... sum weights:',sum(wtf)*p%volume
            t0=t1
        endif

        ! Set the weights.
        wtf=wtf/sum(wtf)
        wts=0.0_r8
        do ipt=1,qp%n_full_point
            jpt=qp%ap(ipt)%irreducible_index
            wts(jpt)=wts(jpt)+wtf(ipt)
        enddo
        do ipt=1,qp%n_irr_point
            qp%ip(ipt)%integration_weight=wts(ipt)
        enddo
        do ipt=1,qp%n_full_point
            qp%ap(ipt)%integration_weight=wtf(ipt)
        enddo

        ! And cleanup
        call mem%deallocate(wts,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(wtf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(multiplicity,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block vorocells

    ! And final cleanup.
    call mem%deallocate(map_to_unique,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    if ( dt%np .gt. 0 ) then
        call mem%deallocate(ar,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ir,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    endif
end subroutine

!> polylogarithm guy, defined for x<1
elemental function lo_polylog_Li2(x) result(val)
    implicit none
    !> argument
    real(r8), intent(in) :: x
    !> value
    real(r8) :: val

    ! some constants
    real(r8), parameter :: zero= 0.0_r8
    real(r8), parameter :: one = 1.0_r8
    real(r8), parameter :: half= 0.5_r8
    real(r8), parameter :: malf= -0.5_r8
    real(r8), parameter :: mone= -1.0_r8
    real(r8), parameter :: mtwo= -2.0_r8
    real(r8), parameter :: pi3 = 3.289868133696453_r8
    real(r8), parameter :: pi6 = 1.644934066848226_r8
    real(r8), dimension(19), parameter :: ci=[&
    -0.0000000000000001_r8,&
     0.0000000000000009_r8,&
    -0.0000000000000061_r8,&
     0.0000000000000404_r8,&
    -0.0000000000002701_r8,&
     0.0000000000018226_r8,&
    -0.0000000000124433_r8,&
     0.0000000000861210_r8,&
    -0.0000000006057848_r8,&
     0.0000000043454506_r8,&
    -0.0000000319334127_r8,&
     0.0000002419518085_r8,&
    -0.0000019078495939_r8,&
     0.0000158841554188_r8,&
    -0.0001430418444234_r8,&
     0.0014575108406227_r8,&
    -0.0185884366501460_r8,&
     0.4097598753307711_r8,&
     0.4299669356081370_r8]

    real(r8) :: y,t,s,a
    real(r8) :: h,alfa,b0,b1,b2
    integer :: i

    ! function values known for certain arguments
    if( x .eq. one ) then
        val=pi6
        return
    else if( x .eq. mone ) then
        val=malf*pi6
        return
    end if

    ! evaluate properly
    t=-x
    if(t .le. mtwo) then
        y=mone/(one+t)
        s=one
        a=-pi3+half*(log(-t)**2-log(one+one/t)**2)
    else if(t .lt. mone) then
        y=mone-t
        s=mone
        a=log(-t)
        a=-pi6+a*(a+log(one+one/t))
    else if(t .le. malf) then
        y=(mone-t)/t
        s=one
        a=log(-t)
        a=-pi6+a*(malf*a+log(one+t))
    else if(t .lt. zero) then
        y=-t/(one+t)
        s=mone
        a=half*log(one+t)**2
    else if(t .le. one) then
        y=t
        s=one
        a=zero
    else
        y=one/t
        s=mone
        a=pi6+half*log(t)**2
    end if

    h=y+y-one
    alfa=h+h
    b1=zero
    b2=zero
    ! do i = 18,0,-1
    !     b0=c(i)+alfa*b1-b2
    !     b2=b1
    !     b1=b0
    ! enddo
    do i = 1,19 !18,0,-1
        b0=ci(i)+alfa*b1-b2
        b2=b1
        b1=b0
    enddo

    val=-(s*(b0-h*b2)+a)
end function

!> polylogarithm guy, defined for x<1 I think
elemental function lo_polylog_Li3(x) result(val)
    !> argument
    real(r8), intent(in) :: x
    !> value
    real(r8) :: val

    real(r8), parameter :: PI    = 3.141592653589793_r8
    real(r8), parameter :: PI2   = PI*PI
    real(r8), parameter :: zeta2 = 1.644934066848226_r8
    real(r8), parameter :: zeta3 = 1.202056903159594_r8
    real(r8), parameter, dimension(18) :: bf=[&
        1.0_r8,&
        -3.0_r8/8.0_r8,&
        17.0_r8/216.0_r8,&
        -5.0_r8/576.0_r8,&
        1.296296296296296e-04_r8,&
        8.101851851851851e-05_r8,&
       -3.419357160853759e-06_r8,&
       -1.328656462585034e-06_r8,&
        8.660871756109851e-08_r8,&
        2.526087595532039e-08_r8,&
       -2.144694468364064e-09_r8,&
       -5.140110622012978e-10_r8,&
        5.249582114600829e-11_r8,&
        1.088775440663631e-11_r8,&
       -1.277939609449369e-12_r8,&
       -2.369824177308745e-13_r8,&
        3.104357887965462e-14_r8,&
        5.261758629912506e-15_r8]
    real(r8), parameter, dimension(7) :: cs=[&
         -3.472222222222222e-03_r8,&
          1.157407407407407e-05_r8,&
         -9.841899722852104e-08_r8,&
          1.148221634332745e-09_r8,&
         -1.581572499080917e-11_r8,&
          2.419500979252515e-13_r8,&
         -3.982897776989488e-15_r8]

    complex(r8) :: z,u,u2,u3,c0,c1,rest,lmz,cval
    real(r8) :: az,pz,lnz

    ! Check known values?
    if ( abs(x) .lt. lo_tiny ) then
        val=0.0_r8
        return
    elseif ( abs(x-1.0_r8) .lt. lo_tiny ) then
        val=zeta3
        return
    elseif ( abs(x+1.0_r8) .lt. lo_tiny ) then
        val=-0.75_r8*zeta3
        return
    elseif ( abs(x-0.5_r8) .lt. lo_tiny ) then
        val=0.5372131936080402_r8
    endif

    ! promote to complex? Just a little bit complex.
    z=cmplx(x,lo_tiny,r8)
    az=abs(z)
    pz=atan2(x,lo_tiny)
    lnz=log(az)

    if ( lnz*lnz+pz*pz .lt. 1.0_r8 ) then
        u  = log(z)
        u2 = u*u
        u3 = u*u2
        c0 = zeta3 + zeta2*u - u3/12.0_r8
        c1 = 0.25_r8 * (3.0_r8 - 2.0_r8*log(-u));
        ! evaluate
        cval= c0 +&
             u2 * (c1 +&
             u2 * (cs(1) +&
             u2 * (cs(2) +&
             u2 * (cs(3) +&
             u2 * (cs(4) +&
             u2 * (cs(5) +&
             u2 * (cs(6) +&
             u2 * (cs(7)))))))))
        val=real(cval,r8)
        return
    endif

    if ( az .le. 1 ) then
        u = -log(1.0_r8 - z)
        rest=0.0_r8
    else
        lmz = log(-z)
        u = -log(1.0_r8 - 1.0_r8/z)
        rest = -lmz*(lmz*lmz/6.0_r8 + PI2/6.0_r8)
    endif

    ! I wonder if this sums in the right order. Hmmm.
    cval = rest + &
          u * (bf(1) + &
          u * (bf(2) + &
          u * (bf(3) + &
          u * (bf(4) + &
          u * (bf(5) + &
          u * (bf(6) + &
          u * (bf(7) + &
          u * (bf(8) + &
          u * (bf(9) + &
          u * (bf(10) + &
          u * (bf(11) + &
          u * (bf(12) + &
          u * (bf(13) + &
          u * (bf(14) + &
          u * (bf(15) + &
          u * (bf(16) + &
          u * (bf(17) + &
          u * (bf(18)))))))))))))))))))
    val=real(cval,r8)
end function

end submodule
