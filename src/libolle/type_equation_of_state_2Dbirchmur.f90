#include "precompilerdefinitions"
submodule (type_equation_of_state) type_equation_of_state_2Dbirchmur
implicit none
contains

!> fit BM to data. Very robust.
module subroutine fit_2D_birch_murnaghan_parameters(eos,volume,energy,eta,verbosity)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(out) :: eos
    !> volumes
    real(flyt), dimension(:), intent(in) :: volume
    !> energy
    real(flyt), dimension(:), intent(in) :: energy
    !> c/a ratio
    real(flyt), dimension(:), intent(in) :: eta
    !> how much to talk
    integer, intent(in) :: verbosity

    integer :: n
    n=size(energy)
    if ( verbosity .gt. 0 ) then
        write(*,*) ''
        write(*,*) 'FITTING 2-D BIRCH-MURNAGHAN EQUATION OF STATE USING ',tochar(n),' POINTS'
    endif

    ! do a least squares fit
    getpar: block
        real(flyt), dimension(:,:), allocatable :: wA,wB
        real(flyt) :: x,y
        integer :: i

        n=size(volume)
        lo_allocate(wA(n,9))
        lo_allocate(wB(n,1))
        ! create the coefficient matrix and energy vector
        do i=1,n
            x=volume(i)**(-2.0_flyt/3.0_flyt)
            y=eta(i)
            wA(i,:)=[1.0_flyt,x,x**2,x**3,y,x*y,y**2,x*y**2,y*x**2]
            wB(i,1)=energy(i)
        enddo
        ! fit
        call lo_dgels(wA,wB)
        ! store values
        eos%A0=wB(1,1)
        eos%A1=wB(2,1)
        eos%A2=wB(3,1)
        eos%A3=wB(4,1)
        eos%B1=wB(5,1)
        eos%B2=wB(6,1)
        eos%B3=wB(7,1)
        eos%B4=wB(8,1)
        eos%B5=wB(9,1)
        if ( verbosity .gt. 0 ) write(*,*) '... got raw fit'
    end block getpar

    ! find a reasonable V0 and eta0
    findv0: block
        integer :: i,j,iter,ctr
        real(flyt) :: eta0,E0,V0,B0,B0p,x,y,E1
        real(flyt) :: C0,C1,C2,C3
        real(flyt) :: gradx,grady,h11,h12,h22,g111,g112,g122,f0

        eta0=lo_huge
        E0=lo_huge
        V0=lo_huge
        j=0
        do i=1,n
            if ( energy(i) .lt. E0 ) then
                E0=min(E0,energy(i))
                V0=volume(i)
                eta0=eta(i)
                j=i
            endif
        enddo

        f0=lo_huge
        f0=min(f0,abs(V0-minval(volume)))
        f0=min(f0,abs(V0-maxval(volume)))
        f0=min(f0,abs(eta0-minval(eta)))
        f0=min(f0,abs(eta0-maxval(eta)))
        if ( f0 .gt. lo_sqtol ) then
            ! Newton minimize the gradient from this guess
            E1=lo_huge
            ctr=0
            do iter=1,1000
                ! my volume parameter
                x=v0**(-2.0_flyt/3.0_flyt)
                y=eta0
                ! current energy
                E0  = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3 + eos%B1*y + eos%B2*x*y + eos%B3*y**2 + eos%B4*x*y**2 + eos%B5*y*x**2
                ! Gradient ( derivative in V and eta )
                gradx = (-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*y + 2*eos%B5*x*y + eos%B4*y**2))/3
                grady = eos%B1 + eos%B2*x + eos%B5*x**2 + 2*eos%B3*y + 2*eos%B4*x*y
                ! partial Hessian ( second derivative in V and eta )
                h11=(2*x**4*(5*eos%A1 + 14*eos%A2*x + 27*eos%A3*x**2 + 5*eos%B2*y + 14*eos%B5*x*y + 5*eos%B4*y**2))/9
                h22=2*(eos%B3 + eos%B4*x)
                h12=(-2*x**(2.5_flyt)*(eos%B2 + 2*eos%B5*x + 2*eos%B4*y))/3
                ! partial third derivatives
                g111=(-8*x**(5.5_flyt)*(10*eos%A1 + 35*eos%A2*x + 81*eos%A3*x**2 + 10*eos%B2*y + 35*eos%B5*x*y + 10*eos%B4*y**2))/27
                g112=(2*x**4*(5*eos%B2 + 14*eos%B5*x + 10*eos%B4*y))/9
                g122=(-4*eos%B4*x**(2.5_flyt))/3
                ! usual parameters
                B0=h11*V0
                B0p=-g111*V0**2/B0-1
                ! unusual
                C0=-h22/E0      ! Dimensionless eta-eta derivative
                C1=h12/B0       ! Dimensionless eta-V derivative
                C2=-g122/B0     ! Dimensionless V-eta-eta derivative
                C3=-V0*g112/B0  ! Dimensionless V-V-eta derivative
                ! The usual Birch-Murnaghan parameters, and the new ones
                if ( verbosity .gt. 0 ) write(*,"(1X,I3,9(2X,F12.6))") iter,E0,V0,eta0,B0,B0p,C0,C1,C2,C3
                if ( abs(E0-E1) .lt. 1E-20_flyt ) then
                    ! we are done
                    ctr=ctr+1
                    if ( ctr .gt. 2 ) exit
                else
                    ! iterate a little closer
                    V0=V0-gradx/h11
                    eta0=eta0-grady/h22
                    E1=E0
                    ctr=0
                endif
            enddo
        else
            if ( verbosity .gt. 0 ) write(*,*) 'WARNING: found the minimum energy at the edge of the data'
            ! just use what I got.
            x=v0**(-2.0_flyt/3.0_flyt)
            y=eta0
            ! current energy
            E0  = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3 + eos%B1*y + eos%B2*x*y + eos%B3*y**2 + eos%B4*x*y**2 + eos%B5*y*x**2
            ! Gradient ( derivative in V and eta )
            gradx = (-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*y + 2*eos%B5*x*y + eos%B4*y**2))/3
            grady = eos%B1 + eos%B2*x + eos%B5*x**2 + 2*eos%B3*y + 2*eos%B4*x*y
            ! partial Hessian ( second derivative in V and eta )
            h11=(2*x**4*(5*eos%A1 + 14*eos%A2*x + 27*eos%A3*x**2 + 5*eos%B2*y + 14*eos%B5*x*y + 5*eos%B4*y**2))/9
            h22=2*(eos%B3 + eos%B4*x)
            h12=(-2*x**(2.5_flyt)*(eos%B2 + 2*eos%B5*x + 2*eos%B4*y))/3
            ! third derivatives
            g111=(-8*x**(5.5_flyt)*(10*eos%A1 + 35*eos%A2*x + 81*eos%A3*x**2 + 10*eos%B2*y + 35*eos%B5*x*y + 10*eos%B4*y**2))/27
            g112=(2*x**4*(5*eos%B2 + 14*eos%B5*x + 10*eos%B4*y))/9
            g122=(-4*eos%B4*x**(2.5_flyt))/3
            ! usual parameters
            B0=h11*V0
            B0p=-g111*V0**2/B0-1
            ! unusual
            C0=-h22/E0      ! Dimensionless eta-eta derivative
            C1=h12/B0       ! Dimensionless eta-V derivative
            C2=-g122/B0     ! Dimensionless V-eta-eta derivative
            C3=-V0*g112/B0  ! Dimensionless V-V-eta derivative
        endif
        ! Calculate the fit parameters back again
        ! A0=(-8*E0*(-2 + C0*eta0**2) - 3*B0*(-18 + 3*B0p + 18*C1*eta0 - 6*C3*eta0 + 4*C2*eta0**2)*V0)/16
        ! A1=(3*B0*(9*B0p + 4*(-12 + 7*C1*eta0 - 3*C3*eta0 + C2*eta0**2))*V0**1.6666666666666667_flyt)/16
        ! A2=(-3*B0*(9*B0p + 10*C1*eta0 - 6*(7 + C3*eta0))*V0**2.3333333333333335_flyt)/16
        ! A3=(9*B0*(-4 + B0p)*V0**3)/16
        ! B1=C0*E0*eta0 + (3*B0*(9*C1 - 3*C3 + 4*C2*eta0)*V0)/8
        ! B2=(-3*B0*(7*C1 - 3*C3 + 2*C2*eta0)*V0**1.6666666666666667_flyt)/4
        ! B3=(-2*C0*E0 - 3*B0*C2*V0)/4
        ! B4=(3*B0*C2*V0**1.6666666666666667_flyt)/4
        ! B5=(3*B0*(5*C1 - 3*C3)*V0**2.3333333333333335_flyt)/8
        ! write(*,*) A0,eos%A0
        ! write(*,*) A1,eos%A1
        ! write(*,*) A2,eos%A2
        ! write(*,*) A3,eos%A3
        ! write(*,*) B1,eos%B1
        ! write(*,*) B2,eos%B2
        ! write(*,*) B3,eos%B3
        ! write(*,*) B4,eos%B4
        ! write(*,*) B5,eos%B5
        ! store data
        eos%E0=E0
        eos%V0=V0
        eos%B0=B0
        eos%B0p=B0p
        eos%eta0=eta0
        eos%C0=C0
        eos%C1=C1
        eos%C2=C2
        eos%C3=C3
    end block findv0
end subroutine

!> set the parameters
module subroutine set_2D_birch_murnaghan_parameters(eos,E0,V0,eta0,B0,B0p,C0,C1,C2,C3)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(out) :: eos
    real(flyt), intent(in) :: E0,V0,eta0,B0,B0p,C0,C1,C2,C3

    eos%E0=E0*lo_eV_to_Hartree               ! input in eV/atom, convert to Hartree/atom
    eos%V0=V0*lo_volume_A_to_bohr            ! input in A^3/atom, onvert to bohr^3/atom
    eos%B0=B0*lo_pressure_GPa_to_HartreeBohr ! input in GPa, converted to Hartree/bohr^3
    eos%B0p=B0p                              ! dimensionless
    eos%eta0=eta0                            ! dimensionless
    eos%C0=C0                                ! dimensionless
    eos%C1=C1                                ! dimensionless
    eos%C2=C2                                ! dimensionless
    eos%C3=C3                                ! dimensionless

    eos%A0=(-8*eos%E0*(-2 + C0*eta0**2) - 3*eos%B0*(-18 + 3*eos%B0p + 18*C1*eos%eta0 - 6*C3*eos%eta0 + 4*C2*eos%eta0**2)*eos%V0)/16
    eos%A1=(3*eos%B0*(9*eos%B0p + 4*(-12 + 7*C1*eos%eta0 - 3*C3*eos%eta0 + C2*eos%eta0**2))*eos%V0**1.6666666666666667_flyt)/16
    eos%A2=(-3*eos%B0*(9*eos%B0p + 10*C1*eos%eta0 - 6*(7 + C3*eos%eta0))*eos%V0**2.3333333333333335_flyt)/16
    eos%A3=(9*eos%B0*(-4 + eos%B0p)*eos%V0**3)/16
    eos%B1=C0*eos%E0*eos%eta0 + (3*eos%B0*(9*C1 - 3*C3 + 4*C2*eos%eta0)*eos%V0)/8
    eos%B2=(-3*eos%B0*(7*C1 - 3*C3 + 2*C2*eos%eta0)*eos%V0**1.6666666666666667_flyt)/4
    eos%B3=(-2*C0*eos%E0 - 3*eos%B0*C2*eos%V0)/4
    eos%B4=(3*eos%B0*C2*eos%V0**1.6666666666666667_flyt)/4
    eos%B5=(3*eos%B0*(5*C1 - 3*C3)*eos%V0**2.3333333333333335_flyt)/8
end subroutine

!> 2D Birch-Munaghan energy
module elemental function twod_birch_murnaghan_energy_from_volume_eta(eos,V,eta) result(E)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> structural parameter
    real(flyt), intent(in) :: eta
    !> Energy
    real(flyt) :: E

    real(flyt) :: x
    x=v**(-2.0_flyt/3.0_flyt)
    E = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3 + eos%B1*eta + eos%B2*x*eta + eos%B3*eta**2 + eos%B4*x*eta**2 + eos%B5*eta*x**2
end function
!> 2D Birch-Munaghan energy with eta contraction
module elemental function twod_birch_murnaghan_energy_from_volume(eos,V) result(E)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Energy
    real(flyt) :: E

    real(flyt) :: x,eta
    x=v**(-2.0_flyt/3.0_flyt)
    eta=-( eos%B1 + eos%B2*x + eos%B5*x**2 )/( 2*eos%B3 + 2*eos%B4*x )
    E = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3 + eos%B1*eta + eos%B2*x*eta + eos%B3*eta**2 + eos%B4*x*eta**2 + eos%B5*eta*x**2
end function

!> 2D Birch-Munaghan pressure
module elemental function twod_birch_murnaghan_pressure_from_volume_eta(eos,V,eta) result(P)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> structural parameter
    real(flyt), intent(in) :: eta
    !> Energy
    real(flyt) :: P

    real(flyt) :: x
    x=v**(-2.0_flyt/3.0_flyt)
    P = (-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*eta + 2*eos%B5*x*eta + eos%B4*eta**2))/3.0_flyt
    P = -P
end function
!> 2D Birch-Munaghan pressure with eta contractoin
module elemental function twod_birch_murnaghan_pressure_from_volume(eos,V) result(P)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Pressure
    real(flyt) :: P

    real(flyt) :: x,eta
    x=v**(-2.0_flyt/3.0_flyt)
    eta=-( eos%B1 + eos%B2*x + eos%B5*x**2 )/( 2*eos%B3 + 2*eos%B4*x )
    P = (-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*eta + 2*eos%B5*x*eta + eos%B4*eta**2))/3.0_flyt
    P = -P
end function
!!> 2D Birch-Munaghan bulk modulus from volume eta
!elemental function twod_birch_bulkmodulus_from_volume_eta(eos,V,eta) result(P)
!    !> Birch-Murnaghan equation of state
!    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
!    !> Volume
!    real(flyt), intent(in) :: V
!    !> structural parameter
!    real(flyt), intent(in) :: eta
!    !> Energy
!    real(flyt) :: P
!
!    real(flyt) :: x
!    x=v**(-2.0_flyt/3.0_flyt)
!    P = (-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*eta + 2*eos%B5*x*eta + eos%B4*eta**2))/3.0_flyt
!    P = -P
!end function
!!> 2D Birch-Munaghan bulk modulus from volume
!elemental function twod_birch_bulkmodulus_from_volume(eos,V) result(P)
!    !> Birch-Murnaghan equation of state
!    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
!    !> Volume
!    real(flyt), intent(in) :: V
!    !> structural parameter
!    !> Energy
!    real(flyt) :: P
!
!    real(flyt) :: x,eta
!    x=v**(-2.0_flyt/3.0_flyt)
!    P = (-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*eta + 2*eos%B5*x*eta + eos%B4*eta**2))/3.0_flyt
!    P = -P
!end function

!> eta from volume
module elemental function twod_birch_murnaghan_eta_from_volume(eos,V) result(eta)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Structural parameter
    real(flyt) :: eta

    real(flyt) :: x
    x=V**(-2.0_flyt/3.0_flyt)
    eta=-( eos%B1 + eos%B2*x + eos%B5*x**2 )/( 2*eos%B3 + 2*eos%B4*x )
end function

!> volume from pressure, with eta contracted
module elemental function twod_birch_murnaghan_volume_from_pressure(eos,P) result(V)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
    !> Pressure
    real(flyt), intent(in) :: P
    !> Volume
    real(flyt) :: V

    real(flyt) :: Pd,pw,pwref,x,eta
    integer :: iter

    V=eos%V0 ! starting volume
    pwref=P
    do iter=1,150
        x=V**(-2.0_flyt/3.0_flyt)
        eta=-( eos%B1 + eos%B2*x + eos%B5*x**2 )/( 2*eos%B3 + 2*eos%B4*x )
        pw = -(-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*eta + 2*eos%B5*x*eta + eos%B4*eta**2))/3.0_flyt
        pw=pw-pwref
        if ( abs(pw) .lt. 1E-15_flyt ) then
            ! we are done
            exit
        else
            ! iterate with the pressure derivative
            Pd = (2*(x**4)*(5*eos%A1 + 14*eos%A2*x + 27*eos%A3*x**2 + 5*eos%B2*eta + 14*eos%B5*x*eta + 5*eos%B4*eta**2))/9.0_flyt
            if ( abs(pw/Pd) .gt. V*0.25_flyt ) then
                ! Too far for a Newton step, do some stupid steps
                if ( pw .lt. 0.0_flyt ) then
                    V=V*0.95_flyt
                else
                    V=V*1.05_flyt
                endif
            else
                ! iterate with the pressure derivative
                V=V+pw/Pd
            endif
        endif
    enddo
end function

!> energy from pressure, with eta contracted
module elemental function twod_birch_murnaghan_energy_from_pressure(eos,P) result(E)
    !> Birch-Murnaghan equation of state
    class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
    !> Pressure
    real(flyt), intent(in) :: P
    !> Volume
    real(flyt) :: E

    real(flyt) :: Pd,pw,pwref,x,eta,V
    integer :: iter

    V=eos%V0 ! starting volume
    pwref=P
    do iter=1,150
        x=V**(-2.0_flyt/3.0_flyt)
        eta=-( eos%B1 + eos%B2*x + eos%B5*x**2 )/( 2*eos%B3 + 2*eos%B4*x )
        pw = -(-2*x**(2.5_flyt)*(eos%A1 + 2*eos%A2*x + 3*eos%A3*x**2 + eos%B2*eta + 2*eos%B5*x*eta + eos%B4*eta**2))/3.0_flyt
        pw=pw-pwref
        if ( abs(pw) .lt. 1E-15_flyt ) then
            ! we are done
            exit
        else
            Pd = (2*(x**4)*(5*eos%A1 + 14*eos%A2*x + 27*eos%A3*x**2 + 5*eos%B2*eta + 14*eos%B5*x*eta + 5*eos%B4*eta**2))/9.0_flyt
            if ( abs(pw/Pd) .gt. V*0.25_flyt ) then
                ! Too far for a Newton step, do some stupid steps
                if ( pw .lt. 0.0_flyt ) then
                    V=V*0.95_flyt
                else
                    V=V*1.05_flyt
                endif
            else
                ! iterate with the pressure derivative
                V=V+pw/Pd
            endif
        endif
    enddo
    E = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3 + eos%B1*eta + eos%B2*x*eta + eos%B3*eta**2 + eos%B4*x*eta**2 + eos%B5*eta*x**2
end function

end submodule
