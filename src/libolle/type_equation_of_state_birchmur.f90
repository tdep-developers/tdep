#include "precompilerdefinitions"
submodule (type_equation_of_state) type_equation_of_state_birchmur
implicit none
contains

!> fit BM to data. Very robust.
module subroutine fit_birch_murnaghan_parameters(eos,volume,energy,verbosity)
    !> Birch-Murnaghan equation of state
    class(lo_eos_birch_murnaghan), intent(out) :: eos
    !> volumes
    real(flyt), dimension(:), intent(in) :: volume
    !> energy
    real(flyt), dimension(:), intent(in) :: energy
    !> how much to talk
    integer, intent(in) :: verbosity

    integer :: n

    n=size(energy)
    if ( verbosity .gt. 0 ) then
        write(*,*) ''
        write(*,*) 'FITTING BIRCH-MURNAGHAN EQUATION OF STATE USING ',tochar(n),' POINTS'
    endif

    ! do a least squares fit to get the 'A' parameters
    getpar: block
        real(flyt), dimension(:,:), allocatable :: wA,wB
        real(flyt) :: x
        integer :: i

        n=size(volume)
        lo_allocate(wA(n,4))
        lo_allocate(wB(n,1))
        ! create the coefficient matrix and energy vector
        do i=1,n
            x=volume(i)**(-2.0_flyt/3.0_flyt)
            wA(i,:)=[1.0_flyt,x,x**2,x**3]
            wB(i,1)=energy(i)
        enddo
        ! fit
        call lo_dgels(wA,wB)
        ! store values
        eos%A0=wB(1,1)
        eos%A1=wB(2,1)
        eos%A2=wB(3,1)
        eos%A3=wB(4,1)
        if ( verbosity .gt. 0 ) write(*,*) '... got raw fit: ',tochar([eos%A0,eos%A1,eos%A2,eos%A3])
    end block getpar

    ! find a reasonable V0
    findv0: block
        integer :: i,j,iter,ctr
        real(flyt) :: E0,E1,V0,P0,B0,B0p,Pd,Pdd,x

        ! Check that the minimum was not at the edges:
        E0=lo_huge
        V0=lo_huge
        j=0
        do i=1,n
            if ( energy(i) .lt. E0 ) then
                E0=min(E0,energy(i))
                V0=volume(i)
                j=i
            endif
        enddo

        if ( j .ne. 1 .and. j .ne. n ) then
            ! Newton-step to find V0
            E1=lo_huge
            ctr=0
            do iter=1,100
                x=v0**(-2.0_flyt/3.0_flyt)                                                 ! convert V to x
                E0  = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3                              ! energy
                P0  = -(-2*x**(2.5_flyt)*(eos%A1 + x*(2*eos%A2 + 3*eos%A3*x)))/3.0_flyt    ! pressure
                Pd  =  (2*x**4*(5*eos%A1 + x*(14*eos%A2 + 27*eos%A3*x)))/9.0_flyt          ! second derivative energy
                Pdd = (-8*x**5.5_flyt*(10*eos%A1 + x*(35*eos%A2 + 81*eos%A3*x)))/27.0_flyt ! third derivative energy
                !
                B0=Pd*V0
                B0p=-Pdd*V0*V0/(B0)-1
                if ( verbosity .gt. 0 ) write(*,"(1X,'   finding V0',1X,I3,4(1X,F18.7),1X,E19.12)") iter,E0,V0,B0,B0p,P0
                if ( abs(E1-E0) .lt. lo_sqtol ) then
                    ! we are done
                    ctr=ctr+1
                    if ( ctr .gt. 3 ) exit
                    V0=V0+P0/PD
                    E1=E0
                else
                    ! not done, iterate
                    V0=V0+P0/PD
                    E1=E0
                    ctr=0
                endif
            enddo
            ! Calculate the parameters back again, to make sure it makes sense
            !A0 = (16*E0 + 54*B0*V0 - 9*B0*B0p*V0)/16
            !A1 = (9*B0*(-16 + 3*B0p))/(16*(V0**(-5.0_flyt/3.0_flyt)))
            !A2 = (-9*B0*(-14 + 3*B0p))/(16*(V0**(-7.0_flyt/3.0_flyt)))
            !A3 = (9*(-4*B0*V0**3 + B0*B0p*V0**3))/16
            !write(*,*) 'A0',A0,eos%A0,abs(A0/eos%A0-1.0_flyt)
            !write(*,*) 'A1',A1,eos%A1,abs(A1/eos%A1-1.0_flyt)
            !write(*,*) 'A2',A2,eos%A2,abs(A2/eos%A2-1.0_flyt)
            !write(*,*) 'A3',A3,eos%A3,abs(A3/eos%A3-1.0_flyt)
        else
            x=v0**(-2.0_flyt/3.0_flyt)                                                 ! convert V to x
            E0  = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3                              ! energy
            P0  = -(-2*x**(2.5_flyt)*(eos%A1 + x*(2*eos%A2 + 3*eos%A3*x)))/3.0_flyt    ! pressure
            Pd  =  (2*x**4*(5*eos%A1 + x*(14*eos%A2 + 27*eos%A3*x)))/9.0_flyt          ! second derivative energy
            Pdd = (-8*x**5.5_flyt*(10*eos%A1 + x*(35*eos%A2 + 81*eos%A3*x)))/27.0_flyt ! third derivative energy
            !
            B0=Pd*V0
            B0p=-Pdd*V0*V0/(B0)-1
        !    ! Just make do with what I have, the fit is still good.
        !    if ( verbosity .gt. 0 ) write(*,*) 'WARNING: found the minimum energy at the edge of the data'
        !    !A0 = (16*E0 + 54*B0*V0 - 9*B0*B0p*V0)/16
        !    !A1 = (9*B0*(-16 + 3*B0p))/(16*(V0**(-5.0_flyt/3.0_flyt)))
        !    !A2 = (-9*B0*(-14 + 3*B0p))/(16*(V0**(-7.0_flyt/3.0_flyt)))
        !    !A3 = (9*(-4*B0*V0**3 + B0*B0p*V0**3))/16
        !    !write(*,*) 'A0',A0,eos%A0,abs(A0/eos%A0-1.0_flyt)
        !    !write(*,*) 'A1',A1,eos%A1,abs(A1/eos%A1-1.0_flyt)
        !    !write(*,*) 'A2',A2,eos%A2,abs(A2/eos%A2-1.0_flyt)
        !    !write(*,*) 'A3',A3,eos%A3,abs(A3/eos%A3-1.0_flyt)
        endif
        ! store data
        eos%E0=E0
        eos%V0=V0
        eos%B0=B0
        eos%B0p=B0p
    end block findv0
end subroutine
!> set the BM parameters
module subroutine set_birch_murnaghan_parameters(eos,E0,V0,B0,B0p,B0pp)
    !> Birch-Murnaghan equation of state
    class(lo_eos_birch_murnaghan), intent(out) :: eos
    !> Parameters
    real(flyt), intent(in) :: E0,V0,B0,B0p
    real(flyt), intent(in), optional :: B0pp

    eos%E0=E0*lo_eV_to_Hartree               ! input in eV/atom, convert to Hartree/atom
    eos%V0=V0*lo_volume_A_to_bohr            ! input in A^3/atom, onvert to bohr^3/atom
    eos%B0=B0*lo_pressure_GPa_to_HartreeBohr ! input in GPa, converted to Hartree/bohr^3
    eos%B0p=B0p                              ! dimensionless
    if ( present(B0pp) ) then
        ! Throw an error?
    endif
    eos%A0 = (16*eos%E0 + 54*eos%B0*eos%V0 - 9*eos%B0*eos%B0p*eos%V0)/16
    eos%A1 = (9*eos%B0*(-16 + 3*eos%B0p))/(16*(eos%V0**(-5.0_flyt/3.0_flyt)))
    eos%A2 = (-9*eos%B0*(-14 + 3*eos%B0p))/(16*(eos%V0**(-7.0_flyt/3.0_flyt)))
    eos%A3 = (9*(-4*eos%B0*eos%V0**3 + eos%B0*eos%B0p*eos%V0**3))/16
end subroutine
!> BM energy
module elemental function birch_murnaghan_energy_from_volume(eos,V) result(E)
    !> Birch-Murnaghan equation of state
    class(lo_eos_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Energy
    real(flyt) :: E

    real(flyt) :: x
    x=V**(-2.0_flyt/3.0_flyt)
    E = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3
end function
!> BM pressure
module elemental function birch_murnaghan_pressure_from_volume(eos,V) result(P)
    !> Birch-Murnaghan equation of state
    class(lo_eos_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Pressure
    real(flyt) :: P

    real(flyt) :: x
    x=V**(-2.0_flyt/3.0_flyt)
    P = -(-2*x**(2.5_flyt)*(eos%A1 + x*(2*eos%A2 + 3*eos%A3*x)))/3.0_flyt
end function
!> BM bulk modulues
module elemental function birch_murnaghan_bulkmodulus_from_volume(eos,V) result(B)
    !> Birch-Murnaghan equation of state
    class(lo_eos_birch_murnaghan), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Bulk modulus
    real(flyt) :: B

    real(flyt) :: x
    x=V**(-2.0_flyt/3.0_flyt)
    B=( 2*x**(2.5_flyt)*(5*eos%A1 + x*(14*eos%A2 + 27*eos%A3*x)) )/9.0_flyt
end function
!> BM vol from pressure
module elemental function birch_murnaghan_volume_from_pressure(eos,P) result(V)
    !> Birch-Murnaghan equation of state
    class(lo_eos_birch_murnaghan), intent(in) :: eos
    !> Pressure (GPa)
    real(flyt), intent(in) :: P
    !> Volume
    real(flyt) :: V

    real(flyt) :: Pd,pw,x,pwref
    integer :: iter

    V=eos%V0
    pwref=P
    ! make quick rough guesses
    do iter=1,20
        x=V**(-2.0_flyt/3.0_flyt)  ! convert to x
        pw = -(-2*x**(2.5_flyt)*(eos%A1 + x*(2*eos%A2 + 3*eos%A3*x)))/3.0_flyt ! current pressure
        pw = pw-pwref
        if ( abs(pw) .lt. lo_tol ) then
            exit
        elseif ( pw .gt. 0.0_flyt ) then
            V=V*1.03_flyt
        else
            V=V*0.95_flyt
        endif
    enddo

    ! Newton-step to get it right
    do iter=1,1000
        x=V**(-2.0_flyt/3.0_flyt)  ! convert to x
        pw = -(-2*x**(2.5_flyt)*(eos%A1 + x*(2*eos%A2 + 3*eos%A3*x)))/3.0_flyt ! current pressure
        pw = pw-pwref
        if ( abs(pw) .lt. 1E-20_flyt ) then
            ! we are done
            exit
        else
            ! iterate with the pressure derivative
            Pd = (2*x**4*(5*eos%A1 + x*(14*eos%A2 + 27*eos%A3*x)))/9.0_flyt
            V=V+pw/Pd
        endif
    enddo
end function
!> BM energy from pressure
module elemental function birch_murnaghan_energy_from_pressure(eos,P) result(E)
    !> Birch-Murnaghan equation of state
    class(lo_eos_birch_murnaghan), intent(in) :: eos
    !> Pressure (GPa)
    real(flyt), intent(in) :: P
    !> Energy
    real(flyt) :: E

    real(flyt) :: Pd,pw,x,pwref,V
    integer :: iter

    ! copy of the above, to locate the correct volume
    V=eos%V0
    pwref=P
    ! stupid rough guesses first
    do iter=1,20
        x=V**(-2.0_flyt/3.0_flyt)  ! convert to x
        pw = -(-2*x**(2.5_flyt)*(eos%A1 + x*(2*eos%A2 + 3*eos%A3*x)))/3.0_flyt ! current pressure
        pw = pw-pwref
        if ( pw .gt. lo_tol ) then
            V=V*1.03_flyt
        else
            V=V*0.95_flyt
        endif
    enddo
    do iter=1,1000
        x=V**(-2.0_flyt/3.0_flyt)  ! convert to x
        pw = -(-2*x**(2.5_flyt)*(eos%A1 + x*(2*eos%A2 + 3*eos%A3*x)))/3.0_flyt ! current pressure
        pw = pw-pwref
        if ( abs(pw) .lt. 1E-20_flyt ) then
            ! we are done
            exit
        else
            ! iterate with the pressure derivative
            Pd = (2*x**4*(5*eos%A1 + x*(14*eos%A2 + 27*eos%A3*x)))/9.0_flyt
            V=V+pw/Pd
        endif
    enddo
    ! evaluate energy at this volume
    E = eos%A0+eos%A1*x+eos%A2*x**2+eos%A3*x**3
end function

end submodule
