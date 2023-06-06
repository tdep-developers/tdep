#include "precompilerdefinitions"
submodule (type_equation_of_state) type_equation_of_state_vinet
implicit none
contains

!> Fit Vinet to volume-energy data. A little robust.
module subroutine fit_vinet_parameters(eos,volume,energy,verbosity)
    !> Birch-Murnaghan equation of state
    class(lo_eos_vinet), intent(out) :: eos
    !> volumes
    real(flyt), dimension(:), intent(in) :: volume
    !> energy
    real(flyt), dimension(:), intent(in) :: energy
    !> how much to talk
    integer, intent(in) :: verbosity

    type(lo_eos_birch_murnaghan) :: bm
    integer :: n

    n=size(energy)
    if ( verbosity .gt. 0 ) then
        write(*,*) ''
        write(*,*) 'FITTING VINET EQUATION OF STATE USING ',tochar(n),' POINTS'
    endif

    ! get a rough guess for the nonlinear solver
    rough: block
        real(flyt) :: V1,V2,z1,z2,e1,e2,deltaV,f1,f2,y0,y1,odiff
        real(flyt) :: aE0,aV0,aB0,aB0p,adiff
        real(flyt) :: bE0,bV0,bB0,bB0p,bdiff
        complex(flyt) :: cb01,cb02,x01,x02

        ! Get a Birch-Murnaghan to get something to start from.
        call bm%fit(volume,energy,verbosity=0)
        ! By calculating Bmod/P I can get equations for B0p and V0, if I do it for two volumes.
        ! The equations have two solutions, I pick the best one out of the two. I have to make sure
        ! the two volumes have nonzero pressure.
        V1=minval(volume)
        V2=maxval(volume)
        deltaV=V2-V1
        V1=V1+0.333*deltaV
        V2=V2-0.3333*deltaV

        do
            if ( abs(bm%pressure_from_volume(V1)) .lt. 1E-3_flyt*lo_pressure_GPa_to_HartreeBohr ) then
                V1=V1+1E-5_flyt
            else
                exit
            endif
        enddo
        do
            if ( abs(bm%pressure_from_volume(V2)) .lt. 1E-3_flyt*lo_pressure_GPa_to_HartreeBohr ) then
                V2=V2-1E-5_flyt
            else
                exit
            endif
        enddo

        ! get Bmod/P
        e1=-bm%bulkmodulus_from_volume(V1)/bm%pressure_from_volume(V1)
        e2=-bm%bulkmodulus_from_volume(V2)/bm%pressure_from_volume(V2)
        ! And the strange intermediate parameters
        z1=V1**(1.0_flyt/3.0_flyt)
        z2=V2**(1.0_flyt/3.0_flyt)
        cb01=(-((2 + 3*e2)**2*z1**4) + (2 + 3*e2)*(5 + 3*e1 + 3*e2)*z1**3*z2 + (2 + 3*e1)*(5 + 3*e1 + 3*e2)*z1*z2**3 - &
                (2 + 3*e1)*z2**2*((2 + 3*e1)*z2**2 + Sqrt((z1 - z2)*(2*z1 + 3*e2*z1 - 2*z2 - 3*e1*z2)*&
                ((2 + 3*e2)*z1**2 - 3*(e1 + e2)*z1*z2 + (2 + 3*e1)*z2**2))) + &
                z1**2*(-3*(4 + 5*e1 + 5*e2 + 6*e1*e2)*z2**2 + (2 + 3*e2)*&
                Sqrt((z1 - z2)*(2*z1 + 3*e2*z1 - 2*z2 - 3*e1*z2)*((2 + 3*e2)*z1**2 - 3*(e1 + e2)*z1*z2 + (2 + 3*e1)*z2**2))))/&
                (3*z1*(z1 - z2)*z2*((2 + 3*e2)*z1 - (2 + 3*e1)*z2))

        cb02=(-((2 + 3*e2)**2*z1**4) + (2 + 3*e2)*(5 + 3*e1 + 3*e2)*z1**3*z2 +(2+3*e1)*(5+3*e1+3*e2)*z1*z2**3 +(2 + 3*e1)*z2**2*&
                (-((2 + 3*e1)*z2**2) + Sqrt((z1 - z2)*(2*z1 + 3*e2*z1 - 2*z2 - 3*e1*z2)*((2 + 3*e2)*z1**2 - 3*(e1 + e2)*z1*z2 + &
                (2 + 3*e1)*z2**2))) + z1**2*(-3*(4 + 5*e1 + 5*e2 + 6*e1*e2)*z2**2 - (2 + 3*e2)* &
                Sqrt((z1 - z2)*(2*z1 + 3*e2*z1 - 2*z2 - 3*e1*z2)*((2 + 3*e2)*z1**2 - 3*(e1 + e2)*z1*z2 + (2 + 3*e1)*z2**2))))/&
                (3.*z1*(z1 - z2)*z2*((2 + 3*e2)*z1 - (2 + 3*e1)*z2))

        x01=(z1 + z2 + Sqrt((z1 - z2)*((2 + 3*e2)*z1 - (2 + 3*e1)*z2)*((2 + 3*e2)*z1**2 - 3*(e1 + e2)*z1*z2 + (2 + 3*e1)*z2**2))/&
            (-((2 + 3*e2)*z1) + (2 + 3*e1)*z2))/2.0_flyt

        x02=(z1 + z2 + Sqrt((z1 - z2)*(2*z1 + 3*e2*z1 - 2*z2 - 3*e1*z2)*&
                ((2 + 3*e2)*z1**2 - 3*(e1 + e2)*z1*z2 + (2 + 3*e1)*z2**2))/(2*z1 + 3*e2*z1 - 2*z2 - 3*e1*z2))/2.0_flyt

        e1=bm%energy_from_volume(V1)
        e2=bm%energy_from_volume(V2)
        ! Convert these to reasonable numbers
        aB0p=real(cb01)
        aV0=real(x01)**3
            y0=(V1/aV0)**(1.0_flyt/3.0_flyt)
            y1=exp( -(1.5_flyt)*(aB0p-1)*(y0-1) )
            f1=(aV0*(-2 + (5 + 3*aB0p*(y0 - 1) - 3*y0)*y1))/(aB0p - 1)**2
            y0=(V2/aV0)**(1.0_flyt/3.0_flyt)
            y1=exp( -(1.5_flyt)*(aB0p-1)*(y0-1) )
            f2=(aV0*(-2 + (5 + 3*aB0p*(y0 - 1) - 3*y0)*y1))/(aB0p - 1)**2
        aE0=-(-e2*f1 + e1*f2)/(f1 - f2)
        aB0=-(e1 - e2)/(2*(f1 - f2))

        bB0p=real(cb02)
        bV0=real(x02)**3
            y0=(V1/bV0)**(1.0_flyt/3.0_flyt)
            y1=exp( -(1.5_flyt)*(bB0p-1)*(y0-1) )
            f1=(bV0*(-2 + (5 + 3*bB0p*(y0 - 1) - 3*y0)*y1))/(bB0p - 1)**2
            y0=(V2/bV0)**(1.0_flyt/3.0_flyt)
            y1=exp( -(1.5_flyt)*(bB0p-1)*(y0-1) )
            f2=(bV0*(-2 + (5 + 3*bB0p*(y0 - 1) - 3*y0)*y1))/(bB0p - 1)**2
        bE0=-((-(e2*f1) + e1*f2)/(f1 - f2))
        bB0=-(e1 - e2)/(2*(f1 - f2))

        eos%E0 =bm%E0
        eos%V0 =bm%V0
        eos%B0 =bm%B0
        eos%B0p=bm%B0p
        odiff=sum( (eos%energy_from_volume(volume)-energy)**2 )/size(volume)
        eos%E0 =aE0
        eos%V0 =aV0
        eos%B0 =aB0
        eos%B0p=aB0p
        if ( aB0 .gt. 0.0_flyt .and. aB0p .gt. 0.0_flyt .and. aV0 .gt. 0.0_flyt ) then
            adiff=sum( (eos%energy_from_volume(volume)-energy)**2 )/size(volume)
        else
            adiff=1E10_flyt
        endif
        eos%E0 =bE0
        eos%V0 =bV0
        eos%B0 =bB0
        eos%B0p=bB0p
        bdiff=sum( (eos%energy_from_volume(volume)-energy)**2 )/size(volume)
        if ( bB0 .gt. 0.0_flyt .and. bB0p .gt. 0.0_flyt .and. bV0 .gt. 0.0_flyt ) then
            bdiff=sum( (eos%energy_from_volume(volume)-energy)**2 )/size(volume)
        else
            bdiff=1E11_flyt
        endif

        ! Decide on my starting guess based on energy differences
        if ( odiff .lt. adiff .and. odiff .lt. bdiff ) then
            eos%E0 =bm%E0
            eos%V0 =bm%V0
            eos%B0 =bm%B0
            eos%B0p=bm%B0p
        elseif ( adiff .lt. odiff .and. adiff .lt. bdiff ) then
            eos%E0 =aE0
            eos%V0 =aV0
            eos%B0 =aB0
            eos%B0p=aB0p
        else
            eos%E0 =bE0
            eos%V0 =bV0
            eos%B0 =bB0
            eos%B0p=bB0p
        endif
    end block rough

    findpar: block
        integer :: outer,inner,ctr,pct
        real(flyt), dimension(4) :: grad,hess,oldgrad
        real(flyt) :: step,diff0,diff1,diff2,x

        ! start with a little steepest descent
        call gradient(eos,energy,volume,oldgrad,hess)
        diff2=lo_huge
        ctr=0
        pct=1
        do outer=1,20000
            call gradient(eos,energy,volume,grad,hess)
            x=norm2(grad)
            grad=(grad*0.5_flyt+0.5_flyt*oldgrad)
            oldgrad=grad
            grad=grad/hess
            !
            step=1E-5_flyt/norm2(grad)
            diff0=sum( ( energy-eos%energy_from_volume(volume) )**2 )/n
            inloop: do inner=1,10000
                eos%E0 =eos%E0  -grad(1)*step
                eos%V0 =eos%V0  -grad(2)*step
                eos%B0 =eos%B0  -grad(3)*step
                eos%B0p=eos%B0p -grad(4)*step
                diff1=sum( ( energy-eos%energy_from_volume(volume) )**2 )/n
                if ( diff1 .lt. diff0 ) then
                    step=step*1.1_flyt
                    diff0=diff1
                else
                    eos%E0=eos%E0  +grad(1)*step
                    eos%V0=eos%V0  +grad(2)*step
                    eos%B0=eos%B0  +grad(3)*step
                    eos%B0p=eos%B0p+grad(4)*step
                    step=step*0.5_flyt
                    if ( step .lt. 1E-20_flyt ) then
                        exit inloop
                    endif
                endif
            enddo inloop
            if ( verbosity .gt. 0 ) then
                if ( mod(outer,pct) .eq. 0 ) then
                    write(*,*) eos%E0,sqrt(diff0),outer,inner,x,pct
                    pct=pct+outer
                endif
            endif

            if ( abs(diff2-diff0) .lt. 1E-20_flyt ) then
                ctr=ctr+1
                if ( ctr .ge. 16 ) exit
            else
                diff2=diff0
                ctr=0
            endif
        enddo
    end block findpar

    contains
    !> gradient of the energy difference squared
    subroutine gradient(eos,energy,volume,grad,hess)
        type(lo_eos_vinet), intent(inout) :: eos
        real(flyt), dimension(:), intent(in) :: energy,volume
        real(flyt), dimension(4), intent(out) :: grad,hess

        real(flyt) :: E0,V0,B0,B0p,invsz
        real(flyt), dimension(9) :: sc_wt1,sc_wt2,sc_dl,sc_val
        real(flyt), parameter :: delta=5E-5_flyt
        integer :: j
        ! Stencil weights. 9-point derivatives.
        sc_dl=[-4,-3,-2,-1,0,1,2,3,4]*delta
        sc_wt1=[3,-32,168,-672,0,672,-168,32,-3]/(840.0_flyt*delta)
        sc_wt2=[-1.0_flyt/560.0_flyt, 8.0_flyt/315.0_flyt, -1.0_flyt/5.0_flyt, 8.0_flyt/5.0_flyt,-205.0_flyt/72.0_flyt,&
                8.0_flyt/5.0_flyt,-1.0_flyt/5.0_flyt,8.0_flyt/315.0_flyt,-1.0_flyt/560.0_flyt]/(delta**2)
        grad=0.0_flyt
        hess=0.0_flyt
        ! Reset the parameters
        E0=eos%E0
        V0=eos%V0
        B0=eos%B0
        B0p=eos%B0p
        invsz=1.0_flyt/size(energy)
        ! Start in E0-direction.
        do j=1,9
            ! reset parameters
            eos%E0 =E0+sc_dl(j)
            eos%V0 =V0
            eos%B0 =B0
            eos%B0p=B0p
            sc_val(j)=sum( ( energy-eos%energy_from_volume(volume) )**2 )*invsz
        enddo
        grad(1)=sum(sc_wt1*sc_val)
        hess(1)=sum(sc_wt2*sc_val)
        do j=1,9
            ! reset parameters
            eos%E0 =E0
            eos%V0 =V0+sc_dl(j)
            eos%B0 =B0
            eos%B0p=B0p
            sc_val(j)=sum( ( energy-eos%energy_from_volume(volume) )**2 )*invsz
        enddo
        grad(2)=sum(sc_wt1*sc_val)
        hess(2)=sum(sc_wt2*sc_val)
        do j=1,9
            ! reset parameters
            eos%E0 =E0
            eos%V0 =V0
            eos%B0 =B0+sc_dl(j)
            eos%B0p=B0p
            sc_val(j)=sum( ( energy-eos%energy_from_volume(volume) )**2 )*invsz
        enddo
        grad(3)=sum(sc_wt1*sc_val)
        hess(3)=sum(sc_wt2*sc_val)
        do j=1,9
            ! reset parameters
            eos%E0 =E0
            eos%V0 =V0
            eos%B0 =B0
            eos%B0p=B0p+sc_dl(j)
            sc_val(j)=sum( ( energy-eos%energy_from_volume(volume) )**2 )*invsz
        enddo
        grad(4)=sum(sc_wt1*sc_val)
        hess(4)=sum(sc_wt2*sc_val)
        ! and set the parameters back to the original
        eos%E0 =E0
        eos%V0 =V0
        eos%B0 =B0
        eos%B0p=B0p
    end subroutine
end subroutine

!> set the Vinet parameters
module subroutine set_vinet_parameters(eos,E0,V0,B0,B0p,B0pp)
    !> Vinet equation of state
    class(lo_eos_vinet), intent(out) :: eos
    !> Parameters
    real(flyt), intent(in) :: E0,V0,B0,B0p
    real(flyt), intent(in), optional :: B0pp

    eos%E0=E0*lo_eV_to_Hartree               ! input in eV/atom, convert to Hartree/atom
    eos%V0=V0*lo_volume_A_to_bohr            ! input in A^3/atom, onvert to bohr^3/atom
    eos%B0=B0*lo_pressure_GPa_to_HartreeBohr ! input in GPa, converted to eV/A^3
    eos%B0p=B0p                              ! dimensionless
    if ( present(B0pp) ) then
        ! Do something about this.
    endif
end subroutine

!> Vinet energy
module elemental function vinet_energy_from_volume(eos,V) result(E)
    !> Vinet equation of state
    class(lo_eos_vinet), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Energy
    real(flyt) :: E

    real(flyt) :: y0,y1
    E=0.0_flyt
    y0=(V/eos%V0)**(1.0_flyt/3.0_flyt)
    y1=exp( -(1.5_flyt)*(eos%B0p-1)*(y0-1) )
    E=eos%E0-(2*eos%B0*eos%V0*(-2 + (5 + 3*eos%B0p*(y0-1) - 3*y0)*y1)) / (eos%B0p-1)**2
end function
!> Vinet pressure
module elemental function vinet_pressure_from_volume(eos,V) result(P)
    !> Vinet equation of state
    class(lo_eos_vinet), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Pressure (GPa)
    real(flyt) :: P

    real(flyt) :: y0,y1
    y0=(V/eos%V0)**(1.0_flyt/3.0_flyt)
    y1=exp( (1.5_flyt)*(eos%B0p-1)*(1-y0) )
    P=-(3.0_flyt*eos%B0*(y0-1)*y1)/y0**2
    P=P
end function
!> Vinet bulk modulus
module elemental function vinet_bulkmodulus_from_volume(eos,V) result(B)
    !> Vinet equation of state
    class(lo_eos_vinet), intent(in) :: eos
    !> Volume
    real(flyt), intent(in) :: V
    !> Bulk modulus (GPa)
    real(flyt) :: B

    real(flyt) :: y0,y1

    y0=(V/eos%V0)**(1.0_flyt/3.0_flyt)
    y1=exp( (1.5_flyt)*(eos%B0p-1)*(1-y0) )
    B=-(eos%B0*(-4 + (5 - 3*eos%B0p)*y0 + 3*(-1 + eos%B0p)*y0**2)*y1)/(2*y0**2)
    B=B
end function
!> Vinet vol from pressure
module elemental function vinet_volume_from_pressure(eos,P) result(V)
    !> Vinet equation of state
    class(lo_eos_vinet), intent(in) :: eos
    !> Pressure (GPa)
    real(flyt), intent(in) :: P
    !> Volume
    real(flyt) :: V

    real(flyt) :: Pd,pw,pwref,y0,y1
    integer :: iter

    V=eos%V0
    pwref=P
    ! stupid rough guesses first
    do iter=1,30
        y0=(V/eos%V0)**(1.0_flyt/3.0_flyt)
        y1=exp( (1.5_flyt)*(eos%B0p-1)*(1-y0) )
        pw=-(3.0_flyt*eos%B0*(y0-1)*y1)/y0**2
        pw = pw-pwref
        if ( pw .gt. lo_tol ) then
            V=V*1.03_flyt
        else
            V=V*0.95_flyt
        endif
    enddo

    do iter=1,1000
        y0=(V/eos%V0)**(1.0_flyt/3.0_flyt)
        y1=exp( (1.5_flyt)*(eos%B0p-1)*(1-y0) )
        pw=-(3.0_flyt*eos%B0*(y0-1)*y1)/y0**2
        pw = pw-pwref
        if ( abs(pw) .lt. 1E-20_flyt ) then
            ! we are done
            exit
        else
            ! iterate with the pressure derivative
            Pd =-(eos%B0*(-4 + (5 - 3*eos%B0p)*y0 + 3*(-1 + eos%B0p)*y0**2)*y1)/(2*y0**2)/eos%V0
            V=V+pw/Pd
        endif
    enddo
end function
!> BM energy from pressure
module elemental function vinet_energy_from_pressure(eos,P) result(E)
    !> Vinet equation of state
    class(lo_eos_vinet), intent(in) :: eos
    !> Pressure (GPa)
    real(flyt), intent(in) :: P
    !> Energy
    real(flyt) :: E

    real(flyt) :: Pd,pw,y0,y1,pwref,V
    integer :: iter

    E=0.0_flyt
    V=eos%V0
    pwref=P
    ! stupid rough guesses first
    do iter=1,20
        y0=(V/eos%V0)**(1.0_flyt/3.0_flyt)
        y1=exp( (1.5_flyt)*(eos%B0p-1)*(1-y0) )
        pw=-(3.0_flyt*eos%B0*(y0-1)*y1)/y0**2
        pw = pw-pwref
        if ( pw .gt. lo_tol ) then
            V=V*1.03_flyt
        else
            V=V*0.95_flyt
        endif
    enddo

    do iter=1,1000
        y0=(V/eos%V0)**(1.0_flyt/3.0_flyt)
        y1=exp( (1.5_flyt)*(eos%B0p-1)*(1-y0) )
        pw=-(3.0_flyt*eos%B0*(y0-1)*y1)/y0**2
        pw = pw-pwref
        if ( abs(pw) .lt. 1E-20_flyt ) then
            ! we are done
            exit
        else
            ! iterate with the pressure derivative
            Pd =-(eos%B0*(-4 + (5 - 3*eos%B0p)*y0 + 3*(-1 + eos%B0p)*y0**2)*y1)/(2*y0**2)/eos%V0
            V=V+pw/Pd
        endif
    enddo
    E=eos%E0-(2*eos%B0*eos%V0*(-2 + (5 + 3*eos%B0p*(y0-1) - 3*y0)*y1)) / (eos%B0p-1)**2
end function

end submodule
