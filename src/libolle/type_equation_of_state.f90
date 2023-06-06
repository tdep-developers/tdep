#include "precompilerdefinitions"
module type_equation_of_state
!! All manner of different equations of state
!@TODO add validity ranges
!@TODO add Tait equation of state to aid in inversion
use konstanter, only : flyt,lo_huge,lo_tol,lo_sqtol,lo_pressure_GPa_to_HartreeBohr,lo_pressure_HartreeBohr_to_GPa,&
                       lo_volume_A_to_bohr,lo_eV_to_Hartree
use gottochblandat, only: tochar,lo_linspace
use type_blas_lapack_wrappers, only: lo_dgels
implicit none

private
! generic equation of state
public :: lo_eos
public :: lo_eos_1d
public :: lo_eos_2d
! specific one-dimensional
public :: lo_eos_birch_murnaghan
public :: lo_eos_vinet
! specific two-dimensional
public :: lo_eos_2d_birch_murnaghan

!> generalized isothermal equation of state
type, abstract :: lo_eos
    !> energy at equilibrium volume
    real(flyt) :: E0=lo_huge
    !> equilibrium volume
    real(flyt) :: V0=lo_huge
    !> bulk modulus at equilibrium volume
    real(flyt) :: B0=lo_huge
    !> dB/dP at equilibrium volume
    real(flyt) :: B0p=lo_huge
    !> d^2B/dP^2 at equilibrium volume
    real(flyt) :: B0pp=lo_huge
    contains
        procedure(energy_from_volume),      deferred :: energy_from_volume
        procedure(pressure_from_volume),    deferred :: pressure_from_volume
        procedure(volume_from_pressure),    deferred :: volume_from_pressure
        procedure(energy_from_pressure),    deferred :: energy_from_pressure
end type

!> for 1-D fits
type, abstract, extends(lo_eos) :: lo_eos_1d
    contains
        procedure(bulkmodulus_from_volume), deferred :: bulkmodulus_from_volume
        procedure(generate),                deferred :: generate
end type

!> for 2-D fits
type, abstract, extends(lo_eos) :: lo_eos_2d
    !> equilibrium
    real(flyt) :: eta0=lo_huge
    !> Dimensionless eta-eta derivative
    real(flyt) :: C0=lo_huge
    !> Dimensionless eta-V derivative
    real(flyt) :: C1=lo_huge
    !> Dimensionless V-eta-eta derivative
    real(flyt) :: C2=lo_huge
    !> Dimensionless V-V-eta derivative
    real(flyt) :: C3=lo_huge
    contains
        procedure(generate_2d),              deferred :: generate_2d
        procedure(energy_from_volume_eta),   deferred :: energy_from_volume_eta
        procedure(pressure_from_volume_eta), deferred :: pressure_from_volume_eta
        procedure(eta_from_volume),          deferred :: eta_from_volume
end type

! A null type - the nothing equation of state. Make no sense right now, but makes life a lot easier later.
!type, extends(lo_eos) :: lo_eos_null
!end type

! abstract interfaces common to all equations of state
interface
    elemental function energy_from_volume(eos,V) result(E)
        import :: lo_eos,flyt
        class(lo_eos), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: E
    end function
    elemental function pressure_from_volume(eos,V) result(P)
        import :: lo_eos,flyt
        class(lo_eos), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: P
    end function
    elemental function volume_from_pressure(eos,P) result(V)
        import :: lo_eos,flyt
        class(lo_eos), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: V
    end function
    elemental function energy_from_pressure(eos,P) result(V)
        import :: lo_eos,flyt
        class(lo_eos), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: V
    end function
end interface

! abstract interfaces, 1D
interface
    subroutine generate(eos,E0,V0,B0,B0p,B0pp)
        import :: lo_eos_1d,flyt
        class(lo_eos_1d), intent(out) :: eos
        real(flyt), intent(in) :: E0,V0,B0,B0p
        real(flyt), intent(in), optional :: B0pp
    end subroutine
    elemental function bulkmodulus_from_volume(eos,V) result(B)
        import :: lo_eos_1d,flyt
        class(lo_eos_1d), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: B
    end function
end interface
! abstract interfaces, 2D
interface
    subroutine generate_2D(eos,E0,V0,eta0,B0,B0p,C0,C1,C2,C3)
        import :: lo_eos_2d,flyt
        class(lo_eos_2d), intent(out) :: eos
        real(flyt), intent(in) :: E0,V0,eta0,B0,B0p,C0,C1,C2,C3
    end subroutine
    elemental function energy_from_volume_eta(eos,V,eta) result(E)
        import :: lo_eos_2d,flyt
        class(lo_eos_2d), intent(in) :: eos
        real(flyt), intent(in) :: V,eta
        real(flyt) :: E
    end function
    elemental function pressure_from_volume_eta(eos,V,eta) result(E)
        import :: lo_eos_2d,flyt
        class(lo_eos_2d), intent(in) :: eos
        real(flyt), intent(in) :: V,eta
        real(flyt) :: E
    end function
    elemental function eta_from_volume(eos,V) result(eta)
        import :: lo_eos_2d,flyt
        class(lo_eos_2d), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: eta
    end function
end interface

! Now I can define the explicit equations of state

!> third order Birch-Murnaghan equation of state
type, extends(lo_eos_1d) :: lo_eos_birch_murnaghan
    !> precalculated factors
    real(flyt) :: A0=lo_huge,A1=lo_huge,A2=lo_huge,A3=lo_huge
    contains
        ! common to all
        procedure :: energy_from_volume=>birch_murnaghan_energy_from_volume
        procedure :: pressure_from_volume=>birch_murnaghan_pressure_from_volume
        procedure :: energy_from_pressure=>birch_murnaghan_energy_from_pressure
        procedure :: volume_from_pressure=>birch_murnaghan_volume_from_pressure
        ! for this one only
        procedure :: generate=>set_birch_murnaghan_parameters
        procedure :: fit=>fit_birch_murnaghan_parameters
        procedure :: bulkmodulus_from_volume=>birch_murnaghan_bulkmodulus_from_volume
end type
interface
    module subroutine fit_birch_murnaghan_parameters(eos,volume,energy,verbosity)
        class(lo_eos_birch_murnaghan), intent(out) :: eos
        real(flyt), dimension(:), intent(in) :: volume
        real(flyt), dimension(:), intent(in) :: energy
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine set_birch_murnaghan_parameters(eos,E0,V0,B0,B0p,B0pp)
        class(lo_eos_birch_murnaghan), intent(out) :: eos
        real(flyt), intent(in) :: E0,V0,B0,B0p
        real(flyt), intent(in), optional :: B0pp
    end subroutine
    module elemental function birch_murnaghan_energy_from_volume(eos,V) result(E)
        class(lo_eos_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: E
    end function
    module elemental function birch_murnaghan_pressure_from_volume(eos,V) result(P)
        class(lo_eos_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: P
    end function
    module elemental function birch_murnaghan_bulkmodulus_from_volume(eos,V) result(B)
        class(lo_eos_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: B
    end function
    module elemental function birch_murnaghan_volume_from_pressure(eos,P) result(V)
        class(lo_eos_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: V
    end function
    module elemental function birch_murnaghan_energy_from_pressure(eos,P) result(E)
        class(lo_eos_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: E
    end function
end interface

!> Vinet equation of state
type, extends(lo_eos_1d) :: lo_eos_vinet
    contains
        ! common to all EOS
        procedure :: energy_from_volume=>vinet_energy_from_volume
        procedure :: pressure_from_volume=>vinet_pressure_from_volume
        procedure :: energy_from_pressure=>vinet_energy_from_pressure
        procedure :: volume_from_pressure=>vinet_volume_from_pressure
        ! only for this
        procedure :: generate=>set_vinet_parameters
        procedure :: fit=>fit_vinet_parameters
        procedure :: bulkmodulus_from_volume=>vinet_bulkmodulus_from_volume
end type
interface
    module subroutine fit_vinet_parameters(eos,volume,energy,verbosity)
        !> Birch-Murnaghan equation of state
        class(lo_eos_vinet), intent(out) :: eos
        !> volumes
        real(flyt), dimension(:), intent(in) :: volume
        !> energy
        real(flyt), dimension(:), intent(in) :: energy
        !> how much to talk
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine set_vinet_parameters(eos,E0,V0,B0,B0p,B0pp)
        class(lo_eos_vinet), intent(out) :: eos
        real(flyt), intent(in) :: E0,V0,B0,B0p
        real(flyt), intent(in), optional :: B0pp
    end subroutine
    module elemental function vinet_energy_from_volume(eos,V) result(E)
        class(lo_eos_vinet), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: E
    end function
    module elemental function vinet_pressure_from_volume(eos,V) result(P)
        class(lo_eos_vinet), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: P
    end function
    module elemental function vinet_bulkmodulus_from_volume(eos,V) result(B)
        class(lo_eos_vinet), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: B
    end function
    module elemental function vinet_volume_from_pressure(eos,P) result(V)
        class(lo_eos_vinet), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: V
    end function
    module elemental function vinet_energy_from_pressure(eos,P) result(E)
        class(lo_eos_vinet), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: E
    end function
end interface

!> third order 2-D Birch-Murnaghan equation of state
type, extends(lo_eos_2d) :: lo_eos_2D_birch_murnaghan
    !> raw fit parameters
    real(flyt) :: A0=lo_huge,A1=lo_huge,A2=lo_huge,A3=lo_huge
    real(flyt) :: B1=lo_huge,B2=lo_huge,B3=lo_huge,B4=lo_huge,B5=lo_huge
    contains
        ! common to all EOS
        procedure :: energy_from_volume=>twod_birch_murnaghan_energy_from_volume
        procedure :: pressure_from_volume=>twod_birch_murnaghan_pressure_from_volume
        procedure :: energy_from_volume_eta=>twod_birch_murnaghan_energy_from_volume_eta
        procedure :: pressure_from_volume_eta=>twod_birch_murnaghan_pressure_from_volume_eta
        ! this one only
        procedure :: volume_from_pressure=>twod_birch_murnaghan_volume_from_pressure
        procedure :: energy_from_pressure=>twod_birch_murnaghan_energy_from_pressure
        procedure :: eta_from_volume=>twod_birch_murnaghan_eta_from_volume

        procedure :: generate_2D=>set_2D_birch_murnaghan_parameters
        procedure :: fit=>fit_2D_birch_murnaghan_parameters
end type
interface
    module subroutine fit_2D_birch_murnaghan_parameters(eos,volume,energy,eta,verbosity)
        class(lo_eos_2D_birch_murnaghan), intent(out) :: eos
        real(flyt), dimension(:), intent(in) :: volume
        real(flyt), dimension(:), intent(in) :: energy
        real(flyt), dimension(:), intent(in) :: eta
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine set_2D_birch_murnaghan_parameters(eos,E0,V0,eta0,B0,B0p,C0,C1,C2,C3)
        class(lo_eos_2D_birch_murnaghan), intent(out) :: eos
        real(flyt), intent(in) :: E0,V0,eta0,B0,B0p,C0,C1,C2,C3
    end subroutine
    module elemental function twod_birch_murnaghan_energy_from_volume_eta(eos,V,eta) result(E)
        class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt), intent(in) :: eta
        real(flyt) :: E
    end function
    module elemental function twod_birch_murnaghan_energy_from_volume(eos,V) result(E)
        class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: E
    end function
    module elemental function twod_birch_murnaghan_pressure_from_volume_eta(eos,V,eta) result(P)
        class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt), intent(in) :: eta
        real(flyt) :: P
    end function
    module elemental function twod_birch_murnaghan_pressure_from_volume(eos,V) result(P)
        class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: P
    end function
    module elemental function twod_birch_murnaghan_eta_from_volume(eos,V) result(eta)
        class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: V
        real(flyt) :: eta
    end function
    module elemental function twod_birch_murnaghan_volume_from_pressure(eos,P) result(V)
        class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: V
    end function
    module elemental function twod_birch_murnaghan_energy_from_pressure(eos,P) result(E)
        class(lo_eos_2D_birch_murnaghan), intent(in) :: eos
        real(flyt), intent(in) :: P
        real(flyt) :: E
    end function
end interface

!!> Tait equation of state, fourth order
!type, extends(lo_eos) :: lo_eos_tait
!end type

!contains

end module
