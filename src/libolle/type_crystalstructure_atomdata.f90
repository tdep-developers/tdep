#include "precompilerdefinitions"
submodule (type_crystalstructure) type_crystalstructure_atomdata
implicit none
contains

!> Returns the atomic symbol from atomic number
module function z_to_symbol(z_nucleus) result(symbol)
    !> Nuclear charge, e.g. 14
    integer, intent(in) :: z_nucleus
    !> Atomic symbol, e.g. "Si"
    character(len=2) :: symbol

    character(len=2), parameter :: symbols_of_z(103) = [&
        &'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
        &'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
        &'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
        &'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
        &'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
        &'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
        &'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
        &'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
        &'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
        &'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
        &'Md', 'No', 'Lr']
    symbol = symbols_of_z(z_nucleus)
end function

!> Translates an atomic symbol into atomic number
module function symbol_to_z(symbol) result(z_nucleus)
    !> Atomic symbol, e.g. "Si"
    character(len=*), intent(in) :: symbol
    !> Nuclear charge, e.g. 14
    integer :: z_nucleus
    !
    select case(trim(adjustl(symbol)))
        case('H' )  ;  z_nucleus=1
        case('He')  ;  z_nucleus=2
        case('Li')  ;  z_nucleus=3
        case('Be')  ;  z_nucleus=4
        case('B' )  ;  z_nucleus=5
        case('C' )  ;  z_nucleus=6
        case('N' )  ;  z_nucleus=7
        case('O' )  ;  z_nucleus=8
        case('F' )  ;  z_nucleus=9

        case('Ne')  ;  z_nucleus=10
        case('Na')  ;  z_nucleus=11
        case('Mg')  ;  z_nucleus=12
        case('Al')  ;  z_nucleus=13
        case('Si')  ;  z_nucleus=14
        case('P' )  ;  z_nucleus=15
        case('S' )  ;  z_nucleus=16
        case('Cl')  ;  z_nucleus=17
        case('Ar')  ;  z_nucleus=18
        case('K' )  ;  z_nucleus=19

        case('Ca')  ;  z_nucleus=20
        case('Sc')  ;  z_nucleus=21
        case('Ti')  ;  z_nucleus=22
        case('V' )  ;  z_nucleus=23
        case('Cr')  ;  z_nucleus=24
        case('Mn')  ;  z_nucleus=25
        case('Fe')  ;  z_nucleus=26
        case('Co')  ;  z_nucleus=27
        case('Ni')  ;  z_nucleus=28
        case('Cu')  ;  z_nucleus=29

        case('Zn')  ;  z_nucleus=30
        case('Ga')  ;  z_nucleus=31
        case('Ge')  ;  z_nucleus=32
        case('As')  ;  z_nucleus=33
        case('Se')  ;  z_nucleus=34
        case('Br')  ;  z_nucleus=35
        case('Kr')  ;  z_nucleus=36
        case('Rb')  ;  z_nucleus=37
        case('Sr')  ;  z_nucleus=38
        case('Y' )  ;  z_nucleus=39

        case('Zr')  ;  z_nucleus=40
        case('Nb')  ;  z_nucleus=41
        case('Mo')  ;  z_nucleus=42
        case('Tc')  ;  z_nucleus=43
        case('Ru')  ;  z_nucleus=44
        case('Rh')  ;  z_nucleus=45
        case('Pd')  ;  z_nucleus=46
        case('Ag')  ;  z_nucleus=47
        case('Cd')  ;  z_nucleus=48
        case('In')  ;  z_nucleus=49

        case('Sn')  ;  z_nucleus=50
        case('Sb')  ;  z_nucleus=51
        case('Te')  ;  z_nucleus=52
        case('I' )  ;  z_nucleus=53
        case('Xe')  ;  z_nucleus=54
        case('Cs')  ;  z_nucleus=55
        case('Ba')  ;  z_nucleus=56
        case('La')  ;  z_nucleus=57
        case('Ce')  ;  z_nucleus=58
        case('Pr')  ;  z_nucleus=59

        case('Nd')  ;  z_nucleus=60
        case('Pm')  ;  z_nucleus=61
        case('Sm')  ;  z_nucleus=62
        case('Eu')  ;  z_nucleus=63
        case('Gd')  ;  z_nucleus=64
        case('Tb')  ;  z_nucleus=65
        case('Dy')  ;  z_nucleus=66
        case('Ho')  ;  z_nucleus=67
        case('Er')  ;  z_nucleus=68
        case('Tm')  ;  z_nucleus=69

        case('Yb')  ;  z_nucleus=70
        case('Lu')  ;  z_nucleus=71
        case('Hf')  ;  z_nucleus=72
        case('Ta')  ;  z_nucleus=73
        case('W' )  ;  z_nucleus=74
        case('Re')  ;  z_nucleus=75
        case('Os')  ;  z_nucleus=76
        case('Ir')  ;  z_nucleus=77
        case('Pt')  ;  z_nucleus=78
        case('Au')  ;  z_nucleus=79

        case('Hg')  ;  z_nucleus=80
        case('Tl')  ;  z_nucleus=81
        case('Pb')  ;  z_nucleus=82
        case('Bi')  ;  z_nucleus=83
        case('Po')  ;  z_nucleus=84
        case('At')  ;  z_nucleus=85
        case('Rn')  ;  z_nucleus=86
        case('Fr')  ;  z_nucleus=87
        case('Ra')  ;  z_nucleus=88
        case('Ac')  ;  z_nucleus=89

        case('Th')  ;  z_nucleus=90
        case('Pa')  ;  z_nucleus=91
        case('U' )  ;  z_nucleus=92
        case('Np')  ;  z_nucleus=93
        case('Pu')  ;  z_nucleus=94
        case('Am')  ;  z_nucleus=95
        case('Cm')  ;  z_nucleus=96
        case('Bk')  ;  z_nucleus=97
        case('Cf')  ;  z_nucleus=98
        case('Es')  ;  z_nucleus=89

        case('Fm')  ;  z_nucleus=100
        case('Md')  ;  z_nucleus=101
        case('No')  ;  z_nucleus=102
        case('Lr')  ;  z_nucleus=103
        case default
            write(*,*) 'TROUBLE!!!, no atomic number available for '//trim(adjustl(symbol))
            write(*,*) 'Are you really sure that is an element?'
            stop
    end select
    !
end function symbol_to_z

! !> inelastic neutron cross-section
module function neutron_cross_section(atomic_number) result(xs)
    !> Nuclear charge, e.g. 14
    integer, intent(in) :: atomic_number
    !> inelastic cross section
    real(flyt) :: xs

    select case(trim(z_to_symbol(atomic_number)))
        case("H")
            xs=1.756800_flyt
        case("He")
            xs=1.340000_flyt
        case("Li")
            xs=0.454000_flyt
        case("Be")
            xs=7.630000_flyt
        case("B")
            xs=3.540000_flyt
        case("C")
            xs=5.551000_flyt
        case("N")
            xs=11.010000_flyt
        case("O")
            xs=4.232000_flyt
        case("F")
            xs=4.017000_flyt
        case("Ne")
            xs=2.620000_flyt
        case("Na")
            xs=1.660000_flyt
        case("Mg")
            xs=3.631000_flyt
        case("Al")
            xs=1.495000_flyt
        case("Si")
            xs=2.163000_flyt
        case("P")
            xs=3.307000_flyt
        case("S")
            xs=1.018600_flyt
        case("Cl")
            xs=11.525700_flyt
        case("Ar")
            xs=0.458000_flyt
        case("K")
            xs=1.690000_flyt
        case("Ca")
            xs=2.780000_flyt
        case("Sc")
            xs=19.000000_flyt
        case("Ti")
            xs=1.485000_flyt
        case("V")
            xs=0.018400_flyt
        case("Cr")
            xs=1.660000_flyt
        case("Mn")
            xs=1.750000_flyt
        case("Fe")
            xs=11.220000_flyt
        case("Co")
            xs=0.779000_flyt
        case("Ni")
            xs=13.300000_flyt
        case("Cu")
            xs=7.485000_flyt
        case("Zn")
            xs=4.054000_flyt
        case("Ga")
            xs=6.675000_flyt
        case("Ge")
            xs=8.420000_flyt
        case("As")
            xs=5.440000_flyt
        case("Se")
            xs=7.980000_flyt
        case("Br")
            xs=5.800000_flyt
        case("Kr")
            xs=7.670000_flyt
        case("Rb")
            xs=6.320000_flyt
        case("Sr")
            xs=6.190000_flyt
        case("Y")
            xs=7.550000_flyt
        case("Zr")
            xs=6.440000_flyt
        case("Nb")
            xs=6.253000_flyt
        case("Mo")
            xs=5.670000_flyt
        case("Tc")
            xs=5.800000_flyt
        case("Ru")
            xs=6.210000_flyt
        case("Rh")
            xs=4.340000_flyt
        case("Pd")
            xs=4.390000_flyt
        case("Ag")
            xs=4.407000_flyt
        case("Cd")
            xs=3.040000_flyt
        case("In")
            xs=2.080000_flyt
        case("Sn")
            xs=4.871000_flyt
        case("Sb")
            xs=3.900000_flyt
        case("Te")
            xs=4.230000_flyt
        case("I")
            xs=3.500000_flyt
        case("Xe")
            xs=2.960000_flyt
        case("Cs")
            xs=3.690000_flyt
        case("Ba")
            xs=3.230000_flyt
        case("La")
            xs=8.530000_flyt
        case("Ce")
            xs=2.940000_flyt
        case("Pr")
            xs=2.640000_flyt
        case("Nd")
            xs=7.430000_flyt
        case("Pm")
            xs=1.0_flyt
        case("Sm")
            xs=0.422000_flyt
        case("Eu")
            xs=6.570000_flyt
        case("Gd")
            xs=29.300000_flyt
        case("Tb")
            xs=6.840000_flyt
        case("Dy")
            xs=35.900000_flyt
        case("Ho")
            xs=8.060000_flyt
        case("Er")
            xs=7.630000_flyt
        case("Tm")
            xs=6.280000_flyt
        case("Yb")
            xs=19.420000_flyt
        case("Lu")
            xs=6.530000_flyt
        case("Hf")
            xs=7.600000_flyt
        case("Ta")
            xs=6.000000_flyt
        case("W")
            xs=2.970000_flyt
        case("Re")
            xs=10.600000_flyt
        case("Os")
            xs=14.400000_flyt
        case("Ir")
            xs=14.100000_flyt
        case("Pt")
            xs=11.580000_flyt
        case("Au")
            xs=7.320000_flyt
        case("Hg")
            xs=20.240000_flyt
        case("Tl")
            xs=9.678000_flyt
        case("Pb")
            xs=11.115000_flyt
        case("Bi")
            xs=9.148000_flyt
        case("Po")
            xs=0.000000_flyt
        case("At")
            xs=0.000000_flyt
        case("Rn")
            xs=0.000000_flyt
        case("Fr")
            xs=0.000000_flyt
        case("Ra")
            xs=1.0_flyt
        case("Ac")
            xs=0.000000_flyt
        case("Th")
            xs=13.360000_flyt
        case("Pa")
            xs=10.400000_flyt
        case("U")
            xs=8.903000_flyt
        case("Np")
            xs=14.000000_flyt
        case("Pu")
            xs=1.0_flyt
        case("Am")
            xs=8.700000_flyt
        case("Cm")
            xs=0.000000_flyt
    end select
end function

end submodule
