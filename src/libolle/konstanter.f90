#include "precompilerdefinitions"
#include "gitinformation"
module konstanter
use, intrinsic :: iso_c_binding, only: c_float,c_double,c_int32_t,c_int64_t
use, intrinsic :: iso_fortran_env, only : output_unit
!!
!! Container for physical constants, tolerances and things like that. Note that I changed the name to swedish since when you try to link with other stuff, it's rather common to have a module named "constants".
!!
implicit none

!> precisions, consistent with the definitions in ELSI. Unclear if this is different than the iso_fortran_env.
integer, parameter :: r4 = c_float
integer, parameter :: r8 = c_double
integer, parameter :: i4 = c_int32_t
integer, parameter :: i8 = c_int64_t

!> Default output unit. I make it a 'saved' variable in case someone wants to redirect it.
integer, save :: lo_iou=output_unit

!> global floating point precision
integer,parameter :: flyt=r8 !dubbel

!> imaginary i
complex(flyt), parameter :: lo_imag = (0.0_flyt,1.0_flyt)
!> pi
real(flyt),parameter :: lo_pi=3.141592653589793_flyt
!> 2*pi
real(flyt),parameter :: lo_twopi=6.283185307179586_flyt

! Ok it's a little weird that I have redefined my constants lately, but they
! where not completely consistent. Now I went to NIST and have made sure it actually
! make sense, and is consistent with each other and so on. These are the
! recommended values from NIST as of November 2017

!> eV to Joule
real(flyt), parameter :: lo_eV_to_Joule=1.6021766208E-19_flyt
real(flyt), parameter :: lo_Joule_to_eV=1.0_flyt/lo_eV_to_Joule
!> eV to Hartree
real(flyt), parameter :: lo_Hartree_to_eV=27.21138602_flyt
real(flyt), parameter :: lo_eV_to_Hartree=1.0_flyt/lo_Hartree_to_eV
!> Hartree to Joule
real(flyt), parameter :: lo_Hartree_to_Joule=lo_Hartree_to_eV*lo_eV_to_Joule
real(flyt), parameter :: lo_Joule_to_Hartree=1.0_flyt/lo_Hartree_to_Joule
!> Bohr radius to m
real(flyt), parameter :: lo_bohr_to_m=0.52917721067E-10_flyt
real(flyt), parameter :: lo_m_to_bohr=1.0_flyt/lo_bohr_to_m
!> Bohr radius to angstrom
real(flyt), parameter :: lo_bohr_to_A=0.52917721067_flyt
real(flyt), parameter :: lo_A_to_bohr=1.0_flyt/lo_bohr_to_A
!> Planck constant
real(flyt), parameter :: lo_hbar_eV=6.582119514E-16_flyt
real(flyt), parameter :: lo_hbar_Joule=lo_hbar_eV*lo_eV_to_Joule
real(flyt), parameter :: lo_hbar_Hartree=lo_hbar_eV*lo_eV_to_Hartree
!> Boltzmann constant
real(flyt), parameter :: lo_kb_eV=8.6173303E-5_flyt
real(flyt), parameter :: lo_kb_Hartree=lo_kb_eV*lo_eV_to_Hartree
real(flyt), parameter :: lo_kb_Joule=lo_kb_eV*lo_eV_to_Joule
!> Speed of light in m/s
real(flyt), parameter :: lo_c_ms=299792458_flyt
!> Atomic and electron mass unit
real(flyt), parameter :: lo_amu_to_kg=1.660539040E-27_flyt
real(flyt), parameter :: lo_emu_to_kg=9.10938356E-31_flyt
real(flyt), parameter :: lo_amu_to_emu=lo_amu_to_kg/lo_emu_to_kg
real(flyt), parameter :: lo_emu_to_amu=1.0_flyt/lo_amu_to_emu

!> convert forces in eV/A to atomic units, Hartree/bohr
real(flyt), parameter :: lo_force_eVA_to_HartreeBohr=lo_eV_to_Hartree/lo_A_to_bohr
real(flyt), parameter :: lo_force_HartreeBohr_to_eVA=1.0_flyt/lo_force_eVA_to_HartreeBohr
!> convert volumes
real(flyt), parameter :: lo_volume_A_to_bohr=lo_A_to_bohr**3
real(flyt), parameter :: lo_volume_bohr_to_A=1.0_flyt/lo_volume_A_to_bohr
!> convert pressure in GPa to Hartree/bohr^3
real(flyt), parameter :: lo_pressure_GPa_to_HartreeBohr=1E9_flyt*(lo_bohr_to_m**3)/lo_Hartree_to_Joule
real(flyt), parameter :: lo_pressure_HartreeBohr_to_GPa=1.0_flyt/lo_pressure_GPa_to_HartreeBohr
!> convert pressure in eV/A to atomic units, Hartree/bohr
real(flyt), parameter :: lo_pressure_eVA_to_HartreeBohr=lo_eV_to_Hartree/(lo_A_to_bohr)**3
real(flyt), parameter :: lo_pressure_HartreeBohr_to_eVA=1.0_flyt/lo_pressure_eVA_to_HartreeBohr
!> convert time to atomic units
real(flyt), parameter :: lo_time_au_to_s=lo_hbar_Joule/lo_Hartree_to_Joule
real(flyt), parameter :: lo_time_s_to_au=1.0_flyt/lo_time_au_to_s
real(flyt), parameter :: lo_time_au_to_fs=1E15_flyt*lo_time_au_to_s
real(flyt), parameter :: lo_time_fs_to_au=1.0_flyt/lo_time_au_to_fs
!> convert velocities
real(flyt), parameter :: lo_velocity_au_to_ms=lo_bohr_to_m/lo_time_au_to_s
real(flyt), parameter :: lo_velocity_ms_to_au=1.0_flyt/lo_velocity_au_to_ms
real(flyt), parameter :: lo_velocity_au_to_Afs=lo_velocity_au_to_ms*1E-5_flyt
real(flyt), parameter :: lo_velocity_Afs_to_au=1.0_flyt/lo_velocity_au_to_Afs
!> convert forceconstants
real(flyt), parameter :: lo_forceconstant_1st_eVA_to_HartreeBohr=lo_eV_to_Hartree/(lo_A_to_bohr)
real(flyt), parameter :: lo_forceconstant_2nd_eVA_to_HartreeBohr=lo_eV_to_Hartree/(lo_A_to_bohr**2)
real(flyt), parameter :: lo_forceconstant_3rd_eVA_to_HartreeBohr=lo_eV_to_Hartree/(lo_A_to_bohr**3)
real(flyt), parameter :: lo_forceconstant_4th_eVA_to_HartreeBohr=lo_eV_to_Hartree/(lo_A_to_bohr**4)
real(flyt), parameter :: lo_forceconstant_1st_HartreeBohr_to_eVA=1.0_flyt/lo_forceconstant_1st_eVA_to_HartreeBohr
real(flyt), parameter :: lo_forceconstant_2nd_HartreeBohr_to_eVA=1.0_flyt/lo_forceconstant_2nd_eVA_to_HartreeBohr
real(flyt), parameter :: lo_forceconstant_3rd_HartreeBohr_to_eVA=1.0_flyt/lo_forceconstant_3rd_eVA_to_HartreeBohr
real(flyt), parameter :: lo_forceconstant_4th_HartreeBohr_to_eVA=1.0_flyt/lo_forceconstant_4th_eVA_to_HartreeBohr
!> convert phonon frequencies
real(flyt), parameter :: lo_frequency_Hartree_to_Hz=1.0_flyt/lo_hbar_Hartree
real(flyt), parameter :: lo_frequency_Hartree_to_THz=1E-12_flyt/lo_twopi/lo_hbar_Hartree
real(flyt), parameter :: lo_frequency_Hartree_to_meV=1000.0_flyt*lo_Hartree_to_eV
real(flyt), parameter :: lo_frequency_Hartree_to_icm=1.0_flyt/lo_hbar_Hartree/lo_twopi/lo_c_ms/100.0_flyt
real(flyt), parameter :: lo_frequency_THz_to_Hartree=1.0_flyt/lo_frequency_Hartree_to_THz
real(flyt), parameter :: lo_frequency_meV_to_Hartree=1.0_flyt/lo_frequency_Hartree_to_meV
!> convert group velocities
real(flyt), parameter :: lo_groupvel_Hartreebohr_to_ms=lo_bohr_to_m/lo_time_au_to_s
real(flyt), parameter :: lo_groupvel_ms_to_Hartreebohr=1.0_flyt/lo_groupvel_Hartreebohr_to_ms
!> thermal conductivity, atomic units to SI(W/mK)
real(flyt), parameter :: lo_kappa_au_to_SI=lo_Hartree_to_Joule/lo_time_au_to_s/lo_bohr_to_m
real(flyt), parameter :: lo_kappa_SI_to_au=1.0_flyt/lo_kappa_au_to_SI

! Some default tolerances
!> Tolerance for realspace distances to be 0
real(flyt), parameter :: lo_tol=1E-5_flyt
!> Tolerance for realspace squared distances to be 0
real(flyt), parameter :: lo_sqtol=lo_tol**2
!> Tolerance for reciprocal distances to be 0
real(flyt), parameter :: lo_rectol=1E-6_flyt
!> Tolerance for reciprocal squared distances
real(flyt), parameter :: lo_sqrectol=lo_rectol**2
!> Tolerance for angles in degrees to be 0
real(flyt), parameter :: lo_degreetol=1E-4_flyt
!> Tolerance for angles in radians to be 0
real(flyt), parameter :: lo_radiantol=lo_degreetol*180.0_flyt/lo_pi
!> Tolerance for phonon frequencies to be 0.
real(flyt), parameter :: lo_freqtol=lo_tol*1E-4_flyt
!> Tolerance for phonon group velocities to be zero
real(flyt), parameter :: lo_phonongroupveltol=lo_tol*1E-5_flyt
!> Tolerance for temperatures to be 0, in K
real(flyt), parameter :: lo_temperaturetol=1E-3_flyt
!> large number
real(flyt), parameter :: lo_huge=huge(1.0_flyt)
!> small number
real(flyt), parameter :: lo_tiny=tiny(1.0_flyt)
!> large integer
real(flyt), parameter :: lo_hugeint=huge(1)

!> Variable that holds exit status, for catching exceptions
integer, save :: lo_status=0
!> A strange vector I use to make some things more deterministic. Don't ask.
real(flyt), dimension(3), parameter :: lo_degenvector=[1.0_flyt,1.33_flyt,1.45624623_flyt]*1E-8_flyt

! The default gnuplot terminal, decided by precompiler flags.
#if GPwxt
!> default gnuplot terminal on linux
character(len=3), parameter :: lo_gnuplotterminal='wxt'
#elif GPaqua
!> default gnuplot terminal on osx
character(len=4), parameter :: lo_gnuplotterminal='aqua'
#elif GPqt
!> default gnuplot terminal on osx
character(len=4), parameter :: lo_gnuplotterminal='qt'
#else
!> fallback gnuplot terminal
character(len=3), parameter :: lo_gnuplotterminal='qt'
#endif

! Some exit codes for when things go horribly wrong.

!> Something has an unexpected dimension
integer, parameter :: lo_exitcode_baddim=1
!> BLAS or LAPACK fails
integer, parameter :: lo_exitcode_blaslapack=2
!> Unphysical value, i.e. negative temperature or something like that.
integer, parameter :: lo_exitcode_physical=3
!> Bad symmetry
integer, parameter :: lo_exitcode_symmetry=4
!> Something off with the arguments sent to the routine
integer, parameter :: lo_exitcode_param=5
!> IO error
integer, parameter :: lo_exitcode_io=6
!> MPI problem
integer, parameter :: lo_exitcode_mpi=7

! Some constant strings to use in the documentation
character(len=12), parameter :: lo_author='Olle Hellman'
character(len=3), parameter :: lo_version='1.2'
character(len=3), parameter :: lo_licence='MIT'
character(len=256), parameter :: lo_gitbranch=gitbranch
character(len=256), parameter :: lo_gitrevision=gitrevision
end module
