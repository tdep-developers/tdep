#include "precompilerdefinitions"
program convert_phonopy_to_forceconstant
!!{!src/convert_phonopy_to_forceconstant/manual.md!}
use konstanter, only: flyt, lo_status, lo_tol, lo_pressure_HartreeBohr_to_GPa
use gottochblandat, only: lo_does_file_exist, lo_linear_least_squares, lo_chop
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_symmetrylist, only: lo_symlist
use type_forcemap, only: lo_forcemap
use ifc_solvers, only: lo_irreducible_forceconstant_from_supercell_dynmat, lo_solve_for_borncharges

use options, only: lo_opts
use io

implicit none

type(lo_mpi_helper) :: mw
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc, ss
type(lo_forceconstant_secondorder) :: fc
real(flyt), dimension(:, :, :, :), allocatable :: phonopy_fc
real(flyt), dimension(:, :, :), allocatable :: born_charges
real(flyt), dimension(3, 3) :: dielectric_tensor
integer, dimension(:), allocatable :: ssind
logical :: polar

! Grab some options
init: block
    ! init MPI
    call mw%init()
    ! Parse options
    call opts%parse()
    if (mw%talk .eqv. .false.) opts%verbosity = -100
    ! Get the unitcell
    call uc%readfromfile(trim(opts%uc_filename), verbosity=0)
    ! Get the symmetry
    call uc%classify('wedge', timereversal=.true.)
    if (mw%talk) write (*, *) '... read unitcell'
    ! And the supercell
    call ss%readfromfile(trim(opts%ss_filename), verbosity=0)
    ! Sort out the relation between unitcell and supercell
    call ss%classify('supercell', uc)
    if (mw%talk) write (*, *) '... read supercell'
    ! Grab the phonopy forceconstant
    call read_new_phonopy_forceconstant(ss, trim(opts%fc_filename), phonopy_fc, ssind)
    if (mw%talk) write (*, *) '... read phonopy forceconstant'
    ! And grab born charges and dielectric tensor, if applicable
    if (lo_does_file_exist(trim(opts%born_filename))) then
        call read_phonopy_born_file(uc, trim(opts%born_filename), born_charges, dielectric_tensor)
        polar = .true.
        if (mw%talk) write (*, *) '... read born effective charges'
    else
        polar = .false.
    end if
end block init

! Another attempt at solving it backwards
bwards: block
    type(lo_symlist) :: sl
    type(lo_forcemap) :: map
    real(flyt), dimension(:, :, :, :), allocatable :: polar_fc, contracted_polar_fc
    real(flyt) :: f0
    integer :: i

    ! First figure out all the symmetries
    if (opts%truncate) then
        f0 = ss%maxcutoff()
        call sl%generate(uc, ss, f0, -1.0_flyt, -1.0_flyt, polar=polar, transposition=.true., spacegroup=.true., &
                         verbosity=opts%verbosity + 2, wzdim=[-1, -1, -1], nj2=-1, nj3=-1, nj4=-1, firstorder=.false., &
                         magcutoff2=-1.0_flyt, magsinglet=.false.)
    else
        call sl%generate(uc, ss, -1.0_flyt, -1.0_flyt, -1.0_flyt, polar=polar, transposition=.true., spacegroup=.true., &
                         verbosity=opts%verbosity + 2, wzdim=[-2, -2, -2], nj2=-1, nj3=-1, nj4=-1, firstorder=.false., &
                         magcutoff2=-1.0_flyt, magsinglet=.false.)
    end if
    call map%generate(sl, uc, ss, polarcorrectiontype=3, verbosity=opts%verbosity + 2)
    if (opts%truncate) then
        call map%forceconstant_constraints(uc, .true., .true., .true., opts%verbosity + 10)
    else
        call map%forceconstant_constraints(uc, .true., .true., .true., opts%verbosity + 10)
        !call map%forceconstant_constraints(uc,.false.,.false.,.false.,opts%verbosity+10)
    end if

    if (mw%talk) write (*, *) '... figured out symmetry things'

    ! Get longrange forceconstant, if necessary
    lo_allocate(contracted_polar_fc(3, 3, ss%na, size(ssind)))
    contracted_polar_fc = 0.0_flyt
    if (polar) then
        lo_allocate(polar_fc(3, 3, ss%na, ss%na))
        polar_fc = 0.0_flyt
        call lo_solve_for_borncharges(map, Z=born_charges, eps=dielectric_tensor, verbosity=opts%verbosity)
        map%ifc_pair = 1.0_flyt
        call map%get_secondorder_forceconstant(uc, fc)
        if (mw%r .eq. mw%n - 1) then
            call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, polar_fc, 1E-15_flyt)
        end if
        call mw%bcast(polar_fc, mw%n - 1, __FILE__, __LINE__)
        do i = 1, size(ssind)
            contracted_polar_fc(:, :, :, i) = polar_fc(:, :, :, ssind(i))
        end do
        lo_deallocate(polar_fc)
        if (mw%talk) write (*, *) '... got longrange forceconstant'
    end if

    ! Solve for the irreducible
    call lo_irreducible_forceconstant_from_supercell_dynmat(map, ss, dynmat_T=phonopy_fc, lrdynmat_T=contracted_polar_fc, &
                                                            enforce_constraints=opts%truncate, fullhavemass=.false., longrangehavemass=.false., &
                                                            subset=ssind, mw=mw, verbosity=opts%verbosity)

    if (mw%talk) write (*, *) '... solved for forceconstants'
    if (mw%talk) then
        call map%get_secondorder_forceconstant(uc, fc)
        call fc%get_elastic_constants(uc)
        write (*, *) 'elastic constants (GPa):'
        do i = 1, 6
            write (*, "(6(3X,F15.5))") lo_chop(fc%elastic_constants_voigt(:, i)*lo_pressure_HartreeBohr_to_GPa, lo_tol)
        end do

        if (map%constraints%neq2 .gt. 0) then
            f0 = sum(abs(matmul(map%constraints%eq2, map%ifc_pair)))
        else
            f0 = 0.0_flyt
        end if
        write (*, *) '    viol', f0, '(violation of symmetry)'
    end if
    ! Finally, dump it to file
    call fc%writetofile(uc, 'outfile.converted_forceconstant')

end block bwards

if (mw%talk) write (*, *) 'All done!'
call mw%destroy()

end program
