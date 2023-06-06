#include "precompilerdefinitions"
program convert_abinit_ddb_to_forceconstant
!!{!src/convert_abinit_ddb_to_forceconstant/manual.md!}
use konstanter, only: flyt, lo_huge, lo_tol, lo_sqtol, lo_pressure_HartreeBohr_to_GPa
use gottochblandat, only: lo_chop, lo_clean_fractional_coordinates, tochar, open_file
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
!use type_symmetrylist, only: lo_symlist
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_forcemap, only: lo_forcemap
use type_qpointmesh, only: lo_qpoint

use options, only: lo_opts
use io, only: read_ddb_file
!use tofc, only: forceconstant_from_dynmat
use ifc_solvers, only: lo_irreducible_forceconstant_from_qmesh_dynmat, lo_solve_for_borncharges

implicit none
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc
complex(flyt), dimension(:, :, :), allocatable :: dynmat
real(flyt), dimension(:, :, :), allocatable :: born_effective_charges
real(flyt), dimension(:, :), allocatable :: qvecs
real(flyt), dimension(3, 3) :: dielectric_tensor

! Grab stuff from DDB files
init: block
    ! Read options
    call mw%init()
    call mem%init()
    call opts%parse()
    if (mw%talk .eqv. .false.) opts%verbosity = -100
    ! Read everything from DDB files
    call read_ddb_file(opts%filename, uc, dynmat, qvecs, born_effective_charges, dielectric_tensor, opts%verbosity)
end block init

! Get a forceconstant from this
solve: block
    type(lo_interaction_tensors) :: sl
    type(lo_forcemap) :: map
    type(lo_crystalstructure) :: ss
    type(lo_forceconstant_secondorder) :: fc
    type(lo_qpoint), dimension(:), allocatable :: qp
    complex(flyt), dimension(:, :, :, :), allocatable :: wdm
    complex(flyt), dimension(:, :, :), allocatable :: lrdynmat, ldynmat
    real(flyt), dimension(3) :: v0
    real(flyt) :: f0
    integer, dimension(3) :: qdim
    integer :: i, j, a1, a2, ii, jj, q, qq, nq, u
    logical :: polar

    polar = .false.
    ! Now ... try to figure out the q-dimensions. Not too shabby. I'm fairly certain it's quite stable.
    v0 = lo_huge
    do j = 1, 3
    do i = 1, size(qvecs, 2)
        if (abs(qvecs(j, i)) .gt. lo_sqtol) v0(j) = min(v0(j), qvecs(j, i))
    end do
    end do
    qdim = int(anint(1.0_flyt/v0))
    ! Only override if no qgrid is specified.
    if ( opts%qgrid(1) .eq. -1 ) opts%qgrid=qdim
    write(*,*) "Supercell Size : "
    write(*,*) opts%qgrid
    ! Get a supercell
    call uc%build_supercell(ss, opts%qgrid)
    ! Figure out if it's polar or not
    if (opts%forcenopolar) then
        ! nope, no polar corrections
        polar = .false.
    elseif (sum(abs(dielectric_tensor)) .gt. lo_sqtol .and. sum(abs(born_effective_charges)) .gt. lo_sqtol) then
        polar = .true.
        if (mw%talk) then
            write (*, *) '... found born charges and dielectric tensor'
            if (opts%verbosity .gt. 0) then
                write (*, *) 'raw dielectric tensor:'
                do i = 1, 3
                    write (*, *) dielectric_tensor(:, i)
                end do
                do j = 1, uc%na
                    write (*, *) 'raw Born effective charge atom ', tochar(j)
                    do i = 1, 3
                        write (*, *) born_effective_charges(:, i, j)
                    end do
                end do
            end if
        end if
    else
        polar = .false.
    end if
    ! Get the symmetry stuff
    if (opts%cutoff > 0.0_flyt) then
        f0 = opts%cutoff
        call sl%generate(uc, ss, cutoff2=f0, cutoff3=-1.0_flyt, cutoff4=-1.0_flyt, polar=polar, transposition=.true., &
                         spacegroup=.true., verbosity=opts%verbosity + 1, firstorder=.false., wzdim=[-1, -1, -1], mw=mw, mem=mem)
    elseif (opts%truncate) then
        f0 = ss%maxcutoff() !*0.25
        call sl%generate(uc, ss, cutoff2=f0, cutoff3=-1.0_flyt, cutoff4=-1.0_flyt, polar=polar, transposition=.true., &
                         spacegroup=.true., verbosity=opts%verbosity + 1, firstorder=.false., wzdim=[-1, -1, -1], mw=mw, mem=mem)
    else
        call sl%generate(uc, ss, cutoff2=-1.0_flyt, cutoff3=-1.0_flyt, cutoff4=-1.0_flyt, polar=polar, transposition=.true., &
                         spacegroup=.true., verbosity=opts%verbosity + 1, firstorder=.false., wzdim=[-2, -2, -2], mw=mw, mem=mem)
    end if

    ! And the forcemap
    call map%generate(uc, ss, polarcorrectiontype=3, st=sl, mw=mw, mem=mem, verbosity=opts%verbosity + 1)
    !call map%generate(sl,uc,ss,polarcorrectiontype=3,verbosity=opts%verbosity+1)
    call map%forceconstant_constraints(uc, .true., .true., .true., opts%verbosity + 1)

    ! Might as well do this thing in parallel. Could get slow for many q-points otherwise.
    nq = 0
    do i = 1, size(qvecs, 2)
        if (mod(i, mw%n) .eq. mw%r) nq = nq + 1
    end do

    ! Pack the qpoints and dynamical matrices the way I like it, per rank
    if (nq .gt. 0) then
        lo_allocate(qp(nq))
        lo_allocate(ldynmat(uc%na*3, uc%na*3, nq))
        lo_allocate(lrdynmat(uc%na*3, uc%na*3, nq))
        if (map%polar .gt. 0) then
            call lo_solve_for_borncharges(map, p=uc, Z=born_effective_charges, eps=dielectric_tensor, verbosity=opts%verbosity, mw=mw, mem=mem)
            call map%get_secondorder_forceconstant(uc, fc, verbosity=opts%verbosity, mem=mem)
            lo_allocate(wdm(3, 3, uc%na, uc%na))
            wdm = 0.0_flyt
        end if
        ldynmat = 0.0_flyt
        lrdynmat = 0.0_flyt
        q = 0
        do qq = 1, size(qvecs, 2)
            if (mod(qq, mw%n) .ne. mw%r) cycle
            ! Grab the q-points
            q = q + 1
            v0 = lo_clean_fractional_coordinates(qvecs(:, qq))
            v0 = lo_chop(v0, 1E-12_flyt)
            v0 = matmul(uc%reciprocal_latticevectors, v0)
            qp(q)%r = v0 - uc%bz%gshift(v0)
            qp(q)%n_invariant_operation = 0
            ! Store the relevant dynamical matrices for this rank, first the shortrange
            ldynmat(:, :, q) = dynmat(:, :, qq)
            ! Then maybe longrange
            if (map%polar .gt. 0) then
                call fc%longrange_dynamical_matrix(wdm, uc, qp(q)%r) !,reconly=.true.)
                do a1 = 1, uc%na
                do a2 = 1, uc%na
                    do i = 1, 3
                    do j = 1, 3
                        ii = (a1 - 1)*3 + i
                        jj = (a2 - 1)*3 + j
                        lrdynmat(jj, ii, q) = wdm(i, j, a1, a2)
                    end do
                    end do
                end do
                end do
            end if
        end do
    end if

    if (mw%talk) write (*, *) '... rearranged dynamical matrices'

    ! Solve backwards
    if (opts%truncate) then
        call lo_irreducible_forceconstant_from_qmesh_dynmat(map, uc, qp, nq, ldynmat, lrdynmat, &
                                                            fullhavemasses=.false., longrangehavemasses=.false., mw=mw, verbosity=opts%verbosity + 1, enforce=opts%truncate)
    elseif (opts%cutoff > 0) then
        call lo_irreducible_forceconstant_from_qmesh_dynmat(map, uc, qp, nq, ldynmat, lrdynmat, &
                                                            fullhavemasses=.false., longrangehavemasses=.false., mw=mw, verbosity=opts%verbosity + 1, enforce=opts%enforcesym)
    else
        call lo_irreducible_forceconstant_from_qmesh_dynmat(map, uc, qp, nq, ldynmat, lrdynmat, &
                                                            fullhavemasses=.false., longrangehavemasses=.false., mw=mw, verbosity=opts%verbosity + 1, enforce=opts%enforcesym)
    end if

    ! Then dump
    if (mw%talk) then
        call map%get_secondorder_forceconstant(uc, fc, mem, 0)
        call fc%get_elastic_constants(uc)
        write (*, *) ''
        write (*, *) 'elastic constants (GPa):'
        do i = 1, 6
            write (*, "(6(3X,F15.5))") lo_chop(fc%elastic_constants_voigt(:, i)*lo_pressure_HartreeBohr_to_GPa, lo_tol)
        end do
        ! Dump forceconstant
        call uc%writetofile('outfile.ucposcar', 1)
        write (*, *) '... wrote unitcell'
        call fc%writetofile(uc, 'outfile.forceconstant')
        write (*, *) '... wrote forceconstant'
        ! Also dump the Born effective charges+dielectric tensor
        if (map%polar .gt. 0) then
            u = open_file('out', 'outfile.lotosplitting')
            do i = 1, 3
                write (u, *) fc%loto%eps(:, i)
            end do
            do j = 1, fc%na
            do i = 1, 3
                write (u, *) fc%loto%born_effective_charges(:, i, j)
            end do
            end do
            close (u)
            write (*, *) '... wrote dielectric tensor and Born charges'
        end if

        if (opts%verbosity .gt. 0 .and. map%polar .gt. 0) then
            write (*, *) ''
            write (*, *) 'Dielectric tensor (raw, symmetric)'
            do i = 1, 3
                write (*, "(1X,3(1X,F12.6),2X,3(1X,F12.6))") dielectric_tensor(:, i), fc%loto%eps(:, i)
            end do
            do j = 1, fc%na
                write (*, *) 'Born effective charge atom ', tochar(j), ' (raw,symmetric)'
                do i = 1, 3
                    write (*, "(1X,3(1X,F12.6),2X,3(1X,F12.6))") born_effective_charges(:, i, j), fc%loto%born_effective_charges(:, i, j)
                end do
            end do

        end if
    end if
end block solve

if (mw%talk) write (*, *) 'All done!'
call mw%destroy()

end program
