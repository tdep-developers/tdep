#include "precompilerdefinitions"
program convert_XX_to_forceconstant
!!{!src/convert_XX_to_forceconstant/manual.md!}
use konstanter, only: r8, lo_huge, lo_tol, lo_sqtol, lo_pressure_HartreeBohr_to_GPa
use gottochblandat, only: lo_chop, lo_clean_fractional_coordinates, tochar, open_file, lo_unflatten_2tensor, lo_flattentensor
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
!use type_symmetrylist, only: lo_symlist
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_forcemap, only: lo_forcemap,lo_coeffmatrix_unitcell_Z_singlet
use type_qpointmesh, only: lo_qpoint

use options, only: lo_opts
use io, only: read_ddb_file
!use tofc, only: forceconstant_from_dynmat
use ifc_solvers, only: lo_irreducible_forceconstant_from_qmesh_dynmat, lo_solve_for_borncharges
use lo_longrange_electrostatics, only: lo_ewald_parameters

implicit none
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc
complex(r8), dimension(:, :, :), allocatable :: dynmat
real(r8), dimension(:, :, :), allocatable :: born_effective_charges
real(r8), dimension(:, :), allocatable :: qvecs
real(r8), dimension(3, 3) :: dielectric_tensor

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
    complex(r8), dimension(:, :, :, :), allocatable :: wdm
    complex(r8), dimension(:, :, :), allocatable :: lrdynmat, ldynmat
    real(r8), dimension(3) :: vv
    real(r8) :: f0
    integer, dimension(3) :: qdim
    integer :: i, j, a1, a2, ii, jj, q, qq, nq, u
    logical :: polar

    polar = .false.
    ! Now ... try to figure out the q-dimensions. Not too shabby. I'm fairly certain it's quite stable.
    vv = lo_huge
    do j = 1, 3
    do i = 1, size(qvecs, 2)
        if (abs(qvecs(j, i)) .gt. lo_sqtol) vv(j) = min(vv(j), qvecs(j, i))
    end do
    end do
    qdim = int(anint(1.0_r8/vv))
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
    if (opts%cutoff > 0.0_r8) then
        f0 = opts%cutoff
        call sl%generate(uc, ss, cutoff2=f0, cutoff3=-1.0_r8, cutoff4=-1.0_r8, polar=polar, transposition=.true., &
                         spacegroup=.true., verbosity=opts%verbosity + 1, firstorder=.false., wzdim=[-1, -1, -1], mw=mw, mem=mem)
    elseif (opts%truncate) then
        f0 = ss%maxcutoff() !*0.25
        call sl%generate(uc, ss, cutoff2=f0, cutoff3=-1.0_r8, cutoff4=-1.0_r8, polar=polar, transposition=.true., &
                         spacegroup=.true., verbosity=opts%verbosity + 1, firstorder=.false., wzdim=[-1, -1, -1], mw=mw, mem=mem)
    else
        call sl%generate(uc, ss, cutoff2=-1.0_r8, cutoff3=-1.0_r8, cutoff4=-1.0_r8, polar=polar, transposition=.true., &
                         spacegroup=.true., verbosity=opts%verbosity + 1, firstorder=.false., wzdim=[-2, -2, -2], mw=mw, mem=mem)
    end if

    ! And the forcemap
    call map%generate(uc, ss, polarcorrectiontype=3, st=sl, mw=mw, mem=mem, verbosity=opts%verbosity + 1)
    !call map%generate(sl,uc,ss,polarcorrectiontype=3,verbosity=opts%verbosity+1)

    ! As soon as we know the Born charges we can construct the IFC constraints
    buildconstraints: block
        type(lo_ewald_parameters) :: ew
        complex(r8), dimension(:,:,:,:), allocatable :: D0
        real(r8), dimension(:,:,:,:), allocatable :: rottensor
        real(r8), dimension(:,:,:), allocatable :: wZ,wD
        real(r8), dimension(:,:), allocatable :: coeff_Z
        real(r8), dimension(:), allocatable :: v0
        real(r8), dimension(:), allocatable :: hermitian_rhs, huang_rhs, rotational_rhs
        real(r8), dimension(3,3,3,3) :: bracket
        real(r8), dimension(3,3) :: eps,m0,m1
        integer :: k,l

        allocate(hermitian_rhs(9*uc%na))
        allocate(rotational_rhs(27*uc%na))
        allocate(huang_rhs(81))
        hermitian_rhs=0.0_r8
        rotational_rhs=0.0_r8
        huang_rhs=0.0_r8

        if ( map%polar .gt. 0 ) then
            ! First thing we need is the polar dynamical matrix to know
            ! how we should build the Hermiticity constraints.
            !eps=lo_unflatten_2tensor(matmul(map%eps_global_shell%coeff, map%xuc%x_eps_global))
            eps = dielectric_tensor

            allocate(wZ(3,3,uc%na))
            allocate(wD(3,3,uc%na))
            allocate(v0(9*uc%na))
            allocate(coeff_Z(uc%na*9,map%xuc%nx_Z_singlet))
            allocate(D0(3,3,uc%na,uc%na))
            wZ=0.0_r8
            wD=0.0_r8
            v0=0.0_r8
            coeff_Z=0.0_r8
            D0=0.0_r8

            call lo_coeffmatrix_unitcell_Z_singlet(map, coeff_Z)
            v0=matmul(coeff_Z,map%xuc%x_Z_singlet)
            do a1 = 1, uc%na
                wZ(:, :, a1) = lo_unflatten_2tensor(v0((a1 - 1)*9 + 1:a1*9))
            end do

            call ew%set(uc, eps, 2, 1E-20_r8, verbosity=-1)
            call ew%longrange_dynamical_matrix(uc, [0.0_r8, 0.0_r8, 0.0_r8], wZ, wD, eps, D0, reconly=.true., chgmult=.true.)

            do a1=1,uc%na
                m0=0.0_r8
                m1=0.0_r8
                do a2=1,uc%na
                    m0=m0+real(d0(:,:,a1,a2),r8)
                    m1=m1+real(d0(:,:,a2,a1),r8)
                enddo
                m0=m0-m1
                hermitian_rhs( (a1-1)*9+1:a1*9 ) = -lo_flattentensor(m0)
            enddo
            hermitian_rhs = lo_chop(hermitian_rhs,1E-12_r8)

            ! Next up we do the huang + rotational invariances
            allocate(rottensor(3,3,3,uc%na))
            rottensor=0.0_r8
            call ew%longrange_elastic_constant_bracket( uc, wZ, eps, bracket, rottensor, reconly=.true., mw=mw, verbosity=opts%verbosity)

            huang_rhs=0.0_r8
            ii=0
            do i=1,3
            do j=1,3
            do k=1,3
            do l=1,3
                ii=ii+1
                huang_rhs(ii) = bracket(k,l,i,j) - bracket(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            huang_rhs = lo_chop(huang_rhs,1E-12_r8)

            ii=0
            do a1=1,uc%na
            do i=1,3
            do j=1,3
            do k=1,3
                ii=ii+1
                rotational_rhs(ii) = rottensor(i,j,k,a1)-rottensor(i,k,j,a1)
            enddo
            enddo
            enddo
            enddo
            rotational_rhs = lo_chop(rotational_rhs,1E-12_r8)

        else
            ! We just have zeros
            huang_rhs=0.0_r8
            hermitian_rhs=0.0_r8
            rotational_rhs=0.0_r8
        endif

        call map%forceconstant_constraints(uc,rotational=.true.,huanginvariances=.true.,hermitian=.true.,hermitian_rhs=hermitian_rhs,huang_rhs=huang_rhs,rotational_rhs=rotational_rhs,verbosity=opts%verbosity)
    end block buildconstraints

!    call map%forceconstant_constraints(uc, rotational=.true., huanginvariances=.true., hermitian=.true., &
!&        hermitian_rhs=hermitian_rhs,huang_rhs=huang_rhs,rotational_rhs=rotational_rhs, verbosity=opts%verbosity + 1)


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
            wdm = 0.0_r8
        end if
        ldynmat = 0.0_r8
        lrdynmat = 0.0_r8
        q = 0
        do qq = 1, size(qvecs, 2)
            if (mod(qq, mw%n) .ne. mw%r) cycle
            ! Grab the q-points
            q = q + 1
            vv = lo_clean_fractional_coordinates(qvecs(:, qq))
            vv = lo_chop(vv, 1E-12_r8)
            vv = matmul(uc%reciprocal_latticevectors, vv)
            qp(q)%r = vv - uc%bz%gshift(vv)
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
        call fc%get_elastic_constants(uc,mw,-1)
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
