#include "precompilerdefinitions"
program generate_structure
!!{!src/generate_structure/manual.md!}
use konstanter, only: flyt, lo_pi, lo_sqtol, lo_A_to_bohr
use gottochblandat, only: open_file, tochar, lo_determ, lo_chop, lo_get_axis_angles
use geometryfunctions, only: lo_inscribed_sphere_in_box
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_forcemap, only: lo_forcemap
use type_crystalstructure, only: lo_crystalstructure
! use type_sqs, only: lo_sqs

use options, only: lo_opts
! use magneticdisorder, only: lo_magdisorder
use autocell, only: return_supercellmatrix
implicit none

type(lo_mpi_helper) :: mw
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc, ss, p
type(lo_forcemap) :: map
! type(lo_magdisorder) :: mag
! type(lo_sqs) :: sqs

! Set some options and parse input file
init: block
    call opts%parse()
    call mw%init()
    if (mw%talk .eqv. .false.) opts%verbosity = -100
    ! read positions
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    if (mw%talk) write (*, *) '... read unitcell'
end block init

! Create the supercell
getsupercell: block
    real(flyt), dimension(3, 3) :: tm
    real(flyt), parameter :: perfect_fill = 0.523598775598299_flyt
    real(flyt) :: a, b, c, al, be, gm, fillratio, r0
    integer, dimension(3, 3) :: supercellmatrix
    integer :: i, j, u
    character(len=2000) :: dumstr

    ! decide how to build a supercell, a few different routines are provided.
    if (sum(abs(opts%ndssdim)) .gt. 0.0_flyt) then
        ! non-diagonal supercell
        call uc%build_supercell(ss, nondiagdimensions=opts%ndssdim)
        supercellmatrix = opts%ndssdim
    elseif (opts%desired_na .gt. 0) then
        ! number of atoms instead of dimensions are defined
        call return_supercellmatrix(uc, opts%desired_na, supercellmatrix)
        call uc%build_supercell(ss, nondiagdimensions=supercellmatrix)
    else
        ! normal, diagonal supercell
        call uc%build_supercell(ss, opts%ssdim)
        supercellmatrix = 0
        do i = 1, 3
            supercellmatrix(i, i) = opts%ssdim(i)
        end do
    end if

    ! build a supercell
    if (mw%talk) write (*, *) '... built supercell'

    r0 = lo_inscribed_sphere_in_box(ss%latticevectors)
    fillratio = 100*(4*lo_pi*(r0**3)/(ss%volume*3.0_flyt))/perfect_fill
    call lo_get_axis_angles(ss%latticevectors, a, b, c, al, be, gm)

    if (mw%talk) then
        write (*, *) ' Supercellmatrix:'
        do i = 1, 3
            write (*, *) tochar(supercellmatrix(:, i))
        end do
        write (*, "(1X,'         Filling ratio:',(2X,F14.7,A))") fillratio, ' % of the ideal cube'
        write (*, "(1X,'           a,b,c (Ang):',3(2X,F14.7))") a/lo_A_to_bohr, b/lo_A_to_bohr, c/lo_A_to_bohr
        write (*, "(1X,'alpha,beta,gamma (deg):',3(2X,F14.7))") al*180/lo_pi, be*180/lo_pi, gm*180/lo_pi
        write (*, *) '       number of atoms: ', lo_determ(supercellmatrix)*uc%na
        write (*, *) ''

        select case (opts%outputformat)
        case (1) ! VASP
            call ss%writetofile('outfile.ssposcar', opts%outputformat)
            write (*, *) '... wrote supercell in VASP format'
        case (2) ! Abinit
            call ss%writetofile('outfile.supercell_abinit', opts%outputformat)
            write (*, *) '... wrote supercell in Abinit format'
        case (3) ! LAMMPS
            call lo_stop_gracefully(['Native LAMMPS IO was removed, please use external converters.'], 8)
        case (4) ! FHI-Aims
            call ss%writetofile('outfile.supercell_aims', opts%outputformat, transformationmatrix=tm)
            write (*, *) '... wrote supercell in FHI-Aims format'
        case (5) ! xyz format, for i-pi
            ! create that weird comment line that I-PI likes:
            write (*, *) 'FIXME ATOMIC UNITS IPI'
            stop
            dumstr = "# CELL{H}: "
            do i = 1, 3
            do j = 1, 3
                dumstr = trim(dumstr)//" "//tochar(ss%latticevectors(j, i), ndecimals=10)
            end do
            end do
            dumstr = trim(dumstr)//" cell{angstrom} positions{angstrom}"
            u = open_file('out', 'outfile.cell_ipi.xyz')
            write (u, *) ss%na
            write (u, *) trim(dumstr)
            do i = 1, ss%na
                write (u, "(2X,A4,4X,3(1X,E19.12))") trim(ss%atomic_symbol(ss%species(i))), lo_chop(ss%rcart(:, i), lo_sqtol)
            end do
            close (u)
            ! and print the normal output stuff
            call uc%writetofile('outfile.uc_ipi', 1)
            call ss%writetofile('outfile.ss_ipi', 1)
            ! and the lammps file
            ss%latticevectors = ss%latticevectors
            ss%inv_latticevectors = ss%inv_latticevectors
            call ss%writetofile('outfile.supercell_lammps_ipi', 3, transformationmatrix=tm)
        case (6) ! QE
            call ss%writetofile('outfile.supercell_qe', opts%outputformat)
            write (*, *) '... wrote supercell in Quantum ESPRESSO format'
        case (7) ! Parsec
            call uc%writetofile('outfile.unitcell_parsec', opts%outputformat)
            call ss%writetofile('outfile.supercell_parsec', opts%outputformat)
            write (*, *) '... wrote supercell in PARSEC format'
        end select
    end if
end block getsupercell

! ! If I want to do something disordered, this is where I generate that.
! disorderthings: block
!     integer :: i
!     ! Now stop if there is no alloy stuff going on.
!     if ( uc%info%alloy .or. uc%info%collmag .or. uc%info%noncollmag ) then
!         ! Best to generate a symmetry table
!         if ( opts%cutoff2 .lt. 0 ) opts%cutoff2=ss%maxcutoff()
!         call sl%generate(uc,ss,opts%cutoff2,-1.0_flyt,-1.0_flyt,verbosity=opts%verbosity,polar=.false.,&
!                          wraparound=.true.,magcutoff2=opts%cutoff2)
!         call map%generate(sl,uc,ss,polarcorrectiontype=0,verbosity=opts%verbosity)
!     else
!         ! If neither an alloy or magnetically disordered, stop here.
!         if ( mw%talk ) write(*,*) '... done'
!         call mw%destroy()
!         stop
!     endif
!
!     ! Generate magnetically disordered things
!     if ( mw%talk ) then
!         if ( uc%info%collmag .or. uc%info%noncollmag ) then
!             ! Get a bunch of magnetic configurations
!             call mag%generate(sl,uc,ss)
!             call mag%optimize(ss,opts%magnconf,opts%magnbin)
!             ! Dump them
!             call mag%dump_configurations(ss)
!         endif
!     endif
!
!     if ( uc%info%alloy ) then
!         ! Generate an sqs, optimize on randomness and so on
!         call sqs%generate(uc,ss,sl,max(opts%nsqs,5),opts%verbosity,mw)
!         ! Dump a configuration
!         do i=1,opts%nsqs
!             call sqs%returnstructure(ss,p,i,.true.,opts%verbosity)
!             call p%classify('bravais')
!             call p%writetofile('outfile.sqs_'//tochar(i,3),1)
!         enddo
!     endif
!
! end block disorderthings

call mw%destroy()

end program
