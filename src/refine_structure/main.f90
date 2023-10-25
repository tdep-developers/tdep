program refine_structure
!!{!src/refine_structure/manual.md!}
use konstanter, only: r8
use gottochblandat, only: tochar, walltime
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_symset
use lo_memtracker, only: lo_mem_helper

use options, only: lo_opts
use refine !, only: refine_unitcell,refine_supercell,refine_lattice,refine_one_cell
use lo_spacegroup, only: lo_symmetry_group
implicit none
! for normal phonon things
type(lo_opts) :: opts
type(lo_mem_helper) :: mem
type(lo_crystalstructure) :: p, uc, ss
type(lo_symset) :: sym
type(lo_symmetry_group) :: sg

logical :: supercell
real(r8) :: timer

! Get the command line arguments
timer = walltime()
call opts%parse()
call mem%init()

write (*, *) ''
write (*, *) 'REFINING STRUCTURE'
write (*, *) '    unitcell: ', trim(opts%unitcell_filename)
write (*, *) '   supercell: ', trim(opts%supercell_filename)
write (*, *) '   prototype: ', trim(opts%prototype_unitcell)

! Read the unitcell with broken symmetry
call uc%readfromfile(trim(opts%unitcell_filename), verbosity=opts%verbosity)

write (*, *) 'Read unit cell'

! Refining just a cell
call refine_one_cell(uc, p, opts%tolerance_lattice, opts%tolerance_internal, mem, 1)
call p%writetofile('outfile.refined_cell', 1)
write (*, *) '... wrote refined      cell to "outfile.refined_cell"'

tst1: block
    real(r8), dimension(:, :, :), allocatable :: ops
    real(r8), dimension(:, :), allocatable :: tr

    !call sg%generate_pool_of_point_operations(uc%latticevectors,ops,mem,1)
    !call sg%generate_pool_of_translations(uc%latticevectors,uc%r,uc%species,ops,tr,mem,1)

end block tst1

! Read the supercell, if applicable
if (trim(opts%supercell_filename) .ne. 'none') then
    call ss%readfromfile(trim(opts%supercell_filename))
    supercell = .true.
else
    supercell = .false.
end if

! Maybe fix the supercell?
if (supercell) then
    call refine_supercell(p, ss)
    call ss%writetofile('outfile.refined_supercell', 1)
    write (*, *) '... wrote refined supercell to "outfile.refined_supercell"'
end if

!write(*,*) 'FIXME PRESERVE HEADER'

! Now, is the proposed unitcell actually a supercell of some other cell? That is a little annoying, but fixable.
! call find_true_unitcell(uc)
! call uc%classify('bravais')
! stop
! Is this a supercell?

!
! ! Read the prototype?
! if ( trim(opts%prototype_unitcell) .ne. 'none' ) then
!     ! Read the prototype, and grab the symmetry operations from there
!     call p%readfromfile(trim(opts%prototype_unitcell))
!     p%info%verbosity=2
!     call p%classify('bravais',tolerance=opts%tolerance,refine=.false.)
!     call sym%generate(p%latticevectors,.true.,p%r,p%atomic_number,tolerance=opts%tolerance)
!     p%info%verbosity=0
!     ! And insert this basis into the unitcell
!     uc%latticevectors=p%latticevectors
!     uc%inv_latticevectors=p%inv_latticevectors
! else
!     uc%info%verbosity=2
!     call uc%classify('bravais',tolerance=opts%tolerance,refine=.false.)
!     call sym%generate(uc%latticevectors,.true.,uc%r,uc%atomic_number,tolerance=opts%tolerance)
!     uc%info%verbosity=0
! endif
!
! write(*,*) '... found ',tochar(sym%n),' symmetry operations'
!
! ! Now I have neat latticevectors, and some symmetryoperations. Make sure
! ! the unitcell is properly symmetric
! call refine_unitcell(uc,sym)
! ! Dump it
! call uc%writetofile('outfile.refined_cell',1)
! write(*,*) '... wrote refined unitcell to "outfile.refined_cell"'
!
! ! Maybe fix the supercell?
! if ( supercell ) then
!     call refine_supercell(uc,ss)
!     call ss%writetofile('outfile.refined_supercell',1)
!     write(*,*) '... wrote refined supercell to "outfile.refined_supercell"'
! endif
!
! ! Get the primitive cell, might be of interest?
! ! findprim: block
! !    call return_primitive_cell(uc)
! ! end block findprim
!
! write(*,*) 'All done! (',tochar(walltime()-timer),'s)'

end program
