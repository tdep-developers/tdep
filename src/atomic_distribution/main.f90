program atomic_distribution
!!{!src/atomic_distribution/manual.md!}
use konstanter, only: r8
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: lo_does_file_exist
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use lo_symmetry_of_interactions, only: lo_interaction_tensors

use options, only: lo_opts

use pairmapping
use pair_distribution
use mean_square_displacement
use diffraction

!use vectordist
!use timedistance_correlation
!use correlationfunction
implicit none

type(lo_opts) :: opts
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_crystalstructure) :: uc, ss
type(lo_mdsim) :: sim
type(lo_pairmapping) :: pm
type(lo_pair_distribution) :: pdf
type(lo_powderdiffraction) :: df
type(lo_mean_square_displacement) :: msd

! Fetch the normal things
init: block
    type(lo_interaction_tensors) :: sl

    ! get command line arguments
    call opts%parse()
    call mw%init()
    if (.not. mw%talk) opts%verbosity = -100
    call mem%init()
    ! read positions
    call uc%readfromfile('infile.ucposcar')
    call uc%classify('wedge', timereversal=.true.)
    if (mw%talk) write (*, *) '... read unitcell'
    call ss%readfromfile('infile.ssposcar')
    if (mw%talk) write (*, *) '... read supercell'

    ! Get all kinds of symmetry stuff
    if (opts%cutoff .lt. 0.0_r8) then
        opts%cutoff = ss%maxcutoff()
    end if

    ! Get some symmetry things
    call sl%generate(uc, ss, opts%cutoff, -1.0_r8, -1.0_r8, .false., mw, mem, opts%verbosity)

    call pm%setup_symmetry(sl, uc, ss, mw, mem, opts%verbosity + 1)

    ! Also need the actual simulation
    if (lo_does_file_exist('infile.sim.hdf5')) then
        call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=opts%stride) !,mw=mw)
        call sim%get_nonpbc_positions()
    else
        call sim%read_from_file(verbosity=opts%verbosity + 2, stride=opts%stride, dynamics=.true.)
    end if

    !call sim%remove_force_and_center_of_mass_drift()

end block init

! Then we calculate the mean square displacement, could be useful
call msd%generate(uc, ss, sim, mw)
if (mw%talk) then
    call msd%write_to_hdf5()
    call msd%write_to_plaintext(uc)
    write (*, *) 'Wrote mean square displacement to file'
end if

! Powder diffraction pattern
call df%generate(uc, ss, sim, mw, opts%verbosity + 1)
if (mw%talk) then
    call df%write_to_hdf5('outfile.powder_diffraction.hdf5')
end if

! First we calculate the symmetry-projected radial pair distribution function
call pdf%bin(uc, ss, pm, sim, opts%nbin, opts%promise_no_diffusion, mw, mem, opts%verbosity + 1)
! Write this to file
if (mw%talk) then
    call pdf%write_to_hdf5(pm, uc, sim, 'outfile.pair_distribution.hdf5')
    write (*, *) 'Wrote pair distribution function to file'
end if

! ! Now revive the old-school atomic distribution. It was more useful than I initially thought.
!
! ! ! generate the vector distribution (only in case of no diffusion, this makes little sense if things are melted)
! ! if ( pdf%diffusion .eqv. .false. ) then
! !     call vd%generate(pm,sim,opts%bintype,opts%transform,opts%nbin)
! !     call vd%write_to_hdf5(pm)
! ! else
! !     ! if it's melted, go with the probability density instead!
! !     write(*,*) 'Probability density function, but not yet. Nag on me to fix this.'
! !     stop
! ! endif

if (opts%verbosity .gt. 0) write (*, *) 'All done!'
call mw%destroy()

end program
