program phasespace_surface
use konstanter, only: flyt,lo_status
use gottochblandat, only: walltime,tochar
use mpi_wrappers, only: lo_mpi_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh,lo_wedge_mesh,lo_qpoint,lo_generate_qmesh,lo_read_qmesh_from_file
use type_phonon_dispersions, only: lo_phonon_dispersions

use options, only: lo_opts
use type_phasespacesurface, only: lo_phasespacesurface
use lo_memtracker, only: lo_mem_helper
! 
implicit none
!
type(lo_opts) :: opts
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_crystalstructure) :: uc
class(lo_qpoint_mesh), allocatable :: qp
type(lo_mpi_helper) :: mw
!
type(lo_phasespacesurface) :: ps
type(lo_qpoint) :: qpoint
real(flyt) :: t0
type(lo_mem_helper) :: mem

call mw%init()
t0=walltime()
call opts%parse()
call mem%init()


! Read unitcell
write(*,*) '... read unitcell poscar'
call uc%readfromfile('infile.ucposcar')
call uc%classify('wedge', timereversal=.true.)
if ( opts%readiso ) then
    write(*,*) '... reading isotope distribution from file'
    call uc%readisotopefromfile()
endif

! Read the force constant
call fc%readfromfile(uc,'infile.forceconstant',mem,opts%verbosity)
write(*,*) '... read second order forceconstant'
!if ( opts%thirdorder ) then
if ( opts%intensities ) then
    call fct%readfromfile(uc,'infile.forceconstant_thirdorder')
    write(*,*) '... read third order forceconstant'
endif

! Get a q-mesh
if ( opts%readqmesh ) then
    call lo_read_qmesh_from_file(qp,uc,'infile.qgrid.hdf5',mem,verbosity=opts%verbosity)
    write(*,*) '... read mesh from file'
else
    write(*,*) '... generate mesh'
    call lo_generate_qmesh(qp,uc,opts%qgrid,'wedge',verbosity=opts%verbosity,timereversal=opts%timereversal,headrankonly=.true.,mw=mw,mem=mem)
endif

! Get the actual q-point we are interested in
qpoint%n_invariant_operation=0
if ( trim(opts%highsymmetrypoint) .ne. 'none' ) then
    qpoint%r=uc%coordinate_from_high_symmetry_point_label(opts%highsymmetrypoint)
else
    qpoint%r=uc%fractional_to_cartesian(opts%qpoint,reciprocal=.true.)
endif
qpoint%r=qpoint%r-uc%bz%gshift(qpoint%r)

! Build the surfaces
select type(qp)
type is(lo_wedge_mesh)
    if ( opts%povray ) then
        call ps%generate(qp,uc,fc,fct,qpoint,opts%verbosity,opts%modespec,&
                         calcintens=opts%intensities,mem=mem)
    else
        call ps%generate(qp,uc,fc,fct,qpoint,opts%verbosity,calcintens=opts%intensities,mem=mem)
    endif
class default
    write(*,*) 'Unknown mesh type'
    stop
end select


! Also to a pov-ray file, if necessary!
if ( opts%povray ) then
    call ps%write_to_povray(uc,'outfile.phasespacesurface',opts%thetaphi(1),opts%thetaphi(2),opts%modespec,opts%povrayquality)
else
    ! Dump it to file
cnt: block
    integer :: b1,b2,b3

    write(*,*) 'plus'
    do b1=1,uc%na*3
    do b2=1,uc%na*3
    do b3=1,uc%na*3
        write(*,*) b1,b2,b3,ps%plus(b1,b2,b3)%ntri
    enddo
    enddo
    enddo
    write(*,*) 'minus'
    do b1=1,uc%na*3
    do b2=1,uc%na*3
    do b3=1,uc%na*3
        write(*,*) b1,b2,b3,ps%minus(b1,b2,b3)%ntri
    enddo
    enddo
    enddo
end block cnt    
    !call ps%write_to_hdf5(uc,'outfile.phasespacesurface.hdf5')
endif


write(*,*) 'All done in ',tochar(walltime()-t0),'s'
call mpi_finalize(lo_status)

end program
