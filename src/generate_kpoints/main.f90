program generate_kpoints
!!{!src/generate_kpoints/manual.md!}
use konstanter, only: r8, lo_iou, lo_sqtol
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use gottochblandat, only: walltime, lo_chop, open_file, tochar, lo_determ, lo_kmesh_density
use hdf5_wrappers, only: lo_hdf5_helper
use dump_data

use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_bandstructure, lo_qpoint_mesh, lo_generate_qmesh, lo_wedge_mesh, lo_monkhorst_pack_mesh, lo_fft_mesh
use options, only: lo_opts
use refine, only: get_histograms, fake_integration, massage_mesh ! massage_mesh,get_volume_histogram,fake_integration
implicit none

type(lo_opts) :: opts
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_hdf5_helper) :: h5

type(lo_crystalstructure) :: p
class(lo_qpoint_mesh), allocatable :: qp
real(r8) :: timer

! Grab options, structure and symmetry
init: block
    timer = walltime()
    call opts%parse()
    call mw%init()
    call h5%initialize()
    if (mw%talk .eqv. .false.) then
        opts%verbosity = -100
    end if
    call mem%init()

    ! Structure
    call p%readfromfile('infile.ucposcar', opts%verbosity)
    call p%classify('wedge', timereversal=.true.)

    call mem%assertzero(__FILE__, __LINE__, mw%comm)
end block init

! Get a mesh
genmesh: block
    integer, parameter :: nhist = 100
    real(r8), dimension(nhist) :: hist_vx, hist_sx
    real(r8), dimension(nhist) :: hist_vy, hist_sy
    real(r8), dimension(3, 3) :: dm0, dm1, I3
    real(r8) :: idealvolume, idealside, nqpre, nqpost
    integer :: i, u

    ! Maybe adjust density?
    if (opts%autodens) then
        i = int(anint(sum(opts%qgrid)/3.0_r8))
        i = max(i, 1)
        opts%qgrid = lo_kmesh_density(p%reciprocal_latticevectors, p%na, i)
    end if

    ! First build the mesh
    select case (opts%meshtype)
    case (1)
        call lo_generate_qmesh(qp, p, &
                               griddensity=opts%qgrid, &
                               meshtype='monkhorst', &
                               timereversal=.true., &
                               headrankonly=.false., &
                               mw=mw, mem=mem, verbosity=opts%verbosity + 1)
    case (2)
        call lo_generate_qmesh(qp, p, &
                               griddensity=opts%qgrid, &
                               meshtype='fft', &
                               timereversal=.true., &
                               headrankonly=.false., &
                               mw=mw, mem=mem, verbosity=opts%verbosity + 1)
    case (3)
        call lo_generate_qmesh(qp, p, &
                               griddensity=opts%qgrid, &
                               meshtype='wedge', &
                               timereversal=.true., &
                               headrankonly=.false., &
                               mw=mw, mem=mem, verbosity=opts%verbosity + 1)
    case (4)
        call lo_generate_qmesh(qp, p, &
                               griddensity=opts%qgrid, &
                               meshtype='commensurate', &
                               timereversal=.true., &
                               headrankonly=.false., &
                               mw=mw, mem=mem, verbosity=opts%verbosity + 1)
    end select

    if (opts%refinemesh .and. opts%meshtype .eq. 3) then
        select type (qp); type is (lo_wedge_mesh)
            I3 = 0.0_r8
            do i = 1, 3
                I3(i, i) = 1.0_r8
            end do
            ! Before I refine, decide on what a sensible
            ! range for the volume histograms are.
            idealvolume = (1.0_r8/p%volume)/qp%n_full_tet                     ! ideal in the original division, that is.
            idealside = (idealvolume*6.0_r8*sqrt(2.0_r8))**(1.0_r8/3.0_r8)    ! and the perfect side from that
            ! Get the pre-refined histograms
            call get_histograms(qp, hist_vx, hist_vy, hist_sx, hist_sy, idealvolume, idealside, mw)
            call fake_integration(p, qp, dm0)
            nqpre = qp%n_full_point**(1.0_r8/3.0_r8)
            ! Dump the pre-refined
            if (mw%talk) then
                u = open_file('out', 'outfile.qmesh_histograms_pre_refinement')
                do i = 1, nhist
                    write (u, "(4(1X,E19.12))") hist_vx(i), hist_vy(i), hist_sx(i), hist_sy(i)
                end do
                close (u)
            end if

            ! Now, try to refine the mesh somewhat. Do it in steps.
            call massage_mesh(qp, p, idealside, idealvolume, 1000, opts%forcetolerance*100, opts%sizetolerance + 2.0_r8, .true., mw, mem, opts%verbosity + 1)
            ! call massage_mesh(qp, p, idealside, idealvolume, 1000, opts%forcetolerance*100, opts%sizetolerance + 1_r8, .false., mw, mem, opts%verbosity + 1)
            ! call massage_mesh(qp, p, idealside, idealvolume, 1000, opts%forcetolerance*50, opts%sizetolerance + 0.5_r8, .false., mw, mem, opts%verbosity + 1)
            ! call massage_mesh(qp, p, idealside, idealvolume, 1000, opts%forcetolerance*50, opts%sizetolerance, .false., mw, mem, opts%verbosity + 1)
            ! call massage_mesh(qp, p, idealside, idealvolume, 1000, opts%forcetolerance*10, -1.0_r8, .false., mw, mem, opts%verbosity + 1)
            ! call massage_mesh(qp, p, idealside, idealvolume, 1000, opts%forcetolerance, -1.0_r8, .false., mw, mem, opts%verbosity + 1)
            ! call massage_mesh(qp, p, idealside, idealvolume, 1000, opts%forcetolerance/10, -1.0_r8, .false., mw, mem, opts%verbosity + 1)
            ! do i=1,2
            !     call massage_mesh(qp,p,idealside,idealvolume,1000,opts%forcetolerance,opts%sizetolerance,mw,mem,opts%verbosity+1)
            ! enddo
            ! call massage_mesh(qp,p,idealside,idealvolume,1000,opts%forcetolerance/10,opts%sizetolerance,mw,mem,opts%verbosity+1)

            ! Dump the post-refined
            call get_histograms(qp, hist_vx, hist_vy, hist_sx, hist_sy, idealvolume, idealside, mw)
            call fake_integration(p, qp, dm1)
            nqpost = qp%n_full_point**(1.0_r8/3.0_r8)
            ! Dump the pre-refined
            if (mw%talk) then
                u = open_file('out', 'outfile.qmesh_histograms_post_refinement')
                do i = 1, nhist
                    write (u, "(4(1X,E19.12))") hist_vx(i), hist_vy(i), hist_sx(i), hist_sy(i)
                end do
                close (u)
            end if

            if (mw%talk) then
                write (*, *) '... if the mesh is perfect this should be an identity matrix:'
                write (*, *) '... pre refinement:', norm2(I3 - dm0), nqpre
                do i = 1, 3
                    write (*, *) lo_chop(dm0(:, i), 1E-14_r8)
                end do
                write (*, *) '... post refinement:', norm2(I3 - dm1), nqpost
                do i = 1, 3
                    write (*, *) lo_chop(dm1(:, i), 1E-14_r8)
                end do
            end if
        end select
    else
        if (mw%talk) then
            call fake_integration(p, qp, dm0)
            write (lo_iou, *) '... if the mesh is perfect this should be an identity matrix:'
            do i = 1, 3
                write (lo_iou, *) lo_chop(dm0(:, i), 1E-14_r8)
            end do
        end if
    end if

    ! And finally, write the q-mesh to file
    if (mw%talk) then
        call qp%write_to_file(p, 'outfile.qgrid.hdf5', mem, opts%verbosity + 1)
    end if

    call mem%assertzero(__FILE__, __LINE__, mw%comm)
end block genmesh

! Dump a k-point path
dumppath: block
    type(lo_bandstructure) :: bs
    real(r8), dimension(3) :: v0, v1
    integer :: i, u
    character(len=1000) :: opf

    ! A standard path
    if (opts%readpathfromfile) then
        call bs%read_path_from_file(p, mw, opts%verbosity)
    else
        call bs%standardpath(p, mw, opts%verbosity)
        bs%n_point_per_path = opts%nqpath
    end if

    ! Dump a path file
    u = open_file('out', 'outfile.qpoints_dispersion')
    write (u, *) 'CUSTOM'
    write (u, *) opts%nqpath, '! number of points on each path'
    write (u, *) bs%n_path, '! number of paths'
    opf = "(3(1X,F16.13),2X,3(1X,F16.13),1X,2(1X,A3))"
    do i = 1, bs%n_path
        v0 = lo_chop(matmul(p%inv_reciprocal_latticevectors, bs%segment(i)%r1), lo_sqtol)
        v1 = lo_chop(matmul(p%inv_reciprocal_latticevectors, bs%segment(i)%r2), lo_sqtol)
        write (u, opf) v0, v1, bs%symb_q_start(i), trim(bs%symb_q_end(i))
    end do
    write (u, *) ''
    write (u, '(A)') '# below are all available high-symmetry points. It is possible that all might not have names.'
    do i = 1, p%irrw%nnodes
        v1 = lo_chop(matmul(p%inv_reciprocal_latticevectors, p%irrw%r(:, i)), lo_sqtol)
        write (u, "('#',3(1X,F16.13),2X,A)") v1, p%irrw%label(i)
    end do
    close (u)
end block dumppath

if (mw%talk) then
    write (lo_iou, *) ''
    write (lo_iou, *) 'All done! (', tochar(walltime() - timer), 's)'
end if
call h5%finalize()
call mw%destroy()

end program
