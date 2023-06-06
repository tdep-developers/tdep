module lo_electron_bandstructure_on_path
!!
!! Generate electron energies along paths in reciprocal space
!!
use konstanter, only: r8,lo_iou,lo_huge,lo_hugeint,lo_status,lo_gnuplotterminal,lo_Hartree_to_eV,lo_bohr_to_A
use gottochblandat, only: open_file,tochar,walltime,lo_progressbar_init,lo_progressbar,lo_trueNtimes
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_bandstructure
use lo_electron_dispersion_relations, only: lo_electron_dispersions_kpoint
use type_effective_hamiltonian, only: lo_effective_hamiltonian
use hdf5_wrappers, only: lo_hdf5_helper

implicit none

private
public :: lo_electron_bandstructure

!> Electron bandstructure along a path
type, extends(lo_bandstructure) :: lo_electron_bandstructure
    !> number of bands
    integer :: n_band=-lo_hugeint
    !> number of spin channels
    integer :: n_spin=-lo_hugeint
    !> number of basis functions
    integer :: n_basis=-lo_hugeint
    !> k-points that hold energies and similar things
    type(lo_electron_dispersions_kpoint), dimension(:), allocatable :: k
    contains
        !> write dispersions to file
        procedure :: write_to_file
        ! !> wite phasespace to file
        ! procedure :: write_phasespace_to_file
        ! !> wite phasespace to file
        ! procedure :: write_linewidth_to_file
        !> dump everything to hdf5
        procedure :: write_to_hdf5
        !> generate from effective Hamiltonian
        procedure :: generate
end type

contains

!> Solve the KS eigenproblem on a path
subroutine generate(bs,eh,p,npts,readpath,eigenvectors,mw,mem,verbosity)
    !> electron bandstructure
    class(lo_electron_bandstructure), intent(out) :: bs
    !> aims Hamiltonian
    type(lo_effective_hamiltonian), intent(in) :: eh
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> number of points on path
    integer, intent(in) :: npts
    !> read path from file?
    logical, intent(in) :: readpath
    !> keep the eigenvectors?
    logical, intent(in) :: eigenvectors
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk?
    integer, intent(in) :: verbosity

    real(r8) :: timer

    ! Start timer and tick
    timer=walltime()
    call mem%tick()

    ! Set some basic things
    init: block
        integer :: i

        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'Calculating ks eigenvalues on k-point path (',tochar(mw%n),' ranks)'
        endif
        bs%n_band=eh%n_band
        bs%n_basis=eh%n_basis
        bs%n_spin=eh%n_spin
        bs%n_point_per_path=npts

        if ( readpath ) then
            call bs%read_path_from_file(p,mw,verbosity)
        else
            call bs%standardpath(p,mw,verbosity)
        endif

        ! make some space
        allocate(bs%k(bs%n_point))
        do i=1,bs%n_point
            allocate(bs%k(i)%eigenvalue(bs%n_band,bs%n_spin))
            bs%k(i)%eigenvalue=0.0_r8
            if ( eigenvectors ) then
                allocate(bs%k(i)%eigenvector(bs%n_basis,bs%n_band,bs%n_spin))
                bs%k(i)%eigenvector=0.0_r8
            endif
        enddo
    end block init

    solve: block
        complex(r8), dimension(:,:,:,:), allocatable :: egvbuf
        real(r8), dimension(:,:,:), allocatable :: evbuf
        integer :: i

        call mem%allocate(evbuf,[bs%n_band,bs%n_spin,bs%n_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        evbuf=0.0_r8

        if ( eigenvectors ) then
            call mem%allocate(egvbuf,[bs%n_basis,bs%n_band,bs%n_spin,bs%n_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            egvbuf=0.0_r8
        endif

        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        do i=1,bs%n_point
            if ( mod(i,mw%n) .ne. mw%r ) cycle

            if ( eigenvectors )  then
                call eh%solve_eigenproblem(bs%q(i),evbuf(:,:,i),mem,egvbuf(:,:,:,i))
            else
                call eh%solve_eigenproblem(bs%q(i),evbuf(:,:,i),mem)
            endif

            ! dump progress?
            if ( verbosity .gt. 0 .and. lo_trueNtimes(i,25,bs%n_point) ) then
                call lo_progressbar(' ... solving KS eigenproblem',i,bs%n_point,walltime()-timer)
            endif
        enddo

        ! Collect from all ranks
        call mw%allreduce('sum',evbuf)
        if ( eigenvectors ) call mw%allreduce('sum',egvbuf)

        do i=1,bs%n_point
            bs%k(i)%eigenvalue=evbuf(:,:,i)
            if ( eigenvectors ) then
                bs%k(i)%eigenvector=egvbuf(:,:,:,i)
            endif
        enddo

        call mem%deallocate(evbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( eigenvectors ) call mem%deallocate(egvbuf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) call lo_progressbar(' ... solving KS eigenproblem',bs%n_point,bs%n_point,walltime()-timer)
    end block solve

    call mem%tock(__FILE__,__LINE__,mw%comm)

    if ( verbosity .gt. 0 ) then
        write(*,*) 'Done calculating eigenvalues on path (',tochar(walltime()-timer),'s)'
    endif
end subroutine

! !> Read the electronic band structure from a VASP EIGENVAL file
! subroutine read_from_file(bs,filename,uc,verbosity)
!     class(lo_electron_bandstructure), intent(inout) :: bs
!     type(lo_crystalstructure), intent(inout) :: uc
!     character(len=*), intent(in) :: filename
!     integer, intent(in), optional :: verbosity
!     !
!     integer :: i,j,k,u
!     integer :: nk,nb,verb
!     character(len=1) :: trams
!     real(r8), dimension(3) :: v0
!     real(r8) :: f0
!
!     if ( present(verbosity) ) then
!         verb=verbosity
!     else
!         verb=0
!     endif
!
!     if ( verb .gt. 0 ) then
!         write(*,*) ''
!         write(*,*) 'Reading electron bandstructure from VASP eigenval file'
!     endif
!
!     ! Grab the q-points from file. Mandatory.
!     call bs%readpathfromfile(uc)
!     allocate(bs%k(bs%n_point))
!
!     u=open_file('in',trim(filename))
!         ! skip the header
!         read(u,*) trams
!         read(u,*) trams
!         read(u,*) trams
!         read(u,*) trams
!         read(u,*) trams
!         read(u,*) i,nk,nb
!         bs%n_band=nb
!         bs%n_spin=1
!         ! now it's the real blocks.
!         do i=1,nk
!             read(u,*) v0,f0
!             allocate(bs%k(i)%eigenvalue(nb,1))
!             do j=1,nb
!                 read(u,*) k,bs%k(i)%eigenvalue(j,1)
!             enddo
!         enddo
!     close(u)
! end subroutine

!> Dump the electron bandstructure to hdf5
subroutine write_to_hdf5(bs,filename,efermi,bandgap,mem)
    !> bandstructure
    class(lo_electron_bandstructure), intent(in) :: bs
    !> filename
    character(len=*), intent(in) :: filename
    !> fermi level
    real(r8), intent(in) :: efermi
    !> band gap
    real(r8), intent(in) :: bandgap
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:,:,:), allocatable :: d3
    character(len=1000) :: lblstr
    character(len=10), dimension(:), allocatable :: dumlbl
    integer :: i,j,k

    if ( bs%n_spin .gt. 1 ) then
        write(*,*) 'FIXME SPIN_POLARIZED WRITE HDF5'
        stop
    endif

    ! First dump it in hdf5, create the file
    call h5%init(__FILE__,__LINE__)
    call h5%open_file('write',trim(filename))

    ! Write the x-axis, the ticks for the x-axis, and labels for the ticks
    call h5%store_data(bs%q_axis/lo_bohr_to_A       ,h5%file_id,'k_values',enhet='A^-1')
    call h5%store_data(bs%q_axis_ticks/lo_bohr_to_A ,h5%file_id,'k_ticks' ,enhet='A^-1')

    allocate(dumlbl(size(bs%q_axis_tick_labels,1)))
    do i=1,size(dumlbl,1)
        k=0
        dumlbl(i)='         '
        do j=1,len(dumlbl(i))
            if ( bs%q_axis_tick_labels(i)(j:j) .ne. '|' ) then
                k=k+1
                dumlbl(i)(k:k)=trim(bs%q_axis_tick_labels(i)(j:j))
            endif
        enddo
        ! hdf5 gets angry if I try to write unicode things.
        if ( trim(dumlbl(i)) .eq. 'Γ' ) then
            dumlbl(i)='G'
        endif
    enddo
    lblstr=''
    do i=1,size(dumlbl,1)
        lblstr=trim(lblstr)//' '//trim(dumlbl(i))
    enddo
    call h5%store_attribute(trim(adjustl(lblstr)),h5%file_id,'k_tick_labels')
    deallocate(dumlbl)

    ! Write fermi-level and eventual bandgap
    call h5%store_attribute(efermi*lo_Hartree_to_eV ,h5%file_id,'fermi_level')
    call h5%store_attribute(bandgap*lo_Hartree_to_eV,h5%file_id,'bandgap')

    ! Write electron energies?
    call mem%allocate(d3,[bs%n_band,bs%n_point,bs%n_spin],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    d3=0.0_r8

    if ( allocated(bs%k(1)%eigenvalue) ) then
        do j=1,bs%n_spin
        do i=1,bs%n_point
            d3(:,i,j)=bs%k(i)%eigenvalue(:,j)*lo_Hartree_to_eV
        enddo
        enddo
        call h5%store_data(d3,h5%file_id,'electron_energies',enhet='eV',dimensions='band,kpoint,ispin')
    endif

    if ( allocated(bs%k(1)%electronphononlinewidth) ) then
        do j=1,bs%n_spin
        do i=1,bs%n_point
            d3(:,i,j)=bs%k(i)%electronphononlinewidth(:,j)*lo_Hartree_to_eV
        enddo
        enddo
        call h5%store_data(d3,h5%file_id,'electron_linewidth',enhet='eV')
    endif

    if ( allocated(bs%k(1)%electronphononphasespace) ) then
        do j=1,bs%n_spin
        do i=1,bs%n_point
            d3(:,i,j)=bs%k(i)%electronphononphasespace(:,j)
        enddo
        enddo
        call h5%store_data(d3,h5%file_id,'electron_phasespace')
    endif

    call mem%deallocate(d3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! close the file
    call h5%close_file()
    call h5%destroy(__FILE__,__LINE__)
end subroutine

! !> Dump the electron phasespace to file
! subroutine write_phasespace_to_file(bs,filename)
!     !> bandstructure
!     class(lo_electron_bandstructure), intent(in) :: bs
!     !> filename
!     character(len=*), intent(in) :: filename
!
!     character(len=100) :: opf
!     integer :: i,j,u
!
!     ! Dump the raw data
!     opf="("//tochar(bs%n_band*bs%n_spin+1)//"(1X,E18.12))"
!     u=open_file('out',trim(filename))
!         do i=1,bs%n_point
!             ! Write
!             select case(bs%n_spin)
!             case(1)
!                 write(u,opf) bs%q_axis(i),bs%k(i)%electronphononphasespace(:,1)*lo_Hartree_to_eV
!             case(2)
!                 write(u,opf) bs%q_axis(i),bs%k(i)%electronphononphasespace(:,1)*lo_Hartree_to_eV,&
!                     bs%k(i)%electronphononphasespace(:,2)*lo_Hartree_to_eV
!             end select
!         enddo
!     close(u)
!
!     ! Nice gnuplot file for the dispersion relations
!     u=open_file('out',trim(filename)//'.gnuplot')
!         write(u,*) 'set terminal '//lo_gnuplotterminal//' size 500,350 enhanced font "CMU Serif,10"'
!         write(u,"('set xrange [',F12.5,':',F12.5,']')") minval(bs%q_axis),maxval(bs%q_axis)
!         ! write(u,"('set yrange [',F12.5,':',F12.5,']')") emin,emax
!         ! set the ticks
!         write(u,*) 'unset xtics'
!         write(u,*) 'set xtics ( "'//trim(bs%q_axis_tick_labels(1))//'" 0.0 ) '
!         do i=2,bs%n_path+1
!             write(u,'(A)',advance='no') 'set xtics add ('
!             write(u,'(A)',advance='no') '"'//trim(bs%q_axis_tick_labels(i))//'" '
!             write(u,'(F9.6)',advance='no') real(bs%q_axis_ticks(i))
!             write(u,*) ' )'
!         enddo
!         ! set gridlines at tics
!         write(u,*) 'set grid xtics lc rgb "#888888" lw 1 lt 0'
!         write(u,*) 'set xzeroaxis linewidth 0.1 linecolor 0 linetype 1'
!         write(u,*) 'set ytics scale 0.5'
!         write(u,*) 'set xtics scale 0.5'
!         write(u,*) 'set mytics 10'
!         write(u,*) 'unset key'
!         write(u,*) ' set ylabel "Energy (eV)" '
!         ! Plot
!         write(u,'(A)',advance='no') 'plot'
!         do j=2,bs%n_band+1
!             write(u,'(A)',advance='no') ' "'//trim(filename)//'" u 1:'//tochar(j)//' w line lc rgb "#318712"'
!             if( j .lt. bs%n_band+1 ) then
!                 write(u,'(A)') ',\'
!             endif
!         enddo
!     close(u)
! end subroutine

! !> Dump the electron linewidth to file
! subroutine write_linewidth_to_file(bs,filename)
!     !> bandstructure
!     class(lo_electron_bandstructure), intent(in) :: bs
!     !> filename
!     character(len=*), intent(in) :: filename
!
!     character(len=100) :: opf
!     integer :: i,j,u
!
!     ! Dump the raw data
!     opf="("//tochar(bs%n_band*bs%n_spin+1)//"(1X,E18.12))"
!     u=open_file('out',trim(filename))
!         do i=1,bs%n_point
!             ! Write
!             select case(bs%n_spin)
!             case(1)
!                 write(u,opf) bs%q_axis(i),bs%k(i)%electronphononlinewidth(:,1)*lo_Hartree_to_eV
!             case(2)
!                 write(u,opf) bs%q_axis(i),bs%k(i)%electronphononlinewidth(:,1)*lo_Hartree_to_eV,&
!                     bs%k(i)%electronphononlinewidth(:,2)*lo_Hartree_to_eV
!             end select
!         enddo
!     close(u)
!
!     ! Nice gnuplot file for the dispersion relations
!     u=open_file('out',trim(filename)//'.gnuplot')
!         write(u,*) 'set terminal '//lo_gnuplotterminal//' size 500,350 enhanced font "CMU Serif,10"'
!         write(u,"('set xrange [',F12.5,':',F12.5,']')") minval(bs%q_axis),maxval(bs%q_axis)
!         ! write(u,"('set yrange [',F12.5,':',F12.5,']')") emin,emax
!         ! set the ticks
!         write(u,*) 'unset xtics'
!         write(u,*) 'set xtics ( "'//trim(bs%q_axis_tick_labels(1))//'" 0.0 ) '
!         do i=2,bs%n_path+1
!             write(u,'(A)',advance='no') 'set xtics add ('
!             write(u,'(A)',advance='no') '"'//trim(bs%q_axis_tick_labels(i))//'" '
!             write(u,'(F9.6)',advance='no') real(bs%q_axis_ticks(i))
!             write(u,*) ' )'
!         enddo
!         ! set gridlines at tics
!         write(u,*) 'set grid xtics lc rgb "#888888" lw 1 lt 0'
!         write(u,*) 'set xzeroaxis linewidth 0.1 linecolor 0 linetype 1'
!         write(u,*) 'set ytics scale 0.5'
!         write(u,*) 'set xtics scale 0.5'
!         write(u,*) 'set mytics 10'
!         write(u,*) 'unset key'
!         write(u,*) ' set ylabel "Energy (eV)" '
!         ! Plot
!         write(u,'(A)',advance='no') 'plot'
!         do j=2,bs%n_band+1
!             write(u,'(A)',advance='no') ' "'//trim(filename)//'" u 1:'//tochar(j)//' w line lc rgb "#318712"'
!             if( j .lt. bs%n_band+1 ) then
!                 write(u,'(A)') ',\'
!             endif
!         enddo
!     close(u)
! end subroutine

!> Dump the electron bandstructure to file
subroutine write_to_file(bs,filename,efermi,emin,emax)
    !> bandstructure
    class(lo_electron_bandstructure), intent(in) :: bs
    !> filename
    character(len=*), intent(in) :: filename
    !> some energies for neater plots
    real(r8), intent(in) :: efermi,emin,emax

    character(len=100) :: opf
    integer :: i,j,u

    ! Dump the raw data
    opf="("//tochar(bs%n_band*bs%n_spin+1)//"(1X,E18.12))"
    u=open_file('out',trim(filename))
        do i=1,bs%n_point
            ! Write
            select case(bs%n_spin)
            case(1)
                write(u,opf) bs%q_axis(i),(bs%k(i)%eigenvalue(:,1)-efermi)*lo_Hartree_to_eV
            case(2)
                write(u,opf) bs%q_axis(i),(bs%k(i)%eigenvalue(:,1)-efermi)*lo_Hartree_to_eV,&
                    (bs%k(i)%eigenvalue(:,2)-efermi)*lo_Hartree_to_eV
            end select
        enddo
    close(u)

    ! Nice gnuplot file for the dispersion relations
    u=open_file('out',trim(filename)//'.gnuplot')
        write(u,*) 'set terminal '//lo_gnuplotterminal//' size 500,350 enhanced font "CMU Serif,10"'
        write(u,"('set xrange [',F12.5,':',F12.5,']')") minval(bs%q_axis),maxval(bs%q_axis)
        write(u,"('set yrange [',F12.5,':',F12.5,']')") emin,emax
        ! set the ticks
        write(u,*) 'unset xtics'
        write(u,*) 'set xtics ( "'//trim(bs%q_axis_tick_labels(1))//'" 0.0 ) '
        do i=2,bs%n_path+1
            write(u,'(A)',advance='no') 'set xtics add ('
            write(u,'(A)',advance='no') '"'//trim(bs%q_axis_tick_labels(i))//'" '
            write(u,'(F9.6)',advance='no') real(bs%q_axis_ticks(i))
            write(u,*) ' )'
        enddo
        ! set gridlines at tics
        write(u,*) 'set grid xtics lc rgb "#888888" lw 1 lt 0'
        write(u,*) 'set xzeroaxis linewidth 0.1 linecolor 0 linetype 1'
        write(u,*) 'set ytics scale 0.5'
        write(u,*) 'set xtics scale 0.5'
        write(u,*) 'set mytics 10'
        write(u,*) 'unset key'
        write(u,*) ' set ylabel "Energy (eV)" '
        ! Plot
        write(u,'(A)',advance='no') 'plot'
        do j=2,bs%n_band+1
            write(u,'(A)',advance='no') ' "outfile.electron_bandstructure" u 1:'//tochar(j)//' w line lc rgb "#318712"'
            if( j .lt. bs%n_band+1 ) then
                write(u,'(A)') ',\'
            endif
        enddo
    close(u)
end subroutine

!subroutine write_intensity(bs,efermi,filename,logscale)
!    class(lo_electron_bandstructure), intent(inout) :: bs
!    real(r8), intent(in) :: efermi
!    character(len=*), intent(in), optional :: filename
!    logical, intent(in), optional :: logscale
!    !
!    character(len=500) :: fn,gpfn,hdffn,lblstr
!    character(len=10), dimension(:), allocatable :: dumlbl
!    integer :: i,j,k
!    logical :: ls
!    integer(HID_T) :: file_id
!
!    ! Another filename perhaps?
!    if ( present(filename) ) then
!        fn=trim(filename)
!        gpfn=trim(filename)//'.gnuplot'
!        hdffn=trim(filename)//'.hdf5'
!    else
!        fn='outfile.spectralfunction'
!        gpfn='outfile.spectralfunction.gnuplot'
!        hdffn='outfile.spectralfunction.hdf5'
!    endif
!
!    ! Logscale?
!    if ( present(logscale) ) then
!        ls=logscale
!    else
!        ls=.true.
!    endif
!    !
!    ! Write files
!    !
!!    ne=size(bs%intensity,2)
!!    u=open_file('out',trim(fn))
!!        do i=1,bs%n_point
!!        do j=1,ne
!!            write(u,*) bs%q_axis(i),bs%energyaxis(j)-efermi,bs%intensity(i,j)
!!        enddo
!!        enddo
!!    close(u)
!!    !
!!    u=open_file('out',trim(gpfn))
!!        write(u,*) 'set terminal '//lo_gnuplotterminal//' size 700,500 enhanced font "CMU Serif,8"'
!!        write(u,*) 'set border lw 0.5'
!!        write(u,*) 'set pm3d map'
!!        write(u,"('set xrange [',F12.5,':',F12.5,']')") minval(bs%q_axis),maxval(bs%q_axis)
!!        write(u,"('set yrange [',F12.5,':',F12.5,']')") emin,emax
!!        !
!!        ! set the ticks
!!        !
!!        write(u,*) 'unset xtics'
!!        write(u,*) 'set xtics ( "'//trim(bs%q_axis_tick_labels(1))//'" 0.0 ) '
!!        do i=2,bs%npath+1
!!            write(u,'(A)',advance='no') 'set xtics add ('
!!            write(u,'(A)',advance='no') '"'//trim(bs%q_axis_tick_labels(i))//'" '
!!            write(u,'(F9.6)',advance='no') (i-1)*maxval(bs%q_axis)/bs%npath
!!            write(u,*) ' )'
!!        enddo
!!        write(u,*) 'set tics out'
!!        write(u,*) 'set ytics scale 0.5'
!!        write(u,*) 'set xtics scale 0.5'
!!        write(u,*) 'unset key'
!!        write(u,*) 'unset colorbox'
!!        write(u,*) 'set ylabel "Energy (eV)"'
!!        call lo_dump_palette_to_gnuplot(u,logscale=ls)
!!        write(u,*) 'plot "'//trim(fn)//'" u 1:2:3 w image'
!!    close(u)
!    !
!    ! Write it as HDF5
!    !
!    call h5open_f(lo_status)
!    call h5fcreate_f(trim(hdffn), H5F_ACC_TRUNC_F, file_id, lo_status)
!write(*,*) 'opened hdf5'
!    ! store the axes and intensity
!    call lo_h5_store_data(bs%q_axis,file_id,'q_values')
!write(*,*) 'wrote q-axis'
!    call lo_h5_store_data(bs%energyaxis-efermi,file_id,'energy_values')
!write(*,*) 'wrote energy-axis'
!    call lo_h5_store_data(bs%intensity,file_id,'intensity')
!write(*,*) 'wrote intensity'
!    ! the x-axis ticks
!    call lo_h5_store_data(bs%q_axis_ticks,file_id,'q_ticks')
!    ! and the tick labels, horribly annoying.
!    allocate(dumlbl(size(bs%q_axis_tick_labels,1)))
!    do i=1,size(dumlbl,1)
!        k=0
!        dumlbl(i)='         '
!        do j=1,len(dumlbl(i))
!            if ( bs%q_axis_tick_labels(i)(j:j) .ne. '|' ) then
!                k=k+1
!                dumlbl(i)(k:k)=bs%q_axis_tick_labels(i)(j:j)
!            endif
!        enddo
!        if ( trim(dumlbl(i)) .eq. 'Γ' ) then
!            dumlbl(i)='G'
!        endif
!    enddo
!    lblstr=''
!    do i=1,size(dumlbl,1)
!        lblstr=trim(lblstr)//' '//trim(dumlbl(i))
!    enddo
!    ! and some extra data
!    call lo_h5_store_attribute(trim(adjustl(lblstr)),file_id,'q_tick_labels')
!    call lo_h5_store_attribute(bs%ne,file_id,'number_of_energies')
!    call lo_h5_store_attribute(bs%n_point,file_id,'number_of_kpoints')
!    call lo_h5_store_attribute('eV',file_id,'energy_unit')
!    call h5fclose_f(file_id, lo_status)
!    call h5close_f(lo_status)
!    !
!end subroutine

end module
