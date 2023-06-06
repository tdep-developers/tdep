#include "precompilerdefinitions"
program crystal_structure_info
!!{!src/crystal_structure_info/manual.md!}
use konstanter, only: flyt, lo_tol, lo_sqtol, lo_bohr_to_A, lo_gitbranch, lo_gitrevision
use gottochblandat, only: lo_frobnorm, lo_chop, open_file, tochar, lo_clean_fractional_coordinates, lo_determ, &
                          lo_fetch_tolerance, walltime, lo_return_unique, lo_invert3x3matrix, qsort
use mpi_wrappers, only: lo_mpi_helper
use options, only: lo_opts
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_bandstructure
use type_symmetryoperation, only: lo_spacegroup_operation, lo_symset
use type_voronoi, only: lo_voronoi_diagram
implicit none

type(lo_opts) :: opts
type(lo_mpi_helper) :: mw
type(lo_crystalstructure) :: p
real(flyt) :: timer

! Read and classify structure
init: block
    ! Start timer
    timer = walltime()
    ! Get the command line arguments
    call opts%parse()
    ! Init MPI
    call mw%init()
    ! Sanity test, don't run this in parallel.
    if (mw%n .gt. 1) then
        write (*, *) 'Do not run this in parallel'
        call mw%destroy()
        stop
    end if

    ! Read unitcell
    call p%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    !call p%classify('irrep', timereversal=opts%timereversal)
    call p%classify('wedge', timereversal=opts%timereversal)
end block init

! Report on structure
report1: block
    real(flyt) :: f0
    integer :: i

    write (*, *) '  git branch: ', trim(lo_gitbranch)
    write (*, *) 'git revision: ', trim(lo_gitrevision)

    write (*, *) 'Info about the crystal structure:'
    write (*, *) '... I believe the basis is ', trim(p%info%bravaislatticelongname), ' (', trim(p%info%bravaislattice), ') with the basis vectors'
    write (*, "(1X,'a1:',3(2X,F18.12))") p%latticevectors(:, 1)*lo_bohr_to_A
    write (*, "(1X,'a2:',3(2X,F18.12))") p%latticevectors(:, 2)*lo_bohr_to_A
    write (*, "(1X,'a3:',3(2X,F18.12))") p%latticevectors(:, 3)*lo_bohr_to_A
    write (*, *) '... and the atoms at these positions (in fractional coordinates)'
    do i = 1, p%na
        write (*, '(1X,A8,3(2x,F15.13))') trim(p%atomic_symbol(p%species(i))), p%r(:, i)
    end do

    ! Print what the conventional basis is
    f0 = lo_frobnorm(p%latticevectors - p%info%unique_conventional_basis)
    if (f0 .gt. lo_tol) then
        write (*, *) 'I believe this is the conventional lattice:'
        write (*, "(1X,'a1:',3(2X,F18.12))") p%info%unique_conventional_basis(:, 1)*lo_bohr_to_A
        write (*, "(1X,'a2:',3(2X,F18.12))") p%info%unique_conventional_basis(:, 2)*lo_bohr_to_A
        write (*, "(1X,'a3:',3(2X,F18.12))") p%info%unique_conventional_basis(:, 3)*lo_bohr_to_A
    end if
end block report1

! Report on high-symmetry points of the BZ
report2: block
    type(lo_bandstructure) :: bs
    real(flyt), dimension(3) :: v0, v1
    integer :: i, u
    character(len=1000) :: opf

    write (*, *) ' '
    write (*, *) 'Available high symmetry points:'
    write (*, *) '(the ones called NP are in fact high symmetry points, but noone has bothered to name them)'
    write (*, *) '      Cartesian coordinates                        Fractional coordinates'
    do i = 1, p%irrw%nnodes
        v0 = lo_chop(p%irrw%r(:, i), lo_tol)
        v1 = lo_chop(matmul(p%inv_reciprocal_latticevectors, v0), lo_tol)
        write (*, "(1X,I2,1X,A4,3(1X,F13.10),3X,3(1X,F13.10))") i, p%irrw%label(i), v0/lo_bohr_to_A, v1
    end do

    call p%write_bz_to_hdf5('outfile.brillouin_zone.hdf5')
    call bs%standardpath(p, mw, verbosity=0)
    bs%n_point_per_path = 100
    write (*, *) ' '
    write (*, *) 'Per default, the following path will be used:'
    do i = 1, bs%n_path
        write (*, *) trim(bs%symb_q_start(i)), ' -> ', trim(bs%symb_q_end(i))
    end do
    write (*, *) 'It is written in "outfile.qpoints_dispersion", modify and'
    write (*, *) 'copy it to "infile.qpoints_dispersion" if you want'

    ! Dump a path file
    u = open_file('out', 'outfile.qpoints_dispersion')
    opf = "(A5,22X,' ! Bravais lattice type')"
    write (u, opf) p%info%bravaislattice
    opf = "(I5,22X,' ! Number of points on each path')"
    write (u, opf) 100
    opf = "(I5,22X,' ! Number paths between special points')"
    write (u, opf) bs%n_path
    opf = "(A3,1X,A3,20X,' ! Starting and ending special point')"
    write (u, opf) bs%symb_q_start(1), bs%symb_q_end(1)
    opf = "(A3,1X,A3,20X,' !')"
    do i = 2, bs%n_path
        write (u, opf) bs%symb_q_start(i), bs%symb_q_end(i)
    end do
    close (u)

end block report2

if (opts%printsymmetry) then
    report3: block
        type(lo_symset) :: sym
        integer, dimension(:), allocatable :: di
        integer :: i, j, u
        character(len=1000) :: opf, str

        write (*, *) ' Symmetry operations:'
        do i = 1, p%sym%n
            write (*, *) 'Operation ', tochar(i), ' ', p%sym%op(i)%whatkind
            write (*, *) '        Axis:', lo_chop(p%sym%op(i)%axis, lo_sqtol)
            if (p%sym%op(i)%fold .gt. 0) write (*, *) '       Angle: ', tochar(int(anint(360.0_flyt/p%sym%op(i)%fold)))
            opf = '        FMAP:'
            do j = 1, size(p%sym%op(i)%fmap)
                opf = trim(opf)//' '//tochar(j)//'->'//tochar(p%sym%op(i)%fmap(j))
            end do
            write (*, *) trim(opf)
            !write(*,*) '        FMAP:',p%sym%op(i)%fmap
            write (*, '(1X,A,3(1X,F9.6))') ' translation:', p%sym%op(i)%ftr
            write (*, *) ' ... in Cartesian coordinates'
            do j = 1, 3
                write (*, "(6(3X,F15.11))") p%sym%op(i)%m(j, :)
                !write(*,"(3(3X,A10))") lo_fancy_real2char(p%sym%op(i)%m(j,1)),lo_fancy_real2char(p%sym%op(i)%m(j,2)),lo_fancy_real2char(p%sym%op(i)%m(j,3))
            end do
            write (*, *) ' ... in fractional coordinates'
            do j = 1, 3
                write (*, "(6(3X,F15.11))") p%sym%op(i)%fm(j, :)
                !write(*,"(6(3X,A3))") lo_fancy_real2char(p%sym%op(i)%fm(j,1)),lo_fancy_real2char(p%sym%op(i)%fm(j,2)),lo_fancy_real2char(p%sym%op(i)%fm(j,3))
            end do
            write (*, *) ' '
        end do

        ! Dump symmetry operations to Mathematica-format, in case you need to play with it
        u = open_file('out', 'outfile.symops_mathematica')
        write (u, *) 'ops=ConstantArray[0,{'//tochar(p%sym%n)//',3,3}];'
        do i = 1, p%sym%n
            do j = 1, 3
                write (u, *) 'ops[['//tochar(i)//','//tochar(j)// &
                    ']]={ '//lo_fancy_real2char(p%sym%op(i)%m(1, j))//','//lo_fancy_real2char(p%sym%op(i)%m(2, j))//','//lo_fancy_real2char(p%sym%op(i)%m(3, j))//'};'
            end do
        end do

        write (u, *) 'characters=ConstantArray[0,{'//tochar(p%sym%n_irrep)//','//tochar(p%sym%n)//'}];'
        do i = 1, p%sym%n_irrep
        do j = 1, p%sym%n
            write (u, *) 'characters[['//tochar(i)//','//tochar(j)// &
                ']]='//tochar(int(p%sym%ir(i)%character_element(j)))//';'
        end do
        end do

        write (u, *) 'fmap=ConstantArray[0,{'//tochar(p%sym%n)//','//tochar(p%na)//'}];'
        do i = 1, p%sym%n
            do j = 1, p%na
                write (u, *) 'fmap[['//tochar(i)//','//tochar(j)//']]='//tochar(p%sym%op(i)%fmap(j))//';'
            end do
        end do
        close (u)

        ! write (*, *) ' Character table:'

        ! ! ! Fancy output
        ! allocate (di(p%sym%n_class))
        ! di = 0
        ! str = "+-"
        ! do i = 1, p%sym%n_class
        !     di(i) = i
        !     str = trim(adjustl(str))//'----'
        ! end do
        ! write (*, *) ''
        ! write (*, "(6X,"//tochar(p%sym%n_class)//"(2X,I2))") di
        ! write (*, "(5X,A)") trim(str)
        ! do i = 1, p%sym%n_class
        !     write (*, "(2X,I2,' |',"//tochar(p%sym%n_class)//"(2X,I2))") i, int(p%sym%ir(i)%character_class) !int(dm0(:,i))
        ! end do
        ! deallocate (di)
        ! write (*, *) ''

        ! ! New symmetry set
        ! call sym%generate(p%latticevectors, .true., p%r, p%species, verbosity=1, tolerance=1E-6_flyt, symmorphic=.true.)
        !call sym%get_character_table(1)

        ! allocate(di(sym%n_class))
        ! di=0
        ! str="+-"
        ! do i=1,sym%n_class
        !     di(i)=i
        !     str=trim(adjustl(str))//'----'
        ! enddo
        ! write(*,*) ''
        ! write(*,"(6X,"//tochar(sym%n_class)//"(2X,I2))") di
        ! write(*,"(5X,A)") trim(str)
        ! do i=1,sym%n_class
        !     write(*,"(2X,I2,' |',"//tochar(sym%n_class)//"(2X,I2))") i,int(sym%ir(i)%character_class) !int(dm0(:,i))
        ! enddo
        ! deallocate(di)
        ! write(*,*) ''

    end block report3
end if

!! Dump the character table
!report4: block
!    character(len=1000) :: str
!    integer, dimension(p%sym%n_class) :: di
!    integer :: i
!
!    str="+-"
!    do i=1,p%sym%n_class
!        di(i)=i
!        str=trim(adjustl(str))//'----'
!    enddo
!    write(*,*) ''
!    write(*,*) 'Character table:'
!    write(*,"(6X,"//tochar(p%sym%n_class)//"(2X,I2))") di
!    write(*,"(5X,A)") trim(str)
!    do i=1,p%sym%n_irrep
!        write(*,"(2X,I2,' |',"//tochar(p%sym%n_class)//"(2X,I2))") i,int(p%sym%ir(i)%character_class)
!    enddo
!    write(*,*) ''
!end block report4

write (*, *) 'All done!'

contains

! function distance_point_box(point,r0,box) result(d)
!     real(flyt), dimension(3), intent(in) :: point,r0
!     real(flyt), dimension(3,3), intent(in) :: box
!     real(flyt) :: d
!
!     real(flyt), dimension(3,8) :: corners
!     real(flyt), dimension(3) :: va,vb,vc
!     ! So, we have a weird box thing. That means
!     ! 8 corners
!     ! 12 edges
!     ! 6 faces
!     corners(:,1)=r0
!     corners(:,2)=r0 + box(:,1)
!     corners(:,3)=r0 + box(:,2)
!     corners(:,4)=r0 + box(:,3)
!     corners(:,5)=r0 + box(:,1)+box(:,2)
!     corners(:,6)=r0 + box(:,2)+box(:,3)
!     corners(:,7)=r0 + box(:,1)+box(:,3)
!     corners(:,8)=r0 + box(:,1)+box(:,2)+box(:,3)
!
!
!     ! First face
!     va=box(:,1)
!     vb=box(:,2)
!
!
!
!
!
! end function

!subroutine heapsort(a)
!   real(flyt), intent(inout) :: a(0:)
!   integer :: start, n, bottom
!   real(flyt) :: temp
!
!   n = size(a)
!   do start = (n - 2) / 2, 0, -1
!     call siftdown(a, start, n);
!   end do
!
!   do bottom = n - 1, 1, -1
!     temp = a(0)
!     a(0) = a(bottom)
!     a(bottom) = temp;
!     call siftdown(a, 0, bottom)
!   end do
!end subroutine heapsort
!
!subroutine siftdown(a, start, bottom)
!  real(flyt), intent(in out) :: a(0:)
!  integer, intent(in) :: start, bottom
!  integer :: child, root
!  real(flyt) :: temp
!
!  root = start
!  do while(root*2 + 1 < bottom)
!    child = root * 2 + 1
!
!    if (child + 1 < bottom) then
!      if (a(child) < a(child+1)) child = child + 1
!    end if
!
!    if (a(root) < a(child)) then
!      temp = a(child)
!      a(child) = a (root)
!      a(root) = temp
!      root = child
!    else
!      return
!    end if
!  end do
!end subroutine siftdown

!> converts a double to a character, and tries to find well defined numbers
function lo_fancy_real2char(x, tolerance) result(s)
    !> number to be converted
    real(flyt), intent(in) :: x
    !> optional tolerance
    real(flyt), intent(in), optional :: tolerance
    !> the string
    character(len=:), allocatable :: s

    real(flyt) :: tol
    real(flyt), dimension(67), parameter :: welldefined = [ &
                                            0.0000000000000000_flyt, &
                                            0.1000000000000000_flyt, &
                                            -0.1000000000000000_flyt, &
                                            0.1111111111111111_flyt, &
                                            -0.1111111111111111_flyt, &
                                            0.1250000000000000_flyt, &
                                            -0.1250000000000000_flyt, &
                                            0.1428571428571428_flyt, &
                                            -0.1428571428571428_flyt, &
                                            0.1666666666666667_flyt, &
                                            -0.1666666666666667_flyt, &
                                            0.2000000000000000_flyt, &
                                            -0.2000000000000000_flyt, &
                                            0.2222222222222222_flyt, &
                                            -0.2222222222222222_flyt, &
                                            0.2500000000000000_flyt, &
                                            -0.2500000000000000_flyt, &
                                            0.2857142857142857_flyt, &
                                            -0.2857142857142857_flyt, &
                                            0.3000000000000000_flyt, &
                                            -0.3000000000000000_flyt, &
                                            0.3333333333333333_flyt, &
                                            -0.3333333333333333_flyt, &
                                            0.3750000000000000_flyt, &
                                            -0.3750000000000000_flyt, &
                                            0.4000000000000000_flyt, &
                                            -0.4000000000000000_flyt, &
                                            0.4285714285714285_flyt, &
                                            -0.4285714285714285_flyt, &
                                            0.4444444444444444_flyt, &
                                            -0.4444444444444444_flyt, &
                                            0.5000000000000000_flyt, &
                                            -0.5000000000000000_flyt, &
                                            0.5555555555555556_flyt, &
                                            -0.5555555555555556_flyt, &
                                            0.5714285714285714_flyt, &
                                            -0.5714285714285714_flyt, &
                                            0.6000000000000000_flyt, &
                                            -0.6000000000000000_flyt, &
                                            0.6250000000000000_flyt, &
                                            -0.6250000000000000_flyt, &
                                            0.6666666666666667_flyt, &
                                            -0.6666666666666667_flyt, &
                                            0.7000000000000000_flyt, &
                                            -0.7000000000000000_flyt, &
                                            0.7142857142857143_flyt, &
                                            -0.7142857142857143_flyt, &
                                            0.7500000000000000_flyt, &
                                            -0.7500000000000000_flyt, &
                                            0.7777777777777778_flyt, &
                                            -0.7777777777777778_flyt, &
                                            0.8000000000000000_flyt, &
                                            -0.8000000000000000_flyt, &
                                            0.8333333333333334_flyt, &
                                            -0.8333333333333334_flyt, &
                                            0.8571428571428571_flyt, &
                                            -0.8571428571428571_flyt, &
                                            0.8660254037844386_flyt, &
                                            -0.8660254037844386_flyt, &
                                            0.8750000000000000_flyt, &
                                            -0.8750000000000000_flyt, &
                                            0.8888888888888888_flyt, &
                                            -0.8888888888888888_flyt, &
                                            0.9000000000000000_flyt, &
                                            -0.9000000000000000_flyt, &
                                            1.0000000000000000_flyt, &
                                            -1.0000000000000000_flyt]
    character(len=12), dimension(67), parameter :: welldefined_string = [ &
                                                   '0          ', & !    0.0000000000000000_flyt,&
                                                   '1/10       ', & !    0.1000000000000000_flyt,&
                                                   '-1/10      ', & !   -0.1000000000000000_flyt,&
                                                   'x          ', & !    0.1111111111111111_flyt,&
                                                   'x          ', & !   -0.1111111111111111_flyt,&
                                                   '1/8        ', & !    0.1250000000000000_flyt,&
                                                   '/1/8       ', & !   -0.1250000000000000_flyt,&
                                                   'x          ', & !    0.1428571428571428_flyt,&
                                                   'x          ', & !   -0.1428571428571428_flyt,&
                                                   '1/6        ', & !    0.1666666666666667_flyt,&
                                                   '-1/6       ', & !   -0.1666666666666667_flyt,&
                                                   '1/5        ', & !    0.2000000000000000_flyt,&
                                                   '-1/5       ', & !   -0.2000000000000000_flyt,&
                                                   'x          ', & !    0.2222222222222222_flyt,&
                                                   'x          ', & !   -0.2222222222222222_flyt,&
                                                   '1/4        ', & !    0.2500000000000000_flyt,&
                                                   '-1/4       ', & !   -0.2500000000000000_flyt,&
                                                   'x          ', & !    0.2857142857142857_flyt,&
                                                   'x          ', & !   -0.2857142857142857_flyt,&
                                                   'x          ', & !    0.3000000000000000_flyt,&
                                                   'x          ', & !   -0.3000000000000000_flyt,&
                                                   'x          ', & !    0.3333333333333333_flyt,&
                                                   'x          ', & !   -0.3333333333333333_flyt,&
                                                   'x          ', & !    0.3750000000000000_flyt,&
                                                   'x          ', & !   -0.3750000000000000_flyt,&
                                                   'x          ', & !    0.4000000000000000_flyt,&
                                                   'x          ', & !   -0.4000000000000000_flyt,&
                                                   'x          ', & !    0.4285714285714285_flyt,&
                                                   'x          ', & !   -0.4285714285714285_flyt,&
                                                   'x          ', & !    0.4444444444444444_flyt,&
                                                   'x          ', & !   -0.4444444444444444_flyt,&
                                                   '1/2        ', & !    0.5000000000000000_flyt,&
                                                   '-1/2       ', & !   -0.5000000000000000_flyt,&
                                                   'x          ', & !    0.5555555555555556_flyt,&
                                                   'x          ', & !   -0.5555555555555556_flyt,&
                                                   'x          ', & !    0.5714285714285714_flyt,&
                                                   'x          ', & !   -0.5714285714285714_flyt,&
                                                   'x          ', & !    0.6000000000000000_flyt,&
                                                   'x          ', & !   -0.6000000000000000_flyt,&
                                                   'x          ', & !    0.6250000000000000_flyt,&
                                                   'x          ', & !   -0.6250000000000000_flyt,&
                                                   'x          ', & !    0.6666666666666667_flyt,&
                                                   'x          ', & !   -0.6666666666666667_flyt,&
                                                   'x          ', & !    0.7000000000000000_flyt,&
                                                   'x          ', & !   -0.7000000000000000_flyt,&
                                                   'x          ', & !    0.7142857142857143_flyt,&
                                                   'x          ', & !   -0.7142857142857143_flyt,&
                                                   'x          ', & !    0.7500000000000000_flyt,&
                                                   'x          ', & !   -0.7500000000000000_flyt,&
                                                   'x          ', & !    0.7777777777777778_flyt,&
                                                   'x          ', & !   -0.7777777777777778_flyt,&
                                                   'x          ', & !    0.8000000000000000_flyt,&
                                                   'x          ', & !   -0.8000000000000000_flyt,&
                                                   'x          ', & !    0.8333333333333334_flyt,&
                                                   'x          ', & !   -0.8333333333333334_flyt,&
                                                   'x          ', & !    0.8571428571428571_flyt,&
                                                   'x          ', & !   -0.8571428571428571_flyt,&
                                                   'Sqrt[3]/2  ', & !    0.8660254037844386_flyt,&
                                                   '-Sqrt[3]/2 ', & !   -0.8660254037844386_flyt,&
                                                   'x          ', & !    0.8750000000000000_flyt,&
                                                   'x          ', & !   -0.8750000000000000_flyt,&
                                                   'x          ', & !    0.8888888888888888_flyt,&
                                                   'x          ', & !   -0.8888888888888888_flyt,&
                                                   '9/10       ', & !    0.9000000000000000_flyt,&
                                                   '-9/10      ', & !   -0.9000000000000000_flyt,&
                                                   '1          ', & !    1.0000000000000000_flyt,&
                                                   '-1         '] !   -1.0000000000000000_flyt]
    integer :: i

    if (present(tolerance)) then
        tol = tolerance
    else
        tol = lo_sqtol
    end if

    do i = 1, 67
        if (abs(x - welldefined(i)) .lt. tol) then
            s = trim(welldefined_string(i))
            return
        end if
    end do
    s = 'dunno'

end function

end program
