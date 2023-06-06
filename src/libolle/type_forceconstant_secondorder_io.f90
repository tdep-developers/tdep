submodule(type_forceconstant_secondorder) type_forceconstant_secondorder_io
use konstanter, only: lo_forceconstant_2nd_HartreeBohr_to_eVA, lo_forceconstant_2nd_eVA_to_HartreeBohr, lo_emu_to_amu, &
                      lo_exitcode_io, lo_bohr_to_A
use gottochblandat, only: open_file
use type_qpointmesh, only: lo_qpoint_mesh,lo_generate_qmesh
use type_voronoi_distancetable, only: lo_voronoi_distancetable
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper

!use type_forceconstant_secondorder, only: does_forceconstant_fit
implicit none
contains

!> write the forceconstant to file.
module subroutine writetofile(fc, p, fn)
    !> second order force constant
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: fn

    integer :: i, j, a1, a2, u
    real(r8) :: f0
    real(r8), dimension(3) :: v
    character(len=1000) :: opf

    ! Get the actual cutoff.
    f0 = 0.0_r8
    do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%pair(i)%i2
            v = p%rcart(:, a2) - p%rcart(:, a1) + fc%atom(a1)%pair(i)%lv2
            f0 = max(f0, norm2(v))
        end do
    end do

    ! Dump it
    u = open_file('out', trim(fn))
    ! Print the necessary stuff
    write (u, "(1X,I10,15X,'How many atoms per unit cell')") fc%na
    write (u, "(1X,F20.15,5X,'Realspace cutoff (A)')") (f0 + lo_sqtol)*lo_bohr_to_A
    do a1 = 1, fc%na
        write (u, "(1X,I10,15X,'How many neighbours does atom ',I3,' have')") fc%atom(a1)%n, a1
        do i = 1, fc%atom(a1)%n
            opf = "(1X,I10,15X,'In the unit cell, what is the index of neighbour "//tochar(i)//" of atom "//tochar(a1)//"')"
            write (u, opf) fc%atom(a1)%pair(i)%i2
            ! The lattice vector needs to be in reduced coordinates
            v = matmul(p%inv_latticevectors, fc%atom(a1)%pair(i)%lv2)
            do j = 1, 3
                v(j) = anint(v(j))*1.0_r8
            end do
            write (u, *) v
            ! And the actual forceconstant
            do j = 1, 3
                v = fc%atom(a1)%pair(i)%m(j, :)*lo_forceconstant_2nd_HartreeBohr_to_eVA
                write (u, *) v
            end do
        end do
    end do
    ! And all the Born-charge stuff, if necessary
    if (fc%polar .eqv. .false.) then
        write (u, *) '0   # This contains no information about Born charges or dielectric tensors'
    else
        write (u, *) fc%loto%correctiontype, '# This is a forceconstant for a polar material.'
        write (u, "(3(1X,F20.12),' Dielectric tensor xx xy xz')") fc%loto%eps(:, 1)
        write (u, "(3(1X,F20.12),' Dielectric tensor yx yy yz')") fc%loto%eps(:, 2)
        write (u, "(3(1X,F20.12),' Dielectric tensor zx zy zz')") fc%loto%eps(:, 3)
        write (u, *) fc%ew%lambda, ' # Coupling parameter in Ewald summation'
        write (u, *) fc%loto%nx_Z, ' # number of irreducible components in the Born charges'
        do i = 1, fc%loto%nx_Z
            if (i .eq. 1) then
                write (u, "(1X,F20.12,1X,A)") fc%loto%x_Z(i), '# irrep of Born charge'
            else
                write (u, "(1X,F20.12)") fc%loto%x_Z(i)
            end if
        end do
        opf = "(1X,"//tochar(fc%loto%nx_Z)//"(1X,F18.12))"
        do i = 1, fc%na*9
            write (u, opf) fc%loto%coeff_Z(i, :)
        end do
        ! Write the actual Born effective charges. Not really used, but could come handy
        do a1 = 1, fc%na
            write (u, "(3(1X,F20.12),1X,A)") fc%loto%born_effective_charges(:, 1, a1), 'Born effective charge atom '//tochar(a1)
            do i = 2, 3
                write (u, "(3(1X,F20.12))") fc%loto%born_effective_charges(:, i, a1)
            end do
        end do
    end if
    ! Print some auxiliary information, if it is there. Such as norm of
    ! forceconstant per shell, which shells there are and so on.
    if (fc%npairshells .gt. 0 .and. allocated(fc%pairshell)) then
        write (u, "(1X,I10,15X,'Number of irreducible coordination shells')") fc%npairshells
        do i = 1, fc%npairshells
            write (u, "(1X,I10,1X,F16.10,1X,F16.10,15X,'number atoms in shell, radius, norm of forceconstant',I0)") &
                fc%pairshell(i)%n, fc%pairshell(i)%rad*lo_bohr_to_A, fc%pairshell(i)%norm*lo_forceconstant_2nd_HartreeBohr_to_eVA, i
        end do
        do i = 1, fc%npairshells
            do j = 1, fc%pairshell(i)%n
                write (u, "(1X,3(1X,F18.12),2(1X,I0))") lo_chop(matmul(p%inv_latticevectors, fc%pairshell(i)%vec(:, j)), lo_sqtol), fc%pairshell(i)%atind(j), fc%pairshell(i)%pairind(j)
            end do
        end do
    end if
    close (u)
end subroutine

!> Read a force constant from file
module subroutine readfromfile(fc, p, fn, mem, verbosity)
    !> the force constant
    class(lo_forceconstant_secondorder), intent(out) :: fc
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> the filename
    character(len=*), intent(in) :: fn
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> How much to talk
    integer, intent(in) :: verbosity

    integer :: i, j, a1, a2, u
    real(r8), dimension(3) :: v1, v2, lv1, lv2, ucv1, ucv2, r
    real(r8) :: lambda

    ! Set everything to nothing
    fc%na = -1
    fc%cutoff = -1
    fc%elastic_constants_voigt = -1
    fc%elastic_constants_tensor = -1
    fc%npairshells = -1
    fc%npairop = -1
    fc%nifc = -1
    fc%nconstraints = -1
    fc%polar = .false.
    fc%loto%correctiontype = -1
    fc%loto%eps = -1
    fc%loto%nx_Z = -1

    ! Get the stuff from file
    u = open_file('in', trim(fn))
    read (u, *) fc%na
    if (fc%na .ne. p%na) then
        call lo_stop_gracefully(['different number of atoms in "'//trim(fn)//'" and the structure'], lo_exitcode_io, __FILE__, __LINE__)
    end if
    read (u, *) fc%cutoff
    allocate (fc%atom(fc%na))
    do a1 = 1, fc%na
        read (u, *) fc%atom(a1)%n
        allocate (fc%atom(a1)%pair(fc%atom(a1)%n))
        do i = 1, fc%atom(a1)%n
            fc%atom(a1)%pair(i)%i1 = -1
            fc%atom(a1)%pair(i)%i2 = -1
            fc%atom(a1)%pair(i)%lv1 = -1
            fc%atom(a1)%pair(i)%lv2 = -1
            fc%atom(a1)%pair(i)%r = -1
            fc%atom(a1)%pair(i)%m = -1
            fc%atom(a1)%pair(i)%weight = -1
            fc%atom(a1)%pair(i)%irreducible_shell = -1
            fc%atom(a1)%pair(i)%irreducible_operation = -1
        end do
        ! read stuff and set some stuff to something
        do i = 1, fc%atom(a1)%n
            fc%atom(a1)%pair(i)%i1 = a1
            read (u, *) fc%atom(a1)%pair(i)%i2
            read (u, *) fc%atom(a1)%pair(i)%lv2
            do j = 1, 3
                read (u, *) v1
                fc%atom(a1)%pair(i)%m(j, :) = v1*lo_forceconstant_2nd_eVA_to_HartreeBohr
            end do
        end do
    end do
    ! that was all the basic stuff. Now check if it's a polar material
    read (u, *) i
    if (i .eq. 0) then
        ! not a polar material, don't bother
        fc%polar = .false.
        fc%loto%correctiontype = 0
    else
        ! Polar material, have to read a bunch more things
        fc%polar = .true.
        fc%loto%correctiontype = i
        do i = 1, 3
            read (u, *) fc%loto%eps(:, i)
        end do
        read (u, *) lambda
        read (u, *) fc%loto%nx_Z
        allocate (fc%loto%x_Z(fc%loto%nx_Z))
        allocate (fc%loto%coeff_Z(fc%na*9, fc%loto%nx_Z))
        fc%loto%x_Z = 0.0_r8
        fc%loto%coeff_Z = 0.0_r8
        do i = 1, fc%loto%nx_Z
            read (u, *) fc%loto%x_Z(i)
        end do
        do i = 1, fc%na*9
            read (u, *) fc%loto%coeff_Z(i, :)
        end do
        ! and space for the born charges
        allocate (fc%loto%born_effective_charges(3, 3, fc%na))
        allocate (fc%loto%born_onsite_correction(3, 3, fc%na))
        fc%loto%born_effective_charges = 0.0_r8
        fc%loto%born_onsite_correction = 0.0_r8
    end if
    close (u)

    ! Convert from fractional to Cartesian, also calculate the actual cutoff in Cartesian coordinates
    fc%cutoff = 0.0_r8
    do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%pair(i)%i2
            ! get all the vectors right
            ucv1 = p%r(:, a1)
            ucv2 = p%r(:, a2)
            lv1 = 0.0_r8
            lv2 = fc%atom(a1)%pair(i)%lv2
            v1 = lv1 + ucv1
            v2 = lv2 + ucv2
            r = v2 - v1
            fc%atom(a1)%pair(i)%lv1 = lo_chop(p%fractional_to_cartesian(lv1), 1E-13_r8)
            fc%atom(a1)%pair(i)%lv2 = lo_chop(p%fractional_to_cartesian(lv2), 1E-13_r8)
            fc%atom(a1)%pair(i)%r = lo_chop(p%fractional_to_cartesian(r), 1E-13_r8)
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%pair(i)%r))
        end do
    end do
    ! and a tiny tolerance to the cutoff.
    fc%cutoff = fc%cutoff + lo_sqtol

    ! Fix more things if it's polar
    if (fc%polar) then
        ! Set Ewald parameters, charges, Hermiticity and ASR corrections
        call fc%set_ewald_and_enforce_borncharge_hermiticity(p, mem, verbosity, fixlambda=lambda)
        ! Maybe pad and weight it. No, deprecated.
        ! select case(fc%loto%correctiontype)
        ! case(1)
        !     ! Just pad it
        !     call pad_and_weight(fc,p,shells=.false.)
        ! case(2)
        !     ! Pad it
        !     call pad_and_weight(fc,p,shells=.true.)
        !     ! Remove longrange part
        !     call remove_longrange_part_from_forceconstant(fc,p,verb)
        ! case(3)
        !     ! Don't think I should do anything
        ! case default
        !     ! do nothing
        ! end select
    end if
end subroutine

!> dump the forceconstant as an ddb-file, to make abinit happy
module subroutine write_to_anaddb(fc,uc,qgrid,mw,mem)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(inout) :: fc
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> list of q-points
    integer, dimension(3),intent(in) :: qgrid
    type(lo_mpi_helper), intent(inout) :: mw
    type(lo_mem_helper), intent(inout) :: mem

    integer :: unddb, nq
    complex(r8), dimension(:,:,:), allocatable :: dynmat
    real(r8), dimension(:), allocatable :: zion
    real(r8), dimension(:,:), allocatable :: qvectors
    class(lo_qpoint_mesh), allocatable :: qp


    ! Abinit wants dynamical matrices, I have forceconstants. I need to generate a q-mesh
    ! tight enough to ensure that there is no information loss. Will probably pad it a little
    ! to be on the safe side.
!    findqmesh: block
!        type(lo_voronoi_distancetable) :: dt
!        real(r8), dimension(3) :: v0
!        real(r8) :: cutoff
!        integer, dimension(:,:), allocatable :: supercells
!        integer, dimension(3) :: dims
!        integer :: i,j,k,l,nq
!
!        ! I always screw up the cutoff, do it again here to be on the safe side
!        cutoff=0.0_r8
!        do i=1,fc%na
!        do j=1,fc%atom(i)%n
!            cutoff=max(cutoff,norm2(fc%atom(i)%pair(j)%r))
!        enddo
!        enddo
!        ! Get a bunch of reasonable supercells
!!        call generate_possible_supercells(uc,cutoff,supercells)
!!        ! Find one large enough to fit the forceconstant
!!        ssl1: do i=1,size(supercells,2)
!!            call dt%generate_wzcutoff( uc%r,uc%latticevectors,supercells(:,i) )
!!            if ( does_forceconstant_fit(fc,dt) ) then
!!                dims=supercells(:,i)
!!                exit ssl1
!!            endif
!!        enddo ssl1
!        dims=(/6,6,6/)
!        ! add a little safety margin to the dimensions, not like it actually costs anything
!        dims=dims+2
!
!        ! build the list of qvectors
!        nq=product(dims)
!        allocate(qvectors(3,nq))
!        l=0
!        do i=1,dims(1)
!        do j=1,dims(2)
!        do k=1,dims(3)
!            l=l+1
!            v0(1)=(i-1.0_r8)/(1.0_r8*dims(1))
!            v0(2)=(j-1.0_r8)/(1.0_r8*dims(2))
!            v0(3)=(k-1.0_r8)/(1.0_r8*dims(3))
!            qvectors(:,l)=uc%fractional_to_cartesian(v0,reciprocal=.true.)
!        enddo
!        enddo
!        enddo
!    end block findqmesh

    ! Now I have a decent list of q-vectors. For each of those, I need a dynamical matrix.
    getdynmat: block
        integer :: i,nb
        complex(r8), dimension(:,:,:), allocatable :: Dq
        real(r8) :: abiqpts(3,12)

        nb=uc%na*3
        call lo_generate_qmesh(qp,uc,qgrid,'fft',timereversal=.true.,headrankonly=.false.,&
           mw=mw,mem=mem,verbosity=1)
        nq = qp%n_irr_point

!nq = 12
!abiqpts =reshape ((/ & 
!  0.00000000E+00,  0.00000000E+00,  0.00000000E+00,  &
!  2.50000000E-01,  0.00000000E+00,  0.00000000E+00, & 
!  5.00000000E-01,  0.00000000E+00,  0.00000000E+00, & 
!  2.50000000E-01,  2.50000000E-01,  0.00000000E+00, & 
!  0.00000000E+00,  0.00000000E+00,  2.50000000E-01, & 
!  2.50000000E-01,  0.00000000E+00,  2.50000000E-01, & 
!  5.00000000E-01,  0.00000000E+00,  2.50000000E-01, & 
!  2.50000000E-01,  2.50000000E-01,  2.50000000E-01, & 
!  0.00000000E+00,  0.00000000E+00,  5.00000000E-01, & 
!  2.50000000E-01,  0.00000000E+00,  5.00000000E-01, & 
!  5.00000000E-01,  0.00000000E+00,  5.00000000E-01, & 
!  2.50000000E-01,  2.50000000E-01,  5.00000000E-01  /), (/3,12/) )

        allocate(dynmat(nb,nb,nq))
        allocate(Dq(nb,nb,3))
        allocate(qvectors(3,nq))
        do i=1,nq
            qvectors(:,i) = uc%cartesian_to_fractional(qp%ip(i)%r,reciprocal=.true.)

!qvectors(:,i) = abiqpts(:,i)
!qp%ip(i)%r = uc%fractional_to_cartesian(abiqpts(:,i), reciprocal=.true.)
!print *, 'iq qpt ', i, qvectors(:,i)
            call fc%dynamicalmatrix(uc,qp%ip(i),dynmat(:,:,i),mem,Dq,qdirection=(/1.0d0,0.0d0,0.0d0/),skipnonanalytical=.true.)
        enddo
        deallocate(Dq)
    end block getdynmat

    ! Now let's try to print this. Not pretty, but seems to work.
    printstuff: block
        real(r8), dimension(:,:), allocatable :: tnons
        real(r8), dimension(:,:), allocatable :: spinat
        real(r8), dimension(:, :), allocatable :: rotmat
        real(r8), dimension(:), allocatable :: amu,znucl ! dumna,dummass
        integer, dimension(:,:,:), allocatable :: symrel
        integer, dimension(:), allocatable :: symafm,typat
        integer, dimension(103), parameter :: valence_charge = [&
         &1,   2,&
         &3,   4,    3,   4,   5,   6,   7,   8,&
         &9,   2,    3,   4,   5,   6,   7,   8,&
         &9,   10,  11,  12,  13,  14,  15,  16, 17, 18, 19, 12,  3,  4,  5,  6,  7,  8,&
         &9,   10,  11,  12,  13,  14,  15,  16, 17, 18, 19, 12,  3,  4,  5,  6,  7,  8,&
         &9,   10,  11,  12,  13,  14,  15,  16, 17, 18, 19, 20, 21, 22, 23, 24, 11,&
              &12,  13,  14,  15,  16,  17,  18, 19, 12, 13, 14, 15, 16,  7,  8,&
         &9,   10,  11,  12,  13,  14,  15,  16, 17, 18, 19, 20, 21, 22, 23, 24, 11]
!       &'H ',                                                                                                   'He',&
!       &'Li', 'Be',                                                               'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
!       &'Na', 'Mg',                                                               'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',&
!       &'K ', 'Ca',  'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',&
!       &'Rb', 'Sr',  'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',  'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe',&
!       &'Cs', 'Ba',  'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',  'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
!       &                   'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',&
!       &'Fr', 'Ra',  'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',  'Es', 'Fm', 'Md', 'No', 'Lr']

        integer, dimension(18) :: ngfft
        integer :: i,j,a1,a2,ii,jj,x,y,nb
        integer :: u, nq
        integer :: ntypat
        character(len=500) :: string
        character(len=500) :: filnam

        if (uc%info%havespacegroup .eqv. .false.) then
            write (*, *) 'NO SPACEGROUP'
            stop
        end if
        ! some constants
        ntypat=size(uc%element_counter,1)   ! number of atom types
        nb=uc%na*3
        nq = size(qvectors, 2)
        string = "generated by TDEP with lots of dummy variables"
        filnam = 'outfile.many_dynamical_matrices_DDB'
        ngfft = 0
        ngfft(1:3) = [10, 10, 10]
        ! abinit wants many things per atom type
        allocate(amu(ntypat))
        allocate(typat(uc%na))
        allocate(zion(ntypat))
        allocate(znucl(ntypat))
        do i=1,ntypat
            do j=1,uc%na
                if ( uc%species(j) .eq. i ) then
                    zion(i)=valence_charge(uc%atomic_number(j))
                    znucl(i)=uc%atomic_number(j)
                    amu(i)=uc%mass(j)*lo_emu_to_amu
                    typat(j)=i
                endif
            enddo
        enddo
        ! some stuff per atom
        allocate (spinat(3, uc%na))
        spinat = 0.0_r8
        ! stuff per symmetry operation
        allocate (symrel(3, 3, uc%sym%n)) ! some space for the symops
        allocate (tnons(3, uc%sym%n))    ! some space for the reduced translations
        allocate (symafm(uc%sym%n))     ! some space for the antiferro operation flags
        symafm = 1
        do i = 1, uc%sym%n
            symrel(:, :, i) = int(uc%sym%op(i)%fm)
            tnons(:, i) = uc%sym%op(i)%ftr
        end do

        ! Dump the first stuffs
        call ddb_io_out (string, &
            filnam, &
            uc%na,&         ! matom
            1,&             ! mband
            1,&             ! mkpt
            uc%sym%n,&      ! msym
            ntypat,&        ! ntypat (number of atom types)
            unddb,&             ! unddb, the output unit (hopefully returned)
            100401,&        ! vrsddb (no idea)
            [1.0_r8, 1.0_r8, 1.0_r8],& ! acell
            amu,&   ! amu (masses, I think)
            1.0_r8,&      ! dilatmx (no idea)
            1.0_r8,&      ! ecut    (no idea)
            0.0_r8,&      ! ecutsm  (no idea)
            1,&       ! intxc (no idea)
            7,&       ! iscf  (no idea)
            1,&       ! ixc   (no idea)
            [0.0_r8, 0.0_r8, 0.0_r8],& ! kpt
            1.0_r8, &      ! kptnrm
            uc%na,&          ! natom
            [1],&            ! nband
            ngfft,&      ! ngfft
            1,&          ! nkpt
            1,&          ! nspden
            1,&          ! nspinor
            1,&          ! nsppol
            uc%sym%n,&   ! nsym
            maxval(uc%species),&     ! ntypat
            [2.0_r8],&    ! occ
            1, &            ! occopt
            2.0_r8,&      ! pawecutdg
            uc%latticevectors,&       ! rprim
            0.0_r8,&   ! dfpt_sciss
            spinat,&     ! spinat
            symafm,&     ! symafm
            symrel,&     ! symrel
            tnons,&      ! tnons
            1.e-20_r8,&       ! tolwfr
            0.0_r8,&          ! tphysel
            0.001_r8, &       ! tsmear
            uc%species(:),&     ! typat
            0,&          ! usepaw
            [0.0_r8],& ! wtk
            uc%r(:,:),&  ! xred
            zion,&       ! zion
            znucl)       ! znucl
            !TODO: a few of these should be updated to the physical values for the current system: znucl kpt and so on

    end block printstuff

    printdm : block 
        complex(r8), dimension(:,:), allocatable :: dmt
        real(r8), dimension(:,:), allocatable :: rotmat
        real(r8), dimension(3,3) :: bec
        real(r8), dimension(3,3) :: eps
        real(r8), dimension(3,3) :: m1,m0
        integer :: i,j,a1,a2,ii,jj,x,y,nb
        integer :: npert

        m1=0.0d0
        do i=1,3
            m1(i,i)=1.0d0
        enddo

        ! Now print the actual stuff you need:
        write (unddb,*) ""
        write (unddb,*) "No information on the potentials yet"
        write (unddb,*) ""
        write (unddb,*) "**** Database of total energy derivatives ****"
        write (unddb,'(a,i4)') "Number of data blocks=  ", nq

        ! big rotation matrix from cartesian to reduced coordinates
        nb=uc%na*3
        allocate(dmt(nb,nb))
        allocate(rotmat(nb,nb))
        rotmat = 0.0d0
        do a1=1,uc%na
          do x=1,3
            do y=1,3
               ii=(a1-1)*3+x
               jj=(a1-1)*3+y
               !rotmat (ii,jj) = uc%latticevectors(x,y)
               rotmat (ii,jj) = uc%inv_reciprocal_latticevectors(x,y)
            end do
          end do
        end do
        do i=1,nq
            npert =  uc%na*uc%na*3*3
            if (norm2(qvectors(:,i)) < 1.e-10 .and. fc%polar .eqv. .true.) then
               npert = npert + 2*uc%na*3*3 + 3*3
            end if
            write (unddb,*) ""
            write (unddb,'(a32,12x,i8)') " 2nd derivatives (non-stat.)  - # elements : ", npert
            write (unddb,'(a,3es16.8,f6.1)') " qpt", qvectors(:,i), 1.0

            dmt=dynmat(:,:,i)
            ! TODO: replace by a BLAS call?
            !dmt = matmul( transpose(rotmat), matmul(dmt, rotmat) )
            dmt = matmul( rotmat, matmul(dmt, transpose(rotmat)) )
            ! this is in eV/reduced coordinate - convert to Hartree and units of electron mass
            !dmt = dmt / lo_Hartree * lo_u2me
            do a2=1,uc%na
            do y=1,3
                jj=(a2-1)*3+y
                do a1=1,uc%na
                do x=1,3
                    ii=(a1-1)*3+x
                    dmt(ii,jj)=dmt(ii,jj)*sqrt(uc%mass(a1)*uc%mass(a2)) !*lo_emu_to_amu
                    write (unddb,'(4I4,2E23.15)') x,a1,y,a2,&
                    !lo_chop(real(dmt(jj,ii)),lo_sqtol),lo_chop(aimag(dmt(jj,ii)),lo_sqtol)
                    lo_chop(real(dmt(ii,jj)),lo_sqtol),lo_chop(-aimag(dmt(ii,jj)),lo_sqtol)
                enddo
                enddo
            enddo
            enddo
            !Born charges and epsilon, if we are at Gamma
            if (norm2(qvectors(:,i)) < 1.e-10 .and. fc%polar .eqv. .true.) then
               do a2=1,uc%na
                   ! Born Charges for E field=a1 direction x, and atoms a2 direction y + symmetric
                   a1 = uc%na+2
                   bec=fc%loto%born_effective_charges(:,:,a2)
                   m0=bec-m1*zion(uc%species(a2))
                   m0=matmul(uc%inv_reciprocal_latticevectors,matmul(m0,transpose(uc%inv_latticevectors)))*lo_twopi
                   do y=1,3
                   do x=1,3
                       write (unddb,'(4I4,2E23.15)') x,a1,y,a2,&
                         lo_chop(m0(x,y),lo_sqtol), 0.0d0
                       write (unddb,'(4I4,2E23.15)') y,a2,x,a1,&
                         lo_chop(m0(y,x),lo_sqtol), 0.0d0
                   enddo
                   enddo
               enddo

               ! epsilon
               a2 = uc%na+2
               a1 = uc%na+2
               eps = fc%loto%eps(:,:)
               m0 = (m1-eps)*uc%volume/(2.0*lo_twopi)
               m0 = m0*(lo_twopi**2)
               m0 = matmul(uc%inv_latticevectors,matmul(m0,transpose(uc%inv_latticevectors)))
               do y=1,3
               do x=1,3
                   write (unddb,'(4I4,2E23.15)') x,a1,y,a2,&
                     lo_chop(m0(x,y),lo_sqtol), 0.0d0
               enddo
               enddo
            endif
        enddo

        close(unddb)
        deallocate(dmt)
        deallocate(rotmat)
    end block printdm
    deallocate(dynmat)
    deallocate(qvectors)
end subroutine

!> mini-version of my dynamical matrix calculator, repeated here to avoid circular dependencies
subroutine simple_dynmat(fc, uc, qv, dynmat)
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(in) :: fc
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> q-vector
    real(r8), dimension(3), intent(in) :: qv
    !> dynamical matric
    complex(r8), dimension(:, :), intent(out) :: dynmat

    complex(r8), dimension(3, 3) :: m0
    complex(r8) :: exp_ikr
    real(r8) :: k_dot_r
    integer :: a1, a2, i, ii, jj, mu, nu
    !
    dynmat = 0.0_r8
    do a1 = 1, fc%na
    do i = 1, fc%atom(a1)%n
        a2 = fc%atom(a1)%pair(i)%i2
        k_dot_r = dot_product(fc%atom(a1)%pair(i)%lv2, qv)*lo_twopi
        exp_ikr = cmplx(cos(k_dot_r), sin(k_dot_r), r8)
        m0 = exp_ikr*fc%atom(a1)%pair(i)%m*uc%invsqrtmass(a1)*uc%invsqrtmass(a2)
        do mu = 1, 3
        do nu = 1, 3
            ii = (a2 - 1)*3 + nu
            jj = (a1 - 1)*3 + mu
            dynmat(ii, jj) = dynmat(ii, jj) + m0(mu, nu)
        end do
        end do
    end do
    end do
end subroutine

! The following is a routine yanked from Abinit to get a reasonable header in the anaddb-file.
! Since this module has everything marked private by default, it should still work fine to link
! this with Abinit.

!!****f* ABINIT/ddb_io_out
!!
!! NAME
!! ddb_io_out
!!
!! FUNCTION
!! Open Derivative DataBase, then
!! reads or write Derivative DataBase preliminary information.
!! Note: only one processor read or write the DDB.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! amu(mtypat)=mass of the atoms (atomic mass unit)
!! dilatmx=the maximal dilatation factor
!! character(len=500) dscrpt:string that describe the output database
!! ecut=kinetic energy planewave cutoff (hartree)
!! ecutsm=smearing energy for plane wave kinetic energy (Ha)
!! character(len=500) filnam: name of output file
!! intxc=control xc quadrature
!! iscf=parameter controlling scf or non-scf choice
!! ixc=exchange-correlation choice parameter
!! kpt(3,mkpt)=k point set (reduced coordinates)
!! kptnrm=normalisation of k points
!! matom=maximum number of atoms
!! mband=maximum number of bands
!! mkpt=maximum number of special points
!! msym=maximum number of symetries
!! mtypat=maximum number of atom types
!! natom=number of atoms in the unit cell
!! nband(mkpt)=number of bands at each k point, for each polarization
!! ngfft(18)=contain all needed information about 3D FFT,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkpt=number of k points
!! nspden=number of spin-density components
!! nspinor=number of spinorial components of the wavefunctions
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements in space group
!! ntypat=number of atom types
!! occ(mband*mkpt)=occupation number for each band and k
!! occopt=option for occupancies
!! pawecutdg=cut-off for fine "double grid" used in PAW calculations (unused for NCPP)
!! rprim(3,3)=dimensionless primitive translations in real space
!! dfpt_sciss=scissor shift (Ha)
!! spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym)=symmetry operations in real space
!! tnons(3,msym)=nonsymmorphic translations for symmetry operations
!! tolwfr=tolerance on largest wf residual
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature) in Hartree
!! typat(matom)=type of each atom
!! unddb=unit number for output
!! usepaw=flag for PAW
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.
!! wtk(mkpt)=weight assigned to each k point
!! xred(3,matom)=reduced atomic coordinates
!! zion(mtypat)=valence charge of each type of atom
!! znucl(mtypat)=atomic number of atom type
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,gstate,mblktyp1,mblktyp5,nonlinear,respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddb_io_out(dscrpt, filnam, matom, mband,&
&  mkpt, msym, mtypat, unddb, vrsddb,&
&  acell, amu, dilatmx, ecut, ecutsm, intxc, iscf, ixc, kpt, kptnrm,&
&  natom, nband, ngfft, nkpt, nspden, nspinor, nsppol, nsym, ntypat, occ, occopt,&
&  pawecutdg, rprim, dfpt_sciss, spinat, symafm, symrel, tnons, tolwfr, tphysel, tsmear,&
&  typat, usepaw, wtk, xred, zion, znucl)

    implicit none

!Arguments -------------------------------
!scalars
    integer, intent(in) :: matom, mband, mkpt, msym, mtypat, unddb, vrsddb
    integer, intent(in) :: intxc, iscf, ixc, natom, nkpt, nspden, nspinor, nsppol, nsym
    integer, intent(in) :: ntypat, occopt, usepaw
    real(r8), intent(in) :: dilatmx, ecut, ecutsm, kptnrm, pawecutdg, dfpt_sciss, tolwfr, tphysel
    real(r8), intent(in) :: tsmear
    character(len=500), intent(in) :: dscrpt, filnam
!arrays
    integer, intent(in) :: nband(mkpt), ngfft(18), symafm(msym), symrel(3, 3, msym)
    integer, intent(in) :: typat(matom)
    real(r8), intent(in) :: acell(3), amu(mtypat), kpt(3, mkpt), occ(mband*mkpt)
    real(r8), intent(in) :: rprim(3, 3), spinat(3, matom), tnons(3, msym), wtk(mkpt)
    real(r8), intent(in) :: xred(3, matom), zion(mtypat), znucl(mtypat)

!Local variables -------------------------
!Set routine version number here:
!scalars
    integer, parameter :: vrsio8 = 100401
    integer :: bantot, ii, ij, ikpt, iline, im
    integer :: ioerr
    character(len=500) :: message
!arrays
    character(len=9) :: name(9)

! *********************************************************************

!Check ioddb8 version number (vrsio8) against mkddb version number (vrsddb)
    if (vrsio8 /= vrsddb) then
        write (message, '(a,a,a,i10,a,a,i10,a)')&
     &   ' ddb_io_out: WARNING -', char(10),&
     &   '  The input/output DDB version number=', vrsio8, char(10),&
     &   '  is not equal to the DDB version number=', vrsddb, '.'
        write (6, *) message
    end if

!Open the output derivative database.
!(version 2.1. : changed because of a bug in a Perl script
!should set up a name checking procedure, with change of name
!like for the output file)
    open (file=filnam, unit=unddb, status='unknown', form='formatted', iostat=ioerr)
    if (ioerr /= 0) then
        write (*, *) message
        stop
    end if

!Write the heading
    write (unddb, '(/,a,/,a,i10,/,/,a,a,/)') &
   & ' **** DERIVATIVE DATABASE ****    ',&
   & '+DDB, Version number', vrsddb, ' ', trim(dscrpt)

!Write the descriptive data
    !1. usepaw
    write (unddb, '(1x,a9,i10)') '   usepaw', usepaw
    !2. natom
    write (unddb, '(1x,a9,i10)') '    natom', natom
    !3. nkpt
    write (unddb, '(1x,a9,i10)') '     nkpt', nkpt
    !4. nsppol
    write (unddb, '(1x,a9,i10)') '   nsppol', nsppol
    !5. nsym
    write (unddb, '(1x,a9,i10)') '     nsym', nsym
    !6. ntypat
    write (unddb, '(1x,a9,i10)') '   ntypat', ntypat
    !7. occopt
    write (unddb, '(1x,a9,i10)') '   occopt', occopt
    !8. nband
    if (occopt == 2) then
        im = 12
        name(1) = '    nband'
        do iline = 1, (nkpt + 11)/12
            if (iline == (nkpt + 11)/12) im = nkpt - 12*(iline - 1)
            write (unddb, '(1x,a9,5x,12i5)') name(1), (nband((iline - 1)*12 + ii), ii=1, im)
            name(1) = '         '
        end do
        bantot = 0
        do ikpt = 1, nkpt
            bantot = bantot + nband(ikpt)
        end do
    else
        write (unddb, '(1x,a9,i10)') '    nband', nband(1)
        bantot = nkpt*nband(1)
    end if

    !9. acell
    write (unddb, '(1x,a9,3d22.14)') '    acell', acell
    !10. amu
    im = 3
    name(1) = '      amu'
    do iline = 1, (ntypat + 2)/3
        if (iline == (ntypat + 2)/3) im = ntypat - 3*(iline - 1)
        write (unddb, '(1x,a9,3d22.14)') name(1), (amu((iline - 1)*3 + ii), ii=1, im)
        name(1) = '         '
    end do
    !11. dilatmx
    write (unddb, '(1x,a9,d22.14)') '  dilatmx', dilatmx
    !12. ecut
    write (unddb, '(1x,a9,d22.14)') '     ecut', ecut
    !12b. pawecutdg (PAW)
    if (usepaw == 1) then
        write (unddb, '(1x,a9,d22.14)') 'pawecutdg', pawecutdg
    end if
    !13. ecutsm
    write (unddb, '(1x,a9,d22.14)') '   ecutsm', ecutsm
    !14. intxc
    write (unddb, '(1x,a9,i10)') '    intxc', intxc
    !15. iscf
    write (unddb, '(1x,a9,i10)') '     iscf', iscf
    !16. ixc
    write (unddb, '(1x,a9,i10)') '      ixc', ixc
    !17. kpt
    name(1) = '      kpt'
    do iline = 1, nkpt
        write (unddb, '(1x,a9,3d22.14)') name(1), (kpt(ii, iline), ii=1, 3)
        name(1) = '      '
    end do
    !18. kptnrm
    write (unddb, '(1x,a9,d22.14)') '   kptnrm', kptnrm
    !19. ngfft
    write (unddb, '(1x,a9,5x,3i5)') '    ngfft', ngfft(1:3)
    !20. nspden
    write (unddb, '(1x,a9,i10)') '   nspden', nspden
    !21. nspinor
    write (unddb, '(1x,a9,i10)') '  nspinor', nspinor
    !22. occ
    if (occopt == 2) then
        im = 3
        name(1) = '      occ'
        do iline = 1, (bantot + 2)/3
            if (iline == (bantot + 2)/3) im = bantot - 3*(iline - 1)
            write (unddb, '(1x,a9,3d22.14)') name(1), (occ((iline - 1)*3 + ii), ii=1, im)
            name(1) = '         '
        end do
    else
        im = 3
        name(1) = '      occ'
        do iline = 1, (nband(1) + 2)/3
            if (iline == (nband(1) + 2)/3) im = nband(1) - 3*(iline - 1)
            write (unddb, '(1x,a9,3d22.14)') name(1), (occ((iline - 1)*3 + ii), ii=1, im)
            name(1) = '         '
        end do
    end if
    !23. rprim
    name(1) = '    rprim'
    do iline = 1, 3
        write (unddb, '(1x,a9,3d22.14)') name(1), (rprim(ii, iline), ii=1, 3)
        name(1) = '      '
    end do
    !24. dfpt_sciss
    write (unddb, '(1x,a11,d22.14)') ' dfpt_sciss', dfpt_sciss
    !25. spinat
    name(1) = '   spinat'
    do iline = 1, natom
        write (unddb, '(1x,a9,3d22.14)') name(1), (spinat(ii, iline), ii=1, 3)
        name(1) = '         '
    end do
    !26. symafm
    im = 12
    name(1) = '   symafm'
    do iline = 1, (nsym + 11)/12
        if (iline == (nsym + 11)/12) im = nsym - 12*(iline - 1)
        write (unddb, '(1x,a9,5x,12i5)') name(1), (symafm((iline - 1)*12 + ii), ii=1, im)
        name(1) = '         '
    end do
    !27. symrel
    name(1) = '   symrel'
    do iline = 1, nsym
        write (unddb, '(1x,a9,5x,9i5)') name(1), ((symrel(ii, ij, iline), ii=1, 3), ij=1, 3)
        name(1) = '         '
    end do
    !28. tnons
    name(1) = '    tnons'
    do iline = 1, nsym
        write (unddb, '(1x,a9,3d22.14)') name(1), (tnons(ii, iline), ii=1, 3)
        name(1) = '         '
    end do
    !29. tolwfr
    write (unddb, '(1x,a9,d22.14)') '   tolwfr', tolwfr
    !30. tphysel
    write (unddb, '(1x,a9,d22.14)') '  tphysel', tphysel
    !31. tsmear
    write (unddb, '(1x,a9,d22.14)') '   tsmear', tsmear
    !32. typat
    im = 12
    name(1) = '    typat'
    do iline = 1, (natom + 11)/12
        if (iline == (natom + 11)/12) im = natom - 12*(iline - 1)
        write (unddb, '(1x,a9,5x,12i5)') name(1), (typat((iline - 1)*12 + ii), ii=1, im)
        name(1) = '         '
    end do
    !33. wtk
    name(1) = '      wtk'
    im = 3
    do iline = 1, (nkpt + 2)/3
        if (iline == (nkpt + 2)/3) im = nkpt - 3*(iline - 1)
        write (unddb, '(1x,a9,3d22.14)') name(1), (wtk((iline - 1)*3 + ii), ii=1, im)
        name(1) = '         '
    end do
    !34. xred
    name(1) = '     xred'
    do iline = 1, natom
        write (unddb, '(1x,a9,3d22.14)') name(1), (xred(ii, iline), ii=1, 3)
        name(1) = '         '
    end do
    !35. znucl
    name(1) = '    znucl'
    write (*, *) 'znucl', znucl, size(znucl, 1)
    im = 3
    do iline = 1, (ntypat + 2)/3
        if (iline == (ntypat + 2)/3) im = ntypat - 3*(iline - 1)
        write (unddb, '(1x,a9,3d22.14)') name(1), (znucl((iline - 1)*3 + ii), ii=1, im)
        name(1) = '         '
    end do
    !36. zion
    name(1) = '     zion # NB: may not match original pseudopotential, but will give the right Born charges'
    write (*, *) 'zion', zion
    im = 3
    do iline = 1, (ntypat + 2)/3
        if (iline == (ntypat + 2)/3) im = ntypat - 3*(iline - 1)
        write (unddb, '(1x,a9,3d22.14)') name(1), (zion((iline - 1)*3 + ii), ii=1, im)
        name(1) = '         '
    end do
end subroutine ddb_io_out

end submodule
