#include "precompilerdefinitions"
module type_jij_secondorder
    use konstanter, only: flyt,lo_twopi,lo_tol,lo_sqtol,lo_huge,lo_hugeint,lo_exitcode_io,lo_Hartree_to_eV,lo_eV_to_Hartree,&
                          lo_Bohr_to_A,lo_A_to_Bohr
    use gottochblandat, only: open_file,walltime,lo_chop,tochar,lo_verletbox,lo_trueNtimes,&
                       lo_progressbar_init,lo_progressbar,lo_stop_gracefully,lo_trace
    use type_crystalstructure, only: lo_crystalstructure
    use type_distancetable, only: lo_distancetable
    implicit none

    private
    public :: lo_jij_secondorder

    !> information for a single pair
    type lo_jij2_pair
        !> index in the unit cell to atom 1 and atom 2
        integer :: i1=-lo_hugeint,i2=-lo_hugeint
        !> lattice vectors positioning the unit cell
        real(flyt), dimension(3) :: lv1=lo_huge,lv2=lo_huge
        !> vector between atom 1 and 2
        real(flyt), dimension(3) :: r=lo_huge
        !> tensorial Jij
        real(flyt), dimension(3,3) :: J=lo_huge
        !> lattice-spin term
        real(flyt), dimension(3,3) :: bird=lo_huge
        !> spin-lattice term
        real(flyt), dimension(3,3) :: dude=lo_huge
        !> longfluct term
        real(flyt), dimension(3) :: Q=lo_huge
        !> which irreducible shell does it belong to?
        integer :: irreducible_shell=-lo_hugeint
        !> what operation takes the shell to here
        integer :: irreducible_operation=-lo_hugeint
    end type

    !> list of pairs per atom
    type lo_jij2_atom
        !> how many pairs
        integer :: n=-lo_hugeint
        !> mean magnetic moment
        real(flyt) :: m0=-lo_huge
        !> curvature around magnetic moments
        real(flyt) :: m0dev=-lo_huge
        !> self-term for bird & dude
        real(flyt), dimension(3,3) :: selfterm_bird=lo_huge
        real(flyt), dimension(3,3) :: selfterm_dude=lo_huge
        !> information about each pair
        type(lo_jij2_pair), dimension(:), allocatable :: pair
    end type

    !> Secondorder magnetic interactions
    type lo_jij_secondorder
        !> How much info to write?
        integer :: verbosity=-lo_hugeint
        !> Number of atoms in the cell
        integer :: na=-lo_hugeint
        !> Length of the longest pair + a tiny margin
        real(flyt) :: cutoff=lo_huge
        !> List of atoms
        type(lo_jij2_atom), allocatable, dimension(:) :: atom
        contains
            !> write to file
            procedure :: writetofile
            !> read from file
            procedure :: readfromfile
            !> remap it to another force constant
            procedure :: remap
            !> get magnon energies
            procedure :: magnonenergy
            !> create a fake prototype
            procedure :: fake_jij
            !> evaluate the onsite energy
            !procedure :: longitudinal_energy
    end type
contains

!pure function longitudinal_energy(jij,m,m0,sigma) result(e)
!    !> jij thingy
!    class(lo_jij_secondorder), intent(in) :: jij
!    !> current magnetic moment
!    real(flyt), intent(in) :: m
!    !> current magnetic reference
!    real(flyt), intent(in) :: m0
!    !> curvature
!    real(flyt), intent(in) :: sigma
!    !> energy
!    real(flyt) :: e
!
!
!    e=sigma*(m**2-m0**2)**2
!end function

!> Create fake exchange interactions to seed calculations with
subroutine fake_jij(jij,uc,J,m0)
    !> Heisenberg pair interactions
    class(lo_jij_secondorder), intent(out) :: jij
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> nearest neighbour interaction
    real(flyt), intent(in) :: J
    !> mean magnetic moment
    real(flyt), intent(in) :: m0

    type(lo_distancetable) :: dt
    integer :: a1,i,k

    jij%verbosity=0
    jij%na=uc%na
    jij%cutoff=uc%nearest_neighbour_distance()+lo_tol
    call dt%generate(uc%r,uc%latticevectors,jij%cutoff,verbosity=0)
    !
    lo_allocate(jij%atom(jij%na))
    do a1=1,jij%na
        jij%atom(a1)%m0=m0
        jij%atom(a1)%m0dev=0.0_flyt
        jij%atom(a1)%selfterm_bird=0.0_flyt
        jij%atom(a1)%selfterm_dude=0.0_flyt
        jij%atom(a1)%n=dt%particle(a1)%n-1
        lo_allocate(jij%atom(a1)%pair( jij%atom(a1)%n ))
        do i=1,jij%atom(a1)%n
            jij%atom(a1)%pair(i)%i1=a1
            jij%atom(a1)%pair(i)%i2=dt%particle(a1)%ind(i+1)
            jij%atom(a1)%pair(i)%lv1=0.0_flyt
            jij%atom(a1)%pair(i)%lv2=dt%particle(a1)%lv(:,i+1)
            jij%atom(a1)%pair(i)%r=dt%particle(a1)%v(:,i+1)
            jij%atom(a1)%pair(i)%J=0.0_flyt
            jij%atom(a1)%pair(i)%bird=0.0_flyt
            jij%atom(a1)%pair(i)%dude=0.0_flyt
            jij%atom(a1)%pair(i)%Q=0.0_flyt
            jij%atom(a1)%pair(i)%irreducible_shell=0
            jij%atom(a1)%pair(i)%irreducible_operation=0
            do k=1,3
                jij%atom(a1)%pair(i)%J(k,k)=J
            enddo
        enddo
    enddo
end subroutine

!> Calculate the magnon energy at a given q-point
subroutine magnonenergy(jij,uc,qvec,energies)
    !> Heisenberg pair interactions
    class(lo_jij_secondorder), intent(in) :: jij
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> q-vector
    real(flyt), dimension(3), intent(in) :: qvec
    !> magnon energies
    real(flyt), dimension(:), intent(out) :: energies

    complex(flyt), dimension(:,:), allocatable :: DM,DG
    complex(flyt) :: expiqr
    real(flyt), parameter :: third=1.0_flyt/3.0_flyt
    real(flyt) :: qdotr,f0,im1,im2
    integer :: a1,a2,i

    lo_allocate(DM(jij%na,jij%na))
    lo_allocate(DG(jij%na,jij%na))

    DM=0.0_flyt
    DG=0.0_flyt
    do a1=1,uc%na !jij%na
    do i=1,jij%atom(a1)%n
        a2=jij%atom(a1)%pair(i)%i2
        im1=1.0_flyt/jij%atom(a1)%m0
        im2=1.0_flyt/jij%atom(a2)%m0
        !
        qdotr=dot_product(jij%atom(a1)%pair(i)%lv2,qvec)*lo_twopi
        expiqr=cmplx(cos(qdotr),sin(qdotr),flyt)
        !
        f0=lo_trace(jij%atom(a1)%pair(i)%J)*third
        if ( a1 .eq. a2 ) then
            DM(a1,a2)=DM(a1,a2)-f0*expiqr*im1
            DG(a1,a2)=DG(a1,a2)-f0*im1
        elseif ( a2 .gt. a1 ) then
            DM(a1,a2)=DM(a1,a2)-f0*expiqr*im1
            DG(a1,a2)=DG(a1,a2)-f0*im1
        else
            DM(a1,a2)=DM(a1,a2)-f0*conjg(expiqr)*im2
            DG(a1,a2)=DG(a1,a2)-f0*im2
        endif
    enddo
    enddo

    ! Acoustic sum rule thing
    do a1=1,jij%na
        f0=0.0_flyt
        do a2=1,jij%na
            f0=f0+real(DG(a2,a1))
        enddo
        DG(:,a1)=0.0_flyt
        DG(a1,a1)=f0
    enddo
    ! Remove acoustic sum thing
    DM=-(DM-DG)*4

    ! Eigenvalue equation?
    if ( jij%na .gt. 1 ) then
        write(*,*) 'FIXME EIGENVALUE SPIN DISPERSIONS'
    else
        ! Only one atom, no eigenvalue equation to bother with.
        energies(1)=lo_chop(real(DM(1,1)),1E-10_flyt)
    endif

end subroutine

!> Map the force constant to a larger cell
subroutine remap(jij,uc,ss,jijss)
    !> the unitcell forceconstant
    class(lo_jij_secondorder), intent(in) :: jij
    !> the unit cell
    type(lo_crystalstructure), intent(in) :: uc
    !> the supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> the supercell forceconstant
    type(lo_jij_secondorder), intent(out) :: jijss
    !
    type(lo_verletbox) :: vb
    real(flyt), dimension(3) :: r0,r1,r2
    real(flyt) :: t0,t1
    integer :: uca,a1,a2,i
    logical :: largecell

    if ( jij%verbosity .gt. 0 ) then
        t0=walltime()
    endif

    ! Just a small sanity check
    if ( ss%info%supercell .eqv. .false. ) then
        write(*,*) 'The supercell needs to be related to the unitcell'
        stop
    endif
    if ( uc%na .ne. jij%na ) then
        write(*,*) 'Different number of atoms in unitcell and forceconstant.'
        stop
    endif

    if ( ss%na .gt. 1000 ) then
        largecell=.true.
        t1=walltime()
        if ( jij%verbosity .gt. 0 ) write(*,*) '... found a large cell: ',tochar(ss%na), 'atoms'
        ! In this case, build Verlet boxes
        call vb%generate(ss%r,[20,20,20])
    else
        largecell=.false.
    endif

    if ( jij%verbosity .gt. 0 ) call lo_progressbar_init()

    ! Build empty forceconstants
    jijss%na=ss%na
    jijss%cutoff=jij%cutoff
    lo_allocate(jijss%atom(jijss%na))
    do a1=1,jijss%na
        r0=ss%rcart(:,a1)
        uca=ss%info%index_in_unitcell(a1)
        jijss%atom(a1)%n=jij%atom(uca)%n
        jijss%atom(a1)%m0=jij%atom(uca)%m0
        jijss%atom(a1)%m0dev=jij%atom(uca)%m0dev
        jijss%atom(a1)%selfterm_bird=jij%atom(uca)%selfterm_bird
        jijss%atom(a1)%selfterm_dude=jij%atom(uca)%selfterm_dude
        lo_allocate(jijss%atom(a1)%pair( jijss%atom(a1)%n ))
        do i=1,jijss%atom(a1)%n
            jijss%atom(a1)%pair(i)%irreducible_shell=jij%atom(uca)%pair(i)%irreducible_shell
            jijss%atom(a1)%pair(i)%irreducible_operation=jij%atom(uca)%pair(i)%irreducible_operation
            jijss%atom(a1)%pair(i)%J=jij%atom(uca)%pair(i)%J
            jijss%atom(a1)%pair(i)%bird=jij%atom(uca)%pair(i)%bird
            jijss%atom(a1)%pair(i)%dude=jij%atom(uca)%pair(i)%dude
            jijss%atom(a1)%pair(i)%Q=jij%atom(uca)%pair(i)%Q
            jijss%atom(a1)%pair(i)%r=jij%atom(uca)%pair(i)%r
            jijss%atom(a1)%pair(i)%lv1=0.0_flyt
            jijss%atom(a1)%pair(i)%lv2=0.0_flyt
            jijss%atom(a1)%pair(i)%i1=a1
            jijss%atom(a1)%pair(i)%i2=0
            ! Now, the only thing missing are the indices. Try and fix that.
            r1=jijss%atom(a1)%pair(i)%r+r0
            r1=ss%cartesian_to_fractional(r1)
            if ( largecell ) then
                a2=vb%locate(ss%r,r1)
                if ( a2 .gt. 0 ) jijss%atom(a1)%pair(i)%i2=a2
            else
                do a2=1,ss%na
                    if ( sum(abs(ss%r(:,a2)-r1)) .lt. lo_tol ) then
                        jijss%atom(a1)%pair(i)%i2=a2
                        exit
                    endif
                enddo
            endif
            ! sanity check
            if ( jijss%atom(a1)%pair(i)%i2 .eq. 0 ) then
                write(*,*) 'ERROR: failed mapping secondorder Jij for atom',a1
                stop
            endif
            ! fix the a lattice vector?
            ! r = lv2 + r2 - r1
            ! lv2 = r - r2 + r1
            r2=jijss%atom(a1)%pair(i)%r - ss%rcart(:,jijss%atom(a1)%pair(i)%i2)+ss%rcart(:,jijss%atom(a1)%pair(i)%i1)
            r2=anint(matmul(ss%inv_latticevectors,r2))
            jijss%atom(a1)%pair(i)%lv2=matmul(ss%latticevectors,r2)
        enddo
        if ( jij%verbosity .gt. 0 ) then
            if ( lo_trueNtimes(a1,200,jijss%na) ) call lo_progressbar(' ... remapping pair Jij',a1,jijss%na,walltime()-t0)
        endif
    enddo
    if ( jij%verbosity .gt. 0 ) call lo_progressbar(' ... remapping pair Jij',jijss%na,jijss%na,walltime()-t0)
end subroutine

!> write the forceconstant to file.
subroutine writetofile(jij,p,fn)
    !> the force constant
    class(lo_jij_secondorder), intent(in) :: jij
    !> the associated crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> the filename
    character(len=*), intent(in) :: fn
    !
    integer :: i,j,a1,a2,u
    real(flyt) :: f0
    real(flyt), dimension(3) :: v0,v1,v2
    character(len=1000) :: opf

    ! Get the actual cutoff.
    f0=0.0_flyt
    do a1=1,jij%na
        do i=1,jij%atom(a1)%n
            a2=jij%atom(a1)%pair(i)%i2
            v0=p%rcart(:,a2)-p%rcart(:,a1)+jij%atom(a1)%pair(i)%lv2
            f0=max(f0,norm2(v0))
        enddo
    enddo

    ! Dump it
    u=open_file('out',trim(fn))
        ! Print the necessary stuff
        write(u,"(1X,I10,15X,'How many atoms per unit cell')") jij%na
        write(u,"(1X,F20.15,5X,'Realspace cutoff (A)')") (f0+lo_sqtol)*lo_Bohr_to_A
        do a1=1,jij%na
            write(u,"(1X,I10,15X,'Number of pairs of atom ',I3)") jij%atom(a1)%n,a1
            write(u,"(1X,E19.12,1X,E19.12,7X,'Average magnetic moment and standard deviation atom',I3)") jij%atom(a1)%m0,jij%atom(a1)%m0dev,a1
            do i=1,jij%atom(a1)%n
                opf="(1X,I6,1X,I4,1X,I2,10X,'Index of neighbour, irreducible shell, irreducible operation')"
                write(u,opf) jij%atom(a1)%pair(i)%i2,jij%atom(a1)%pair(i)%irreducible_shell,jij%atom(a1)%pair(i)%irreducible_operation
                ! The lattice vector needs to be in reduced coordinates
                v0=matmul(p%inv_latticevectors,jij%atom(a1)%pair(i)%lv2)
                do j=1,3
                    v0(j)=anint(v0(j))*1.0_flyt
                enddo
                write(u,"(1X,3(1X,F19.14))") lo_chop(v0,lo_sqtol)
                ! And the actual forceconstant
                do j=1,3
                    v0=jij%atom(a1)%pair(i)%J(j,:)*lo_Hartree_to_eV
                    v1=jij%atom(a1)%pair(i)%bird(j,:)*lo_Hartree_to_eV
                    v2=jij%atom(a1)%pair(i)%dude(j,:)*lo_Hartree_to_eV
                    f0=jij%atom(a1)%pair(i)%Q(j)*lo_Hartree_to_eV
                    write(u,"(10(1X,E19.12))") v0,v1,v2,f0
                    !v3=jij%atom(a1)%pair(i)%Q(j,:)*lo_Hartree_to_eV
                    !write(u,"(12(1X,E19.12))") v0,v1,v2,v3
                enddo
            enddo
        enddo
    close(u)
end subroutine

!> Read exchange parameters from file
subroutine readfromfile(jij,p,fn,verbosity)
    !> the force constant
    class(lo_jij_secondorder), intent(out) :: jij
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> the filename
    character(len=*), intent(in) :: fn
    !> How much to talk
    integer, intent(in), optional :: verbosity
    !
    integer :: i,j,a1,a2,u
    real(flyt), dimension(3) :: v0,v1,v2,lv1,lv2,ucv1,ucv2,r
    real(flyt) :: f0

    ! Set everything to nothing
    jij%na=-1
    jij%cutoff=-1
    if ( present(verbosity) ) then
        jij%verbosity=verbosity
    else
        jij%verbosity=0
    endif

    ! Get the stuff from file
    u=open_file('in',trim(fn))
        read(u,*) jij%na
        if ( jij%na .ne. p%na ) then
            call lo_stop_gracefully(['Different number of atoms in "'//trim(adjustl(fn))//'" and the structure'],lo_exitcode_io,__FILE__,__LINE__)
        endif
        read(u,*) jij%cutoff
        jij%cutoff=jij%cutoff*lo_A_to_Bohr
        lo_allocate(jij%atom(jij%na))
        do a1=1,jij%na
            read(u,*) jij%atom(a1)%n
            read(u,*) jij%atom(a1)%m0,jij%atom(a1)%m0dev
            lo_allocate(jij%atom(a1)%pair( jij%atom(a1)%n ) )
            ! set everything to nothing
            do i=1,jij%atom(a1)%n
                jij%atom(a1)%pair(i)%i1=-1
                jij%atom(a1)%pair(i)%i2=-1
                jij%atom(a1)%pair(i)%lv1=-1
                jij%atom(a1)%pair(i)%lv2=-1
                jij%atom(a1)%pair(i)%r=-1
                jij%atom(a1)%pair(i)%J=-1
                jij%atom(a1)%pair(i)%bird=-1
                jij%atom(a1)%pair(i)%dude=-1
                jij%atom(a1)%pair(i)%Q=-1
                jij%atom(a1)%pair(i)%irreducible_shell=-1
                jij%atom(a1)%pair(i)%irreducible_operation=-1
            enddo
            ! read stuff and set some stuff to something
            do i=1,jij%atom(a1)%n
                jij%atom(a1)%pair(i)%i1=a1
                read(u,*) jij%atom(a1)%pair(i)%i2,jij%atom(a1)%pair(i)%irreducible_shell,jij%atom(a1)%pair(i)%irreducible_operation
                read(u,*) jij%atom(a1)%pair(i)%lv2
                do j=1,3
                    read(u,*) v0,v1,v2,f0
                    jij%atom(a1)%pair(i)%J(j,:)=lo_chop(v0*lo_eV_to_Hartree,1E-14_flyt)
                    jij%atom(a1)%pair(i)%bird(j,:)=lo_chop(v1*lo_eV_to_Hartree,1E-14_flyt)
                    jij%atom(a1)%pair(i)%dude(j,:)=lo_chop(v2*lo_eV_to_Hartree,1E-14_flyt)
                    jij%atom(a1)%pair(i)%Q(j)=lo_chop(f0*lo_eV_to_Hartree,1E-14_flyt)
                enddo
            enddo

            ! Set the self-terms
            jij%atom(a1)%selfterm_bird=0.0_flyt
            jij%atom(a1)%selfterm_dude=0.0_flyt
            do i=1,jij%atom(a1)%n
                jij%atom(a1)%selfterm_bird=jij%atom(a1)%selfterm_bird-jij%atom(a1)%pair(i)%bird
                jij%atom(a1)%selfterm_dude=jij%atom(a1)%selfterm_dude-jij%atom(a1)%pair(i)%dude
            enddo
            jij%atom(a1)%selfterm_bird=lo_chop( jij%atom(a1)%selfterm_bird, 1E-15_flyt )
            jij%atom(a1)%selfterm_dude=lo_chop( jij%atom(a1)%selfterm_dude, 1E-15_flyt )
        enddo
    close(u)

    ! Convert from fractional to Cartesian, also calculate the actual cutoff in Cartesian coordinates
    jij%cutoff=0.0_flyt
    do a1=1,jij%na
        do i=1,jij%atom(a1)%n
            a2=jij%atom(a1)%pair(i)%i2
            ! get all the vectors right
            ucv1=p%r(:,a1)
            ucv2=p%r(:,a2)
            lv1=0.0_flyt
            lv2=jij%atom(a1)%pair(i)%lv2
            v1=lv1+ucv1
            v2=lv2+ucv2
            r=v2-v1
            jij%atom(a1)%pair(i)%lv1=lo_chop(p%fractional_to_cartesian(lv1),lo_sqtol)
            jij%atom(a1)%pair(i)%lv2=lo_chop(p%fractional_to_cartesian(lv2),lo_sqtol)
            jij%atom(a1)%pair(i)%r=lo_chop(p%fractional_to_cartesian(r),lo_sqtol)
            jij%cutoff=max(jij%cutoff,norm2(jij%atom(a1)%pair(i)%r))
        enddo
    enddo
    ! and a tiny tolerance to the cutoff.
    jij%cutoff=jij%cutoff+lo_sqtol

end subroutine

end module
