#include "precompilerdefinitions"
submodule (type_crystalstructure) type_crystalstructure_symmetry
implicit none
contains

!> Classify the structure in different ways. Will erase any classification that might have existed earlier.
module subroutine classify(p,how,uc,tolerance,refine,timereversal)
    !> crystal structure
    class(lo_crystalstructure), intent(inout) :: p
    !> classify how?
    character(len=*), intent(in) :: how
    !> perhaps a unitcell
    type(lo_crystalstructure), intent(in), optional :: uc
    !> non-default tolerance?
    real(flyt), intent(in), optional :: tolerance
    !> while classifying, refine the cell?
    logical, intent(in), optional :: refine
    !> with the symmetry stuff, consider time-reveral symmetry as well?
    logical, intent(in), optional :: timereversal

    real(flyt) :: tol
    logical :: refinecell

    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_tol
    endif

    if ( present(refine) ) then
        refinecell=refine
    else
        refinecell=.false.
    endif

    select case(trim(how))
    case('bravais')
        ! figure out the name of the Bravais lattice
        call find_bravais_lattice(p,tolerance=tol,refine=refinecell)
    case('bz')
        ! calculate the Brillouin zone
        if ( p%info%havebravais .eqv. .false. ) call find_bravais_lattice(p)
        call get_brillouin_zone(p)
    case('spacegroup')
        ! Get the spacegroup
        if ( present(timereversal) .eqv. .false. ) then
            call lo_stop_gracefully(['Time reversal symmetry needs to be specified for finding the space group!'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        call p%sym%generate(p%latticevectors,timereversal,p%r,p%species,verbosity=p%info%verbosity,tolerance=tol)
        p%info%decidedtimereversal=.true.
        p%info%havespacegroup=.true.
    case('wedge')
        if ( present(timereversal) .eqv. .false. ) then
            call lo_stop_gracefully(['Time reversal symmetry needs to be specified defining the irreducible wedge.'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Get the irreducible wedge. Need quite some stuff for this.
        if ( p%info%havebravais .eqv. .false. ) call find_bravais_lattice(p)
        if ( p%info%havebz .eqv. .false. ) call get_brillouin_zone(p)
        if ( p%info%havespacegroup .eqv. .false. ) then
            call p%sym%generate(p%latticevectors,timereversal,p%r,p%species,verbosity=p%info%verbosity)
            p%info%decidedtimereversal=.true.
            p%info%havespacegroup=.true.
        endif
        if ( p%info%havewedge .eqv. .false. ) then
            call get_irreducible_wedge(p,timereversal)
            call label_highsymmetry_points(p,timereversal)
        endif
    case('supercell')
        ! Figure out how a supercell and unitcell might be related.
        call classify_unitcell_supercell_pair(uc,p)
    case('irrep')
        ! Get the character table for the spacegroup
        if ( present(timereversal) .eqv. .false. ) then
            call lo_stop_gracefully(['Time reversal symmetry needs to be specified to get the irrep'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( p%info%havebravais .eqv. .false. ) call find_bravais_lattice(p)
        if ( p%info%havebz .eqv. .false. ) call get_brillouin_zone(p)
        if ( p%info%havespacegroup .eqv. .false. ) then
            call p%sym%generate(p%latticevectors,timereversal,p%r,p%species,verbosity=p%info%verbosity)
            p%info%decidedtimereversal=.true.
            p%info%havespacegroup=.true.
        endif
        if ( p%info%havewedge .eqv. .false. ) then
            call get_irreducible_wedge(p,timereversal)
            call label_highsymmetry_points(p,timereversal)
        endif
        call p%sym%get_character_table(p%info%verbosity)
    case default
        call lo_stop_gracefully(['Not a known way to classify things.'],lo_exitcode_param,__FILE__,__LINE__)
    end select
end subroutine

!> Returns the Cartesian coordinates for a high symmetry point from its label
module function coordinate_from_high_symmetry_point_label(p,label,previous) result(qpoint)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the label
    character(len=*), intent(in) :: label
    !> The previous point
    real(flyt), dimension(3), intent(in), optional :: previous
    !> the q-point
    real(flyt), dimension(3) :: qpoint
    !
    integer :: i,j
    real(flyt) :: f0,f1

    ! First make sure that all things are initialized properly.
    if ( p%info%havewedge .eqv. .false. ) then
        call lo_stop_gracefully(['Need the wedge to fetch coordinates from labels.'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%info%pointslabelled .eqv. .false. ) then
        call lo_stop_gracefully(['Points need to be labelled before they can be fetched'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    if ( present(previous) ) then
        ! find the closest point that matches
        j=0
        f0=lo_huge
        do i=1,p%bz%nhighsymmetrypoints
            if ( trim(label) .eq. trim(p%bz%label(i)) ) then
                ! found a point matching.
                f1=norm2(previous-p%bz%highsymmetrypoints(:,i))
                if ( f1 .lt. f0 ) then
                    ! I want the closest matching point.
                    j=j+1
                    qpoint=p%bz%highsymmetrypoints(:,i)
                    f0=f1
                endif
            endif
        enddo
        !
        if ( j .eq. 0 ) then
            ! I failed
            call lo_stop_gracefully(['Could not find a tabulated value for the point '//trim(label)],lo_exitcode_param,__FILE__,__LINE__)
        else
            ! I'm done
            return
        endif
    else
        ! just get the first point that matches
        do i=1,p% irrw%nnodes
            if ( trim(label) .eq. trim(p%irrw%label(i)) ) then
                ! found the point
                qpoint=p%irrw%r(:,i)
                return
            endif
        enddo
        do i=1,p%bz%nhighsymmetrypoints
            if ( trim(label) .eq. trim(p%bz%label(i)) ) then
                qpoint=p%bz%highsymmetrypoints(:,i)
                return
            endif
        enddo
    endif
    ! if I made it here, I failed
    call lo_stop_gracefully(['Could not find a tabulated value for the point '//trim(label)],lo_exitcode_param,__FILE__,__LINE__)
end function

! Below are not exposed outside this submodule.

!> figure out how a unitcell and a supercell is related
subroutine classify_unitcell_supercell_pair(uc,ss)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss

    real(flyt), dimension(3,3) :: ssm
    integer :: multiple

    ! Maybe I already did this? Does not matter, do it again, it's fast anyway.
    if ( allocated(ss%info%index_in_unitcell) ) deallocate(ss%info%index_in_unitcell)
    if ( allocated(ss%info%cellindex) ) deallocate(ss%info%cellindex)

    ! Start by assuming it is not a supercell
    ss%info%supercell=.false.

    if ( ss%info%verbosity .gt. 0 ) then
        write(*,*) ' '
        write(*,*) 'Classifying crystal structure: matching unit and supercell '
    endif

    ! A series of sanity checks, to make really sure they actually are a unit-supercell pair
    checklatticevectors: block
        real(flyt) :: f0
        integer :: i,j
        logical :: dl

        ! First that the volume of the supercell is an integer multiple of the unitcell
        f0=ss%volume/uc%volume
        if ( abs(f0-anint(f0)) .gt. lo_tol ) then
            call lo_stop_gracefully(['Supercell/unitcell volume ratio is non-integer.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! The number of atoms must also have the same multiple as the volume
        multiple=int(anint(ss%volume/uc%volume))
        if ( uc%na*multiple .ne. ss%na ) then
            call lo_stop_gracefully(['Inconsistent ratio of atoms between supercell and unitcell.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! Get the supercell latticevectors in terms of the unitcell latticevectors, these
        ! must also have integer values
        ssm=matmul(uc%inv_latticevectors,ss%latticevectors)
        dl=.false.
        do j=1,3
        do i=1,3
            if ( abs( ssm(i,j)-anint(ssm(i,j)) ) .gt. lo_tol ) dl=.true.
        enddo
        enddo
        if ( dl ) then
            call lo_stop_gracefully(['Supercell matrix with non-integer values.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        ! Clean it a little
        ssm=lo_chop(anint(ssm),lo_sqtol)
    end block checklatticevectors

    if ( ss%info%verbosity .gt. 0 ) then
        write(*,*) '... passed the basic sanity checks'
    endif

    ! So at this point I am pretty sure it's a valid supercell.
    ! Now try to match it per atom
    matchatoms: block
        real(flyt), dimension(:,:), allocatable :: ucv,ssv
        integer, dimension(:,:), allocatable :: ssi
        integer, dimension(:), allocatable :: ucctr,ssctr,ssuci
        integer :: i,j
        logical :: dl
        !
        lo_allocate(ucv(5,uc%na))
        lo_allocate(ssv(5,ss%na))
        lo_allocate(ssi(3,ss%na))
        lo_allocate(ucctr(uc%na))
        lo_allocate(ssctr(ss%na))
        lo_allocate(ssuci(ss%na))

        ! Build some vectors that are good for comparisons!
        do i=1,uc%na
            ucv(1:3,i)=lo_chop(lo_clean_fractional_coordinates(uc%r(:,i),lo_sqtol),lo_sqtol)
            ucv(4,i)=uc%atomic_number(i)
            if ( allocated(uc%flavor) ) then
                ucv(5,i)=uc%flavor(i)
            else
                ucv(5,i)=0
            endif
        enddo

        ! Same for the supercell, but with some more stuff
        do i=1,ss%na
            ssi(:,i)=floor(matmul(ssm,ss%r(:,i))+lo_sqtol)
            ssv(1:3,i)=lo_chop( lo_clean_fractional_coordinates( matmul(ssm,ss%r(:,i))-ssi(:,i)*1.0_flyt, lo_sqtol) ,lo_sqtol)
            ssv(4,i)=ss%atomic_number(i)
            if ( allocated(ss%flavor) ) then
                ssv(5,i)=uc%flavor(i)
            else
                ssv(5,i)=0
            endif
        enddo

        ! Now match it
        ucctr=0
        ssctr=0
        ssuci=0
        do i=1,ss%na
            do j=1,uc%na
                if ( sum(abs(ucv(:,j)-ssv(:,i))) .lt. lo_tol ) then
                    ucctr(j)=ucctr(j)+1
                    ssuci(i)=j
                    ssctr(i)=ssctr(i)+1
                endif
            enddo
        enddo

        ! Plenty of sanity checks
        dl=.false.
        do i=1,ss%na
            if ( ssctr(i) .eq. 0 ) then
                write(*,*) 'Atom ',tochar(i),' in the supercell was not matched'
                dl=.true.
            elseif ( ssctr(i) .gt. 1 ) then
                write(*,*) 'Atom ',tochar(i),' in the supercell was matched ',tochar(ssctr(i)),' times'
                dl=.true.
            endif
        enddo
        do i=1,uc%na
            if ( ucctr(i) .ne. multiple ) then
                write(*,*) 'Atom ',tochar(i),' in the unitcell was matched ',tochar(ucctr(i)),' times, I expected ',tochar(multiple)
                dl=.true.
            endif
        enddo
        if ( dl ) stop
        if ( ss%info%verbosity .gt. 0 ) then
            write(*,*) '... matched all the atoms'
        endif

        ! If I made it here, everything is ok I hope.
        ss%info%supercell=.true.
        ss%info%supercellmatrix=int(anint(ssm))
        lo_allocate(ss%info%cellindex(3,ss%na))
        lo_allocate(ss%info%index_in_unitcell(ss%na))
        ss%info%cellindex=ssi+1
        ss%info%index_in_unitcell=ssuci

        ! And some cleanup
        lo_deallocate(ucv)
        lo_deallocate(ssv)
        lo_deallocate(ssi)
        lo_deallocate(ucctr)
        lo_deallocate(ssctr)
        lo_deallocate(ssuci)
    end block matchatoms

    if ( ss%info%verbosity .gt. 0 ) then
        write(*,*) 'Successfully matched unit and supercell'
    endif
end subroutine

!> Generates the Brillouin zone. I just build the Voronoi diagram of the reciprocal lattice.
subroutine get_brillouin_zone(p)
    !> The crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !
    type(lo_voronoi_diagram) :: voro
    real(flyt), dimension(3,1) :: dumpts
    real(flyt), dimension(3) :: v0
    real(flyt) :: rc
    integer :: i,j,l

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) ' '
        write(*,*) 'Generating Brillouin zone.'
    endif

    ! get the Voronoi diagram for the reciprocal lattice
    rc=lo_bounding_sphere_of_box(p%reciprocal_latticevectors)*3.0_flyt ! just a little safety margin
    dumpts=0.0_flyt

    call voro%generate(dumpts,p%reciprocal_latticevectors,cutoff=rc,verbosity=p%info%verbosity)
    ! stuff this into the BZ structure
    p%bz%nnodes=voro%cell(1)%nnodes
    p%bz%nfaces=voro%cell(1)%nfaces
    p%bz%nedges=voro%cell(1)%nedges
    lo_allocate(p%bz%node( p%bz%nnodes ))
    lo_allocate(p%bz%face( p%bz%nfaces ))
    lo_allocate(p%bz%edge( p%bz%nedges ))
    lo_allocate(p%bz%r( 3,p%bz%nnodes ))

    do i=1,p%bz%nfaces
        p%bz%face(i)%n=voro%cell(1)%face(i)%n
        lo_allocate(p%bz%face(i)%ind( p%bz%face(i)%n ))
    enddo
    p%bz%r=voro%cell(1)%r
    p%bz%node=voro%cell(1)%node
    p%bz%face=voro%cell(1)%face
    p%bz%edge=voro%cell(1)%edge
    p%bz%rmin=voro%cell(1)%rmin
    p%bz%rmax=voro%cell(1)%rmax
    ! I'm not 100% sure that this assignment is valid, but it seems to work.
    ! It should according to the f2008 standard, if I got it right.
    p%bz%polyhedron=voro%cell(1)%polyhedron

    ! How many high-symmetry points are there in total? (+1 is gamma)
    p%bz%nhighsymmetrypoints=p%bz%nnodes+p%bz%nedges+p%bz%nfaces+1
    ! Get a list of all high-symmetry points
    lo_allocate(p%bz%highsymmetrypoints(3,p%bz%nhighsymmetrypoints))
    ! add all nodes
    l=0
    do i=1,p%bz%nnodes
        l=l+1
        p%bz%highsymmetrypoints(:,l)=p%bz%r(:,i)
    enddo
    ! add centers of all edges
    do i=1,p%bz%nedges
        l=l+1
        p%bz%highsymmetrypoints(:,l)=(p%bz%r( :, p%bz%edge(i)%i1)+p%bz%r( :, p%bz%edge(i)%i2 ))*0.5_flyt
    enddo
    ! and centers of all faces
    do i=1,p%bz%nfaces
        l=l+1
        v0=0.0_flyt
        do j=1,p%bz%face(i)%n
            v0=v0+p%bz%r(:, p%bz%face(i)%ind(j) )
        enddo
        v0=v0/(p%bz%face(i)%n*1.0_flyt)
        p%bz%highsymmetrypoints(:,l)=v0
    enddo
    ! and gamma
    p%bz%highsymmetrypoints(:,l+1)=0.0_flyt

    ! and tell the world we did it
    p%info%havebz=.true.

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) 'Got Brillouin zone with '//tochar(p%bz%nnodes)//' nodes, '//&
                   tochar(p%bz%nfaces)//' faces, '//tochar(p%bz%nedges)//' edges and '//&
                   tochar(p%bz%nhighsymmetrypoints)//' high symmetry points.'
    endif
end subroutine

!> Will return the irreducible wedge of the BZ. This might be unnecesarily complicated, but it seems to work quite reliably.
subroutine get_irreducible_wedge(p,timereversal)
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> consider timereversal symmetry
    logical, intent(in) :: timereversal
    !
    real(flyt), dimension(:,:), allocatable :: tetr,untetr
    real(flyt) :: t0
    integer, dimension(:,:), allocatable :: teti,unteti
    integer :: ntet,nuntet

    if ( p%info%verbosity .gt. 0 ) then
        t0=walltime()
        write(*,*) ''
        write(*,*) 'Generating the irreducible wedge'
    endif

    ! Check the prereqs
    if ( p%info%havebravais .eqv. .false. ) then
        call lo_stop_gracefully(['Need the Bravais lattice to get the wedge'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%info%havespacegroup .eqv. .false. ) then
        call lo_stop_gracefully(['Need the spacegroup to get the wedge'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%info%havebz .eqv. .false. ) then
        call lo_stop_gracefully(['Need the BZ to get the wedge'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%sym%timereversal .neqv. timereversal ) then
        call lo_stop_gracefully(['Conflicting information regarding time-reversal symmetry'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! First thing to do is to chop the BZ into tetrahedrons
    buildtetrahedrons: block
        real(flyt), dimension(:,:), allocatable :: points
        real(flyt), dimension(3) :: facecentroid
        real(flyt) :: f0
        integer, dimension(:), allocatable :: di
        real(flyt), dimension(3,4) :: t1,t2
        real(flyt), dimension(3) :: v0
        integer :: np,i_gamma,tetcount,i_1,i_2,i_centroid
        integer :: i,j,k,l,ii,jj,kk

        ! fetch the points
        np=p%bz%nhighsymmetrypoints
        lo_allocate(tetr(3,np))
        tetr=p%bz%highsymmetrypoints
        if ( p%info%verbosity .gt. 0 ) write(*,*) '... found '//tochar(np)//' high symmetry points.'

        ! count number of tetrahedrons
        ntet=0
        do i=1,p%bz%nfaces
            ntet=ntet+p%bz%face(i)%n*2
        enddo
        lo_allocate(teti(4,ntet))

        ! need to keep track of where gamma is
        i_gamma=np

        ! Construct the tetrahedrons
        tetcount=0
        do i=1,p%bz%nfaces
            ! Get the points on this face that are along the edge: the nodes and midpoints of edges.
            lo_allocate(points( 3, p%bz%face(i)%n*2 ))
            l=0
            facecentroid=0.0_flyt
            do j=1,p%bz%face(i)%n
                ! add the point
                l=l+1
                points(:,l)=p%bz%r( : , p%bz%face(i)%ind(j) )
                ! and the one on the middle of the path
                l=l+1
                k=lo_index_in_periodic_array(j+1,p%bz%face(i)%n)
                points(:,l)=( p%bz%r( : , p%bz%face(i)%ind(j) ) + p%bz%r( : , p%bz%face(i)%ind(k) ) )*0.5_flyt
                ! and add to the centroid of the face
                facecentroid=facecentroid+p%bz%r(:, p%bz%face(i)%ind(j) )/(p%bz%face(i)%n*1.0_flyt)
            enddo
            ! Anglesort the points. Should not be necessary, but you never know.
            call p%bz%face(i)%plane%anglesort(points)

            ! Find these points in the list of tetrahedron nodes.
            lo_allocate(di( p%bz%face(i)%n*2 ))
            i_centroid=0
            do j=1,np
                ! find the centroid
                v0=facecentroid-tetr(:,j)
                if ( lo_sqnorm(v0) .lt. lo_sqtol ) i_centroid=j
                ! find the others
                do k=1,size(points,2)
                    v0=tetr(:,j)-points(:,k)
                    if ( lo_sqnorm(v0) .lt. lo_sqtol ) di(k)=j
                enddo
            enddo

            ! Build some tetrahedrons, walk along the perimeter of the face
            do j=1,size(di,1)
                tetcount=tetcount+1
                ! get the two points not gamma and the centroid
                i_1=di( j )
                i_2=di( lo_index_in_periodic_array(j+1, size(di,1) ) )
                ! store the tetrahedron
                teti(:,tetcount)=[i_gamma,i_centroid,i_1,i_2]
            enddo

            lo_deallocate(points)
            lo_deallocate(di)
        enddo
        ! sanity checks, first check the number of tetrahedrons
        if ( tetcount .ne. ntet ) then
            call lo_stop_gracefully(['Bad number of tetrahedrons, expected '//tochar(ntet)//' found '//tochar(tetcount)],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        ! check that they make up the entire BZ
        f0=0.0_flyt
        do i=1,ntet
            t1=tetr(:,teti(:,i))
            f0=f0+lo_unsigned_tetrahedron_volume(t1)
        enddo
        if ( abs(f0-1.0_flyt/p%volume) .gt. lo_tol ) then
            call lo_stop_gracefully(['Volume of tetrahedrons does not add up to the volume of the BZ'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        ! They can not be degenerate
        do i=1,ntet
        do j=i+1,ntet
            t1=tetr(:,teti(:,i))
            t2=tetr(:,teti(:,j))
            ! To avoid annoying things with the Cray compiler I had to inline this function.
            kk=0
            do ii=1,4
                f0=lo_huge
                do jj=1,4
                    v0=t1(:,ii)-t2(:,jj)
                    f0=min(f0,lo_sqnorm(v0))
                enddo
                if ( f0 .lt. lo_sqtol ) kk=kk+1
            enddo
            if ( kk .eq. 4 ) then
                call lo_stop_gracefully(['Tetrahedron '//tochar(i)//' and '//tochar(j)//' are degenerate'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        enddo
        ! The only thing left to check if if there is overlap, but that has never happened yet, so I don't bother.
    end block buildtetrahedrons
    if ( p%info%verbosity .gt. 0 ) write(*,*) '... chopped the BZ into '//tochar(ntet)//' tetrahedrons'

    ! Reduce the tetrahedrons by symmetry
    symmtet: block
        real(flyt), dimension(:,:), allocatable :: tetcentroids,dr1
        real(flyt), dimension(3,6) :: strangepts
        real(flyt), dimension(3,4) :: t1,t2,t3
        real(flyt), dimension(3) :: v0
        real(flyt) :: f0,f1,bsrad
        integer, dimension(:), allocatable :: di1,di2,tetclass,untet
        integer :: i,j,l,o,pt,ii,jj,kk

        lo_allocate(di1(ntet))
        di1=1
        ! Start by findin the unique tetrahedrons
        tetloop1: do i=1,ntet
            t1=tetr(:,teti(:,i))
            do j=i+1,ntet
                do o=1,p%sym%n
                    t2=tetr(:,teti(:,j))
                    ! Rotate the tetrahedron, this used to be a function but the Cray compilers got angry
                    do ii=1,4
                        t2(:,ii)=lo_operate_on_vector(p%sym%op(o),t2(:,ii),reciprocal=.true.)
                    enddo
                    ! Inlined comparison function
                    kk=0
                    do ii=1,4
                        f0=lo_huge
                        do jj=1,4
                            v0=t1(:,ii)-t2(:,jj)
                            f0=min(f0,lo_sqnorm(v0))
                        enddo
                        if ( f0 .lt. lo_sqtol ) kk=kk+1
                    enddo
                    if ( kk .eq. 4 ) then
                        di1(i)=0
                        cycle tetloop1
                    endif
                    if ( timereversal ) then
                        ! Inlined comparison function
                        t3=-t2
                        kk=0
                        do ii=1,4
                            f0=lo_huge
                            do jj=1,4
                                v0=t1(:,ii)-t3(:,jj)
                                f0=min(f0,lo_sqnorm(v0))
                            enddo
                            if ( f0 .lt. lo_sqtol ) kk=kk+1
                        enddo
                        if ( kk .eq. 4 ) then
                            di1(i)=0
                            cycle tetloop1
                        endif
                    endif
                enddo
            enddo
        enddo tetloop1
        ! Store the unique tetrahedrons
        nuntet=sum(di1)
        lo_allocate(untet(nuntet))
        l=0
        do i=1,ntet
            if ( di1(i) .eq. 1 ) then
                l=l+1
                untet(l)=i
            endif
        enddo
        lo_deallocate(di1)
        ! Figure out what kind each tetrahedron is
        lo_allocate(tetclass(ntet))
        tetclass=0
        tetloop2: do i=1,ntet
            t1=tetr(:,teti(:,i))
            do j=1,nuntet
                do o=1,p%sym%n
                    t2=tetr(:,teti(:,untet(j)))
                    do ii=1,4
                        t2(:,ii)=lo_operate_on_vector(p%sym%op(o),t2(:,ii),reciprocal=.true.)
                    enddo
                    t3=-t2

                    ! First comparison, inlined because Cray
                    kk=0
                    do ii=1,4
                        f0=lo_huge
                        do jj=1,4
                            v0=t1(:,ii)-t2(:,jj)
                            f0=min(f0,lo_sqnorm(v0))
                        enddo
                        if ( f0 .lt. lo_sqtol ) kk=kk+1
                    enddo
                    if ( kk .eq. 4 ) then
                        tetclass(i)=j
                        cycle tetloop2
                    endif
                    if ( timereversal ) then
                        ! Second comparison, inlined because Cray
                        kk=0
                        do ii=1,4
                            f0=lo_huge
                            do jj=1,4
                                v0=t1(:,ii)-t3(:,jj)
                                f0=min(f0,lo_sqnorm(v0))
                            enddo
                            if ( f0 .lt. lo_sqtol ) kk=kk+1
                        enddo
                        if ( kk .eq. 4 ) then
                            tetclass(i)=j
                            cycle tetloop2
                        endif
                    endif
                enddo
            enddo
            if ( tetclass(i) .eq. 0 ) then
                write(*,*) 'Failed classifying tetrahedron',i
                stop
            endif
        enddo tetloop2
        if ( p%info%verbosity .gt. 0 ) write(*,*) '... reduced to '//tochar(nuntet)//' tetrahedrons'

        ! The wedge is built by picking one of each of the irreducible tetrahedrons. There are many many ways to
        ! do this, but I only need one, and one that generates a polyhedron that is as convex as possible.
        ! My algorithm seems dumber than it is, but it actually has not failed (so far).

        ! The centroids of each tetrahedron
        lo_allocate(tetcentroids(3,ntet))
        do i=1,ntet
            do j=1,3
                tetcentroids(j,i)=lo_mean( tetr(j,teti(:,i)) )
            enddo
        enddo
        ! some random points
        strangepts(:,1)=[1010, 123  ,34  ]*1.0_flyt
        strangepts(:,2)=[10,   1023 ,34  ]*1.0_flyt
        strangepts(:,3)=[10,   123  ,1034]*1.0_flyt
        strangepts(:,4)=[1010, 1023 ,34  ]*1.0_flyt
        strangepts(:,5)=[1010, 123  ,1034]*1.0_flyt
        strangepts(:,6)=[1110, 1023 ,1036]*1.0_flyt

        ! So, the algorithm is that I pick one random point, and then the irreducible tetrahedrons closest
        ! to this point. Then I measure the bounding sphere of the centroids of the irreducible. The smallest
        ! bounding sphere is the most compact representation of the wedge, at last I hope so.
        lo_allocate(di1(nuntet))
        lo_allocate(di2(nuntet))
        di1=0
        di2=0
        bsrad=lo_huge
        do pt=1,6
            ! For each unique tetrahedron, pick the one closest to the strange point
            do i=1,nuntet
                f0=lo_huge
                do j=1,ntet
                    if ( tetclass(j) .ne. i ) cycle
                    v0=tetcentroids(:,j)-strangepts(:,pt)
                    f1=lo_sqnorm(v0)
                    if ( f1 .lt. f0 ) then
                        di1(i)=j
                        f0=f1
                    endif
                enddo
            enddo
            ! Now calculate the bounding sphere
            f0=0.0_flyt
            do i=1,nuntet
            do j=i+1,nuntet
                v0=tetcentroids(:,di1(i))-tetcentroids(:,di1(j))
                f0=max(f0,lo_sqnorm( v0 ))
            enddo
            enddo
            if ( f0+1E-8_flyt .lt. bsrad ) then
                bsrad=f0
                di2=di1
            endif
        enddo

        ! There is no need to keep all the points anymore, the irreducible tetrahedra are all we need
        ! First reduce the points to the irreducible
        lo_allocate(dr1(3,nuntet*4))
        l=0
        do i=1,nuntet
        do j=1,4
            l=l+1
            dr1(:,l)=tetr(:,teti(j,di2(i)))
        enddo
        enddo
        call lo_return_unique(dr1,untetr,lo_tol)
        ! and fix the tetrahedron indices
        lo_allocate(unteti(4,nuntet))
        unteti=0
        do i=1,nuntet
        do j=1,4
            do l=1,size(untetr,2)
                v0=untetr(:,l)-tetr(:, teti(j,di2(i)) )
                if ( lo_sqnorm(v0) .lt. lo_sqtol ) then
                    unteti(j,i)=l
                    cycle
                endif
            enddo
            if ( unteti(j,i) .eq. 0 ) then
                call lo_stop_gracefully(['Failed finding irreducible tetrahedron index'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        enddo
        ! And now I should have the wedge as a reduced list of tetrahedrons
    end block symmtet

    ! With the short list of tetrahedrons, build the actual wedge!
    buildwedge: block
        type(lo_plane) :: plane
        real(flyt), dimension(:,:), allocatable :: dum,dumun,points
        real(flyt), dimension(3,3) :: planepts
        real(flyt), dimension(3) :: v0,v1,v2
        integer, dimension(:), allocatable :: di
        integer :: i,j,k,l,ii,jj

        ! First, store the points
        p%irrw%nnodes=size(untetr,2)
        lo_allocate(p%irrw%r(3,p%irrw%nnodes))
        p%irrw%r=untetr

        ! Now I want to find the planes that enclose these points. This is faster than I thought it would be.
        ! I just construct all possible planes from the unique points, and keep those planes where all the
        ! unique points are on the negative side, or on the plane.
        lo_allocate(dum(4,p%irrw%nnodes**3))
        l=0
        do i=1,p%irrw%nnodes
        do j=1,p%irrw%nnodes
        do k=1,p%irrw%nnodes
            if ( i .eq. j ) cycle
            if ( j .eq. k ) cycle
            if ( k .eq. i ) cycle
            planepts(:,1)=p%irrw%r(:,i)
            planepts(:,2)=p%irrw%r(:,j)
            planepts(:,3)=p%irrw%r(:,k)
            ! first a dummy check that the points are not coplanaer:
            v0=planepts(:,1)-planepts(:,3)
            v1=planepts(:,2)-planepts(:,3)
            v2=lo_cross(v0,v1)
            if ( lo_sqnorm( v2 ) .lt. lo_sqtol ) cycle
            ! construct this plane
            call plane%generate(planepts)
            ! check that all points are on the negative side, or on the plane
            jj=0
            do ii=1,p%irrw%nnodes
                if ( plane%distance_to_point( p%irrw%r(:,ii) ) .lt. lo_tol ) jj=jj+1
            enddo
            ! if all points are below, it is a keeper!
            if ( jj .eq. p%irrw%nnodes ) then
                l=l+1
                dum(1:3,l)=plane%normal
                dum(4,l)=plane%p
            endif
        enddo
        enddo
        enddo

        ! it's enough to have the unique planes
        call lo_return_unique(dum(:,1:l),dumun)
        lo_deallocate(dum)

        ! Now I know the number of faces!
        p%irrw%nfaces=size(dumun,2)
        lo_allocate(p%irrw%face( p%irrw%nfaces ))

        lo_allocate(points(3,p%irrw%nnodes))
        ! Just have to reconstruct the planes again, from the Hessian normal form. With that
        ! I build the polygons that enclose the irreducible wedge.
        do i=1,p%irrw%nfaces
            ! Define the enclosing planes
            v0=dumun(1:3,i)
            v1=dumun(1:3,i)*(-dumun(4,i))
            call p%irrw%face(i)%plane%generate( normal=v0, point=v1 )
            ! count points on this plane
            points=0.0_flyt
            l=0
            do j=1,p%irrw%nnodes
                if ( abs(p%irrw%face(i)%plane%distance_to_point( p%irrw%r(:,j) )) .lt. lo_tol ) then
                    l=l+1
                    points(:,l)=p%irrw%r(:,j)
                endif
            enddo
            ! find the centroid
            v0=0.0_flyt
            do j=1,l
                v0=v0+points(:,j)/l
            enddo
            ! See if any of the points is the centroid
            k=0
            do j=1,l
                if ( lo_sqnorm(v0-points(:,j)) .lt. lo_sqtol ) k=k+1
            enddo
            ! the centroid should not be added.
            l=l-k

            ! make some space
            if ( l .eq. 0 ) then
                write(*,*) 'Warning, bad plane in the irreducible wedge. Something is strange.'
                ! I don't have to kill it since this really only matters for plots and stuff.
                ! but maybe I should. Dunno.
                p%irrw%face(i)%n=0
                cycle
            endif

            ! Store them
            p%irrw%face(i)%n=l
            lo_allocate(p%irrw%face(i)%ind(l))
            lo_allocate(di(l))
            lo_allocate(dum(3,l))
            l=0
            do j=1,p%irrw%nnodes
                if ( abs(p%irrw%face(i)%plane%distance_to_point( p%irrw%r(:,j) )) .lt. lo_tol ) then
                if ( lo_sqnorm(v0-p%irrw%r(:,j)) .gt. lo_sqtol ) then
                    l=l+1
                    p%irrw%face(i)%ind(l)=j
                    dum(:,l)=p%irrw%r(:,j)
                endif
                endif
            enddo
            ! Anglesort to make plots nice
            call p%irrw%face(i)%plane%anglesort(dum,di)
            p%irrw%face(i)%ind=p%irrw%face(i)%ind( di )
            lo_deallocate(di)
            lo_deallocate(dum)
        enddo
    end block buildwedge
    ! And I have the irreducible wedge!
    p%info%havewedge=.true.
    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) 'Generated the irreducible wedge in '//tochar(walltime()-t0)//'s'
    endif
end subroutine

!> Figures out the Bravais family of the lattice and specific member of that family. Converts the basis to an a,b,c,alpha,beta,gamma representation. From this the Bravais family is deduced. Using that information, a standard orientation of that lattice is generated to figure out the space group and the nomenclature of the high symmetry points. Not stable at all.
subroutine find_bravais_lattice(p,tolerance,refine)
    !> The crystal structure
    class(lo_crystalstructure), intent(inout) :: p
    !> Tolerance, if not default
    real(flyt), intent(in), optional :: tolerance
    !> Should the lattice vectors be refined to match their Bravais family?
    logical, intent(in), optional :: refine
    !
    real(flyt) :: tol,breakpoint
    integer :: i
    !
    real(flyt), dimension(3,3) :: pri,con
    real(flyt) :: f0
    character(len=10) :: family,familymember,oldfamily,newfamily
    character(len=50) :: familymemberlongname
    logical :: refinevectors

    ! default tolerance or something else?
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_tol
    endif

    ! Should I adjust the lattice vectors to match the Bravais lattice, to high precision?
    if ( present(refine) ) then
        refinevectors=refine
    else
        refinevectors=.false.
    endif

    ! Get the family and stuff
    call find_bravais_family(p%latticevectors,tol,family,familymember,familymemberlongname)


    ! Perhaps refine it properly?
    if ( refinevectors ) then

        ! First check where the break point is, as in what tolerance changes the determined Bravais lattice
        f0=1E-4_flyt ! Very rough tolerance to start with
        call find_bravais_family(p%latticevectors,f0,oldfamily,familymember,familymemberlongname)
        breakpoint=0.0_flyt
        do i=1,50
            f0=f0*0.5_flyt
            call find_bravais_family(p%latticevectors,f0,newfamily,familymember,familymemberlongname)
            if ( oldfamily .ne. newfamily ) then
                breakpoint=f0
                exit
            else
                oldfamily=newfamily
            endif
        enddo

        ! Set the tolerance above the breakpoint, I guess
        ! Not sure what to do about the break point. Could be useful for some heuristics in the future.
        if ( breakpoint .gt. lo_tiny ) then
            if ( p%info%verbosity .gt. 0 ) write(*,*) '... Bravais lattice changes from ',trim(oldfamily),' to ',trim(newfamily), ' at tolerance:',breakpoint
            f0=breakpoint*2
        else
            f0=lo_tol
        endif

        ! Temporarily get the canonical basis, and the transformation to it
        call get_conventional_abc_and_standard_basis(p%latticevectors,family,pri,con,p%info%unitcell_lattice_parameter)
        p%info%unique_primitive_basis=lo_chop(pri,lo_sqtol)
        call find_transformation_between_current_and_conventional(p)
        pri=matmul(pri,p%info%permutation_to_unique)

        ! Write the old lattice vectors
        if ( p%info%verbosity .gt. 0 ) then
             write(*,*) '... input lattice vectors (scaled with lattice parameter):'
             do i=1,3
                write(*,*) p%latticevectors(:,i)/p%info%unitcell_lattice_parameter
             enddo
        endif

        ! Now adjust the latticevectors to match this exactly
        call adjust_latticevectors(p%latticevectors,pri)

        ! Reclassify with the adjusted lattice vectors
        call find_bravais_family(p%latticevectors,f0,newfamily,familymember,familymemberlongname)
        call get_conventional_abc_and_standard_basis(p%latticevectors,family,pri,con,p%info%unitcell_lattice_parameter)
        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) '... slightly refined lattice vectors:'
            do i=1,3
                write(*,*) p%latticevectors(:,i)/p%info%unitcell_lattice_parameter
            enddo
        endif
        ! Add a final touch, to get e.g. sqrt(3)/2 to all digits
        p%latticevectors=lo_chop(p%latticevectors/p%info%unitcell_lattice_parameter,1E-8_flyt)*p%info%unitcell_lattice_parameter

        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) '... completely refined lattice vectors:'
            do i=1,3
                write(*,*) lo_chop(p%latticevectors(:,i)/p%info%unitcell_lattice_parameter,1E-13_flyt)
            enddo
        endif

        ! Now update all the coordinates and things that depend on the latticevectors
        p%inv_latticevectors=lo_chop(lo_invert3x3matrix(p%latticevectors),1E-13_flyt)
        p%reciprocal_latticevectors=lo_chop(lo_reciprocal_basis(p%latticevectors),1E-13_flyt)
        p%inv_reciprocal_latticevectors=lo_chop(lo_invert3x3matrix(p%reciprocal_latticevectors),1E-13_flyt)
    endif

    ! It's neat to store the standardized conventional and primitive basis
    call get_conventional_abc_and_standard_basis(p%latticevectors,family,pri,con,p%info%unitcell_lattice_parameter)
    p%info%bravaislattice=trim(adjustl(familymember))
    p%info%bravaislatticelongname=trim(adjustl(familymemberlongname))
    p%info%unique_primitive_basis=lo_chop(pri,lo_sqtol)
    p%info%unique_conventional_basis=lo_chop(con,lo_sqtol)
    ! It's also convenient to have the transformation from the current format to the
    ! standardized one. It is an annoying transformation+permutation.
    call find_transformation_between_current_and_conventional(p)
    ! Now I have the Bravais lattice
    p%info%havebravais=.true.
    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) 'Found Bravais lattice: ',trim(p%info%bravaislatticelongname),' (',trim(p%info%bravaislattice),')'
    endif
end subroutine

!> Make sure the distances and angles make sense to really high precision
subroutine adjust_latticevectors(basis,reference)
    real(flyt), dimension(3,3), intent(inout) :: basis
    real(flyt), dimension(3,3), intent(in) :: reference
    !
    real(flyt), dimension(3,3) :: grad,m1
    real(flyt) :: step,f0,f1
    integer :: iter,swctr

    ! Set some starting parameters
    m1=basis
    f0=pardist(basis,reference)
    grad=gradient(basis,reference,basis)
    step=1E-12_flyt/lo_frobnorm(grad)
    swctr=1
    ! Actual minimization
    do iter=1,10000
        if ( f0 .lt. 1E-14_flyt ) exit
        ! update the basis
        m1=lo_chop(m1-step*grad,1E-14_flyt)
        ! new distance
        f1=pardist(m1,reference)
        ! Might be satisified here
        if ( f1 .lt. f0 ) then
            step=step*1.1_flyt
        else
            m1=lo_chop(m1+step*grad,1E-14_flyt)
            swctr=swctr+1
            step=1E-12_flyt/lo_frobnorm(grad)
            step=step/sqrt(swctr*1.0_flyt)
            grad=gradient(m1,reference,basis)
        endif
        f0=f1
    enddo
    ! And make sure the volume does not change
    f0=lo_determ(basis)
    f1=lo_determ(m1)
    basis=m1*f0/f1
    contains

    ! Calculate the gradient of the distance difference with a 9-point stencil
    function gradient(m0,m1,m2) result(grad)
        real(flyt), dimension(3,3) :: m0,m1,m2,grad
        !
        integer :: i,j,k
        real(flyt) :: delta
        real(flyt), dimension(3,3) :: dm
        real(flyt), dimension(8) :: sc_wt,sc_dl,sc_val
        !
        delta=1E-8_flyt
        !
        sc_wt=[3,-32,168,-672,672,-168,32,-3]/(840.0_flyt*delta)
        sc_dl=[-4,-3,-2,-1,1,2,3,4]*delta
        !
        do i=1,3
        do j=1,3
            ! Only provide a gradient for non-zero entries in the original
            if ( abs(m2(j,i)) .gt. 1E-10_flyt ) then
                do k=1,8
                    dm=0.0_flyt
                    dm(j,i)=sc_dl(k)
                    sc_val(k)=pardist(m0+dm,m1)
                enddo
                grad(j,i)=sum(sc_wt*sc_val)
            else
                grad(j,i)=0.0_flyt
            endif
        enddo
        enddo
    end function

    !> distance in parameters such as a,b,c,alpha,beta,gamma
    function pardist(m0,m1) result(d)
        real(flyt), dimension(3,3), intent(in) :: m0,m1
        real(flyt) :: d
        !
        real(flyt), dimension(6) :: p0,p1

        ! Calculate the two a,b,c,alpha,beta,gamma
        call lo_get_axis_angles(m0,p0(1),p0(2),p0(3),p0(4),p0(5),p0(6))
        call lo_get_axis_angles(m1,p1(1),p1(2),p1(3),p1(4),p1(5),p1(6))
        ! Squared difference
        p0=(p0-p1)**2
        d=sqrt(sum(p0))
    end function
end subroutine

!> Find the transformation that takes the current definition of the lattice to a standardized one.
subroutine find_transformation_between_current_and_conventional(p)
    type(lo_crystalstructure), intent(inout) :: p

    integer :: i
    real(flyt), dimension(3,3) :: m0,m1
    real(flyt), dimension(3,3,6) :: permutationmatrices
    real(flyt), dimension(6) :: y0,y1
    real(flyt) :: a,b,c,al,be,gm,f0,f1

    ! Premutation matrices
    permutationmatrices(1,:,1)=[1,0,0]*1.0_flyt
    permutationmatrices(2,:,1)=[0,1,0]*1.0_flyt
    permutationmatrices(3,:,1)=[0,0,1]*1.0_flyt

    permutationmatrices(1,:,2)=[1,0,0]*1.0_flyt
    permutationmatrices(2,:,2)=[0,0,1]*1.0_flyt
    permutationmatrices(3,:,2)=[0,1,0]*1.0_flyt

    permutationmatrices(1,:,3)=[0,1,0]*1.0_flyt
    permutationmatrices(2,:,3)=[1,0,0]*1.0_flyt
    permutationmatrices(3,:,3)=[0,0,1]*1.0_flyt

    permutationmatrices(1,:,4)=[0,0,1]*1.0_flyt
    permutationmatrices(2,:,4)=[1,0,0]*1.0_flyt
    permutationmatrices(3,:,4)=[0,1,0]*1.0_flyt

    permutationmatrices(1,:,5)=[0,1,0]*1.0_flyt
    permutationmatrices(2,:,5)=[0,0,1]*1.0_flyt
    permutationmatrices(3,:,5)=[1,0,0]*1.0_flyt

    permutationmatrices(1,:,6)=[0,0,1]*1.0_flyt
    permutationmatrices(2,:,6)=[0,1,0]*1.0_flyt
    permutationmatrices(3,:,6)=[1,0,0]*1.0_flyt

    ! Find the correct permutation.
    m0=p%info%unique_primitive_basis
    call lo_get_axis_angles(m0,a,b,c,al,be,gm)
    y0=[a,b,c,al,be,gm]
    f0=lo_huge
    do i=1,6
        m1=matmul(p%latticevectors,permutationmatrices(:,:,i))
        call lo_get_axis_angles(m1,a,b,c,al,be,gm)
        y1=[a,b,c,al,be,gm]
        f1=sum(abs(y0-y1))
        if ( f1 .lt. f0 ) then
            p%info%permutation_to_unique=permutationmatrices(:,:,i)
            f0=f1
        endif
    enddo
    ! And the transformation
    p%info%transformation_to_unique=matmul(p%inv_latticevectors,p%info%unique_primitive_basis)
end subroutine

!> Given three latticevectors and a tolerance, figure out the Bravais family.
subroutine find_bravais_family(latticevectors,tolerance,family,member,memberlongname)
    !> latticevectors
    real(flyt), dimension(3,3), intent(in) :: latticevectors
    !> the tolerance, in A
    real(flyt), intent(in) :: tolerance
    !> Bravais family
    character(len=10) :: family
    !> Specific member in that family
    character(len=10) :: member
    !> Long name for this familymemeber
    character(len=50) :: memberlongname
    !
    real(flyt), dimension(3,3) :: m0
    real(flyt), dimension(3,3) :: pri,con
    real(flyt) :: a,b,c,al,be,gm
    real(flyt) :: ka,kb,kc,kal,kbe,kgm
    real(flyt) :: ap,bp,cp
    real(flyt) :: a0,b0,c0,al0,be0,gm0
    real(flyt) :: d90,d60,d120,dbcc
    real(flyt) :: f0,f1
    real(flyt) :: tol_dist,tol_angle

    ! some angles to test against
    d90=1.570796326794897_flyt ! 90 degrees in radians
    d60=1.047197551196598_flyt ! 60 degrees in radians
    d120=2.094395102393195_flyt ! 120 degrees in radians
    dbcc=1.910633236249019_flyt ! the bcc angle

    ! Scale the distance and angle tolerances together somehow
    tol_dist=tolerance
    tol_angle=tolerance*180.0_flyt/lo_pi

    ! Get a,b,c,alpha,beta,gamma
    call get_a_lt_b_lt_c_order(latticevectors,m0)
    call lo_get_axis_angles(m0,a0,b0,c0,al0,be0,gm0)

    ! Guess the family first. The trick was to move these around to catch things in the correct order.
    if ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 1 CUB Simple cubic 1.000000 1.000000 1.000000 ang 90.000000 90.000000 90.000000
        family='CUB'
        member='CUB'
        memberlongname='simple cubic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d60,tol_angle) ) then
        ! 2 FCC Face centered cubic 0.707107 0.707107 0.707107 ang 60.000000 60.000000 60.000000
        family='FCC'
        member='FCC'
        memberlongname='face centered cubic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,dbcc,tol_angle) ) then
        ! 3 BCC Body centered cubic 0.866025 0.866025 0.866025 ang 109.471221 109.471221 109.471221
        family='BCC'
        member='BCC'
        memberlongname=''
    elseif ( lo_two_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 4 TET Tetragonal 1.000000 1.000000 1.100000 ang 90.000000 90.000000 90.000000
        family='TET'
        member='TET'
        memberlongname='tetragonal'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_two_equal(al0,be0,gm0,tol_angle)) then ! .and. lo_three_larger_than_num(a0,be0,gm0,d90,tol_angle) ) then
        ! 5 BCT1 Body-centered tetragonal 0.838153 0.838153 0.838153 ang 106.753588 106.753588 115.054967
        family='BCT'
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(con,a,b,c,al,be,gm)
        if ( a .gt. c ) then
            member='BCT1'
        else
            member='BCT2'
        endif
        memberlongname='body-centered tetragonal'
    elseif ( lo_none_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 7 ORC Orthorhombic 1.000000 1.100000 1.200000 ang 90.000000 90.000000 90.000000
        family='ORC'
        member='ORC'
        memberlongname='orthorhombic'
    elseif ( lo_one_equal_to_num(al0,be0,gm0,d120,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) .and. lo_two_equal(a0,b0,c0,tol_dist) ) then
        ! 13 HEX Hexagonal 1.000000 1.000000 1.200000 ang 90.000000 90.000000 120.000000
        family='HEX'
        member='HEX'
        memberlongname='hexagonal'
    elseif ( lo_one_equal_to_num(al0,be0,gm0,d60,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) .and. lo_two_equal(a0,b0,c0,tol_dist) ) then
        ! 13 HEX Hexagonal 1.000000 1.000000 1.200000 ang 90.000000 90.000000 60.000000
        family='HEX2'
        member='HEX2'
        memberlongname='hexagonal'
    elseif ( lo_two_equal(a0,b0,c0,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 12 ORCC C-centered orthorhombic 0.743303 0.743303 1.300000 ang 90.000000 90.000000 95.452622
        family='ORCC'
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(pri,a,b,c,al,be,gm)
        if ( gm .lt. lo_pi*0.5_flyt ) then
            member='ORCC1'
        else
            member='ORCC2'
        endif
        memberlongname='c-centered orthorhombic'
    elseif ( lo_none_equal(a0,b0,c0,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 16 MCL Monoclinic 1.000000 1.100000 1.200000 ang 80.000000 90.000000 90.000000
        family='MCL'
        member='MCL'
        memberlongname='monoclinic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_none_equal(al0,be0,gm0,tol_angle) ) then
        ! 11 ORCI Body-centered orthorhombic 0.955249 0.955249 0.955249 ang 116.875594 109.693369 102.178552
        family='ORCI'
        member='ORCI'
        memberlongname='body-centered orthorhombic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal(al0,be0,gm0,tol_angle) ) then
        ! 14 RHL1 Rhombohedral 1.000000 1.000000 1.000000 ang 80.000000 80.000000 80.000000
        family='RHL'
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(con,a,b,c,al,be,gm)
        if ( al .lt. d90 ) then
            member='RHL1'
        else
            member='RHL2'
        endif
        memberlongname='rhombohedral'
    elseif ( lo_two_equal(a0,b0,c0,tol_dist) .and. lo_two_equal(al0,be0,gm0,tol_angle) ) then
        ! 17 MCLC1 C-centered monoclinic 0.743303 0.743303 1.200000 ang 82.617700 82.617700 84.547378
        family='MCLC'
        ! this one is annoying to differentiate
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(con,a,b,c,al,be,gm)
        ! get the reciprocal basis as well
        m0=lo_reciprocal_basis(pri)
        call lo_get_axis_angles(m0,ka,kb,kc,kal,kbe,kgm)
        ! sort stuff out
        if ( abs(kgm-d90) .lt. tol_angle ) then
            member='MCLC2'
        elseif ( kgm .gt. d90 ) then
            member='MCLC1'
        else
            f0=b*cos(al)/c+(b*sin(al)/a)**2
            if ( abs(f0-1.0_flyt) .lt. tol_dist ) then
                member='MCLC4'
            elseif ( f0 .lt. 1.0_flyt ) then
                member='MCLC3'
            else
                member='MCLC5'
            endif
        endif
        memberlongname='c-centered monoclinic'
    else
        ! final annoying one, differentiate ORCF and TRI.
        call get_c_lt_b_lt_a_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,al,be,gm)
        ! get the ORCF conventional lattice parameters
        a=2*Sqrt(bp*cp*Cos(al))
        b=2*Sqrt(ap*cp*Cos(be))
        c=2*Sqrt(ap*bp*Cos(gm))
        ! build the primitive again
        m0(:,1)=(/0.0_flyt, b/2, c/2/)
        m0(:,2)=(/a/2, 0.0_flyt, c/2/)
        m0(:,3)=(/a/2, b/2, 0.0_flyt/)
        ! this should not match for a trigonal lattice
        if ( abs(ap-norm2(m0(:,1))) .lt. tol_dist .and. &
             abs(bp-norm2(m0(:,2))) .lt. tol_dist .and. &
             abs(cp-norm2(m0(:,3))) .lt. tol_dist ) then
            ! 8 ORCF1 Face-centered orthorhombic 1.250000 1.118034 0.901388 ang 75.636697 60.050918 44.312385
            family='ORCF'
            ! check stuff
            call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
            call lo_get_axis_angles(con,a,b,c,al,be,gm)
            f0=1.0_flyt/(a**2)
            f1=1.0_flyt/(b**2)+1.0_flyt/(c**2)
            if ( abs(f0-f1) .lt. tol_dist ) then
                member='ORCF3'
            elseif ( f0 .gt. f1 ) then
                member='ORCF1'
            else
                member='ORCF2'
            endif
            memberlongname='face-centered orthorhombic'
        else
            ! 22 TRI1a Triclinic 1.000000 1.100000 1.200000 ang 60.000000 70.000000 75.000000
            family='TRI'
            call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
            call lo_get_axis_angles(con,a,b,c,al,be,gm)
            m0=lo_reciprocal_basis(con)
            call lo_get_axis_angles(m0,ka,kb,kc,kal,kbe,kgm)
            if ( lo_three_smaller_than_num(kal,kbe,kgm,d90,tol_angle) ) then
                member='TRI1b'
            elseif ( lo_three_larger_than_num(kal,kbe,kgm,d90,tol_angle) ) then
                member='TRI1a'
            elseif ( lo_two_smaller_than_num(kal,kbe,kgm,d90,tol_angle) ) then
                member='TRI2b'
            else
                member='TRI2a'
            endif
            memberlongname='triclinic'
        endif
    endif
end subroutine

!> Given the Bravais family, I can figure out the lattice parameter of the conventional cell, as well as a standardized primitive basis.
subroutine get_conventional_abc_and_standard_basis(latticevectors,family,pri,con,alat)
    !> the lattice vectors
    real(flyt), dimension(3,3), intent(in) :: latticevectors
    !> Bravais family name
    character(len=*), intent(in) :: family
    !> standardized primitive basis
    real(flyt), dimension(3,3), intent(out) :: pri
    !> standardized conventional basis
    real(flyt), dimension(3,3), intent(out) :: con
    !> the "a" lattice parameter, in some standard sense
    real(flyt), intent(out), optional :: alat
    !
    real(flyt), dimension(3,3) :: m0
    real(flyt) :: a,b,c,ap,bp,cp,ra,rb,rg

    select case(trim(family))
    case('CUB')
    ! CUB  a=b=c, al=be=gm=90
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a,0.0_flyt,0.0_flyt]
        pri(:,2)=[0.0_flyt,a,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,a]
        con=pri
    case('FCC')
    ! FCC  a=b=c, al=be=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*ap
        b = sqrt(2.0_flyt)*bp
        c = sqrt(2.0_flyt)*cp
        pri(:,1)=[0.0_flyt,a/2,a/2]
        pri(:,2)=[a/2,0.0_flyt,a/2]
        pri(:,3)=[a/2,a/2,0.0_flyt]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,a,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,a]
    case('BCC')
    ! BCC  a=b=c, al=be=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = (2*ap)/sqrt(3.0_flyt)
        b = (2*bp)/sqrt(3.0_flyt)
        c = (2*cp)/sqrt(3.0_flyt)
        pri(:,1)=[-a/2,a/2,a/2]
        pri(:,2)=[a/2,-a/2,a/2]
        pri(:,3)=[a/2,a/2,-a/2]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,a,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,a]
    case('TET')
    ! TET  a=b, al=be=gm=90
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a,0.0_flyt,0.0_flyt]
        pri(:,2)=[0.0_flyt,a,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('BCT')
    ! BCT  a=b=c, al=be
        call get_al_eq_be_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*ap*Sqrt(-Cos(ra) - Cos(rg))
        b = sqrt(2.0_flyt)*bp*Sqrt(-Cos(ra) - Cos(rg))
        c = 2*cp*Sqrt(-Cos(rb))
        pri(:,1)=[-a/2,a/2,c/2]
        pri(:,2)=[a/2,-a/2,c/2]
        pri(:,3)=[a/2,a/2,-c/2]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,a,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('ORC')
    ! ORC  a<b<c al=be=gm=90
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a,0.0_flyt,0.0_flyt]
        pri(:,2)=[0.0_flyt,b,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('ORCF')
    ! ORCF  c<b<a al<be<gm
        call get_c_lt_b_lt_a_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = 2*Sqrt(bp*cp*Cos(ra))
        b = 2*Sqrt(ap*cp*Cos(rb))
        c = 2*Sqrt(ap*bp*Cos(rg))
        pri(:,1)=[0.0_flyt,b/2,c/2]
        pri(:,2)=[a/2,0.0_flyt,c/2]
        pri(:,3)=[a/2,b/2,0.0_flyt]
        !
        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,b,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('ORCI')
    ! ORCI  a=b=c gm<be<al
        call get_gm_lt_be_lt_al_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*Sqrt(ap**2*(-Cos(rb) - Cos(rg)))
        b = sqrt(2.0_flyt)*Sqrt(ap**2*(-Cos(ra) - Cos(rg)))
        c = sqrt(2.0_flyt)*Sqrt(ap**2*(-Cos(ra) - Cos(rb)))
        pri(:,1)=[-a/2,b/2,c/2]
        pri(:,2)=[a/2,-b/2,c/2]
        pri(:,3)=[a/2,b/2,-c/2]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,b,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('ORCC')
        call get_al_eq_be_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*Sqrt(dot_product(m0(1,:),m0(1,:)) + dot_product(m0(1,:),m0(2,:)))
        b = sqrt(2.0_flyt)*Sqrt(-dot_product(m0(1,:),m0(2,:)) + dot_product(m0(2,:),m0(2,:)))
        c = cp

        pri(:,1)=[a/2,b/2,0.0_flyt]
        pri(:,2)=[a/2,-b/2,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,b,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('HEX')
    ! HEX  a=b al=be
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[ a/2,-(sqrt(3.0_flyt)*a)/2 ,0.0_flyt]
        pri(:,2)=[ a/2, (sqrt(3.0_flyt)*a)/2 ,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('HEX2')
    ! HEX  a=b al=be
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[ a, 0.0_flyt ,0.0_flyt]
        pri(:,2)=[ a/2, (sqrt(3.0_flyt)*a)/2 ,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('RHL')
    ! RHL  a=b=c al=be=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a*Cos(ra/2), -(a*Sin(ra/2)), 0.0_flyt ]
        pri(:,2)=[a*Cos(ra/2),  a*Sin(ra/2),   0.0_flyt ]
        pri(:,3)=[a*Cos(ra)*(1/Cos(ra/2)), 0.0_flyt, a*Sqrt(1 - Cos(ra)**2*(1/Cos(ra/2))**2) ]
        con=pri
    case('MCL')
    ! MCL  a<b<c be=gm
        call get_mcl_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a, 0.0_flyt, 0.0_flyt]
        pri(:,2)=[0.0_flyt, b, 0.0_flyt]
        pri(:,3)=[0.0_flyt, c*Cos(ra), c*Sin(ra)]
        con=pri
    case('MCLC')
    ! MCLC  a=b al=be
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = norm2(m0(:,1)-m0(:,2))
        b = norm2(m0(:,1)+m0(:,2))
        c = cp
        ra = acos(dot_product(m0(:,1)+m0(:,2),m0(:,3))/(b*c))
        pri(:,1)=[a/2, b/2, 0.0_flyt]
        pri(:,2)=[-a/2, b/2, 0.0_flyt]
        pri(:,3)=[0.0_flyt, c*Cos(ra), c*Sin(ra)]
        !
        con(:,1)=[a, 0.0_flyt, 0.0_flyt]
        con(:,2)=[0.0_flyt, b, 0.0_flyt]
        con(:,3)=[0.0_flyt, c*Cos(ra), c*Sin(ra)]
    case('TRI')
    ! TRI  a<b<c, al/=be/=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp

        pri(:,1)=[a, 0.0_flyt, 0.0_flyt]
        pri(:,2)=[b*Cos(rg), b*Sin(rg), 0.0_flyt]
        pri(:,3)=[c*Cos(rb), c*(Cos(ra) - Cos(rb)*Cos(rg))*(1/Sin(rg)),&
                 c*(1/Sin(rg))*Sqrt(-Cos(ra)**2 - Cos(rb)**2 + 2*Cos(ra)*Cos(rb)*Cos(rg) + Sin(rg)**2)]
        con=pri
    end select
    ! keep the "a" lattice parameter
    if ( present(alat) ) then
        alat=a
    endif
end subroutine

!> Sort the lattice vectors so that a<b<c
subroutine get_a_lt_b_lt_c_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    real(flyt), dimension(3) :: abc
    integer, dimension(3) :: ind
    !
    call lo_get_axis_angles(original_basis(:,:),a,b,c,al,be,gm)
    abc(1)=a
    abc(2)=b
    abc(3)=c
    call qsort(abc,ind)
    reordered_basis=original_basis(:,ind)
end subroutine

!> sort the lattice vectors in MCL order
subroutine get_mcl_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm,f0,f1
    integer, dimension(:,:), allocatable :: perms
    integer :: i
    !
    reordered_basis=0.0_flyt
    call lo_permutations(perms,3)
    f0=lo_huge
    do i=1,size(perms,2)
        call lo_get_axis_angles(original_basis(:,perms(:,i)),a,b,c,al,be,gm)
        f1=abs(a-c)+abs(b-c)+abs(be-0.5_flyt*lo_pi)+abs(gm-0.5_flyt*lo_pi)
        if ( f1 .lt. f0 ) then
            reordered_basis=original_basis(:,perms(:,i))
            f0=f1
        endif
    enddo
end subroutine

!> sort the lattice vectors so that c<b<a
subroutine get_c_lt_b_lt_a_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    real(flyt), dimension(3) :: abc
    integer, dimension(3) :: ind
    !
    call lo_get_axis_angles(original_basis,a,b,c,al,be,gm)
    abc(1)=a
    abc(2)=b
    abc(3)=c
    call qsort(abc,ind)
    reordered_basis(:,1)=original_basis(:,ind(3))
    reordered_basis(:,2)=original_basis(:,ind(2))
    reordered_basis(:,3)=original_basis(:,ind(1))
end subroutine

!> sort the lattice vectors so that gamma<beta<allpha
subroutine get_gm_lt_be_lt_al_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    real(flyt), dimension(3) :: abc
    integer, dimension(3) :: ind
    !
    call lo_get_axis_angles(original_basis(:,:),a,b,c,al,be,gm)
    abc(1)=al
    abc(2)=be
    abc(3)=gm
    call qsort(abc,ind)
    reordered_basis(:,1)=original_basis(:,ind(3))
    reordered_basis(:,2)=original_basis(:,ind(2))
    reordered_basis(:,3)=original_basis(:,ind(1))
end subroutine

!> sort the lattice vectors so that a=b
subroutine get_a_eq_b_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    integer, dimension(:,:), allocatable :: perms
    integer :: i
    !
    reordered_basis=0.0_flyt
    call lo_permutations(perms,3)
    do i=1,size(perms,2)
        call lo_get_axis_angles(original_basis(:,perms(:,i)),a,b,c,al,be,gm)
        if ( abs(a-b) .lt. lo_tol ) then
            reordered_basis=original_basis(:,perms(:,i))
            exit
        endif
    enddo
end subroutine

!> sort the lattice vectors so that alpha=beta
subroutine get_al_eq_be_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    integer, dimension(:,:), allocatable :: perms
    integer :: i
    !
    reordered_basis=0.0_flyt
    call lo_permutations(perms,3)
    do i=1,size(perms,2)
        call lo_get_axis_angles(original_basis(:,perms(:,i)),a,b,c,al,be,gm)
        if ( abs(al-be) .lt. lo_radiantol ) then
            reordered_basis=original_basis(:,perms(:,i))
            exit
        endif
    enddo
end subroutine

!> Tedious routine to label the high symmetry points of the BZ. There is no logic to it whatsoever. I want to retire it and just use Miller indices instead, at least they are logical.
subroutine label_highsymmetry_points(p,timereversal)
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> consider time reversal symmetry in the labelling?
    logical, intent(in) :: timereversal

    type(lo_symset) :: sym
    real(flyt), dimension(:,:), allocatable :: pts
    real(flyt), dimension(:,:), allocatable :: r,rw
    real(flyt), dimension(3,3) :: bas,ibas
    real(flyt), dimension(3) :: v0
    real(flyt) :: zero,one,half,fourth,fiveovereight,threeoverfour,threeovereight,twooverthree,oneoverthree
    real(flyt) :: m,g,s,f,d,x,w,l,q
    real(flyt) :: a,b,c,al,be,gm
    integer :: i,j,o,np,nhit,nmiss
    character(len=10), dimension(:), allocatable :: lbl

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) ''
        write(*,*) 'Trying to label high-symmetry points'
    endif

    ! First block is to fetch the tabulated values
    zero=0.0_flyt
    one=1.0_flyt
    half=0.5_flyt
    fourth=0.25_flyt
    fiveovereight=5.0_flyt/8.0_flyt
    threeoverfour=3.0_flyt/4.0_flyt
    threeovereight=3.0_flyt/8.0_flyt
    twooverthree=2.0_flyt/3.0_flyt
    oneoverthree=1.0_flyt/3.0_flyt

    ! Structural parameters from the conventional cell in the correct order
    call lo_get_axis_angles(p%info%unique_conventional_basis,a,b,c,al,be,gm)
    select case(trim(p%info%bravaislattice))
        case('CUB ')
            ! Table 2 3 CUB
            lo_allocate(pts(3,3))
            lo_allocate(lbl(3))
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/   half,   half,   zero/)
            pts(:, 3)=(/   zero,   half,   zero/)
            lbl( 1)='R'
            lbl( 2)='M'
            lbl( 3)='X'
        case('FCC ')
            ! Table 3 5 FCC
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/fiveovereight, fourth,fiveovereight/)
            pts(:, 2)=(/threeovereight,threeovereight,threeoverfour/)
            pts(:, 3)=(/   half, fourth,threeoverfour/)
            pts(:, 4)=(/   half,   half,   half/)
            pts(:, 5)=(/   half,   zero,   half/)
            lbl( 1)='U'
            lbl( 2)='K'
            lbl( 3)='W'
            lbl( 4)='L'
            lbl( 5)='X'
        case('BCC ')
        ! Table 4 3 BCC
            lo_allocate(pts(3,3))
            lo_allocate(lbl(3))
            pts(:, 1)=(/ fourth, fourth, fourth/)
            pts(:, 2)=(/   half,  -half,   half/)
            pts(:, 3)=(/   zero,   zero,   half/)
            lbl( 1)='P'
            lbl( 2)='H'
            lbl( 3)='N'
        case('TET ')
        ! Table 5 5 TET
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/   zero,   half,   half/)
            pts(:, 2)=(/   half,   half,   half/)
            pts(:, 3)=(/   zero,   half,   zero/)
            pts(:, 4)=(/   half,   half,   zero/)
            pts(:, 5)=(/   zero,   zero,   half/)
            lbl( 1)='R'
            lbl( 2)='A'
            lbl( 3)='X'
            lbl( 4)='M'
            lbl( 5)='Z'
        case('BCT1 ')
        ! Table 6 6 BCT1
            lo_allocate(pts(3,6))
            lo_allocate(lbl(6))
            !#g=(1+c^2/a^2)/4
            g=(one+(c/a)**2)*0.25_flyt
            pts(:, 1)=(/   zero,   zero,   half/)
            pts(:, 2)=(/  -half,   half,   half/)
            pts(:, 3)=(/      g,      g,     -g/)
            pts(:, 4)=(/   zero,   half,   zero/)
            pts(:, 5)=(/     -g,  one-g,      g/)
            pts(:, 6)=(/ fourth, fourth, fourth/)
            lbl( 1)='X'
            lbl( 2)='M'
            lbl( 3)='Z'
            lbl( 4)='N'
            lbl( 5)='Z1'
            lbl( 6)='P'
        case('BCT2 ')
        ! Table 7 8 BCT2
            lo_allocate(pts(3,8))
            lo_allocate(lbl(8))
            !#g=(1+a^2/c^2)/4
            !#f=a^2/(2*c^2)
            g=(one+(a/c)**2)*0.25_flyt
            f=(a**2)/(2.0_flyt*c**2)
            pts(:, 1)=(/   zero,   zero,   half/)
            pts(:, 2)=(/   zero,   half,   zero/)
            pts(:, 3)=(/     -f,      f,   half/)
            pts(:, 4)=(/ fourth, fourth, fourth/)
            pts(:, 5)=(/   half,   half,      f/)
            pts(:, 6)=(/     -g,      g,      g/)
            pts(:, 7)=(/   half,   half,  -half/)
            pts(:, 8)=(/      g,  one-g,      g/)
            lbl( 1)='X'
            lbl( 2)='N'
            lbl( 3)='Y'
            lbl( 4)='P'
            lbl( 5)='Y1'
            lbl( 6)='R'
            lbl( 7)='Z'
            lbl( 8)='R1'
        case('ORC ')
        ! Table 8 7 ORC
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,   zero,   half/)
            pts(:, 2)=(/   half,   half,   half/)
            pts(:, 3)=(/   half,   zero,   zero/)
            pts(:, 4)=(/   half,   half,   zero/)
            pts(:, 5)=(/   zero,   half,   zero/)
            pts(:, 6)=(/   zero,   half,   half/)
            pts(:, 7)=(/   zero,   zero,   half/)
            lbl( 1)='U'
            lbl( 2)='R'
            lbl( 3)='X'
            lbl( 4)='S'
            lbl( 5)='Y'
            lbl( 6)='T'
            lbl( 7)='Z'
        case('ORCF1')
        ! Table 9 8 ORCF1 ORCF3
            lo_allocate(pts(3,8))
            lo_allocate(lbl(8))
            !#f=(1+a^2/b^2-a^2/c^2)/4
            !#g=(1+a^2/b^2+a^2/c^2)/4
            f=(one+a**2/b**2-a**2/c**2)*0.25_flyt
            g=(one+a**2/b**2+a**2/c**2)*0.25_flyt
            pts(:, 1)=(/   zero,      g,      g/)
            pts(:, 2)=(/   half, half+f,      f/)
            pts(:, 3)=(/    one,  one-g,  one-g/)
            pts(:, 4)=(/   half, half-f,  one-f/)
            pts(:, 5)=(/   half,   zero,   half/)
            pts(:, 6)=(/   half,   half,   half/)
            pts(:, 7)=(/   half,   half,   zero/)
            pts(:, 8)=(/    one,   half,   half/)
            lbl( 1)='X'
            lbl( 2)='A'
            lbl( 3)='X1'
            lbl( 4)='A1'
            lbl( 5)='Y'
            lbl( 6)='L'
            lbl( 7)='Z'
            lbl( 8)='T'
        case('ORCF3')
        ! Table 9 8 ORCF1 ORCF3
            lo_allocate(pts(3,8))
            lo_allocate(lbl(8))
            !#f=(1+a^2/b^2-a^2/c^2)/4
            !#g=(1+a^2/b^2+a^2/c^2)/4
            f=(one+a**2/b**2-a**2/c**2)*0.25_flyt
            g=(one+a**2/b**2+a**2/c**2)*0.25_flyt
            pts(:, 1)=(/   zero,      g,      g/)
            pts(:, 2)=(/   half, half+f,      f/)
            pts(:, 3)=(/    one,  one-g,  one-g/)
            pts(:, 4)=(/   half, half-f,  one-f/)
            pts(:, 5)=(/   half,   zero,   half/)
            pts(:, 6)=(/   half,   half,   half/)
            pts(:, 7)=(/   half,   half,   zero/)
            pts(:, 8)=(/    one,   half,   half/)
            lbl( 1)='X'
            lbl( 2)='A'
            lbl( 3)='X1'
            lbl( 4)='A1'
            lbl( 5)='Y'
            lbl( 6)='L'
            lbl( 7)='Z'
            lbl( 8)='T'
        case('ORCF2 ')
        ! Table 10 10 ORCF2
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            !#g=(1+a^2/b^2-a^2/c^2)/4
            !#d=(1+b^2/a^2-b^2/c^2)/4
            !#s=(1+c^2/b^2-c^2/a^2)/4
            g=(one+a**2/b**2-a**2/c**2)*0.25_flyt
            d=(one+b**2/a**2-b**2/c**2)*0.25_flyt
            s=(one+c**2/b**2-c**2/a**2)*0.25_flyt
            pts(:, 1)=(/  one-s, half-s,   half/)
            pts(:, 2)=(/   half, half-g,  one-g/)
            pts(:, 3)=(/      s, half+s,   half/)
            pts(:, 4)=(/   half, half+g,      g/)
            pts(:, 5)=(/   zero,   half,   half/)
            pts(:, 6)=(/ half-d,   half,  one-d/)
            pts(:, 7)=(/   half,   zero,   half/)
            pts(:, 8)=(/ half+d,   half,      d/)
            pts(:, 9)=(/   half,   half,   zero/)
            pts(:,10)=(/   half,   half,   half/)
            lbl( 1)='H'
            lbl( 2)='C'
            lbl( 3)='H1'
            lbl( 4)='C1'
            lbl( 5)='X'
            lbl( 6)='D'
            lbl( 7)='Y'
            lbl( 8)='D1'
            lbl( 9)='Z'
            lbl(10)='L'
        case('ORCI ')
        ! Table 11 12 ORCI
            lo_allocate(pts(3,12))
            lo_allocate(lbl(12))
            !#f=(1+a^2/c^2)/4
            !#g=(1+b^2/c^2)/4
            !#d=(b^2-a^2)/(4*c^2)
            !#l=(a^2+b^2)/(4*c^2)
            f=(one+a**2/c**2)*0.25_flyt
            g=(one+b**2/c**2)*0.25_flyt
            d=(b**2-a**2)/(4*c**2)
            l=(a**2+b**2)/(4*c**2)
            pts(:, 1)=(/ fourth, fourth, fourth/)
            pts(:, 2)=(/     -l,      l, half-d/)
            pts(:, 3)=(/     -f,      f,      f/)
            pts(:, 4)=(/      l,     -l, half+d/)
            pts(:, 5)=(/      f,  one-f,      f/)
            pts(:, 6)=(/ half-d, half+d,     -l/)
            pts(:, 7)=(/      g,     -g,      g/)
            pts(:, 8)=(/   zero,   half,   zero/)
            pts(:, 9)=(/  one-g,      g,     -g/)
            pts(:,10)=(/   half,   zero,   zero/)
            pts(:,11)=(/   half,   half,  -half/)
            pts(:,12)=(/   zero,   zero,   half/)
            lbl( 1)='W'
            lbl( 2)='L'
            lbl( 3)='X'
            lbl( 4)='L1'
            lbl( 5)='X1'
            lbl( 6)='L2'
            lbl( 7)='Y'
            lbl( 8)='R'
            lbl( 9)='Y1'
            lbl(10)='S'
            lbl(11)='Z'
            lbl(12)='T'
        case('ORCC1')
        ! Table 12 9 ORCC
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            f=(one+b**2/a**2)*0.25_flyt
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/      f,  one-f,   half/)
            pts(:, 3)=(/     -f,      f,   zero/)
            pts(:, 4)=(/     -f,      f,   half/)
            pts(:, 5)=(/      f,  one-f,   zero/)
            pts(:, 6)=(/   zero,   half,   half/)
            pts(:, 7)=(/  -half,   half,   zero/)
            pts(:, 8)=(/   zero,   half,   zero/)
            pts(:, 9)=(/   zero,   zero,   half/)
            pts(:,10)=(/   half,   half,   zero/)
            lbl( 1)='T'
            lbl( 2)='A'
            lbl( 3)='X'
            lbl( 4)='A1'
            lbl( 5)='X1'
            lbl( 6)='R'
            lbl( 7)='Y'
            lbl( 8)='S'
            lbl( 9)='Z'
            lbl(10)='Y'
        case('ORCC2')
            ! Table 12 9 ORCC
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            f=(one+a**2/b**2)*0.25_flyt
            pts(:, 1)=(/  -half,   half,   half/)
            pts(:, 2)=(/      f,      f,   half/)
            pts(:, 3)=(/      f,      f,   zero/)
            pts(:, 4)=(/     -f,  one-f,   half/)
            pts(:, 5)=(/     -f,  one-f,   zero/)
            pts(:, 6)=(/   zero,   half,   half/)
            pts(:, 7)=(/  -half,   half,   zero/)
            pts(:, 8)=(/   zero,   half,   zero/)
            pts(:, 9)=(/   zero,   zero,   half/)
            pts(:,10)=(/  -half,   half,   zero/)
            lbl( 1)='T'
            lbl( 2)='A'
            lbl( 3)='X'
            lbl( 4)='A1'
            lbl( 5)='X1'
            lbl( 6)='R'
            lbl( 7)='Y'
            lbl( 8)='S'
            lbl( 9)='Z'
            lbl(10)='Y'
        case('HEX ')
        ! Table 13 5 HEX
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/  oneoverthree,  oneoverthree,   zero/)
            pts(:, 2)=(/   zero,   zero,   half/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/ oneoverthree,  oneoverthree,   half/)
            pts(:, 5)=(/   half,   zero,   zero/)
            lbl( 1)='K'
            lbl( 2)='A'
            lbl( 3)='L'
            lbl( 4)='H'
            lbl( 5)='M'
        case('HEX2 ')
        ! Table 13 5 HEX
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/   zero,   zero,   half/)
            pts(:, 2)=(/  oneoverthree,  twooverthree,   half/)
            pts(:, 3)=(/  half,  half,  half/)
            pts(:, 4)=(/   half,   half,   zero/)
            pts(:, 5)=(/  oneoverthree,   twooverthree,   zero/)

            lbl( 1)='A'
            lbl( 2)='L'
            lbl( 3)='H'
            lbl( 4)='K'
            lbl( 5)='M'
        case('RHL1 ')
        ! Table 14 11 RHL1
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            !#g=(1+4*cos(alpha))/(2+4*cos(alpha))
            !#m=3/4-g/2
            g=(1.0_flyt+4.0_flyt*cos(al))/(2.0_flyt+4.0_flyt*cos(al))
            m=3.0_flyt/4.0_flyt-g/2.0_flyt
            pts(:, 1)=(/      g,      m,      m/)
            pts(:, 2)=(/      g,   half,  one-g/)
            pts(:, 3)=(/  one-m,  one-m,  one-g/)
            pts(:, 4)=(/   half,  one-g,  g-one/)
            pts(:, 5)=(/      m,      m,  g-one/)
            pts(:, 6)=(/   half,   half,   zero/)
            pts(:, 7)=(/  one-m,      m,   zero/)
            pts(:, 8)=(/   half,   zero,   zero/)
            pts(:, 9)=(/      m,   zero,     -m/)
            pts(:,10)=(/   half,   half,   half/)
            !pts(:,11)=(/   zero,   zero,  -half/)
            lbl( 1)='P'
            lbl( 2)='B'
            lbl( 3)='P1'
            lbl( 4)='B1'
            lbl( 5)='P2'
            lbl( 6)='F'
            lbl( 7)='Q'
            lbl( 8)='L'
            lbl( 9)='X'
            lbl(10)='Z'
            !lbl(11)='L1'
        case('RHL2 ')
        ! Table 15 7 RHL2
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            !#g=1/( 2*tan(alpha/2)^2 )
            !#m=3/4-g/2
            g=1.0_flyt/( 2*tan(al/2)**2 )
            m=3.0_flyt/4.0_flyt-g/2.0_flyt
            pts(:, 1)=(/      m,  m-one,  m-one/)
            pts(:, 2)=(/   half,  -half,   zero/)
            pts(:, 3)=(/      g,      g,      g/)
            pts(:, 4)=(/   half,   zero,   zero/)
            pts(:, 5)=(/  one-g,     -g,     -g/)
            pts(:, 6)=(/  one-m,     -m,  one-m/)
            pts(:, 7)=(/   half,  -half,   half/)
            lbl( 1)='P1'
            lbl( 2)='F'
            lbl( 3)='Q'
            lbl( 4)='L'
            lbl( 5)='Q1'
            lbl( 6)='P'
            lbl( 7)='Z'
        case('MCL ')
        ! Table 16 15 MCL
            lo_allocate(pts(3,15))
            lo_allocate(lbl(15))
            !#g=(1-b*cos(alpha)/c)/(2*sin(alpha)^2)
            !#m=1/2-g*c*cos(alpha)/b
            ! I think the curtarolo paper is a little wrong, had to modify this
            if ( al .lt. lo_pi*0.5_flyt ) then
                g=(1-b*cos(al)/c)/(2*sin(al)**2)
                m=0.5_flyt-g*c*cos(al)/b
               ! write(*,*) ' eta',g
               ! write(*,*) '  nu',m
            else
                !g=(1+c*cos(alpha)/b)/(2*sin(alpha)^2)
                !m=1/2+g*b*cos(alpha)/c
                m=(1.0_flyt+c*cos(al)/b)/(2*sin(al)*sin(al))
                g=0.5_flyt+m*b*cos(al)/c
                m=1.0_flyt-m
                !g=1.0_flyt-g
                !write(*,*) ' eta',g
                !write(*,*) '  nu',m
            endif
            pts(:, 1)=(/   zero,      g,     -m/)
            pts(:, 2)=(/   half,   half,   zero/)
            pts(:, 3)=(/   half,      g,  one-m/)
            pts(:, 4)=(/   zero,   half,   half/)
            pts(:, 5)=(/   half,  one-g,      m/)
            pts(:, 6)=(/   half,   zero,   half/)
            pts(:, 7)=(/   half,      g,     -m/)
            pts(:, 8)=(/   half,   zero,  -half/)
            pts(:, 9)=(/   zero,   half,   zero/)
            pts(:,10)=(/   half,   half,   half/)
            pts(:,11)=(/   zero,   zero,   half/)
            pts(:,12)=(/   zero,      g,  one-m/)
            pts(:,13)=(/   zero,   zero,  -half/)
            pts(:,14)=(/   zero,  one-g,      m/)
            pts(:,15)=(/   half,   zero,   zero/)
            lbl( 1)='H2'
            lbl( 2)='A'
            lbl( 3)='M'
            lbl( 4)='C'
            lbl( 5)='M1'
            lbl( 6)='D'
            lbl( 7)='M2'
            lbl( 8)='D1'
            lbl( 9)='X'
            lbl(10)='E'
            lbl(11)='Y'
            lbl(12)='H'
            lbl(13)='Y1:q'
            lbl(14)='H1'
            lbl(15)='Z'
        case('MCLC1')
        ! Table 17 16 MCLC1  MCLC2
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#f=(2-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g=1/2+2*f*c*cos(alpha)/b
            !#w=3/4-(a/(2*b*sin(alpha)))^2
            !#s=w+(3/4-w)*b*cos(alpha)/c
            f=(2.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g=0.5_flyt+2*f*c*cos(al)/b
            w=0.75_flyt-(a/(2*b*sin(al)))**2
            s=w+(0.75_flyt-w)*b*cos(al)/c
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/   half,   zero,   zero/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   zero,  -half,   zero/)
            pts(:, 5)=(/  one-w,  w-one,   zero/)
            pts(:, 6)=(/  one-f,  one-f,  one-g/)
            pts(:, 7)=(/      w,  one-w,   zero/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  w-one,     -w,   zero/)
            pts(:,10)=(/     -f,     -f,  one-g/)
            pts(:,11)=(/   half,   half,   zero/)
            pts(:,12)=(/  one-f,     -f,  one-g/)
            pts(:,13)=(/  -half,  -half,   zero/)
            pts(:,14)=(/      s,  one-s,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/  one-s,  s-one,   half/)
            lbl( 1)='L'
            lbl( 2)='N'
            lbl( 3)='M'
            lbl( 4)='N1'
            lbl( 5)='X'
            lbl( 6)='F'
            lbl( 7)='X1'
            lbl( 8)='F1'
            lbl( 9)='X2'
            lbl(10)='F2'
            lbl(11)='Y'
            lbl(12)='F3'
            lbl(13)='Y1'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='I1'
        case('MCLC2')
        ! Table 17 16 MCLC1  MCLC2
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#f=(2-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g=1/2+2*f*c*cos(alpha)/b
            !#w=3/4-(a/(2*b*sin(alpha)))^2
            !#s=w+(3/4-w)*b*cos(alpha)/c
            f=(2.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g=0.5_flyt+2*f*c*cos(al)/b
            w=0.75_flyt-(a/(2*b*sin(al)))**2
            s=w+(0.75_flyt-w)*b*cos(al)/c
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/   half,   zero,   zero/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   zero,  -half,   zero/)
            pts(:, 5)=(/  one-w,  w-one,   zero/)
            pts(:, 6)=(/  one-f,  one-f,  one-g/)
            pts(:, 7)=(/      w,  one-w,   zero/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  w-one,     -w,   zero/)
            pts(:,10)=(/     -f,     -f,  one-g/)
            pts(:,11)=(/   half,   half,   zero/)
            pts(:,12)=(/  one-f,     -f,  one-g/)
            pts(:,13)=(/  -half,  -half,   zero/)
            pts(:,14)=(/      s,  one-s,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/  one-s,  s-one,   half/)
            lbl( 1)='L'
            lbl( 2)='N'
            lbl( 3)='M'
            lbl( 4)='N1'
            lbl( 5)='X'
            lbl( 6)='F'
            lbl( 7)='X1'
            lbl( 8)='F1'
            lbl( 9)='X2'
            lbl(10)='F2'
            lbl(11)='Y'
            lbl(12)='F3'
            lbl(13)='Y1'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='I1'
        case('MCLC3')
        ! Table 18 16 MCLC3 MCLC4
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#l  =  (1 + b^2/a^2)/4
            !#d  =  b*c*cos(alpha)/(2*a^2)
            !#f  =  l - 1/4 + (1-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g  =  1/2 + 2*f*c*cos(alpha)/b
            !#s  =  1 + f - 2*l
            !#w  =  g - 2d
            l  =  (1.0_flyt + b**2/a**2)*0.25_flyt
            d  =  b*c*cos(al)/(2*a**2)
            f  =  l - 0.25_flyt + (1.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g  =  0.5_flyt + 2*f*c*cos(al)/b
            s  =  1.0_flyt + f - 2*l
            w  =  g - 2*d
            pts(:, 1)=(/   half,   zero,   zero/)
            pts(:, 2)=(/  one-s,  one-s,  one-w/)
            pts(:, 3)=(/   zero,   half,   zero/)
            pts(:, 4)=(/      s,      s,  one-w/)
            pts(:, 5)=(/   half,   half,   zero/)
            pts(:, 6)=(/  one-s,      s,  one-w/)
            pts(:, 7)=(/      l,      l,      d/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  one-l,     -l,      d/)
            pts(:,10)=(/  one-f,     -f,  one-g/)
            pts(:,11)=(/      l,      l,      d/)
            pts(:,12)=(/      f,      f,  one-g/)
            pts(:,13)=(/      l,  l-one,      d/)
            pts(:,14)=(/   half,   half,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/   half,   zero,   half/)
            lbl( 1)='N'
            lbl( 2)='F'
            lbl( 3)='N1'
            lbl( 4)='F1'
            lbl( 5)='X'
            lbl( 6)='F2'
            lbl( 7)='Y'
            lbl( 8)='H'
            lbl( 9)='Y1'
            lbl(10)='H1'
            lbl(11)='Y2'
            lbl(12)='H2'
            lbl(13)='Y3'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='M'
        case('MCLC4')
        ! Table 18 16 MCLC3 MCLC4
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#l  =  (1 + b^2/a^2)/4
            !#d  =  b*c*cos(alpha)/(2*a^2)
            !#f  =  l - 1/4 + (1-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g  =  1/2 + 2*f*c*cos(alpha)/b
            !#s  =  1 + f - 2*l
            !#w  =  g - 2d
            l  =  (1.0_flyt + b**2/a**2)*0.25_flyt
            d  =  b*c*cos(al)/(2*a**2)
            f  =  l - 0.25_flyt + (1.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g  =  0.5_flyt + 2*f*c*cos(al)/b
            s  =  1.0_flyt + f - 2*l
            w  =  g - 2*d
            pts(:, 1)=(/   half,   zero,   zero/)
            pts(:, 2)=(/  one-s,  one-s,  one-w/)
            pts(:, 3)=(/   zero,   half,   zero/)
            pts(:, 4)=(/      s,      s,  one-w/)
            pts(:, 5)=(/   half,   half,   zero/)
            pts(:, 6)=(/  one-s,      s,  one-w/)
            pts(:, 7)=(/      l,      l,      d/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  one-l,     -l,      d/)
            pts(:,10)=(/  one-f,     -f,  one-g/)
            pts(:,11)=(/      l,      l,      d/)
            pts(:,12)=(/      f,      f,  one-g/)
            pts(:,13)=(/      l,  l-one,      d/)
            pts(:,14)=(/   half,   half,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/   half,   zero,   half/)
            lbl( 1)='N'
            lbl( 2)='F'
            lbl( 3)='N1'
            lbl( 4)='F1'
            lbl( 5)='X'
            lbl( 6)='F2'
            lbl( 7)='Y'
            lbl( 8)='H'
            lbl( 9)='Y1'
            lbl(10)='H1'
            lbl(11)='Y2'
            lbl(12)='H2'
            lbl(13)='Y3'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='M'
        case('MCLC5 ')
        ! Table 19 16 MCLC5
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#f  = ( (b^2/a^2)  + (1-b*cos(alpha)/c)/(sin(alpha)^2 )/4
            !#l  =  g/2 + (b/(2*a))^2 - b*c*cos(alpha)/(2*a^2)
            !#x  =  (4*m - 1 - (b*sin(alpha)/a)^2) * c/(2*b*cos(alpha))
            !#g  =  1/2 + 2*f*c*cos(a)/b
            !#m  =  2*l - f
            !#d  =  f*c*cos(alpha)/b + x/2 - 1/4
            !#q  =  1 - f*(a/b)^2
            f  = ( (b**2/a**2)  + (1.0_flyt-b*cos(al)/c)/(sin(al)**2) )*0.25_flyt
            g  =  0.5_flyt + 2*f*c*cos(al)/b
            l  =  g*0.5_flyt + (b/(2*a))**2 - b*c*cos(al)/(2*a**2)
            m  =  2*l - f
            x  =  (4*m - 1.0_flyt - (b*sin(al)/a)**2)*c/(2*b*cos(al))
            d  =  f*c*cos(al)/b + x/2 - 0.25_flyt
            q  =  1.0_flyt - f*(a/b)**2
            !
            pts(:, 1)=(/   half,   zero,   half/)
            pts(:, 2)=(/      m,      m,      x/)
            pts(:, 3)=(/   half,   zero,   zero/)
            pts(:, 4)=(/  one-m,  one-m,  one-x/)
            pts(:, 5)=(/   zero,  -half,   zero/)
            pts(:, 6)=(/      m,  m-one,      x/)
            pts(:, 7)=(/   half,  -half,   zero/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/      l,      l,      d/)
            pts(:,10)=(/  one-f,     -f,  one-g/)
            pts(:,11)=(/  one-l,     -l,     -d/)
            pts(:,12)=(/     -f,     -f,  one-g/)
            pts(:,13)=(/     -l,     -l,     -d/)
            pts(:,14)=(/      q,  one-q,   half/)
            pts(:,15)=(/      l,  l-one,      d/)
            pts(:,16)=(/  one-q,  q-one,   half/)
            lbl( 1)='M'
            lbl( 2)='F'
            lbl( 3)='N'
            lbl( 4)='F1'
            lbl( 5)='N1'
            lbl( 6)='F2'
            lbl( 7)='X'
            lbl( 8)='H'
            lbl( 9)='Y'
            lbl(10)='H1'
            lbl(11)='Y1'
            lbl(12)='H2'
            lbl(13)='Y2'
            lbl(14)='I'
            lbl(15)='Y3'
            lbl(16)='I1'
        case('TRI1a')
        ! Table 20 7 TRI1a TRI2a
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,   half,   zero/)
            pts(:, 2)=(/   zero,   half,   half/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   half,   half,   half/)
            pts(:, 5)=(/   half,   zero,   zero/)
            pts(:, 6)=(/   zero,   half,   zero/)
            pts(:, 7)=(/   zero,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case('TRI2a')
        ! Table 20 7 TRI1a TRI2a
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,   half,   zero/)
            pts(:, 2)=(/   zero,   half,   half/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   half,   half,   half/)
            pts(:, 5)=(/   half,   zero,   zero/)
            pts(:, 6)=(/   zero,   half,   zero/)
            pts(:, 7)=(/   zero,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case('TRI1b')
        ! Table 21 7 TRI1b TRI2b
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,  -half,   zero/)
            pts(:, 2)=(/   zero,   zero,   half/)
            pts(:, 3)=(/  -half,  -half,   half/)
            pts(:, 4)=(/   zero,  -half,   half/)
            pts(:, 5)=(/   zero,  -half,   zero/)
            pts(:, 6)=(/   half,   zero,   zero/)
            pts(:, 7)=(/  -half,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case('TRI2b')
        ! Table 21 7 TRI1b TRI2b
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,  -half,   zero/)
            pts(:, 2)=(/   zero,   zero,   half/)
            pts(:, 3)=(/  -half,  -half,   half/)
            pts(:, 4)=(/   zero,  -half,   half/)
            pts(:, 5)=(/   zero,  -half,   zero/)
            pts(:, 6)=(/   half,   zero,   zero/)
            pts(:, 7)=(/  -half,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case default
            write(*,*) 'Unknown Bravais lattice type'
            lo_allocate(pts(3,1))
            lo_allocate(lbl(1))
            pts=0.0_flyt
            lbl(1)='GM'
    end select

    ! First try to label the full BZ, set the labels to nothing:
    if ( allocated(p%bz%label) ) deallocate(p%bz%label)
    lo_allocate(p%bz%label( p%bz%nhighsymmetrypoints ))
    do i=1,p%bz%nhighsymmetrypoints
        p%bz%label(i)='NP'
    enddo

    ! get the tabulated version of the reciprocal lattice
    bas=lo_reciprocal_basis(p%info%unique_primitive_basis)
    ibas=lo_invert3x3matrix( bas )

    ! fetch the coordinates for the bz nodes
    lo_allocate(r(3,p%bz%nhighsymmetrypoints))
    do i=1,p%bz%nhighsymmetrypoints
        ! to fractional
        r(:,i)=p%cartesian_to_fractional( p%bz%highsymmetrypoints(:,i), reciprocal=.true., pbc=.false. )
        ! then strange permutation
        r(:,i)=matmul(r(:,i),p%info%permutation_to_unique)
        ! to Cartesian again, in the new coordinate system
        r(:,i)=matmul(bas,r(:,i))
    enddo

    ! Also get the tabulated points to cartesian
    np=size(pts,2)
    do i=1,np
        pts(:,i)=matmul(bas,pts(:,i))
    enddo
    ! Get all possible symmetry operations in this represenation
    call sym%generate(bas,timereversal)
    if ( p%info%verbosity .gt. 0 ) write(*,*) '... generated symmorphic group'

    ! Try to label things:
    nhit=1
    nmiss=0
    ! minus 1 because I know that the last point is gamma
    bzptloop: do i=1,p%bz%nhighsymmetrypoints-1
        do o=1,sym%n
            ! rotate the BZ points
            v0=lo_operate_on_vector(sym%op(o),r(:,i),reciprocal=.true.)
            do j=1,np
                if ( lo_sqnorm(v0-pts(:,j)) .lt. lo_sqtol ) then
                    ! I found a point, I guess
                    nhit=nhit+1
                    p%bz%label(i)=trim(lbl(j))
                    cycle bzptloop
                endif
                if ( sym%timereversal ) then
                    if ( lo_sqnorm(v0+pts(:,j)) .lt. lo_sqtol ) then
                        ! I found a point, I guess
                        nhit=nhit+1
                        p%bz%label(i)=trim(lbl(j))
                        cycle bzptloop
                    endif
                endif
            enddo
        enddo
        ! if I made it here, I guess I did not find a label.
        nmiss=nmiss+1
    enddo bzptloop
    p%bz%label(p%bz%nhighsymmetrypoints)='GM'

    if ( p%info%verbosity .gt. 0 ) then
        !write(*,*) '... tried labelling: '//trim(int2char(nhit))//' hits, '//trim(int2char(nmiss))//' misses'
        write(*,*) '... tried labelling: ',nhit,' hits, ',nmiss,' misses'
        if ( nmiss .eq. 0 ) then
            write(*,*) '... highly successful, all points found'
        elseif ( nhit .eq. 0 ) then
            write(*,*) '... did not go very well, no points found'
        else
            write(*,*) '... moderate success, some points found'
        endif
    endif

    ! Test backward, see if all points got assigned
    nhit=0
    nmiss=0
    bwloop: do i=1,np
        do j=1,p%bz%nhighsymmetrypoints
            if ( trim(lbl(i)) .eq. trim(p%bz%label(j)) ) then
                nhit=nhit+1
                cycle bwloop
            endif
        enddo
        nmiss=nmiss+1
    enddo bwloop

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) '... testing assignment: '//tochar(nhit)//' hits, '//tochar(nmiss)//' misses'
        if ( nmiss .eq. 0 ) then
            write(*,*) '... highly successful, all points got assigned'
        elseif ( nhit .eq. 0 ) then
            write(*,*) '... did not go very well, nothing assigned'
        else
            write(*,*) '... moderate success, some assigned'
        endif
    endif

    ! Now that the full BZ is somewhat labelled, try to label the wedge as well.
    ! fetch the coordinates for the wedge
    if ( allocated(p%irrw%label) ) deallocate(p%irrw%label)
    lo_allocate(p%irrw%label(p%irrw%nnodes))
    lo_allocate(rw(3,p%irrw%nnodes))
    do i=1,p%irrw%nnodes
        p%irrw%label(i)='NP'
    enddo
    do i=1,p%irrw%nnodes
        ! to fractional
        rw(:,i)=p%cartesian_to_fractional( p%irrw%r(:,i), reciprocal=.true., pbc=.false. )
        ! then strange permutation
        rw(:,i)=matmul(rw(:,i),p%info%permutation_to_unique)
        ! to Cartesian again, in the new coordinate system
        rw(:,i)=matmul(bas,rw(:,i))
    enddo

    ! Figure out labels for the wedge nodes
    nhit=0
    nmiss=0
    ! first try without any operations, might have success with that.
    wedgeptloop1: do i=1,p%irrw%nnodes
        do j=1,p%bz%nhighsymmetrypoints
            if ( lo_sqnorm(rw(:,i)-r(:,j)) .lt. lo_sqtol ) then
                ! I found a point, if it is labelled
                if ( trim(p%bz%label(j)) .ne. 'NP' ) then
                    nhit=nhit+1
                    p%irrw%label(i)=trim(p%bz%label(j))
                    cycle wedgeptloop1
                endif
            endif
        enddo
        ! if I made it here, I guess I did not find a label.
        nmiss=nmiss+1
    enddo wedgeptloop1

    ! If not successful, try again with all operations. It might work.
    if ( nmiss .ne. 0 ) then
        !
        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) '... not complete success first pass, trying again'
        endif
        ! New attempt, with all operations
        wedgeptloop2: do i=1,p%irrw%nnodes
            ! skip if already fixed
            if ( trim(p%irrw%label(i)) .ne. 'NP' ) cycle
            ! test vs all BZ points
            do j=1,p%bz%nhighsymmetrypoints
                ! skip if this one also failed
                if ( trim(p%bz%label(j)) .eq. 'NP' ) cycle
                ! test vs all operations
                do o=1,sym%n
                    ! rotate the wedge point
                    v0=lo_operate_on_vector(sym%op(o),rw(:,i),reciprocal=.true.)
                    if ( lo_sqnorm(v0-r(:,j)) .lt. lo_sqtol ) then
                        ! I found a point, I guess
                        nhit=nhit+1
                        p%irrw%label(i)=trim(p%bz%label(j))
                        cycle wedgeptloop2
                    endif
                    if ( sym%timereversal ) then
                        if ( lo_sqnorm(v0+r(:,j)) .lt. lo_sqtol ) then
                            nhit=nhit+1
                            p%irrw%label(i)=trim(p%bz%label(j))
                            cycle wedgeptloop2
                        endif
                    endif
                enddo
            enddo
            ! if I made it here, I guess I did not find a label.
            nmiss=nmiss+1
        enddo wedgeptloop2
    endif

    ! Report on success
    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) '... labelled wedge: '//tochar(nhit)//' hits, '//tochar(nmiss)//' misses'
    endif

    ! Name the new points in the wedge, maybe
    if ( nmiss .ne. 0 ) then
        ! name the points in the wedge that don't have real names
        j=0
        do i=1,p%irrw%nnodes
            if ( trim(p%irrw%label(i)) .eq. 'NP' ) then
                j=j+1
                p%irrw%label(i)='NP'//tochar(j)
            endif
        enddo

        ! Now go over the wedge and spread these names:
        bzmissloop: do i=1,p%bz%nhighsymmetrypoints
            if ( trim(p%bz%label(i)) .ne. 'NP' ) cycle
            do j=1,p%irrw%nnodes
                do o=1,sym%n
                    v0=lo_operate_on_vector(sym%op(o),rw(:,j),reciprocal=.true.)
                    if ( lo_sqnorm(v0-r(:,i)) .lt. lo_sqtol ) then
                        p%bz%label(i)=trim(p%irrw%label(j))
                        cycle bzmissloop
                    endif
                    if ( sym%timereversal ) then
                        if ( lo_sqnorm(v0+r(:,i)) .lt. lo_sqtol ) then
                            p%bz%label(i)=trim(p%irrw%label(j))
                            cycle bzmissloop
                        endif
                    endif
                enddo
            enddo
            ! Hopefully, I will never make it here
            write(*,*) 'Completely failed labelling point ',tochar(i)
            stop
        enddo bzmissloop
    endif
    p%info%pointslabelled=.true.
    if ( p%info%verbosity .gt. 0 ) write(*,*) 'Did my best at labelling points'
    !
end subroutine


!> Check if three numbers are equal within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_equal(a,b,c,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-b) .lt. thres ) i=i+1
    if ( abs(a-c) .lt. thres ) i=i+1
    if ( abs(b-c) .lt. thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_equal=.true.
    else
        lo_three_equal=.false.
    endif
end function

!> Check if exactly one number is equal to some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_one_equal_to_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-num) .lt. thres ) i=i+1
    if ( abs(b-num) .lt. thres ) i=i+1
    if ( abs(c-num) .lt. thres ) i=i+1
    !
    if ( i .eq. 1 ) then
        lo_one_equal_to_num=.true.
    else
        lo_one_equal_to_num=.false.
    endif
end function

!> Check if exactly two numbers are equal to some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_two_equal_to_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-num) .lt. thres ) i=i+1
    if ( abs(b-num) .lt. thres ) i=i+1
    if ( abs(c-num) .lt. thres ) i=i+1
    !
    if ( i .eq. 2 ) then
        lo_two_equal_to_num=.true.
    else
        lo_two_equal_to_num=.false.
    endif
end function

!> Check if three numbers are equal to some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_equal_to_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-num) .lt. thres ) i=i+1
    if ( abs(b-num) .lt. thres ) i=i+1
    if ( abs(c-num) .lt. thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_equal_to_num=.true.
    else
        lo_three_equal_to_num=.false.
    endif
end function

!> Check if three numbers are less than some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_smaller_than_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check how many equal
    i=0
    if ( a .lt. num-thres ) i=i+1
    if ( b .lt. num-thres ) i=i+1
    if ( c .lt. num-thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_smaller_than_num=.true.
    else
        lo_three_smaller_than_num=.false.
    endif
end function

!> Check if three numbers are greater than some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_larger_than_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check how many equal
    i=0
    if ( a .gt. num-thres ) i=i+1
    if ( b .gt. num-thres ) i=i+1
    if ( c .gt. num-thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_larger_than_num=.true.
    else
        lo_three_larger_than_num=.false.
    endif
end function

!> Check if three numbers are less than some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_two_smaller_than_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check how many equal
    i=0
    if ( a .lt. num-thres ) i=i+1
    if ( b .lt. num-thres ) i=i+1
    if ( c .lt. num-thres ) i=i+1
    !
    if ( i .eq. 2 ) then
        lo_two_smaller_than_num=.true.
    else
        lo_two_smaller_than_num=.false.
    endif
end function

!> Check if two out of three numbers are equal within some tolerance. If no tolerance is specified, I guess one. Will return false for three equal numbers.
logical pure function lo_two_equal(a,b,c,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    !
    i=0
    if ( abs(a-b) .lt. thres ) i=i+1
    if ( abs(a-c) .lt. thres ) i=i+1
    if ( abs(b-c) .lt. thres ) i=i+1
    !
    if ( i .eq. 1 ) then
        lo_two_equal=.true.
    else
        lo_two_equal=.false.
    endif
end function

!> Check that none of three numbers are equal within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_none_equal(a,b,c,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    !
    i=0
    if ( abs(a-b) .lt. thres ) i=i+1
    if ( abs(a-c) .lt. thres ) i=i+1
    if ( abs(b-c) .lt. thres ) i=i+1
    !
    if ( i .eq. 0 ) then
        lo_none_equal=.true.
    else
        lo_none_equal=.false.
    endif
end function

end submodule
