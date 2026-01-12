#include "precompilerdefinitions"
!! Get the symmetry operations for a structure in a neat way.
module type_symmetryoperation
use konstanter, only: i8,r8, lo_huge, lo_hugeint, lo_pi, lo_twopi, lo_tol, lo_sqtol, lo_radiantol, &
                      lo_exitcode_symmetry, lo_exitcode_param, lo_exitcode_baddim, lo_status
use gottochblandat, only: tochar,walltime,lo_clean_fractional_coordinates,lo_chop,lo_frobnorm,lo_invert3x3matrix,lo_determ,&
                          lo_sqnorm,qsort,lo_fetch_tolerance, lo_return_unique, lo_reciprocal_basis, lo_stop_gracefully,&
                          lo_flattentensor,lo_unflatten_2tensor,lo_invert_real_matrix,lo_real_pseudoinverse,&
                          lo_general_real_eigenvalues_eigenvectors,lo_triplegemm
use geometryfunctions, only: lo_rotation_matrix_from_axis_and_angle,lo_improper_rotation_matrix_from_axis_and_angle
implicit none

private
public :: lo_spacegroup_operation
public :: lo_symset
public :: lo_operate_on_vector
public :: lo_operate_on_secondorder_tensor
public :: lo_eigenvector_transformation_matrix
public :: lo_expandoperation_pair
public :: lo_expandoperation_triplet
public :: lo_expandoperation_quartet

!> symmetry operation
type lo_spacegroup_operation
    !> rotation in cartesian coordinates
    real(r8), dimension(3,3) :: m=lo_huge
    !> inverse rotation in cartesian coordinates
    real(r8), dimension(3,3) :: im=lo_huge
    !> operation in realspace fractional coordinates
    real(r8), dimension(3,3) :: fm=lo_huge
    !> inverse operation in realspace fractional coordinates
    real(r8), dimension(3,3) :: ifm=lo_huge
    !> operation in reciprocal fractional coordinates
    real(r8), dimension(3,3) :: rfm=lo_huge
    !> inverse operation in reciprocal fractional coordinates
    real(r8), dimension(3,3) :: irfm=lo_huge

    !> translation, cartesian coordinates in realspace
    real(r8), dimension(3) :: tr=lo_huge
    !>  translation, fractional coordinates in realspace
    real(r8), dimension(3) :: ftr=lo_huge

    !> what kind of operation
    character(len=20) :: whatkind='dunno'
    !> if it's a rotation, what kind?
    integer :: fold=-lo_hugeint
    !> axis of rotation or rotoinversion
    real(r8), dimension(3) :: axis=lo_huge
    !> K=F(k,S) function per Maradudin and Vosko, what atom gets transformed to what atom.
    integer, dimension(:), allocatable :: fmap

    ! !> in terms of euler angles in the y-convention, what is the alpha angle
    ! real(r8) :: fwd_alpha=-lo_huge,bwd_alpha=-lo_huge
    ! !> beta angle
    ! real(r8) :: fwd_beta=-lo_huge,bwd_beta=-lo_huge
    ! !> gamma angle
    ! real(r8) :: fwd_gamma=-lo_huge,bwd_gamma=-lo_huge
end type

!> conjugacy classes of the symmetry operations
type lo_symset_conjugacyclass
    !> number of members in class
    integer :: n_member=-lo_hugeint
    !> members
    integer, dimension(:), allocatable :: member
end type

!> irreducible representations
type lo_symset_irreducible_representation
    !> dimension of representation
    integer :: dimension=-lo_hugeint
    !> character per class
    real(r8), dimension(:), allocatable :: character_class
    !> character per element
    real(r8), dimension(:), allocatable :: character_element
end type

!> Set of symmetry operations
type lo_symset
    !> number of operations
    integer :: n=-lo_hugeint
    !> operations
    type(lo_spacegroup_operation), dimension(:), allocatable :: op
    !> Should time-reversal symmetry be enforced?
    logical :: timereversal=.false.
    !> What is the degeneracy of each atom? As in which atom is equivalent to which?
    integer, dimension(:), allocatable :: degeneracy
    !> Which atoms is it degenerate with
    integer, dimension(:,:), allocatable :: degenerate_atom

    !> how many conjugacy classes are there
    integer :: n_class=-lo_hugeint
    !> how many irreducible representations are there
    integer :: n_irrep=-lo_hugeint
    !> conjugacy classes
    type(lo_symset_conjugacyclass), dimension(:), allocatable :: cc
    !> irreducible represenations
    type(lo_symset_irreducible_representation), dimension(:), allocatable :: ir

    !> how many unique atoms are there
    integer :: n_irreducible_atom=-lo_hugeint
    !> index to the unique atoms (n_irreducible_atom)
    integer, dimension(:), allocatable :: irr_to_all
    !> for all the atoms, which unique atom is it? (n_atom)
    integer, dimension(:), allocatable :: all_to_irr
    !> for all the atoms, which operations takes the irreducible atom to this atom
    integer, dimension(:), allocatable :: op_all_from_irr
    !> for each irreducible atom, how many atoms does it unfold to?
    integer, dimension(:), allocatable :: irr_unfold_ctr
    !> for each irreducible atom, which atoms does it unfold to?
    integer, dimension(:,:), allocatable :: irr_unfold_index
    !> for each irreducible atom, which operation does it unfold with?
    integer, dimension(:,:), allocatable :: irr_unfold_operation

    !> multiplication table, for forward-forward operations
    integer, dimension(:,:), allocatable :: multiplication_table

    contains
        !> get the reduction arrays, op(i)*op(j)=op(k)
        procedure :: reductionlist
        !> get spacegroup
        procedure :: generate
        !> calculate the character table
        procedure :: get_character_table
        !> size in memory
        procedure :: size_in_mem
end type

!> Helper box to speed up symmetry detection
type lo_speciesbox
    !> how many of this species
    integer :: n=-lo_hugeint
    !> the points
    real(r8), dimension(:,:), allocatable :: r
    !> the original index of these points
    integer, dimension(:), allocatable :: ind
end type

contains

!> Apply an operation to a vector
pure function lo_operate_on_vector(op,v,reciprocal,inverse,fractional) result(w)
    !> operation
    type(lo_spacegroup_operation), intent(in) :: op
    !> vector
    real(r8), dimension(3), intent(in) :: v
    !> should the inverse operation be applied?
    logical, intent(in), optional :: inverse
    !> is it an operation in reciprocal space
    logical, intent(in), optional :: reciprocal
    !> should it be done in fractional coordinates?
    logical, intent(in), optional :: fractional
    !> rotated vector
    real(r8), dimension(3) :: w

    logical :: rs,inv,fr
    real(r8), dimension(3,3) :: m
    real(r8), dimension(3) :: t

    ! Ok, eight possibilites, forward or inverse, realspace or reciprocal, cartesian or fractional
    if ( present(inverse) ) then
        inv=inverse
    else
        inv=.false.
    endif
    if ( present(reciprocal) ) then
        rs=reciprocal
    else
        rs=.false.
    endif
    if ( present(fractional) ) then
        fr=fractional
    else
        fr=.false.
    endif

    ! Get the rotation part.
    if ( fr ) then
        if ( inv ) then
            if ( rs ) then
                m=op%irfm ! inverse, reciprocal,fractional
            else
                m=op%ifm  ! inverse, realspace, fractional
            endif
        else
            if ( rs ) then
                m=op%rfm  ! forward, reciprocal, fractional
            else
                m=op%fm   ! forward, realspace, fractional
            endif
        endif
    else
        if ( inv ) then
            m=op%im       ! Cartesian, real/reciprocal, inverse
        else
            m=op%m        ! Cartesian, real/reciprocal, forward
        endif
    endif
    ! get the translation part
    if ( rs ) then
        t=0.0_r8
    else
        if ( fr ) then
            t=op%ftr
        else
            t=op%tr
        endif
    endif
    ! Begin operation "actual operation"
    if ( inv ) then
        w=v-t
        w=matmul(m,w)
    else
        w=matmul(m,v)
        w=w+t
    endif
end function

!> Apply operation to a 3x3 tensor
pure function lo_operate_on_secondorder_tensor(op,m) result(n)
    !> The operation
    type(lo_spacegroup_operation), intent(in) :: op
    !> The original tensor
    real(r8), dimension(3,3), intent(in) :: m
    !> The rotated tensor
    real(r8), dimension(3,3) :: n

    integer :: i,j,ii,jj
    n=0.0_r8
    do j=1,3
    do i=1,3
        do jj=1,3
        do ii=1,3
            n(i,j)=n(i,j)+m(ii,jj)*op%m(i,ii)*op%m(j,jj)
        enddo
        enddo
    enddo
    enddo
end function

!> Transformation matrix for phonon eigenvectors corresponding to a symmetry operation
pure subroutine lo_eigenvector_transformation_matrix(gm,rcart,qv,op,inverseoperation)
    !> transformation matrix
    complex(r8), dimension(:,:), intent(out) :: gm
    !> position vectors in Cartesian coordinates, generally crystalstructure%rcart
    real(r8), dimension(:,:), intent(in) :: rcart
    !> the q-vector
    real(r8), dimension(3), intent(in) :: qv
    !> the point operation
    type(lo_spacegroup_operation), intent(in) :: op
    !> should the inverse transformation be used?
    logical, intent(in), optional :: inverseoperation
    !
    complex(r8) :: expiqr
    real(r8), dimension(3) :: v0
    real(r8) :: iqvv
    integer :: a1,a2,na
    logical :: inv

    ! Forward or backwards operation
    if ( present(inverseoperation) ) then
        inv=inverseoperation
    else
        inv=.false.
    endif

    na=size(rcart,2)
    gm=0.0_r8
    do a2=1,na
        a1=op%fmap(a2)
        v0=lo_operate_on_vector(op,rcart(:,a1),inverse=.true.)-rcart(:,a2)
        iqvv=-dot_product(v0,qv)*lo_twopi
        expiqr=cmplx(cos(iqvv),sin(iqvv),r8)
        if ( inv ) then
            gm( (a1-1)*3+1:a1*3 , (a2-1)*3+1:a2*3 )=op%im*expiqr
        else
            gm( (a1-1)*3+1:a1*3 , (a2-1)*3+1:a2*3 )=op%m*expiqr
        endif
    enddo
end subroutine

!> Returns the symmetry operations of the space group of a crystal structure. If the set of points is omitted, the allowed operations of the parent lattice is returned.
subroutine generate(sym,basis,timereversal,points,pointclassification,pointflavor,possible,symmorphic,verbosity,tolerance)
    !> set
    class(lo_symset), intent(out) :: sym
    !> basis
    real(r8), dimension(3,3), intent(in) :: basis
    !> should time-reversal symmetry be enforced?
    logical, intent(in) :: timereversal
    !> points in the lattice
    real(r8), dimension(:,:), intent(in), optional :: points
    !> classification of the points if they are not indistinguishable
    integer, dimension(:), intent(in), optional :: pointclassification
    !> perhaps they are classified in more than one way
    integer, dimension(:), intent(in), optional :: pointflavor
    !> return the possible operations, not the valid
    logical, intent(in), optional :: possible
    !> generate only the symmorphic operations
    logical, intent(in), optional :: symmorphic
    !> how much to talk
    integer, intent(in), optional :: verbosity
    !> non-default tolerance?
    real(r8), intent(in), optional :: tolerance

    ! Dummy points, species and flavors
    type(lo_speciesbox), dimension(:), allocatable :: sbox
    real(r8), dimension(:,:), allocatable :: dumpts
    integer, dimension(:), allocatable :: dumclass,dumflavor
    ! Dummy operations
    real(r8), dimension(:,:,:), allocatable :: op,fop
    real(r8), dimension(:,:), allocatable :: tr,ftr
    ! Local stuff
    real(r8), dimension(3,3) :: inversebasis,reciprocalbasis,inversereciprocalbasis
    real(r8) :: timer,tol,t0,t1
    integer :: verb,np
    logical :: only_possible,only_symmorphic

    ! set some things and figure out exactly what to do.
    init: block
        integer, dimension(:), allocatable :: dumi
        integer :: nspecies,i,j
        ! How much to talk?
        if ( present(verbosity) ) then
            timer=walltime()
            t0=timer
            t1=0.0_r8
            verb=verbosity
        else
            verb=0
        endif

        if ( verb .gt. 0 ) then
            write(*,*) ' '
            write(*,*) 'Generating spacegroup operations'
        endif

        ! Initialize everything to nothing.
        sym%n=-1
        sym%n_class=-1
        sym%n_irrep=-1

        ! Some tolerances
        if ( present(tolerance) ) then
            tol=tolerance
        else
            tol=lo_tol
        endif
        ! consider timereversal symmetries?
        sym%timereversal=timereversal

        ! I have to sort out how much to care when checking the symmetries.
        ! First, are there any points provided?
        if ( present(points) ) then
            np=size(points,2)
            lo_allocate(dumpts(3,np))
            dumpts=points
        else
            lo_allocate(dumpts(3,1))
            dumpts=0.0_r8
            np=1
        endif
        ! Ok, other things that could be provided is that the points have a species
        ! and a flavor. Not using the flavor at the moment, but to be safe for the
        ! future.
        lo_allocate(dumclass(np))
        lo_allocate(dumflavor(np))
        dumclass=1
        dumflavor=0

        ! Any species provided?
        if ( present(pointclassification) ) then
            if ( size(pointclassification) .ne. np ) then
                call lo_stop_gracefully(['You need to have a classification for every point'],lo_exitcode_param)
            endif
            dumclass=pointclassification
        endif
        ! Any extra flavor on the points?
        if ( present(pointflavor) ) then
            if ( size(pointflavor,1) .ne. np ) then
                call lo_stop_gracefully(['You need to have a flavor for every point'],lo_exitcode_param)
            endif
            dumflavor=pointflavor
        endif

        ! Perhaps just return the possible operations, not the valid?
        only_symmorphic=.false.
        if ( present(possible) ) then
            only_possible=possible
        else
            only_possible=.false.
            only_symmorphic=.true.
        endif
        ! Only the symmorphic part of the space group?
        if ( present(symmorphic) ) then
            only_symmorphic=symmorphic
        else
            only_symmorphic=.false.
        endif

        ! Get all the versions of basis and stuff
        reciprocalbasis=lo_chop( lo_reciprocal_basis(basis) , 1E-11_r8 )
        inversereciprocalbasis=lo_invert3x3matrix( reciprocalbasis )
        inversebasis=lo_invert3x3matrix( basis )

        ! We might need the speciesboxes:
        if ( only_possible .eqv. .false. ) then
            ! count the number of different species
            call lo_return_unique(dumclass,dumi)
            nspecies=size(dumi,1)
            lo_allocate(sbox(nspecies))
            ! count species of each kind
            do i=1,nspecies
                sbox(i)%n=0
            enddo
            do i=1,np
                do j=1,nspecies
                    if ( dumi(j) .eq. dumclass(i) ) sbox(j)%n=sbox(j)%n+1
                enddo
            enddo
            ! make space
            do i=1,nspecies
                lo_allocate(sbox(i)%r(3,sbox(i)%n))
                lo_allocate(sbox(i)%ind(sbox(i)%n))
                sbox(i)%n=0
            enddo
            ! stuff the boxes
            do i=1,np
                do j=1,nspecies
                    if ( dumi(j) .eq. dumclass(i) ) then
                        sbox(j)%n=sbox(j)%n+1
                        sbox(j)%r(:,sbox(j)%n)=dumpts(:,i)
                        sbox(j)%ind(sbox(j)%n)=i
                    endif
                enddo
            enddo
            lo_deallocate(dumi)
            if ( verb .gt. 0 ) then
                t1=walltime()
                write(*,*) '... rearranged atoms per species in '//tochar(t1-t0)//'s'
                t0=t1
            endif
        endif

        ! Report progress, so far
        if ( verb .gt. 0 ) then
            if ( present(points) ) write(*,*) '... considering '//tochar(np)//' points in the cell'
            if ( present(pointclassification) ) write(*,*) '... found '//tochar(maxval(dumclass))//' different species'
            if ( present(pointflavor) ) write(*,*) '... found '//tochar(maxval(dumflavor))//' different flavors'
        endif
    end block init

    ! If we only want the possible operations, just pack them up and we are done
    if ( only_possible ) then
    symmorphicops: block
        ! Dummy operations
        real(r8), dimension(3,3) :: primlat
        integer :: o

        ! Find lattice vectors of the smallest cell
        call find_primitive_cell(basis,dumpts,dumclass,dumflavor,primlat,tol)
        ! A list of the possible operations
        call return_all_operations(op,fop,primlat,basis,tol)
        if ( verb .gt. 0 ) write(*,*) '... found '//tochar(size(op,3))//' possible proper/improper rotations in '&
                                      //tochar(walltime()-timer)//'s'

        sym%n=size(op,3)
        lo_allocate(sym%op(sym%n))
        do o=1,sym%n
            ! point operations
            sym%op(o)%m=op(:,:,o)
            sym%op(o)%im=transpose(sym%op(o)%m)
            sym%op(o)%fm=fop(:,:,o)
            sym%op(o)%ifm=anint( lo_invert3x3matrix( sym%op(o)%fm ) )
            sym%op(o)%rfm=lo_chop( matmul(matmul(inversereciprocalbasis,sym%op(o)%m),reciprocalbasis),1E-10_r8 )
            sym%op(o)%irfm=lo_chop( lo_invert3x3matrix(sym%op(o)%rfm),1E-10_r8 )
            ! translations
            sym%op(o)%tr=0.0_r8
            sym%op(o)%ftr=0.0_r8
            ! the fmap is trivial
            lo_allocate(sym%op(o)%fmap(1))
            ! Better to set it to something stupid, I think
            sym%op(o)%fmap=-lo_hugeint
        enddo
        ! classify the operations
        call figure_out_what_operations_are(sym)
        ! Some cleanup
        lo_deallocate(op)
        lo_deallocate(fop)
        lo_deallocate(dumpts)
        lo_deallocate(dumclass)
        lo_deallocate(dumflavor)

        if ( verb .gt. 0 ) write(*,*) 'Returned the possible operations in '//tochar(walltime()-timer)//'s'
        ! Done here!
        return
    end block symmorphicops
    endif

    ! Ok, if we made it here, we want the full spacegroup.
    spaceops: block
        real(r8), dimension(:,:,:), allocatable :: dumop1
        real(r8), dimension(:,:), allocatable :: dumtr1,dr0
        real(r8), dimension(:), allocatable :: dv
        real(r8), dimension(3,3) :: primlat,m0
        real(r8), dimension(3) :: v0
        integer, dimension(:,:), allocatable :: di
        integer, dimension(:), allocatable :: fmp
        integer :: i,j,k,l,o,ii
        logical :: opvalid

        t0=walltime()
        ! Find lattice vectors of the smallest cell
        call find_primitive_cell(basis,dumpts,dumclass,dumflavor,primlat,tol)
        ! A list of the possible operations
        call return_all_operations(op,fop,primlat,basis,tol)
        if ( verb .gt. 0 ) write(*,*) '... found '//tochar(size(op,3))//' possible proper/improper rotations in '&
                                      //tochar(walltime()-timer)//'s'
        ! Get the possible translations
        if ( only_symmorphic .eqv. .true. ) then
            ! This means I only want operations with no translations.
            lo_allocate(tr(3,1))
            lo_allocate(ftr(3,1))
            tr=0.0_r8
            ftr=0.0_r8
        else
            call return_possible_translations(tr,ftr,op,basis,inversebasis,dumpts,dumclass,tol)
        endif
        if ( verb .gt. 0 ) then
            t1=walltime()
            write(*,*) '... found '//tochar(size(tr,2))//' possible translations in '//tochar(t1-t0)//'s'
            t0=t1
        endif

        ! Get a first rough count
        l=size(op,3)*size(tr,2)
        lo_allocate(dumop1(3,3,l))
        lo_allocate(dumtr1(3,l))
        lo_allocate(fmp(np))
        fmp=0
        dumop1=0.0_r8
        dumtr1=0.0_r8
        l=0
        do i=1,size(op,3)
        do j=1,size(tr,2)
            ! test this operation
            call check_operation(fop(:,:,i),ftr(:,j),sbox,basis,fmp,opvalid,tol)
            ! maybe store it
            if ( opvalid ) then
                l=l+1
                dumop1(:,:,l)=fop(:,:,i)
                dumtr1(:,l)=ftr(:,j)
            endif
        enddo
        enddo
        ! Kill the big arrays to replace with more conveniently sized ones
        lo_deallocate(op)
        lo_deallocate(fop)
        lo_deallocate(tr)
        lo_deallocate(ftr)

        ! Now fill it out with combinations. Only needed far black phosphorous
        ! so far, but keeping it in all cases, might be important at some point.
        !@TODO probably revisit this and make sure it does not get slow.
        lo_allocate(dr0(12,l**2))
        lo_allocate(dv(12))
        dr0=0.0_r8
        k=0
        do i=1,l
        jl: do j=1,l
            m0=matmul(dumop1(:,:,i),dumop1(:,:,j))
            v0=matmul(dumop1(:,:,i),dumtr1(:,j))+dumtr1(:,i)
            v0=lo_clean_fractional_coordinates(v0)
            dv(1:3)=v0
            dv(4:12)=lo_flattentensor(m0)
            do ii=1,k
                if ( sum(abs(dv-dr0(:,ii))) .lt. lo_sqtol*10 ) cycle jl
            enddo
            call check_operation(m0,v0,sbox,basis,fmp,opvalid,tol)
            ! Store only unique @todo Check if I should have test for opvalid here
            k=k+1
            dr0(:,k)=dv
        enddo jl
        enddo
        lo_deallocate(dumop1)
        lo_deallocate(dumtr1)
        lo_deallocate(dv)

        if ( verb .gt. 0 ) then
            t1=walltime()
            write(*,*) '... found '//tochar(k)//' valid operations in '//tochar(t1-t0)//'s'
            t0=t1
        endif

        ! Now store this
        sym%n=k
        lo_allocate(sym%op(sym%n))
        do o=1,sym%n
            ! Actual operation, fractional coordinates. Check that it is not stupid.
            m0=lo_unflatten_2tensor(dr0(4:12,o))
            ! Actual translation, in fractional coordinates
            v0=dr0(1:3,o)

            ! Not sure if this is actually a good check. Have to think.
            if ( sum(abs(m0-anint(m0))) .gt. lo_tol ) then
                call lo_stop_gracefully(['Invalid operation. Could never happen, unless you are not using the primitive cell. Which you should.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif

            ! Check I did not do something stupid, and get the fmap
            call check_operation(m0,v0,sbox,basis,fmp,opvalid,tol)
            if ( opvalid .eqv. .false. ) then
                call lo_stop_gracefully(['Invalid operation. Could never happen.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif

            ! Store fractional operation
            sym%op(o)%fm=anint(m0)
            sym%op(o)%ifm=anint(lo_invert3x3matrix(m0))
            ! Then make Cartesian operations, and check that they are not stupid
            m0=matmul(basis,sym%op(o)%fm)
            sym%op(o)%m=lo_chop( matmul(m0,inversebasis), 1E-10_r8 )
            sym%op(o)%im=transpose(sym%op(o)%m)
            if ( sum(abs(sym%op(o)%im-lo_invert3x3matrix(sym%op(o)%m))) .gt. lo_tol ) then
                call lo_stop_gracefully(['Invalid operation. Could never happen.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
            ! Reciprocal fractional operation
            sym%op(o)%rfm=lo_chop(matmul(matmul(inversereciprocalbasis,sym%op(o)%m),reciprocalbasis),1E-10_r8)
            sym%op(o)%irfm=lo_chop( lo_invert3x3matrix(sym%op(o)%rfm),1E-10_r8 )
            ! translations
            sym%op(o)%ftr=v0
            sym%op(o)%tr=lo_chop(matmul(basis,v0),1E-11_r8)
            ! the K=F(k,S) map
            lo_allocate(sym%op(o)%fmap( np ))
            sym%op(o)%fmap=fmp
        enddo

        ! Classify the operations? Could be meaningful.
        call figure_out_what_operations_are(sym)

        ! Figure out degeneracies
        lo_allocate(sym%degeneracy(np))
        lo_allocate(di(np,np))
        di=0
        do o=1,sym%n
            do i=1,np
                j=sym%op(o)%fmap(i)
                di(i,j)=1
                di(j,i)=1
            enddo
        enddo
        ! Store the degeneracy
        sym%degeneracy=0
        do i=1,np
            sym%degeneracy(i)=sum(di(:,i))
        enddo
        lo_allocate(sym%degenerate_atom(maxval(sym%degeneracy),np))
        sym%degenerate_atom=0
        do i=1,np
            l=0
            do j=1,np
                if ( di(j,i) .eq. 1 ) then
                    l=l+1
                    sym%degenerate_atom(l,i)=j
                endif
            enddo
        enddo

        ! Some cleanup
        lo_deallocate(fmp)
        do i=1,size(sbox,1)
            lo_deallocate(sbox(i)%r)
            lo_deallocate(sbox(i)%ind)
        enddo
        lo_deallocate(sbox)
        lo_deallocate(dumpts)
        lo_deallocate(dumflavor)
        lo_deallocate(di)
    end block spaceops

    ! Get the unique atoms
    uniqueatoms: block
        real(r8), dimension(3,3) :: m0
        real(r8), dimension(3) :: v0,v1
        integer, dimension(:), allocatable :: di
        integer :: i,j,k,l,o,a1,a2,ctr

        ! Use the list of operations to check which atom transforms to which
        lo_allocate(di(np))
        di=1
        do o=1,sym%n
            do a1=1,np
                a2=sym%op(o)%fmap(a1)
                if ( a2 .gt. a1 ) di(a2)=0
            enddo
        enddo

        ! Now I know how many irreducible atoms there are.
        sym%n_irreducible_atom=sum(di)

        ! Make some space
        lo_allocate(sym%irr_to_all(sym%n_irreducible_atom))
        lo_allocate(sym%all_to_irr(np))
        lo_allocate(sym%op_all_from_irr(np))
        sym%irr_to_all=0
        sym%all_to_irr=0
        sym%op_all_from_irr=0

        ! Store the indices to the irreducible
        i=0
        do a1=1,np
            if ( di(a1) .eq. 0 ) cycle
            i=i+1
            sym%irr_to_all(i)=a1
        enddo

        ! Slightly wasteful way of mapping backwards, but I think it's fine in
        ! terms of speed. Can make it O(N) if necessary.
        atl1: do a1=1,np
            do o=1,sym%n
            do i=1,sym%n_irreducible_atom
                j=sym%irr_to_all(i)
                a2=sym%op(o)%fmap(j)
                if ( a2 .eq. a1) then
                    sym%all_to_irr(a1)=i
                    sym%op_all_from_irr(a1)=o
                    cycle atl1
                endif
            enddo
            enddo
        enddo atl1

        ! Get the unfolding map as well, while I'm at it.
        ! First count the degeneracy
        allocate(sym%irr_unfold_ctr(np))
        sym%irr_unfold_ctr=0
        do a1=1,np
            j=sym%all_to_irr(a1)
            sym%irr_unfold_ctr(j)=sym%irr_unfold_ctr(j)+1
        enddo
        ! Make some space
        j=maxval(sym%irr_unfold_ctr)
        allocate(sym%irr_unfold_index(j,np))
        allocate(sym%irr_unfold_operation(j,np))
        sym%irr_unfold_index=-1
        sym%irr_unfold_operation=-1
        ! Store the data
        sym%irr_unfold_ctr=0
        do a1=1,np
            j=sym%all_to_irr(a1)
            o=sym%op_all_from_irr(a1)
            sym%irr_unfold_ctr(j)=sym%irr_unfold_ctr(j)+1
            sym%irr_unfold_index(sym%irr_unfold_ctr(j),j)=a1
            sym%irr_unfold_operation(sym%irr_unfold_ctr(j),j)=o
        enddo

        ! And get the multiplication table
        allocate(sym%multiplication_table(sym%n,sym%n))
        sym%multiplication_table=0
        do i=1,sym%n
        do j=1,sym%n
            ! How do you chain operations again?
            m0=matmul(sym%op(i)%fm,sym%op(j)%fm)
            ! And the translation
            v0=matmul(sym%op(i)%fm,sym%op(j)%ftr)+sym%op(i)%ftr
            v0=lo_clean_fractional_coordinates(v0)
            ctr=0
            l=0
            do k=1,sym%n
                if ( sum(abs(m0-sym%op(k)%fm)) .lt. lo_tol ) then
                    v1=lo_clean_fractional_coordinates(v0-sym%op(k)%ftr+0.5_r8)-0.5_r8
                    if ( sum(abs(v1)) .lt. lo_tol ) then
                        l=k
                        ctr=ctr+1
                    endif
                endif
            enddo
            if ( ctr .eq. 1 ) then
                sym%multiplication_table(i,j)=l
            else
                call lo_stop_gracefully(['The group is not a group, but I thought it was.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        enddo

        ! And cleanup
        lo_deallocate(dumclass)
        lo_deallocate(di)
    end block uniqueatoms

    ! Final report
    if ( verb .gt. 0 ) then
        write(*,*) '... identified and classified operations in '//tochar(walltime()-timer)//'s'
        write(*,*) 'Generated the relevant spacegroup operations.'
    endif
end subroutine

!> Return the valid point operations for this basis
pure subroutine return_all_operations(ops,fracops,primbasis,basis,tolerance)
    !> basis
    real(r8), dimension(3,3), intent(in) :: primbasis
    !> basis for the lattice I actually want
    real(r8), dimension(3,3), intent(in) :: basis
    !> point operations
    real(r8), dimension(:,:,:), allocatable, intent(out) :: ops
    !> point operations in fractional coordinates
    real(r8), dimension(:,:,:), allocatable, intent(out) :: fracops
    !> tolerance
    real(r8), intent(in) :: tolerance

    real(r8), dimension(3,3,96) :: bflist
    real(r8), dimension(:,:,:), allocatable :: list
    real(r8), dimension(3,3,48) :: fraclist
    real(r8), dimension(3,3) :: inversebasis,m,ident,canbasis,caninvbasis
    integer, dimension(3) :: fracvals
    integer, dimension(3,3) :: op
    integer :: i,j,k,l
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! do a quick Delaunay reduction of the basis to canonical form, in case it's not done already
    call reduce_basis_to_something_decent(primbasis,canbasis)

    ! I need the inverse basis of both the normal basis and the reduced.
    caninvbasis=lo_invert3x3matrix(canbasis)
    inversebasis=lo_invert3x3matrix(basis)

    ! And an identity matrix
    ident=0.0_r8
    ident(1,1)=1.0_r8
    ident(2,2)=1.0_r8
    ident(3,3)=1.0_r8

    ! Build all possible operations, all 3^9 if them.
    fracvals=[-1,0,1]
    k=0
    il1: do i1=1,3
    do i2=1,3
    do i3=1,3
    do i4=1,3
    do i5=1,3
    do i6=1,3
    do i7=1,3
    do i8=1,3
    do i9=1,3
        ! possible operation in fractional coordinates
        op(:,1)=[fracvals(i1),fracvals(i2),fracvals(i3)]
        op(:,2)=[fracvals(i4),fracvals(i5),fracvals(i6)]
        op(:,3)=[fracvals(i7),fracvals(i8),fracvals(i9)]
        ! Convert it to Cartesian coordinates
        m=matmul(canbasis,op)
        m=matmul(m,caninvbasis)
        ! Check that op*op^T=op*op^-1
        if ( lo_frobnorm( matmul( m,transpose(m) )-ident ) .lt. tolerance ) then
            k=k+1
            bflist(:,:,k)=m
        endif
        ! done when I hit 48
        if ( k .eq. 48 ) exit il1
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo il1

    ! This is to catch some edge-cases I find really annoying, tolerance thing.
    if ( k .lt. 48 ) then
        l=0
        il2: do i1=1,3
        do i2=1,3
        do i3=1,3
        do i4=1,3
        do i5=1,3
        do i6=1,3
        do i7=1,3
        do i8=1,3
        do i9=1,3
            ! possible operation in fractional coordinates
            op(:,1)=[fracvals(i1),fracvals(i2),fracvals(i3)]
            op(:,2)=[fracvals(i4),fracvals(i5),fracvals(i6)]
            op(:,3)=[fracvals(i7),fracvals(i8),fracvals(i9)]
            ! Convert it to Cartesian coordinates
            m=matmul(basis,op)
            m=matmul(m,inversebasis)
            ! Check that op*op^T=op*op^-1
            if ( lo_frobnorm( matmul( m,transpose(m) )-ident ) .lt. tolerance ) then
                l=l+1
                k=k+1
                bflist(:,:,k)=m
            endif
            ! done when I hit 48
            if ( l .eq. 48 ) exit il2
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo il2
    endif

    ! Keep the unique operations
    call lo_return_unique(bflist(:,:,1:k),list,tolerance)
    k=size(list,3) ! number of operations

    ! Make the values clean and nice
    list=lo_chop(list,tolerance*1E-2_r8)
    ! Get the fractional list
    fraclist=0.0_r8
    do i=1,k
        fraclist(:,:,i)=matmul(matmul(inversebasis,list(:,:,i)),basis)
    enddo
    fraclist=lo_chop(fraclist,tolerance*1E-2_r8)

    ! Store the ones that are ok, but put the identity first!
    lo_allocate(ops(3,3,k))
    lo_allocate(fracops(3,3,k))
    j=0
    do i=1,k
        if ( lo_frobnorm(list(:,:,i)-ident) .lt. tolerance ) then
            j=i
            ops(:,:,1)=list(:,:,i)
            fracops(:,:,1)=fraclist(:,:,i)
        endif
    enddo
    ! and then the others
    l=1
    do i=1,k
        if ( i .ne. j ) then
            l=l+1
            ops(:,:,l)=list(:,:,i)
            fracops(:,:,l)=fraclist(:,:,i)
        endif
    enddo
    ! try to put inversion in second place?
    do i=1,k
        if ( norm2(ops(:,:,i)+ident) .lt. lo_sqtol ) then
            m=ops(:,:,2)
            ops(:,:,2)=-ident
            ops(:,:,i)=m
        endif
    enddo
    deallocate(list)
end subroutine

!> Classify the operations, so that I know what is a rotation, a reflection and so on. Right now only classifies the rotation part.
subroutine figure_out_what_operations_are(sym)
    !> the operations
    type(lo_symset), intent(inout) :: sym
    !
    integer :: i,j,k,o
    integer :: fold
    real(r8), dimension(3,3) :: m,m1,identity
    real(r8), dimension(3) :: vec,v0,v1
    real(r8), dimension(3) :: welldefinedvalues
    real(r8) :: tr,eta,theta,sintheta,costheta,f0
    character(len=100) :: whatkind

    ! an identity matrix
    identity=0.0_r8
    identity(1,1)=1.0_r8
    identity(2,2)=1.0_r8
    identity(3,3)=1.0_r8

    ! Should add some more here.
    welldefinedvalues(1)=0.0_r8
    welldefinedvalues(2)=1.0_r8
    welldefinedvalues(3)=0.5_r8

    do o=1,sym%n
        ! get some basic things
        m=sym%op(o)%m           ! the operation
        eta=lo_determ(m)        ! the determinant
        tr=m(1,1)+m(2,2)+m(3,3) ! the trace
        ! The angle can only be pi,pi/2,pi/3,pi/4,pi/6 or zero
        f0=(tr-eta)*0.5_r8
        if ( abs(f0+1.0_r8) .lt. lo_tol ) f0=-1.0_r8
        if ( abs(f0-1.0_r8) .lt. lo_tol ) f0=1.0_r8
        f0=acos( f0 )
        if ( abs(f0) .lt. lo_radiantol ) then
            theta=0.0_r8
            costheta=cos(theta)
            sintheta=sin(theta)
            fold=0
        elseif ( abs(f0-lo_twopi/2.0_r8) .lt. lo_radiantol ) then
            theta=lo_twopi/2.0_r8
            costheta=cos(theta)
            sintheta=sin(theta)
            fold=2
        elseif ( abs(f0-lo_twopi/3.0_r8) .lt. lo_radiantol ) then
            theta=lo_twopi/3.0_r8
            costheta=cos(theta)
            sintheta=sin(theta)
            fold=3
        elseif ( abs(f0-lo_twopi/4.0_r8) .lt. lo_radiantol ) then
            theta=lo_twopi/4.0_r8
            costheta=cos(theta)
            sintheta=sin(theta)
            fold=4
        elseif ( abs(f0-lo_twopi/6.0_r8) .lt. lo_radiantol ) then
            theta=lo_twopi/6.0_r8
            costheta=cos(theta)
            sintheta=sin(theta)
            fold=6
        else
            ! small sanity check
            call lo_stop_gracefully(['Strange angle for a rotation matrix: '//tochar(theta*180/lo_pi)],lo_exitcode_symmetry)
        endif

        vec=0.0_r8
        whatkind='dunno'
        ! Ok, got nice angles, try to classify the operations. First the trivial.
        if ( lo_frobnorm(m-identity) .lt. lo_tol ) then
            whatkind='identity'
        elseif ( lo_frobnorm(m+identity) .lt. lo_tol ) then
            whatkind='inversion'
        elseif ( abs(eta*tr+1.0_r8) .gt. lo_tol .and. abs(eta*tr-3.0_r8) .gt. lo_tol ) then
            ! Easy cases:
            f0=sqrt( (3.0_r8-eta*tr)*(1.0_r8+eta*tr) )
            vec(1)=( m(3,2)-m(2,3) )/f0
            vec(2)=( m(1,3)-m(3,1) )/f0
            vec(3)=( m(2,1)-m(1,2) )/f0
            !
            if ( eta .gt. 0.0_r8 ) then
                ! easy, pure rotation
                whatkind='rotation'
            else
                ! pure reflection or rotoreflection
                if ( abs(theta) .lt. lo_radiantol ) then
                    ! I don't think this can happen
                    whatkind='reflection'
                else
                    ! I think this always happens.
                    whatkind='rotoreflection'
                endif
            endif
        else
            ! Trickier cases
            do i=1,3
                f0=m(i,i)+m(i,i)+eta-tr
                v0(i)=f0/( 3.0_r8*eta-tr )
                ! check for well defined values. Just to be on the safe side.
                do j=1,size(welldefinedvalues,1)
                    if ( abs(v0(i)-welldefinedvalues(j)) .lt. lo_sqtol ) then
                        v0(i)=welldefinedvalues(j)
                    endif
                enddo
                v0(i)=sqrt(v0(i))
            enddo
            ! Figure out the correct signs.
            signloop: do i=-1,1,2
            do j=-1,1,2
            do k=-1,1,2
                !
                v1(1)=i*v0(1)
                v1(2)=j*v0(2)
                v1(3)=k*v0(3)
                !
                if ( eta .gt. 0.0_r8 ) then
                    m1=lo_rotation_matrix_from_axis_and_angle(v1,theta)
                else
                    m1=lo_improper_rotation_matrix_from_axis_and_angle(v1,theta)
                endif
                !
                f0=lo_frobnorm(m1-m)
                if ( f0 .lt. lo_tol ) then
                    ! got it
                    vec=v1
                    exit signloop
                endif
            enddo
            enddo
            enddo signloop
            !
            if ( eta .gt. 0.0_r8 ) then
                whatkind='rotation'
            else
                if ( abs(theta) .lt. lo_radiantol ) then
                    ! no angle, then it's a pure reflection I suppose
                    whatkind='reflection'
                else
                    whatkind='rotoreflection'
                endif
            endif
            !
        endif

        ! Summarize the info:
        sym%op(o)%whatkind=trim(whatkind)
        sym%op(o)%fold=fold
        sym%op(o)%axis=vec
    enddo
end subroutine

!> Find the primitive cell
subroutine find_primitive_cell(basis,r,species,flavor,primbasis,tolerance)
    !> cell to be reduced
    real(r8), dimension(3,3), intent(in) :: basis
    !> atoms in the cell
    real(r8), dimension(:,:), intent(in) :: r
    !> species classification
    integer, dimension(:), intent(in) :: species
    !> flavor classification
    integer, dimension(:), intent(in) :: flavor
    !> primitive lattice vectors
    real(r8), dimension(3,3), intent(out) :: primbasis
    !> tolerance
    real(r8), intent(in) :: tolerance
    !
    real(r8), dimension(:,:), allocatable :: r0,r1
    real(r8), dimension(3,3) :: basis0,basis1
    integer, dimension(:), allocatable :: species0,flavor0,species1,flavor1
    integer :: i,na
    !
    na=size(r,2)
    allocate(r0(3,na))
    allocate(species0(na))
    allocate(flavor0(na))
    basis0=basis
    r0=r
    species0=species
    flavor0=flavor

    ! Recursively find a smaller cell
    redloop: do
        ! Try to find a smaller cell
        call reduce_cell(basis0,r0,species0,flavor0,basis1,r1,species1,flavor1,tolerance)
        if ( abs(abs(lo_determ(basis0))-abs(lo_determ(basis1))) .lt. tolerance ) then
            ! this means no new cell, we are converged.
            exit redloop
        else
            ! this means a new cell, update and try again!
            if ( allocated(r0) )       deallocate(r0)
            if ( allocated(species0) ) deallocate(species0)
            if ( allocated(flavor0) )  deallocate(flavor0)
            i=size(r1,2)
            allocate(r0(3,i))
            allocate(species0(i))
            allocate(flavor0(i))
            basis0=basis1
            r0=r1
            species0=species1
            flavor0=flavor1
            if ( allocated(r1) )       deallocate(r1)
            if ( allocated(species1) ) deallocate(species1)
            if ( allocated(flavor1) )  deallocate(flavor1)
        endif
    enddo redloop

    ! Here we should be done.
    primbasis=basis0

    if ( allocated(r0) )       deallocate(r0)
    if ( allocated(r1) )       deallocate(r1)
    if ( allocated(species0) ) deallocate(species0)
    if ( allocated(species1) ) deallocate(species1)
    if ( allocated(flavor0) )  deallocate(flavor0)
    if ( allocated(flavor1) )  deallocate(flavor1)
end subroutine

!> single iteration of the cell reduction
subroutine reduce_cell(latticevectors,r,species,flavor,newlatticevectors,newr,newspecies,newflavor,tolerance)
    !> original lattice vectors
    real(r8), dimension(3,3), intent(in) :: latticevectors
    !> original positions
    real(r8), dimension(:,:), intent(in) :: r
    !> original species
    integer, dimension(:), intent(in) :: species
    !> additional flavor for symmetry thing
    integer, dimension(:), intent(in) :: flavor
    !> tolerance
    real(r8), intent(in) :: tolerance

    real(r8), dimension(3,3), intent(out) :: newlatticevectors
    real(r8), dimension(:,:), allocatable, intent(out) :: newr
    integer, dimension(:), allocatable, intent(out) :: newspecies
    integer, dimension(:), allocatable,intent(out) :: newflavor

    real(r8), dimension(:,:), allocatable :: vecs,dum
    real(r8), dimension(3,3) :: m0,m1
    real(r8) :: vol,f0,tol,reltol,frtol
    integer :: i,j,k,l,nvecs,na
    logical :: cellvalid

    ! set the tolerance
    call lo_fetch_tolerance(tolerance,latticevectors,realspace_cart_tol=tol,realspace_fractional_tol=frtol,relative_tol=reltol)

    na=size(r,2)
    ! Get a list of possible lattice vectors: All the vectors in the cell + old latticevectors
    l=na+3
    allocate(dum(3,na+3))
    dum(:,1:na)=lo_clean_fractional_coordinates(r,frtol**2)
    do i=1,na
        dum(:,i)=matmul(latticevectors,dum(:,i))
    enddo
    do i=1,3
        dum(:,i+na)=latticevectors(:,i)
    enddo
    ! Reduce these to the unique. Should not be necessary, but never hurts.
    call lo_return_unique(dum,vecs,frtol)
    nvecs=size(vecs,2)

    ! volume of original cell
    vol=abs(lo_determ(latticevectors))
    ! Loop over all vectors to try to find a new cell.
    vl1: do i=1,nvecs
    vl2: do j=i+1,nvecs
    vl3: do k=j+1,nvecs
        ! possible new set of lattice vectors
        m0(:,1)=vecs(:,i)
        m0(:,2)=vecs(:,j)
        m0(:,3)=vecs(:,k)
        ! some quick tests to see if there is any point in continuing:
        f0=abs(lo_determ(m0))
        if ( f0 .lt. tol ) cycle vl3                           ! stupid cell if volume is zero
        if ( abs(mod(vol/f0,1.0_r8)) .gt. reltol ) cycle vl3 ! has to be an integer division
        if ( int(anint(vol/f0)) .eq. 1 ) cycle vl3             ! has to make the cell smaller, not the same
        ! I can get really retarded cells from this, do the Delaunay thing on it
        call reduce_basis_to_something_decent(m0,m1)
        ! Now do an actual test
        call test_cell(latticevectors,r,species,flavor,m1,cellvalid,newr,newspecies,newflavor,tolerance)
        if ( cellvalid ) then
            ! If I find a new cell, we will exit this routine and start over. Makes things way way faster.
            newlatticevectors=m1
            deallocate(vecs)
            return
        endif
    enddo vl3
    enddo vl2
    enddo vl1
    ! If we made it here, the cell is already the smallest. Return the same basis again.
    newlatticevectors=latticevectors
end subroutine

!> reduce a set of lattice vectors to canonical form
pure subroutine reduce_basis_to_something_decent(basis,canonical_basis)
    !> original lattice vectors
    real(r8), dimension(3,3), intent(in) :: basis
    !> canonical lattice vectors
    real(r8), dimension(3,3), intent(out) :: canonical_basis

    real(r8), dimension(3,4) :: ext_basis
    real(r8), dimension(3,7) :: dum
    real(r8), dimension(7) :: d
    real(r8) :: f0
    integer, dimension(7) :: ind
    integer :: i,j,k

    ! First do the Delaunay reduction thing
    ext_basis=0.0_r8
    ext_basis(:,1:3)=basis
    do i=1,3
        ext_basis(:,4)=ext_basis(:,4)-basis(:,i)
    enddo
    ! Iterate and stuff
    iterloop: do
        do i=1,4
        do j=i+1,4
            ! check scalar product
            f0=dot_product(ext_basis(:,i),ext_basis(:,j))
            if ( f0 .gt. lo_sqtol ) then
                ! subtract i from all not i and j
                do k=1,4
                    if ( k .ne. i .and. k .ne. j ) then
                        ext_basis(:,k)=ext_basis(:,k)+ext_basis(:,i)
                    endif
                enddo
                ! flip sign of i
                ext_basis(:,i)=-ext_basis(:,i)
                ! start over!
                cycle iterloop
            endif
        enddo
        enddo
        ! If we make it here, we are done!
        exit iterloop
    enddo iterloop

    ! Look through all variants to get shortest vectors
    dum=0.0_r8
    dum(:,1:4)=ext_basis
    dum(:,5)=ext_basis(:,1)+ext_basis(:,2)
    dum(:,6)=ext_basis(:,2)+ext_basis(:,3)
    dum(:,7)=ext_basis(:,1)+ext_basis(:,3)
    do i=1,7
        d(i)=norm2(dum(:,i))
    enddo
    call qsort(d,ind)
    !
    do i=3,7
        canonical_basis(:,1)=dum(:,ind(1))
        canonical_basis(:,2)=dum(:,ind(2))
        canonical_basis(:,3)=dum(:,ind(i))
        if ( abs(lo_determ(canonical_basis)) .gt. lo_tol ) exit
    enddo
    if ( lo_determ(canonical_basis) .lt. 0.0_r8 ) canonical_basis=-canonical_basis
end subroutine

!> Arrays to reduce to operations to one, the product of two symmetry operations in a group is a third operation in the same group.
subroutine reductionlist(sym,a,tr1,tr2)
    !> the set of operations
    class(lo_symset), intent(in) :: sym
    !> the reduction arraj
    integer, dimension(:,:), allocatable, intent(out) :: a
    !> should the first operation be inverse?
    logical, intent(in), optional :: tr1
    !> should the second operation be inverse?
    logical, intent(in), optional :: tr2

    integer :: i,j,k
    real(r8), dimension(3,3) :: m,m0,m1

    lo_allocate(a(sym%n,sym%n))
    a=0
    do i=1,sym%n
    do j=1,sym%n

        if ( present(tr1) .and. tr1 ) then
            m0=transpose(sym%op(i)%m)
        else
            m0=sym%op(i)%m
        endif

        if ( present(tr2) .and. tr2 ) then
            m1=transpose(sym%op(j)%m)
        else
            m1=sym%op(j)%m
        endif

        m=matmul(m0,m1)
        do k=1,sym%n
            if ( lo_frobnorm(m-sym%op(k)%m) .lt. lo_tol ) then
                a(i,j)=k
                exit
            endif
        enddo
    enddo
    enddo
end subroutine

!> expand a 3x3 symmetry operation to a 9x9 operation, to operate on flattened tensors
pure function lo_expandoperation_pair(o) result(bigo)
    !> the original 3x3 operation
    real(r8), dimension(3,3), intent(in) :: o
    !> the resulting 9x9 operation
    real(r8), dimension(9,9) :: bigo
    !
    !integer :: i,j,ii,jj,k,l
    !
    !bigo=0.0_r8
    !k=0
    !do i=1,3
    !do ii=1,3
    !    k=k+1
    !    l=0
    !    do j=1,3
    !    do jj=1,3
    !        l=l+1
    !        bigo(k,l)=o(ii,jj)*o(i,j)
    !    enddo
    !    enddo
    !enddo
    !enddo
    bigo(1,1)=o(1,1)*o(1,1)
    bigo(1,2)=o(1,2)*o(1,1)
    bigo(1,3)=o(1,3)*o(1,1)
    bigo(1,4)=o(1,1)*o(1,2)
    bigo(1,5)=o(1,2)*o(1,2)
    bigo(1,6)=o(1,3)*o(1,2)
    bigo(1,7)=o(1,1)*o(1,3)
    bigo(1,8)=o(1,2)*o(1,3)
    bigo(1,9)=o(1,3)*o(1,3)
    bigo(2,1)=o(2,1)*o(1,1)
    bigo(2,2)=o(2,2)*o(1,1)
    bigo(2,3)=o(2,3)*o(1,1)
    bigo(2,4)=o(2,1)*o(1,2)
    bigo(2,5)=o(2,2)*o(1,2)
    bigo(2,6)=o(2,3)*o(1,2)
    bigo(2,7)=o(2,1)*o(1,3)
    bigo(2,8)=o(2,2)*o(1,3)
    bigo(2,9)=o(2,3)*o(1,3)
    bigo(3,1)=o(3,1)*o(1,1)
    bigo(3,2)=o(3,2)*o(1,1)
    bigo(3,3)=o(3,3)*o(1,1)
    bigo(3,4)=o(3,1)*o(1,2)
    bigo(3,5)=o(3,2)*o(1,2)
    bigo(3,6)=o(3,3)*o(1,2)
    bigo(3,7)=o(3,1)*o(1,3)
    bigo(3,8)=o(3,2)*o(1,3)
    bigo(3,9)=o(3,3)*o(1,3)
    bigo(4,1)=o(1,1)*o(2,1)
    bigo(4,2)=o(1,2)*o(2,1)
    bigo(4,3)=o(1,3)*o(2,1)
    bigo(4,4)=o(1,1)*o(2,2)
    bigo(4,5)=o(1,2)*o(2,2)
    bigo(4,6)=o(1,3)*o(2,2)
    bigo(4,7)=o(1,1)*o(2,3)
    bigo(4,8)=o(1,2)*o(2,3)
    bigo(4,9)=o(1,3)*o(2,3)
    bigo(5,1)=o(2,1)*o(2,1)
    bigo(5,2)=o(2,2)*o(2,1)
    bigo(5,3)=o(2,3)*o(2,1)
    bigo(5,4)=o(2,1)*o(2,2)
    bigo(5,5)=o(2,2)*o(2,2)
    bigo(5,6)=o(2,3)*o(2,2)
    bigo(5,7)=o(2,1)*o(2,3)
    bigo(5,8)=o(2,2)*o(2,3)
    bigo(5,9)=o(2,3)*o(2,3)
    bigo(6,1)=o(3,1)*o(2,1)
    bigo(6,2)=o(3,2)*o(2,1)
    bigo(6,3)=o(3,3)*o(2,1)
    bigo(6,4)=o(3,1)*o(2,2)
    bigo(6,5)=o(3,2)*o(2,2)
    bigo(6,6)=o(3,3)*o(2,2)
    bigo(6,7)=o(3,1)*o(2,3)
    bigo(6,8)=o(3,2)*o(2,3)
    bigo(6,9)=o(3,3)*o(2,3)
    bigo(7,1)=o(1,1)*o(3,1)
    bigo(7,2)=o(1,2)*o(3,1)
    bigo(7,3)=o(1,3)*o(3,1)
    bigo(7,4)=o(1,1)*o(3,2)
    bigo(7,5)=o(1,2)*o(3,2)
    bigo(7,6)=o(1,3)*o(3,2)
    bigo(7,7)=o(1,1)*o(3,3)
    bigo(7,8)=o(1,2)*o(3,3)
    bigo(7,9)=o(1,3)*o(3,3)
    bigo(8,1)=o(2,1)*o(3,1)
    bigo(8,2)=o(2,2)*o(3,1)
    bigo(8,3)=o(2,3)*o(3,1)
    bigo(8,4)=o(2,1)*o(3,2)
    bigo(8,5)=o(2,2)*o(3,2)
    bigo(8,6)=o(2,3)*o(3,2)
    bigo(8,7)=o(2,1)*o(3,3)
    bigo(8,8)=o(2,2)*o(3,3)
    bigo(8,9)=o(2,3)*o(3,3)
    bigo(9,1)=o(3,1)*o(3,1)
    bigo(9,2)=o(3,2)*o(3,1)
    bigo(9,3)=o(3,3)*o(3,1)
    bigo(9,4)=o(3,1)*o(3,2)
    bigo(9,5)=o(3,2)*o(3,2)
    bigo(9,6)=o(3,3)*o(3,2)
    bigo(9,7)=o(3,1)*o(3,3)
    bigo(9,8)=o(3,2)*o(3,3)
    bigo(9,9)=o(3,3)*o(3,3)
    bigo=lo_chop(bigo,1E-11_r8)
end function

!> expand a 3x3 symmetry operation to a 27x27 operation, to operate on vector format of tensors
pure function lo_expandoperation_triplet(o) result(bigo)
    !> the original 3x3 operation
    real(r8), dimension(3,3), intent(in) :: o
    !> the resulting 27x27 operation
    real(r8), dimension(27,27) :: bigo

    integer :: i,ii,iii,j,jj,jjj,l,m

    bigo=0.0_r8
    l=0
    do i=1,3
    do ii=1,3
    do iii=1,3
        l=l+1
        m=0
        do j=1,3
        do jj=1,3
        do jjj=1,3
            m=m+1
            bigo(l,m)=o(iii,jjj)*o(ii,jj)*o(i,j)
        enddo
        enddo
        enddo
    enddo
    enddo
    enddo
    bigo=lo_chop(bigo,1E-11_r8)
end function

!> expand a 3x3 symmetry operation to a 81x81 operation, to operate on vector format of tensors
pure function lo_expandoperation_quartet(o) result(bigo)
    !> the original 3x3 operation
    real(r8), dimension(3,3), intent(in) :: o
    !> the resulting 27x27 operation
    real(r8), dimension(81,81) :: bigo

    integer :: i,ii,iii,iiii,j,jj,jjj,jjjj,l,m
    bigo=0.0_r8
    l=0
    do i=1,3
    do ii=1,3
    do iii=1,3
    do iiii=1,3
        l=l+1
        m=0
        do j=1,3
        do jj=1,3
        do jjj=1,3
        do jjjj=1,3
            m=m+1
            bigo(l,m)=o(iiii,jjjj)*o(iii,jjj)*o(ii,jj)*o(i,j)
        enddo
        enddo
        enddo
        enddo
    enddo
    enddo
    enddo
    enddo
    bigo=lo_chop(bigo,1E-11_r8)
end function

!> Calculate the character table for the space group
subroutine get_character_table(sym,verbosity)
    !> symmetry operations
    class(lo_symset), intent(inout) :: sym
    !> talk a lot?
    integer, intent(in) :: verbosity

    !integer :: n_class
    !type(lo_symset_conjugacyclass), dimension(:), allocatable :: cc
    !type(lo_symset_irreducible_representation), dimension(:), allocatable :: ir
    integer, dimension(sym%n,sym%n) :: multiplication_table
    real(r8) :: timer,t0,t1

    init: block
        if ( verbosity .gt. 0 ) then
            t0=walltime()
            timer=t0
            write(*,*) ''
            write(*,*) 'CALCULATING CHARACTER TABLE'
        endif
    end block init

    ! Divide the operations into conjugacy classes
    getclasses: block
        real(r8), dimension(3,3) :: m0
        real(r8), dimension(3) :: v0,v1
        integer, dimension(:,:), allocatable :: conjugacies,hi
        integer, dimension(:), allocatable :: dc1,di,dj
        integer, dimension(sym%n,sym%n) :: dc0
        integer :: i,j,k,l,ctr

        ! First get the multiplication table
        multiplication_table=0
        do i=1,sym%n
        do j=1,sym%n
            ! How do you chain operations again?
            m0=matmul(sym%op(i)%fm,sym%op(j)%fm)
            ! And the translation
            v0=matmul(sym%op(i)%fm,sym%op(j)%ftr)+sym%op(i)%ftr
            v0=lo_clean_fractional_coordinates(v0)
            ctr=0
            do k=1,sym%n
                if ( sum(abs(m0-sym%op(k)%fm)) .lt. lo_tol ) then
                    v1=lo_clean_fractional_coordinates(v0-sym%op(k)%ftr+0.5_r8)-0.5_r8
                    if ( sum(abs(v1)) .lt. lo_tol ) then
                        l=k
                        ctr=ctr+1
                    endif
                endif
            enddo
            if ( ctr .eq. 1 ) then
                multiplication_table(i,j)=l
            else
                call lo_stop_gracefully(['The group is not a group, but I thought it was.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        enddo
        ! Now get the conjugacy guys
        dc0=0
        do i=1,sym%n
            do j=1,sym%n
                ! Build the conjugation: the operations, from right to left do
                ! v' = im2*(v-t2)
                ! v'' = m1*v'+t1
                ! v''' = m2*v''+t2
                ! Combined it becomes
                ! v''' = m2*(m1*im2*(v-t2)+t1)+t2
                ! v''' = m2*m1*im2*v - m2*m1*im2*t2 + m2*t1 + t2
                m0 = matmul(sym%op(i)%fm,sym%op(j)%ifm)
                m0 = matmul(sym%op(j)%fm,m0)
                v0 = -matmul(m0,sym%op(j)%ftr) + matmul(sym%op(j)%fm,sym%op(i)%ftr) + sym%op(j)%ftr
                v0 = lo_clean_fractional_coordinates(v0)
                l=0
                ctr=0
                do k=1,sym%n
                    if ( sum(abs(m0-sym%op(k)%fm)) .lt. lo_tol ) then
                        v1=lo_clean_fractional_coordinates(v0-sym%op(k)%ftr+0.5_r8)-0.5_r8
                        if ( sum(abs(v1)) .lt. lo_tol ) then
                            l=k
                            ctr=ctr+1
                        endif
                    endif
                enddo
                if ( ctr .eq. 1 ) then
                    dc0(j,i)=l
                else
                    call lo_stop_gracefully(['The group is not a group, but I thought it was.'],lo_exitcode_symmetry,__FILE__,__LINE__)
                endif
            enddo
            ! Sort it to make it easier to find the unique
            call qsort(dc0(:,i))
        enddo
        ! Then return the unique
        call lo_return_unique(dc0,conjugacies)
        ! I want them sorted by number of members, makes it neater I think.
        ! Just a few extra steps, no big deal.
        sym%n_class=size(conjugacies,2)
        allocate(di(sym%n_class))
        allocate(dj(sym%n_class))
        allocate(hi(sym%n,sym%n_class))
        di=0
        dj=0
        hi=0
        do i=1,sym%n_class
            call lo_return_unique(conjugacies(:,i),dc1)
            di(i)=size(dc1)     ! size of class
            hi(1:di(i),i)=dc1   ! members of class
            deallocate(dc1)
        enddo
        ! sort by size. I think that the identity element stays
        ! put as class 1, always, but better make sure.
        di(1)=0
        call qsort(di,dj)
        di(1)=1

        ! Maybe sort these by size in the future, but who cares
        allocate(sym%cc(sym%n_class))
        do i=1,sym%n_class
            sym%cc(i)%n_member = di(i)
            allocate(sym%cc(i)%member( sym%cc(i)%n_member ))
            sym%cc(i)%member = hi(1:sym%cc(i)%n_member,dj(i))
        enddo

        deallocate(di,dj)

        ! First sanity check?
        allocate(di(sym%n))
        di=-1
        do i=1,sym%n_class
        do j=1,sym%cc(i)%n_member
            k=sym%cc(i)%member(j)
            di(k)=di(k)+1
        enddo
        enddo
        if ( sum(abs(di)) .ne. 0 ) then
            call lo_stop_gracefully(['Inclomplete conjugacy classes.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        deallocate(di)
        ! Think of more sanity checks later.

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... got '//tochar(sym%n_class)//' conjugacy classes ('//tochar(t1-t0)//'s)'
            t0=t1
        endif
    end block getclasses

    getchar: block
        real(r8), dimension(:,:,:), allocatable :: cnrs0
        real(r8), dimension(:,:), allocatable :: left_eigenvectors,right_eigenvectors
        real(r8), dimension(:,:), allocatable :: inv_right_eigenvectors,dm0,dm1
        real(r8), dimension(:), allocatable :: v0
        real(r8) :: f0
        integer, dimension(:), allocatable :: di
        integer :: i,j,k,l,ii,jj,iop,jop,kop,ctr

        ! This builds the magic CNRS parameter that I found in a few papers.
        allocate(cnrs0(sym%n_class,sym%n_class,sym%n_class))
        cnrs0=0.0_r8
        do i=1,sym%n_class
            do j=1,sym%n_class
            do k=1,sym%n_class
                ctr=0
                do ii=1,sym%cc(i)%n_member
                do jj=1,sym%cc(j)%n_member
                    iop=sym%cc(i)%member(ii)
                    jop=sym%cc(j)%member(jj)
                    kop=sym%cc(k)%member(1)
                    if ( multiplication_table(iop,jop) .eq. kop ) ctr=ctr+1
                enddo
                enddo
                cnrs0(k,j,i)=ctr
            enddo
            enddo
        enddo

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... got CNRS ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! Get the joint eigenvectors for these matrices
        allocate(left_eigenvectors(sym%n_class,sym%n_class))
        allocate(right_eigenvectors(sym%n_class,sym%n_class))
        call fancy_joint_eigenvector_thing(cnrs0,left_eigenvectors,right_eigenvectors)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... got joint eigenvectors ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! Some temporary space for figuring stuff out.
        allocate(inv_right_eigenvectors(sym%n_class,sym%n_class))
        allocate(dm0(sym%n_class,sym%n_class))
        allocate(dm1(sym%n_class,sym%n_class))
        allocate(v0(sym%n_class))
        allocate(di(sym%n_class))
        dm0=0.0_r8
        dm1=0.0_r8
        v0=0.0_r8
        di=0
        call lo_invert_real_matrix(right_eigenvectors,inv_right_eigenvectors)

        ! Now I can start fiddling with the irreducible representations.
        ! I will not do it all at once right away, first it should be sorted
        ! by dimension. Makes it prettier, I think. Or something. First I
        ! get the class characters, always a good idea.
        do i=1,sym%n_class
            call lo_triplegemm(inv_right_eigenvectors,cnrs0(:,:,i),right_eigenvectors,dm0)
            do j=1,sym%n_class
                dm1(j,i)=dm0(j,j)
                dm0(j,j)=0.0_r8
            enddo
            ! Sanity check
            if ( sum(abs(dm0)) .gt. 1E-12_r8*sym%n_class*sym%n_class ) then
                call lo_stop_gracefully(['Failed making similarity transform'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        ! And it was transposed. Who knew. Not me at least.
        dm0=transpose(dm1)
        ! Ok, so now dm0 holds the class characters per irrep as (:,i)

        ! sym%n_irrep is always the same as sym%n_class, but it gets confusing when
        ! I loop which is which, so I try to clarify a little.
        sym%n_irrep=sym%n_class
        ! First get the dimesion thingy per irrep.
        do i=1,sym%n_irrep
            f0=0.0_r8
            do j=1,sym%n_class
                f0=f0+dm0(j,i)*dm0(j,i)/real(sym%cc(j)%n_member,r8)
            enddo
            ! Several sanity checks
            if ( abs(f0) .lt. 1E-10_r8 ) then
                call lo_stop_gracefully(['Dimension of irrep has to be nonzero'],lo_exitcode_symmetry,__FILE__,__LINE__)
            else
                f0=sqrt(real(sym%n,r8)/f0)
            endif
            if ( abs(f0-anint(f0)) .gt. 1E-9_r8 ) then
                call lo_stop_gracefully(['Dimension of irrep has to be an integer'],lo_exitcode_symmetry,__FILE__,__LINE__)
            else
                f0=anint(f0)
            endif
            ! If nothing failed, here is the dimension stored in v0
            v0(i)=f0
        enddo

        ! Next is to convert the class constant to character.
        do i=1,sym%n_irrep
            do j=1,sym%n_class
                f0=v0(i)/real(sym%cc(j)%n_member)
                dm1(j,i)=f0*dm0(j,i)
            enddo
        enddo

        ! Now find the one with just 1 in it, and make sure it ends up first
        do i=1,sym%n_irrep
            if ( sum(abs(dm1(:,i)-1)) .lt. 1E-10_r8 ) then
                v0(i)=0
            endif
        enddo
        ! Now sort it by dimension, and also round to integers if applicable.
        call qsort(v0,di)
        v0(1)=1
        do i=1,sym%n_irrep
            dm0(:,i)=dm1(:,di(i))
            do j=1,sym%n_class
                if ( abs(dm0(j,i)-anint(dm0(j,i))) .lt. 1E-8_r8 ) dm0(j,i)=anint(dm0(j,i))
            enddo
        enddo

        ! Store this.
        allocate(sym%ir(sym%n_irrep))
        do i=1,sym%n_irrep
            sym%ir(i)%dimension=int(anint(v0(i)))
            allocate(sym%ir(i)%character_class(sym%n_class))
            allocate(sym%ir(i)%character_element(sym%n))
            sym%ir(i)%character_class=dm0(:,i)
            do j=1,sym%n_class
            do k=1,sym%cc(j)%n_member
                l=sym%cc(j)%member(k)
                sym%ir(i)%character_element(l)=dm0(j,i)
            enddo
            enddo
        enddo

        ! Test orthoginality? Maybe a good idea.
        do i=1,sym%n_irrep
            do j=1,sym%n_irrep
                dm1(j,i)=dot_product(sym%ir(i)%character_element,sym%ir(j)%character_element)
            enddo
            dm1(i,i)=dm1(i,i)-sym%n
        enddo
        if ( sum(abs(dm1)) .gt. 1E-12_r8*sym%n_class*sym%n_class ) then
            call lo_stop_gracefully(['irreps not orthogonal'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... got character table ('//tochar(t1-t0)//'s)'
            ! t0=t1
            ! ! Fancy output
            ! str="+-"
            ! do i=1,sym%n_class
            !     di(i)=i
            !     str=trim(adjustl(str))//'----'
            ! enddo
            ! write(*,*) ''
            ! write(*,"(6X,"//tochar(sym%n_class)//"(2X,I2))") di
            ! write(*,"(5X,A)") trim(str)
            ! do i=1,sym%n_class
            !     write(*,"(2X,I2,' |',"//tochar(sym%n_class)//"(2X,I2))") i,int(dm0(:,i))
            ! enddo
            ! write(*,*) ''
            write(*,*) 'Done calculting character table ('//tochar(walltime()-timer)//'s)'
        endif
    end block getchar
end subroutine

!> This one is starting to turn a little nasty.
subroutine fancy_joint_eigenvector_thing(M,left_eigenvectors,right_eigenvectors)
    !> many matrices to simultaneuously diagonalize
    real(r8), dimension(:,:,:), intent(in) :: M
    !> left eigenvectors
    real(r8), dimension(:,:), intent(out) :: left_eigenvectors
    !> right eigenvectors
    real(r8), dimension(:,:), intent(out) :: right_eigenvectors

    complex(r8), dimension(:), allocatable :: eigenvalues
    real(r8) :: not_an_eigenvalue
    integer, dimension(:,:), allocatable :: group_member
    integer, dimension(:), allocatable :: group_ctr
    integer :: n,n_group

    !@TODO remove left eigenvectors? Don't really do anything at the moment.

    ! Make sure it's a solveable problem. I kind of assume it is not too large
    ! if too large it will take forever.
    init: block
        real(r8), dimension(:,:), allocatable :: mA,mB
        integer :: i,j
        ! Size of problem:
        n=size(M,1)

        ! Some sanity checks. This first one is actually not a problem, but it means
        ! something stupid was sent to this routine, and that is never good.
        if ( size(M,1) .ne. size(M,2) ) then
            call lo_stop_gracefully(['Can only solve n x n x n problems'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( size(M,1) .ne. size(M,3) ) then
            call lo_stop_gracefully(['Can only solve n x n x n problems'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! Check that all matrices commute.
        allocate(mA(n,n))
        allocate(mB(n,n))
        do i=1,n
        do j=1,n
            call dgemm('N','N',n,n,n,1.0_r8,M(:,:,i),n,M(:,:,j),n,0.0_r8,mA,n)
            call dgemm('N','N',n,n,n,1.0_r8,M(:,:,j),n,M(:,:,i),n,0.0_r8,mB,n)
            if ( sum(abs(mA-mB)) .gt. 1E-13_r8*n*n ) then
                call lo_stop_gracefully(['Matrices have to commute, all of them.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        enddo
        deallocate(mA)
        deallocate(mB)

        ! Space for eigenvalues
        allocate(eigenvalues(n))
        eigenvalues=0.0_r8
    end block init

    ! First thing to do is to pick a sensible place to start from
    ! which is the matrix with as low degeneracy as possible, I think.
    ! makes for an easier start and fewer iterations below.
    preiter: block
        integer :: i,j,mi
        mi=0
        j=lo_hugeint
        do i=1,n
            call lo_general_real_eigenvalues_eigenvectors(M(:,:,i),eigenvalues,left_eigenvectors,right_eigenvectors)
            call sort_and_reorder_eigenvectors(eigenvalues,left_eigenvectors,right_eigenvectors,n_group,group_ctr,group_member)
            if ( maxval(group_ctr) .lt. j ) then
                mi=i
                j=maxval(group_ctr)
            endif
        enddo
        call lo_general_real_eigenvalues_eigenvectors(M(:,:,mi),eigenvalues,left_eigenvectors,right_eigenvectors)
        call sort_and_reorder_eigenvectors(eigenvalues,left_eigenvectors,right_eigenvectors,n_group,group_ctr,group_member)

        ! Get a number that is definitely not an eigenvalue
        not_an_eigenvalue=ceiling(maxval(abs(eigenvalues)))+1
    end block preiter

    ! Now try to actually solve things
    slv: block
        integer, parameter :: maxiter=1000 ! ??? Should be plenty, I think.
        complex(r8), dimension(:), allocatable :: sm_egv
        real(r8), dimension(:,:), allocatable :: V0,V1,V0i,dm0,sm_left,sm_right
        integer, dimension(:,:), allocatable :: sm_member
        integer, dimension(:), allocatable :: sm_ctr
        integer :: iter,gi,n_vecs,i,j,ii,jj,mj
        integer :: n_sm_group

        iterloop: do iter=1,maxiter
            ! Check for convergence
            if ( n_group .eq. n ) then
                ! This means all degeneracies are resolved
                exit iterloop
            endif

            if ( iter .eq. maxiter ) then
                ! This should never ever happen. In theory.
                call lo_stop_gracefully(['Recursing too many times, not good.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif

            grloop: do gi=1,n_group
                ! skip non-degenerate eigenvalues, they don't need splitting
                if ( group_ctr(gi) .eq. 1 ) cycle

                ! Some temporary space
                n_vecs=group_ctr(gi)
                allocate(V0(n,n_vecs))
                allocate(V1(n,n_vecs))
                allocate(V0i(n_vecs,n))
                allocate(dm0(n_vecs,n_vecs))
                allocate(sm_egv(n_vecs))
                allocate(sm_left(n_vecs,n_vecs))
                allocate(sm_right(n_vecs,n_vecs))
                do i=1,n_vecs
                    V0(:,i)=right_eigenvectors(:,group_member(i,gi))
                    V1(:,i)=left_eigenvectors(:,group_member(i,gi))
                enddo
                call lo_real_pseudoinverse(V0,V0i)
                mj=0
                j=lo_hugeint
                ! Find a good matrix to project on.
                do i=1,n
                    call lo_triplegemm(V0i,M(:,:,i),V0,dm0)
                    ! empty matrices don't help.
                    if ( sum(abs(dm0)) .lt. lo_sqtol ) cycle
                    call lo_general_real_eigenvalues_eigenvectors(dm0,sm_egv,sm_left,sm_right)
                    call sort_and_reorder_eigenvectors(sm_egv,sm_left,sm_right,n_sm_group,sm_ctr,sm_member)
                    if ( maxval(sm_ctr) .lt. j ) then
                        mj=i
                        j=maxval(sm_ctr)
                    endif
                enddo
                if ( mj .eq. 0 ) then
                    ! This should be impossible.
                    call lo_stop_gracefully(['Could not split eigenspace, should never happen.'],lo_exitcode_symmetry,__FILE__,__LINE__)
                else
                    ! This should always happen.
                    call lo_triplegemm(V0i,M(:,:,mj),V0,dm0)
                    call lo_general_real_eigenvalues_eigenvectors(dm0,sm_egv,sm_left,sm_right)
                    call sort_and_reorder_eigenvectors(sm_egv,sm_left,sm_right,n_sm_group,sm_ctr,sm_member)
                endif
                ! Now we have a new, tiny group. First apply the rotation to the currect subspace
                ! to fix things about life in general, and the eigenvectors in particular.
                V0 = matmul(V0,sm_right)
                V1 = matmul(V1,sm_left)
                do i=1,n_vecs
                    right_eigenvectors(:,group_member(i,gi))=V0(:,i)
                    left_eigenvectors(:,group_member(i,gi))=V1(:,i)
                enddo
                ! Here is a magical thing. For each distinct eigenvalue, change the corresponding
                ! global eigenvalue to a new number. Then re-calculate the groups for the full
                ! thing, then cycle the outer loop. This will have the effect that every iterloop,
                ! one degenerate subspace is split. Eventually all of them get split.
                ! Kind of a recursive thing but I never get actual recursive things to work.
                do i=1,n_sm_group
                    not_an_eigenvalue=not_an_eigenvalue+1
                    do j=1,sm_ctr(i)
                        ii=sm_member(j,i)       ! Index in small group
                        jj=group_member(ii,gi)  ! Global index
                        eigenvalues(jj)=not_an_eigenvalue
                    enddo
                enddo

                ! Now recalculate the global groups
                call sort_and_reorder_eigenvectors(eigenvalues,left_eigenvectors,right_eigenvectors,n_group,group_ctr,group_member)
                ! Some cleanup
                deallocate(V0,V0i,V1,dm0,sm_egv,sm_left,sm_right)
                ! And start over!
                cycle iterloop
            enddo grloop
        enddo iterloop
    end block slv
end subroutine

!> Takes a set of eigenvectors and eigenvalues and sorts and groups them according to multiplicity.
subroutine sort_and_reorder_eigenvectors(eigenvalues,left_eigenvectors,right_eigenvectors,n_group,group_ctr,group_member)
    !> eigenvalues
    complex(r8), dimension(:), intent(inout) :: eigenvalues
    !> left and right eigenvectors. left are (:,i), right are also (:,i)
    real(r8), dimension(:,:), intent(inout) :: left_eigenvectors,right_eigenvectors
    !> how many groups of degenerate eigenvalues are there
    integer, intent(out) :: n_group
    !> counter for each group
    integer, dimension(:), allocatable, intent(out) :: group_ctr
    !> members of each group
    integer, dimension(:,:), allocatable, intent(out) :: group_member

    real(r8), dimension(size(eigenvalues),size(eigenvalues)) :: dma,dmb
    complex(r8), dimension(size(eigenvalues)) :: dc
    integer, dimension(size(eigenvalues),size(eigenvalues)) :: dj
    integer, dimension(size(eigenvalues)) :: mult,di
    integer :: i,j,n

    ! First sort everything by multiplicity to make life a little easier.
    ! Not completely sure this step is necessary, but cost almost nothing
    ! anyway, so I might as well.
    n=size(eigenvalues)
    mult=0
    do i=1,n
    do j=1,n
        if ( abs(eigenvalues(i)-eigenvalues(j)) .lt. lo_sqtol ) mult(i)=mult(i)+1
    enddo
    enddo
    call qsort(mult,di)
    dc=eigenvalues
    dma=left_eigenvectors
    dmb=right_eigenvectors
    do i=1,n
        eigenvalues(i)=dc(di(i))
        left_eigenvectors(:,i)=dma(:,di(i))
        right_eigenvectors(:,i)=dmb(:,di(i))
    enddo

    ! Ok, good start. Now sort all of this into groups of constant eigenvalue.
    ! Slightly ninja loop, build the groups and multiplicites on-the-fly, NlogN
    dc=maxval(abs(eigenvalues))*10
    n_group=0
    di=0
    dj=0
    ul: do i=1,n
        do j=1,n_group
            if ( abs(eigenvalues(i)-dc(j)) .lt. lo_sqtol ) then
                di(j)=di(j)+1
                dj( di(j),j )=i
                cycle ul
            endif
        enddo
        n_group=n_group+1
        dc(n_group)=eigenvalues(i)
        di(n_group)=di(n_group)+1
        dj( di(n_group),n_group )=i
    enddo ul

    allocate(group_ctr(n_group))
    allocate(group_member(maxval(di),n_group))
    group_ctr=di(1:n_group)
    do i=1,n_group
        group_member(1:group_ctr(i),i)=dj(1:group_ctr(i),i)
    enddo
end subroutine

!> Check if an operation is valid
pure subroutine check_operation(op,tr,sbox,basis,fmap,valid,tolerance)
    !> The point operation
    real(r8), dimension(3,3), intent(in) :: op
    !> The translation
    real(r8), dimension(3), intent(in) :: tr
    !> the points, decomposed into boxes per species
    type(lo_speciesbox), dimension(:), intent(in) :: sbox
    !> the basis
    real(r8), dimension(3,3), intent(in) :: basis
    !> the mappy stuff
    integer, dimension(:), intent(inout) :: fmap
    !> is it valid?
    logical, intent(out) :: valid
    !> tolerance
    real(r8), intent(in) :: tolerance
    !
    real(r8), dimension(:,:), allocatable :: dum
    real(r8), dimension(3) :: v0
    real(r8) :: sqtol
    integer :: i,j,sp,ii,jj
    logical :: match

    ! Set the tolerances
    call lo_fetch_tolerance(tolerance,basis,realspace_fractional_tol=sqtol,squared=.true.)

    valid=.false.
    fmap=0
    ! check per species separately, makes it faster
    jj=0
    do sp=1,size(sbox,1)
        ! get the rotated lattice
        allocate(dum(3,sbox(sp)%n))
        do i=1,sbox(sp)%n
            ! choose a particle
            v0=sbox(sp)%r(:,i)
            ! rotate it
            v0=matmul(op,v0)
            ! add translation
            v0=v0+tr
            ! move back to first unitcell
            dum(:,i)=lo_clean_fractional_coordinates(v0,sqtol)
        enddo

        ! now I have the rotated, compare with all
        ii=0
        do i=1,sbox(sp)%n
            match=.false.
            do j=1,sbox(sp)%n
                ! compare vectors
                v0=sbox(sp)%r(:,i)-dum(:,j)
                v0=lo_clean_fractional_coordinates(v0+0.5_r8,sqtol)-0.5_r8
                ! compare
                if ( lo_sqnorm(v0) .lt. sqtol ) then
                    match=.true.
                    fmap( sbox(sp)%ind(j) )=sbox(sp)%ind(i)
                    exit
                endif
            enddo
            ! If there was no match, the operation is invalid
            if ( match .eqv. .false. ) then
                valid=.false.
                fmap=0
                return
            else
                ! if there was a match, continue
                ii=ii+1
            endif
        enddo
        ! Where all points matched?
        if ( ii .eq. sbox(sp)%n ) then
            jj=jj+1
        endif
        deallocate(dum)
    enddo

    ! It might be valid!
    if ( jj .eq. size(sbox,1) ) then
        valid=.true.
    else
        valid=.false.
        fmap=0
    endif
end subroutine

!> test if a new cell is a smaller unit cell
subroutine test_cell(oldcell,oldr,oldspecies,oldflavor,cell,valid,newr,newspecies,newflavor,tolerance)
    !> original cell
    real(r8), dimension(3,3), intent(in) :: oldcell
    !> original positions
    real(r8), dimension(:,:), intent(in) :: oldr
    !> original species
    integer, dimension(:), intent(in) :: oldspecies
    !> original flavor
    integer, dimension(:), intent(in) :: oldflavor
    !> new cell
    real(r8), dimension(3,3), intent(in) :: cell
    !> is it a new cell?
    logical, intent(out) :: valid
    !> new positions in case of a valid cell
    real(r8), dimension(:,:), allocatable, intent(out) :: newr
    !> new species in case of a valid cell
    integer, dimension(:), allocatable, intent(out) :: newspecies
    !> new flavor in case of a valid cell
    integer, dimension(:), allocatable, intent(out) :: newflavor
    !> tolerance
    real(r8), intent(in) :: tolerance
    !
    real(r8), dimension(3,3) :: icell,m0
    real(r8), dimension(:,:), allocatable :: dumr1,dumr2
    real(r8) :: tol,sqtol,f0
    integer :: i,j,l,mult,oldna,newna

    ! set the tolerance
    call lo_fetch_tolerance(tolerance,oldcell,realspace_fractional_tol=tol)
    sqtol=tol**2

    valid=.false.
    ! How many atoms should I get of each?
    f0=abs(lo_determ(oldcell)/lo_determ(cell))
    mult=int(anint(abs(lo_determ(oldcell)/lo_determ(cell))))

    ! m0 is a matrix that converts to fractional coordinates in the new cell
    icell=lo_invert3x3matrix(cell)
    m0=matmul(icell,oldcell)

    ! convert coordinates to fractional with respect to the new cell. Also
    ! I decorate the positions with the species and flavor, to keep track of that.
    oldna=size(oldr,2)
    allocate(dumr1(5,oldna))
    do i=1,oldna
        dumr1(1:3,i)=lo_clean_fractional_coordinates(matmul(m0,oldr(:,i)),sqtol)
        dumr1(4:5,i)=[oldspecies(i),oldflavor(i)]
    enddo
    ! Reduce this to the union of the list, the unique
    call lo_return_unique(dumr1,dumr2,tol)

    ! Did it become the correct number of atoms?
    if ( abs(size(dumr2,2)*mult-oldna) .ne. 0 ) then
        valid=.false.
        return
    endif

    ! Slightly better check: Each of the atoms in the reduced cell must appear
    ! exactly mult times in the old cell.
    do i=1,size(dumr2,2)
        l=0
        do j=1,oldna
            if ( sum(abs(dumr1(:,j)-dumr2(:,i))) .lt. tol ) l=l+1
        enddo
        ! kill if not true
        if ( l .ne. mult ) then
            valid=.false.
            return
        endif
    enddo

    ! ok, now I say this is valid, return the new cell
    newna=size(dumr2,2)
    allocate(newr(3,newna))
    allocate(newspecies(newna))
    allocate(newflavor(newna))
    newr=dumr2(1:3,:)
    newspecies=int(anint(dumr2(4,:)))
    newflavor=int(anint(dumr2(5,:)))
    valid=.true.
end subroutine

!> Return the possible translations for this lattice
pure subroutine return_possible_translations(translations,fractranslations,ops,basis,inversebasis,points,species,tolerance)
    !> basis
    real(r8), dimension(3,3), intent(in) :: basis
    !> inverse of basis
    real(r8), dimension(3,3), intent(in) :: inversebasis
    !> translpoint operations
    real(r8), dimension(:,:), allocatable, intent(out) :: translations
    !> translations in fractional coordinates
    real(r8), dimension(:,:), allocatable, intent(out) :: fractranslations
    !> point operations
    real(r8), dimension(:,:,:), allocatable, intent(in) :: ops
    !> positions of atoms
    real(r8), dimension(:,:), intent(in) :: points
    !> species of the atoms
    integer, dimension(:), intent(in) :: species
    !> tolerance
    real(r8), intent(in) :: tolerance

    integer :: i,j,k,l,n,o
    real(r8), dimension(3) :: v0,v1
    real(r8), dimension(:,:), allocatable :: dum,r
    real(r8) :: fractol

    ! Set the tolerance
    call lo_fetch_tolerance(tolerance,basis,realspace_fractional_tol=fractol)
    fractol=fractol*100.0_r8

    ! I operate under the assumption that a translation has to move one
    ! atom to another. Why do I know that? Well, I am allowed to rigidly translate the
    ! lattice within the cell without changing the symmetry. So I can put any atom at
    ! the origin. It will not move anywhere by the rotation part, so the translation
    ! part has to move it to another site with the same kind of atom. This can be repeated
    ! for all kinds of atoms. So the possible translations are then all the pair vectors in
    ! the cell, that move an atom of species A.

    ! get cartesian positions
    n=size(points,2)
    lo_allocate(r(3,n))
    do i=1,n
        r(:,i)=matmul(basis,points(:,i))
    enddo
    ! First rough count, I only want translations between the same species
    l=0
    do i=1,n
    do j=1,n
        if ( species(i) .eq. species(j) ) l=l+1
    enddo
    enddo

    ! I consider the action of the rotations as well, so multiply with that.
    l=l*size(ops,3)
    lo_allocate(dum(3,l))
    dum=0.0_r8
    l=0
    do i=1,n
    do j=i,n
        if ( species(i) .eq. species(j) ) then
            do o=1,size(ops,3)
                ! possible vector
                v0=r(:,j)-matmul(ops(:,:,o),r(:,i))
                v0=lo_clean_fractional_coordinates(matmul(inversebasis,v0))
                do k=1,3
                    if ( v0(k) .gt. fractol ) then
                        v1(k)=1.0_r8/v0(k)
                    else
                        v1(k)=0.0_r8
                    endif
                enddo
                ! vector has to be a 1/n-thingy, where n is an integer
                if ( sum(abs(v1-anint(v1))) .lt. fractol ) then
                    l=l+1
                    dum(:,l)=r(:,j)-matmul(ops(:,:,o),r(:,i))
                endif

            enddo
        endif
    enddo
    enddo

    ! back to fractional
    do i=1,l
        dum(:,i)=matmul(inversebasis,dum(:,i))
    enddo

    ! Chop to my well-defined values, and make sure the translations are 0-1, not including 1
    dum(:,1:l)=lo_clean_fractional_coordinates(lo_chop(lo_clean_fractional_coordinates(dum(:,1:l)),fractol))

    ! Now keep only the unique translations, no point in having millions of them.
    call lo_return_unique(dum(:,1:l),fractranslations,fractol)
    ! sanity check

    lo_allocate(translations(3,size(fractranslations,2)))
    do i=1,size(fractranslations,2)
        translations(:,i)=matmul(basis,fractranslations(:,i))
    enddo
end subroutine

!> size in memory, approximately.
pure function size_in_mem(sym) result(mem)
    !> symmetry operations
    class(lo_symset), intent(in) :: sym
    !> size in memory, in bytes
    integer(i8) :: mem

    integer :: i
    mem=0
    mem=mem+storage_size(sym)
    if ( allocated(sym%degeneracy) ) mem=mem+size(sym%degeneracy)*storage_size(sym%degeneracy)
    if ( allocated(sym%degenerate_atom) ) mem=mem+size(sym%degenerate_atom)*storage_size(sym%degenerate_atom)
    if ( allocated(sym%irr_to_all) ) mem=mem+size(sym%irr_to_all)*storage_size(sym%irr_to_all)
    if ( allocated(sym%all_to_irr) ) mem=mem+size(sym%all_to_irr)*storage_size(sym%all_to_irr)
    if ( allocated(sym%op_all_from_irr) ) mem=mem+size(sym%op_all_from_irr)*storage_size(sym%op_all_from_irr)
    if ( allocated(sym%irr_unfold_ctr) ) mem=mem+size(sym%irr_unfold_ctr)*storage_size(sym%irr_unfold_ctr)
    if ( allocated(sym%irr_unfold_index) ) mem=mem+size(sym%irr_unfold_index)*storage_size(sym%irr_unfold_index)
    if ( allocated(sym%irr_unfold_operation) ) mem=mem+size(sym%irr_unfold_operation)*storage_size(sym%irr_unfold_operation)
    if ( allocated(sym%multiplication_table) ) mem=mem+size(sym%multiplication_table)*storage_size(sym%multiplication_table)
    if ( allocated(sym%op) ) then
        mem=mem+size(sym%op)*storage_size(sym%op)
        do i=1,size(sym%op)
            mem=mem+size(sym%op(i)%fmap)*storage_size(sym%op(i)%fmap)
        enddo
    endif
    if ( allocated(sym%cc) ) then
        mem=mem+size(sym%cc)*storage_size(sym%cc)
        do i=1,size(sym%cc)
            if ( allocated(sym%cc(i)%member) ) mem=mem+size(sym%cc(i)%member)*storage_size(sym%cc(i)%member)
        enddo
    endif
    if ( allocated(sym%ir) ) then
        mem=mem+sym%n*storage_size(sym%ir)
        do i=1,size(sym%ir)
            if ( allocated(sym%ir(i)%character_class) ) mem=mem+size(sym%ir(i)%character_class)*storage_size(sym%ir(i)%character_class)
            if ( allocated(sym%ir(i)%character_element) ) mem=mem+size(sym%ir(i)%character_element)*storage_size(sym%ir(i)%character_element)
        enddo
    endif
    mem=mem/8
end function

! Stuff below is not quite done, but might come handy some day.
!
! !> multiply operations, op3 = op1*op2. Dangerous to use for now, does not update everything to make op3 a complete operation.
! subroutine multiply_operations(op1,op2,op3,inv1,inv2)
!     type(lo_spacegroup_operation), intent(in) :: op1,op2
!     type(lo_spacegroup_operation), intent(out) :: op3
!     logical, intent(in) :: inv1,inv2
!
!     real(r8), dimension(3,3) :: m1,m2
!     real(r8), dimension(3) :: v0,v1,v2
!     integer :: whatkind
!
!     ! So, forward operation is
!     ! v' = m*v+t
!     ! Backward operation is
!     ! v' = im*(v-t)
!     ! So we get for ways of combining operations:
!     if ( inv1 ) then
!         if ( inv2 ) then
!             !whatkind=1
!             ! im2*(v-t2)
!             ! im1*(v-t1)
!             ! im1*(im2*(v-t2)-t1) = im1*im2*v - im1*im2*t2 - im1*t1
!             op3%fm = matmul(op1%ifm,op2%ifm)
!             op3%ftr = -matmul(op3%fm,op2%ftr)-matmul(op1%ifm,op1%ftr)
!         else
!             ! m2*v + t2
!             ! im1*(v-t1)
!             ! im1*(m2*v+t2-t1) = im1*m2*v + im1*(t2-t1)
!             op3%fm = matmul(op1%ifm,op2%fm)
!             op3%ftr = matmul(op1%ifm,op2%ftr-op1%ftr)
!         endif
!     else
!         if ( inv2 ) then
!             ! im2*(v - t2)
!             ! m1*im2*(v -t2) + t1 = m1*im2*v - m1*im2*t2 + t1
!             op3%fm = matmul(op1%fm,op2%ifm)
!             op3%ftr = -matmul(op3%fm,op2%ftr)+op1%tr
!
!         else
!             ! m2*v + t2
!             ! m1*m2*(v +t2) + t1 = m1*m2*v - m1*t2 + t1
!             op3%fm = matmul(op1%fm,op2%fm)
!             op3%ftr = matmul(op1%fm,op2%ftr)+op1%tr
!         endif
!     endif
!     ! Clean the translation part
!     op3%ftr = lo_clean_fractional_coordinates(op3%ftr)
! end subroutine
!
! !> conjugate operations op3 = op2 * op1 * op2^-1
! subroutine conjugate_operations(op1,op2,op3)
!     type(lo_spacegroup_operation), intent(in) :: op1,op2
!     type(lo_spacegroup_operation), intent(out) :: op3
!
!     ! The operations, from right to left
!     ! v = im2*(v-t2)
!     ! v = m1*v+t1
!     ! v = m2*v+t2
!     ! Combined it becomes
!     ! v = m2*(m1*im2*(v-t2)+t1)+t2
!     ! v = m2*m1*im2*v - m2*m1*im2*t2 + m2*t1 + t2
!     op3%fm = matmul(op1%fm,op2%ifm)
!     op3%fm = matmul(op2%fm,op3%fm)
!     op3%ftr = -matmul(op3%fm,op2%ftr) + matmul(op2%fm,op1%ftr) + op2%ftr
!     op3%ftr = lo_clean_fractional_coordinates(op3%ftr)
! end subroutine
!
! function compare_operations(op1,op2) result(match)
!     type(lo_spacegroup_operation), intent(in) :: op1,op2
!     logical :: match
!
!     real(r8), dimension(3) :: v0
!     if ( sum(abs(op1%fm-op2%fm)) .gt. lo_sqtol ) then
!         match=.false.
!         return
!     endif
!
!     v0=lo_clean_fractional_coordinates(op1%ftr-op2%ftr+0.5_r8)-0.5_r8
!     if ( sum(abs(v0)) .gt. lo_tol ) then
!         match=.false.
!         return
!     endif
!
!     match=.true.
! end function



end module
