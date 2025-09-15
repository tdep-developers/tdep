submodule (lo_symmetry_of_interactions) lo_symmetry_of_interactions_helpers
use gottochblandat, only: lo_chop,lo_determ,lo_sqnorm,lo_return_tensor_transpositions
use type_blas_lapack_wrappers, only: lo_gemm
implicit none
contains

!> check that cutoffs are sane
module subroutine check_cutoff(ss,cutoff,tol,dt,wraparound,mw)
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> cutoff to check
    real(r8), intent(inout) :: cutoff
    !> tolerance for degeneracy
    real(r8), intent(in) :: tol
    !> unitcell distance table
    type(lo_distancetable), intent(in) :: dt
    !> are wraparound distances allowed?
    logical, intent(in) :: wraparound
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw

    integer :: iatom,ipair,ctr

    ! If the cutoff is too small, do nothing
    if ( cutoff .lt. 0.0_r8 ) then
        return
    endif
    ! make sure it's not too long
    if ( .not.wraparound ) cutoff=min(ss%maxcutoff(),cutoff)
    ! and not too short
    cutoff=max(ss%mincutoff(),cutoff)
    ! Check if the cutoffs coincide with pair distances, that causes instabilities.
    do
        ctr=0
        do iatom=1,dt%np
        do ipair=1,dt%particle(iatom)%n
            if ( abs(cutoff-dt%particle(iatom)%d(ipair)) .lt. tol ) ctr=ctr+1
        enddo
        enddo
        if ( ctr .eq. 0 ) then
            exit
        else
            cutoff=cutoff+4*tol
            ! if ( mw%talk ) then
            !     write(lo_iou,*) 'WARNING: cutoff coincides with a coordination shell, increasing it:'
            !     write(lo_iou,*) '    from:',cutoff
            !     cutoff=cutoff+4*tol
            !     write(lo_iou,*) '      to:',cutoff
            ! endif
        endif
    enddo
end subroutine

!> Align distance tables
module subroutine align_distance_tables(dtuc,dtss,ss,tol,mw,mem)
    !> unitcell distance table
    class(lo_distancetable), intent(in) :: dtuc
    !> supercell distance table
    class(lo_distancetable), intent(inout) :: dtss
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> tolerance
    real(r8), intent(in) :: tol
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(3) :: v0
    real(r8) :: sqtol
    integer, dimension(:), allocatable :: di
    integer :: iatom,jatom,i,j

    ! A smart thing to do is to make sure that the distance tables are sorted identically, such that
    ! for each atom in the supercell (i) that correspond to an atom in the unitcell (j), then neighbour
    ! n from both (i) and (j) are the same. This makes matching the tuplets from the prototypes to the
    ! supercell fast, and way easier. Also, this is an early and good place for some sanity checks.
    sqtol=tol**2
    do iatom=1,dtss%np
        jatom=ss%info%index_in_unitcell(iatom)
        if ( dtuc%particle(jatom)%n .ne. dtss%particle(iatom)%n ) then
            call lo_stop_gracefully(['Inconsistent number of neighbours'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        call mem%allocate(di,dtuc%particle(jatom)%n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=0
        p1l: do i=1,dtss%particle(iatom)%n
            do j=1,dtuc%particle(jatom)%n
                ! same distance
                if ( abs(dtuc%particle(jatom)%d(j)-dtss%particle(iatom)%d(i)) .gt. tol ) cycle
                ! same vectors
                v0=dtuc%particle(jatom)%v(:,j)-dtss%particle(iatom)%v(:,i)
                if ( lo_sqnorm( v0 ) .lt. sqtol ) then
                    di(j)=i
                    cycle p1l
                endif
            enddo
            ! if we make it here, there is a discrepancy in the distance tables
            call lo_stop_gracefully(['Inconsistency in distance tables.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        enddo p1l
        ! now I just reorder everything
        dtss%particle(iatom)%ind=dtss%particle(iatom)%ind(di)
        dtss%particle(iatom)%d  =dtss%particle(iatom)%d(di)
        dtss%particle(iatom)%v  =dtss%particle(iatom)%v(:,di)
        dtss%particle(iatom)%lv =dtss%particle(iatom)%lv(:,di)
        if ( allocated(dtss%particle(iatom)%weight) ) then
            dtss%particle(iatom)%weight=dtss%particle(iatom)%weight(di)
        endif
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    enddo
end subroutine

!> expand operations to include transpositions as well
module subroutine expandoperations(sl,uc,spacegroup,transposition,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> consider spacegroup symmetries
    logical, intent(in) :: spacegroup
    !> consider transposition symmetries
    logical, intent(in) :: transposition
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:,:), allocatable :: trm_pair,trm_triplet,trm_quartet
    integer, dimension(:,:), allocatable :: prm_pair,prm_triplet,prm_quartet
    integer :: nop

    init: block
        integer :: o
        ! Make space for the transpositions. I had to make these allocatable since I ran inte weird weird bugs
        call mem%allocate(trm_pair   ,[9,9,2],   persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(trm_triplet,[27,27,6], persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(trm_quartet,[81,81,24],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(prm_pair   ,[2,2],     persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(prm_triplet,[3,6],     persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(prm_quartet,[4,24],    persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        trm_pair=0.0_r8
        trm_triplet=0.0_r8
        trm_quartet=0.0_r8
        prm_pair=0
        prm_triplet=0
        prm_quartet=0

        ! Get all the transpositions
        call lo_return_tensor_transpositions(trm_pair,trm_triplet,trm_quartet,prm_pair,prm_triplet,prm_quartet)
        ! start by figuring out if I should consider spacegroup/transpositions or not.
        if ( spacegroup ) then
            nop=uc%sym%n ! all operations
        else
            nop=1 ! only identity
        endif

        ! Store the singlet operations, never hurts
        sl%n_singlet_operation=nop
        allocate(sl%singletop(nop))
        do o=1,nop
            sl%singletop(o)%permind=-1
            sl%singletop(o)%opind=o
            sl%singletop(o)%detop=int(anint(lo_determ(uc%sym%op(o)%m)))
            sl%singletop(o)%m3=uc%sym%op(o)%m
            sl%singletop(o)%im3=uc%sym%op(o)%im
            allocate(sl%singletop(o)%fmap(uc%na))
            sl%singletop(o)%fmap=uc%sym%op(o)%fmap
        enddo
    end block init

    fixpair: block
        real(r8), dimension(:,:,:), allocatable :: so9
        integer :: l,t,o,ntr

        if ( transposition ) then
            ntr=size(prm_pair,2) ! all permutations
        else
            ntr=1 ! only identity
        endif

        call mem%allocate(so9,[9,9,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        so9=0.0_r8
        do o=1,nop
            so9(:,:,o)=expandoperation_pair(uc%sym%op(o)%m)
        enddo
        ! count and make space
        sl%n_pair_operation=nop*ntr
        allocate(sl%pairop(sl%n_pair_operation))
        ! populate
        l=0
        do t=1,ntr
        do o=1,nop
            l=l+1
            ! Unclear why minus sign. Put past Olle put it there, and
            ! I won't change it to cause trouble for future Olle.
            if ( t .eq. 1 ) then
                sl%pairop(l)%m3=uc%sym%op(o)%m
                sl%pairop(l)%im3=uc%sym%op(o)%im
            else
                sl%pairop(l)%m3=-uc%sym%op(o)%m
                sl%pairop(l)%im3=-uc%sym%op(o)%im
            endif
            call lo_gemm(so9(:,:,o),trm_pair(:,:,t),sl%pairop(l)%sotr)
            sl%pairop(l)%perm=prm_pair(:,t)
            sl%pairop(l)%permind=t
            sl%pairop(l)%opind=o
            sl%pairop(l)%detop=int(anint(lo_determ(uc%sym%op(o)%m)))
            allocate(sl%pairop(l)%fmap(uc%na))
            sl%pairop(l)%fmap=uc%sym%op(o)%fmap
        enddo
        enddo
        call mem%deallocate(so9,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block fixpair

    ! Now the same thing for the triplets
    fixtriplet: block
        real(r8), dimension(:,:,:), allocatable :: so27
        integer :: l,t,o,ntr

        if ( transposition ) then
            ntr=size(prm_triplet,2) ! all permutations
        else
            ntr=1 ! only identity
        endif

        call mem%allocate(so27,[27,27,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        so27=0.0_r8
        do o=1,nop
            so27(:,:,o)=expandoperation_triplet(uc%sym%op(o)%m)
        enddo
        ! count and make space
        sl%n_triplet_operation=nop*ntr
        allocate(sl%tripletop(sl%n_triplet_operation))
        ! populate
        l=0
        do t=1,ntr
        do o=1,nop
            l=l+1
            sl%tripletop(l)%m3=uc%sym%op(o)%m
            sl%tripletop(l)%im3=uc%sym%op(o)%im
            call lo_gemm(so27(:,:,o),trm_triplet(:,:,t),sl%tripletop(l)%sotr)
            sl%tripletop(l)%perm=prm_triplet(:,t)
            sl%tripletop(l)%permind=t
            sl%tripletop(l)%opind=o
            sl%tripletop(l)%detop=int(anint(lo_determ(uc%sym%op(o)%m)))
            allocate(sl%tripletop(l)%fmap(uc%na))
            sl%tripletop(l)%fmap=uc%sym%op(o)%fmap
        enddo
        enddo
        call mem%deallocate(so27,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block fixtriplet

    fixquartet: block
        real(r8), dimension(:,:,:), allocatable :: so81
        integer :: l,t,o,ntr

        if ( transposition ) then
            ntr=size(prm_quartet,2) ! all permutations
        else
            ntr=1 ! only identity
        endif

        ! And for quartets, eventually
        call mem%allocate(so81,[81,81,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        so81=0.0_r8
        do o=1,nop
            so81(:,:,o)=expandoperation_quartet(uc%sym%op(o)%m)
        enddo
        ! count and make space
        sl%n_quartet_operation=nop*ntr
        allocate(sl%quartetop(sl%n_quartet_operation))
        ! populate
        l=0
        do t=1,ntr
        do o=1,nop
            l=l+1
            sl%quartetop(l)%m3=uc%sym%op(o)%m
            sl%quartetop(l)%im3=uc%sym%op(o)%im
            call lo_gemm(so81(:,:,o),trm_quartet(:,:,t),sl%quartetop(l)%sotr)
            sl%quartetop(l)%perm=prm_quartet(:,t)
            sl%quartetop(l)%permind=t
            sl%quartetop(l)%opind=o
            sl%quartetop(l)%detop=int(anint(lo_determ(uc%sym%op(o)%m)))
            allocate(sl%quartetop(l)%fmap(uc%na))
            sl%quartetop(l)%fmap=uc%sym%op(o)%fmap
        enddo
        enddo
        call mem%deallocate(so81,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block fixquartet

    call mem%deallocate(trm_pair   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(trm_triplet,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(trm_quartet,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(prm_pair   ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(prm_triplet,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(prm_quartet,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> expand a 3x3 symmetry operation to a 9x9 operation, to operate on flattened tensors
module pure function expandoperation_pair(o) result(bigo)
    !> the original 3x3 operation
    real(r8), dimension(3,3), intent(in) :: o
    !> the resulting 9x9 operation
    real(r8), dimension(9,9) :: bigo
    !
    integer :: i,j,ii,jj,k,l
    !
    bigo=0.0_r8
    k=0
    do i=1,3
    do ii=1,3
        k=k+1
        l=0
        do j=1,3
        do jj=1,3
            l=l+1
            bigo(k,l)=o(ii,jj)*o(i,j)
        enddo
        enddo
    enddo
    enddo
    bigo=lo_chop(bigo,1E-11_r8)
end function

!> expand a 3x3 symmetry operation to a 27x27 operation, to operate on vector format of tensors
module pure function expandoperation_triplet(o) result(bigo)
    !> the original 3x3 operation
    real(r8), dimension(3,3), intent(in) :: o
    !> the resulting 27x27 operation
    real(r8), dimension(27,27) :: bigo
    !
    integer :: i,ii,iii,j,jj,jjj,l,m
    !
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
module pure function expandoperation_quartet(o) result(bigo)
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

end submodule
