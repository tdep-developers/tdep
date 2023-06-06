submodule (lo_symmetry_of_interactions) lo_symmetry_of_interactions_tuplets
use gottochblandat, only: lo_sqnorm,lo_progressbar_init,lo_progressbar
use lo_sorting, only: lo_qsort
implicit none
contains

!> build the coordination shells from tuplets
module subroutine tuplets_to_shells(sl,uc,ss,sh,mw,mem,verbosity)
    !> symmetry table
    class(lo_interaction_tensors), intent(inout) :: sl
    !> structures
    type(lo_crystalstructure), intent(in) :: uc,ss
    !> helper
    type(lo_symtabhelper), intent(inout) :: sh
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    !real(r8) :: t0,t1

    call mem%tick()

    ! This is a little tricky. Hmmm. Think the trick is to make two versions,
    ! one set uf tuplets that are transposeable, and one set that is not?
    ! Later problem I suppose. Start with the singlets, might as well always
    ! do those, cost nothing and we have the information alredy.
    singlet: block
        integer :: iatom,i

        ! This is easy enough, already taken care of in the spacegroup thing, no
        ! need to do anything fancy.
        allocate(sh%singlet(uc%sym%n_irreducible_atom))
        do i=1,uc%sym%n_irreducible_atom
            iatom=uc%sym%irr_to_all(i)
            sh%singlet(i)%i1=sl%uc_singlet(iatom)%i1
            sh%singlet(i)%n_unfold=uc%sym%irr_unfold_ctr(i)
            allocate(sh%singlet(i)%unfold_index(sh%singlet(i)%n_unfold))
            allocate(sh%singlet(i)%unfold_operation(sh%singlet(i)%n_unfold))
            sh%singlet(i)%unfold_index=uc%sym%irr_unfold_index(1:sh%singlet(i)%n_unfold,i)
            sh%singlet(i)%unfold_operation=uc%sym%irr_unfold_operation(1:sh%singlet(i)%n_unfold,i)
        enddo
        do iatom=1,uc%na
            sl%uc_singlet%unique=uc%sym%all_to_irr( sl%uc_singlet(iatom)%i1 )
            sl%uc_singlet%operation=uc%sym%op_all_from_irr( sl%uc_singlet(iatom)%i1 )
        enddo
        do iatom=1,ss%na
            sl%ss_singlet%unique=uc%sym%all_to_irr( sl%ss_singlet(iatom)%i1 )
            sl%uc_singlet%operation=uc%sym%op_all_from_irr( sl%ss_singlet(iatom)%i1 )
        enddo
    end block singlet

    ! Get unique pairs
    pair: block
        type(lo_tuplet_ucpair), dimension(:), allocatable :: pair
        integer, dimension(:), allocatable :: ind
        integer :: iatom,i,j,k,l

        ! count candidate pairs
        l=0
        do i=1,uc%sym%n_irreducible_atom
            iatom=uc%sym%irr_to_all(i)
            l=l+sh%uc_ctr_pair(iatom)
        enddo
        allocate(pair(l))
        ! store candidate pairs
        l=0
        do i=1,uc%sym%n_irreducible_atom
            iatom=uc%sym%irr_to_all(i)
            do j=1,sh%uc_ctr_pair(iatom)
                k=sh%uc_offset_pair(iatom)+j
                l=l+1
                pair(l)=sl%uc_pair(k)
            enddo
        enddo
        ! return the indices to the unique
        call return_unique_tuplets(sl,sh,pair,sl%pairop,ss,ind,mem,verbosity)
        ! store the unique
        allocate(sh%pair(size(ind)))
        do i=1,size(ind)
            j=ind(i)
            sh%pair(i)%i1             =pair(j)%i1
            sh%pair(i)%i2             =pair(j)%i2
            sh%pair(i)%ui1            =pair(j)%ui1
            sh%pair(i)%ui2            =pair(j)%ui2
            sh%pair(i)%v              =pair(j)%v
            sh%pair(i)%ind_ss_pair    =pair(j)%ind_ss_pair
            sh%pair(i)%n_unfold       =0
        enddo
        ! cleanup
        deallocate(ind)
        deallocate(pair)
        ! ensure the symmetry connections are constructed every which way
        call build_tuplet_connection(sh%pair,sl%uc_pair,sl%ss_pair,sl%pairop,sl,sh,ss,&
            sh%ss_ctr_pair,sh%uc_offset_pair,sh%ss_offset_pair,mw,mem,verbosity)
    end block pair

    ! Get unique triplets
    if ( sl%consider_triplet ) then
    triplet: block
        type(lo_tuplet_uctriplet), dimension(:), allocatable :: triplet
        integer, dimension(:), allocatable :: ind
        integer :: iatom,i,j,k,l

        ! count candidate triplets
        l=0
        do i=1,uc%sym%n_irreducible_atom
            iatom=uc%sym%irr_to_all(i)
            l=l+sh%uc_ctr_triplet(iatom)
        enddo
        allocate(triplet(l))
        ! store candidate triplets
        l=0
        do i=1,uc%sym%n_irreducible_atom
            iatom=uc%sym%irr_to_all(i)
            do j=1,sh%uc_ctr_triplet(iatom)
                k=sh%uc_offset_triplet(iatom)+j
                l=l+1
                triplet(l)=sl%uc_triplet(k)
            enddo
        enddo
        ! return the indices to the unique
        call return_unique_tuplets(sl,sh,triplet,sl%tripletop,ss,ind,mem,verbosity)
        ! store the unique
        allocate(sh%triplet(size(ind)))
        do i=1,size(ind)
            j=ind(i)
            sh%triplet(i)%i1             = triplet(j)%i1
            sh%triplet(i)%i2             = triplet(j)%i2
            sh%triplet(i)%i3             = triplet(j)%i3
            sh%triplet(i)%ui1            = triplet(j)%ui1
            sh%triplet(i)%ui2            = triplet(j)%ui2
            sh%triplet(i)%ui3            = triplet(j)%ui3
            sh%triplet(i)%v1             = triplet(j)%v1
            sh%triplet(i)%v2             = triplet(j)%v2
            sh%triplet(i)%v3             = triplet(j)%v3
            sh%triplet(i)%ind_ss_triplet = triplet(j)%ind_ss_triplet
            sh%triplet(i)%n_unfold       = 0
        enddo
        ! cleanup
        deallocate(ind)
        deallocate(triplet)
        ! ensure the symmetry connections are constructed every which way
        call build_tuplet_connection(sh%triplet,sl%uc_triplet,sl%ss_triplet,sl%tripletop,sl,sh,ss,&
            sh%ss_ctr_triplet,sh%uc_offset_triplet,sh%ss_offset_triplet,mw,mem,verbosity)
    end block triplet
    endif

    ! Get unique quartets
    if ( sl%consider_quartet ) then
    quartet: block
        type(lo_tuplet_ucquartet), dimension(:), allocatable :: quartet
        integer, dimension(:), allocatable :: ind
        integer :: iatom,i,j,k,l

        ! count candidate quartets
        l=0
        do i=1,uc%sym%n_irreducible_atom
            iatom=uc%sym%irr_to_all(i)
            l=l+sh%uc_ctr_quartet(iatom)
        enddo
        allocate(quartet(l))
        ! store candidate quartets
        l=0
        do i=1,uc%sym%n_irreducible_atom
            iatom=uc%sym%irr_to_all(i)
            do j=1,sh%uc_ctr_quartet(iatom)
                k=sh%uc_offset_quartet(iatom)+j
                l=l+1
                quartet(l)=sl%uc_quartet(k)
            enddo
        enddo
        ! return the indices to the unique
        call return_unique_tuplets(sl,sh,quartet,sl%quartetop,ss,ind,mem,verbosity)
        ! store the unique
        allocate(sh%quartet(size(ind)))
        do i=1,size(ind)
            j=ind(i)
            sh%quartet(i)%i1             = quartet(j)%i1
            sh%quartet(i)%i2             = quartet(j)%i2
            sh%quartet(i)%i3             = quartet(j)%i3
            sh%quartet(i)%i4             = quartet(j)%i4
            sh%quartet(i)%ui1            = quartet(j)%ui1
            sh%quartet(i)%ui2            = quartet(j)%ui2
            sh%quartet(i)%ui3            = quartet(j)%ui3
            sh%quartet(i)%ui4            = quartet(j)%ui4
            sh%quartet(i)%v1             = quartet(j)%v1
            sh%quartet(i)%v2             = quartet(j)%v2
            sh%quartet(i)%v3             = quartet(j)%v3
            sh%quartet(i)%v4             = quartet(j)%v4
            sh%quartet(i)%ind_ss_quartet = quartet(j)%ind_ss_quartet
            sh%quartet(i)%n_unfold       = 0
        enddo
        ! cleanup
        deallocate(ind)
        deallocate(quartet)
        ! ensure the symmetry connections are constructed every which way
        call build_tuplet_connection(sh%quartet,sl%uc_quartet,sl%ss_quartet,sl%quartetop,sl,sh,ss,&
            sh%ss_ctr_quartet,sh%uc_offset_quartet,sh%ss_offset_quartet,mw,mem,verbosity)
    end block quartet
    endif

    call mem%tock(__FILE__,__LINE__)
end subroutine

!> construct all the tuplets
module subroutine construct_tuplets(sl,uc,ss,sh,rc3,rc4,dc3,verbosity)
    !> symmetry table
    class(lo_interaction_tensors), intent(inout) :: sl
    !> structures
    type(lo_crystalstructure), intent(in) :: uc,ss
    !> helper
    type(lo_symtabhelper), intent(inout) :: sh
    !> cutoffs for all the things
    real(r8), intent(in) :: rc3,rc4,dc3
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8) :: rc

    ! build the singlets, easiest. Might as well always do that.
    singlet: block
        integer :: iatom
        allocate(sl%uc_singlet(uc%na))
        allocate(sl%ss_singlet(ss%na))
        do iatom=1,uc%na
            sl%uc_singlet(iatom)%i1=iatom
        enddo
        do iatom=1,ss%na
            sl%ss_singlet(iatom)%j1=iatom
            sl%ss_singlet(iatom)%i1=ss%info%index_in_unitcell( sl%ss_singlet(iatom)%j1 )
        enddo
        do iatom=1,uc%na
            call locate_unitcell_tuplet_in_supercell(sl%uc_singlet(iatom),sl,ss,sh,sl%uc_singlet(iatom)%ind_ss_singlet)
        enddo

        if ( verbosity .gt. 0 ) write(lo_iou,*) '... built singlets'
    end block singlet

    ! I have also concluded that I always have pairs
    pair: block
        integer :: iatom,ipair,l

        ! Count pairs
        sl%n_uc_pair=0
        sl%n_ss_pair=0
        do iatom=1,uc%na
            sl%n_uc_pair=sl%n_uc_pair+sh%dtuc%particle(iatom)%n
        enddo
        do iatom=1,ss%na
            sl%n_ss_pair=sl%n_ss_pair+sh%dtss%particle(iatom)%n
        enddo
        ! Space for pairs
        allocate(sl%uc_pair(sl%n_uc_pair))
        allocate(sl%ss_pair(sl%n_ss_pair))
        ! Store pairs
        l=0
        do iatom=1,uc%na
        do ipair=1,sh%dtuc%particle(iatom)%n
            l=l+1
            sl%uc_pair(l)%v=sh%dtuc%particle(iatom)%v(:,ipair)
            sl%uc_pair(l)%i1=iatom
            sl%uc_pair(l)%i2=sh%dtuc%particle(iatom)%ind(ipair)
            sl%uc_pair(l)%ui1=uc%sym%all_to_irr(sl%uc_pair(l)%i1)
            sl%uc_pair(l)%ui2=uc%sym%all_to_irr(sl%uc_pair(l)%i2)
        enddo
        enddo
        l=0
        do iatom=1,ss%na
        do ipair=1,sh%dtss%particle(iatom)%n
            l=l+1
            sl%ss_pair(l)%v=sh%dtss%particle(iatom)%v(:,ipair)
            sl%ss_pair(l)%j1=iatom
            sl%ss_pair(l)%j2=sh%dtss%particle(iatom)%ind(ipair)
            sl%ss_pair(l)%i1=ss%info%index_in_unitcell(sl%ss_pair(l)%j1)
            sl%ss_pair(l)%i2=ss%info%index_in_unitcell(sl%ss_pair(l)%j2)
            sl%ss_pair(l)%ui1=uc%sym%all_to_irr(sl%ss_pair(l)%i1)
            sl%ss_pair(l)%ui2=uc%sym%all_to_irr(sl%ss_pair(l)%i2)
        enddo
        enddo

        ! Count and get some offsets for later.
        allocate(sh%uc_ctr_pair(uc%na))
        allocate(sh%ss_ctr_pair(ss%na))
        allocate(sh%uc_offset_pair(uc%na))
        allocate(sh%ss_offset_pair(ss%na))
        sh%uc_ctr_pair=0
        sh%ss_ctr_pair=0
        sh%uc_offset_pair=0
        sh%ss_offset_pair=0
        do ipair=1,sl%n_uc_pair
            iatom=sl%uc_pair(ipair)%i1
            sh%uc_ctr_pair(iatom)=sh%uc_ctr_pair(iatom)+1
        enddo
        do ipair=1,sl%n_ss_pair
            iatom=sl%ss_pair(ipair)%j1
            sh%ss_ctr_pair(iatom)=sh%ss_ctr_pair(iatom)+1
        enddo
        l=0
        do iatom=1,uc%na
            sh%uc_offset_pair(iatom)=l
            l=l+sh%uc_ctr_pair(iatom)
        enddo
        l=0
        do iatom=1,ss%na
            sh%ss_offset_pair(iatom)=l
            l=l+sh%ss_ctr_pair(iatom)
        enddo

        if ( verbosity .gt. 0 ) write(lo_iou,*) '... built pairs'
    end block pair

    ! check wether triplets are relevant?
    rc=maxval([rc3,dc3])
    if ( rc .gt. 0.0_r8 ) then
    triplet: block
        real(r8), dimension(3) :: v0
        real(r8) :: rcsq
        integer :: iatom,itriplet,i,j,l

        rcsq=rc**2
        ! Count triplets
        sl%n_uc_triplet=0
        do iatom=1,uc%na
            do i=1,sh%dtuc%particle(iatom)%n
            if ( sh%dtuc%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtuc%particle(iatom)%n
                    if ( sh%dtuc%particle(iatom)%d(j) .gt. rc ) cycle
                    v0=sh%dtuc%particle(iatom)%v(:,j)-sh%dtuc%particle(iatom)%v(:,i)
                    if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                    sl%n_uc_triplet=sl%n_uc_triplet+1
                enddo
            enddo
        enddo
        sl%n_ss_triplet=0
        do iatom=1,ss%na
            do i=1,sh%dtss%particle(iatom)%n
            if ( sh%dtss%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtss%particle(iatom)%n
                    if ( sh%dtss%particle(iatom)%d(j) .gt. rc ) cycle
                    v0=sh%dtss%particle(iatom)%v(:,j)-sh%dtss%particle(iatom)%v(:,i)
                    if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                    sl%n_ss_triplet=sl%n_ss_triplet+1
                enddo
            enddo
        enddo

        ! Store triplets
        allocate(sl%uc_triplet(sl%n_uc_triplet))
        allocate(sl%ss_triplet(sl%n_ss_triplet))
        l=0
        do iatom=1,uc%na
            do i=1,sh%dtuc%particle(iatom)%n
            if ( sh%dtuc%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtuc%particle(iatom)%n
                    if ( sh%dtuc%particle(iatom)%d(j) .gt. rc ) cycle
                    v0=sh%dtuc%particle(iatom)%v(:,j)-sh%dtuc%particle(iatom)%v(:,i)
                    if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                    l=l+1
                    sl%uc_triplet(l)%i1=iatom
                    sl%uc_triplet(l)%i2=sh%dtuc%particle(iatom)%ind(i)
                    sl%uc_triplet(l)%i3=sh%dtuc%particle(iatom)%ind(j)
                    sl%uc_triplet(l)%v1=0.0_r8
                    sl%uc_triplet(l)%v2=sh%dtuc%particle(iatom)%v(:,i)
                    sl%uc_triplet(l)%v3=sh%dtuc%particle(iatom)%v(:,j)
                    sl%uc_triplet(l)%ui1=uc%sym%all_to_irr(sl%uc_triplet(l)%i1)
                    sl%uc_triplet(l)%ui2=uc%sym%all_to_irr(sl%uc_triplet(l)%i2)
                    sl%uc_triplet(l)%ui3=uc%sym%all_to_irr(sl%uc_triplet(l)%i3)
                enddo
            enddo
        enddo
        l=0
        do iatom=1,ss%na
            do i=1,sh%dtss%particle(iatom)%n
            if ( sh%dtss%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtss%particle(iatom)%n
                    if ( sh%dtss%particle(iatom)%d(j) .gt. rc ) cycle
                    v0=sh%dtss%particle(iatom)%v(:,j)-sh%dtss%particle(iatom)%v(:,i)
                    if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                    l=l+1
                    sl%ss_triplet(l)%j1=iatom
                    sl%ss_triplet(l)%j2=sh%dtss%particle(iatom)%ind(i)
                    sl%ss_triplet(l)%j3=sh%dtss%particle(iatom)%ind(j)
                    sl%ss_triplet(l)%i1=ss%info%index_in_unitcell(sl%ss_triplet(l)%j1)
                    sl%ss_triplet(l)%i2=ss%info%index_in_unitcell(sl%ss_triplet(l)%j2)
                    sl%ss_triplet(l)%i3=ss%info%index_in_unitcell(sl%ss_triplet(l)%j3)
                    sl%ss_triplet(l)%v1=0.0_r8
                    sl%ss_triplet(l)%v2=sh%dtss%particle(iatom)%v(:,i)
                    sl%ss_triplet(l)%v3=sh%dtss%particle(iatom)%v(:,j)
                    sl%ss_triplet(l)%ui1=uc%sym%all_to_irr(sl%ss_triplet(l)%i1)
                    sl%ss_triplet(l)%ui2=uc%sym%all_to_irr(sl%ss_triplet(l)%i2)
                    sl%ss_triplet(l)%ui3=uc%sym%all_to_irr(sl%ss_triplet(l)%i3)
                enddo
            enddo
        enddo

        ! Count and get some offsets for later.
        allocate(sh%uc_ctr_triplet(uc%na))
        allocate(sh%ss_ctr_triplet(ss%na))
        allocate(sh%uc_offset_triplet(uc%na))
        allocate(sh%ss_offset_triplet(ss%na))
        sh%uc_ctr_triplet=0
        sh%ss_ctr_triplet=0
        sh%uc_offset_triplet=0
        sh%ss_offset_triplet=0
        do itriplet=1,sl%n_uc_triplet
            iatom=sl%uc_triplet(itriplet)%i1
            sh%uc_ctr_triplet(iatom)=sh%uc_ctr_triplet(iatom)+1
        enddo
        do itriplet=1,sl%n_ss_triplet
            iatom=sl%ss_triplet(itriplet)%j1
            sh%ss_ctr_triplet(iatom)=sh%ss_ctr_triplet(iatom)+1
        enddo
        l=0
        do iatom=1,uc%na
            sh%uc_offset_triplet(iatom)=l
            l=l+sh%uc_ctr_triplet(iatom)
        enddo
        l=0
        do iatom=1,ss%na
            sh%ss_offset_triplet(iatom)=l
            l=l+sh%ss_ctr_triplet(iatom)
        enddo

        if ( verbosity .gt. 0 ) write(lo_iou,*) '... built triplets'
    end block triplet
    else
        ! no relevant triplets
        sl%n_uc_triplet=0
        sl%n_ss_triplet=0
    endif

    ! check wether quartets are relevant?
    rc=maxval([rc4])
    if ( rc .gt. 0.0_r8 ) then
    quartet: block
        real(r8), dimension(3) :: v0
        real(r8) :: rcsq
        integer :: iatom,i,j,k,l,iquartet

        rcsq=rc**2
        ! Count triplets
        sl%n_uc_quartet=0
        do iatom=1,uc%na
            do i=1,sh%dtuc%particle(iatom)%n
                if ( sh%dtuc%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtuc%particle(iatom)%n
                    if ( sh%dtuc%particle(iatom)%d(j) .gt. rc ) cycle
                    do k=1,sh%dtuc%particle(iatom)%n
                        if ( sh%dtuc%particle(iatom)%d(k) .gt. rc ) cycle
                        v0=sh%dtuc%particle(iatom)%v(:,j)-sh%dtuc%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtuc%particle(iatom)%v(:,k)-sh%dtuc%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtuc%particle(iatom)%v(:,k)-sh%dtuc%particle(iatom)%v(:,j)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        sl%n_uc_quartet=sl%n_uc_quartet+1
                    enddo
                enddo
            enddo
        enddo
        sl%n_ss_quartet=0
        do iatom=1,ss%na
            do i=1,sh%dtss%particle(iatom)%n
                if ( sh%dtss%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtss%particle(iatom)%n
                    if ( sh%dtss%particle(iatom)%d(j) .gt. rc ) cycle
                    do k=1,sh%dtss%particle(iatom)%n
                        if ( sh%dtss%particle(iatom)%d(k) .gt. rc ) cycle
                        v0=sh%dtss%particle(iatom)%v(:,j)-sh%dtss%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtss%particle(iatom)%v(:,k)-sh%dtss%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtss%particle(iatom)%v(:,k)-sh%dtss%particle(iatom)%v(:,j)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        sl%n_ss_quartet=sl%n_ss_quartet+1
                    enddo
                enddo
            enddo
        enddo

        ! Store quartets
        allocate(sl%uc_quartet(sl%n_uc_quartet))
        allocate(sl%ss_quartet(sl%n_ss_quartet))

        l=0
        do iatom=1,uc%na
            do i=1,sh%dtuc%particle(iatom)%n
                if ( sh%dtuc%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtuc%particle(iatom)%n
                    if ( sh%dtuc%particle(iatom)%d(j) .gt. rc ) cycle
                    do k=1,sh%dtuc%particle(iatom)%n
                        if ( sh%dtuc%particle(iatom)%d(k) .gt. rc ) cycle
                        v0=sh%dtuc%particle(iatom)%v(:,j)-sh%dtuc%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtuc%particle(iatom)%v(:,k)-sh%dtuc%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtuc%particle(iatom)%v(:,k)-sh%dtuc%particle(iatom)%v(:,j)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        l=l+1
                        sl%uc_quartet(l)%i1=iatom
                        sl%uc_quartet(l)%i2=sh%dtuc%particle(iatom)%ind(i)
                        sl%uc_quartet(l)%i3=sh%dtuc%particle(iatom)%ind(j)
                        sl%uc_quartet(l)%i4=sh%dtuc%particle(iatom)%ind(k)
                        sl%uc_quartet(l)%v1=0.0_r8
                        sl%uc_quartet(l)%v2=sh%dtuc%particle(iatom)%v(:,i)
                        sl%uc_quartet(l)%v3=sh%dtuc%particle(iatom)%v(:,j)
                        sl%uc_quartet(l)%v4=sh%dtuc%particle(iatom)%v(:,k)
                        sl%uc_quartet(l)%ui1=uc%sym%all_to_irr(sl%uc_quartet(l)%i1)
                        sl%uc_quartet(l)%ui2=uc%sym%all_to_irr(sl%uc_quartet(l)%i2)
                        sl%uc_quartet(l)%ui3=uc%sym%all_to_irr(sl%uc_quartet(l)%i3)
                        sl%uc_quartet(l)%ui4=uc%sym%all_to_irr(sl%uc_quartet(l)%i4)
                    enddo
                enddo
            enddo
        enddo
        l=0
        do iatom=1,ss%na
            do i=1,sh%dtss%particle(iatom)%n
                if ( sh%dtss%particle(iatom)%d(i) .gt. rc ) cycle
                do j=1,sh%dtss%particle(iatom)%n
                    if ( sh%dtss%particle(iatom)%d(j) .gt. rc ) cycle
                    do k=1,sh%dtss%particle(iatom)%n
                        if ( sh%dtss%particle(iatom)%d(k) .gt. rc ) cycle
                        v0=sh%dtss%particle(iatom)%v(:,j)-sh%dtss%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtss%particle(iatom)%v(:,k)-sh%dtss%particle(iatom)%v(:,i)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        v0=sh%dtss%particle(iatom)%v(:,k)-sh%dtss%particle(iatom)%v(:,j)
                        if ( lo_sqnorm(v0) .gt. rcsq ) cycle
                        l=l+1
                        sl%ss_quartet(l)%j1=iatom
                        sl%ss_quartet(l)%j2=sh%dtss%particle(iatom)%ind(i)
                        sl%ss_quartet(l)%j3=sh%dtss%particle(iatom)%ind(j)
                        sl%ss_quartet(l)%j4=sh%dtss%particle(iatom)%ind(k)
                        sl%ss_quartet(l)%i1=ss%info%index_in_unitcell(sl%ss_quartet(l)%j1)
                        sl%ss_quartet(l)%i2=ss%info%index_in_unitcell(sl%ss_quartet(l)%j2)
                        sl%ss_quartet(l)%i3=ss%info%index_in_unitcell(sl%ss_quartet(l)%j3)
                        sl%ss_quartet(l)%i4=ss%info%index_in_unitcell(sl%ss_quartet(l)%j4)
                        sl%ss_quartet(l)%v1=0.0_r8
                        sl%ss_quartet(l)%v2=sh%dtss%particle(iatom)%v(:,i)
                        sl%ss_quartet(l)%v3=sh%dtss%particle(iatom)%v(:,j)
                        sl%ss_quartet(l)%v4=sh%dtss%particle(iatom)%v(:,k)
                        sl%ss_quartet(l)%ui1=uc%sym%all_to_irr(sl%ss_quartet(l)%i1)
                        sl%ss_quartet(l)%ui2=uc%sym%all_to_irr(sl%ss_quartet(l)%i2)
                        sl%ss_quartet(l)%ui3=uc%sym%all_to_irr(sl%ss_quartet(l)%i3)
                        sl%ss_quartet(l)%ui4=uc%sym%all_to_irr(sl%ss_quartet(l)%i4)
                    enddo
                enddo
            enddo
        enddo

        ! Count and get some offsets for later.
        allocate(sh%uc_ctr_quartet(uc%na))
        allocate(sh%ss_ctr_quartet(ss%na))
        allocate(sh%uc_offset_quartet(uc%na))
        allocate(sh%ss_offset_quartet(ss%na))
        sh%uc_ctr_quartet=0
        sh%ss_ctr_quartet=0
        sh%uc_offset_quartet=0
        sh%ss_offset_quartet=0
        do iquartet=1,sl%n_uc_quartet
            iatom=sl%uc_quartet(iquartet)%i1
            sh%uc_ctr_quartet(iatom)=sh%uc_ctr_quartet(iatom)+1
        enddo
        do iquartet=1,sl%n_ss_quartet
            iatom=sl%ss_quartet(iquartet)%j1
            sh%ss_ctr_quartet(iatom)=sh%ss_ctr_quartet(iatom)+1
        enddo
        l=0
        do iatom=1,uc%na
            sh%uc_offset_quartet(iatom)=l
            l=l+sh%uc_ctr_quartet(iatom)
        enddo
        l=0
        do iatom=1,ss%na
            sh%ss_offset_quartet(iatom)=l
            l=l+sh%ss_ctr_quartet(iatom)
        enddo

        if ( verbosity .gt. 0 ) write(lo_iou,*) '... built quartets'
    end block quartet
    else
        ! no relevant triplets
        sl%n_uc_quartet=0
        sl%n_ss_quartet=0
    endif

    ! Identify unitcell tuplets in supercell
    locatetuplets: block
        integer :: i

        do i=1,sl%n_uc_pair
            call locate_unitcell_tuplet_in_supercell(sl%uc_pair(i),sl,ss,sh,sl%uc_pair(i)%ind_ss_pair)
        enddo
        do i=1,sl%n_uc_triplet
            call locate_unitcell_tuplet_in_supercell(sl%uc_triplet(i),sl,ss,sh,sl%uc_triplet(i)%ind_ss_triplet)
        enddo
        do i=1,sl%n_uc_quartet
            call locate_unitcell_tuplet_in_supercell(sl%uc_quartet(i),sl,ss,sh,sl%uc_quartet(i)%ind_ss_quartet)
        enddo
    end block locatetuplets
end subroutine

!> create hash from tuplet
function hash_tuplet(tuplet) result(hash)
    !> tuplet
    class(lo_tuplet), intent(in) :: tuplet
    !> hashed thing
    real(r8) :: hash

    real(r8), dimension(3) :: v0,v1,v2,v3,v4

    hash=0.0_r8
    select type(t=>tuplet)
    class is(lo_tuplet_pair)
        ! pairs, easy enough
        hash=hash+t%ui1+t%ui2
        hash=hash+lo_sqnorm(t%v)
    class is(lo_tuplet_triplet)
        ! triplets, also not tricky
        hash=hash+t%ui1+t%ui2+t%ui3
        v1=t%v2
        v2=t%v3
        v3=t%v3-t%v2
        hash=hash+lo_sqnorm(v1)+lo_sqnorm(v2)+lo_sqnorm(v3)
        hash=hash+dot_product(v1,v2)**2+dot_product(v1,v3)**2+dot_product(v2,v3)**2
    class is(lo_tuplet_quartet)
        ! quartets, a little annoying
        hash=hash+t%ui1+t%ui2+t%ui3+t%ui4
        ! center of mass
        v0=(t%v1+t%v2+t%v3+t%v4)*0.25_r8
        ! distances to tetrahedron corners from the center of mass
        v1=t%v1-v0
        v2=t%v2-v0
        v3=t%v3-v0
        v4=t%v4-v0
        hash=hash+lo_sqnorm(v1)+lo_sqnorm(v2)+lo_sqnorm(v3)+lo_sqnorm(v4)
    end select
end function

!> work out connections between tuplets
subroutine build_tuplet_connection(unique_tuplet,uc_tuplet,ss_tuplet,operations,sl,sh,ss,ss_ctr,uc_offset,ss_offset,mw,mem,verbosity)
    !> tuplets to reduce
    class(lo_tuplet), dimension(:), intent(inout) :: unique_tuplet
    class(lo_tuplet), dimension(:), intent(inout) :: uc_tuplet
    class(lo_tuplet), dimension(:), intent(inout) :: ss_tuplet
    !> set of operations to reduce with
    class(lo_tensor_tupletop), dimension(:), intent(in) :: operations
    !> symmetry table
    type(lo_interaction_tensors), intent(in) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> counters and offsets
    integer, dimension(:), intent(in) :: ss_ctr,uc_offset,ss_offset
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:), allocatable :: hash_unique,hash_uc,hash_ss
    real(r8) :: t0
    integer, dimension(:,:), allocatable :: di
    integer :: n_unique,n_uc,n_ss,n_op
    integer :: itup,iop,iun,iuc,iss,i,j,k
    integer :: a1,a2

    t0=walltime()
    if ( verbosity .gt. 0 ) call lo_progressbar_init()

    ! Get the sizes
    n_op=size(operations)
    n_unique=size(unique_tuplet)
    n_uc=size(uc_tuplet)
    n_ss=size(ss_tuplet)

    call mem%allocate(hash_unique,n_unique,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(hash_uc,n_uc,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(hash_ss,n_ss,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    hash_unique=0.0_r8
    hash_uc=0.0_r8
    hash_ss=0.0_r8

    do i=1,n_unique
        if ( mod(i,mw%n) .ne. mw%r ) cycle
        hash_unique(i)=hash_tuplet(unique_tuplet(i))
    enddo
    do i=1,n_uc
        if ( mod(i,mw%n) .ne. mw%r ) cycle
        hash_uc(i)=hash_tuplet(uc_tuplet(i))
    enddo
    do i=1,n_ss
        if ( mod(i,mw%n) .ne. mw%r ) cycle
        hash_ss(i)=hash_tuplet(ss_tuplet(i))
    enddo
    call mw%allreduce('sum',hash_unique)
    call mw%allreduce('sum',hash_uc)
    call mw%allreduce('sum',hash_ss)

    ! Then start mapping the unitcell
    call mem%allocate(di,[2,n_uc],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    di=0
    tuploop1: do iuc=1,n_uc
        ! make it parallel?
        if ( mod(iuc,mw%n) .ne. mw%r ) cycle
        ! we don't know what we are
        itup=0
        iop=0
        unloop1: do iun=1,n_unique
            ! hash has to match
            if ( abs(hash_uc(iuc)-hash_unique(iun)) .gt. lo_tol ) cycle unloop1
            ! check w.r.t all operations
            oploop1: do i=1,n_op
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,unique_tuplet(iun),uc_tuplet(iuc),operations(i),lo_tol) ) then
                    itup=iun
                    iop=i
                    exit unloop1
                endif
            enddo oploop1
        enddo unloop1
        if ( itup .gt. 0 ) then
            di(:,iuc)=[itup,iop]
        else
            call lo_stop_gracefully(['Could not locate tuplet'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        ! report
        if ( verbosity .gt. 0 ) then
            select type(unique_tuplet)
            class is(lo_tuplet_pair)
                call lo_progressbar(' ... mapping pairs',iuc,n_uc+n_ss,walltime()-t0)
            class is(lo_tuplet_triplet)
                call lo_progressbar(' ... mapping triplets',iuc,n_uc+n_ss,walltime()-t0)
            class is(lo_tuplet_quartet)
                call lo_progressbar(' ... mapping quartets',iuc,n_uc+n_ss,walltime()-t0)
            class default
                call lo_stop_gracefully(['Unknown tuplet type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            end select
        endif
    enddo tuploop1

    ! Sync and store
    call mw%allreduce('sum',di)
    do i=1,n_uc
        uc_tuplet(i)%unique=di(1,i)
        uc_tuplet(i)%operation=di(2,i)
    enddo
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! Now the supercell. Slightly differently, we can reuse some information here.
    call mem%allocate(di,[2,n_ss],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    di=0
    do a1=1,ss%na
    do i=1,ss_ctr(a1)
        iss=ss_offset(a1)+i
        if ( mod(iss,mw%n) .ne. mw%r ) cycle
        ! which atom in the unit cell is this
        a2=ss%info%index_in_unitcell(a1)
        iuc=uc_offset(a2)+i
        ! If the stars are aligned, unitcell pair k and supercell pair j are the
        ! same, so it should be real easy to match
        itup=0
        iop=0
        j=uc_tuplet(iuc)%unique
        k=uc_tuplet(iuc)%operation
        if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,unique_tuplet(j),ss_tuplet(iss),operations(k),lo_tol) ) then
            itup=j
            iop=k
        else
            call lo_stop_gracefully(['Thinking is hard!'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        if ( itup .gt. 0 ) then
            di(:,iss)=[itup,iop]
        else
            call lo_stop_gracefully(['Could not locate tuplet.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        ! report
        if ( verbosity .gt. 0 .and. iss .lt. n_ss ) then
            select type(unique_tuplet)
            class is(lo_tuplet_pair)
                call lo_progressbar(' ... mapping pairs',n_uc+iss,n_uc+n_ss,walltime()-t0)
            class is(lo_tuplet_triplet)
                call lo_progressbar(' ... mapping triplets',n_uc+iss,n_uc+n_ss,walltime()-t0)
            class is(lo_tuplet_quartet)
                call lo_progressbar(' ... mapping quartets',n_uc+iss,n_uc+n_ss,walltime()-t0)
            class default
                call lo_stop_gracefully(['Unknown tuplet type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
            end select
        endif
    enddo
    enddo

    ! sync and store
    call mw%allreduce('sum',di)
    do i=1,n_ss
        ss_tuplet(i)%unique=di(1,i)
        ss_tuplet(i)%operation=di(2,i)
    enddo
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! And cleanup
    call mem%deallocate(hash_unique,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(hash_uc,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(hash_ss,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! While I'm at it, create the unfolding pattern, always useful.
    ! At times I miss multiple inheritance, but what can you do.
    do i=1,n_unique
        select type(t=>unique_tuplet(i))
        type is(lo_tuplet_protpair)
            t%n_unfold=0
        type is(lo_tuplet_prottriplet)
            t%n_unfold=0
        type is(lo_tuplet_protquartet)
            t%n_unfold=0
        class default
            call lo_stop_gracefully(['Unknown tuplet type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        end select
    enddo

    ! First count
    do i=1,n_uc
        itup=uc_tuplet(i)%unique
        select type(t=>unique_tuplet(itup))
        type is(lo_tuplet_protpair)
            t%n_unfold=t%n_unfold+1
        type is(lo_tuplet_prottriplet)
            t%n_unfold=t%n_unfold+1
        type is(lo_tuplet_protquartet)
            t%n_unfold=t%n_unfold+1
        end select
    enddo
    ! Make space
    do i=1,n_unique
        select type(t=>unique_tuplet(i))
        type is(lo_tuplet_protpair)
            allocate(t%unfold_index(t%n_unfold))
            allocate(t%unfold_operation(t%n_unfold))
            t%unfold_index=0
            t%unfold_operation=0
            t%n_unfold=0
        type is(lo_tuplet_prottriplet)
            allocate(t%unfold_index(t%n_unfold))
            allocate(t%unfold_operation(t%n_unfold))
            t%unfold_index=0
            t%unfold_operation=0
            t%n_unfold=0
        type is(lo_tuplet_protquartet)
            allocate(t%unfold_index(t%n_unfold))
            allocate(t%unfold_operation(t%n_unfold))
            t%unfold_index=0
            t%unfold_operation=0
            t%n_unfold=0
        end select
    enddo
    ! Store
    do i=1,n_uc
        itup=uc_tuplet(i)%unique
        iop=uc_tuplet(i)%operation
        select type(t=>unique_tuplet(itup))
        type is(lo_tuplet_protpair)
            t%n_unfold=t%n_unfold+1
            t%unfold_index(t%n_unfold)=i
            t%unfold_operation(t%n_unfold)=iop
        type is(lo_tuplet_prottriplet)
            t%n_unfold=t%n_unfold+1
            t%unfold_index(t%n_unfold)=i
            t%unfold_operation(t%n_unfold)=iop
        type is(lo_tuplet_protquartet)
            t%n_unfold=t%n_unfold+1
            t%unfold_index(t%n_unfold)=i
            t%unfold_operation(t%n_unfold)=iop
        end select
    enddo

    ! final report?
    if ( verbosity .gt. 0 ) then
        select type(unique_tuplet)
        class is(lo_tuplet_pair)
            call lo_progressbar(' ... mapping pairs',n_uc+n_ss,n_uc+n_ss,walltime()-t0)
        class is(lo_tuplet_triplet)
            call lo_progressbar(' ... mapping triplets',n_uc+n_ss,n_uc+n_ss,walltime()-t0)
        class is(lo_tuplet_quartet)
            call lo_progressbar(' ... mapping quartets',n_uc+n_ss,n_uc+n_ss,walltime()-t0)
        class default
            call lo_stop_gracefully(['Unknown tuplet type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        end select
    endif
    ! Could be worth building the outwards spinny thing?
end subroutine

!> return the unique tuplets from a long list of tuplets
subroutine return_unique_tuplets(sl,sh,tuplets,operations,ss,index_to_unique,mem,verbosity)
    !> settings and stuff
    type(lo_interaction_tensors), intent(in) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> tuplets to reduce
    class(lo_tuplet), dimension(:), intent(in) :: tuplets
    !> set of operations to reduce with
    class(lo_tensor_tupletop), dimension(:), intent(in) :: operations
    !> structure
    type(lo_crystalstructure), intent(in) :: ss
    !> indices to unique tuplets
    integer, dimension(:), allocatable, intent(out) :: index_to_unique
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:), allocatable :: uhash,dr
    real(r8) :: t0,h1
    integer, dimension(:), allocatable :: di,dj
    integer :: i,t1,t2,op,unique_counter,ntup,nop
    logical :: newtuplet

    ! start the timer
    t0=walltime()
    ! number of tuplets
    ntup=size(tuplets)
    ! number of operations
    nop=size(operations)

    ! some temporary space
    call mem%allocate(di,ntup,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(uhash,ntup,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    di=0
    uhash=-lo_huge

    if ( verbosity .gt. 0 ) then
        call lo_progressbar_init()
    endif

    ! Start reducing
    unique_counter=0
    do t1=1,ntup
        ! I assume this is a brand new tuplet
        newtuplet=.true.
        ! hash it
        h1=hash_tuplet(tuplets(t1))
        ! now we take the next tuplet, and compare with those that are already there
        t2l: do t2=1,unique_counter
            ! hash has to match
            if ( abs(h1-uhash(t2)) .gt. lo_tol ) cycle
            ! Compare with all of them
            do op=1,nop
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,tuplets( t1 ),tuplets( di(t2) ),operations(op),lo_tol) ) then
                    newtuplet=.false.
                    exit t2l
                endif
            enddo
        enddo t2l

        ! Now, I might have a new unique tuplet!
        if ( newtuplet ) then
            ! update the counter
            unique_counter=unique_counter+1
            ! add the hash to the list of unique hashes
            uhash(unique_counter)=h1
            ! make a note of the unique tuplets
            di(unique_counter)=t1
        endif
        ! and tell the world how it's going
        if ( verbosity .gt. 0 ) then
            select type(tuplets)
            class is(lo_tuplet_pair)
                call lo_progressbar(' ... reducing pairs',t1,ntup,walltime()-t0)
            class is(lo_tuplet_triplet)
                call lo_progressbar(' ... reducing triplets',t1,ntup,walltime()-t0)
            class is(lo_tuplet_quartet)
                call lo_progressbar(' ... reducing quartets',t1,ntup,walltime()-t0)
            class default
                call lo_stop_gracefully(['Unknown tuplet type'],lo_exitcode_param,__FILE__,__LINE__)
            end select
        endif
    enddo

    ! It's always neat to return the tuplets sorted by size. Not that
    ! it matters much, but is not a bad idea.
    call mem%allocate(dj,unique_counter,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(dr,unique_counter,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dj=0
    dr=0.0_r8
    do i=1,unique_counter
        t1=di(i)
        select type(tuplets)
        class is(lo_tuplet_pair)
            dr(i)=norm2( tuplets(t1)%v )
        class is(lo_tuplet_triplet)
            dr(i)=norm2(tuplets(t1)%v2) + norm2(tuplets(t1)%v3) + norm2(tuplets(t1)%v2-tuplets(t1)%v3)
        class is(lo_tuplet_quartet)
            dr(i)=norm2(tuplets(t1)%v2) + norm2(tuplets(t1)%v3) + norm2(tuplets(t1)%v4) + &
            norm2(tuplets(t1)%v2-tuplets(t1)%v3) + &
            norm2(tuplets(t1)%v2-tuplets(t1)%v4) + &
            norm2(tuplets(t1)%v3-tuplets(t1)%v4)
        class default
            call lo_stop_gracefully(['Unknown tuplet type'],lo_exitcode_param,__FILE__,__LINE__)
        end select
    enddo

    ! Return them sorted
    call lo_qsort(dr,dj)
    allocate(index_to_unique(unique_counter))
    !index_to_unique=di(1:unique_counter)
    do i=1,unique_counter
       index_to_unique(i)=di(dj(i))
    enddo
    call mem%deallocate(dj,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(dr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(uhash,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> horribly convoluted way of comparing two tuplets
module function compare_tuplets_after_one_is_rotated(sl,sh,ss,t1,t2,op,tol) result(match)
    !> symmetry list thing
    type(lo_interaction_tensors), intent(in) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> two tuplets to compare
    class(lo_tuplet), intent(in) :: t1,t2
    !> operation to apply
    class(lo_tensor_tupletop), intent(in) :: op
    !> tolerance
    real(r8), intent(in) :: tol
    !> do they match?
    logical :: match

    integer :: ssi,a1,a2
    real(r8), dimension(3) :: v0
    real(r8) :: sqtol
    integer, dimension(2) :: pj1,pj2,pi1,pi2
    integer, dimension(3) :: tj1,ti1,ti2
    integer, dimension(4) :: qj1,qi1,qi2
    logical :: donesomething

    match=.false.
    donesomething=.false.
    sqtol=tol**2

    ! Stupid polymorphism, I hate it. Now we are at pairs, perhaps.
    select type(t1); class is(lo_tuplet_pair)
    select type(t2); class is(lo_tuplet_pair)
    select type(op); type is(lo_tensor_pairop)
        ! at least if I make it here, my overly ambitious polymorphic stuff is working.
        donesomething=.true.
        ! these are unitcell indices
        pi1(1)=t1%i1
        pi1(2)=t1%i2
        pi2(1)=t2%i1
        pi2(2)=t2%i2

        pj1(1)=t1%ui1
        pj1(2)=t1%ui2
        pj2(1)=t2%ui1
        pj2(2)=t2%ui2

        pi1=pi1(op%perm)
        pj1=pj1(op%perm)
        ! So, with this operation is basically
        ! op*t1%i1 = t2%i2, or op*pi1(1)=pi2(2)
        ! so that operation needs to be valid
        if ( op%fmap(pi1(1)) .ne. pi2(1) ) then
            match=.false.
            return
        endif
        ! species stuff should match
        if ( sum(abs(pj1-pj2)) .ne. 0 ) then
            match=.false.
            return
        endif
        ! apply the transposition + operation and compare
        v0=matmul(op%m3,t1%v)-t2%v
        if ( lo_sqnorm(v0) .lt. sqtol ) then
            match=.true.
            return
        else
            match=.false.
            return
        endif
        ! and now the pair is compared! great success
    end select
    end select
    end select

    ! now, perhaps, I want to compare two triplets instead!
    select type(t1); class is(lo_tuplet_uctriplet)
    select type(t2)
        class is(lo_tuplet_uctriplet)
        select type(op); type is(lo_tensor_tripletop)
            donesomething=.true.
            ! this version is comparing two unitcell triplets.
            ! fetch where in the supercell t1 is.
            ssi=t1%ind_ss_triplet
            ! get the supercell indices for the first one
            tj1=[sl%ss_triplet(ssi)%j1,sl%ss_triplet(ssi)%j2,sl%ss_triplet(ssi)%j3]
            ! permute
            tj1=tj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell( tj1(1) )
            a2=t2%i1
            ! is the operation valid between these?
            if ( op%fmap(a1) .ne. a2 ) then
                match=.false.
                return
            endif
            ! test if all the atoms are the same
            ti1=[t1%ui1,t1%ui2,t1%ui3]
            ti2=[t2%ui1,t2%ui2,t2%ui3]
            if ( sum(abs(ti1(op%perm)-ti2)) .ne. 0 ) then
                match=.false.
                return
            endif
            ! Fetch indices for the permuted
            call locate_unitcell_tuplet_via_indices(t1,tj1,sl,sh,ssi)
            ! now compare with rotations
            if ( lo_sqnorm(matmul(op%m3,sl%ss_triplet(ssi)%v2)-t2%v2) .lt. lo_sqtol ) then
                if ( lo_sqnorm(matmul(op%m3,sl%ss_triplet(ssi)%v3)-t2%v3) .lt. lo_sqtol ) then
                    match=.true.
                    return
                else
                    match=.false.
                    return
                endif
            else
                match=.false.
                return
            endif
        end select
        class is(lo_tuplet_sstriplet)
        select type(op); type is(lo_tensor_tripletop)
            donesomething=.true.
            ! now I am comparing the first triplet from the unitcell, and the second triplet is from the
            ! supercell. Fetch where in the supercell t1 is.
            ssi=t1%ind_ss_triplet
            ! get the supercell indices for the first one
            tj1=[sl%ss_triplet(ssi)%j1,sl%ss_triplet(ssi)%j2,sl%ss_triplet(ssi)%j3]
            ! permute
            tj1=tj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell(tj1(1))
            a2=ss%info%index_in_unitcell(t2%j1)
            ! are there valid operations betwee these?
            if ( op%fmap(a1) .ne. a2 ) then
                match=.false.
                return
            endif
            ! test if all the atoms are the same
            ti1=[t1%ui1,t1%ui2,t1%ui3]
            ti2=[t2%ui1,t2%ui2,t2%ui3]
            if ( sum(abs(ti1(op%perm)-ti2)) .ne. 0 ) then
                match=.false.
                return
            endif
            ! Fetch indices for the permuted
            call locate_unitcell_tuplet_via_indices(t1,tj1,sl,sh,ssi)
            ! now compare with rotations
            if ( lo_sqnorm(matmul(op%m3,sl%ss_triplet(ssi)%v2)-t2%v2) .lt. lo_sqtol ) then
                if ( lo_sqnorm(matmul(op%m3,sl%ss_triplet(ssi)%v3)-t2%v3) .lt. lo_sqtol ) then
                    match=.true.
                    return
                else
                    match=.false.
                    return
                endif
            else
                match=.false.
                return
            endif
        end select
    end select
    end select

    ! now, perhaps, I want to compare two quartets instead!
    select type(t1); class is(lo_tuplet_ucquartet)
    select type(t2)
        class is(lo_tuplet_ucquartet)
        select type(op); type is(lo_tensor_quartetop)
            donesomething=.true.
            ! this version is comparing two unitcell triplets.
            ! fetch where in the supercell t1 is.
            ssi=t1%ind_ss_quartet
            ! get the supercell indices for the first one
            qj1=[sl%ss_quartet(ssi)%j1,sl%ss_quartet(ssi)%j2,sl%ss_quartet(ssi)%j3,sl%ss_quartet(ssi)%j4]
            ! permute
            qj1=qj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell(qj1(1))
            a2=t2%i1
            ! are there valid operations between these?
            if ( op%fmap(a1) .ne. a2 ) then
                match=.false.
                return
            endif
            ! test if all the atoms are the same
            qi1=[t1%ui1,t1%ui2,t1%ui3,t1%ui4]
            qi2=[t2%ui1,t2%ui2,t2%ui3,t2%ui4]
            if ( sum(abs(qi1(op%perm)-qi2)) .ne. 0 ) then
                match=.false.
                return
            endif
            call locate_unitcell_tuplet_via_indices(t1,qj1,sl,sh,ssi)
            ! now compare with rotations
            match=.false.
            if ( lo_sqnorm(matmul(op%m3,sl%ss_quartet(ssi)%v2)-t2%v2) .lt. sqtol ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss_quartet(ssi)%v3)-t2%v3) .lt. sqtol ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss_quartet(ssi)%v4)-t2%v4) .lt. sqtol ) then
                match=.true.
            endif
            endif
            endif
        end select
        class is(lo_tuplet_ssquartet)
        select type(op); type is(lo_tensor_quartetop)
            donesomething=.true.
            ! now I am comparing the first quartet from the unitcell, and the second
            ! quartet is from the supercell.
            ! fetch where in the supercell t1 is.
            ssi=t1%ind_ss_quartet
            ! get the supercell indices for the first one
            qj1=[sl%ss_quartet(ssi)%j1,sl%ss_quartet(ssi)%j2,sl%ss_quartet(ssi)%j3,sl%ss_quartet(ssi)%j4]
            ! permute
            qj1=qj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell(qj1(1))
            a2=ss%info%index_in_unitcell(t2%j1)
            if ( op%fmap(a1) .ne. a2 ) then
                match=.false.
                return
            endif
            ! test if all the atoms are the same
            qi1=[t1%ui1,t1%ui2,t1%ui3,t1%ui4]
            qi2=[t2%ui1,t2%ui2,t2%ui3,t2%ui4]
            if ( sum(abs(qi1(op%perm)-qi2)) .ne. 0 ) then
                match=.false.
                return
            endif
            call locate_unitcell_tuplet_via_indices(t1,qj1,sl,sh,ssi)
            ! now compare with rotations
            match=.false.
            if ( lo_sqnorm(matmul(op%m3,sl%ss_quartet(ssi)%v2)-t2%v2) .lt. sqtol ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss_quartet(ssi)%v3)-t2%v3) .lt. sqtol ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss_quartet(ssi)%v4)-t2%v4) .lt. sqtol ) then
                match=.true.
            endif
            endif
            endif
        end select
    end select
    end select

    ! we should never make it here, in that case I have done something wrong.
    if ( donesomething .eqv. .false. ) then
        call lo_stop_gracefully(['Thinking is hard'],lo_exitcode_param,__FILE__,__LINE__)
    endif
end function

!> locate unitcell tuplet in supercell
subroutine locate_unitcell_tuplet_in_supercell(tuplet,sl,ss,sh,ind)
    !> tuplet
    class(lo_tuplet), intent(in) :: tuplet
    !> symmetry table
    class(lo_interaction_tensors), intent(in) :: sl
    !> structures
    type(lo_crystalstructure), intent(in) :: ss
    !> helper
    type(lo_symtabhelper), intent(inout) :: sh
    !> index to tuplet
    integer, intent(out) :: ind

    integer :: i,j,k

    select type(t=>tuplet)
    type is(lo_tuplet_ucsinglet)
        ind=0
        ssl1: do i=1,ss%na
            j=sl%ss_singlet(i)%i1
            if ( j .eq. t%i1 ) then
                ind=i
                exit ssl1
            endif
        enddo ssl1
        if ( ind .eq. 0 ) then
            call lo_stop_gracefully(['Could not locate singlet'],lo_exitcode_symmetry,__FILE__,__LINE__)
        else
            ! went fine!
            return
        endif
    type is(lo_tuplet_ucpair)
        ind=0
        ssl2: do i=1,ss%na
            if ( ss%info%index_in_unitcell(i) .ne. t%i1 ) cycle
            do j=1,sh%ss_ctr_pair(i)
                k=sh%ss_offset_pair(i)+j
                if ( sl%ss_pair(k)%i2 .ne. t%i2 ) cycle
                if ( lo_sqnorm(sl%ss_pair(k)%v-t%v) .gt. lo_sqtol ) cycle
                ind=k
                exit ssl2
            enddo
        enddo ssl2
        if ( ind .eq. 0 ) then
            call lo_stop_gracefully(['Could not locate pair'],lo_exitcode_symmetry,__FILE__,__LINE__)
        else
            ! went fine!
            return
        endif
    type is(lo_tuplet_uctriplet)
        ind=0
        ssl3: do i=1,ss%na
            if ( ss%info%index_in_unitcell(i) .ne. t%i1 ) cycle
            do j=1,sh%ss_ctr_triplet(i)
                k=sh%ss_offset_triplet(i)+j
                if ( sl%ss_triplet(k)%i2 .ne. t%i2 ) cycle
                if ( sl%ss_triplet(k)%i3 .ne. t%i3 ) cycle
                if ( lo_sqnorm(sl%ss_triplet(k)%v2-t%v2) .gt. lo_sqtol ) cycle
                if ( lo_sqnorm(sl%ss_triplet(k)%v3-t%v3) .gt. lo_sqtol ) cycle
                ind=k
                exit ssl3
            enddo
        enddo ssl3
        if ( ind .eq. 0 ) then
            call lo_stop_gracefully(['Could not locate triplet'],lo_exitcode_symmetry,__FILE__,__LINE__)
        else
            ! went fine!
            return
        endif
    type is(lo_tuplet_ucquartet)
        ind=0
        ssl4: do i=1,ss%na
            if ( ss%info%index_in_unitcell(i) .ne. t%i1 ) cycle
            do j=1,sh%ss_ctr_quartet(i)
                k=sh%ss_offset_quartet(i)+j
                if ( sl%ss_quartet(k)%i2 .ne. t%i2 ) cycle
                if ( sl%ss_quartet(k)%i3 .ne. t%i3 ) cycle
                if ( sl%ss_quartet(k)%i4 .ne. t%i4 ) cycle
                if ( lo_sqnorm(sl%ss_quartet(k)%v2-t%v2) .gt. lo_sqtol ) cycle
                if ( lo_sqnorm(sl%ss_quartet(k)%v3-t%v3) .gt. lo_sqtol ) cycle
                if ( lo_sqnorm(sl%ss_quartet(k)%v4-t%v4) .gt. lo_sqtol ) cycle
                ind=k
                exit ssl4
            enddo
        enddo ssl4
        if ( ind .eq. 0 ) then
            call lo_stop_gracefully(['Could not locate quartet'],lo_exitcode_symmetry,__FILE__,__LINE__)
        else
            ! went fine!
            return
        endif
    end select
end subroutine

!> locate unitcell tuplet in supercell
subroutine locate_unitcell_tuplet_via_indices(tuplet,atomind,sl,sh,ind)
    !> tuplet (to determine what kind it is)
    class(lo_tuplet), intent(in) :: tuplet
    !> index to atoms
    integer, dimension(:), intent(in) :: atomind
    !> symmetry table
    class(lo_interaction_tensors), intent(in) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> index to tuplet
    integer, intent(out) :: ind

    integer :: j1,j2,j3,j4,i,j

    select type(t=>tuplet)
    class is(lo_tuplet_triplet)
        ind=0
        j1=atomind(1)
        j2=atomind(2)
        j3=atomind(3)
        trloop: do i=1,sh%ss_ctr_triplet(j1)
            j=sh%ss_offset_triplet(j1)+i
            if ( sl%ss_triplet(j)%j2 .ne. j2 ) cycle
            if ( sl%ss_triplet(j)%j3 .ne. j3 ) cycle
            ! got it, I think!
            ind=j
            exit trloop
        enddo trloop
        if ( ind .eq. 0 ) then
            call lo_stop_gracefully(['Could not locate triplet'],lo_exitcode_symmetry,__FILE__,__LINE__)
        else
            ! went fine!
            return
        endif
    class is(lo_tuplet_quartet)
        ind=0
        j1=atomind(1)
        j2=atomind(2)
        j3=atomind(3)
        j4=atomind(4)
        qtloop: do i=1,sh%ss_ctr_quartet(j1)
            j=sh%ss_offset_quartet(j1)+i
            if ( sl%ss_quartet(j)%j2 .ne. j2 ) cycle
            if ( sl%ss_quartet(j)%j3 .ne. j3 ) cycle
            if ( sl%ss_quartet(j)%j4 .ne. j4 ) cycle
            ! got it, I think!
            ind=j
            exit qtloop
        enddo qtloop
        if ( ind .eq. 0 ) then
            call lo_stop_gracefully(['Could not locate triplet'],lo_exitcode_symmetry,__FILE__,__LINE__)
        else
            ! went fine!
            return
        endif
    class default
        call lo_stop_gracefully(['Unknown tuplet type'],lo_exitcode_param,__FILE__,__LINE__)
    end select
end subroutine

end submodule
