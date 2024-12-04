submodule (lo_symmetry_of_interactions) lo_symmetry_of_interactions_nullspace
!! gets the nullspace for all types of tuplets
use gottochblandat, only: lo_real_nullspace_coefficient_matrix,lo_sqnorm,lo_identitymatrix,lo_transpositionmatrix,lo_unflatten_3tensor
use type_blas_lapack_wrappers, only: lo_gemm
implicit none
contains

!> nullspace for on-site forces
module subroutine nullspace_fc_singlet(sl,sh,uc,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:), allocatable :: iM,rotM,invarM,coeffA
    real(r8), dimension(:,:), allocatable :: m0,m1,m2
    real(r8) :: f0
    integer :: iop,ne,nx,a1,a2,i,j,l

    ! First get the nullspace for the full thing:
    ne=uc%na*3
    call mem%allocate(iM    ,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(rotM  ,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(invarM,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    iM=0.0_r8
    rotM=0.0_r8
    invarM=0.0_r8
    call lo_identitymatrix(iM)
    ! Get the rotational invariance things
    do iop=1,sl%n_singlet_operation
        rotM=0.0_r8
        do a1=1,uc%na
            a2=sl%singletop(iop)%fmap(a1)
            rotM( (a2-1)*3+1:a2*3, (a1-1)*3+1:a1*3 )=sl%singletop(iop)%m3
        enddo
        invarM=invarM+rotM-IM
    enddo

    ! Then make sure it adds up to zero
    rotM=0.0_r8
    do a1=1,uc%na
        do i=1,3
            rotM( (a1-1)*3+i , i )=1.0_r8
        enddo
    enddo
    invarM=invarM+rotM
    ! Get the nullspace
    call lo_real_nullspace_coefficient_matrix(&
        invarM=invarM,&
        coeff=coeffA,&
        nvar=nx,&
        tolerance=lo_sqtol)

    ! Intermediate cleanup
    call mem%deallocate(iM    ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(rotM  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(invarM,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! Could be easy:
    if ( nx .eq. 0 ) then
        sl%nx_fc_singlet=0
        sl%n_fc_singlet_shell=0
        sl%have_fc_singlet=.false.
        return
    endif

    ! Apparently not, now we have to store things
    sl%nx_fc_singlet=nx
    sl%n_fc_singlet_shell=size(sh%singlet)
    sl%have_fc_singlet=.true.
    allocate(sl%fc_singlet_shell(sl%n_fc_singlet_shell))
    ! Do some sanity checks while sorting this out, never a bad idea
    call mem%allocate(m0,[3,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(m1,[3,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(m2,[3,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    shloop: do i=1,sl%n_fc_singlet_shell
        a1=sh%singlet(i)%i1             ! index in the unit cell
        m0=coeffA( (a1-1)*3+1:a1*3,: )  ! partial coefficient matrix
        ! Test that my brain is not turning to mush of old age:
        do j=1,sh%singlet(i)%n_unfold
            iop=sh%singlet(i)%unfold_operation(j)
            a2=sh%singlet(i)%unfold_index(j)
            m1=matmul(sl%singletop(iop)%m3,m0)
            m2=coeffA( (a2-1)*3+1:a2*3,: )
            m2=m2-m1
            m2=abs(m2)
            f0=sum(m2)
            if ( f0 .gt. lo_sqtol ) then
               call lo_stop_gracefully(['Bad operation singlets'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo
        ! Seems fine, count relevant for this shell:
        l=0
        do j=1,nx
            f0=sum(abs(m0(:,j)))
            if ( f0 .gt. lo_sqtol ) then
                l=l+1
            endif
        enddo
        ! In case of nothing, cheap exit
        if ( l .eq. 0 ) then
            sl%fc_singlet_shell(i)%nx=0
            cycle shloop
        endif
        ! In case of something, store:
        sl%fc_singlet_shell(i)%nx=l
        allocate(sl%fc_singlet_shell(i)%ind_local(l))
        allocate(sl%fc_singlet_shell(i)%ind_global(l))
        allocate(sl%fc_singlet_shell(i)%coeff(3,l))
        sl%fc_singlet_shell(i)%ind_local=0
        sl%fc_singlet_shell(i)%ind_global=0
        sl%fc_singlet_shell(i)%coeff=0.0_r8
        l=0
        do j=1,nx
            f0=sum(abs(m0(:,j)))
            if ( f0 .gt. lo_sqtol ) then
                l=l+1
                sl%fc_singlet_shell(i)%ind_global(l)=j
                sl%fc_singlet_shell(i)%ind_local(l)=j
                sl%fc_singlet_shell(i)%coeff(1:3,l)=m0(1:3,j)
            endif
        enddo
    enddo shloop

    ! And store the unfolding
    do i=1,sl%n_fc_singlet_shell
        sl%fc_singlet_shell(i)%n_unfold = sh%singlet(i)%n_unfold
        allocate( sl%fc_singlet_shell(i)%unfold_index( sl%fc_singlet_shell(i)%n_unfold ) )
        allocate( sl%fc_singlet_shell(i)%unfold_operation( sl%fc_singlet_shell(i)%n_unfold ) )
        sl%fc_singlet_shell(i)%unfold_index=sh%singlet(i)%unfold_index
        sl%fc_singlet_shell(i)%unfold_operation=sh%singlet(i)%unfold_operation
    enddo

    ! And a little cleanup
    call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(m1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(m2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> grab the nullspace for forceconstant pairs
module subroutine nullspace_fc_pair(sl,sh,ss,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> supercell (for some odd reason?)
    type(lo_crystalstructure), intent(in) :: ss
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:,:), allocatable :: op
    real(r8), dimension(:,:), allocatable :: m0
    integer, dimension(:), allocatable :: di
    integer :: shell,iop,ctr,nop
    integer :: global_ctr,i

    ! Space for the shells
    sl%n_fc_pair_shell=size(sh%pair)
    allocate(sl%fc_pair_shell(sl%n_fc_pair_shell))

    !@TODO populate shells with things? Will decide later what is needed.

    ! Some temporary space
    call mem%allocate(di,sl%n_pair_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    di=0
    global_ctr=0
    shloop: do shell=1,sl%n_fc_pair_shell
        ! If it's a self-term, we deal with it separately
        if ( lo_sqnorm(sh%pair(shell)%v) .lt. lo_sqtol ) then
            sl%fc_pair_shell(shell)%nx=0
            cycle shloop
        endif
        ! Figure out which the invariant operations are? Count first, I suppose
        di=0
        ctr=0
        do iop=1,sl%n_pair_operation
            if ( mod(iop,mw%n) .ne. mw%r ) cycle
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sh%pair(shell),sh%pair(shell),sl%pairop(iop),lo_tol) ) then
                ctr=ctr+1
                di(iop)=1
            endif
        enddo
        call mw%allreduce('max',di)
        call mw%size_and_offset(ctr,nop)

        if ( nop .gt. 0 ) then
            call mem%allocate(op,[9,9,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            op=0.0_r8
        else
            ! should be impossible to have zero operations
            call lo_stop_gracefully(['No invariant operations for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        ! Then store the operations and sync across ranks
        ctr=0
        do iop=1,sl%n_pair_operation
            if ( mod(iop,mw%n) .ne. mw%r ) cycle
            if ( di(iop) .eq. 0 ) cycle
            ctr=ctr+1
            op(:,:,ctr+mw%ctr_offset_per_rank(mw%r+1))=sl%pairop(iop)%sotr
        enddo
        call mw%allreduce('sum',op)
        ! Get the nullspace
        call lo_real_nullspace_coefficient_matrix(&
            invariant_operations=op,&
            coeff=sl%fc_pair_shell(shell)%coeff,&
            nvar=sl%fc_pair_shell(shell)%nx,&
            varind=sl%fc_pair_shell(shell)%ind_local,&
            tolerance=lo_sqtol)

        ! Store the global indices
        if ( sl%fc_pair_shell(shell)%nx .gt. 0 ) then
            allocate(sl%fc_pair_shell(shell)%ind_global( sl%fc_pair_shell(shell)%nx ))
            sl%fc_pair_shell(shell)%ind_global=0
            do i=1,sl%fc_pair_shell(shell)%nx
                global_ctr=global_ctr+1
                sl%fc_pair_shell(shell)%ind_global(i)=global_ctr
            enddo
        else
            ! Should always be some in this case
            call lo_stop_gracefully(['No independent for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif

        ! This is a sensible place to figure out wether I am worthy or not.
        call mem%allocate(m0,[9,sl%fc_pair_shell(shell)%nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        m0=0.0_r8
        do iop=1,size(op,3)
            call lo_gemm(op(:,:,iop),sl%fc_pair_shell(shell)%coeff,m0)
            if ( norm2(m0-sl%fc_pair_shell(shell)%coeff) .gt. lo_sqtol*size(m0) ) then
                call lo_stop_gracefully(['Not as invariant as it should be.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo

        ! temporary cleanup
        call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    enddo shloop
    ! And make a note of the total number of independent
    sl%nx_fc_pair=global_ctr

    ! This is also a sensible place to generate the unfolding thing, i.e. how to construct
    ! all the pairs from the irreducible and not the other way around.
    do shell=1,sl%n_fc_pair_shell
        sl%fc_pair_shell(shell)%n_unfold = sh%pair(shell)%n_unfold
        allocate(sl%fc_pair_shell(shell)%unfold_index( sh%pair(shell)%n_unfold ))
        allocate(sl%fc_pair_shell(shell)%unfold_operation( sh%pair(shell)%n_unfold ))
        sl%fc_pair_shell(shell)%unfold_index    =sh%pair(shell)%unfold_index
        sl%fc_pair_shell(shell)%unfold_operation=sh%pair(shell)%unfold_operation
    enddo

    ! Cleanup
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> grab the nullspace for forceconstant triplets
module subroutine nullspace_fc_triplet(sl,sh,ss,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> supercell (for some odd reason?)
    type(lo_crystalstructure), intent(in) :: ss
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:,:), allocatable :: op
    real(r8), dimension(:,:), allocatable :: m0
    integer, dimension(:), allocatable :: di
    integer :: shell,iop,ctr,nop
    integer :: global_ctr,i

    ! Space for the shells
    sl%n_fc_triplet_shell=size(sh%triplet)
    allocate(sl%fc_triplet_shell(sl%n_fc_triplet_shell))

    !@TODO populate shells with things? Will decide later what is needed.

    ! Some temporary space
    call mem%allocate(di,sl%n_triplet_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    di=0
    global_ctr=0
    shloop: do shell=1,sl%n_fc_triplet_shell
        ! If it's a self-term, we deal with it separately? Yes no maybe.
        if ( lo_sqnorm(sh%triplet(shell)%v2) .lt. lo_sqtol ) then
        if ( lo_sqnorm(sh%triplet(shell)%v3) .lt. lo_sqtol ) then
            sl%fc_triplet_shell(shell)%nx=0
            cycle shloop
        endif
        endif
        ! Figure out which the invariant operations are? Count first, I suppose
        di=0
        ctr=0
        do iop=1,sl%n_triplet_operation
            if ( mod(iop,mw%n) .ne. mw%r ) cycle
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sh%triplet(shell),sh%triplet(shell),sl%tripletop(iop),lo_tol) ) then
                ctr=ctr+1
                di(iop)=1
            endif
        enddo
        call mw%allreduce('max',di)
        call mw%size_and_offset(ctr,nop)

        if ( nop .gt. 0 ) then
            call mem%allocate(op,[27,27,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            op=0.0_r8
        else
            ! should be impossible to have zero operations
            call lo_stop_gracefully(['No invariant operations for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        ! Then store the operations and sync across ranks
        ctr=0
        do iop=1,sl%n_triplet_operation
            if ( mod(iop,mw%n) .ne. mw%r ) cycle
            if ( di(iop) .eq. 0 ) cycle
            ctr=ctr+1
            op(:,:,ctr+mw%ctr_offset_per_rank(mw%r+1))=sl%tripletop(iop)%sotr
        enddo
        call mw%allreduce('sum',op)
        ! Get the nullspace
        call lo_real_nullspace_coefficient_matrix(&
            invariant_operations=op,&
            coeff=sl%fc_triplet_shell(shell)%coeff,&
            nvar=sl%fc_triplet_shell(shell)%nx,&
            varind=sl%fc_triplet_shell(shell)%ind_local,&
            tolerance=lo_sqtol)

        ! Store the global indices
        if ( sl%fc_triplet_shell(shell)%nx .gt. 0 ) then
            allocate(sl%fc_triplet_shell(shell)%ind_global( sl%fc_triplet_shell(shell)%nx ))
            sl%fc_triplet_shell(shell)%ind_global=0
            do i=1,sl%fc_triplet_shell(shell)%nx
                global_ctr=global_ctr+1
                sl%fc_triplet_shell(shell)%ind_global(i)=global_ctr
            enddo
        else
            ! Should always be some in this case
            call lo_stop_gracefully(['No independent for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif

        ! This is a sensible place to figure out wether I am worthy or not.
        call mem%allocate(m0,[27,sl%fc_triplet_shell(shell)%nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        m0=0.0_r8
        do iop=1,size(op,3)
            call lo_gemm(op(:,:,iop),sl%fc_triplet_shell(shell)%coeff,m0)
            if ( norm2(m0-sl%fc_triplet_shell(shell)%coeff) .gt. lo_sqtol*size(m0) ) then
                call lo_stop_gracefully(['Not as invariant as it should be.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo

        ! temporary cleanup
        call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    enddo shloop
    ! And make a note of the total number of independent
    sl%nx_fc_triplet=global_ctr

    ! This is also a sensible place to generate the unfolding thing, i.e. how to construct
    ! all the triplets from the irreducible and not the other way around.
    do shell=1,sl%n_fc_triplet_shell
        sl%fc_triplet_shell(shell)%n_unfold = sh%triplet(shell)%n_unfold
        allocate(sl%fc_triplet_shell(shell)%unfold_index( sh%triplet(shell)%n_unfold ))
        allocate(sl%fc_triplet_shell(shell)%unfold_operation( sh%triplet(shell)%n_unfold ))
        sl%fc_triplet_shell(shell)%unfold_index    =sh%triplet(shell)%unfold_index
        sl%fc_triplet_shell(shell)%unfold_operation=sh%triplet(shell)%unfold_operation
    enddo

    ! Cleanup
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> grab the nullspace for forceconstant quartets
module subroutine nullspace_fc_quartet(sl,sh,ss,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> supercell (for some odd reason?)
    type(lo_crystalstructure), intent(in) :: ss
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:,:), allocatable :: op
    real(r8), dimension(:,:), allocatable :: m0
    integer, dimension(:), allocatable :: di
    integer :: shell,iop,ctr,nop
    integer :: global_ctr,i

    ! Space for the shells
    sl%n_fc_quartet_shell=size(sh%quartet)
    allocate(sl%fc_quartet_shell(sl%n_fc_quartet_shell))

    !@TODO populate shells with things? Will decide later what is needed.

    ! Some temporary space
    call mem%allocate(di,sl%n_quartet_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    di=0
    global_ctr=0
    shloop: do shell=1,sl%n_fc_quartet_shell
        ! If it's a self-term, we deal with it separately? Yes no maybe.
        if ( lo_sqnorm(sh%quartet(shell)%v2) .lt. lo_sqtol ) then
        if ( lo_sqnorm(sh%quartet(shell)%v3) .lt. lo_sqtol ) then
        if ( lo_sqnorm(sh%quartet(shell)%v4) .lt. lo_sqtol ) then
            sl%fc_quartet_shell(shell)%nx=0
            cycle shloop
        endif
        endif
        endif
        ! Figure out which the invariant operations are? Count first, I suppose
        di=0
        ctr=0
        do iop=1,sl%n_quartet_operation
            if ( mod(iop,mw%n) .ne. mw%r ) cycle
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sh%quartet(shell),sh%quartet(shell),sl%quartetop(iop),lo_tol) ) then
                ctr=ctr+1
                di(iop)=1
            endif
        enddo
        call mw%allreduce('max',di)
        call mw%size_and_offset(ctr,nop)

        if ( nop .gt. 0 ) then
            call mem%allocate(op,[81,81,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            op=0.0_r8
        else
            ! should be impossible to have zero operations
            call lo_stop_gracefully(['No invariant operations for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        ! Then store the operations and sync across ranks
        ctr=0
        do iop=1,sl%n_quartet_operation
            if ( mod(iop,mw%n) .ne. mw%r ) cycle
            if ( di(iop) .eq. 0 ) cycle
            ctr=ctr+1
            op(:,:,ctr+mw%ctr_offset_per_rank(mw%r+1))=sl%quartetop(iop)%sotr
        enddo
        call mw%allreduce('sum',op)
        ! Get the nullspace
        call lo_real_nullspace_coefficient_matrix(&
            invariant_operations=op,&
            coeff=sl%fc_quartet_shell(shell)%coeff,&
            nvar=sl%fc_quartet_shell(shell)%nx,&
            varind=sl%fc_quartet_shell(shell)%ind_local,&
            tolerance=lo_sqtol)

        ! Store the global indices
        if ( sl%fc_quartet_shell(shell)%nx .gt. 0 ) then
            allocate(sl%fc_quartet_shell(shell)%ind_global( sl%fc_quartet_shell(shell)%nx ))
            sl%fc_quartet_shell(shell)%ind_global=0
            do i=1,sl%fc_quartet_shell(shell)%nx
                global_ctr=global_ctr+1
                sl%fc_quartet_shell(shell)%ind_global(i)=global_ctr
            enddo
        else
            ! Should always be some in this case
            call lo_stop_gracefully(['No independent for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif

        ! This is a sensible place to figure out wether I am worthy or not.
        call mem%allocate(m0,[81,sl%fc_quartet_shell(shell)%nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        m0=0.0_r8
        do iop=1,size(op,3)
            call lo_gemm(op(:,:,iop),sl%fc_quartet_shell(shell)%coeff,m0)
            if ( norm2(m0-sl%fc_quartet_shell(shell)%coeff) .gt. lo_tol*size(m0) ) then
                call lo_stop_gracefully(['Not as invariant as it should be.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo

        ! temporary cleanup
        call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    enddo shloop
    ! And make a note of the total number of independent
    sl%nx_fc_quartet=global_ctr

    ! This is also a sensible place to generate the unfolding thing, i.e. how to construct
    ! all the quartets from the irreducible and not the other way around.
    do shell=1,sl%n_fc_quartet_shell
        sl%fc_quartet_shell(shell)%n_unfold = sh%quartet(shell)%n_unfold
        allocate(sl%fc_quartet_shell(shell)%unfold_index( sh%quartet(shell)%n_unfold ))
        allocate(sl%fc_quartet_shell(shell)%unfold_operation( sh%quartet(shell)%n_unfold ))
        sl%fc_quartet_shell(shell)%unfold_index    =sh%quartet(shell)%unfold_index
        sl%fc_quartet_shell(shell)%unfold_operation=sh%quartet(shell)%unfold_operation
    enddo

    ! Cleanup
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> nullspace for dielectric tensor
module subroutine nullspace_eps_global(sl,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:), allocatable :: iM,rotM,invarM
    integer :: iop

    ! First get the nullspace for the full thing:
    call mem%allocate(iM    ,[9,9],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(rotM  ,[9,9],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(invarM,[9,9],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    iM=0.0_r8
    rotM=0.0_r8
    invarM=0.0_r8
    call lo_identitymatrix(iM)
    ! Get the rotational invariance things
    do iop=1,sl%n_singlet_operation
        rotM=expandoperation_pair(sl%singletop(iop)%m3)
        invarM=invarM+rotM-IM
    enddo
    ! And the transposition symmetry
    call lo_transpositionmatrix(rotM)
    invarM=invarM+rotM-IM
    ! Get the nullspace
    call lo_real_nullspace_coefficient_matrix(&
        invarM=invarM,&
        coeff=sl%eps_global_shell%coeff,&
        nvar=sl%eps_global_shell%nx,&
        tolerance=lo_sqtol)
    if ( sl%eps_global_shell%nx .lt. 1 .or. sl%eps_global_shell%nx .gt. 6 ) then
        call lo_stop_gracefully(['Unphysical number of components in the dielectric tensor'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
    else
        sl%have_eps_global=.true.
        sl%nx_eps_global=sl%eps_global_shell%nx
    endif
    ! Intermediate cleanup
    call mem%deallocate(iM    ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(rotM  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(invarM,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> nullspace for dielectric on-site response, i.e. Raman tensor
module subroutine nullspace_eps_singlet(sl,sh,uc,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !real(r8), dimension(27) :: v27
    real(r8), dimension(27,27) :: m27
    real(r8), dimension(:,:), allocatable :: iM,rotM,invarM,coeffA
    real(r8), dimension(:,:), allocatable :: m0,m1,m2
    real(r8) :: f0
    integer :: iop,ne,nx,a1,a2,i,j,l

    ! First get the nullspace for the full thing:
    ne=uc%na*27
    call mem%allocate(iM    ,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(rotM  ,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(invarM,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    iM=0.0_r8
    rotM=0.0_r8
    invarM=0.0_r8
    call lo_identitymatrix(iM)

    ! Spacegroup things should hold
    do iop=1,sl%n_triplet_operation
        ! Can not permute the position index (index 3) to any of the electric field indices.
        if ( sl%tripletop(iop)%perm(3) .ne. 3 ) cycle
        rotM=0.0_r8
        do a1=1,uc%na
            a2=sl%tripletop(iop)%fmap(a1)
            rotM( (a2-1)*27+1:a2*27, (a1-1)*27+1:a1*27 )=sl%tripletop(iop)%sotr
        enddo
        invarM=invarM+rotM-IM
    enddo

    ! Then make sure it adds up to zero
    rotM=0.0_r8
    do a1=1,uc%na
        do i=1,27
            rotM( (a1-1)*27+i , i )=1.0_r8
        enddo
    enddo
    invarM=invarM+rotM
    ! Get the nullspace
    call lo_real_nullspace_coefficient_matrix(&
        invarM=invarM,&
        coeff=coeffA,&
        nvar=nx,&
        tolerance=lo_sqtol)
    ! Intermediate cleanup
    call mem%deallocate(iM    ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(rotM  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(invarM,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! Could be easy:
    if ( nx .eq. 0 ) then
        sl%nx_eps_singlet=0
        sl%n_eps_singlet_shell=0
        sl%have_eps_singlet=.false.
        return
    endif

    ! Apparently not, now we have to store things
    sl%nx_eps_singlet=nx
    sl%n_eps_singlet_shell=size(sh%singlet)
    sl%have_eps_singlet=.true.
    allocate(sl%eps_singlet_shell(sl%n_eps_singlet_shell))
    ! Do some sanity checks while sorting this out, never a bad idea
    call mem%allocate(m0,[27,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(m1,[27,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(m2,[27,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    m0=0.0_r8
    m1=0.0_r8
    m2=0.0_r8
    shloop: do i=1,sl%n_eps_singlet_shell
        a1=sh%singlet(i)%i1               ! index in the unit cell
        m0=coeffA( (a1-1)*27+1:a1*27,: )  ! partial coefficient matrix
        ! Test that my brain is not turning to mush of old age:
        do j=1,sh%singlet(i)%n_unfold
            iop=sh%singlet(i)%unfold_operation(j)
            a2=sh%singlet(i)%unfold_index(j)
            m27=expandoperation_triplet(sl%singletop(iop)%m3)
            call lo_gemm(m27,m0,m1)
            m2=coeffA( (a2-1)*27+1:a2*27,: )
            m2=m2-m1
            m2=abs(m2)
            f0=sum(m2)
            if ( f0 .gt. lo_sqtol ) then
                call lo_stop_gracefully(['Bad operation singlets'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo
        ! Seems fine, count relevant for this shell:
        l=0
        do j=1,nx
            f0=sum(abs(m0(:,j)))
            if ( f0 .gt. lo_sqtol ) then
                l=l+1
            endif
        enddo
        ! In case of nothing, cheap exit
        if ( l .eq. 0 ) then
            sl%eps_singlet_shell(i)%nx=0
            cycle shloop
        endif
        ! In case of something, store:
        sl%eps_singlet_shell(i)%nx=l
        allocate(sl%eps_singlet_shell(i)%ind_local(l))
        allocate(sl%eps_singlet_shell(i)%ind_global(l))
        allocate(sl%eps_singlet_shell(i)%coeff(27,l))
        sl%eps_singlet_shell(i)%ind_local=0
        sl%eps_singlet_shell(i)%ind_global=0
        sl%eps_singlet_shell(i)%coeff=0.0_r8
        l=0
        do j=1,nx
            f0=sum(abs(m0(:,j)))
            if ( f0 .gt. lo_sqtol ) then
                l=l+1
                sl%eps_singlet_shell(i)%ind_global(l)=j
                sl%eps_singlet_shell(i)%ind_local(l)=j
                sl%eps_singlet_shell(i)%coeff(:,l)=m0(:,j)
            endif
        enddo
    enddo shloop

    ! And store the unfolding
    do i=1,sl%n_eps_singlet_shell
        sl%eps_singlet_shell(i)%n_unfold = sh%singlet(i)%n_unfold
        allocate( sl%eps_singlet_shell(i)%unfold_index( sl%eps_singlet_shell(i)%n_unfold ) )
        allocate( sl%eps_singlet_shell(i)%unfold_operation( sl%eps_singlet_shell(i)%n_unfold ) )
        sl%eps_singlet_shell(i)%unfold_index=sh%singlet(i)%unfold_index
        sl%eps_singlet_shell(i)%unfold_operation=sh%singlet(i)%unfold_operation
    enddo

    ! And a little cleanup
    call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(m1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(m2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> nullspace for dielectric multi-site response, i.e. second order Raman tensor
module subroutine nullspace_eps_pair(sl,sh,uc,ss,cutoff,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> cutoff
    real(r8), intent(in) :: cutoff
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_tuplet_protpair), dimension(:), allocatable :: pair
    integer, dimension(:), allocatable :: ind_to_quartet_op
    integer :: ind_to_transpose

    ! This is a little tricky. The tensor I'm looking at is
    ! eps( ea,eb,uc,ud)
    ! that is, two electric fields and two displacements. The pure pair part should
    ! transform pretty normally, so it's enough to keep the prototypes that are within
    ! the cutoff.
    getpairs: block
        integer :: i,j,k,l,ctr

        ! Get a copy of the prototype pairs that matter. Not sure if it's
        ! really necessary to make this copy, but perhaps. When I want to
        ! generalize later it's a decent idea.
        ctr=0
        do i=1,size(sh%pair)
            if ( norm2(sh%pair(i)%v) .lt. cutoff ) ctr=ctr+1
        enddo
        allocate(pair(ctr))
        ctr=0
        do i=1,size(sh%pair)
            if ( norm2(sh%pair(i)%v) .lt. cutoff ) then
                ctr=ctr+1
                pair(ctr)=sh%pair(ctr)
            endif
        enddo

        ! Figure out which symmetry operations I'm allowed to use. Trick is to find the
        ! quartet operations that correspond to pair operations, if that makes any sense.
        ! Perhaps not, but it does in fact make sense. It's just that a couple of permutations
        ! are not allowed, the pair operations can permute sites, which means I only
        ! want the quartet operations that permute the two last indices. Slightly brain damaging.
        call mem%allocate(ind_to_quartet_op,sl%n_pair_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ind_to_quartet_op=0
        l=0
        do i=1,sl%n_pair_operation
            k=0
            ctr=0
            do j=1,sl%n_quartet_operation
                ! have to be the same rotation
                if ( sl%quartetop(j)%opind .ne. sl%pairop(i)%opind ) cycle
                ! can not permute the electric field indices
                if ( sl%quartetop(j)%perm(1) .ne. 1 ) cycle
                if ( sl%quartetop(j)%perm(2) .ne. 2 ) cycle
                ! next one is annoying.
                if ( sum(abs(sl%pairop(i)%perm-[1,2])) .eq. 0 ) then
                    if ( sum(abs(sl%quartetop(j)%perm(3:4)-[3,4])) .ne. 0 ) cycle
                else
                    if ( sum(abs(sl%quartetop(j)%perm(3:4)-[4,3])) .ne. 0 ) cycle
                endif
                ctr=ctr+1
                k=j
            enddo
            if ( ctr .gt. 1 ) call lo_stop_gracefully(['More than one operation matches, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            if ( k .eq. 0 ) then
                call lo_stop_gracefully(['Could not find operation, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            else
                ind_to_quartet_op(i)=k
            endif
        enddo

        ! Also find the index to the operation that makes a simple transpose of the
        ! electric field indices.
        ind_to_transpose=0
        ctr=0
        do i=1,sl%n_quartet_operation
            if ( sl%quartetop(i)%opind .ne. 1 ) cycle
            if ( sum(abs(sl%quartetop(i)%perm-[2,1,3,4])) .ne. 0 ) cycle
            ind_to_transpose=i
            ctr=ctr+1
        enddo
        if ( ctr .gt. 1 ) call lo_stop_gracefully(['More than one operation matches, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        if ( ind_to_transpose .eq. 0 ) then
            call lo_stop_gracefully(['Could not find operation, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
    end block getpairs

    ! Make a little space for the shells and set some things.
    sl%n_eps_pair_shell=size(pair)
    allocate(sl%eps_pair_shell(sl%n_eps_pair_shell))
    sl%have_eps_pair=.true.
    sl%nx_eps_pair=0

    ! Work out the nullspace of each shell.
    getnspace: block
        real(r8), dimension(:,:,:), allocatable :: op
        !real(r8), dimension(:,:), allocatable :: m0
        integer, dimension(:), allocatable :: di
        integer :: shell,iop,jop,ctr,nop,global_ctr,i

        call mem%allocate(di,sl%n_pair_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=0
        global_ctr=0
        shloop: do shell=1,sl%n_eps_pair_shell
            ! If it's a self-term, we deal with it separately? Probably.
            if ( lo_sqnorm(pair(shell)%v) .lt. lo_sqtol ) then
                sl%eps_pair_shell(shell)%nx=0
                cycle shloop
            endif

            ! Figure out which the invariant operations are?
            di=0
            ctr=0
            do iop=1,sl%n_pair_operation
                if ( mod(iop,mw%n) .ne. mw%r ) cycle
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,pair(shell),pair(shell),sl%pairop(iop),lo_tol) ) then
                    ctr=ctr+1
                    di(iop)=1
                endif
            enddo
            call mw%allreduce('max',di)
            call mw%size_and_offset(ctr,nop)
            if ( nop .gt. 0 ) then
                call mem%allocate(op,[81,81,nop+1],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                op=0.0_r8
            else
                call lo_stop_gracefully(['No invariant operations for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
            ! Then store the operations and sync across ranks
            ctr=0
            do iop=1,sl%n_pair_operation
                if ( mod(iop,mw%n) .ne. mw%r ) cycle
                if ( di(iop) .eq. 0 ) cycle
                ctr=ctr+1
                ! this is the difference from normal pairs, I use triplet operations!
                jop=ind_to_quartet_op(iop)
                op(:,:,ctr+mw%ctr_offset_per_rank(mw%r+1))=sl%quartetop(jop)%sotr
            enddo
            call mw%allreduce('sum',op)
            op(:,:,nop+1)=sl%quartetop(ind_to_transpose)%sotr

            ! Get the nullspace
            call lo_real_nullspace_coefficient_matrix(&
                invariant_operations=op,&
                coeff=sl%eps_pair_shell(shell)%coeff,&
                nvar=sl%eps_pair_shell(shell)%nx,&
                varind=sl%eps_pair_shell(shell)%ind_local,&
                tolerance=lo_sqtol)
            ! Store the global indices
            if ( sl%eps_pair_shell(shell)%nx .gt. 0 ) then
                allocate(sl%eps_pair_shell(shell)%ind_global( sl%eps_pair_shell(shell)%nx ))
                sl%eps_pair_shell(shell)%ind_global=0
                do i=1,sl%eps_pair_shell(shell)%nx
                    global_ctr=global_ctr+1
                    sl%eps_pair_shell(shell)%ind_global(i)=global_ctr
                enddo
            else
                ! Should always be some in this case? Or, not sure.
                call lo_stop_gracefully(['No independent for shell. Not sure if real problem.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
            ! temporary cleanup
            call mem%deallocate(op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        enddo shloop
        ! cleanup
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! Make a note of the total number
        sl%nx_eps_pair=global_ctr

        ! Get the unfolding information
        do shell=1,sl%n_eps_pair_shell
            sl%eps_pair_shell(shell)%n_unfold = pair(shell)%n_unfold
            allocate(sl%eps_pair_shell(shell)%unfold_index( pair(shell)%n_unfold ))
            allocate(sl%eps_pair_shell(shell)%unfold_operation( pair(shell)%n_unfold ))
            sl%eps_pair_shell(shell)%unfold_index    =pair(shell)%unfold_index
            sl%eps_pair_shell(shell)%unfold_operation=pair(shell)%unfold_operation
        enddo
    end block getnspace

    ! Could be the case that everything disappears
    if ( sl%nx_eps_pair .le. 0 ) then
        sl%have_eps_pair=.false.
    endif

    ! And cleanup
    deallocate(pair)
    call mem%deallocate(ind_to_quartet_op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> nullspace for on-site Born charges
module subroutine nullspace_Z_singlet(sl,sh,uc,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(9,9) :: m9
    real(r8), dimension(:,:), allocatable :: iM,rotM,invarM,coeffA
    real(r8), dimension(:,:), allocatable :: m0,m1,m2
    integer :: iop,ne,nx,a1,a2,i,j,l

    ! First get the nullspace for the full thing:
    ne=uc%na*9
    call mem%allocate(iM    ,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(rotM  ,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(invarM,[ne,ne],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    iM=0.0_r8
    rotM=0.0_r8
    invarM=0.0_r8
    call lo_identitymatrix(iM)
    ! Get spacegroup symmetries in there
    do iop=1,sl%n_singlet_operation
        m9=expandoperation_pair(sl%singletop(iop)%m3)
        rotM=0.0_r8
        do a1=1,uc%na
            a2=sl%singletop(iop)%fmap(a1)
            rotM( (a2-1)*9+1:a2*9, (a1-1)*9+1:a1*9 )=m9
        enddo
        invarM=invarM+rotM-IM
    enddo
    ! Then make sure it adds up to zero
    rotM=0.0_r8
    do a1=1,uc%na
        do i=1,9
            rotM( (a1-1)*9+i , i )=1.0_r8
        enddo
    enddo
    invarM=invarM+rotM

    ! This could be a sensible place to add rotational invariances?

    ! Get the nullspace
    call lo_real_nullspace_coefficient_matrix(&
        invarM=invarM,&
        coeff=coeffA,&
        nvar=nx,&
        tolerance=lo_sqtol)
    ! Intermediate cleanup
    call mem%deallocate(iM    ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(rotM  ,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(invarM,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! Could be easy:
    if ( nx .eq. 0 ) then
        sl%nx_Z_singlet=0
        sl%n_Z_singlet_shell=0
        sl%have_Z_singlet=.false.
        return
    endif

    ! Apparently not, now we have to store things
    sl%nx_Z_singlet=nx
    sl%n_Z_singlet_shell=size(sh%singlet)
    sl%have_Z_singlet=.true.
    allocate(sl%Z_singlet_shell(sl%n_Z_singlet_shell))
    ! Do some sanity checks while sorting this out, never a bad idea
    call mem%allocate(m0,[9,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(m1,[9,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(m2,[9,nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    shloop: do i=1,sl%n_Z_singlet_shell
        a1=sh%singlet(i)%i1             ! index in the unit cell
        m0=coeffA( (a1-1)*9+1:a1*9,: )  ! partial coefficient matrix
        ! Test that my brain is not turning to mush of old age:
        do j=1,sh%singlet(i)%n_unfold
            iop=sh%singlet(i)%unfold_operation(j)
            a2=sh%singlet(i)%unfold_index(j)
            m9=expandoperation_pair(sl%singletop(iop)%m3)
            m1=matmul(m9,m0)
            m2=coeffA( (a2-1)*9+1:a2*9,: )
            if ( norm2(m1-m2) .gt. lo_sqtol ) then
                call lo_stop_gracefully(['Bad operation Z singlets'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo
        ! Seems fine, count relevant for this shell:
        l=0
        do j=1,nx
            if ( sum(abs(m0(:,j))) .gt. lo_sqtol ) then
                l=l+1
            endif
        enddo
        ! In case of nothing, cheap exit
        if ( l .eq. 0 ) then
            sl%Z_singlet_shell(i)%nx=0
            cycle shloop
        endif
        ! In case of something, store:
        sl%Z_singlet_shell(i)%nx=l
        allocate(sl%Z_singlet_shell(i)%ind_local(l))
        allocate(sl%Z_singlet_shell(i)%ind_global(l))
        allocate(sl%Z_singlet_shell(i)%coeff(9,l))
        sl%Z_singlet_shell(i)%ind_local=0
        sl%Z_singlet_shell(i)%ind_global=0
        sl%Z_singlet_shell(i)%coeff=0.0_r8
        l=0
        do j=1,nx
            if ( sum(abs(m0(:,j))) .gt. lo_sqtol ) then
                l=l+1
                sl%Z_singlet_shell(i)%ind_global(l)=j
                sl%Z_singlet_shell(i)%ind_local(l)=j
                sl%Z_singlet_shell(i)%coeff(:,l)=m0(:,j)
            endif
        enddo
    enddo shloop

    ! ! And store the unfolding
    do i=1,sl%n_Z_singlet_shell
        sl%Z_singlet_shell(i)%n_unfold = sh%singlet(i)%n_unfold
        allocate( sl%Z_singlet_shell(i)%unfold_index( sl%Z_singlet_shell(i)%n_unfold ) )
        allocate( sl%Z_singlet_shell(i)%unfold_operation( sl%Z_singlet_shell(i)%n_unfold ) )
        sl%Z_singlet_shell(i)%unfold_index=sh%singlet(i)%unfold_index
        sl%Z_singlet_shell(i)%unfold_operation=sh%singlet(i)%unfold_operation
    enddo

    ! And a little cleanup
    call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(m1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(m2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> nullspace for Z', the pair born charges or whatever you call it.
module subroutine nullspace_Z_pair(sl,sh,uc,ss,cutoff,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> cutoff
    real(r8), intent(in) :: cutoff
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_tuplet_protpair), dimension(:), allocatable :: pair
    integer, dimension(:), allocatable :: ind_to_triplet_op

    ! So, the tensor we are working with is
    !   Z(ea,ua,ub)
    ! That is one electric field and two positions. That means transpositions between
    ! index 1 and the other two are not allowed. The pair part is still the normal pair
    ! part, but when checking the invariance under the small group I should use triplet
    ! operations, but only those with the allowed permutations.
    getpairs: block
        integer :: i,j,k,ctr

        ! Get a copy of the prototype pairs that matter. Not sure if it's
        ! really necessary to make this copy, but perhaps. When I want to
        ! generalize later it's a decent idea.
        ctr=0
        do i=1,size(sh%pair)
            if ( norm2(sh%pair(i)%v) .lt. cutoff ) ctr=ctr+1
        enddo
        allocate(pair(ctr))
        ctr=0
        do i=1,size(sh%pair)
            if ( norm2(sh%pair(i)%v) .lt. cutoff ) then
                ctr=ctr+1
                pair(ctr)=sh%pair(ctr)
            endif
        enddo

        ! Figure out which symmetry operations I'm allowed to use. Trick is to find the triplet
        ! operations that correspond to pure pair operations. It makes sense, I promise.
        call mem%allocate(ind_to_triplet_op,sl%n_pair_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ind_to_triplet_op=0
        do i=1,sl%n_pair_operation
            k=0
            ctr=0
            do j=1,sl%n_triplet_operation
                ! have to be the same rotation
                if ( sl%tripletop(j)%opind .ne. sl%pairop(i)%opind ) cycle
                ! can not permute the electric field indices
                if ( sl%tripletop(j)%perm(1) .ne. 1 ) cycle
                ! next one is annoying.
                if ( sum(abs(sl%pairop(i)%perm-[1,2])) .eq. 0 ) then
                    if ( sum(abs(sl%tripletop(j)%perm(2:3)-[2,3])) .ne. 0 ) cycle
                else
                    if ( sum(abs(sl%tripletop(j)%perm(2:3)-[2,3])) .ne. 0 ) cycle
                endif
                ctr=ctr+1
                k=j
            enddo
            if ( ctr .gt. 1 ) call lo_stop_gracefully(['More than one operation matches, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            if ( k .eq. 0 ) then
                call lo_stop_gracefully(['Could not find operation, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            else
                ind_to_triplet_op(i)=k
            endif
        enddo
    end block getpairs

    ! Make a little space for the shells and set some things.
    ! sl%n_Z_pair_shell=size(pair)
    ! allocate(sl%Z_pair_shell(sl%n_Z_pair_shell))
    ! sl%have_Z_pair=.true.
    ! sl%nx_Z_pair=0

    ! Work out the nullspace of each shell.
    getnspace: block
        type(lo_tensor_shell), dimension(:), allocatable :: dsh
        real(r8), dimension(:,:,:), allocatable :: op
        real(r8), dimension(:,:), allocatable :: m0
        integer, dimension(:), allocatable :: di
        integer :: shell,iop,jop,ctr,nop,global_ctr,i

        ! Create dummy shells
        allocate(dsh(size(pair)))

        call mem%allocate(di,sl%n_pair_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=0
        shloop: do shell=1,size(pair)
            ! If it's a self-term, we deal with it separately? Probably.
            if ( lo_sqnorm(pair(shell)%v) .lt. lo_sqtol ) then
                dsh(shell)%nx=0
                !sl%Z_pair_shell(shell)%nx=0
                cycle shloop
            endif

            ! Figure out which the invariant operations are?
            di=0
            ctr=0
            do iop=1,sl%n_pair_operation
                if ( mod(iop,mw%n) .ne. mw%r ) cycle
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,pair(shell),pair(shell),sl%pairop(iop),lo_tol) ) then
                    ctr=ctr+1
                    di(iop)=1
                endif
            enddo
            call mw%allreduce('max',di)
            call mw%size_and_offset(ctr,nop)

            if ( nop .gt. 0 ) then
                call mem%allocate(op,[27,27,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                op=0.0_r8
            else
                call lo_stop_gracefully(['No invariant operations for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif

            ! Then store the operations and sync across ranks
            ctr=0
            do iop=1,sl%n_pair_operation
                if ( mod(iop,mw%n) .ne. mw%r ) cycle
                if ( di(iop) .eq. 0 ) cycle
                ctr=ctr+1
                ! this is the difference from normal pairs, I use triplet operations!
                jop=ind_to_triplet_op(iop)
                op(:,:,ctr+mw%ctr_offset_per_rank(mw%r+1))=sl%tripletop(jop)%sotr
            enddo
            call mw%allreduce('sum',op)
            ! Get the nullspace
            call lo_real_nullspace_coefficient_matrix(&
                invariant_operations=op,&
                coeff=dsh(shell)%coeff,&
                nvar=dsh(shell)%nx,&
                varind=dsh(shell)%ind_local,&
                tolerance=lo_sqtol)

            if ( dsh(shell)%nx .gt. 0 ) then
                ! Test if I am still worthy
                call mem%allocate(m0,[27,dsh(shell)%nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                m0=0.0_r8
                do iop=1,size(op,3)
                    call lo_gemm(op(:,:,iop),dsh(shell)%coeff,m0)
                    if ( norm2(m0-dsh(shell)%coeff) .gt. lo_sqtol*size(m0) ) then
                        call lo_stop_gracefully(['Not as invariant as it should be. I am not worthy.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                    endif
                enddo
                call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            endif

            ! Make a note that triplets are a thing
            sl%consider_triplet=.true.

            ! temporary cleanup
            call mem%deallocate(op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        enddo shloop

        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! Now make sure to store the shells that have a nonzero number of components, also means
        ! I only store the pairs with nonzero contributions.
        ctr=0
        i=0
        do shell=1,size(pair)
            if ( dsh(shell)%nx .gt. 0 ) then
                ctr=ctr+1
                i=i+1
            endif
            if ( lo_sqnorm(pair(shell)%v) .lt. lo_sqtol ) ctr=ctr+1
        enddo

        ! Now return the correct shells
        if ( i .gt. 0 ) then
            sl%n_Z_pair_shell=ctr
            allocate(sl%Z_pair_shell(sl%n_Z_pair_shell))
            sl%nx_Z_pair=0
            ctr=0
            global_ctr=0
            ! Means I have meaningful contributions, store the relevant ones.
            do shell=1,size(pair)
                if ( lo_sqnorm(pair(shell)%v) .lt. lo_sqtol ) then
                    ctr=ctr+1
                elseif ( dsh(shell)%nx .gt. 0 ) then
                    ctr=ctr+1
                else
                    ! irrelevant shell
                    cycle
                endif
                ! Store the shell
                sl%Z_pair_shell(ctr)=dsh(shell)
                ! Store pairs
                sl%Z_pair_shell(ctr)%n_unfold = pair(shell)%n_unfold
                allocate(sl%Z_pair_shell(ctr)%unfold_index( pair(shell)%n_unfold ))
                allocate(sl%Z_pair_shell(ctr)%unfold_operation( pair(shell)%n_unfold ))
                sl%Z_pair_shell(ctr)%unfold_index    =pair(shell)%unfold_index
                sl%Z_pair_shell(ctr)%unfold_operation=pair(shell)%unfold_operation
                ! Store indices
                if ( sl%Z_pair_shell(ctr)%nx .gt. 0 ) then
                    allocate(sl%Z_pair_shell(ctr)%ind_global( sl%Z_pair_shell(ctr)%nx ))
                    sl%Z_pair_shell(ctr)%ind_global=0
                    do i=1,sl%Z_pair_shell(ctr)%nx
                        global_ctr=global_ctr+1
                        sl%Z_pair_shell(ctr)%ind_global(i)=global_ctr
                    enddo
                endif
            enddo
            sl%nx_Z_pair=global_ctr
            sl%have_Z_pair=.true.
        else
            ! Everything disappeared with symmetry
            sl%nx_Z_pair=0
            sl%n_Z_pair_shell=0
            sl%have_Z_pair=.false.
        endif

        deallocate(dsh)
    end block getnspace

    ! And cleanup
    deallocate(pair)
    call mem%deallocate(ind_to_triplet_op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> nullspace for Z''
module subroutine nullspace_Z_triplet(sl,sh,uc,ss,cutoff,mw,mem)
    !> symmetry table
    type(lo_interaction_tensors), intent(inout) :: sl
    !> helper
    type(lo_symtabhelper), intent(in) :: sh
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> cutoff
    real(r8), intent(in) :: cutoff
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_tuplet_prottriplet), dimension(:), allocatable :: triplet
    integer, dimension(:), allocatable :: ind_to_quartet_op

    ! So, the tensor we are working with is
    !   Z(ea,ua,ub,uc)
    ! So, transpose the last three but not the first index.
    gettriplets: block
        integer :: i,j,k,ctr

        ! Get a copy of the prototype triplets that matter.
        ctr=0
        do i=1,size(sh%triplet)
            if ( norm2(sh%triplet(i)%v2) .gt. cutoff ) cycle
            if ( norm2(sh%triplet(i)%v3) .gt. cutoff ) cycle
            if ( norm2(sh%triplet(i)%v2-sh%triplet(i)%v3) .gt. cutoff ) cycle
            ctr=ctr+1
        enddo
        allocate(triplet(ctr))
        ctr=0
        do i=1,size(sh%triplet)
            if ( norm2(sh%triplet(i)%v2) .gt. cutoff ) cycle
            if ( norm2(sh%triplet(i)%v3) .gt. cutoff ) cycle
            if ( norm2(sh%triplet(i)%v2-sh%triplet(i)%v3) .gt. cutoff ) cycle
            ctr=ctr+1
            triplet(ctr)=sh%triplet(ctr)
        enddo

        ! Figure out which symmetry operations I'm allowed to use. Trick is to find the triplet
        ! operations that correspond to pure triplet operations. It makes sense, I promise.
        call mem%allocate(ind_to_quartet_op,sl%n_triplet_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ind_to_quartet_op=0
        do i=1,sl%n_triplet_operation
            k=0
            ctr=0
            do j=1,sl%n_quartet_operation
                ! have to be the same rotation
                if ( sl%quartetop(j)%opind .ne. sl%tripletop(i)%opind ) cycle
                ! can not permute the electric field index
                if ( sl%quartetop(j)%perm(1) .ne. 1 ) cycle
                ! now it gets annoying. If the triplet permutation is
                ! [1,2,3], then the quartet operation must be [2,3,4]. I think.
                if ( sum(abs(sl%tripletop(i)%perm-(sl%quartetop(j)%perm(2:4)-1))) .eq. 0 ) then
                    ctr=ctr+1
                    k=j
                endif
            enddo

            if ( ctr .gt. 1 ) call lo_stop_gracefully(['More than one operation matches, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            if ( k .eq. 0 ) then
                call lo_stop_gracefully(['Could not find operation, impossible'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            else
                ind_to_quartet_op(i)=k
            endif
        enddo
    end block gettriplets

    ! Work out the nullspace of each shell.
    getnspace: block
        type(lo_tensor_shell), dimension(:), allocatable :: dsh
        real(r8), dimension(:,:,:), allocatable :: op
        real(r8), dimension(:,:), allocatable :: m0
        integer, dimension(:), allocatable :: di
        integer :: shell,iop,jop,ctr,nop,global_ctr,i

        ! Get some dummy shells
        allocate(dsh(size(triplet)))

        call mem%allocate(di,sl%n_triplet_operation,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=0
        global_ctr=0
        shloop: do shell=1,size(triplet)
            ! If it's a self-term, we deal with it separately? Probably.
            if ( lo_sqnorm(triplet(shell)%v2)+lo_sqnorm(triplet(shell)%v3) .lt. lo_sqtol ) then
                dsh(shell)%nx=0
                !sl%Z_triplet_shell(shell)%nx=0
                cycle shloop
            endif

            ! Figure out which the invariant operations are?
            di=0
            ctr=0
            do iop=1,sl%n_triplet_operation
                if ( mod(iop,mw%n) .ne. mw%r ) cycle
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,triplet(shell),triplet(shell),sl%tripletop(iop),lo_tol) ) then
                    ctr=ctr+1
                    di(iop)=1
                endif
            enddo
            call mw%allreduce('max',di)
            call mw%size_and_offset(ctr,nop)

            if ( nop .gt. 0 ) then
                call mem%allocate(op,[81,81,nop],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                op=0.0_r8
            else
                call lo_stop_gracefully(['No invariant operations for shell'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif

            ! Then store the operations and sync across ranks
            ctr=0
            do iop=1,sl%n_triplet_operation
                if ( mod(iop,mw%n) .ne. mw%r ) cycle
                if ( di(iop) .eq. 0 ) cycle
                ctr=ctr+1
                ! this is the difference from normal triplets, I use quartet operations!
                jop=ind_to_quartet_op(iop)
                op(:,:,ctr+mw%ctr_offset_per_rank(mw%r+1))=sl%quartetop(jop)%sotr
            enddo
            call mw%allreduce('sum',op)
            ! Get the nullspace
            call lo_real_nullspace_coefficient_matrix(&
                invariant_operations=op,&
                coeff=dsh(shell)%coeff,&
                nvar=dsh(shell)%nx,&
                varind=dsh(shell)%ind_local,&
                tolerance=lo_sqtol)

            if ( dsh(shell)%nx .gt. 0 ) then
                ! Test if I am still worthy
                call mem%allocate(m0,[81,dsh(shell)%nx],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                m0=0.0_r8
                do iop=1,size(op,3)
                    call lo_gemm(op(:,:,iop),dsh(shell)%coeff,m0)
                    if ( norm2(m0-dsh(shell)%coeff) .gt. lo_sqtol*size(m0) ) then
                        call lo_stop_gracefully(['Not as invariant as it should be. I am not worthy.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                    endif
                enddo
                call mem%deallocate(m0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            endif

            ! Make a note that quartets are a thing if it was not really obvious
            sl%consider_quartet=.true.

            ! temporary cleanup
            call mem%deallocate(op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        enddo shloop

        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! Now, only store the shells with nonzero number of independent.
        ctr=0
        i=0
        do shell=1,size(triplet)
            if ( dsh(shell)%nx .gt. 0 ) then
                ctr=ctr+1
                i=i+1
            endif
            if ( lo_sqnorm(triplet(shell)%v2)+lo_sqnorm(triplet(shell)%v3) .lt. lo_sqtol ) ctr=ctr+1
        enddo

        if ( i .gt. 0 ) then
            ! This means I have significant shells
            sl%n_Z_triplet_shell=ctr
            allocate(sl%Z_triplet_shell(sl%n_Z_triplet_shell))
            sl%nx_Z_triplet=0
            ctr=0
            global_ctr=0
            ! Means I have meaningful contributions, store the relevant ones.
            do shell=1,size(triplet)
                if ( lo_sqnorm(triplet(shell)%v2)+lo_sqnorm(triplet(shell)%v3) .lt. lo_sqtol ) then
                    ctr=ctr+1
                elseif ( dsh(shell)%nx .gt. 0 ) then
                    ctr=ctr+1
                else
                    ! irrelevant shell
                    cycle
                endif
                ! Store the shell
                sl%Z_triplet_shell(ctr)=dsh(shell)
                ! Store triplets
                sl%Z_triplet_shell(ctr)%n_unfold = triplet(shell)%n_unfold
                allocate(sl%Z_triplet_shell(ctr)%unfold_index( triplet(shell)%n_unfold ))
                allocate(sl%Z_triplet_shell(ctr)%unfold_operation( triplet(shell)%n_unfold ))
                sl%Z_triplet_shell(ctr)%unfold_index    =triplet(shell)%unfold_index
                sl%Z_triplet_shell(ctr)%unfold_operation=triplet(shell)%unfold_operation
                ! Store indices
                if ( sl%Z_triplet_shell(ctr)%nx .gt. 0 ) then
                    allocate(sl%Z_triplet_shell(ctr)%ind_global( sl%Z_triplet_shell(ctr)%nx ))
                    sl%Z_triplet_shell(ctr)%ind_global=0
                    do i=1,sl%Z_triplet_shell(ctr)%nx
                        global_ctr=global_ctr+1
                        sl%Z_triplet_shell(ctr)%ind_global(i)=global_ctr
                    enddo
                endif
            enddo
            sl%nx_Z_triplet=global_ctr
            sl%have_Z_triplet=.true.
        else
            ! Everything disappeared with symmetry
            sl%nx_Z_triplet=0
            sl%n_Z_triplet_shell=0
            sl%have_Z_triplet=.false.
        endif
    end block getnspace

    ! And cleanup
    deallocate(triplet)
    call mem%deallocate(ind_to_quartet_op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

! !> nullspace for magnetic exchange pair interactions
! module subroutine nullspace_magnetic_pair(sl,sh,uc,cutoff,mw,mem)
!     !> symmetry table
!     type(lo_interaction_tensors), intent(inout) :: sl
!     !> helper
!     type(lo_symtabhelper), intent(in) :: sh
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!     !> cutoff
!     real(r8), intent(in) :: cutoff
!     !> mpi helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> memory helper
!     type(lo_mem_helper), intent(inout) :: mem
!
!     type(lo_tuplet_protpair), dimension(:), allocatable :: pair
!
!     ! Get the list of pairs to reduce.
!     getpairs: block
!         integer :: i,ctr
!
!         ! Get a copy of the prototype pairs that matter. Not sure if it's
!         ! really necessary to make this copy, but perhaps. When I want to
!         ! generalize later it's a decent idea.
!         ctr=0
!         do i=1,size(sh%pair)
!             if ( norm2(sh%pair(i)%v) .lt. cutoff ) ctr=ctr+1
!         enddo
!         allocate(pair(ctr))
!         ctr=0
!         do i=1,size(sh%pair)
!             if ( norm2(sh%pair(i)%v) .lt. cutoff ) then
!                 ctr=ctr+1
!                 pair(ctr)=sh%pair(ctr)
!             endif
!         enddo
!         write(*,*) 'FIXME MAGNETIC JIJ'
!         stop
!     end block getpairs
! end subroutine

end submodule
