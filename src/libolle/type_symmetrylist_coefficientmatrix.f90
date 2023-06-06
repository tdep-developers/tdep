#include "precompilerdefinitions"
submodule (type_symmetrylist) type_symmetrylist_coefficientmatrix
implicit none
contains

!> Create the nullspace for each irreducible singlet shell @TODO replicate the one from forces here, works better.
module subroutine nullspace_singletshells(sl,sh)
    type(lo_symlist), intent(inout) :: sl
    type(lo_symtabhelper), intent(inout) :: sh

    integer, dimension(:), allocatable :: opind
    integer :: shell,op,opctr
    integer :: i,j,l,a1
    real(flyt), dimension(:,:,:), allocatable :: listofM
    real(flyt), dimension(3,3) :: I3,invarM

    call lo_identitymatrix(I3)
    lo_allocate(listofM(3,3,sl%sym%n))
    lo_allocate(opind(sl%sym%n))

    shloop: do shell=1,sl%nsingletshells
        ! Now I assume it's a rather normal shell. Figure out which operations
        ! leave this atom invariant
        listofM=0.0_flyt
        opctr=0
        a1=sl%singletshell(shell)%protatom
        oploop: do op=1,sl%sym%n
            ! is it ok?
            if ( sh%clist(a1,a1,op) ) then
                ! this might be a new, interesting operation
                do i=1,opctr
                    if ( sum(abs(listofM(:,:,i)-sl%sym%op(op)%m)) .lt. lo_tol ) then
                        ! sorry, not new, cycle
                        cycle oploop
                    endif
                enddo
                ! yup, it was new!
                opctr=opctr+1
                listofM(:,:,opctr)=sl%sym%op(op)%m
            endif
        enddo oploop

        ! Add it up
        invarM=0.0_flyt
        do i=1,opctr
            invarM=invarM+listofM(:,:,i)-I3
        enddo
        call lo_nullspace_coefficient_matrix(invarM,sl%singletshell(shell)%coeff_M,&
                                             sl%singletshell(shell)%nfc,sl%singletshell(shell)%relfc)
        if ( sl%singletshell(shell)%nfc .gt. 0 ) then
            lo_allocate(sl%singletshell(shell)%relfcind( sl%singletshell(shell)%nfc ))
            sl%singletshell(shell)%relfcind=0
        endif
    enddo shloop

    ! count forceconstants and index them
    l=0
    do i=1,sl%nsingletshells
        ! count the others
        do j=1,sl%singletshell(i)%nfc
            l=l+1
            sl%singletshell(i)%relfcind(j)=l
        enddo
    enddo
    sl%nsingletifc=l
end subroutine

!> Create the nullspace for the Born effective charges
module subroutine nullspace_Zshells(sl,uc)
    type(lo_symlist), intent(inout) :: sl
    type(lo_crystalstructure), intent(in) :: uc

    real(flyt), dimension(:,:), allocatable :: rotM,invarM,IM
    real(flyt), dimension(9,9) :: m9
    integer :: ne,i,a1,a2,o

    ne=uc%na*9
    lo_allocate(invarM(ne,ne))
    lo_allocate(rotM(ne,ne))
    lo_allocate(IM(ne,ne))
    invarM=0.0_flyt
    rotM=0.0_flyt
    call lo_identitymatrix(IM)

    ! These are the spacegroup thingies, and the fmap takes care of the transform between the atoms.
    do o=1,uc%sym%n
        ! Expand the operation
        m9=lo_expandoperation_pair(uc%sym%op(o)%m)
        rotM=0.0_flyt
        do a1=1,uc%na
            a2=uc%sym%op(o)%fmap(a1)
            rotM( (a2-1)*9+1:a2*9, (a1-1)*9+1:a1*9 )=m9
        enddo
        invarM=invarM+rotM-IM
    enddo
    ! Then charge neutrality
    rotM=0.0_flyt
    do a1=1,uc%na
        do i=1,9
            rotM( (a1-1)*9+i , i )=1.0_flyt
        enddo
    enddo
    invarM=invarM+rotM
    ! Project out the nullspace
    call lo_nullspace_coefficient_matrix(invarM,sl%Zshell%coeff_M,sl%nZtheta,tolerance=lo_sqtol)
end subroutine

!> Create the nullspace for the singlet eps guys
module subroutine nullspace_epssinglet(sl,uc)
    type(lo_symlist), intent(inout) :: sl
    type(lo_crystalstructure), intent(in) :: uc

    real(flyt), dimension(:,:), allocatable :: rotM,invarM,IM
    real(flyt), dimension(27,27,6) :: trm_triplet
    integer, dimension(3,6) :: prm_triplet
    real(flyt), dimension(27,27) :: m27
    integer :: ne,i,a1,a2,o,ipr

    ! Collect the relevant permutation
    call lo_return_tensor_transpositions(trm_triplet=trm_triplet,prm_triplet=prm_triplet)
    ipr=0
    do i=1,6
        if ( prm_triplet(3,i) .eq. 3 .and. prm_triplet(1,i) .eq. 2 ) ipr=i
    enddo

    ne=uc%na*27
    allocate(invarM(ne,ne))
    allocate(rotM(ne,ne))
    allocate(IM(ne,ne))
    invarM=0.0_flyt
    rotM=0.0_flyt
    call lo_identitymatrix(IM)

    ! These are the spacegroup thingies, and the fmap takes care of the transform between the atoms.
    symloop: do o=1,uc%sym%n
        ! Expand the operation
        m27=expandoperation_triplet(uc%sym%op(o)%m)
        rotM=0.0_flyt
        do a1=1,uc%na
            a2=uc%sym%op(o)%fmap(a1)
            rotM( (a2-1)*27+1:a2*27, (a1-1)*27+1:a1*27 )=m27
        enddo
        invarM=invarM+rotM-IM
    enddo symloop
    ! Add the permutation symmetry
    rotM=0.0_flyt
    do a1=1,uc%na
        rotM( (a1-1)*27+1:a1*27 , (a1-1)*27+1:a1*27 )=trm_triplet(:,:,ipr)
    enddo
    invarM=invarM+rotM-IM
    ! Then charge neutrality
    rotM=0.0_flyt
    do a1=1,uc%na
        do i=1,27
            rotM( (a1-1)*27+i , i )=1.0_flyt
        enddo
    enddo
    invarM=invarM+rotM
    ! Project out the nullspace
    call lo_nullspace_coefficient_matrix(invarM,sl%epssingletshell%coeff_M,sl%nepstheta_singlet,tolerance=lo_sqtol)
end subroutine

!> Create the nullspace for the dielectric tensor
module subroutine nullspace_epsshell(sl,uc)
    !> symmetry stuff
    type(lo_symlist), intent(inout) :: sl
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc

    real(flyt), dimension(:,:,:), allocatable :: ops
    real(flyt), dimension(9,9) :: I9,T,invarM
    integer :: i

    ! Transposition operator
    call lo_transpositionmatrix(T)
    call lo_identitymatrix(I9)
    ! Expanded operations
    lo_allocate(ops(9,9,uc%sym%n))
    ! Expand the operations
    do i=1,uc%sym%n
        ops(:,:,i)=lo_expandoperation_pair( uc%sym%op(i)%m )
    enddo
    ops=lo_chop(ops,1E-11_flyt)

    ! Get the coefficient matrix
    invarM=0.0_flyt
    do i=1,uc%sym%n
        invarM=invarM+ops(:,:,i)-I9
    enddo
    if ( sum(abs(invarM)) .gt. lo_tol ) then
        invarM=invarM+T-I9
    endif
    call lo_nullspace_coefficient_matrix(invarM,sl%epsshell%coeff_M,sl%nepstheta,tolerance=lo_sqtol)
end subroutine

! ! Create the nullspace for each irreducible shell
! module subroutine nullspace_pairshells(sl,sh,ss)
!     type(lo_symlist), intent(inout) :: sl
!     type(lo_symtabhelper), intent(inout) :: sh
!     type(lo_crystalstructure), intent(in) :: ss
!
!     integer, dimension(:), allocatable :: opind
!     integer, dimension(9) :: zeroSV
!     integer :: shell,op,opctr,fcctr,fcctr_global
!     integer :: i,j,l
!     real(flyt), dimension(:,:,:), allocatable :: listofM
!     real(flyt), dimension(9,9) :: I9,BMS,U,V,coeffM,invarM
!     real(flyt), dimension(9) :: S
!
!     call lo_identitymatrix(I9)
!     lo_allocate(listofM(9,9,sl%npairop))
!     lo_allocate(opind(sl%npairop))
!
!     fcctr_global=0
!     shloop: do shell=1,sl%npairshells
!
!         ! it might be a self-term, in that case it's near trivial to fix
!         if ( lo_sqnorm(sl%pairshell(shell)%protpair%v) .lt. sl%sqtol_cart ) then
!             sl%pairshell(shell)%nfc=0
!             sl%pairshell(shell)%coeffM=0.0_flyt
!             sl%pairshell(shell)%invarM=0.0_flyt
!             cycle shloop
!         endif
!
!         ! Now I assume it's a rather normal shell. Figure out which operations
!         ! take the pair back to itself.
!         listofM=0.0_flyt
!         opctr=0
!         oploop: do op=1,sl%npairop
!             if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sl%pairshell(shell)%protpair,sl%pairshell(shell)%protpair,sl%pairop(op)) ) then
!                 ! this might be a new, interesting operation
!                 do i=1,opctr
!                     if ( sum(abs(listofM(:,:,i)-sl%pairop(op)%sotr)) .lt. lo_tol ) then
!                         ! sorry, not new, cycle
!                         cycle oploop
!                     endif
!                 enddo
!                 ! yup, it was new!
!                 opctr=opctr+1
!                 listofM(:,:,opctr)=sl%pairop(op)%sotr
!             endif
!         enddo oploop
!
!         ! Seems reasonable, perhaps. Build the thingy
!         BMS=0.0_flyt
!         S=0.0_flyt
!         U=0.0_flyt
!         V=0.0_flyt
!         coeffM=0.0_flyt
!         zeroSV=0
!         fcctr=0
!         do i=1,opctr
!             BMS=BMS+listofM(:,:,i)-I9
!         enddo
!         invarM=BMS
!         !call lo_dgesvd(BMS,S,U,V)
!         call lo_dgesdd(BMS,S,U,V,jobz='A',info=lo_status)
!         l=0
!         do i=1,9
!             if ( abs(S(i)) .lt. lo_tol ) then
!                 l=l+1
!                 zeroSV(l)=i
!             endif
!         enddo
!         if ( l .eq. 0 .or. l .eq. 9 ) then
!             ! symmetry did nothing for this shell
!             coeffM=I9
!             fcctr=9
!         else
!             ! symmetry did something, generate that and clean it a bit
!             fcctr=l
!             coeffM=matmul(U(:,zeroSV(1:fcctr)),V(zeroSV(1:fcctr),:))
!             !call lo_real_gram_schmidt(coeffM)
!             call lo_make_coeffmatrix_tidy(coeffM)
!         endif
!         ! store stuff in the shell
!         sl%pairshell(shell)%invarM=invarM
!         sl%pairshell(shell)%coeffM=coeffM
!         sl%pairshell(shell)%nfc=fcctr
!         lo_allocate(sl%pairshell(shell)%relfc( fcctr ))
!         lo_allocate(sl%pairshell(shell)%relfcind( fcctr ))
!         l=0
!         do i=1,9
!             if ( sum(abs(coeffM(:,i))) .gt. lo_sqtol ) then
!                 l=l+1
!                 sl%pairshell(shell)%relfc(l)=i
!             endif
!         enddo
!
!         ! small sanity check
!         if ( l .ne. fcctr ) then
!             write(*,*) 'ERROR:'
!             write(*,*) '    my Gram-Schmidt or something failed, I got',l,'ifcs when I expected',fcctr
!         endif
!
!         fcctr_global=fcctr_global+fcctr
!
! #ifdef AGRESSIVE_SANITY
! ! Do some more, agressive sanity tests
! sanity: block
!     real(flyt), dimension(9,9) :: m0,m1
!     real(flyt), dimension(9) :: v0,v1
!     real(flyt) :: f0
!
!     ! Check if the SVD was ok, with orthogonal vectors and all that.
!     m0=transpose(V)
!     f0=0.0_flyt
!     do i=1,9
!     do j=i+1,9
!         f0=f0+abs(dot_product(U(:,i),U(:,j)))+abs(dot_product(m0(:,i),m0(:,j)))
!     enddo
!     enddo
!     ! Check the norm of the singular vectors, have to be 1
!     do i=1,9
!         f0=f0+norm2(U(:,i))+norm2(m0(:,i))-2.0_flyt
!     enddo
!     ! The singular vectors are an orthogonal transformation, have to preserve the norm properly
!     do i=1,9
!         v0(i)=i
!     enddo
!     v1=matmul(U,v0)
!     f0=f0+norm2(v1)-norm2(v0)
!     v1=matmul(m0,v0)
!     f0=f0+norm2(v1)-norm2(v0)
!
!     ! Check that the coefficient matrix is truly invariant to the operations it should
!     ! be invariant under
!     m0=sl%pairshell(shell)%coeffM
!     do i=1,opctr
!         call lo_gemm(listofM(:,:,i),m0,m1)
!         f0=f0+sum(abs(m0-m1))
!     enddo
!     if ( f0 .gt. lo_sqtol ) then
!         call lo_stop_gracefully(['Inconsitent symmetry for pair shell '//tochar(shell)],lo_exitcode_symmetry,__FILE__,__LINE__)
!     endif
! end block sanity
! #endif
!
!     enddo shloop
!
!     ! count forceconstants and index them
!     l=0
!     do i=1,sl%npairshells
!         ! count the others
!         do j=1,sl%pairshell(i)%nfc
!             l=l+1
!             sl%pairshell(i)%relfcind(j)=l
!         enddo
!     enddo
!     sl%npairifc=l
!
! end subroutine

! Create the nullspace for each irreducible shell
module subroutine nullspace_tripletshells(sl,sh,ss)
    type(lo_symlist), intent(inout) :: sl
    type(lo_symtabhelper), intent(inout) :: sh
    type(lo_crystalstructure), intent(in) :: ss

    integer, parameter :: dm=27
    integer, dimension(:), allocatable :: opind
    integer, dimension(dm) :: zeroSV
    integer :: shell,op,opctr,fcctr,fcctr_global
    integer :: i,j,l
    real(flyt), dimension(:,:,:), allocatable :: listofM
    real(flyt), dimension(dm,dm) :: I27,BMS,U,V,coeffM
    real(flyt), dimension(dm) :: S
    real(flyt) :: f0

    call lo_identitymatrix(I27)
    lo_allocate(listofM(dm,dm,sl%ntripletop))
    lo_allocate(opind(sl%ntripletop))

    fcctr_global=0
    shloop: do shell=1,sl%ntripletshells

        ! it might be a self-term, in that case it's near trivial to fix
        f0=lo_sqnorm(sl%tripletshell(shell)%prottriplet%v1)+&
           lo_sqnorm(sl%tripletshell(shell)%prottriplet%v2)+&
           lo_sqnorm(sl%tripletshell(shell)%prottriplet%v3)
        if ( f0 .lt. sl%sqtol_cart ) then
            sl%tripletshell(shell)%nfc=0
            sl%tripletshell(shell)%coeffM=0.0_flyt
            cycle shloop
        endif

        ! Now I assume it's a rather normal shell. Figure out which operations
        ! take the pair back to itself.
        listofM=0.0_flyt
        opctr=0
        oploop: do op=1,sl%ntripletop
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sl%tripletshell(shell)%prottriplet,sl%tripletshell(shell)%prottriplet,sl%tripletop(op)) ) then
                ! this might be a new, interesting operation
                do i=1,opctr
                    if ( sum(abs(listofM(:,:,i)-sl%tripletop(op)%sotr)) .lt. lo_tol ) then
                        ! sorry, not new, cycle
                        cycle oploop
                    endif
                enddo
                ! yup, it was new!
                opctr=opctr+1
                listofM(:,:,opctr)=sl%tripletop(op)%sotr
            endif
        enddo oploop

        ! Seems reasonable, perhaps. Build the thingy
        BMS=0.0_flyt
        S=0.0_flyt
        U=0.0_flyt
        V=0.0_flyt
        coeffM=0.0_flyt
        zeroSV=0
        fcctr=0
        do i=1,opctr
            BMS=BMS+listofM(:,:,i)-I27
        enddo
        !call lo_dgesvd(BMS,S,U,V)
        call lo_dgesdd(BMS,S,U,V,jobz='A',info=lo_status)
        l=0
        do i=1,dm
            if ( abs(S(i)) .lt. lo_tol ) then
                l=l+1
                zeroSV(l)=i
            endif
        enddo
        if ( l .eq. 0 .or. l .eq. dm ) then
            ! symmetry did nothing for this shell
            coeffM=I27
            fcctr=dm
        else
            ! symmetry did something, generate that and clean it a bit
            fcctr=l
            coeffM=matmul(U(:,zeroSV(1:fcctr)),V(zeroSV(1:fcctr),:))
            !call lo_real_gram_schmidt(coeffM)
            call lo_make_coeffmatrix_tidy(coeffM)
        endif
        ! store stuff in the shell
        sl%tripletshell(shell)%coeffM=coeffM
        sl%tripletshell(shell)%nfc=fcctr
        lo_allocate(sl%tripletshell(shell)%relfc( fcctr ))
        lo_allocate(sl%tripletshell(shell)%relfcind( fcctr ))
        l=0
        do i=1,dm
            if ( sum(abs(coeffM(:,i))) .gt. lo_sqtol ) then
                l=l+1
                sl%tripletshell(shell)%relfc(l)=i
            endif
        enddo

        ! small sanity check
        if ( l .ne. fcctr ) then
            write(*,*) 'ERROR:'
            write(*,*) '    my Gram-Schmidt or something failed, I got',l,'ifcs when I expected',fcctr
        endif

        fcctr_global=fcctr_global+fcctr

#ifdef AGRESSIVE_SANITY
! Do some more, agressive sanity tests
sanity: block
    real(flyt), dimension(27,27) :: m0,m1
    real(flyt) :: f0

    ! Check if the SVD was ok, with orthogonal vectors and all that.
    m0=transpose(V)
    f0=0.0_flyt
    do i=1,27
    do j=i+1,27
        f0=f0+abs(dot_product(U(:,i),U(:,j)))+abs(dot_product(m0(:,i),m0(:,j)))
    enddo
    enddo

    ! Check that the coefficient matrix is truly invariant to the operations it should
    ! be invariant under
    m0=sl%tripletshell(shell)%coeffM
    do i=1,opctr
        call lo_gemm(listofM(:,:,i),m0,m1)
        f0=f0+sum(abs(m0-m1))
    enddo
    if ( f0 .gt. lo_sqtol ) then
        call lo_stop_gracefully(['Inconsitent symmetry for triplet shell '//tochar(shell)],lo_exitcode_symmetry,__FILE__,__LINE__)
    endif
end block sanity
#endif

    enddo shloop

    ! count forceconstants and index them
    l=0
    do i=1,sl%ntripletshells
        ! count the others
        do j=1,sl%tripletshell(i)%nfc
            l=l+1
            sl%tripletshell(i)%relfcind(j)=l
        enddo
    enddo
    sl%ntripletifc=l
end subroutine

! Create the nullspace for each irreducible shell
module subroutine nullspace_quartetshells(sl,sh,ss)
    type(lo_symlist), intent(inout) :: sl
    type(lo_symtabhelper), intent(inout) :: sh
    type(lo_crystalstructure), intent(in) :: ss

    integer, parameter :: dm=81
    integer, dimension(:), allocatable :: opind
    integer, dimension(dm) :: zeroSV
    integer :: shell,op,opctr,fcctr,fcctr_global
    integer :: i,j,l
    real(flyt), dimension(:,:,:), allocatable :: listofM
    real(flyt), dimension(dm,dm) :: I81,BMS,U,V,coeffM
    real(flyt), dimension(dm) :: S
    real(flyt) :: f0

    call lo_identitymatrix(I81)
    lo_allocate(listofM(dm,dm,sl%nquartetop))
    lo_allocate(opind(sl%nquartetop))

    fcctr_global=0
    shloop: do shell=1,sl%nquartetshells

        ! it might be a self-term, in that case it's near trivial to fix
        f0=lo_sqnorm(sl%quartetshell(shell)%protquartet%v1)+&
           lo_sqnorm(sl%quartetshell(shell)%protquartet%v2)+&
           lo_sqnorm(sl%quartetshell(shell)%protquartet%v3)+&
           lo_sqnorm(sl%quartetshell(shell)%protquartet%v4)
        if ( f0 .lt. sl%sqtol_cart ) then
            sl%quartetshell(shell)%nfc=0
            sl%quartetshell(shell)%coeffM=0.0_flyt
            cycle shloop
        endif

        ! Now I assume it's a rather normal shell. Figure out which operations take the pair back to itself.
        listofM=0.0_flyt
        opctr=0
        oploop: do op=1,sl%nquartetop
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sl%quartetshell(shell)%protquartet,sl%quartetshell(shell)%protquartet,sl%quartetop(op)) ) then
                ! this might be a new, interesting operation
                do i=1,opctr
                    if ( sum(abs(listofM(:,:,i)-sl%quartetop(op)%sotr)) .lt. lo_tol ) then
                        ! sorry, not new, cycle
                        cycle oploop
                    endif
                enddo
                ! yup, it was new!
                opctr=opctr+1
                listofM(:,:,opctr)=sl%quartetop(op)%sotr
            endif
        enddo oploop
        ! Seems reasonable, perhaps. Build the thingy
        BMS=0.0_flyt
        S=0.0_flyt
        U=0.0_flyt
        V=0.0_flyt
        coeffM=0.0_flyt
        zeroSV=0
        fcctr=0
        do i=1,opctr
            BMS=BMS+listofM(:,:,i)-I81
        enddo
        !call lo_dgesvd(BMS,S,U,V,info=lo_status)
        ! Something weird with MKL made me switch SVD routine here. Don't understand, the glitch seems completely random.
        call lo_dgesdd(BMS,S,U,V,jobz='A',info=lo_status)
        l=0
        do i=1,dm
            if ( abs(S(i)) .lt. lo_tol ) then
                l=l+1
                zeroSV(l)=i
            endif
        enddo
        if ( l .eq. 0 .or. l .eq. dm ) then
            ! symmetry did nothing for this shell
            coeffM=I81
            fcctr=dm
        else
            ! symmetry did something, generate that and clean it a bit
            fcctr=l
            coeffM=matmul(U(:,zeroSV(1:fcctr)),V(zeroSV(1:fcctr),:))
            !call lo_real_gram_schmidt(coeffM)
            call lo_make_coeffmatrix_tidy(coeffM)
        endif

        ! store stuff in the shell
        sl%quartetshell(shell)%coeffM=coeffM
        sl%quartetshell(shell)%nfc=fcctr
        lo_allocate(sl%quartetshell(shell)%relfc( fcctr ))
        lo_allocate(sl%quartetshell(shell)%relfcind( fcctr ))
        l=0
        do i=1,dm
            if ( sum(abs(coeffM(:,i))) .gt. lo_sqtol ) then
                l=l+1
                sl%quartetshell(shell)%relfc(l)=i
            else
                coeffM(:,i)=0.0_flyt
            endif
        enddo

        ! small sanity check
        if ( l .ne. fcctr ) then
            write(*,*) 'ERROR:'
            write(*,*) '    my Gram-Schmidt or something failed, I got',l,'ifcs when I expected',fcctr
        endif

        fcctr_global=fcctr_global+fcctr

#ifdef AGRESSIVE_SANITY
! Do some more, agressive sanity tests
sanity: block
    real(flyt), dimension(81,81) :: m0,m1
    real(flyt) :: f0

    ! Check if the SVD was ok, with orthogonal vectors and all that.
    m0=transpose(V)
    f0=0.0_flyt
    do i=1,81
    do j=i+1,81
        f0=f0+abs(dot_product(U(:,i),U(:,j)))+abs(dot_product(m0(:,i),m0(:,j)))
    enddo
    enddo

    ! Check that the coefficient matrix is truly invariant to the operations it should
    ! be invariant under
    m0=sl%quartetshell(shell)%coeffM
    do i=1,opctr
        call lo_gemm(listofM(:,:,i),m0,m1)
        f0=f0+sum(abs(m0-m1))
    enddo
    if ( f0 .gt. lo_sqtol ) then
        call lo_stop_gracefully(['Inconsitent symmetry for quartet shell '//tochar(shell)],lo_exitcode_symmetry,__FILE__,__LINE__)
    endif
end block sanity
#endif

    enddo shloop

    ! count forceconstants and index them
    l=0
    do i=1,sl%nquartetshells
        ! count the others
        do j=1,sl%quartetshell(i)%nfc
            l=l+1
            sl%quartetshell(i)%relfcind(j)=l
        enddo
    enddo
    sl%nquartetifc=l
end subroutine

!> Create the nullspace for each irreducible magnetic pair interaction
module subroutine nullspace_magpairshells(sl,sh,ss)
    type(lo_symlist), intent(inout) :: sl
    type(lo_symtabhelper), intent(inout) :: sh
    type(lo_crystalstructure), intent(in) :: ss

    real(flyt), dimension(:,:,:), allocatable :: listofM,listofM3
    real(flyt), dimension(9,9) :: I9,invarM
    real(flyt), dimension(3,3) :: I3,invarM3
    integer :: shell,op,opctr
    integer :: i,j,ii,jj,kk

    call lo_identitymatrix(I9)
    call lo_identitymatrix(I3)
    lo_allocate(listofM(9,9,sl%npairop))
    lo_allocate(listofM3(3,3,sl%npairop))

    shloop: do shell=1,sl%nmagpairshells
        ! Here we fix the Jij's.
        ! Since this is magnetic, all operations should be allowed!
        ! I'm not sure it's supposed to be all-all though. Hmm. Maybe a slightly
        ! reduced amount would be better. Problem for future Olle.
        invarM=0.0_flyt
        listofM=0.0_flyt
        opctr=0
        opl1: do op=1,sl%npairop
            ! this might be a new, interesting operation
            do i=1,opctr
                if ( sum(abs(listofM(:,:,i)-sl%pairop(op)%sotr)) .lt. lo_tol ) then
                    ! sorry, not new, cycle
                    cycle opl1
                endif
            enddo
            ! yup, it was new!
            opctr=opctr+1
            listofM(:,:,opctr)=sl%pairop(op)%sotr
            invarM=invarM+listofM(:,:,opctr)-I9
        enddo opl1
        call lo_nullspace_coefficient_matrix(invarM,sl%magpairshell(shell)%coeffM_jij,&
                                             sl%magpairshell(shell)%nfc_jij,sl%magpairshell(shell)%relfc_jij)
        if ( sl%magpairshell(shell)%nfc_jij .gt. 0 ) then
            lo_allocate(sl%magpairshell(shell)%relfcind_jij( sl%magpairshell(shell)%nfc_jij ))
            sl%magpairshell(shell)%relfcind_jij=0
        endif

        ! Now create the spin-lattice.
        listofM=0.0_flyt
        opctr=0
        invarM=0.0_flyt
        opl2: do op=1,sl%npairop
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sl%magpairshell(shell)%protpair,&
                 sl%magpairshell(shell)%protpair,sl%pairop(op)) ) then
                if ( sl%pairop(op)%permind .eq. 2 ) cycle
                ! this might be a new, interesting operation
                do i=1,opctr
                    if ( sum(abs(listofM(:,:,i)-sl%pairop(op)%msotr)) .lt. lo_tol ) cycle opl2
                enddo
                ! yup, it was new!
                opctr=opctr+1
                listofM(:,:,opctr)=sl%pairop(op)%msotr
                invarM=invarM+listofM(:,:,opctr)-I9
            endif
        enddo opl2
        call lo_nullspace_coefficient_matrix(invarM,sl%magpairshell(shell)%coeffM_tij,&
             sl%magpairshell(shell)%nfc_tij,sl%magpairshell(shell)%relfc_tij)
        if ( sl%magpairshell(shell)%nfc_tij .gt. 0 ) then
            lo_allocate(sl%magpairshell(shell)%relfcind_tij( sl%magpairshell(shell)%nfc_tij ))
            sl%magpairshell(shell)%relfcind_tij=0
        endif

        ! Try it again, but a rank-1 tensor instead of rank-2
        listofM3=0.0_flyt
        opctr=0
        invarM3=0.0_flyt
        opl3: do op=1,sl%npairop
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sl%magpairshell(shell)%protpair,&
                 sl%magpairshell(shell)%protpair,sl%pairop(op)) ) then
                ! this might be a new, interesting operation
                do i=1,opctr
                    if ( sum(abs(listofM3(:,:,i)-sl%pairop(op)%m3)) .lt. lo_tol ) cycle opl3
                enddo
                ! yup, it was new!
                opctr=opctr+1
                listofM3(:,:,opctr)=sl%pairop(op)%m3
                invarM3=invarM3+listofM3(:,:,opctr)-I3
            endif
        enddo opl3
        call lo_nullspace_coefficient_matrix(invarM3,sl%magpairshell(shell)%coeffM_qij,&
             sl%magpairshell(shell)%nfc_qij,sl%magpairshell(shell)%relfc_qij)
        if ( sl%magpairshell(shell)%nfc_qij .gt. 0 ) then
            lo_allocate(sl%magpairshell(shell)%relfcind_qij( sl%magpairshell(shell)%nfc_qij ))
            sl%magpairshell(shell)%relfcind_qij=0
        endif

    enddo shloop

    ! count forceconstants and index them
    ii=0
    jj=0
    kk=0
    do i=1,sl%nmagpairshells
        ! count the others
        do j=1,sl%magpairshell(i)%nfc_jij
            ii=ii+1
            sl%magpairshell(i)%relfcind_jij(j)=ii
        enddo
        do j=1,sl%magpairshell(i)%nfc_tij
            jj=jj+1
            sl%magpairshell(i)%relfcind_tij(j)=jj
        enddo
        do j=1,sl%magpairshell(i)%nfc_qij
            kk=kk+1
            sl%magpairshell(i)%relfcind_qij(j)=kk
        enddo
    enddo
    sl%nmagpairtheta_jij=ii
    sl%nmagpairtheta_tij=jj
    sl%nmagpairtheta_qij=kk
end subroutine

! ! Create the nullspace
! module subroutine nullspace_epspairshells(sl,sh,ss)
!     type(lo_symlist), intent(inout) :: sl
!     type(lo_symtabhelper), intent(inout) :: sh
!     type(lo_crystalstructure), intent(in) :: ss
!
!     integer, dimension(:), allocatable :: opind
!     integer, dimension(9) :: zeroSV
!     integer :: shell,op,opctr,fcctr,fcctr_global
!     integer :: i,j,l
!     real(flyt), dimension(:,:,:), allocatable :: listofM
!     real(flyt), dimension(9,9) :: I9,BMS,U,V,coeffM,invarM
!     real(flyt), dimension(9) :: S
!
!     !call lo_identitymatrix(I9)
!
!     allocate(listofM(9,9,sl%npairop))
!     allocate(opind(sl%npairop))
!
!     fcctr_global=0
!     shloop: do shell=1,sl%nepsshellpair
!         ! it might be a self-term, in that case it's near trivial to fix
!         ! if ( lo_sqnorm(sl%epspairshell(shell)%protpair%v) .lt. sl%sqtol_cart ) then
!         !     sl%epspairshell(shell)%nfc=0
!         !     sl%epspairshell(shell)%coeffM=0.0_flyt
!         !     cycle shloop
!         ! endif
!         !
!         ! ! Now I assume it's a rather normal shell. Figure out which operations
!         ! ! take the pair back to itself.
!         ! listofM=0.0_flyt
!         ! opctr=0
!         ! opind=0
!         ! oploop: do op=1,sl%npairop
!         !     if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sl%epspairshell(shell)%protpair,sl%epspairshell(shell)%protpair,sl%pairop(op)) ) then
!         !         ! this might be a new, interesting operation
!         !         do i=1,opctr
!         !             if ( sum(abs(listofM(:,:,i)-sl%pairop(op)%sotr)) .lt. lo_tol ) then
!         !                 ! sorry, not new, cycle
!         !                 cycle oploop
!         !             endif
!         !         enddo
!         !         ! yup, it was new!
!         !         opctr=opctr+1
!         !         opind(opctr)=op
!         !         listofM(:,:,opctr)=sl%pairop(op)%sotr
!         !     endif
!         ! enddo oploop
!
!     !     ! Seems reasonable, perhaps. Build the thingy
!     !     BMS=0.0_flyt
!     !     S=0.0_flyt
!     !     U=0.0_flyt
!     !     V=0.0_flyt
!     !     coeffM=0.0_flyt
!     !     zeroSV=0
!     !     fcctr=0
!     !     do i=1,opctr
!     !         BMS=BMS+listofM(:,:,i)-I9
!     !     enddo
!     !     invarM=BMS
!     !     !call lo_dgesvd(BMS,S,U,V)
!     !     call lo_dgesdd(BMS,S,U,V,jobz='A',info=lo_status)
!     !     l=0
!     !     do i=1,9
!     !         if ( abs(S(i)) .lt. lo_tol ) then
!     !             l=l+1
!     !             zeroSV(l)=i
!     !         endif
!     !     enddo
!     !     if ( l .eq. 0 .or. l .eq. 9 ) then
!     !         ! symmetry did nothing for this shell
!     !         coeffM=I9
!     !         fcctr=9
!     !     else
!     !         ! symmetry did something, generate that and clean it a bit
!     !         fcctr=l
!     !         coeffM=matmul(U(:,zeroSV(1:fcctr)),V(zeroSV(1:fcctr),:))
!     !         !call lo_real_gram_schmidt(coeffM)
!     !         call lo_make_coeffmatrix_tidy(coeffM)
!     !     endif
!     !     ! store stuff in the shell
!     !     sl%pairshell(shell)%invarM=invarM
!     !     sl%pairshell(shell)%coeffM=coeffM
!     !     sl%pairshell(shell)%nfc=fcctr
!     !     lo_allocate(sl%pairshell(shell)%relfc( fcctr ))
!     !     lo_allocate(sl%pairshell(shell)%relfcind( fcctr ))
!     !     l=0
!     !     do i=1,9
!     !         if ( sum(abs(coeffM(:,i))) .gt. lo_sqtol ) then
!     !             l=l+1
!     !             sl%pairshell(shell)%relfc(l)=i
!     !         endif
!     !     enddo
!     !
!     !     ! small sanity check
!     !     if ( l .ne. fcctr ) then
!     !         write(*,*) 'ERROR:'
!     !         write(*,*) '    my Gram-Schmidt or something failed, I got',l,'ifcs when I expected',fcctr
!     !     endif
!     !
!     !     fcctr_global=fcctr_global+fcctr
!     enddo shloop
!
!     ! ! count forceconstants and index them
!     ! l=0
!     ! do i=1,sl%npairshells
!     !     ! count the others
!     !     do j=1,sl%pairshell(i)%nfc
!     !         l=l+1
!     !         sl%pairshell(i)%relfcind(j)=l
!     !     enddo
!     ! enddo
!     ! sl%npairifc=l
!
! end subroutine

end submodule
