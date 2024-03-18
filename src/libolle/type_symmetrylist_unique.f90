#include "precompilerdefinitions"
submodule (type_symmetrylist) type_symmetrylist_unique
use type_distancetable
use type_voronoi_distancetable
use type_symmetryoperation, only: lo_spacegroup_operation
implicit none
contains

!> figure out which the unique tuplets are
module subroutine get_unique_tuplets(sl,sh,uc,ss)
    type(lo_symlist), intent(inout) :: sl
    type(lo_symtabhelper), intent(inout) :: sh
    type(lo_crystalstructure), intent(in) :: uc
    type(lo_crystalstructure), intent(in) :: ss

    class(lo_symlist_tuplet), dimension(:), allocatable :: untup
    real(flyt), dimension(:), allocatable :: shellhash
    real(flyt) :: t0
    integer, dimension(:), allocatable :: di
    integer :: i,j,k,l,ii,op
    integer :: a1,a2,a3

    ! So, how this works: I have previously generated structures that contain, for each atom,
    ! a list of all the tuplets it is involved in, where a tuplet is a pair, triplet or quartet.
    ! Schematically, it looks something like this. I have one list for all the atoms in the unit-
    ! cell, and one for all the atoms in the supercell.
    !
    !                     ┌────────┬────────┬───────┬────────┐
    ! ┌────────────┐  ┌──▶│tuplet 1│tuplet 2│  ...  │tuplet n│
    ! │   atom 1   │──┤   └────────┴────────┴───────┴────────┘
    ! └────────────┘  │   ┌────────┬────────┬───────┬────────┐
    !                 └──▶│tuplet 1│tuplet 2│  ...  │tuplet n│
    !                     └────────┴────────┴───────┴────────┘
    !                     ┌────────┬────────┬───────┬────────┐
    ! ┌────────────┐  ┌──▶│tuplet 1│tuplet 2│  ...  │tuplet n│
    ! │   atom 2   │──┤   └────────┴────────┴───────┴────────┘
    ! └────────────┘  │   ┌────────┬────────┬───────┬────────┐
    !                 └──▶│tuplet 1│tuplet 2│  ...  │tuplet n│
    !                     └────────┴────────┴───────┴────────┘
    !
    ! The task is to find the irreducible tuplets, and store these. I define a shell as all the tuplets
    ! that are related via spacegroup or transposition symmetries. For each of these shells, I store a
    ! prototype tuplet.
    ! This routine finds those prototype tuplets, and also how the tuplets in the above structure are
    ! constructed from the prototypes.

    ! the first order should be reasonably simple
    if ( sl%firstorder ) then
        ! number of shells is the same as the number of unique atoms
        sl%nsingletshells=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) ) sl%nsingletshells=sl%nsingletshells+1
        enddo
        ! build the shells
        lo_allocate(sl%singletshell(sl%nsingletshells))
        lo_allocate(di(sl%nuc))
        di=0
        a1=0
        do l=1,sl%nuc
            if ( sh%is_unitcell_atom_prototype(l) .eqv. .false. ) cycle
            a1=a1+1
            sl%singletshell(a1)%protatom=l
            di(l)=a1
        enddo
        ! match the shells, first unitcell
        do a1=1,sl%nuc
            a2=sh%unique_unitcell_atom(a1)
            sl%uc(a1)%singlet%unique_shell=di(a2)
            ii=0
            do op=1,sl%sym%n
                if ( sh%clist(a2,a1,op) ) then
                    ii=op
                    exit
                endif
            enddo
            if ( ii .eq. 0 ) then
                write(*,*) 'Failed mapping from singlet shells to unitcell'
                stop
            else
                sl%uc(a1)%singlet%operation_from_shell_to_atom=ii
            endif
        enddo
        ! then supercell
        do a1=1,sl%nss
            a2=sh%unique_supercell_atom(a1)    ! unique index in unitcell
            a3=ss%info%index_in_unitcell(a1)   ! actual index in unitcell
            sl%ss(a1)%singlet%unique_shell=di(a2)
            ii=0
            do op=1,sl%sym%n
                if ( sh%clist(a2,a3,op) ) then
                    ii=op
                    exit
                endif
            enddo
            if ( ii .eq. 0 ) then
                write(*,*) 'Failed mapping from singlet shells to supercell'
                stop
            else
                sl%ss(a1)%singlet%operation_from_shell_to_atom=ii
            endif
        enddo
        lo_deallocate(di)
    endif

    ! start with the pairs
    if ( sl%secondorder ) then
        ! Return the unique pairs
        call return_unique_tuplets(sl,sh,sh%allpairs,sh%pair_allhash,ss,untup)
        ! Store these as coordination shells
        select type(untup); type is(lo_symlist_ucpair)
            sl%npairshells=size(untup,1)
            lo_allocate(sl%pairshell( sl%npairshells ))
            do i=1,sl%npairshells
                sl%pairshell(i)%protpair=untup(i)
            enddo
        end select
        lo_deallocate(untup)
        ! This might have gone reasonably well. Now the idea is to match this outwards again.

        ! Now I have to match from the prototypes to the full structure.
        ! I will do this with hashes as well, I think that's slightly faster.
        lo_allocate(shellhash(sl%npairshells))
        do i=1,sl%npairshells
            shellhash(i)=hashtuplet( sl%pairshell(i)%protpair )
        enddo

        ! Match it outwards, first the unitcell
        if ( sl%verbosity .gt. 0 ) then
            t0=walltime()
            call lo_progressbar_init()
        endif
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%npair
            l=l+1
            shl1: do j=1,sl%npairshells
                if ( abs(sh%pair_uchash(l)-shellhash(j)) .gt. lo_tol ) cycle shl1
                ! now try to match it properly: try with all operations to match this shell with the prototype
                do op=1,sl%npairop
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                         sl%pairshell(j)%protpair,sl%uc(a1)%pair(i),sl%pairop(op)) ) then
                        sl%uc(a1)%pair(i)%unique_shell=j
                        sl%uc(a1)%pair(i)%operation_from_shell_to_pair=op
                        exit shl1
                    endif
                enddo
            enddo shl1
            ! I should, hopefully, never make it here
            if ( sl%uc(a1)%pair(i)%unique_shell .eq. 0 ) then
                write(*,*) 'ERROR:'
                write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                write(*,*) ' went wrong.'
                write(*,*) ' The vector I want to match:'
                write(*,*) real(sl%uc(a1)%pair(i)%v)
                write(*,*) ' These are the vectors I tried to match with:'
                do j=1,sl%npairshells
                    write(*,*) j,real(sl%pairshell(j)%protpair%v)
                enddo
                write(*,*) ' Future Olle will write some proper diagnostics for this.'
                stop
            endif
        enddo
        if ( sl%verbosity .gt. 0 ) then
            call lo_progressbar(' ... matching pair shells',a1,uc%na+ss%na,walltime()-t0)
        endif
        enddo
        ! Now match it to the supercell
        k=0
        l=0
        do a1=1,sl%nss
            a2=ss%info%index_in_unitcell(a1) ! for the quick test
            do i=1,sl%ss(a1)%npair
                l=l+1
                ! first we try it the easy way
                j=sl%uc(a2)%pair(i)%unique_shell
                op=sl%uc(a2)%pair(i)%operation_from_shell_to_pair
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                     sl%pairshell(j)%protpair,sl%ss(a1)%pair(i),sl%pairop(op)) ) then
                    sl%ss(a1)%pair(i)%unique_shell=j
                    sl%ss(a1)%pair(i)%operation_from_shell_to_pair=op
                else
                    ! The hard way. Should never happen in theory. k is the counter for the number of misses
                    ! and is reported to stdout for debuging.
                    k=k+1
                    shl2: do j=1,sl%npairshells
                        ! skip if hash does not match
                        if ( abs(sh%pair_sshash(l)-shellhash(j)) .gt. lo_tol ) cycle shl2
                        ! now try to match it properly: try with all operations to match this shell with the prototype
                        do op=1,sl%npairop
                            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                                 sl%pairshell(j)%protpair,sl%ss(a1)%pair(i),sl%pairop(op)) ) then
                                sl%ss(a1)%pair(i)%unique_shell=j
                                sl%ss(a1)%pair(i)%operation_from_shell_to_pair=op
                                exit shl2
                            endif
                        enddo
                    enddo shl2
                endif
                ! I should, hopefully, never make it here
                if ( sl%ss(a1)%pair(i)%unique_shell .eq. 0 ) then
                    write(*,*) 'ERROR:'
                    write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                    write(*,*) ' went wrong.'
                    write(*,*) ' The vector I want to match:'
                    write(*,*) real(sl%ss(a1)%pair(i)%v),&
                    tochar([sl%ss(a1)%pair(i)%i1,sl%ss(a1)%pair(i)%i2,sl%ss(a1)%pair(i)%ui1,sl%ss(a1)%pair(i)%ui2])
                    write(*,*) ' These are the vectors I tried to match with:'
                    do j=1,sl%npairshells
                        write(*,*) j,real(sl%pairshell(j)%protpair%v),&
                        tochar([sl%pairshell(j)%protpair%i1,sl%pairshell(j)%protpair%i2,&
                                sl%pairshell(j)%protpair%ui1,sl%pairshell(j)%protpair%ui2])
                    enddo
                    write(*,*) ' Future Olle will write some proper diagnostics for this.'
                    stop
                endif
            enddo
            if ( sl%verbosity .gt. 0 ) then
                if ( k .eq. 0 ) then
                    call lo_progressbar(' ... matching pair shells',uc%na+a1,uc%na+ss%na,walltime()-t0)
                else
                    call lo_progressbar(' ... matching pair shells, misses: '//tochar(k),uc%na+a1,uc%na+ss%na,walltime()-t0)
                endif
            endif
        enddo

        ! For future reference, it might come handy with all the pairs for each shell
        ! first count everything
        do i=1,sl%npairshells
            sl%pairshell(i)%n=0
        enddo
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%npair
            j=sl%uc(a1)%pair(i)%unique_shell
            sl%pairshell(j)%n=sl%pairshell(j)%n+1
        enddo
        enddo
        ! Make some space
        do i=1,sl%npairshells
            lo_allocate(sl%pairshell(i)%atomind( sl%pairshell(i)%n ))
            lo_allocate(sl%pairshell(i)%pairind( sl%pairshell(i)%n ))
            sl%pairshell(i)%atomind=0
            sl%pairshell(i)%pairind=0
            sl%pairshell(i)%n=0
        enddo
        ! store the pairs in the shells
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%npair
            j=sl%uc(a1)%pair(i)%unique_shell
            sl%pairshell(j)%n=sl%pairshell(j)%n+1
            sl%pairshell(j)%atomind( sl%pairshell(j)%n )=a1
            sl%pairshell(j)%pairind( sl%pairshell(j)%n )=i
        enddo
        enddo
        ! some cleanup
        lo_deallocate(shellhash)
    endif ! done with pairs

    ! And perhaps the third order
    if ( sl%thirdorder ) then
        ! Get the unique
        call return_unique_tuplets(sl,sh,sh%alltriplets,sh%triplet_allhash,ss,untup)
        ! Store these as coordination shells
        select type(untup); type is(lo_symlist_uctriplet)
            sl%ntripletshells=size(untup,1)
            lo_allocate(sl%tripletshell( sl%ntripletshells ))
            do i=1,sl%ntripletshells
                sl%tripletshell(i)%prottriplet=untup(i)
            enddo
        end select
        lo_deallocate(untup)
        ! Now match the prototypes to the full structure. I will do this by hashing too.
        lo_allocate(shellhash(sl%ntripletshells))
        do i=1,sl%ntripletshells
            shellhash(i)=hashtuplet( sl%tripletshell(i)%prottriplet )
        enddo

        if ( sl%verbosity .gt. 0 ) then
            t0=walltime()
            call lo_progressbar_init()
        endif

        ! Now I can try to match this
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%ntriplet
            l=l+1
            call locate_triplet_in_supercell(sl%uc(a1)%triplet(i),sl,ss,&
                 sl%uc(a1)%triplet(i)%ssa,sl%uc(a1)%triplet(i)%ssi)
            shl3: do j=1,sl%ntripletshells
                ! skip if hash does not match
                !if ( abs(sh%triplet_uchash(l)-shellhash(j)) .gt. lo_tol ) cycle shl3
                do op=1,sl%ntripletop
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                        sl%tripletshell(j)%prottriplet,sl%uc(a1)%triplet(i),sl%tripletop(op)) ) then
                        sl%uc(a1)%triplet(i)%unique_shell=j
                        sl%uc(a1)%triplet(i)%operation_from_shell_to_triplet=op
                        exit shl3
                    endif
                enddo
            enddo shl3
            ! angry sanity check
            if ( sl%uc(a1)%triplet(i)%unique_shell .eq. 0 ) then
                write(*,*) 'ERROR:'
                write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                write(*,*) ' went wrong.'
                write(*,*) ' I tried to match triplet'
                write(*,*) '       uci: ',tochar([sl%uc(a1)%triplet(i)%i1,&
                                          sl%uc(a1)%triplet(i)%i2,&
                                          sl%uc(a1)%triplet(i)%i3])
                write(*,*) '       uni: ',tochar([sl%uc(a1)%triplet(i)%ui1,&
                                          sl%uc(a1)%triplet(i)%ui2,&
                                          sl%uc(a1)%triplet(i)%ui3])
                write(*,*) '     ssind: ',tochar([sl%uc(a1)%triplet(i)%ssa,&
                                          sl%uc(a1)%triplet(i)%ssi])
                write(*,*) 'ucv1: ',sl%uc(a1)%triplet(i)%v1
                write(*,*) 'ucv2: ',sl%uc(a1)%triplet(i)%v2
                write(*,*) 'ucv3: ',sl%uc(a1)%triplet(i)%v3
                write(*,*) 'And the shell:'
                do j=1,sl%ntripletshells
                    if ( abs(sh%triplet_uchash(l)-shellhash(j)) .gt. lo_tol ) cycle
                    write(*,*) '       uci: ',tochar([sl%tripletshell(j)%prottriplet%i1,&
                                sl%tripletshell(j)%prottriplet%i2,&
                                sl%tripletshell(j)%prottriplet%i3])
                    write(*,*) '       uni: ',tochar([sl%tripletshell(j)%prottriplet%ui1,&
                                sl%tripletshell(j)%prottriplet%ui2,&
                                sl%tripletshell(j)%prottriplet%ui3])
                    write(*,*) tochar([sl%tripletshell(j)%prottriplet%ssa,&
                                sl%tripletshell(j)%prottriplet%ssi])
                    write(*,*) 'ucv1: ',sl%tripletshell(j)%prottriplet%v1
                    write(*,*) 'ucv2: ',sl%tripletshell(j)%prottriplet%v2
                    write(*,*) 'ucv3: ',sl%tripletshell(j)%prottriplet%v3
                enddo
                stop
            endif
        enddo
        if ( sl%verbosity .gt. 0 ) then
            call lo_progressbar(' ... matching triplet shells',a1,uc%na+ss%na,walltime()-t0)
        endif
        enddo

        ! and match the supercell
        k=0
        l=0
        do a1=1,sl%nss
            a2=ss%info%index_in_unitcell(a1)
            do i=1,sl%ss(a1)%ntriplet
                l=l+1
                ! first easy way
                j=sl%uc(a2)%triplet(i)%unique_shell
                op=sl%uc(a2)%triplet(i)%operation_from_shell_to_triplet
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                    sl%tripletshell(j)%prottriplet,sl%ss(a1)%triplet(i),sl%tripletop(op)) ) then
                    sl%ss(a1)%triplet(i)%unique_shell=j
                    sl%ss(a1)%triplet(i)%operation_from_shell_to_triplet=op
                else
                    ! nope, the hard way
                    k=k+1
                    shl4: do j=1,sl%ntripletshells
                        ! skip if hash does not match
                        if ( abs(sh%triplet_sshash(l)-shellhash(j)) .gt. lo_tol ) cycle shl4
                        do op=1,sl%ntripletop
                            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                                sl%tripletshell(j)%prottriplet,sl%ss(a1)%triplet(i),sl%tripletop(op)) ) then
                                sl%ss(a1)%triplet(i)%unique_shell=j
                                sl%ss(a1)%triplet(i)%operation_from_shell_to_triplet=op
                                exit shl4
                            endif
                        enddo
                    enddo shl4
                endif
                ! angry sanity check
                if ( sl%ss(a1)%triplet(i)%unique_shell .eq. 0 ) then
                    write(*,*) 'ERROR:'
                    write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                    write(*,*) ' went wrong.'
                    write(*,*) ' In the supercell, I tried to match triplet'
                    write(*,*) '       uci: ',tochar([sl%ss(a1)%triplet(i)%i1,&
                                              sl%ss(a1)%triplet(i)%i2,&
                                              sl%ss(a1)%triplet(i)%i3])
                    write(*,*) '       uni: ',tochar([sl%ss(a1)%triplet(i)%ui1,&
                                              sl%ss(a1)%triplet(i)%ui2,&
                                              sl%ss(a1)%triplet(i)%ui3])
                    write(*,*) 'ucv1: ',sl%ss(a1)%triplet(i)%v1
                    write(*,*) 'ucv2: ',sl%ss(a1)%triplet(i)%v2
                    write(*,*) 'ucv3: ',sl%ss(a1)%triplet(i)%v3
                    write(*,*) ' ssi: ',tochar([a1,i])
                    write(*,*) 'And the shell:'
                    do j=1,sl%ntripletshells
                        if ( abs(sh%triplet_sshash(l)-shellhash(j)) .gt. lo_tol ) cycle
                        write(*,*) '       uci: ',tochar([sl%tripletshell(j)%prottriplet%i1,&
                                    sl%tripletshell(j)%prottriplet%i2,&
                                    sl%tripletshell(j)%prottriplet%i3])
                        write(*,*) '       nci: ',tochar([sl%tripletshell(j)%prottriplet%ui1,&
                                    sl%tripletshell(j)%prottriplet%ui2,&
                                    sl%tripletshell(j)%prottriplet%ui3])
                        write(*,*) '       ssi: ',tochar([sl%tripletshell(j)%prottriplet%ssa,&
                                    sl%tripletshell(j)%prottriplet%ssi])
                        write(*,*) 'ucv1: ',sl%tripletshell(j)%prottriplet%v1
                        write(*,*) 'ucv2: ',sl%tripletshell(j)%prottriplet%v2
                        write(*,*) 'ucv3: ',sl%tripletshell(j)%prottriplet%v3
                        a2=sl%tripletshell(j)%prottriplet%ssa
                        k=sl%tripletshell(j)%prottriplet%ssi
                        write(*,*) 'ssv1: ',sl%ss(a2)%triplet(k)%v1
                        write(*,*) 'ssv2: ',sl%ss(a2)%triplet(k)%v2
                        write(*,*) 'ssv3: ',sl%ss(a2)%triplet(k)%v3
                        write(*,*) '       uci: ',tochar([sl%ss(a2)%triplet(k)%i1,&
                                              sl%ss(a2)%triplet(k)%i2,&
                                              sl%ss(a2)%triplet(k)%i3])
                        write(*,*) '       uni: ',tochar([sl%ss(a2)%triplet(k)%ui1,&
                                              sl%ss(a2)%triplet(k)%ui2,&
                                              sl%ss(a2)%triplet(k)%ui3])

                    enddo
                    stop
                endif
            enddo
            if ( sl%verbosity .gt. 0 ) then
                if ( k .gt. 0 ) then
                    call lo_progressbar(' ... matching triplet shells '//tochar(k),a1+uc%na,uc%na+ss%na,walltime()-t0)
                else
                    call lo_progressbar(' ... matching triplet shells',a1+uc%na,uc%na+ss%na,walltime()-t0)
                endif
            endif
        enddo

        ! and a little cleanup
        lo_deallocate(shellhash)
    endif

    ! and finally fourth order
    if ( sl%fourthorder ) then
        ! Get the unique
        call return_unique_tuplets(sl,sh,sh%allquartets,sh%quartet_allhash,ss,untup)
        ! Store these as coordination shells
        select type(untup); type is(lo_symlist_ucquartet)
            sl%nquartetshells=size(untup,1)
            lo_allocate(sl%quartetshell( sl%nquartetshells ))
            do i=1,sl%nquartetshells
                sl%quartetshell(i)%protquartet=untup(i)
            enddo
        end select
        lo_deallocate(untup)
        ! Now match the prototypes to the full structure. I will do this by hashing too.
        lo_allocate(shellhash(sl%nquartetshells))
        do i=1,sl%nquartetshells
            shellhash(i)=hashtuplet( sl%quartetshell(i)%protquartet )
        enddo

        if ( sl%verbosity .gt. 0 ) then
            t0=walltime()
            call lo_progressbar_init()
        endif

        ! Now I can try to match this
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nquartet
            l=l+1
            call locate_quartet_in_supercell(sl%uc(a1)%quartet(i),sl,ss,&
                 sl%uc(a1)%quartet(i)%ssa,sl%uc(a1)%quartet(i)%ssi)
            shl5: do j=1,sl%nquartetshells
                ! skip if hash does not match
                if ( abs(sh%quartet_uchash(l)-shellhash(j)) .gt. lo_tol ) cycle shl5
                do op=1,sl%nquartetop
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                        sl%quartetshell(j)%protquartet,sl%uc(a1)%quartet(i),sl%quartetop(op)) ) then
                        sl%uc(a1)%quartet(i)%unique_shell=j
                        sl%uc(a1)%quartet(i)%operation_from_shell_to_quartet=op
                        exit shl5
                    endif
                enddo
            enddo shl5
            ! angry sanity check
            if ( sl%uc(a1)%quartet(i)%unique_shell .eq. 0 ) then
                write(*,*) 'ERROR:'
                write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                write(*,*) ' went wrong.'
                write(*,*) ' I tried to match quartet'
                stop
            endif
        enddo
        if ( sl%verbosity .gt. 0 ) then
            call lo_progressbar(' ... matching quartet shells',a1,uc%na+ss%na,walltime()-t0)
        endif
        enddo
        ! and the supercell
        k=0
        l=0
        do a1=1,sl%nss
            a2=ss%info%index_in_unitcell(a1)
            do i=1,sl%ss(a1)%nquartet
                l=l+1
                ! first easy way
                j=sl%uc(a2)%quartet(i)%unique_shell
                op=sl%uc(a2)%quartet(i)%operation_from_shell_to_quartet
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                    sl%quartetshell(j)%protquartet,sl%ss(a1)%quartet(i),sl%quartetop(op)) ) then
                    sl%ss(a1)%quartet(i)%unique_shell=j
                    sl%ss(a1)%quartet(i)%operation_from_shell_to_quartet=op
                else
                    ! nope, the hard way
                    k=k+1
                    shl6: do j=1,sl%nquartetshells
                        ! skip if hash does not match
                        if ( abs(sh%quartet_sshash(l)-shellhash(j)) .gt. lo_tol ) cycle shl6
                        do op=1,sl%nquartetop
                            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                                sl%quartetshell(j)%protquartet,sl%ss(a1)%quartet(i),sl%quartetop(op)) ) then
                                sl%ss(a1)%quartet(i)%unique_shell=j
                                sl%ss(a1)%quartet(i)%operation_from_shell_to_quartet=op
                                exit shl6
                            endif
                        enddo
                    enddo shl6
                endif
                ! angry sanity check
                if ( sl%ss(a1)%quartet(i)%unique_shell .eq. 0 ) then
                    write(*,*) 'ERROR:'
                    write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                    write(*,*) ' went wrong with the supercell quartets.'
                    stop
                endif
            enddo
            if ( sl%verbosity .gt. 0 ) then
                if ( k .gt. 0 ) then
                    call lo_progressbar(' ... matching quartet shells '//tochar(k),a1+uc%na,uc%na+ss%na,walltime()-t0)
                else
                    call lo_progressbar(' ... matching quartet shells',a1+uc%na,uc%na+ss%na,walltime()-t0)
                endif
            endif
        enddo

        ! and a little cleanup
        lo_deallocate(shellhash)
    endif

    ! Magnetic pairs
    if ( sl%magnetic_pair_interactions ) then
        ! Return the unique magpairs
        call return_unique_tuplets(sl,sh,sh%allmagpairs,sh%magpair_allhash,ss,untup)
        ! Store these as coordination shells
        select type(untup); type is(lo_symlist_ucmagpair)
            sl%nmagpairshells=size(untup,1)
            lo_allocate(sl%magpairshell( sl%nmagpairshells ))
            do i=1,sl%nmagpairshells
                sl%magpairshell(i)%protpair=untup(i)
            enddo
        end select
        lo_deallocate(untup)
        ! Now I have to match from the prototypes to the full structure.
        ! I will do this with hashes as well, I think that's slightly faster.
        lo_allocate(shellhash(sl%nmagpairshells))
        do i=1,sl%nmagpairshells
            shellhash(i)=hashtuplet( sl%magpairshell(i)%protpair )
        enddo

        ! Match it outwards, first the unitcell
        if ( sl%verbosity .gt. 0 ) then
            t0=walltime()
            call lo_progressbar_init()
        endif
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nmagpair
            l=l+1
            ! First match the normal Jij
            mshl1: do j=1,sl%nmagpairshells
                if ( abs(sh%magpair_uchash(l)-shellhash(j)) .gt. lo_tol ) cycle mshl1
                ! now try to match it properly: try with all operations to match this shell with the prototype
                do op=1,sl%npairop
                    if ( sl%pairop(op)%permind .eq. 2 ) then
                        cycle
                    endif
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                         sl%magpairshell(j)%protpair,sl%uc(a1)%magpair(i),sl%pairop(op)) ) then
                        sl%uc(a1)%magpair(i)%unique_shell=j
                        sl%uc(a1)%magpair(i)%operation_from_shell_to_jij=op
                        sl%uc(a1)%magpair(i)%operation_from_shell_to_bird=op
                        sl%uc(a1)%magpair(i)%operation_from_shell_to_dude=-1
                        exit mshl1
                    endif
                enddo
            enddo mshl1
            ! I should, hopefully, never make it here
            if ( sl%uc(a1)%magpair(i)%unique_shell .eq. 0 ) then
                call lo_stop_gracefully(['Could not locate shell. Futuore olle will fix this.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        if ( sl%verbosity .gt. 0 ) then
            call lo_progressbar(' ... matching magpair shells',a1,uc%na+ss%na,walltime()-t0)
        endif
        enddo
        ! Now match it to the supercell
        k=0
        l=0
        do a1=1,sl%nss
            a2=ss%info%index_in_unitcell(a1) ! for the quick test
            do i=1,sl%ss(a1)%nmagpair
                l=l+1
                ! first we try it the easy way
                j=sl%uc(a2)%magpair(i)%unique_shell
                op=sl%uc(a2)%magpair(i)%operation_from_shell_to_jij
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                     sl%magpairshell(j)%protpair,sl%ss(a1)%magpair(i),sl%pairop(op)) ) then
                    sl%ss(a1)%magpair(i)%unique_shell=j
                    sl%ss(a1)%magpair(i)%operation_from_shell_to_jij=op
                    sl%ss(a1)%magpair(i)%operation_from_shell_to_bird=op
                    sl%ss(a1)%magpair(i)%operation_from_shell_to_dude=-1
                else
                    ! The hard way. Should never happen in theory. k is the counter for the number of misses
                    ! and is reported to stdout for debuging.
                    k=k+1
                    mshl2: do j=1,sl%nmagpairshells
                        ! skip if hash does not match
                        if ( abs(sh%magpair_sshash(l)-shellhash(j)) .gt. lo_tol ) cycle mshl2
                        ! now try to match it properly: try with all operations to match this shell with the prototype
                        do op=1,sl%npairop
                            if ( sl%pairop(op)%permind .eq. 2 ) then
                                cycle
                            endif
                            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                                 sl%magpairshell(j)%protpair,sl%ss(a1)%magpair(i),sl%pairop(op)) ) then
                                sl%ss(a1)%magpair(i)%unique_shell=j
                                sl%ss(a1)%magpair(i)%operation_from_shell_to_jij=op
                                sl%ss(a1)%magpair(i)%operation_from_shell_to_bird=op
                                sl%ss(a1)%magpair(i)%operation_from_shell_to_dude=-1
                                exit mshl2
                            endif
                        enddo
                    enddo mshl2
                endif
                ! I should, hopefully, never make it here
                if ( sl%ss(a1)%magpair(i)%unique_shell .eq. 0 ) then
                    call lo_stop_gracefully(['Could not locate shell. Futuore olle will fix this.'],lo_exitcode_symmetry,__FILE__,__LINE__)
                endif
            enddo
            if ( sl%verbosity .gt. 0 ) then
                if ( k .eq. 0 ) then
                    call lo_progressbar(' ... matching magpair shells',uc%na+a1,uc%na+ss%na,walltime()-t0)
                else
                    call lo_progressbar(' ... matching magpair miss '//tochar(k),uc%na+a1,uc%na+ss%na,walltime()-t0)
                endif
            endif
        enddo

        ! For future reference, it might come handy with all the magpairs for each shell
        ! first count everything
        do i=1,sl%nmagpairshells
            sl%magpairshell(i)%n=0
        enddo
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nmagpair
            j=sl%uc(a1)%magpair(i)%unique_shell
            sl%magpairshell(j)%n=sl%magpairshell(j)%n+1
        enddo
        enddo
        ! Make some space
        do i=1,sl%nmagpairshells
            lo_allocate(sl%magpairshell(i)%atomind( sl%magpairshell(i)%n ))
            lo_allocate(sl%magpairshell(i)%pairind( sl%magpairshell(i)%n ))
            sl%magpairshell(i)%atomind=0
            sl%magpairshell(i)%pairind=0
            sl%magpairshell(i)%n=0
        enddo
        ! store the magpairs in the shells
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nmagpair
            j=sl%uc(a1)%magpair(i)%unique_shell
            sl%magpairshell(j)%n=sl%magpairshell(j)%n+1
            sl%magpairshell(j)%atomind( sl%magpairshell(j)%n )=a1
            sl%magpairshell(j)%pairind( sl%magpairshell(j)%n )=i
        enddo
        enddo
        ! some cleanup
        lo_deallocate(shellhash)

        ! No matter what, I will need the unique magnetic on-site things, if not for anything
        ! else, it's useful for averaging.
        sl%nmagsingletshells=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) ) sl%nmagsingletshells=sl%nmagsingletshells+1
        enddo
        ! build the shells
        lo_allocate(sl%magsingletshell(sl%nmagsingletshells))
        lo_allocate(di(sl%nuc))
        di=0
        a1=0
        do l=1,sl%nuc
            if ( sh%is_unitcell_atom_prototype(l) .eqv. .false. ) cycle
            a1=a1+1
            sl%magsingletshell(a1)%protatom=l
            di(l)=a1
        enddo
        ! match the shells, first unitcell
        do a1=1,sl%nuc
            a2=sh%unique_unitcell_atom(a1)
            sl%uc(a1)%magsinglet%unique_shell=di(a2)
            ii=0
            do op=1,sl%sym%n
                if ( sh%clist(a2,a1,op) ) then
                    ii=op
                    exit
                endif
            enddo
            if ( ii .eq. 0 ) then
                call lo_stop_gracefully(['Could not locate shell. Futuore olle will fix this.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            else
                sl%uc(a1)%magsinglet%operation_from_shell_to_atom=ii
            endif
        enddo
        ! then supercell
        do a1=1,sl%nss
            a2=sh%unique_supercell_atom(a1)    ! unique index in unitcell
            a3=ss%info%index_in_unitcell(a1)   ! actual index in unitcell
            sl%ss(a1)%magsinglet%unique_shell=di(a2)
            ii=0
            do op=1,sl%sym%n
                if ( sh%clist(a2,a3,op) ) then
                    ii=op
                    exit
                endif
            enddo
            if ( ii .eq. 0 ) then
                call lo_stop_gracefully(['Could not locate shell. Futuore olle will fix this.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            else
                sl%ss(a1)%magsinglet%operation_from_shell_to_atom=ii
            endif
        enddo
        lo_deallocate(di)
    endif ! done with pairs

    ! start with the pairs
    if ( sl%dielectric_pair ) then
        ! Return the unique pairs (epspairs and zpairs are equivalent)
        call return_unique_tuplets(sl,sh,sh%allepspairs,sh%epspair_allhash,ss,untup)
        ! Store the unique tuplets as coordination shells
        select type(untup); type is(lo_symlist_ucpair)
            sl%nepsshellpair=size(untup,1)
            sl%nZshellpair=size(untup,1)
            allocate(sl%epspairshell( sl%nepsshellpair ))
            allocate(sl%Zpairshell( sl%nZshellpair ))
            do i=1,sl%nZshellpair
                sl%epspairshell(i)%protpair=untup(i)
                sl%Zpairshell(i)%protpair=untup(i)
            enddo
        end select
        deallocate(untup)
        ! This might have gone reasonably well. Now the idea is to match this outwards again.
        ! I will do this with hashes as well, I think that's slightly faster.
        allocate(shellhash(sl%nZshellpair))
        do i=1,sl%nZshellpair
            shellhash(i)=hashtuplet( sl%Zpairshell(i)%protpair )
        enddo

        ! Match it outwards, first the unitcell
        if ( sl%verbosity .gt. 0 ) then
            t0=walltime()
            call lo_progressbar_init()
        endif
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nZpair
            l=l+1
            dshl1: do j=1,sl%nZshellpair
                if ( abs(sh%Zpair_uchash(l)-shellhash(j)) .gt. lo_tol ) cycle dshl1
                ! now try to match it properly: try with all operations to match this shell with the prototype
                do op=1,sl%npairop
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                         sl%Zpairshell(j)%protpair,sl%uc(a1)%Zpair(i),sl%pairop(op)) ) then
                        sl%uc(a1)%Zpair(i)%unique_shell=j
                        sl%uc(a1)%Zpair(i)%operation_from_shell_to_pair=op
                        sl%uc(a1)%epspair(i)%unique_shell=j
                        sl%uc(a1)%epspair(i)%operation_from_shell_to_pair=op
                        exit dshl1
                    endif
                enddo
            enddo dshl1
            ! I should, hopefully, never make it here
            if ( sl%uc(a1)%Zpair(i)%unique_shell .eq. 0 ) then
                write(*,*) 'ERROR:'
                write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                stop
            endif
        enddo
        if ( sl%verbosity .gt. 0 ) then
            call lo_progressbar(' ... matching pair shells',a1,uc%na+ss%na,walltime()-t0)
        endif
        enddo
        ! Now match it to the supercell
        k=0
        l=0
        do a1=1,sl%nss
            a2=ss%info%index_in_unitcell(a1) ! for the quick test
            do i=1,sl%ss(a1)%nZpair
                l=l+1
                ! first we try it the easy way
                j=sl%uc(a2)%Zpair(i)%unique_shell
                op=sl%uc(a2)%Zpair(i)%operation_from_shell_to_pair
                if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                     sl%Zpairshell(j)%protpair,sl%ss(a1)%Zpair(i),sl%pairop(op)) ) then
                    sl%ss(a1)%Zpair(i)%unique_shell=j
                    sl%ss(a1)%Zpair(i)%operation_from_shell_to_pair=op
                    sl%ss(a1)%epspair(i)%unique_shell=j
                    sl%ss(a1)%epspair(i)%operation_from_shell_to_pair=op
                else
                    ! The hard way. Should never happen in theory. k is the counter for the number of misses
                    ! and is reported to stdout for debuging.
                    k=k+1
                    dshl2: do j=1,sl%nZshellpair
                        ! skip if hash does not match
                        if ( abs(sh%Zpair_sshash(l)-shellhash(j)) .gt. lo_tol ) cycle dshl2
                        ! now try to match it properly: try with all operations to match this shell with the prototype
                        do op=1,sl%npairop
                            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                                 sl%Zpairshell(j)%protpair,sl%ss(a1)%Zpair(i),sl%pairop(op)) ) then
                                sl%ss(a1)%Zpair(i)%unique_shell=j
                                sl%ss(a1)%Zpair(i)%operation_from_shell_to_pair=op
                                sl%ss(a1)%epspair(i)%unique_shell=j
                                sl%ss(a1)%epspair(i)%operation_from_shell_to_pair=op
                                exit dshl2
                            endif
                        enddo
                    enddo dshl2
                endif
                ! I should, hopefully, never make it here
                if ( sl%ss(a1)%Zpair(i)%unique_shell .eq. 0 ) then
                    write(*,*) 'ERROR:'
                    write(*,*) ' the unique shells where perhaps not as unique as I thought. Or something else'
                    stop
                endif
            enddo
            if ( sl%verbosity .gt. 0 ) then
                if ( k .eq. 0 ) then
                    call lo_progressbar(' ... matching pair shells',uc%na+a1,uc%na+ss%na,walltime()-t0)
                else
                    call lo_progressbar(' ... matching pair shells, misses: '//tochar(k),uc%na+a1,uc%na+ss%na,walltime()-t0)
                endif
            endif
        enddo
        ! some cleanup
        lo_deallocate(shellhash)
    endif ! done with pairs
end subroutine

!> build all the empty tuplets
module subroutine setuptuplets(sl,sh,uc,ss)
    !> list of symmetry stuff
    type(lo_symlist), intent(inout) :: sl
    !> helper with more stuff
    type(lo_symtabhelper), intent(inout) :: sh
    !> unit and supercell
    type(lo_crystalstructure), intent(in) :: uc,ss

    type(lo_distancetable) :: dmuc,dmss ! distance tables for magnetic interactions
    real(flyt) :: rcsq
    integer :: a1,i,j,k,l

    ! I will need two more distance tables for magnetic interactions
    if ( sl%magnetic_pair_interactions ) then
    newdt: block
        integer, dimension(:), allocatable :: di1
        integer :: uca

        call dmuc%generate(uc%r,uc%latticevectors,sl%mc2,verbosity=0)
        call dmss%generate(ss%r,ss%latticevectors,sl%mc2,verbosity=0)
        do a1=1,ss%na
            uca=ss%info%index_in_unitcell(a1)
            if ( dmuc%particle(uca)%n .ne. dmss%particle(a1)%n ) then
                write(*,*) 'ERROR:'
                write(*,*) 'atom ',tochar(uca),' in the unitcell and atom ',tochar(a1),'in the supercell should have the same'
                write(*,*) 'number of neighbours within the magnetic cutoff, but they do not.'
                stop
            endif
            lo_allocate(di1(dmuc%particle(uca)%n))
            di1=0
            p1l: do i=1,dmss%particle(a1)%n
                do j=1,dmuc%particle(uca)%n
                    ! same distance
                    if ( abs(dmuc%particle(uca)%d(j)-dmss%particle(a1)%d(i)) .gt. sl%tol_cart ) cycle
                    ! same vectors
                    if ( lo_sqnorm( dmuc%particle(uca)%v(:,j)-dmss%particle(a1)%v(:,i) ) .lt. sl%sqtol_cart ) then
                        di1(j)=i
                        cycle p1l
                    endif
                enddo
                ! if we make it here, there is a discrepancy in the distance tables
                write(*,*) 'ERROR:'
                write(*,*) '    discrepancy in the distance tables between the unit and supercell'
                write(*,*) '    make sure they are generated from consistent cells.'
                stop
            enddo p1l
            ! now I just reorder everything
            dmss%particle(a1)%ind=dmss%particle(a1)%ind(di1)
            dmss%particle(a1)%d=dmss%particle(a1)%d(di1)
            dmss%particle(a1)%v=dmss%particle(a1)%v(:,di1)
            dmss%particle(a1)%lv=dmss%particle(a1)%lv(:,di1)
            if ( allocated(dmss%particle(a1)%weight) ) then
                dmss%particle(a1)%weight=dmss%particle(a1)%weight(di1)
            endif
            lo_deallocate(di1)
        enddo
    end block newdt
    endif

    ! First I generate all the tuplets
    ucatomloop: do a1=1,sl%nuc
        ! singlets
        if ( sl%firstorder ) then
            sl%uc(a1)%singlet%unique_shell=0
        endif

        ! pairs in unitcell
        if ( sl%secondorder ) then
            sl%uc(a1)%npair=sh%dtuc%particle(a1)%n
            lo_allocate(sl%uc(a1)%pair( sl%uc(a1)%npair ))
            do i=1,sl%uc(a1)%npair
                sl%uc(a1)%pair(i)%i1=a1
                sl%uc(a1)%pair(i)%i2=sh%dtuc%particle(a1)%ind(i)
                sl%uc(a1)%pair(i)%ui1=sh%unique_unitcell_atom( sl%uc(a1)%pair(i)%i1 )
                sl%uc(a1)%pair(i)%ui2=sh%unique_unitcell_atom( sl%uc(a1)%pair(i)%i2 )
                sl%uc(a1)%pair(i)%v=sh%dtuc%particle(a1)%v(:,i)
                sl%uc(a1)%pair(i)%unique_shell=0
            enddo
        endif
        ! magnetic pairs in unitcell
        if ( sl%magnetic_pair_interactions ) then
            ! Note that I skip the self-term here
            sl%uc(a1)%nmagpair=dmuc%particle(a1)%n-1
            lo_allocate(sl%uc(a1)%magpair( sl%uc(a1)%nmagpair ))
            do i=1,sl%uc(a1)%nmagpair
                sl%uc(a1)%magpair(i)%i1=a1
                sl%uc(a1)%magpair(i)%i2=dmuc%particle(a1)%ind(i+1)
                sl%uc(a1)%magpair(i)%ui1=sh%unique_unitcell_atom( sl%uc(a1)%magpair(i)%i1 )
                sl%uc(a1)%magpair(i)%ui2=sh%unique_unitcell_atom( sl%uc(a1)%magpair(i)%i2 )
                sl%uc(a1)%magpair(i)%v=dmuc%particle(a1)%v(:,i+1)
                sl%uc(a1)%magpair(i)%unique_shell=0
            enddo
        endif

        ! triplets in unitcell
        if ( sl%thirdorder ) then
            rcsq=sl%rc3**2
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%rc3 ) cycle
                do j=1,sh%dtuc%particle(a1)%n
                    if ( sh%dtuc%particle(a1)%d(j) .gt. sl%rc3 ) cycle
                    if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                    endif
                enddo
            enddo
            sl%uc(a1)%ntriplet=l
            lo_allocate(sl%uc(a1)%triplet( sl%uc(a1)%ntriplet ))
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%rc3 ) cycle
                do j=1,sh%dtuc%particle(a1)%n
                    if ( sh%dtuc%particle(a1)%d(j) .gt. sl%rc3 ) cycle
                    if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                        sl%uc(a1)%triplet(l)%i1=a1
                        sl%uc(a1)%triplet(l)%i2=sh%dtuc%particle(a1)%ind(i)
                        sl%uc(a1)%triplet(l)%i3=sh%dtuc%particle(a1)%ind(j)
                        sl%uc(a1)%triplet(l)%ui1=sh%unique_unitcell_atom( sl%uc(a1)%triplet(l)%i1 )
                        sl%uc(a1)%triplet(l)%ui2=sh%unique_unitcell_atom( sl%uc(a1)%triplet(l)%i2 )
                        sl%uc(a1)%triplet(l)%ui3=sh%unique_unitcell_atom( sl%uc(a1)%triplet(l)%i3 )
                        sl%uc(a1)%triplet(l)%v1=0.0_flyt
                        sl%uc(a1)%triplet(l)%v2=sh%dtuc%particle(a1)%v(:,i)
                        sl%uc(a1)%triplet(l)%v3=sh%dtuc%particle(a1)%v(:,j)
                        sl%uc(a1)%triplet(l)%unique_shell=0
                        sl%uc(a1)%triplet(l)%ssa=0
                        sl%uc(a1)%triplet(l)%ssi=0
                    endif
                enddo
            enddo
        endif
        ! quartets in unitcell
        if ( sl%fourthorder ) then
            rcsq=sl%rc4**2
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%rc4 ) cycle
                do j=1,sh%dtuc%particle(a1)%n
                    if ( sh%dtuc%particle(a1)%d(j) .gt. sl%rc4 ) cycle
                    do k=1,sh%dtuc%particle(a1)%n
                        if ( sh%dtuc%particle(a1)%d(k) .gt. sl%rc4 ) cycle
                        if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,j)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,k)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,j)-sh%dtuc%particle(a1)%v(:,k)) .lt. rcsq ) l=l+1
                    enddo
                enddo
            enddo
            sl%uc(a1)%nquartet=l
            lo_allocate(sl%uc(a1)%quartet( sl%uc(a1)%nquartet ))
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%rc4 ) cycle
                do j=1,sh%dtuc%particle(a1)%n
                    if ( sh%dtuc%particle(a1)%d(j) .gt. sl%rc4 ) cycle
                    do k=1,sh%dtuc%particle(a1)%n
                        if ( sh%dtuc%particle(a1)%d(k) .gt. sl%rc4 ) cycle
                        if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,j)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,k)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,j)-sh%dtuc%particle(a1)%v(:,k)) .lt. rcsq ) then
                            l=l+1
                            sl%uc(a1)%quartet(l)%i1=a1
                            sl%uc(a1)%quartet(l)%i2=sh%dtuc%particle(a1)%ind(i)
                            sl%uc(a1)%quartet(l)%i3=sh%dtuc%particle(a1)%ind(j)
                            sl%uc(a1)%quartet(l)%i4=sh%dtuc%particle(a1)%ind(k)
                            sl%uc(a1)%quartet(l)%ui1=sh%unique_unitcell_atom( sl%uc(a1)%quartet(l)%i1 )
                            sl%uc(a1)%quartet(l)%ui2=sh%unique_unitcell_atom( sl%uc(a1)%quartet(l)%i2 )
                            sl%uc(a1)%quartet(l)%ui3=sh%unique_unitcell_atom( sl%uc(a1)%quartet(l)%i3 )
                            sl%uc(a1)%quartet(l)%ui4=sh%unique_unitcell_atom( sl%uc(a1)%quartet(l)%i4 )
                            sl%uc(a1)%quartet(l)%v1=0.0_flyt
                            sl%uc(a1)%quartet(l)%v2=sh%dtuc%particle(a1)%v(:,i)
                            sl%uc(a1)%quartet(l)%v3=sh%dtuc%particle(a1)%v(:,j)
                            sl%uc(a1)%quartet(l)%v4=sh%dtuc%particle(a1)%v(:,k)
                            sl%uc(a1)%quartet(l)%unique_shell=0
                            sl%uc(a1)%quartet(l)%ssa=0
                            sl%uc(a1)%quartet(l)%ssi=0
                        endif
                    enddo
                enddo
            enddo
        endif

        ! dielectric things
        if ( sl%dielectric_pair ) then
            rcsq=sl%dc2**2
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%dc2 ) cycle
                l=l+1
            enddo
            sl%uc(a1)%nepspair=l
            sl%uc(a1)%nZpair=l
            allocate(sl%uc(a1)%epspair(l))
            allocate(sl%uc(a1)%Zpair(l))
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%dc2 ) cycle
                l=l+1
                sl%uc(a1)%epspair(l)%i1=a1
                sl%uc(a1)%epspair(l)%i2=sh%dtuc%particle(a1)%ind(i)
                sl%uc(a1)%epspair(l)%ui1=sh%unique_unitcell_atom( sl%uc(a1)%pair(i)%i1 )
                sl%uc(a1)%epspair(l)%ui2=sh%unique_unitcell_atom( sl%uc(a1)%pair(i)%i2 )
                sl%uc(a1)%epspair(l)%v=sh%dtuc%particle(a1)%v(:,i)
                sl%uc(a1)%epspair(l)%unique_shell=0
                sl%uc(a1)%Zpair(l)%i1=a1
                sl%uc(a1)%Zpair(l)%i2=sh%dtuc%particle(a1)%ind(i)
                sl%uc(a1)%Zpair(l)%ui1=sh%unique_unitcell_atom( sl%uc(a1)%pair(i)%i1 )
                sl%uc(a1)%Zpair(l)%ui2=sh%unique_unitcell_atom( sl%uc(a1)%pair(i)%i2 )
                sl%uc(a1)%Zpair(l)%v=sh%dtuc%particle(a1)%v(:,i)
                sl%uc(a1)%Zpair(l)%unique_shell=0
            enddo
        endif
        if ( sl%dielectric_triplet ) then
            rcsq=sl%dc3**2
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%dc3 ) cycle
                do j=1,sh%dtuc%particle(a1)%n
                    if ( sh%dtuc%particle(a1)%d(j) .gt. sl%dc3 ) cycle
                    if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                    endif
                enddo
            enddo
            sl%uc(a1)%nZtriplet=l
            allocate(sl%uc(a1)%Ztriplet( sl%uc(a1)%ntriplet ))
            l=0
            do i=1,sh%dtuc%particle(a1)%n
                if ( sh%dtuc%particle(a1)%d(i) .gt. sl%dc3 ) cycle
                do j=1,sh%dtuc%particle(a1)%n
                    if ( sh%dtuc%particle(a1)%d(j) .gt. sl%dc3 ) cycle
                    if ( lo_sqnorm(sh%dtuc%particle(a1)%v(:,i)-sh%dtuc%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                        sl%uc(a1)%Ztriplet(l)%i1=a1
                        sl%uc(a1)%Ztriplet(l)%i2=sh%dtuc%particle(a1)%ind(i)
                        sl%uc(a1)%Ztriplet(l)%i3=sh%dtuc%particle(a1)%ind(j)
                        sl%uc(a1)%Ztriplet(l)%ui1=sh%unique_unitcell_atom( sl%uc(a1)%triplet(l)%i1 )
                        sl%uc(a1)%Ztriplet(l)%ui2=sh%unique_unitcell_atom( sl%uc(a1)%triplet(l)%i2 )
                        sl%uc(a1)%Ztriplet(l)%ui3=sh%unique_unitcell_atom( sl%uc(a1)%triplet(l)%i3 )
                        sl%uc(a1)%Ztriplet(l)%v1=0.0_flyt
                        sl%uc(a1)%Ztriplet(l)%v2=sh%dtuc%particle(a1)%v(:,i)
                        sl%uc(a1)%Ztriplet(l)%v3=sh%dtuc%particle(a1)%v(:,j)
                        sl%uc(a1)%Ztriplet(l)%unique_shell=0
                        sl%uc(a1)%Ztriplet(l)%ssa=0
                        sl%uc(a1)%Ztriplet(l)%ssi=0
                    endif
                enddo
            enddo
        endif
    enddo ucatomloop

    ssatomloop: do a1=1,sl%nss
        ! singlets in supercell
        if ( sl%firstorder ) then
            sl%ss(a1)%singlet%unique_shell=0
        endif
        ! pairs in supercell
        if ( sl%secondorder ) then
            sl%ss(a1)%npair=sh%dtss%particle(a1)%n
            lo_allocate(sl%ss(a1)%pair( sl%ss(a1)%npair ))
            do i=1,sl%ss(a1)%npair
                sl%ss(a1)%pair(i)%j1=a1
                sl%ss(a1)%pair(i)%j2=sh%dtss%particle(a1)%ind(i)
                sl%ss(a1)%pair(i)%i1=ss%info%index_in_unitcell( sl%ss(a1)%pair(i)%j1 )
                sl%ss(a1)%pair(i)%i2=ss%info%index_in_unitcell( sl%ss(a1)%pair(i)%j2 )
                sl%ss(a1)%pair(i)%ui1=sh%unique_unitcell_atom( sl%ss(a1)%pair(i)%i1 )
                sl%ss(a1)%pair(i)%ui2=sh%unique_unitcell_atom( sl%ss(a1)%pair(i)%i2 )
                sl%ss(a1)%pair(i)%v=sh%dtss%particle(a1)%v(:,i)
                sl%ss(a1)%pair(i)%unique_shell=0
            enddo
        endif
        ! magnetic pairs in supercell
        if ( sl%magnetic_pair_interactions ) then
            ! Note that I skip the self-term here as well
            sl%ss(a1)%nmagpair=dmss%particle(a1)%n-1
            lo_allocate(sl%ss(a1)%magpair( sl%ss(a1)%nmagpair ))
            do i=1,sl%ss(a1)%nmagpair
                sl%ss(a1)%magpair(i)%j1=a1
                sl%ss(a1)%magpair(i)%j2=dmss%particle(a1)%ind(i+1)
                sl%ss(a1)%magpair(i)%i1=ss%info%index_in_unitcell( sl%ss(a1)%magpair(i)%j1 )
                sl%ss(a1)%magpair(i)%i2=ss%info%index_in_unitcell( sl%ss(a1)%magpair(i)%j2 )
                sl%ss(a1)%magpair(i)%ui1=sh%unique_unitcell_atom( sl%ss(a1)%magpair(i)%i1 )
                sl%ss(a1)%magpair(i)%ui2=sh%unique_unitcell_atom( sl%ss(a1)%magpair(i)%i2 )
                sl%ss(a1)%magpair(i)%v=dmss%particle(a1)%v(:,i+1)
                sl%ss(a1)%magpair(i)%unique_shell=0
            enddo
        endif
        ! triplets in supercell
        if ( sl%thirdorder ) then
            rcsq=sl%rc3**2
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%rc3 ) cycle
                do j=1,sh%dtss%particle(a1)%n
                    if ( sh%dtss%particle(a1)%d(j) .gt. sl%rc3 ) cycle
                    if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                    endif
                enddo
            enddo
            sl%ss(a1)%ntriplet=l
            lo_allocate(sl%ss(a1)%triplet( sl%ss(a1)%ntriplet ))
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%rc3 ) cycle
                do j=1,sh%dtss%particle(a1)%n
                    if ( sh%dtss%particle(a1)%d(j) .gt. sl%rc3 ) cycle
                    if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                        sl%ss(a1)%triplet(l)%j1=a1
                        sl%ss(a1)%triplet(l)%j2=sh%dtss%particle(a1)%ind(i)
                        sl%ss(a1)%triplet(l)%j3=sh%dtss%particle(a1)%ind(j)
                        sl%ss(a1)%triplet(l)%i1=ss%info%index_in_unitcell( sl%ss(a1)%triplet(i)%j1 )
                        sl%ss(a1)%triplet(l)%i2=ss%info%index_in_unitcell( sl%ss(a1)%triplet(i)%j2 )
                        sl%ss(a1)%triplet(l)%i3=ss%info%index_in_unitcell( sl%ss(a1)%triplet(i)%j3 )
                        sl%ss(a1)%triplet(l)%ui1=sh%unique_supercell_atom( sl%ss(a1)%triplet(l)%j1 )
                        sl%ss(a1)%triplet(l)%ui2=sh%unique_supercell_atom( sl%ss(a1)%triplet(l)%j2 )
                        sl%ss(a1)%triplet(l)%ui3=sh%unique_supercell_atom( sl%ss(a1)%triplet(l)%j3 )
                        sl%ss(a1)%triplet(l)%v1=0.0_flyt
                        sl%ss(a1)%triplet(l)%v2=sh%dtss%particle(a1)%v(:,i)
                        sl%ss(a1)%triplet(l)%v3=sh%dtss%particle(a1)%v(:,j)
                        sl%ss(a1)%triplet(l)%unique_shell=0
                    endif
                enddo
            enddo
        endif
        ! quartets in supercell
        if ( sl%fourthorder ) then
            rcsq=sl%rc4**2
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%rc4 ) cycle
                do j=1,sh%dtss%particle(a1)%n
                    if ( sh%dtss%particle(a1)%d(j) .gt. sl%rc4 ) cycle
                    do k=1,sh%dtss%particle(a1)%n
                        if ( sh%dtss%particle(a1)%d(k) .gt. sl%rc4 ) cycle
                        if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,j)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,k)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,j)-sh%dtss%particle(a1)%v(:,k)) .lt. rcsq ) l=l+1
                    enddo
                enddo
            enddo
            sl%ss(a1)%nquartet=l
            lo_allocate(sl%ss(a1)%quartet( sl%ss(a1)%nquartet ))
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%rc4 ) cycle
                do j=1,sh%dtss%particle(a1)%n
                    if ( sh%dtss%particle(a1)%d(j) .gt. sl%rc4 ) cycle
                    do k=1,sh%dtss%particle(a1)%n
                        if ( sh%dtss%particle(a1)%d(k) .gt. sl%rc4 ) cycle
                        if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,j)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,k)) .gt. rcsq ) cycle
                        if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,j)-sh%dtss%particle(a1)%v(:,k)) .lt. rcsq ) then
                            l=l+1
                            sl%ss(a1)%quartet(l)%j1=a1
                            sl%ss(a1)%quartet(l)%j2=sh%dtss%particle(a1)%ind(i)
                            sl%ss(a1)%quartet(l)%j3=sh%dtss%particle(a1)%ind(j)
                            sl%ss(a1)%quartet(l)%j4=sh%dtss%particle(a1)%ind(k)
                            sl%ss(a1)%quartet(l)%i1=ss%info%index_in_unitcell( sl%ss(a1)%quartet(l)%j1 )
                            sl%ss(a1)%quartet(l)%i2=ss%info%index_in_unitcell( sl%ss(a1)%quartet(l)%j2 )
                            sl%ss(a1)%quartet(l)%i3=ss%info%index_in_unitcell( sl%ss(a1)%quartet(l)%j3 )
                            sl%ss(a1)%quartet(l)%i4=ss%info%index_in_unitcell( sl%ss(a1)%quartet(l)%j4 )
                            sl%ss(a1)%quartet(l)%ui1=sh%unique_unitcell_atom( sl%ss(a1)%quartet(l)%i1 )
                            sl%ss(a1)%quartet(l)%ui2=sh%unique_unitcell_atom( sl%ss(a1)%quartet(l)%i2 )
                            sl%ss(a1)%quartet(l)%ui3=sh%unique_unitcell_atom( sl%ss(a1)%quartet(l)%i3 )
                            sl%ss(a1)%quartet(l)%ui4=sh%unique_unitcell_atom( sl%ss(a1)%quartet(l)%i4 )
                            sl%ss(a1)%quartet(l)%v1=0.0_flyt
                            sl%ss(a1)%quartet(l)%v2=sh%dtss%particle(a1)%v(:,i)
                            sl%ss(a1)%quartet(l)%v3=sh%dtss%particle(a1)%v(:,j)
                            sl%ss(a1)%quartet(l)%v4=sh%dtss%particle(a1)%v(:,k)
                            sl%ss(a1)%quartet(l)%unique_shell=0
                        endif
                    enddo
                enddo
            enddo
        endif
        ! dielectric things
        if ( sl%dielectric_pair ) then
            rcsq=sl%dc2**2
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%dc2 ) cycle
                l=l+1
            enddo
            sl%ss(a1)%nepspair=l
            sl%ss(a1)%nZpair=l
            allocate(sl%ss(a1)%epspair(l))
            allocate(sl%ss(a1)%Zpair(l))
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%dc2 ) cycle
                l=l+1
                sl%ss(a1)%epspair(l)%j1=a1
                sl%ss(a1)%epspair(l)%j2=sh%dtss%particle(a1)%ind(i)
                sl%ss(a1)%epspair(l)%i1=ss%info%index_in_unitcell( sl%ss(a1)%pair(i)%j1 )
                sl%ss(a1)%epspair(l)%i2=ss%info%index_in_unitcell( sl%ss(a1)%pair(i)%j2 )
                sl%ss(a1)%epspair(l)%ui1=sh%unique_unitcell_atom( sl%ss(a1)%pair(i)%i1 )
                sl%ss(a1)%epspair(l)%ui2=sh%unique_unitcell_atom( sl%ss(a1)%pair(i)%i2 )
                sl%ss(a1)%epspair(l)%v=sh%dtss%particle(a1)%v(:,i)
                sl%ss(a1)%epspair(l)%unique_shell=0
                sl%ss(a1)%Zpair(l)%j1=a1
                sl%ss(a1)%Zpair(l)%j2=sh%dtss%particle(a1)%ind(i)
                sl%ss(a1)%Zpair(l)%i1=ss%info%index_in_unitcell( sl%ss(a1)%pair(i)%j1 )
                sl%ss(a1)%Zpair(l)%i2=ss%info%index_in_unitcell( sl%ss(a1)%pair(i)%j2 )
                sl%ss(a1)%Zpair(l)%ui1=sh%unique_unitcell_atom( sl%ss(a1)%pair(i)%i1 )
                sl%ss(a1)%Zpair(l)%ui2=sh%unique_unitcell_atom( sl%ss(a1)%pair(i)%i2 )
                sl%ss(a1)%Zpair(l)%v=sh%dtss%particle(a1)%v(:,i)
                sl%ss(a1)%Zpair(l)%unique_shell=0
            enddo
        endif
        if ( sl%dielectric_triplet ) then
            rcsq=sl%dc3**2
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%dc3 ) cycle
                do j=1,sh%dtss%particle(a1)%n
                    if ( sh%dtss%particle(a1)%d(j) .gt. sl%dc3 ) cycle
                    if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                    endif
                enddo
            enddo
            sl%ss(a1)%nZtriplet=l
            allocate(sl%ss(a1)%Ztriplet( sl%ss(a1)%nZtriplet ))
            l=0
            do i=1,sh%dtss%particle(a1)%n
                if ( sh%dtss%particle(a1)%d(i) .gt. sl%dc3 ) cycle
                do j=1,sh%dtss%particle(a1)%n
                    if ( sh%dtss%particle(a1)%d(j) .gt. sl%dc3 ) cycle
                    if ( lo_sqnorm(sh%dtss%particle(a1)%v(:,i)-sh%dtss%particle(a1)%v(:,j)) .lt. rcsq ) then
                        l=l+1
                        sl%ss(a1)%Ztriplet(l)%j1=a1
                        sl%ss(a1)%Ztriplet(l)%j2=sh%dtss%particle(a1)%ind(i)
                        sl%ss(a1)%Ztriplet(l)%j3=sh%dtss%particle(a1)%ind(j)
                        sl%ss(a1)%Ztriplet(l)%i1=ss%info%index_in_unitcell( sl%ss(a1)%triplet(i)%j1 )
                        sl%ss(a1)%Ztriplet(l)%i2=ss%info%index_in_unitcell( sl%ss(a1)%triplet(i)%j2 )
                        sl%ss(a1)%Ztriplet(l)%i3=ss%info%index_in_unitcell( sl%ss(a1)%triplet(i)%j3 )
                        sl%ss(a1)%Ztriplet(l)%ui1=sh%unique_supercell_atom( sl%ss(a1)%triplet(l)%j1 )
                        sl%ss(a1)%Ztriplet(l)%ui2=sh%unique_supercell_atom( sl%ss(a1)%triplet(l)%j2 )
                        sl%ss(a1)%Ztriplet(l)%ui3=sh%unique_supercell_atom( sl%ss(a1)%triplet(l)%j3 )
                        sl%ss(a1)%Ztriplet(l)%v1=0.0_flyt
                        sl%ss(a1)%Ztriplet(l)%v2=sh%dtss%particle(a1)%v(:,i)
                        sl%ss(a1)%Ztriplet(l)%v3=sh%dtss%particle(a1)%v(:,j)
                        sl%ss(a1)%Ztriplet(l)%unique_shell=0
                    endif
                enddo
            enddo
        endif
    enddo ssatomloop

    ! Done building all the tuplets. Could be a neat idea to locate all the unitcell ones
    ! in the supercell
    do a1=1,sl%nuc
        if ( sl%secondorder ) then
            do i=1,sl%uc(a1)%npair
                call locate_pair_in_supercell(sl%uc(a1)%pair(i),sl,ss,&
                     sl%uc(a1)%pair(i)%ssa,sl%uc(a1)%pair(i)%ssi)
            enddo
        endif
        if ( sl%thirdorder ) then
            do i=1,sl%uc(a1)%ntriplet
                call locate_triplet_in_supercell(sl%uc(a1)%triplet(i),sl,ss,&
                     sl%uc(a1)%triplet(i)%ssa,sl%uc(a1)%triplet(i)%ssi)
            enddo
        endif
        if ( sl%fourthorder ) then
            do i=1,sl%uc(a1)%nquartet
                call locate_quartet_in_supercell(sl%uc(a1)%quartet(i),sl,ss,&
                     sl%uc(a1)%quartet(i)%ssa,sl%uc(a1)%quartet(i)%ssi)
            enddo
        endif
        if ( sl%magnetic_pair_interactions ) then
            do i=1,sl%uc(a1)%nmagpair
                call locate_magpair_in_supercell(sl%uc(a1)%magpair(i),sl,ss,&
                     sl%uc(a1)%magpair(i)%ssa,sl%uc(a1)%magpair(i)%ssi)
            enddo
        endif
        if ( sl%dielectric_pair ) then
            do i=1,sl%uc(a1)%nepspair
                call locate_pair_in_supercell(sl%uc(a1)%epspair(i),sl,ss,&
                     sl%uc(a1)%epspair(i)%ssa,sl%uc(a1)%epspair(i)%ssi)
            enddo
            do i=1,sl%uc(a1)%nZpair
                call locate_pair_in_supercell(sl%uc(a1)%Zpair(i),sl,ss,&
                     sl%uc(a1)%Zpair(i)%ssa,sl%uc(a1)%Zpair(i)%ssi)
            enddo
        endif
        if ( sl%dielectric_triplet ) then
            do i=1,sl%uc(a1)%nZtriplet
                call locate_triplet_in_supercell(sl%uc(a1)%Ztriplet(i),sl,ss,&
                     sl%uc(a1)%Ztriplet(i)%ssa,sl%uc(a1)%Ztriplet(i)%ssi)
            enddo
        endif
    enddo

    ! Also, store a temporary list of all the unique tuplets of all orders
    if ( sl%secondorder ) then
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            l=l+sl%uc(a1)%npair
        enddo
        allocate(sh%allpairs(l))
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            do i=1,sl%uc(a1)%npair
                l=l+1
                sh%allpairs(l)=sl%uc(a1)%pair(i)
            enddo
        enddo
    endif
    if ( sl%thirdorder ) then
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            l=l+sl%uc(a1)%ntriplet
        enddo
        allocate(sh%alltriplets(l))
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            do i=1,sl%uc(a1)%ntriplet
                l=l+1
                sh%alltriplets(l)=sl%uc(a1)%triplet(i)
            enddo
        enddo
    endif
    if ( sl%fourthorder ) then
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            l=l+sl%uc(a1)%nquartet
        enddo
        allocate(sh%allquartets(l))
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            do i=1,sl%uc(a1)%nquartet
                l=l+1
                sh%allquartets(l)=sl%uc(a1)%quartet(i)
            enddo
        enddo
    endif
    if ( sl%magnetic_pair_interactions ) then
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            l=l+sl%uc(a1)%nmagpair
        enddo
        allocate(sh%allmagpairs(l))
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            do i=1,sl%uc(a1)%nmagpair
                l=l+1
                sh%allmagpairs(l)=sl%uc(a1)%magpair(i)
            enddo
        enddo
    endif
    if ( sl%dielectric_pair ) then
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            l=l+sl%uc(a1)%nepspair
        enddo
        allocate(sh%allepspairs(l))
        allocate(sh%allZpairs(l))
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            do i=1,sl%uc(a1)%nepspair
                l=l+1
                sh%allepspairs(l)=sl%uc(a1)%epspair(i)
                sh%allZpairs(l)=sl%uc(a1)%Zpair(i)
            enddo
        enddo
    endif
    if ( sl%thirdorder ) then
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            l=l+sl%uc(a1)%nZtriplet
        enddo
        allocate(sh%allZtriplets(l))
        l=0
        do a1=1,uc%na
            if ( sh%is_unitcell_atom_prototype(a1) .eqv. .false. ) cycle
            do i=1,sl%uc(a1)%nZtriplet
                l=l+1
                sh%allZtriplets(l)=sl%uc(a1)%Ztriplet(i)
            enddo
        enddo
    endif

    ! This is also a reasonable place to create all hashes, that I use to speed up the
    ! symmetry detection.
    if ( sl%secondorder ) then
        ! hash all the tuplets, first the unitcell
        l=0
        do a1=1,sl%nuc
            l=l+sl%uc(a1)%npair
        enddo
        lo_allocate(sh%pair_uchash(l))
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%npair
            l=l+1
            sh%pair_uchash(l)=hashtuplet(sl%uc(a1)%pair(i))
        enddo
        enddo
        ! then the supercell
        l=0
        do a1=1,sl%nss
            l=l+sl%ss(a1)%npair
        enddo
        lo_allocate(sh%pair_sshash(l))
        l=0
        do a1=1,sl%nss
        do i=1,sl%ss(a1)%npair
            l=l+1
            sh%pair_sshash(l)=hashtuplet( sl%ss(a1)%pair(i) )
        enddo
        enddo
        ! and the unique
        lo_allocate(sh%pair_allhash(size(sh%allpairs,1)))
        do i=1,size(sh%allpairs,1)
            sh%pair_allhash(i)=hashtuplet( sh%allpairs(i) )
        enddo
    endif
    if ( sl%thirdorder ) then
        ! hash all the tuplets, first the unitcell
        l=0
        do a1=1,sl%nuc
            l=l+sl%uc(a1)%ntriplet
        enddo
        lo_allocate(sh%triplet_uchash(l))
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%ntriplet
            l=l+1
            sh%triplet_uchash(l)=hashtuplet(sl%uc(a1)%triplet(i))
        enddo
        enddo
        ! then the supercell
        l=0
        do a1=1,sl%nss
            l=l+sl%ss(a1)%ntriplet
        enddo
        lo_allocate(sh%triplet_sshash(l))
        l=0
        do a1=1,sl%nss
        do i=1,sl%ss(a1)%ntriplet
            l=l+1
            sh%triplet_sshash(l)=hashtuplet( sl%ss(a1)%triplet(i) )
        enddo
        enddo
        ! and the unique
        lo_allocate(sh%triplet_allhash(size(sh%alltriplets,1)))
        do i=1,size(sh%alltriplets,1)
            sh%triplet_allhash(i)=hashtuplet( sh%alltriplets(i) )
        enddo
    endif
    if ( sl%fourthorder ) then
        ! hash all the tuplets, first the unitcell
        l=0
        do a1=1,sl%nuc
            l=l+sl%uc(a1)%nquartet
        enddo
        lo_allocate(sh%quartet_uchash(l))
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nquartet
            l=l+1
            sh%quartet_uchash(l)=hashtuplet(sl%uc(a1)%quartet(i))
        enddo
        enddo
        ! then the supercell
        l=0
        do a1=1,sl%nss
            l=l+sl%ss(a1)%nquartet
        enddo
        lo_allocate(sh%quartet_sshash(l))
        l=0
        do a1=1,sl%nss
        do i=1,sl%ss(a1)%nquartet
            l=l+1
            sh%quartet_sshash(l)=hashtuplet( sl%ss(a1)%quartet(i) )
        enddo
        enddo
        ! and the unique
        lo_allocate(sh%quartet_allhash(size(sh%allquartets,1)))
        do i=1,size(sh%allquartets,1)
            sh%quartet_allhash(i)=hashtuplet( sh%allquartets(i) )
        enddo
    endif
    if ( sl%magnetic_pair_interactions ) then
        ! hash all the tuplets, first the unitcell
        l=0
        do a1=1,sl%nuc
            l=l+sl%uc(a1)%nmagpair
        enddo
        lo_allocate(sh%magpair_uchash(l))
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nmagpair
            l=l+1
            sh%magpair_uchash(l)=hashtuplet(sl%uc(a1)%magpair(i))
        enddo
        enddo
        ! then the supercell
        l=0
        do a1=1,sl%nss
            l=l+sl%ss(a1)%nmagpair
        enddo
        lo_allocate(sh%magpair_sshash(l))
        l=0
        do a1=1,sl%nss
        do i=1,sl%ss(a1)%nmagpair
            l=l+1
            sh%magpair_sshash(l)=hashtuplet( sl%ss(a1)%magpair(i) )
        enddo
        enddo
        ! and the unique
        lo_allocate(sh%magpair_allhash(size(sh%allmagpairs,1)))
        do i=1,size(sh%allmagpairs,1)
            sh%magpair_allhash(i)=hashtuplet( sh%allmagpairs(i) )
        enddo
    endif
    if ( sl%dielectric_pair ) then
        ! hash all the tuplets, first the unitcell
        l=0
        do a1=1,sl%nuc
            l=l+sl%uc(a1)%nepspair
        enddo
        allocate(sh%epspair_uchash(l))
        allocate(sh%Zpair_uchash(l))
        l=0
        do a1=1,sl%nuc
        do i=1,sl%uc(a1)%nepspair
            l=l+1
            sh%epspair_uchash(l)=hashtuplet(sl%uc(a1)%epspair(i))
            sh%Zpair_uchash(l)=hashtuplet(sl%uc(a1)%Zpair(i))
        enddo
        enddo
        ! then the supercell
        l=0
        do a1=1,sl%nss
            l=l+sl%ss(a1)%nepspair
        enddo
        allocate(sh%epspair_sshash(l))
        allocate(sh%Zpair_sshash(l))
        l=0
        do a1=1,sl%nss
        do i=1,sl%ss(a1)%nepspair
            l=l+1
            sh%epspair_sshash(l)=hashtuplet( sl%ss(a1)%epspair(i) )
            sh%Zpair_sshash(l)=hashtuplet( sl%ss(a1)%Zpair(i) )
        enddo
        enddo
        ! and the unique
        allocate(sh%epspair_allhash(size(sh%allepspairs,1)))
        allocate(sh%Zpair_allhash(size(sh%allZpairs,1)))
        do i=1,size(sh%allepspairs,1)
            sh%epspair_allhash(i)=hashtuplet( sh%allepspairs(i) )
            sh%Zpair_allhash(i)=hashtuplet( sh%allZpairs(i) )
        enddo
    endif
end subroutine

!> just the operations, per shell, that leave that pair invariant.
module subroutine validops_pairshells(sl,sh,ss)
    type(lo_symlist), intent(inout) :: sl
    type(lo_symtabhelper), intent(inout) :: sh
    type(lo_crystalstructure), intent(in) :: ss

    integer, dimension(:), allocatable :: opind
    integer :: shell,op,opctr

    lo_allocate(opind(sl%npairop))
    do shell=1,sl%npairshells
        ! Figure out which operations take the pair back to itself.
        opind=0
        opctr=0
        oploop: do op=1,sl%npairop
            if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,sl%pairshell(shell)%protpair,sl%pairshell(shell)%protpair,sl%pairop(op)) ) then
                opctr=opctr+1
                opind(opctr)=op
            endif
        enddo oploop
        ! store these
        sl%pairshell(shell)%noperations=opctr
        lo_allocate(sl%pairshell(shell)%operations(opctr))
        sl%pairshell(shell)%operations=opind(1:opctr)
    enddo
    lo_deallocate(opind)
end subroutine

!> check if two tuplets are equivalent after an operation
module function compare_tuplets_after_one_is_rotated(sl,sh,ss,t1,t2,op,testperm) result(match)
    type(lo_symlist), intent(in) :: sl
    type(lo_symtabhelper), intent(in) :: sh
    type(lo_crystalstructure), intent(in) :: ss
    class(lo_symlist_tuplet), intent(in) :: t1,t2
    class(lo_symlist_tupletop), intent(in) :: op
    logical, intent(in), optional :: testperm
    logical :: match

    integer :: ssi1,ssa1,a1,a2
    real(flyt), dimension(3) :: v0
    integer, dimension(2) :: pj1,pj2,pi1,pi2
    integer, dimension(3) :: tj1,ti1,ti2
    integer, dimension(4) :: qj1,qi1,qi2
    logical :: donesomething,checkperm
    donesomething=.false.

    ! test both tr*op and op*tr
    if ( present(testperm) ) then
        checkperm=testperm
    else
        checkperm=.false.
    endif

    ! Stupid polymorphism, I hate it. Now we are at pairs, perhaps.
    select type(t1); class is(lo_symlist_pair)
    select type(t2); class is(lo_symlist_pair)
    select type(op); type is(lo_symlist_pairop)
        ! at least if I make it here, my overly ambitious polymorphic stuff is working.
        donesomething=.true.
        ! these are unitcell indices
        pi1(1)=t1%i1
        pi1(2)=t1%i2
        pi2(1)=t2%i1
        pi2(2)=t2%i2
        pj1(1)=sh%unique_unitcell_atom(pi1(1))
        pj1(2)=sh%unique_unitcell_atom(pi1(2))
        pj2(1)=sh%unique_unitcell_atom(pi2(1))
        pj2(2)=sh%unique_unitcell_atom(pi2(2))
        pi1=pi1(op%perm)
        pj1=pj1(op%perm)
        ! So, with this operation is basically
        ! op*t1%i1 = t2%i2, or op*pi1(1)=pi2(2)
        ! so that operation needs to be valid
        if ( sh%clist(pi1(1),pi2(1),op%opind) .eqv. .false. ) then
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
        if ( lo_sqnorm(v0) .lt. sl%sqtol_cart ) then
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
    select type(t1); class is(lo_symlist_uctriplet)
    select type(t2)
        class is(lo_symlist_uctriplet)
        select type(op); type is(lo_symlist_tripletop)
            donesomething=.true.
            ! this version is comparing two unitcell triplets.
            ! fetch where in the supercell t1 is.
            ssa1=t1%ssa
            ssi1=t1%ssi
            ! get the supercell indices for the first one
            tj1=[sl%ss(ssa1)%triplet(ssi1)%j1,sl%ss(ssa1)%triplet(ssi1)%j2,sl%ss(ssa1)%triplet(ssi1)%j3]
            ! permute
            tj1=tj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell(tj1(1))
            a2=ss%info%index_in_unitcell(t2%ssa)
            ! are there valid operations between these?
            if ( sh%clist(a1,a2,op%opind) .eqv. .false. ) then
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
            call locate_triplet_in_supercell_via_indices(tj1,sl,ssa1,ssi1)

            ! now compare with rotations
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%triplet(ssi1)%v2)-t2%v2) .lt. sl%sqtol_cart ) then
                if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%triplet(ssi1)%v3)-t2%v3) .lt. sl%sqtol_cart ) then
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
        class is(lo_symlist_sstriplet)
        select type(op); type is(lo_symlist_tripletop)
            donesomething=.true.
            ! now I am comparing the first triplet from the unitcell, and the second triplet is from the
            ! supercell.
            ! fetch where in the supercell t1 is.
            ssa1=t1%ssa
            ssi1=t1%ssi
            ! get the supercell indices for the first one
            tj1=[sl%ss(ssa1)%triplet(ssi1)%j1,sl%ss(ssa1)%triplet(ssi1)%j2,sl%ss(ssa1)%triplet(ssi1)%j3]
            ! permute
            tj1=tj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell(tj1(1))
            a2=ss%info%index_in_unitcell(t2%j1)
            ! are there valid operations betwee these?
            if ( sh%clist(a1,a2,op%opind) .eqv. .false. ) then
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
            call locate_triplet_in_supercell_via_indices(tj1,sl,ssa1,ssi1)

            ! now compare with rotations
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%triplet(ssi1)%v2)-t2%v2) .lt. sl%sqtol_cart ) then
                if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%triplet(ssi1)%v3)-t2%v3) .lt. sl%sqtol_cart ) then
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
    select type(t1); class is(lo_symlist_ucquartet)
    select type(t2)
        class is(lo_symlist_ucquartet)
        select type(op); type is(lo_symlist_quartetop)
            donesomething=.true.
            ! this version is comparing two unitcell triplets.
            ! fetch where in the supercell t1 is.
            ssa1=t1%ssa
            ssi1=t1%ssi
            ! get the supercell indices for the first one
            qj1=[sl%ss(ssa1)%quartet(ssi1)%j1,sl%ss(ssa1)%quartet(ssi1)%j2,&
                 sl%ss(ssa1)%quartet(ssi1)%j3,sl%ss(ssa1)%quartet(ssi1)%j4]
            ! permute
            qj1=qj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell(qj1(1))
            a2=ss%info%index_in_unitcell(t2%ssa)
            ! are there valid operations betwee these?
            if ( sh%clist(a1,a2,op%opind) .eqv. .false. ) then
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

            ! Fetch indices for the permuted
            call locate_quartet_in_supercell_via_indices(qj1,sl,ssa1,ssi1)

            ! now compare with rotations
            match=.false.
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v2)-t2%v2) .lt. sl%sqtol_cart ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v3)-t2%v3) .lt. sl%sqtol_cart ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v4)-t2%v4) .lt. sl%sqtol_cart ) then
                match=.true.
            endif
            endif
            endif
        end select
        class is(lo_symlist_ssquartet)
        select type(op); type is(lo_symlist_quartetop)
            donesomething=.true.
            ! now I am comparing the first quartet from the unitcell, and the second
            ! quartet is from the supercell.
            ! fetch where in the supercell t1 is.
            ssa1=t1%ssa
            ssi1=t1%ssi
            ! get the supercell indices for the first one
            qj1=[sl%ss(ssa1)%quartet(ssi1)%j1,sl%ss(ssa1)%quartet(ssi1)%j2,&
                 sl%ss(ssa1)%quartet(ssi1)%j3,sl%ss(ssa1)%quartet(ssi1)%j4]
            ! permute
            qj1=qj1(op%perm)
            ! now I know the two possible starting atoms.
            a1=ss%info%index_in_unitcell(qj1(1))
            a2=ss%info%index_in_unitcell(t2%j1)
            ! are there valid operations between these?
            if ( sh%clist(a1,a2,op%opind) .eqv. .false. ) then
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

            ! Fetch indices for the permuted
            call locate_quartet_in_supercell_via_indices(qj1,sl,ssa1,ssi1)

            ! now compare with rotations
            match=.false.
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v2)-t2%v2) .lt. sl%sqtol_cart ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v3)-t2%v3) .lt. sl%sqtol_cart ) then
            if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v4)-t2%v4) .lt. sl%sqtol_cart ) then
                match=.true.
            endif
            endif
            endif
        end select
    end select
    end select
    ! For debugging, remove later
    !select type(t1); class is(lo_symlist_ssquartet)
    !select type(t2); class is(lo_symlist_ssquartet)
    !select type(op); type is(lo_symlist_quartetop)
    !        donesomething=.true.
    !        ! now I am comparing the first quartet from the supercell, and the second
    !        ! quartet is from the supercell.
    !        qj1=[t1%j1,t1%j2,t1%j3,t1%j4]
    !        ! permute
    !        qj1=qj1(op%perm)
    !        ! now I know the two possible starting atoms.
    !        a1=ss%info%index_in_unitcell(qj1(1))
    !        a2=ss%info%index_in_unitcell(t2%j1)
    !        ! are there valid operations between these?
    !        if ( sh%clist(a1,a2,op%opind) .eqv. .false. ) then
    !            match=.false.
    !            return
    !        endif
    !        ! test if all the atoms are the same
    !        qi1=[t1%ui1,t1%ui2,t1%ui3,t1%ui4]
    !        qi2=[t2%ui1,t2%ui2,t2%ui3,t2%ui4]
    !        if ( sum(abs(qi1(op%perm)-qi2)) .ne. 0 ) then
    !            match=.false.
    !            return
    !        endif

    !        ! Fetch indices for the permuted
    !        call locate_quartet_in_supercell_via_indices(qj1,sl,ssa1,ssi1)

    !        ! now compare with rotations
    !        match=.false.
    !        if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v2)-t2%v2) .lt. sl%sqtol_cart ) then
    !        if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v3)-t2%v3) .lt. sl%sqtol_cart ) then
    !        if ( lo_sqnorm(matmul(op%m3,sl%ss(ssa1)%quartet(ssi1)%v4)-t2%v4) .lt. sl%sqtol_cart ) then
    !            match=.true.
    !        endif
    !        endif
    !        endif
    !end select
    !end select
    !end select

    select type(t1); class is(lo_symlist_magpair)
    select type(t2); class is(lo_symlist_magpair)
    select type(op); type is(lo_symlist_pairop)
        ! at least if I make it here, my overly ambitious polymorphic stuff is working.
        donesomething=.true.
        ! these are unitcell indices
        pi1(1)=t1%i1
        pi1(2)=t1%i2
        pi2(1)=t2%i1
        pi2(2)=t2%i2
        pj1(1)=sh%unique_unitcell_atom(pi1(1))
        pj1(2)=sh%unique_unitcell_atom(pi1(2))
        pj2(1)=sh%unique_unitcell_atom(pi2(1))
        pj2(2)=sh%unique_unitcell_atom(pi2(2))
        pi1=pi1(op%perm)
        pj1=pj1(op%perm)
        ! So, with this operation is basically
        ! op*t1%i1 = t2%i2, or op*pi1(1)=pi2(2)
        ! so that operation needs to be valid
        if ( sh%clist(pi1(1),pi2(1),op%opind) .eqv. .false. ) then
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
        if ( lo_sqnorm(v0) .lt. sl%sqtol_cart ) then
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

    if ( donesomething .eqv. .false. ) then
        write(*,*) 'I failed badly at comparing stuff, read up on polymorphism'
        stop
    endif
end function

!> This maps out the unique atoms
module subroutine findunique(sl,sh,uc,ss)
    type(lo_symlist), intent(inout) :: sl
    type(lo_symtabhelper), intent(inout) :: sh
    type(lo_crystalstructure), intent(in) :: uc
    type(lo_crystalstructure), intent(in) :: ss

    real(flyt), dimension(3) :: v0
    real(flyt) :: f0
    integer, dimension(:), allocatable :: di1,di2,di3
    integer :: i,j,l,o
    integer :: a1,a2,uca

    ! generate distance tables for the unit cell and supercell
    if ( sl%wzcutoff ) then
        ! Use wigner-seitz cutoff for the second order
        allocate(lo_voronoi_distancetable::sh%dtuc)
        allocate(lo_voronoi_distancetable::sh%dtss)
        select type(a=>sh%dtuc); type is(lo_voronoi_distancetable)
            if ( sl%wzdim(1) .eq. -2 ) then
                ! Secret option don't bother thinking about it
                call a%generate_wzcutoff(uc%r,uc%latticevectors,[1,1,1],ss%latticevectors)
            else
                call a%generate_wzcutoff(uc%r,uc%latticevectors,sl%wzdim,uc%latticevectors,vorona=uc%na)
            endif
            ! some sanity checks for this. It will get weird if the WZ cutoff is smaller than
            ! any of the other cutoffs, and makes no sense.
            if ( a%cutoff .lt. sl%rc3 ) then
                write(*,*) 'ERROR: The Wigner-Seitz cutoff is smaller than the third order cutoff'
                write(*,*) '    This might work if I put some work into it, but not today.'
                stop
            endif
            if ( a%cutoff .lt. sl%rc4 ) then
                write(*,*) 'ERROR: The Wigner-Seitz cutoff is smaller than the fourth order cutoff'
                write(*,*) '    This might work if I put some work into it, but not today.'
                stop
            endif
            if ( a%cutoff .lt. sl%mc2 ) then
                write(*,*) 'ERROR: The Wigner-Seitz cutoff is smaller than the magnetic pair cutoff'
                write(*,*) '    This might work if I put some work into it, but not today.'
                stop
            endif
            if ( a%cutoff .lt. sl%dc2 ) then
                write(*,*) 'ERROR: The Wigner-Seitz cutoff is smaller than the dielectric pair cutoff'
                write(*,*) '    This might work if I put some work into it, but not today.'
                stop
            endif
            if ( a%cutoff .lt. sl%dc3 ) then
                write(*,*) 'ERROR: The Wigner-Seitz cutoff is smaller than the dielectric triplet cutoff'
                write(*,*) '    This might work if I put some work into it, but not today.'
                stop
            endif
        end select
        select type(a=>sh%dtss); type is(lo_voronoi_distancetable)
            if ( sl%wzdim(1) .eq. -2 ) then
                call a%generate_wzcutoff(ss%r,ss%latticevectors,[1,1,1],ss%latticevectors)
            else
                call a%generate_wzcutoff(ss%r,ss%latticevectors,sl%wzdim,uc%latticevectors,vorona=uc%na)
            endif
        end select
    elseif ( sl%nj2 .ge. 1 ) then
        call lo_stop_gracefully(['Not done implementing graph distancetables'],lo_exitcode_param,__FILE__,__LINE__)
    else
        ! make perfectly normal distance tables
        allocate(lo_distancetable::sh%dtuc,stat=lo_status)
        allocate(lo_distancetable::sh%dtss,stat=lo_status)
        f0=max(sl%rc2,sl%rc3,sl%rc4,sl%mc2,sl%dc2,sl%dc3)
        select type(a=>sh%dtuc); type is(lo_distancetable)
            call a%generate(uc%r,uc%latticevectors,f0,verbosity=0)
        end select
        select type(a=>sh%dtss); type is(lo_distancetable)
            call a%generate(ss%r,ss%latticevectors,f0,verbosity=0)
        end select
    endif

    ! A smart thing to do is to make sure that the distance tables are sorted identically, such that
    ! for each atom in the supercell (i) that correspond to an atom in the unitcell (j), then neighbour
    ! n from both (i) and (j) are the same. This makes matching the tuplets from the prototypes to the
    ! supercell fast, and way easier. Also, this is an early and good place for some sanity checks.
    do a1=1,sh%dtss%np
        uca=ss%info%index_in_unitcell(a1)
        if ( sh%dtuc%particle(uca)%n .ne. sh%dtss%particle(a1)%n ) then
            write(*,*) 'ERROR:'
            write(*,*) 'atom ',tochar(uca),' in the unitcell and atom ',tochar(a1),'in the supercell should have the same'
            write(*,*) 'number of neighbours within the cutoff, but they do not.'
            stop
        endif
        lo_allocate(di1(sh%dtuc%particle(uca)%n))
        di1=0
        p1l: do i=1,sh%dtss%particle(a1)%n
            do j=1,sh%dtuc%particle(uca)%n
                ! same distance
                if ( abs(sh%dtuc%particle(uca)%d(j)-sh%dtss%particle(a1)%d(i)) .gt. sl%tol_cart ) cycle
                ! same vectors
                v0=sh%dtuc%particle(uca)%v(:,j)-sh%dtss%particle(a1)%v(:,i)
                if ( lo_sqnorm( v0 ) .lt. sl%sqtol_cart ) then
                    di1(j)=i
                    cycle p1l
                endif
            enddo
            ! if we make it here, there is a discrepancy in the distance tables
            write(*,*) 'ERROR:'
            write(*,*) '    discrepancy in the distance tables between the unit and supercell'
            write(*,*) '    make sure they are generated from consistent cells.'
            stop
        enddo p1l
        ! now I just reorder everything
        sh%dtss%particle(a1)%ind=sh%dtss%particle(a1)%ind(di1)
        sh%dtss%particle(a1)%d=sh%dtss%particle(a1)%d(di1)
        sh%dtss%particle(a1)%v=sh%dtss%particle(a1)%v(:,di1)
        sh%dtss%particle(a1)%lv=sh%dtss%particle(a1)%lv(:,di1)
        if ( allocated(sh%dtss%particle(a1)%weight) ) then
            sh%dtss%particle(a1)%weight=sh%dtss%particle(a1)%weight(di1)
        endif
        lo_deallocate(di1)
    enddo

    ! Test that I actually got this right:
    do a1=1,sh%dtss%np
        uca=ss%info%index_in_unitcell(a1)
        ! discrepancy for the entire list:
        f0=sum(abs(sh%dtss%particle(a1)%v-sh%dtuc%particle(uca)%v))/sh%dtss%particle(a1)%n
        if ( f0 .gt. sl%tol_cart/100 ) then
            write(*,*) 'ERROR:'
            write(*,*) '    There was too large discrepancy between the neighbour list generated from'
            write(*,*) '    the unitcell and supercell. Acceptable discrepancy is essentially zero.'
            write(*,*) '    There is no valid reason to have any difference whatsoever, so make sure'
            write(*,*) '    The unit and supercell describe the same lattice, to many many digits.'
            stop
        endif
    enddo

    ! Figure out which the unique atoms are.
    lo_allocate(di1(uc%na))
    lo_allocate(di2(uc%na))
    lo_allocate(di3(uc%na))
    di1=1
    di2=0
    di3=0
    do a1=1,uc%na
        if ( di1(a1) .eq. 0 ) cycle
        di2(a1)=a1
        di3(a1)=1
        a3loop: do a2=a1+1,uc%na
            ! trivial, fast version
            do o=1,sl%sym%n
                ! trivial check
                if ( stars_equal_after_operation(sl,uc,uc,sl%sym%op(o),sh%dtuc%particle(a2),sh%dtuc%particle(a1)) ) then
                    di1(a2)=0
                    di2(a2)=a1
                    di3(a2)=o
                    cycle a3loop
                endif
            enddo
        enddo a3loop
    enddo

    if ( sl%verbosity .gt. 1 ) then
        do a1=1,uc%na
            write(*,*) '     atom ',tochar(a1),' equivalent to atom ',tochar(di2(a1)),' with operation ',tochar(di3(a1))
        enddo
    endif

    ! Store this information in a neat way.
    lo_allocate(sh%is_unitcell_atom_prototype(uc%na))
    lo_allocate(sh%unique_unitcell_atom(uc%na))
    lo_allocate(sh%unique_unitcell_operation(uc%na))
    lo_allocate(sh%unique_supercell_atom(ss%na))
    lo_allocate(sh%unique_supercell_operation(ss%na))
    sh%is_unitcell_atom_prototype=.false.
    sh%unique_unitcell_atom=0
    sh%unique_unitcell_operation=0
    sh%unique_supercell_atom=0
    sh%unique_supercell_operation=0

    ! Now sanity test this, for both unit and supercell. If I was a real ninja, it's going to work.
    l=0
    do a1=1,uc%na
        a2=di2(a1)
        o=di3(a1)
        sh%unique_unitcell_atom(a1)=a2
        sh%unique_unitcell_operation(a1)=o
        if ( di1(a1) .eq. 1 ) sh%is_unitcell_atom_prototype(a1)=.true.
        ! sanity test
        if ( stars_equal_after_operation(sl,uc,uc,sl%sym%op(o),sh%dtuc%particle(a1),sh%dtuc%particle(a2)) ) then
            l=l+1
        endif
    enddo
    if ( uc%na-l .ne. 0 ) then
        write(*,*) 'ERROR:'
        write(*,*) '  failed the symmetry mapping from unitcell to unitcell when trying to find unique atoms.'
        stop
    endif

    l=0
    do a1=1,ss%na
        i=ss%info%index_in_unitcell(a1)
        a2=di2(i)
        o=di3(i)
        sh%unique_supercell_atom(a1)=a2
        sh%unique_supercell_operation(a1)=o
        if ( stars_equal_after_operation(sl,ss,uc,sl%sym%op(o),sh%dtss%particle(a1),sh%dtuc%particle(a2)) ) then
            l=l+1
        endif
    enddo
    if ( ss%na-l .ne. 0 ) then
        write(*,*) 'ERROR:'
        write(*,*) '  failed the symmetry mapping from unitcell to supercell after I mapped out the unique atoms.'
        stop
    endif

    ! Seems like it was a decent matching. Now I want to know what operations are valid between what atoms.
    ! The reson I do it in two steps, is that I only have to test between atoms that are already identified
    ! as equivalent, should make it a little faster.
    lo_allocate(sh%clist(uc%na,uc%na,sl%sym%n))
    sh%clist=.false.
    do a1=1,uc%na
    do a2=1,uc%na
        if ( sh%unique_unitcell_atom(a1) .ne. sh%unique_unitcell_atom(a2) ) cycle
        do o=1,sl%sym%n
            if ( stars_equal_after_operation(sl,uc,uc,sl%sym%op(o),sh%dtuc%particle(a1),sh%dtuc%particle(a2)) ) then
                sh%clist(a1,a2,o)=.true.
            endif
        enddo
    enddo
    enddo
end subroutine

! Below are things that are not exposed outside this submodule

!> compare stars of between two atoms @todo Should add a hashing+sorting guy here to make it close to NlogN
function stars_equal_after_operation(sl,p1,p2,op,s1,s2) result(match)
    !> settings and stuff
    class(lo_symlist), intent(in) :: sl
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p1,p2
    !> operation
    type(lo_spacegroup_operation), intent(in) :: op
    !> the stars of vectors
    type(lo_distancetable_particle), intent(in) :: s1,s2
    !> do they match
    logical :: match

    integer :: a1,a2,i,j,k,l
    integer, dimension(s1%n) :: sp1,sp2
    real(flyt), dimension(3) :: v0,v1

    ! Indices in the unit cell to the atoms
    a1=s1%ind(1)
    a2=s2%ind(1)
    ! Start by checking the simplest things first
    if ( p1%atomic_number(a1) .ne. p2%atomic_number(a2) ) then
        match=.false.
        return
    endif
    if ( s1%n .ne. s2%n ) then
        match=.false.
        return
    endif
    do i=1,s1%n
        sp1(i)=p1%atomic_number( s1%ind(i) )
        sp2(i)=p2%atomic_number( s2%ind(i) )
    enddo
    if ( sl%testcoll ) then
        if ( abs(p1%mag%collinear_moment(a1)-p2%mag%collinear_moment(a2)) .gt. sl%tol_mag ) then
            match=.false.
            return
        endif
    endif
    if ( sl%testnoncoll ) then
        v0=p1%mag%noncollinear_moment(:,a1)-p2%mag%noncollinear_moment(:,a2)
        if ( lo_sqnorm(v0) .gt. sl%sqtol_mag ) then
            match=.false.
            return
        endif
    endif

    ! Then test the annoying stuff
    l=0
    s1l: do i=1,s1%n
        v0=matmul(op%m,s1%v(:,i))
        k=0
        do j=1,s2%n
            v1=v0-s2%v(:,j)
            if ( sp1(i) .ne. sp2(j) ) cycle
            if ( lo_sqnorm(v1) .gt. sl%sqtol_cart ) cycle
            if ( sl%testcoll ) then
                if ( abs(p1%mag%collinear_moment( s1%ind(i) )-p2%mag%collinear_moment( s2%ind(j) )) .gt. sl%tol_mag ) cycle
            endif
            if ( sl%testnoncoll ) then
                v1=p1%mag%noncollinear_moment( :,s1%ind(i) )-p2%mag%noncollinear_moment( :,s2%ind(j) )
                if ( lo_sqnorm(v1) .gt. sl%sqtol_mag ) cycle
            endif
            ! if I made it here, I think it matches!
            k=1
            l=l+1
            cycle s1l
        enddo
        if ( k .eq. 0 ) then
            match=.false.
            return
        endif
    enddo s1l

    if ( l .eq. s1%n ) then
        match=.true.
    else
        match=.false.
    endif
end function

!> return the unique tuplets from a long list of tuplets
subroutine return_unique_tuplets(sl,sh,tuplets,hash,ss,unique_tuplets)
    !> settings and stuff
    type(lo_symlist), intent(in) :: sl
    type(lo_symtabhelper), intent(in) :: sh
    !> tuplets to reduce
    class(lo_symlist_tuplet), dimension(:), intent(inout) :: tuplets
    !> hash of the tuplets
    real(flyt), dimension(:), intent(in) :: hash
    !> structure
    type(lo_crystalstructure), intent(in) :: ss
    !> resulting tuplets
    class(lo_symlist_tuplet), dimension(:), allocatable, intent(out) :: unique_tuplets

    class(lo_symlist_tuplet), dimension(:), allocatable :: dumtup
    real(flyt), dimension(:), allocatable :: unhash
    real(flyt) :: t0
    integer :: i,n,t1,t2,op,unique_counter
    logical :: newtuplet

    ! start the timer
    t0=walltime()
    ! number of tuplets
    n=size(tuplets,1)
    ! fetch the hashes, and create an empty list of tuplets
    select type(tuplets)
        type is(lo_symlist_ucpair)
            lo_allocate(lo_symlist_ucpair::dumtup(n))
        class is(lo_symlist_uctriplet)
            lo_allocate(lo_symlist_uctriplet::dumtup(n))
        class is(lo_symlist_ucquartet)
            lo_allocate(lo_symlist_ucquartet::dumtup(n))
        class is(lo_symlist_ucmagpair)
            lo_allocate(lo_symlist_ucmagpair::dumtup(n))
        class default
            call lo_stop_gracefully(['Undefined kind of tuplet. Thinking is hard.'],lo_exitcode_param,__FILE__,__LINE__)
    end select

    ! Create space for the unique hashes
    lo_allocate(unhash(n))
    unhash=-lo_huge

    if ( sl%verbosity .gt. 0 ) then
        call lo_progressbar_init()
    endif

    ! The first tuplet is already stored there.
    unique_counter=0
    do t1=1,n
        ! I assume this is a brand new tuplet
        newtuplet=.true.
        ! now we take the next tuplet, and compare with those that are already there
        t2l: do t2=1,unique_counter
            ! hash has to match
            if ( abs(hash(t1)-unhash(t2)) .gt. lo_tol ) cycle
            ! Compare with all of them
            select type(tuplets)
            class is(lo_symlist_ucpair)
                opl1: do op=1,sl%npairop
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                         tuplets( t1 ),dumtup( t2 ),sl%pairop(op) ) ) then
                        newtuplet=.false.
                        exit t2l
                    endif
                enddo opl1
            class is(lo_symlist_triplet)
                opl2: do op=1,sl%ntripletop
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                         tuplets( t1 ),dumtup( t2 ),sl%tripletop(op) ) ) then
                        newtuplet=.false.
                        exit t2l
                    endif
                enddo opl2
            class is(lo_symlist_quartet)
                opl3: do op=1,sl%nquartetop
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                         tuplets( t1 ),dumtup( t2 ),sl%quartetop(op) ) ) then
                        newtuplet=.false.
                        exit t2l
                    endif
                enddo opl3
            class is(lo_symlist_ucmagpair)
                opl4: do op=1,sl%npairop
                    ! permutations are confusing. Skip those.
                    if ( sl%pairop(op)%permind .eq. 2 ) cycle
                    if ( compare_tuplets_after_one_is_rotated(sl,sh,ss,&
                         tuplets( t1 ),dumtup( t2 ),sl%pairop(op) ) ) then
                        newtuplet=.false.
                        exit t2l
                    endif
                enddo opl4
            end select
        enddo t2l

        ! Now, I might have a new unique tuplet!
        if ( newtuplet ) then
            ! update the counter
            unique_counter=unique_counter+1
            ! add the hash to the list of unique hashes
            unhash(unique_counter)=hash(t1)
            ! and store the tuplet
            select type(tuplets)
            type is(lo_symlist_ucpair)
                select type(dumtup); type is(lo_symlist_ucpair)
                    dumtup(unique_counter)=tuplets(t1)
                end select
            type is(lo_symlist_uctriplet)
                select type(dumtup); type is(lo_symlist_uctriplet)
                    dumtup(unique_counter)=tuplets(t1)
                end select
            type is(lo_symlist_ucquartet)
                select type(dumtup); type is(lo_symlist_ucquartet)
                    dumtup(unique_counter)=tuplets(t1)
                end select
            type is(lo_symlist_ucmagpair)
                select type(dumtup); type is(lo_symlist_ucmagpair)
                    dumtup(unique_counter)=tuplets(t1)
                end select
            end select
        endif
        ! and tell the world how it's going
        if ( sl%verbosity .gt. 0 ) then
            select type(tuplets)
            class is(lo_symlist_pair)
                call lo_progressbar(' ... reducing pairs',t1,n,walltime()-t0)
            class is(lo_symlist_triplet)
                call lo_progressbar(' ... reducing triplets',t1,n,walltime()-t0)
            class is(lo_symlist_quartet)
                call lo_progressbar(' ... reducing quartets',t1,n,walltime()-t0)
            class is(lo_symlist_magpair)
                call lo_progressbar(' ... reducing magnetic pairs',t1,n,walltime()-t0)
            end select
        endif
    enddo

    ! Now try to return the unique
    select type(dumtup)
    type is(lo_symlist_ucpair)
        allocate(lo_symlist_ucpair::unique_tuplets(unique_counter))
        select type(unique_tuplets); type is(lo_symlist_ucpair)
            do i=1,unique_counter
                unique_tuplets(i)=dumtup(i)
            enddo
        end select
    type is(lo_symlist_uctriplet)
        allocate(lo_symlist_uctriplet::unique_tuplets(unique_counter))
        select type(unique_tuplets); type is(lo_symlist_uctriplet)
            do i=1,unique_counter
                unique_tuplets(i)=dumtup(i)
            enddo
        end select
    type is(lo_symlist_ucquartet)
        allocate(lo_symlist_ucquartet::unique_tuplets(unique_counter))
        select type(unique_tuplets); type is(lo_symlist_ucquartet)
            do i=1,unique_counter
                unique_tuplets(i)=dumtup(i)
            enddo
        end select
    type is(lo_symlist_ucmagpair)
        allocate(lo_symlist_ucmagpair::unique_tuplets(unique_counter))
        select type(unique_tuplets); type is(lo_symlist_ucmagpair)
            do i=1,unique_counter
                unique_tuplets(i)=dumtup(i)
            enddo
        end select
    class default
        write(*,*) 'ERROR: This should not be able to happen, trying to return unique tuplets.'
        stop
    end select

    ! and a little cleanup
    lo_deallocate(dumtup)
    lo_deallocate(unhash)
end subroutine

!> given a unitcell pair, where do I find it in the supercell?
subroutine locate_pair_in_supercell(tr,sl,ss,ssa,ssi)
    !> triplet to locate
    type(lo_symlist_ucpair), intent(in) :: tr
    !> table and stuff
    type(lo_symlist), intent(in) :: sl
    !> unit and supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> index to atom in supercell
    integer, intent(out) :: ssa
    !> index to triplet in supercell
    integer, intent(out) :: ssi

    integer :: a1,i
    ssa=0
    ssi=0
    do a1=1,sl%nss
        ! same unitcell atom
        if ( ss%info%index_in_unitcell(a1) .ne. tr%i1 ) cycle
        ! look through the triplets
        do i=1,sl%ss(a1)%npair
            if ( sl%ss(a1)%pair(i)%ui1 .ne. tr%ui1 ) cycle
            if ( sl%ss(a1)%pair(i)%ui2 .ne. tr%ui2 ) cycle
            if ( lo_sqnorm(sl%ss(a1)%pair(i)%v-tr%v) .lt. sl%sqtol_cart ) then
                ssa=a1
                ssi=i
                return
            endif
        enddo
    enddo
    if ( ssa .eq. 0 ) then
        call lo_stop_gracefully(['Could not find unitcell pair in the supercell'],lo_exitcode_symmetry,__FILE__,__LINE__)
    endif
end subroutine

!> given a unitcell magnetic pair, where do I find it in the supercell?
subroutine locate_magpair_in_supercell(tr,sl,ss,ssa,ssi)
    !> triplet to locate
    type(lo_symlist_ucmagpair), intent(in) :: tr
    !> table and stuff
    type(lo_symlist), intent(in) :: sl
    !> unit and supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> index to atom in supercell
    integer, intent(out) :: ssa
    !> index to triplet in supercell
    integer, intent(out) :: ssi

    integer :: a1,i
    ssa=0
    ssi=0
    ia1l: do a1=1,sl%nss
        ! same unitcell atom
        if ( ss%info%index_in_unitcell(a1) .ne. tr%i1 ) cycle
        ! look through the triplets
        do i=1,sl%ss(a1)%nmagpair
            if ( sl%ss(a1)%magpair(i)%ui1 .ne. tr%ui1 ) cycle
            if ( sl%ss(a1)%magpair(i)%ui2 .ne. tr%ui2 ) cycle
            if ( lo_sqnorm(sl%ss(a1)%magpair(i)%v-tr%v) .lt. sl%sqtol_cart ) then
                ssa=a1
                ssi=i
                exit ia1l
            endif
        enddo
    enddo ia1l
    if ( ssa .eq. 0 ) then
        call lo_stop_gracefully(['Could not find unitcell pair in the supercell'],lo_exitcode_symmetry,__FILE__,__LINE__)
    endif
end subroutine

!> given a unitcell triplet, where do I find it in the supercell?
subroutine locate_triplet_in_supercell(tr,sl,ss,ssa,ssi)
    !> triplet to locate
    type(lo_symlist_uctriplet), intent(in) :: tr
    !> table and stuff
    type(lo_symlist), intent(in) :: sl
    !> unit and supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> index to atom in supercell
    integer, intent(out) :: ssa
    !> index to triplet in supercell
    integer, intent(out) :: ssi

    integer :: a1,i
    ssa=0
    ssi=0
    do a1=1,sl%nss
        ! same unitcell atom
        if ( ss%info%index_in_unitcell(a1) .ne. tr%i1 ) cycle
        ! look through the triplets
        do i=1,sl%ss(a1)%ntriplet
            if ( sl%ss(a1)%triplet(i)%ui1 .ne. tr%ui1 ) cycle
            if ( sl%ss(a1)%triplet(i)%ui2 .ne. tr%ui2 ) cycle
            if ( sl%ss(a1)%triplet(i)%ui3 .ne. tr%ui3 ) cycle
            if ( lo_sqnorm(sl%ss(a1)%triplet(i)%v2-tr%v2) .gt. sl%sqtol_cart ) cycle
            if ( lo_sqnorm(sl%ss(a1)%triplet(i)%v3-tr%v3) .lt. sl%sqtol_cart ) then
                ssa=a1
                ssi=i
                return
            endif
        enddo
    enddo
    if ( ssa .eq. 0 ) then
        write(*,*) 'ERROR: could not locate triplet in supercell'
        stop
    endif
end subroutine

!> given a unitcell triplet, where do I find it in the supercell?
subroutine locate_quartet_in_supercell(tr,sl,ss,ssa,ssi)
    !> triplet to locate
    type(lo_symlist_ucquartet), intent(in) :: tr
    !> table and stuff
    type(lo_symlist), intent(in) :: sl
    !> unit and supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> index to atom in supercell
    integer, intent(out) :: ssa
    !> index to triplet in supercell
    integer, intent(out) :: ssi

    integer :: a1,i
    ssa=0
    ssi=0
    do a1=1,sl%nss
        ! same unitcell atom
        if ( ss%info%index_in_unitcell(a1) .ne. tr%i1 ) cycle
        ! look through the triplets
        do i=1,sl%ss(a1)%nquartet
            if ( sl%ss(a1)%quartet(i)%ui1 .ne. tr%ui1 ) cycle
            if ( sl%ss(a1)%quartet(i)%ui2 .ne. tr%ui2 ) cycle
            if ( sl%ss(a1)%quartet(i)%ui3 .ne. tr%ui3 ) cycle
            if ( sl%ss(a1)%quartet(i)%ui4 .ne. tr%ui4 ) cycle
            if ( lo_sqnorm(sl%ss(a1)%quartet(i)%v2-tr%v2) .gt. sl%sqtol_cart ) cycle
            if ( lo_sqnorm(sl%ss(a1)%quartet(i)%v3-tr%v3) .gt. sl%sqtol_cart ) cycle
            if ( lo_sqnorm(sl%ss(a1)%quartet(i)%v4-tr%v4) .lt. sl%sqtol_cart ) then
                ssa=a1
                ssi=i
                return
            endif
        enddo
    enddo
    if ( ssa .eq. 0 ) then
        write(*,*) 'ERROR: could not locate quartet in supercell'
        write(*,*) real(tr%v1),tr%i1,tr%ui1
        write(*,*) real(tr%v2),tr%i2,tr%ui2
        write(*,*) real(tr%v3),tr%i3,tr%ui3
        write(*,*) real(tr%v4),tr%i4,tr%ui4
        stop
    endif
end subroutine

!> find a supercell triplet given the j1,j2,j3 indices
subroutine locate_triplet_in_supercell_via_indices(ind,sl,ssa,ssi)
    !> indices to atoms
    integer, dimension(3), intent(in) :: ind
    !> table and stuff
    type(lo_symlist), intent(in) :: sl
    !> index to atom in supercell
    integer, intent(out) :: ssa
    !> index to triplet in supercell
    integer, intent(out) :: ssi

    integer :: i

    ssa=ind(1)
    ssi=0
    do i=1,sl%ss(ssa)%ntriplet
        if ( sl%ss(ssa)%triplet(i)%j2 .ne. ind(2) ) cycle
        if ( sl%ss(ssa)%triplet(i)%j3 .eq. ind(3) ) then
            ssi=i
        endif
    enddo
    if ( ssi .eq. 0 ) then
        write(*,*) 'ERROR: could not locate triplet in supercell'
    endif
end subroutine

!> find a supercell triplet given the j1,j2,j3 indices
subroutine locate_quartet_in_supercell_via_indices(ind,sl,ssa,ssi)
    !> indices to atoms
    integer, dimension(4), intent(in) :: ind
    !> table and stuff
    type(lo_symlist), intent(in) :: sl
    !> index to atom in supercell
    integer, intent(out) :: ssa
    !> index to triplet in supercell
    integer, intent(out) :: ssi

    integer :: i
    ssa=ind(1)
    ssi=0
    do i=1,sl%ss(ssa)%nquartet
        if ( sl%ss(ssa)%quartet(i)%j2 .ne. ind(2) ) cycle
        if ( sl%ss(ssa)%quartet(i)%j3 .ne. ind(3) ) cycle
        if ( sl%ss(ssa)%quartet(i)%j4 .eq. ind(4) ) then
            ssi=i
        endif
    enddo
    if ( ssi .eq. 0 ) then
        write(*,*) 'ERROR: could not locate quartet in supercell via indices'
        stop
    endif
end subroutine

end submodule
