#include "precompilerdefinitions"
submodule (type_symmetrylist) type_symmetrylist_helper
implicit none
contains

!> expand operations to include transpositions as well
module subroutine expandoperations(sl)
    !> symmetry table
    type(lo_symlist), intent(inout) :: sl

    real(flyt), dimension(:,:,:), allocatable :: trm_pair,trm_triplet,trm_quartet,so81,so27,so9
    integer, dimension(:,:), allocatable :: prm_pair,prm_triplet,prm_quartet
    integer :: l,t,o
    integer :: ntr,nop

    ! Make space for the transpositions. I had to make these allocatable since I ran inte weird weird bugs
    lo_allocate(trm_pair(9,9,2))
    lo_allocate(trm_triplet(27,27,6))
    lo_allocate(trm_quartet(81,81,24))
    lo_allocate(prm_pair(2,2))
    lo_allocate(prm_triplet(3,6))
    lo_allocate(prm_quartet(4,24))
    trm_pair=0.0_flyt
    trm_triplet=0.0_flyt
    trm_quartet=0.0_flyt
    prm_pair=0
    prm_triplet=0
    prm_quartet=0

    ! Get all the transpositions
    call lo_return_tensor_transpositions(trm_pair,trm_triplet,trm_quartet,prm_pair,prm_triplet,prm_quartet)
    ! start by figuring out if I should consider spacegroup/transpositions or not.
    if ( sl%spacegroup ) then
        nop=sl%sym%n
    else
        nop=1 ! only identity
    endif

    ! First build the pair transpositions/operations
    if ( sl%transposition ) then
        ntr=size(prm_pair,2)
    else
        ntr=1 ! only identity
    endif

    ! expand the pair operations
    if ( sl%secondorder ) then
        lo_allocate(so9(9,9,nop))
        so9=0.0_flyt
        do o=1,nop
            so9(:,:,o)=expandoperation_pair(sl%sym%op(o)%m)
        enddo
        ! count and make space
        sl%npairop=nop*ntr
        lo_allocate(sl%pairop(sl%npairop))
        ! populate
        l=0
        do t=1,ntr
        do o=1,nop
            l=l+1
            if ( t .eq. 1 ) then
                sl%pairop(l)%m3=sl%sym%op(o)%m
                sl%pairop(l)%im3=sl%sym%op(o)%im
            else
                sl%pairop(l)%m3=-sl%sym%op(o)%m
                sl%pairop(l)%im3=-sl%sym%op(o)%im
            endif
            call lo_gemm(so9(:,:,o),trm_pair(:,:,t),sl%pairop(l)%sotr)
            sl%pairop(l)%perm=prm_pair(:,t)
            sl%pairop(l)%permind=t
            sl%pairop(l)%opind=o
            ! Pretty sure this is how it should be done.
            sl%pairop(l)%msotr=anint(lo_determ(sl%sym%op(o)%m))*sl%pairop(l)%sotr
            !sl%pairop(l)%msotr=sl%pairop(l)%sotr
        enddo
        enddo
    else
        sl%npairop=0
    endif

    ! Now the same thing for the triplets
    if ( sl%thirdorder ) then
        ! And for quartets, eventually
        lo_allocate(so27(27,27,nop))
        so27=0.0_flyt
        do o=1,nop
            so27(:,:,o)=expandoperation_triplet(sl%sym%op(o)%m)
        enddo
        if ( sl%transposition ) then
            ntr=size(prm_triplet,2)
        else
            ntr=1
        endif
        sl%ntripletop=nop*ntr
        lo_allocate(sl%tripletop(sl%ntripletop))
        ! populate
        l=0
        do t=1,ntr
        do o=1,nop
            l=l+1
            sl%tripletop(l)%m3=sl%sym%op(o)%m
            sl%tripletop(l)%im3=sl%sym%op(o)%im
            call lo_gemm(so27(:,:,o),trm_triplet(:,:,t),sl%tripletop(l)%sotr)
            sl%tripletop(l)%perm=prm_triplet(:,t)
            sl%tripletop(l)%permind=t
            sl%tripletop(l)%opind=o
        enddo
        enddo
    else
        sl%ntripletop=0
    endif

    if ( sl%fourthorder ) then
        ! And for quartets, eventually
        lo_allocate(so81(81,81,nop))
        so81=0.0_flyt
        do o=1,nop
            so81(:,:,o)=expandoperation_quartet(sl%sym%op(o)%m)
        enddo
        if ( sl%transposition ) then
            ntr=size(prm_quartet,2)
        else
            ntr=1
        endif
        sl%nquartetop=ntr*nop
        lo_allocate(sl%quartetop(sl%nquartetop))
        ! populate
        l=0
        do t=1,ntr
        do o=1,nop
            l=l+1
            sl%quartetop(l)%m3=sl%sym%op(o)%m
            sl%quartetop(l)%im3=sl%sym%op(o)%im
            call lo_gemm(so81(:,:,o),trm_quartet(:,:,t),sl%quartetop(l)%sotr)
            sl%quartetop(l)%perm=prm_quartet(:,t)
            sl%quartetop(l)%permind=t
            sl%quartetop(l)%opind=o
        enddo
        enddo
        lo_deallocate(so81)
    else
        sl%nquartetop=0
    endif
    lo_deallocate(trm_pair)
    lo_deallocate(trm_triplet)
    lo_deallocate(trm_quartet)
    lo_deallocate(prm_pair)
    lo_deallocate(prm_triplet)
    lo_deallocate(prm_quartet)
end subroutine

!> check that cutoffs are sane
module subroutine check_cutoffs(uc,ss,rc2,rc3,rc4,wraparound)
    !> the crystal structures
    type(lo_crystalstructure), intent(in) :: uc,ss
    !> the cutoffs
    real(flyt), intent(inout) :: rc2,rc3,rc4
    !> are wraparound distances allowed?
    logical, intent(in) :: wraparound
    !
    type(lo_distancetable) :: dt
    integer :: i,j,k
    real(flyt) :: f0

    ! Check min and max
    rc2=max(ss%mincutoff(),rc2)
    if ( wraparound .eqv. .false. ) rc2=min(ss%maxcutoff(),rc2)
    if ( rc3 .gt. 0.0_flyt ) then
        rc3=max(ss%mincutoff(),rc3)
        if ( wraparound .eqv. .false. ) rc3=min(ss%maxcutoff(),rc3)
    endif
    if ( rc4 .gt. 0.0_flyt ) then
        rc4=max(ss%mincutoff(),rc4)
        if ( wraparound .eqv. .false. ) rc4=min(ss%maxcutoff(),rc4)
    endif

    ! Get a distance table
    f0=max(max(rc2,rc3),rc4)+0.1_flyt
    call dt%generate(uc%r,uc%latticevectors,f0,verbosity=0)
    ! Check if the cutoffs coincide with pair distances, that causes instabilities.
    do
        k=0
        do i=1,dt%np
        do j=1,dt%particle(i)%n
            if ( abs(rc2-dt%particle(i)%d(j)) .lt. 3*lo_tol ) k=k+1
        enddo
        enddo
        if ( k .eq. 0 ) then
            exit
        else
            write(*,*) 'WARNING: second order cutoff coincides with a coordination shell, increasing it:'
            write(*,*) '    from:',rc2
            rc2=rc2+4.0_flyt*lo_tol
            write(*,*) '      to:',rc2
        endif
    enddo
    !
    if ( rc3 .gt. 0.0_flyt ) then
        do
            k=0
            do i=1,dt%np
            do j=1,dt%particle(i)%n
                if ( abs(rc3-dt%particle(i)%d(j)) .lt. 3*lo_tol ) k=k+1
            enddo
            enddo
            !
            if ( k .eq. 0 ) then
                exit
            else
                write(*,*) 'WARNING: third order cutoff coincides with a coordination shell, increasing it:'
                write(*,*) '    from:',rc3
                rc3=rc3+4.0_flyt*lo_tol
                write(*,*) '      to:',rc3
            endif
        enddo
    endif
    !
    if ( rc4 .gt. 0.0_flyt ) then
        do
            k=0
            do i=1,dt%np
            do j=1,dt%particle(i)%n
                if ( abs(rc4-dt%particle(i)%d(j)) .lt. 3*lo_tol ) k=k+1
            enddo
            enddo
            !
            if ( k .eq. 0 ) then
                exit
            else
                write(*,*) 'WARNING: fourth order cutoff coincides with a coordination shell, increasing it:'
                write(*,*) '    from:',rc4
                rc4=rc4+4.0_flyt*lo_tol
                write(*,*) '      to:',rc4
            endif
        enddo
    endif
end subroutine

!> hashing function to create some sort of identifier for the tuplets
module function hashtuplet(t) result(hash)
    class(lo_symlist_tuplet), intent(in) :: t
    real(flyt) :: hash
    !
    real(flyt), dimension(3) :: v0,v1,v2,v3,v4
    hash=0.0_flyt
    select type(t)
    class is(lo_symlist_pair)
        ! pairs, easy
        hash=hash+t%ui1+t%ui2
        hash=hash+lo_sqnorm(t%v)
    class is(lo_symlist_triplet)
        ! triplets, also not tricky
        hash=hash+t%ui1+t%ui2+t%ui3
        v1=t%v2
        v2=t%v3
        v3=t%v3-t%v2
        hash=hash+lo_sqnorm(v1)+lo_sqnorm(v2)+lo_sqnorm(v3)
        hash=hash+dot_product(v1,v2)**2+dot_product(v1,v3)**2+dot_product(v2,v3)**2
    class is(lo_symlist_quartet)
        ! quartets, a little annoying
        hash=hash+t%ui1+t%ui2+t%ui3+t%ui4
        ! center of mass
        v0=(t%v1+t%v2+t%v3+t%v4)*0.25_flyt
        ! distances to tetrahedron corners from the center of mass
        v1=t%v1-v0
        v2=t%v2-v0
        v3=t%v3-v0
        v4=t%v4-v0
        hash=hash+lo_sqnorm(v1)+lo_sqnorm(v2)+lo_sqnorm(v3)+lo_sqnorm(v4)
    class is(lo_symlist_magpair)
        ! pairs, easy
        hash=hash+t%ui1+t%ui2
        hash=hash+lo_sqnorm(t%v)
    end select
end function

! Below are not exposed

!> expand a 3x3 symmetry operation to a 9x9 operation, to operate on flattened tensors
pure function expandoperation_pair(o) result(bigo)
    !> the original 3x3 operation
    real(flyt), dimension(3,3), intent(in) :: o
    !> the resulting 9x9 operation
    real(flyt), dimension(9,9) :: bigo
    !
    integer :: i,j,ii,jj,k,l
    !
    bigo=0.0_flyt
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
    bigo=lo_chop(bigo,1E-11_flyt)
end function

!> expand a 3x3 symmetry operation to a 27x27 operation, to operate on vector format of tensors
module pure function expandoperation_triplet(o) result(bigo)
    !> the original 3x3 operation
    real(flyt), dimension(3,3), intent(in) :: o
    !> the resulting 27x27 operation
    real(flyt), dimension(27,27) :: bigo
    !
    integer :: i,ii,iii,j,jj,jjj,l,m
    !
    bigo=0.0_flyt
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
    bigo=lo_chop(bigo,1E-11_flyt)
end function

!> expand a 3x3 symmetry operation to a 81x81 operation, to operate on vector format of tensors
pure function expandoperation_quartet(o) result(bigo)
    !> the original 3x3 operation
    real(flyt), dimension(3,3), intent(in) :: o
    !> the resulting 27x27 operation
    real(flyt), dimension(81,81) :: bigo

    integer :: i,ii,iii,iiii,j,jj,jjj,jjjj,l,m
    bigo=0.0_flyt
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
    bigo=lo_chop(bigo,1E-11_flyt)
end function

! #ifdef AGRESSIVE_SANITY
! !> Test that everything is as I think it is.
! subroutine sanitytest(sl,sh,uc,ss)
!     type(lo_symlist), intent(inout) :: sl
!     type(lo_symtabhelper), intent(in) :: sh
!     type(lo_crystalstructure), intent(in) :: uc
!     type(lo_crystalstructure), intent(in) :: ss
!
!     ! Start with really basic things, that never ever should be wrong:
!     testbasics: block
!         logical :: err
!         integer :: a1,a2
!
!         err=.false.
!         ! check number of atoms
!         if ( sl%nuc .ne. uc%na ) err=.true.
!         if ( sl%nss .ne. ss%na ) err=.true.
!         if ( sl%secondorder ) then
!             ! check number of pairs
!             do a1=1,sl%nss
!                 a2=ss%info%index_in_unitcell(a1)
!                 if ( sl%uc(a2)%npair .ne. sl%ss(a1)%npair ) err=.true.
!             enddo
!         endif
!         if ( sl%thirdorder ) then
!             ! check number of triplets
!             do a1=1,sl%nss
!                 a2=ss%info%index_in_unitcell(a1)
!                 if ( sl%uc(a2)%ntriplet .ne. sl%ss(a1)%ntriplet ) err=.true.
!             enddo
!         endif
!         if ( err ) then
!             write(*,*) 'ERROR:'
!             write(*,*) '    symmetry table failed basic sanity test.'
!         endif
!     end block testbasics
!
!     ! test that all the tuplets are connected properly
!     ! I also make sure not to use the same transformation functions
!     ! as before, can perhaps catch some bugs that way.
!     testconnections: block
!         logical :: err
!         integer :: a1,op,i,s
!         integer, dimension(:), allocatable :: prota,proti
!         real(flyt), dimension(3,3) :: m0
!         real(flyt), dimension(3) :: v0,v1,v2
!         real(flyt), dimension(9) :: fpm1,fpm2,fpm3
!
!         err=.false.
!         if ( sl%secondorder ) then
!             ! locate the prototype in the full thing
!             lo_allocate(prota(sl%npairshells))
!             lo_allocate(proti(sl%npairshells))
!             prota=0
!             proti=0
!             shl1: do s=1,sl%npairshells
!                 do a1=1,sl%nuc
!                 do i=1,sl%uc(a1)%npair
!                     if ( a1 .ne. sl%pairshell(s)%protpair%i1 ) cycle
!                     if ( sl%uc(a1)%pair(i)%i2 .ne. sl%pairshell(s)%protpair%i2 ) cycle
!                     if ( lo_sqnorm(sl%uc(a1)%pair(i)%v-sl%pairshell(s)%protpair%v) .lt. lo_sqtol ) then
!                         prota(s)=a1
!                         proti(s)=i
!                         cycle shl1
!                     endif
!                 enddo
!                 enddo
!                 if ( prota(s) .eq. 0 .or. proti(s) .eq. 0 ) then
!                     write(*,*) 'ERROR:'
!                     write(*,*) '    could not locate the prototype pair from the shells in the unitcell.'
!                     write(*,*) '    I was looking for this vector:'
!                     write(*,*) sl%pairshell(s)%protpair%v
!                     write(*,*) 'from atom ',tochar(sl%pairshell(s)%protpair%i1),' to ',tochar(sl%pairshell(s)%protpair%i2)
!                     stop
!                 endif
!             enddo shl1
!             ! test operations from shells to unitcell
!             do a1=1,sl%nuc
!             do i=1,sl%uc(a1)%npair
!                 ! just test the vectors first
!                 s=sl%uc(a1)%pair(i)%unique_shell
!                 op=sl%uc(a1)%pair(i)%operation_from_shell_to_pair
!                 v0=sl%pairshell(s)%protpair%v
!                 v1=sl%uc(a1)%pair(i)%v
!                 v2=matmul(sl%pairop(op)%m3,v0)
!                 if ( lo_sqnorm(v2-v1) .gt. sl%sqtol_cart ) then
!                     write(*,*) 'bad unitcell connection in terms of the vectors'
!                     err=.true.
!                 endif
!                 ! then less fake (but still fake) matrices
!                 m0=fake_eam_forceconstant(sh%dtuc,uc,prota(s),proti(s))
!                 fpm1=lo_flattentensor( m0 )
!                 m0=fake_eam_forceconstant(sh%dtuc,uc,a1,i)
!                 fpm2=lo_flattentensor( m0 )
!                 fpm3=matmul(sl%pairop(op)%sotr,fpm1)
!                 if ( sum(abs(fpm3-fpm2))/sum(abs(fpm2)) .gt. lo_tol ) then
!                     write(*,*) 'bad unitcell non-diagonal tensor connection'
!                     err=.true.
!                 endif
!             enddo
!             enddo
!             ! test operations from shells to supercell
!             do a1=1,sl%nss
!             do i=1,sl%ss(a1)%npair
!                 s=sl%ss(a1)%pair(i)%unique_shell
!                 op=sl%ss(a1)%pair(i)%operation_from_shell_to_pair
!                 v0=sl%pairshell(s)%protpair%v
!                 v1=sl%ss(a1)%pair(i)%v
!                 v2=matmul(sl%pairop(op)%m3,v0)
!                 if ( lo_sqnorm(v2-v1) .gt. sl%sqtol_cart ) then
!                     write(*,*) 'bad supercell connection'
!                     err=.true.
!                 endif
!             enddo
!             enddo
!             lo_deallocate(prota)
!             lo_deallocate(proti)
!         endif
!
!         if ( err ) then
!             stop
!         else
!             if ( sl%verbosity .gt. 1 ) write(*,*) '... passed basic sanity tests'
!         endif
!     end block testconnections
!
! end subroutine
! #endif

! #ifdef AGRESSIVE_SANITY
! !> embedded atom pair forceconstant, but not really
! function fake_eam_forceconstant(dt,uc,atom,pair) result(m)
!     type(lo_distancetable), intent(in) :: dt
!     type(lo_crystalstructure), intent(in) :: uc
!     integer, intent(in) :: atom
!     integer, intent(in) :: pair
!     real(flyt), dimension(3,3) :: m
!
!     integer :: i,j,k,ii,jj,a1,a2
!     real(flyt), dimension(:), allocatable :: rhobar
!     real(flyt), dimension(3,3) :: m0
!     real(flyt), dimension(3) :: rij,rjk,rik,v0,v1
!     real(flyt) :: f0,d
!
!     ! make space for temporary stuff
!     lo_allocate(rhobar(uc%na))
!     rhobar=0.0_flyt
!
!     ! Start by calculating rho
!     do i=1,uc%na
!         f0=0.0_flyt
!         do j=1,dt%particle(i)%n
!             d=dt%particle(i)%d(j)
!             if ( d .lt. lo_tol ) cycle
!             ! I add the masses here to have different interactions between different atoms.
!             f0=f0+embed(d)*sqrt(uc%mass(i))*sqrt(uc%mass(dt%particle(i)%ind(j)))
!         enddo
!         rhobar(i)=embed(f0)
!     enddo
!
!     m0=0.0_flyt
!     a1=atom
!     a2=dt%particle(a1)%ind(pair)
!     rij=dt%particle(a1)%v(:,pair)
!
!     ! loop over all neighbours of all atoms to get the embedding part
!     do ii=1,dt%particle(a1)%n
!     do jj=1,dt%particle(a2)%n
!         ! no self-terms
!         if ( dt%particle(a1)%d(ii) .lt. lo_tol ) cycle
!         if ( dt%particle(a2)%d(jj) .lt. lo_tol ) cycle
!         ! I need to find some sort of overlap here. It's supposed to be a sum of k, where the
!         ! k index is all atoms.
!         v0=dt%particle(a1)%v(:,ii)
!         v1=dt%particle(a2)%v(:,jj)+rij
!         if ( lo_sqnorm(v0-v1) .gt. lo_sqtol ) cycle
!         k=dt%particle(a1)%ind(ii)
!         ! ok, now I guess I have an atom that is within the cutoff from both.
!         rik=-v0
!         rjk=rij-v0
!         f0=sqrt(uc%mass(a1)*uc%mass(a2)*uc%mass(k)**2)
!         m0=m0+f0*embedpp(rhobar(k))*embedp(norm2(rik))*embedp(norm2(rjk))*lo_outerproduct(rik,rjk)/(norm2(rik)*norm2(rjk))
!     enddo
!     enddo
!     ! Add up contributions to the forceconstant
!     m=lo_chop( m0+pairpotfc(rij)*uc%mass(a1)*uc%mass(a2) ,lo_sqtol)
!
!     contains
!
!     function embed(x) result(y)
!         real(flyt), intent(in) :: x
!         real(flyt) :: y
!         y=sqrt(x)
!     end function
!     function embedp(x) result(y)
!         real(flyt), intent(in) :: x
!         real(flyt) :: y
!         y=0.5_flyt/sqrt(x)
!     end function
!     function embedpp(x) result(y)
!         real(flyt), intent(in) :: x
!         real(flyt) :: y
!         y=-0.25_flyt/(sqrt(x)*x)
!     end function
!     function pairpotfc(r) result(m)
!         real(flyt), dimension(3), intent(in) :: r
!         real(flyt), dimension(3,3) :: m
!
!         real(flyt) :: r2x,r2y,r2z
!         if ( lo_sqnorm(r) .lt. lo_sqtol ) then
!             m=0.0_flyt
!         else
!             r2x=r(1)
!             r2y=r(2)
!             r2z=r(3)
!             m(1,1)=(-16*r2x**2 )/(r2x**2 + r2y**2 + r2z**2)**3 + 4/(r2x**2 + r2y**2 + r2z**2)**2
!             m(1,2)=(-16*r2x*r2y)/(r2x**2 + r2y**2 + r2z**2)**3
!             m(1,3)=(-16*r2x*r2z)/(r2x**2 + r2y**2 + r2z**2)**3
!             m(2,1)=(-16*r2x*r2y)/(r2x**2 + r2y**2 + r2z**2)**3
!             m(2,2)=(-16*r2y**2 )/(r2x**2 + r2y**2 + r2z**2)**3 + 4/(r2x**2 + r2y**2 + r2z**2)**2
!             m(2,3)=(-16*r2y*r2z)/(r2x**2 + r2y**2 + r2z**2)**3
!             m(3,1)=(-16*r2x*r2z)/(r2x**2 + r2y**2 + r2z**2)**3
!             m(3,2)=(-16*r2y*r2z)/(r2x**2 + r2y**2 + r2z**2)**3
!             m(3,3)=(-16*r2z**2 )/(r2x**2 + r2y**2 + r2z**2)**3 + 4/(r2x**2 + r2y**2 + r2z**2)**2
!         endif
!     end function
! end function
! #endif

end submodule
