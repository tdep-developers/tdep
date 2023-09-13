module pairmapping
!!
!! Build a minimal version of the symmetry list
!!
use konstanter, only: r8, lo_huge, lo_hugeint, lo_sqtol, lo_tol, lo_bohr_to_A, lo_exitcode_symmetry
use gottochblandat, only: walltime, lo_determ, lo_sqnorm, tochar, lo_cross, lo_chop
use type_crystalstructure, only: lo_crystalstructure
use type_symmetrylist, only: lo_symlist
use type_mdsim, only: lo_mdsim
use type_distancetable, only: lo_distancetable
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully

implicit none

private
public :: lo_pairmapping

!> representation of all the atoms in a single coordination shell
type lo_pm_atom_shell
    !> what atoms does it go between?
    integer :: i1 = -lo_hugeint, i2 = -lo_hugeint
    !> what is the prototype pair vector, ideally?
    real(r8), dimension(3) :: prototype_vector = -lo_huge
    !> the coordinate transformation that makes it flat and nice
    real(r8), dimension(3, 3) :: coordinate_transformation = -lo_huge
end type

!> representation of the pairs originating from this atom
type lo_pm_atom
    !> each atom will have a certain number of coordination shells
    integer :: nshell = -lo_hugeint
    !> but it will also have a certain number of neighbours
    integer :: npair = -lo_hugeint
    !> and the actual shells
    type(lo_pm_atom_shell), dimension(:), allocatable :: shell
    !> how many atoms in the unitcell correspond to this unique type
    integer :: n_uc_atoms = -lo_hugeint
    !> which unitcell atoms correspond to this unique type
    integer, dimension(:), allocatable :: unitcell_indices
end type

type lo_pm_ssatom_pair
    !> which two atoms in the supercell are involved
    integer :: i1 = -lo_hugeint, i2 = -lo_hugeint
    !> which shell for that atom
    integer :: shell = -lo_hugeint
    !> and the matrix you should transform with to take it to the prototype
    real(r8), dimension(3, 3) :: m = -lo_huge
    !> the bond distance
    real(r8) :: rad = -lo_huge
    !> the lattice vector for this pair
    real(r8), dimension(3) :: lv = -lo_huge
end type

!> just some shorthand to make the mapping less of a pain
type lo_pm_ssatom
    !> npair
    integer :: npair = -lo_hugeint
    !> info about each pair
    type(lo_pm_ssatom_pair), dimension(:), allocatable :: pair
    !> which unique atom is it?
    integer :: unique_atom = -lo_hugeint
end type

!> keep track of a coordination shell
type lo_pairmapping_shell
    !> what are the atom indices for this shell?
    integer :: i1 = -lo_hugeint
    integer :: i2 = -lo_hugeint
    !> a neat label for this shell?
    character(len=50) :: label
    !> what is the prototype vector for this shell?
    real(r8), dimension(3) :: r
    !> how many supercell pairs are in this shell
    integer :: n_ss_pair = -lo_hugeint
    !> pair indices
    integer, dimension(:, :), allocatable :: pairind
    !> baseline pair  vectors
    real(r8), dimension(:, :), allocatable :: pairvec
    !> pair operation
    integer, dimension(:), allocatable :: pairop
    !> which species pair is this shell?
    integer :: index_species_pair=-lo_hugeint
end type

!> structure that handles how all pairs are related by symmetry
type lo_pairmapping
    !> number of bins in the histogram
    integer :: nh = -lo_hugeint
    !> number of unique atoms
    integer :: na = -lo_hugeint
    !> number of atoms in the supercell
    integer :: nass = -lo_hugeint
    !> a characteristic bond length
    real(r8) :: nndist = -lo_huge
    !> the cutoff for the neighbour list
    real(r8) :: nncutoff = -lo_huge
    !> one distribution per unique atom
    type(lo_pm_atom), dimension(:), allocatable :: atom
    !> a neat list for all the pairs in the supercell
    type(lo_pm_ssatom), dimension(:), allocatable :: ssatom

    !> how many symmetry-distinct shells are there
    integer :: n_shell = -lo_hugeint
    !> coordination shell
    type(lo_pairmapping_shell), dimension(:), allocatable :: sh

    !> number of species shells
    integer :: n_species_pair =-lo_hugeint
    !> labels for species pairs
    character(len=100), dimension(:), allocatable :: species_pair_label
contains
    !> Figure out the symmetry mapping
    procedure :: setup_symmetry
end type

contains

!> initialize things
subroutine setup_symmetry(pm, sl, uc, ss, mw, mem, verbosity)
    !> pair mapping thing
    class(lo_pairmapping), intent(out) :: pm
    !> symmetry list thing
    type(lo_interaction_tensors), intent(in) :: sl
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    character(len=1000), dimension(:), allocatable :: shell_labels
    integer :: ish, ipair, iop, i,j, s1, s2

    ! Start copying things somehow:
    pm%n_shell = sl%n_fc_pair_shell
    allocate (pm%sh(pm%n_shell))
    ! Populate the shells
    do ish = 1, pm%n_shell
        ! index to prototype pair
        ipair = sl%fc_pair_shell(ish)%unfold_index(1)
        iop = sl%fc_pair_shell(ish)%unfold_operation(1)
        if (iop .ne. 1) call lo_stop_gracefully(['Thinking is hard'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)

        ! store which atoms are involved
        pm%sh(ish)%i1 = sl%uc_pair(ipair)%i1
        pm%sh(ish)%i2 = sl%uc_pair(ipair)%i2
        pm%sh(ish)%r = sl%uc_pair(ipair)%v

        ! neat label for the shell
        s1 = pm%sh(ish)%i1
        s2 = pm%sh(ish)%i2
        s1 = uc%species(s1)
        s2 = uc%species(s2)
        if (s1 .gt. s2) then
            pm%sh(ish)%label = trim(uc%atomic_symbol(s2))//'-'//trim(uc%atomic_symbol(s1))
        else
            pm%sh(ish)%label = trim(uc%atomic_symbol(s1))//'-'//trim(uc%atomic_symbol(s2))
        end if

        ! Set things to nothing
        pm%sh(ish)%n_ss_pair = 0
    end do

    ! Count pairs per shell
    do ipair = 1, sl%n_ss_pair
        ish = sl%ss_pair(ipair)%unique
        pm%sh(ish)%n_ss_pair = pm%sh(ish)%n_ss_pair + 1
    end do
    ! Make space
    do ish = 1, pm%n_shell
        allocate (pm%sh(ish)%pairvec(3, pm%sh(ish)%n_ss_pair))
        allocate (pm%sh(ish)%pairind(2, pm%sh(ish)%n_ss_pair))
        allocate (pm%sh(ish)%pairop(pm%sh(ish)%n_ss_pair))
        pm%sh(ish)%pairvec = 0.0_r8
        pm%sh(ish)%pairind = 0
        pm%sh(ish)%pairop = 0
        pm%sh(ish)%n_ss_pair = 0
    end do
    ! populate
    do ipair = 1, sl%n_ss_pair
        ish = sl%ss_pair(ipair)%unique
        pm%sh(ish)%n_ss_pair = pm%sh(ish)%n_ss_pair + 1
        i = pm%sh(ish)%n_ss_pair
        pm%sh(ish)%pairvec(:, i) = sl%ss_pair(ipair)%v
        pm%sh(ish)%pairop(i) = sl%ss_pair(ipair)%operation
        pm%sh(ish)%pairind(:, i) = [sl%ss_pair(ipair)%j1, sl%ss_pair(ipair)%j2]
    end do

    if (verbosity .gt. 0) then
        write (*, *) ''
        write (*, *) 'Identified ', tochar(pm%n_shell), ' coordination shells:'
        do ish = 1, pm%n_shell
            s1 = pm%sh(ish)%i1
            s2 = pm%sh(ish)%i2
            s1 = uc%species(s1)
            s2 = uc%species(s2)
            write (*, *) ish, norm2(pm%sh(ish)%r)*lo_bohr_to_A, pm%sh(ish)%label
        end do
    end if

    ! We might want to bin things into pairs per pair of species, so to say. It makes
    ! sense to say which species pair each symmetry unique pair is.
    i=0
    do s1=1,uc%nelements
    do s2=s1,uc%nelements
        i=i+1
    enddo
    enddo
    pm%n_species_pair=i
    allocate(pm%species_pair_label(pm%n_species_pair))
    i=0
    do s1=1,uc%nelements
    do s2=s1,uc%nelements
        i=i+1
        pm%species_pair_label(i) = trim(uc%atomic_symbol(s1))//'-'//trim(uc%atomic_symbol(s2))
    enddo
    enddo

    do ish=1,pm%n_shell
        j=-1
        do i=1,pm%n_species_pair
            if ( trim(pm%sh(ish)%label) .eq. trim(pm%species_pair_label(i)) ) then
                j=i
                exit
            endif
        enddo
        if ( j .gt. 0 ) then
            pm%sh(ish)%index_species_pair=j
        else
            call lo_stop_gracefully(['Could not find pair'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
    enddo



    ! type(lo_distancetable) :: dt
    ! real(r8), dimension(3,3) :: m0
    ! real(r8), dimension(3) :: v0,v1,v2
    ! integer, dimension(:,:), allocatable :: sh_per_un
    ! integer, dimension(:), allocatable :: di1
    ! integer :: a1,sh,uca,ucb,una
    ! integer :: i,j,k,o
    ! integer :: ii,jj,kk
    ! real(r8) :: f0,f1,t0
    !
    ! ! Need a reductionlist, from the pairs in the unique representation, to a certain
    ! ! coordination shell, along with the operation that takes it there. First make space.
    ! pm%na=sl%nsingletshells
    ! lo_allocate(pm%atom(pm%na))
    ! lo_allocate(di1(sl%npairshells))
    ! lo_allocate(sh_per_un(sl%npairshells,pm%na))
    ! sh_per_un=0
    ! do a1=1,pm%na
    !     ! I need to know how many coordination shells originate from this atom
    !     uca=sl%singletshell(a1)%protatom
    !     di1=0.0_r8
    !     kk=0
    !     ploop1: do i=1,sl%uc(uca)%npair
    !         j=sl%uc(uca)%pair(i)%unique_shell
    !         do k=1,kk
    !             if ( j .eq. di1(k) ) cycle ploop1
    !         enddo
    !         kk=kk+1
    !         di1(kk)=j
    !     enddo ploop1
    !     ! store these shells
    !     sh_per_un(1:kk,a1)=di1(1:kk)
    !     ! make some space
    !     pm%atom(a1)%nshell=kk
    !     pm%atom(a1)%npair=sl%uc(uca)%npair
    !     lo_allocate(pm%atom(a1)%shell( pm%atom(a1)%nshell ))
    !     pm%atom(a1)%n_uc_atoms=0
    ! enddo
    ! ! A little indexing to find how atoms are mapped here
    ! do i=1,sl%nuc
    !     a1=sl%uc(i)%singlet%unique_shell
    !     pm%atom(a1)%n_uc_atoms=pm%atom(a1)%n_uc_atoms+1
    ! enddo
    ! do a1=1,pm%na
    !     lo_allocate(pm%atom(a1)%unitcell_indices( pm%atom(a1)%n_uc_atoms ))
    !     pm%atom(a1)%n_uc_atoms=0
    ! enddo
    ! do i=1,sl%nuc
    !     a1=sl%uc(i)%singlet%unique_shell
    !     pm%atom(a1)%n_uc_atoms=pm%atom(a1)%n_uc_atoms+1
    !     pm%atom(a1)%unitcell_indices( pm%atom(a1)%n_uc_atoms )=i
    ! enddo
    !
    ! ! Get a prototype vector for each shell. I think I want one close to 100
    ! v0=[1004.0_r8,103.0_r8,12.0_r8]
    ! do a1=1,pm%na
    !     uca=sl%singletshell(a1)%protatom ! unitcell atom
    !     do ii=1,pm%atom(a1)%nshell
    !         sh=sh_per_un(ii,a1) ! the unique shell
    !         f0=lo_huge
    !         do i=1,sl%uc(uca)%npair
    !             if ( sl%uc(uca)%pair(i)%unique_shell .ne. sh ) cycle
    !             f1=lo_sqnorm(v0-sl%uc(uca)%pair(i)%v)
    !             if ( f1 .lt. f0 ) then
    !                 ! I have a prototype
    !                 pm%atom(a1)%shell(ii)%i1=sl%uc(uca)%pair(i)%i1
    !                 pm%atom(a1)%shell(ii)%i2=sl%uc(uca)%pair(i)%i2
    !                 pm%atom(a1)%shell(ii)%prototype_vector=sl%uc(uca)%pair(i)%v
    !                 f0=f1
    !             endif
    !         enddo
    !     enddo
    ! enddo
    !
    ! ! Print these prototypes:
    ! write(*,*) ' '
    ! do a1=1,pm%na
    !     write(*,*) 'From atom '//tochar(a1)//' we have '//tochar(pm%atom(a1)%nshell)//' coordination shells'
    !     do sh=1,pm%atom(a1)%nshell
    !         ii=pm%atom(a1)%shell(sh)%i1
    !         jj=pm%atom(a1)%shell(sh)%i2
    !         write(*,"(1X,I4,1X,3(2X,F10.6),' norm: ',F10.5,' atoms: ',A)") sh,lo_chop( pm%atom(a1)%shell(sh)%prototype_vector*lo_bohr_to_A ,lo_sqtol),&
    !              norm2(pm%atom(a1)%shell(sh)%prototype_vector)*lo_bohr_to_A,&
    !              trim(uc%atomic_symbol(uc%species(ii)))//'-'//trim(uc%atomic_symbol(uc%species(jj)))
    !     enddo
    ! enddo
    !
    ! ! Figure out the mapping from the supercell to these shells
    ! pm%nass=sl%nss
    ! lo_allocate(pm%ssatom( pm%nass ))
    ! do a1=1,pm%nass
    !     ! some basic stuff
    !     pm%ssatom(a1)%npair=sl%ss(a1)%npair
    !     lo_allocate(pm%ssatom(a1)%pair( pm%ssatom(a1)%npair ))
    !     pm%ssatom(a1)%unique_atom=sl%ss(a1)%singlet%unique_shell
    !     ! unique atom
    !     una=pm%ssatom(a1)%unique_atom
    !     uca=sl%singletshell(una)%protatom       ! prototype unitcell atom
    !     ucb=ss%info%index_in_unitcell(a1)       ! index in unitcell
    !     do i=1,sl%ss(a1)%npair
    !         ! first store the indices
    !         pm%ssatom(a1)%pair(i)%i1=sl%ss(a1)%pair(i)%j1
    !         pm%ssatom(a1)%pair(i)%i2=sl%ss(a1)%pair(i)%j2
    !         pm%ssatom(a1)%pair(i)%rad=norm2(sl%ss(a1)%pair(i)%v)
    !         ! Which shell should it correspond to?
    !         sh=0
    !         ii=sl%ss(a1)%pair(i)%unique_shell
    !         do j=1,pm%atom(una)%nshell
    !             if ( sh_per_un(j,una) .eq. ii ) then
    !                 sh=j
    !                 exit
    !             endif
    !         enddo
    !         if ( sh .eq. 0 ) then
    !             write(*,*) 'Failed getting the proper shell for atom ',tochar(a1),' in the supercell'
    !         else
    !             pm%ssatom(a1)%pair(i)%shell=sh
    !         endif
    !
    !         ! now find the proper operation.
    !         v0=pm%atom(una)%shell(sh)%prototype_vector   ! vector I want to match
    !         v1=sl%ss(a1)%pair(i)%v                      ! current vector
    !         ii=0
    !         syml: do o=1,uc%sym%n
    !             if ( uc%sym%op(o)%fmap(ucb) .ne. uca ) cycle
    !             v2=matmul(uc%sym%op(o)%m,v1)
    !             if ( lo_sqnorm(v0-v2) .lt. lo_sqtol ) then
    !                 ii=o
    !                 exit syml
    !             endif
    !         enddo syml
    !
    !         if ( ii .eq. 0 ) then
    !             write(*,*) 'Failed getting the proper operation for atom ',tochar(a1),' in the supercell'
    !             stop
    !         else
    !             ! store the proper operation
    !             pm%ssatom(a1)%pair(i)%m=uc%sym%op(ii)%m
    !         endif
    !
    !         ! I also want the lattice-vector associated with this pair
    !         v0=sl%ss(a1)%pair(i)%v
    !         v0=v0-ss%rcart(:,sl%ss(a1)%pair(i)%j2)+ss%rcart(:,a1)
    !         v0=int(anint(ss%cartesian_to_fractional(v0,pbc=.false.)))
    !         v0=ss%fractional_to_cartesian(v0)
    !         pm%ssatom(a1)%pair(i)%lv=v0
    !     enddo
    ! enddo
    !
    ! ! Check that it works
    ! do a1=1,pm%nass
    !     ! unique atom
    !     ii=pm%ssatom(a1)%unique_atom
    !     do i=1,pm%ssatom(a1)%npair
    !         ! the original vector
    !         v0=sl%ss(a1)%pair(i)%v
    !         ! shell
    !         jj=pm%ssatom(a1)%pair(i)%shell
    !         ! the prototype vector
    !         v1=pm%atom( ii )%shell( jj )%prototype_vector
    !         ! check difference
    !         m0=pm%ssatom(a1)%pair(i)%m
    !         !
    !         f0=lo_sqnorm( matmul(m0,v0)-v1 )
    !         if ( f0 .gt. lo_sqtol ) then
    !             write(*,*) 'TROUBLE!!! BAD MAPPING!!!'
    !             stop
    !         endif
    !     enddo
    ! enddo
    !
    !  ! Some auxiliary values that might be useful
    ! pm%nndist=lo_huge
    ! pm%nncutoff=0.0_r8
    ! do a1=1,pm%na
    ! do i=1,pm%atom(a1)%nshell
    !     f0=norm2(pm%atom(a1)%shell(i)%prototype_vector)
    !     if ( f0 .gt. lo_tol ) then
    !         pm%nndist=min(pm%nndist,f0)
    !         pm%nncutoff=max(pm%nncutoff,f0)
    !     endif
    ! enddo
    ! enddo
    ! ! and a little margin for safety
    ! pm%nncutoff=pm%nncutoff+10*lo_tol
    !
    ! ! get a local coordinate system
    ! t0=walltime()
    ! call dt%generate(uc%r,uc%latticevectors,pm%nncutoff*3,verbosity=0)
    ! write(*,*) 'got distabletable',walltime()-t0
    !
    ! ! Find orthogonal vectors
    ! do a1=1,pm%na
    ! do sh=1,pm%atom(a1)%nshell
    !     ! the given vector
    !     v0=pm%atom(a1)%shell(sh)%prototype_vector
    !     !
    !     if ( lo_sqnorm(v0) .lt. lo_sqtol ) then
    !         m0(1,:)=[1, 0, 0]
    !         m0(2,:)=[0, 1, 0]
    !         m0(3,:)=[0, 0, 1]
    !     else
    !         ! find something orthogonal
    !         uca=sl%uc(a1)%singlet%unique_shell
    !         j=0
    !         do i=2,dt%particle(uca)%n
    !             v1=dt%particle(uca)%v(:,i)
    !             if ( abs(dot_product(v0,v1)) .lt. lo_tol ) then
    !                 j=1
    !                 exit
    !             endif
    !         enddo
    !         !
    !         if ( j .gt. 0 ) then
    !             v2=lo_cross(v0,v1)
    !             v0=v0/norm2(v0)
    !             v1=v1/norm2(v1)
    !             v2=v2/norm2(v2)
    !             ! try one version first
    !             m0(1,:)=v1
    !             m0(2,:)=v2
    !             m0(3,:)=v0
    !             if ( lo_determ(m0) .lt. 0.0_r8 ) then
    !                 m0(1,:)=v2
    !                 m0(2,:)=v1
    !                 m0(3,:)=v0
    !             endif
    !         else
    !             m0(1,:)=[1, 0, 0]
    !             m0(2,:)=[0, 1, 0]
    !             m0(3,:)=[0, 0, 1]
    !         endif
    !     endif
    !     ! store the transformation
    !     pm%atom(a1)%shell(sh)%coordinate_transformation=m0
    ! enddo
    ! enddo
end subroutine

end module
