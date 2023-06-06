#include "precompilerdefinitions"
module type_graph_distancetable
!! Handle distance tables defined by number of NN jumps
!use constants
!use gottochblandat
!use geometryfunctions
use type_distancetable, only: lo_distancetable
!use type_voronoi

implicit none
private
public :: lo_graph_distancetable

!> Specialized distance table that uses number of NN jumps as the cutoff
type, extends(lo_distancetable) :: lo_graph_distancetable
    !> number of jumps that define the table
    integer :: njump
    contains
!        procedure :: generate_graphtable
end type

contains

! !> Get a distance table with number of NN jumps as the cutoff
! subroutine generate_graphtable(dt,particles,basis,njump,cutoff,coordination)
!     !> the distance table
!     class(lo_graph_distancetable), intent(out) :: dt
!     !> the basis the cell
!     real(flyt), dimension(3,3), intent(in) :: basis
!     !> the list of particles
!     real(flyt), dimension(:,:), intent(in) :: particles
!     !> number of jumps
!     integer, intent(in) :: njump
!     !> cutoff used when generating
!     real(flyt), intent(in) :: cutoff
!     !> pre-defined coordination numbers
!     integer, dimension(:), intent(in), optional :: coordination
! 
! 
! 
! !    !> the supercell dimensions of this cell
! !    integer, dimension(3), intent(in) :: dimensions
! !    !> the basis to build the voronoi cell from
! !    real(flyt), dimension(3,3), intent(in), optional :: vorobasis
! !    !> how many atoms are contained in this basis
! !    integer, intent(in), optional :: vorona
! !!    !> the actual voronoi polyhedron
! !!    type(lo_polyhedron), intent(out), optional :: polyhedron
! !    !
! !    type(lo_voronoi_diagram) :: voro
! !    type(lo_distancetable) :: odt
! !    real(flyt), dimension(:,:,:), allocatable :: nnvecs
! !    real(flyt), dimension(:), allocatable :: dr,du,dv
! !    real(flyt), dimension(3) :: v0,v1,v2
! !    real(flyt) :: f0,nndist,nntol,sqrc,t0
! !    integer, dimension(:), allocatable :: crd
! !    integer :: a1,a2,i,j,k,l,ndist,i1,i2,i3,ctr,current_jump
! 
!     dt%np=0
! 
! !    ! First thing to do is to generate the normal distance table. I will use the longest
! !    ! possible cutoff in the box, I think.
! !
! !    t0=walltime()
! !    ! initialize all the jump thingies
! !    call odt%generate(particles,basis,cutoff*2)
! !    do a1=1,odt%np
! !        lo_allocate(odt%particle(a1)%jumps( odt%particle(a1)%n ))
! !        odt%particle(a1)%jumps=-1
! !        odt%particle(a1)%jumps(1)=0
! !    enddo
! !    ! Build the voronoi diagram
! !    call voro%generate(particles,basis,cutoff,verbosity=1)
! !    l=0
! !
! !do a1=1,voro%n
! !    write(*,*) a1,voro%cell(a1)%nfaces
! !    ! start traversing the diagram
! !    do i1=1,voro%cell(a1)%nfaces
! !        v0=voro%cell(a1)%face(i1)%neighbourvector
! !write(*,*) i1,v0
! !        ! locate this in the distance table, define it as the first coordination shell
! !        do i2=1,odt%particle(a1)%n
! !            if ( lo_sqnorm(v0-odt%particle(a1)%v(:,i2)) .lt. lo_sqtol ) then
! !                odt%particle(a1)%jumps(i2)=1
! !                exit
! !                !@TODO Here I can add a penalty for layered materials, depending on the vector!
! !            endif
! !        enddo
! !    enddo
! !enddo
! 
!     stop
!     
! !    ! If I have pre-defined coordination number, just use those
! !    if ( present(coordination) ) then
! !        ! sanity check
! !        if ( coordination(1) .lt. 0 ) then
! !            write(*,*) 'Can not work with negative coordination numbers'
! !            write(*,*) 'coord:',tochar(coordination)
! !            stop
! !        endif
! !        ! set the N first to be coordination 1
! !        do a1=1,odt%np
! !            do i=1,coordination(a1)
! !                odt%particle(a1)%jumps(i+1)=1
! !            enddo
! !            ! Check that this makes sense
! !            f0=odt%particle(a1)%d( coordination(a1)+1 )
! !            l=-1
! !            do i=1,odt%particle(a1)%n
! !                if ( odt%particle(a1)%d(i)-f0 .lt. lo_tol ) l=l+1
! !            enddo
! !            if ( l .ne. coordination(a1) ) then
! !                write(*,*) 'Could not define coordination ',tochar(coordination(a1)),' for atom ',tochar(a1)
! !                write(*,*) 'Found ',tochar(l),' neighbours instead.'
! !                stop
! !            endif
! !        enddo
! !    else
! !        ! If not, I have to guess somehow. Not as robust as I would like.
! !        do a1=1,odt%np
! !            ! what is the nearest neighbour distance from this atom?
! !            nndist=odt%particle(a1)%d(2)
! !            ! pick a slightly longer NN distance to compare with
! !            nndist=nndist+nndist*0.1_flyt
! !            ! some random number that says that if we are within this of the 
! !            ! NN distance, it belongs to the first shell.
! !            nntol=nndist*0.2_flyt
! !
! !            l=0
! !            do i=1,odt%particle(a1)%n
! !                if ( abs(odt%particle(a1)%d(i)-nndist) .lt. nntol ) then
! !                    l=l+1
! !                    odt%particle(a1)%jumps(i)=1
! !                endif
! !            enddo
! !        enddo
! !    endif
! !    write(*,*) 'mellantid: ',tochar(walltime()-t0),'s'
! !    t0=walltime()
! !
! !    ! store shorthand for the coordination of each atom
! !    lo_allocate(crd(odt%np))
! !    crd=0
! !    do a1=1,odt%np
! !        l=0
! !        do i=1,odt%particle(a1)%n
! !            if ( odt%particle(a1)%jumps(i) .eq. 1 ) l=l+1
! !        enddo
! !        crd(a1)=l
! !        write(*,*) 'atom',a1,'coord',l
! !    enddo
! !
! !    ! Ok, now it's time to spread this out.
! !    sqrc=odt%cutoff**2
! !    current_jump=1
! !    do
! !        ctr=0
! !        do a1=1,odt%np
! !            ! check if we are done with this atom
! !            if ( minval(odt%particle(a1)%jumps) .ge. 0 ) then 
! !                cycle
! !            endif
! !            ! nope, not done, spread things out.
! !            do i1=1,odt%particle(a1)%n
! !                ! I want to jump from the current jump
! !                if ( odt%particle(a1)%jumps(i1) .ne. current_jump ) cycle
! !                ! so, from this neighbour I need to go to all neighbours that are 1 jump away from this atom.
! !                ! this is the vector to neighbour i1
! !                v1=odt%particle(a1)%v(:,i1) ! ( a2-a1 )
! !                a2=odt%particle(a1)%ind(i1)
! !                ! so from the neighbouring atom, walk to all neighbours within the coordination
! !                do i2=2,odt%particle(a2)%n
! !                    ! has to be a nearest neighbour jump
! !                    if ( odt%particle(a2)%jumps(i2) .ne. 1 ) cycle
! !                    v2=odt%particle(a2)%v(:,i2) ! (a3-a2)
! !                    v0=v2+v1                    ! (a3-a1)
! !                    ! it might be an impossible vector
! !                    if ( lo_sqnorm(v0) .gt. sqrc ) cycle
! !                    ! now locate this vector from atom 1
! !                    do i3=1,odt%particle(a1)%n
! !                        ! skip if it's already determined
! !                        if ( odt%particle(a1)%jumps(i3) .ge. 0 ) cycle
! !                        ! check if I can find it
! !                        if ( lo_sqnorm(v0-odt%particle(a1)%v(:,i3)) .lt. lo_sqtol ) then
! !                            ! got it!
! !                            odt%particle(a1)%jumps(i3)=current_jump+1
! !                            ctr=ctr+1
! !                            exit
! !                        endif
! !                    enddo
! !                enddo
! !            enddo
! !        enddo
! !        ! check if we are done
! !        if ( ctr .eq. 0 ) then
! !            ! this means nothing was assigned. Probably done
! !            exit
! !        else
! !            ! if not, go one jump further
! !            current_jump=current_jump+1
! !        endif
! !        ! extra sanity check
! !        if ( current_jump .gt. 5*njump ) exit
! !    enddo
! !
! !    ! Ok now most things should be assigned. Store this as a normal distance table.
! !
! !do a1=1,odt%np
! !    write(*,*) 'atom',a1,'nn',odt%particle(a1)%n
! !    do i=1,odt%particle(a1)%n
! !        write(*,*) i,odt%particle(a1)%jumps(i),odt%particle(a1)%d(i)
! !    enddo
! !enddo
! !    write(*,*) 'mellantid: ',tochar(walltime()-t0),'s'
! !    t0=walltime()
! 
! end subroutine

end module

