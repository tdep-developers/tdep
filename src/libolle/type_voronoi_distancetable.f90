#include "precompilerdefinitions"
module type_voronoi_distancetable
!! Handle distance tables with non-spherical cutoffs
use konstanter, only: flyt,lo_tol,lo_sqtol
use gottochblandat, only: lo_determ,lo_sqnorm
use geometryfunctions, only: lo_linesegment
use type_distancetable, only: lo_distancetable
use type_voronoi, only: lo_voronoi_diagram

implicit none
private
public :: lo_voronoi_distancetable

!> Specialized distance table that uses a super-Wigner-Seitz cell as the cutoff.
type, extends(lo_distancetable) :: lo_voronoi_distancetable
    !> the lattice vectors used to generate the voronoi cell
    real(flyt), dimension(3,3) :: voronoibasis    
    contains        
        procedure :: generate_wzcutoff 
end type

contains

!> Get a distance table with voronoi cell as the cutoff
subroutine generate_wzcutoff(dt,particles,basis,dimensions,vorobasis,vorona) 
    !> the distance table
    class(lo_voronoi_distancetable), intent(out) :: dt
    !> the basis the cell
    real(flyt), dimension(3,3), intent(in) :: basis
    !> the list of particles
    real(flyt), dimension(:,:), intent(in) :: particles
    !> the supercell dimensions of this cell
    integer, dimension(3), intent(in) :: dimensions
    !> the basis to build the voronoi cell from
    real(flyt), dimension(3,3), intent(in), optional :: vorobasis
    !> how many atoms are contained in this basis
    integer, intent(in), optional :: vorona
    !
    type(lo_distancetable) :: odt
    type(lo_voronoi_diagram) :: voro
    type(lo_linesegment) :: segment
    real(flyt), dimension(size(particles,1),size(particles,2)) :: rcart
    real(flyt), dimension(3,1) :: r0
    real(flyt), dimension(3) :: v0
    integer :: i,j,k,l,i1,i2,expected_pairs
    
    ! First, generate the Voronoi diagram:
    
    ! Create basis for the supercell
    if ( present(vorobasis) ) then
        do i=1,3
            dt%voronoibasis(:,i)=vorobasis(:,i)*dimensions(i)
        enddo
    else
        do i=1,3
            dt%voronoibasis(:,i)=basis(:,i)*dimensions(i)
        enddo
    endif
    ! Figure out how many pairs are expected
    if ( present(vorona) ) then
        if ( present(vorobasis) ) then
            expected_pairs=int(anint(abs( lo_determ(dt%voronoibasis)/lo_determ(vorobasis) )))*vorona
        else
            expected_pairs=int(anint(abs( lo_determ(dt%voronoibasis)/lo_determ(basis) )))*vorona
        endif
    else
        expected_pairs=int(anint(abs( lo_determ(dt%voronoibasis)/lo_determ(basis) )))*size(particles,2)
    endif
    ! put an atom at the origin
    r0=0.0_flyt
    ! and get the actual super-WZ-cell
    call voro%generate(r0,dt%voronoibasis,verbosity=0)
    ! Get the spherical distance table. Thise will include all the neighbours, and then we will remove points
    call odt%generate(particles,basis,voro%cell(1)%rmax*1.02_flyt,verbosity=0)
    
    ! Truncate the distance table
    do i=1,size(particles,2)
        rcart(:,i)=matmul(basis,particles(:,i))
    enddo
    dt%cutoff=voro%cell(1)%rmin
    dt%np=size(particles,2) 
    lo_allocate(dt%particle(dt%np))
    do i=1,dt%np
        ! Count atoms within the WZ-cell
        l=0
        do j=1,odt%particle(i)%n
            ! test if it's inside the cell, lo_tol=1E-5 is the tolerance if somethings is considered to be inside
            if ( voro%cell(1)%polyhedron%is_point_inside( odt%particle(i)%v(:,j) , lo_tol ) ) then
                l=l+1
            endif
        enddo
        
        ! Make some space, now that I know how many there are
        dt%particle(i)%n=l
        lo_allocate(dt%particle(i)%v(3,l))
        lo_allocate(dt%particle(i)%lv(3,l))
        lo_allocate(dt%particle(i)%ind(l))
        lo_allocate(dt%particle(i)%d(l))
        lo_allocate(dt%particle(i)%weight(l))
        ! store the data
        l=0
        do j=1,odt%particle(i)%n
            if ( voro%cell(1)%polyhedron%is_point_inside( odt%particle(i)%v(:,j) , lo_tol ) ) then
                l=l+1
                dt%particle(i)%v(:,l)=odt%particle(i)%v(:,j)
                dt%particle(i)%ind(l)=odt%particle(i)%ind(j)
                dt%particle(i)%d(l)=odt%particle(i)%d(j)
                ! the lattice vector
                dt%particle(i)%lv(:,l)=dt%particle(i)%v(:,l)
                v0=rcart(:,i)-rcart(:,dt%particle(i)%ind(l))
                dt%particle(i)%lv(:,l)=dt%particle(i)%lv(:,l)+v0
            endif
        enddo
        
        ! Calculate the weights.
        dt%particle(i)%weight=0.0_flyt
        particleloop: do j=1,dt%particle(i)%n
            if ( dt%particle(i)%d(j) .lt. voro%cell(1)%rmin ) then
                ! it's safely inside
                dt%particle(i)%weight(j)=1.0_flyt
            else
                ! check with nodes
                do k=1,voro%cell(1)%nnodes
                    if ( lo_sqnorm( dt%particle(i)%v(:,j)-voro%cell(1)%r(:,k) ) .lt. lo_sqtol ) then
                        dt%particle(i)%weight(j)=1.0_flyt/(voro%cell(1)%node(k)%nsharedcells+1)
                        cycle particleloop
                    endif
                enddo
                ! check with edges
                do k=1,voro%cell(1)%nedges
                    i1=voro%cell(1)%edge(k)%i1
                    i2=voro%cell(1)%edge(k)%i2
                    call segment%generate( voro%cell(1)%r(:, i1 ), voro%cell(1)%r(:, i2 ))
                    if ( segment%distance_to_point(dt%particle(i)%v(:,j)) .lt. lo_tol ) then
                        dt%particle(i)%weight(j)=1.0_flyt/(voro%cell(1)%edge(k)%nshared+1)
                        cycle particleloop
                    endif
                enddo
                ! check with faces
                do k=1,voro%cell(1)%nfaces
                    if ( abs(voro%cell(1)%face(k)%plane%distance_to_point(dt%particle(i)%v(:,j))) .lt. lo_tol ) then
                        dt%particle(i)%weight(j)=0.5_flyt
                        cycle particleloop
                    endif
                enddo
                ! and, finally, if it's nothing it just has weight 1
                dt%particle(i)%weight(j)=1.0_flyt
            endif
        enddo particleloop
        
        ! Make a sanity check
        if ( abs(expected_pairs*1.0_flyt-sum(dt%particle(i)%weight)) .gt. lo_tol ) then
            write(*,*) 'WARNING'
            write(*,*) 'improper number of atoms within the cutoff, considering the weights.'
        endif
        !
    enddo
end subroutine

end module

