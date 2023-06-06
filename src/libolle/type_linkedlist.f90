#include "precompilerdefinitions"
module type_linkedlist
!! Handle linked lists with unlimited polymorphism. Requires some fiddling for each kind of object you want to use it for.
use konstanter, only: lo_status

implicit none
private
public :: lo_linked_list

!> A link in the chain
type link
    !> polymorphic object to store data in. Avoiding a pointer here to prevent memory leaks.
    class(*), allocatable :: obj
    !> move to the next object
    type(link), pointer :: next => null()
    !> move to the previous
    type(link), pointer :: previous => null()
end type link

!> Linked list that can hold anything.
type lo_linked_list
    !> First link in list
    type(link), pointer :: first => null()
    !> Last link in list
    type(link), pointer :: last => null()
    !> The active link in the list
    type(link), pointer :: current => null()
    contains
        !> Add an extra object to the list
        procedure :: add => add_polymorphic_object_to_list
        !> Count how many objects there are
        procedure :: countitems => count_items_in_linked_list
        !> Remove the active object
        procedure :: remove => remove_active_link
        !> Replace the active object
        procedure :: replace => replace_active_link
        !> Destroy the whole list
        procedure :: destroy
end type

contains

!> Count how many links I have in the list
function count_items_in_linked_list(list) result(l)
    !> list
    class(lo_linked_list), intent(inout) :: list
    integer :: l

    if ( associated(list%first) ) then
        ! start from the beginning
        list%current => list%first
        l=0
        ! iterate over the whole list
        do
            l=l+1
            if ( associated(list%current%next) ) then
                list%current=>list%current%next
            else
                exit
            endif
        enddo
        ! reset the list
        list%current=>list%first
    else
        ! if there is no first item, there is no list.
        l=0
    endif
end function

!> Add something to the end of the linked list.
subroutine add_polymorphic_object_to_list(list,obj)
    !> The list
    class(lo_linked_list), intent(inout) :: list
    class(*), intent(in) :: obj
    !
    type(link), pointer :: newlink
    !
    ! allocate space for the new link
    lo_allocate(newlink)
    ! make space for the information
    allocate(newlink%obj,source=obj,stat=lo_status)
    if ( lo_status .ne. 0 ) then
        write(*,*) 'failed allocating link'
        stop
    endif
    ! fix the iterators
    if ( associated(list%first) .eqv. .false. ) then
        ! It's the first link
        list%first=>newlink
        list%last=>newlink
        list%current=>newlink
    else
        ! there is at least one link already.
        newlink%previous=>list%last
        ! add it as the last link
        list%last%next=>newlink
        list%last=>newlink
    endif
end subroutine

!> Replace the object in the current link with a new one.
subroutine replace_active_link(list,object)
    !> The list
    class(lo_linked_list), intent(inout) :: list
    !> THe object
    class(*), intent(in) :: object
    ! A lot easier than my first 20 attempts:
    if ( associated(list%current) ) then
        if ( allocated(list%current%obj) ) lo_deallocate(list%current%obj)
        allocate(list%current%obj,source=object,stat=lo_status)
        if ( lo_status .ne. 0 ) then
            write(*,*) 'failed replacing link'
        endif
    endif
end subroutine

!> Remove the current element from the linked list
subroutine remove_active_link(list,wheredidIgo)
    !> The list
    class(lo_linked_list), intent(inout) :: list
    integer, intent(out) :: wheredidIgo
    !
    type(link), pointer :: l,pre,nex

    ! Keep track of the one I want to delete
    l=>null()
    if ( associated(list%current) ) l=>list%current

    ! There are a few possibilities:
    if ( associated(list%current) .eqv. .false. ) then
        ! no list left
        wheredidIgo=0
        return
    elseif ( (associated(list%current%next) .eqv. .false.) .and. &
         (associated(list%current%previous) .eqv. .false.) ) then
        ! the final link, remove list completely
        list%current => null()
        list%first => null()
        list%last => null()
        wheredidIgo=0
    elseif ( associated(list%current%next) .and. associated(list%current%previous) ) then
        ! I'm in the middle, remove element and join around it.
        ! store endpoints
        pre=>list%current%previous
        nex=>list%current%next
        ! remove the current
        list%current=>null()
        ! join the loose ends
        pre%next=>nex
        nex%previous=>pre
        ! put us at the previous one ( or the next? )
        list%current=>pre
        ! I stepped backwards
        wheredidIgo=-1
    elseif ( associated(list%current%next) .eqv. .false. ) then
        ! I'm at the end
        ! kill last element
        list%last => list%current%previous
        list%last%next=> null()
        ! point to the new last?
        list%current => list%last
        ! I stepped backwards
        wheredidIgo=-1
    elseif ( associated(list%current%previous) .eqv. .false. ) then
        ! I'm at the beginning
        ! Destroy first element
        list%first => list%current%next
        list%first%previous => null()
        ! point to the new first?
        list%current => list%first
        ! Now I stepped forward
        wheredidIgo=1
    endif
    ! Deallocate memory
    if ( associated(l) ) then
        if ( allocated(l%obj) ) lo_deallocate(l%obj)
        lo_deallocate(l)
    endif
end subroutine

!> Deallocate a linked list
subroutine destroy(list)
    !> The list
    class(lo_linked_list), intent(inout) :: list
    integer :: i
    do
        if ( associated(list%first) ) then
            list%current=>list%first
            call list%remove(i)
        else
            ! now it's destroyed
            exit
        endif
    enddo
end subroutine

end module
