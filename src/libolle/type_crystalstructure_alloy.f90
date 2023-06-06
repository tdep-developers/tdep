#include "precompilerdefinitions"
submodule (type_crystalstructure) type_crystalstructure_alloy
implicit none
contains

!> creates the permutation that orders the SQS the same way as the reference
module subroutine alloy_site_permutation(ss,sqs,permutation)
    !> reference ideal structure
    class(lo_crystalstructure), intent(in) :: ss
    !> ideal SQS
    type(lo_crystalstructure), intent(in) :: sqs
    !> permutation, how to order the atoms in the sqs as the atoms in the supercell
    integer, dimension(:), intent(out) :: permutation


    init: block
        ! First some sanity checks
        if ( ss%info%alloy .eqv. .false. ) then
            call lo_stop_gracefully(['The reference structure needs to be an alloy.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( ss%na .ne. sqs%na ) then
            call lo_stop_gracefully(['Need the same number of atoms to generate permutation.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( norm2(ss%latticevectors-sqs%latticevectors) .gt. lo_tol ) then
            call lo_stop_gracefully(['The lattices need to match when generating permutations.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Add more if needed eventually.

    end block init

    ! Naive and slow way of matching. Plenty of ways to make it faster if I have to at some point.
    ! Ideas include, but are not limited to
    ! a) Sort both, then compare quickly. Annoying to sort vectors with PBC though.
    ! b) Verlet box for lookups, O(N) and sensibly fast I think. Maybe doublebox.
    ! Whatever. Just double-loop for now. Should not be time-sensitive in any way, but you never know.
    stupidandslow: block
        real(flyt), dimension(3) :: v0
        real(flyt) :: f0
        integer :: i,j,k
        logical, dimension(:), allocatable :: assigned

        allocate(assigned(ss%na))
        assigned=.false.
        permutation=0

        do i=1,ss%na
            k=0
            do j=1,sqs%na
                if ( assigned(j) ) cycle
                v0=ss%r(:,i)-sqs%r(:,j)
                v0=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
                f0=lo_sqnorm(v0)
                if ( f0 .lt. lo_sqtol ) then
                    k=j
                    exit
                endif
            enddo

            if ( k .eq. 0 ) then
                call lo_stop_gracefully(['Could not find permutation.'],&
                                        lo_exitcode_param,__FILE__,__LINE__)
            else
                permutation(i)=k
            endif
        enddo

        ! Check that nothing is unassigned
        if ( count(permutation==0) .ne. 0 ) then
            call lo_stop_gracefully(['Could not find permutation.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
    end block stupidandslow

end subroutine

!> Change the order of the atoms according to some permutation
module subroutine permute_positions(p,permutation,forward)
    !> structure
    class(lo_crystalstructure), intent(inout) :: p
    !> how to permute them
    integer, dimension(:), intent(in) :: permutation
    !> forward or backwards permutation
    logical, intent(in) :: forward

    if ( forward ) then
        p%atomic_number                     = p%atomic_number                   (permutation)
        p%species                           = p%species                         (permutation)
        p%r                                 = p%r                               (:,permutation)
        p%rcart                             = p%rcart                           (:,permutation)
        p%v                                 = p%v                               (:,permutation)
        p%f                                 = p%f                               (:,permutation)
        p%u                                 = p%u                               (:,permutation)
        p%mass                              = p%mass                            (permutation)
        p%invsqrtmass                       = p%invsqrtmass                     (permutation)
        p%inelastic_neutron_cross_section   = p%inelastic_neutron_cross_section (permutation)
        !p%flavor                            = p%flavor                          (permutation)
        if ( allocated(p%info%index_in_unitcell) ) then
            p%info%index_in_unitcell = p%info%index_in_unitcell (permutation)
        endif
        if ( allocated(p%info%cellindex) ) then
            p%info%cellindex = p%info%cellindex(:,permutation)
        endif
    else
        p%atomic_number                   (permutation) = p%atomic_number
        p%species                         (permutation) = p%species
        p%r                               (:,permutation) = p%r
        p%rcart                           (:,permutation) = p%rcart
        p%v                               (:,permutation) = p%v
        p%f                               (:,permutation) = p%f
        p%u                               (:,permutation) = p%u
        p%mass                            (permutation) = p%mass
        p%invsqrtmass                     (permutation) = p%invsqrtmass
        p%inelastic_neutron_cross_section (permutation) = p%inelastic_neutron_cross_section
        !p%flavor                          (permutation) = p%flavor
        if ( allocated(p%info%index_in_unitcell) ) then
            p%info%index_in_unitcell(permutation) = p%info%index_in_unitcell
        endif
        if ( allocated(p%info%cellindex) ) then
            p%info%cellindex(:,permutation) = p%info%cellindex
        endif
    endif

end subroutine

end submodule
