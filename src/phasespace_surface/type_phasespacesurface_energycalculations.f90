
subroutine energy_values_on_vertices(qp1,qp,uc,fc,energy_plus,energy_minus,mem)
    !> the q-point in question
    type(lo_qpoint), intent(in) :: qp1 
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> energy values
    real(flyt), dimension(:,:,:,:), intent(out) :: energy_plus,energy_minus
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !
    type(lo_qpoint) :: qp2,qp3
    type(lo_phonon_dispersions_qpoint) :: drp1,drp2,drp3
    integer :: i,nb,b1,b2,b3
    !
    nb=fc%na*3
    energy_plus=0.0_flyt
    energy_minus=0.0_flyt
    ! get the starting point
    call drp1%generate(fc,uc,mem,qp1)
    do i=1,qp%n_full_point
        qp2%r=qp%ap(i)%r
        qp2%n_invariant_operation=0
        qp3%r=qp2%r+qp1%r
        qp3%n_invariant_operation=0
        call drp2%generate(fc,uc,mem,qp2)
        call drp3%generate(fc,uc,mem,qp3)
        do b1=1,nb
        do b2=1,nb
        do b3=1,nb
            energy_plus(i,b1,b2,b3)=drp1%omega(b1)+drp2%omega(b2)-drp3%omega(b3)
            energy_minus(i,b1,b2,b3)=drp1%omega(b1)-drp2%omega(b2)-drp3%omega(b3)
        enddo
        enddo
        enddo
    enddo
end subroutine


