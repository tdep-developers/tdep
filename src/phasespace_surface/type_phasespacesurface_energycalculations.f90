
subroutine energy_values_on_vertices(qp1,qp,uc,fc,energy_plus,energy_minus)
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
    !
    type(lo_qpoint) :: qp2,qp3
    type(lo_phonon_dispersions_qpoint) :: drp1,drp2,drp3
    integer :: i,nb,b1,b2,b3
    !
    nb=fc%na*3
    energy_plus=0.0_flyt
    energy_minus=0.0_flyt
    ! get the starting point
    call drp1%generate(fc,uc,qp1)
    do i=1,qp%nq_tot
        qp2%w=qp%ap(i)%w
        qp2%noperations=0
        qp3%w=qp2%w+qp1%w
        qp3%noperations=0
        call drp2%generate(fc,uc,qp2)
        call drp3%generate(fc,uc,qp3)
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


