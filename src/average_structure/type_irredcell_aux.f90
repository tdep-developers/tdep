
!> Get irreducible representation from something
subroutine get_internal_irrep_from_unitcell(ic, r)
    !> irreducible representation
    class(lo_irredcell), intent(inout) :: ic
    !> unitcell positions in fractional coordinates
    real(flyt), dimension(:, :), intent(in) :: r

    real(flyt), dimension(:), allocatable :: flatu
    real(flyt), dimension(3) :: v0
    integer :: i

    lo_allocate(flatu(size(r)))
    flatu = 0.0_flyt
    do i = 1, size(r, 2)
        v0 = lo_clean_fractional_coordinates(r(:, i) - ic%reference_positions(:, i) + 0.5_flyt) - 0.5_flyt
        flatu((i - 1)*3 + 1:i*3) = v0
    end do
    call lo_linear_least_squares(ic%coeffM_internal, flatu, ic%x_internal)
end subroutine

!> Get irreducible representation
subroutine get_internal_irrep_from_supercell(ic, uc, ss, r)
    !> irreducible representation
    class(lo_irredcell), intent(inout) :: ic
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> unitcell positions in fractional coordinates
    real(flyt), dimension(:, :), intent(in) :: r

    real(flyt), dimension(:, :), allocatable :: coeffM
    real(flyt), dimension(:), allocatable :: flatu
    real(flyt), dimension(3, 3) :: m0
    real(flyt), dimension(3) :: v0
    integer :: i, j

    lo_allocate(flatu(size(r)))
    lo_allocate(coeffM(size(r), ic%nx_internal))
    flatu = 0.0_flyt
    coeffM = 0.0_flyt
    m0 = matmul(uc%inv_latticevectors, ss%latticevectors)

    do i = 1, size(r, 2)
        j = ss%info%index_in_unitcell(i)
        v0 = matmul(m0, r(:, i)) - ss%info%cellindex(:, i) + 1 - ic%reference_positions(:, j)
        v0 = lo_clean_fractional_coordinates(v0 + 0.5_flyt) - 0.5_flyt
        flatu((i - 1)*3 + 1:i*3) = v0
        coeffM((i - 1)*3 + 1:i*3, :) = ic%coeffM_internal((j - 1)*3 + 1:j*3, :)
    end do
    call lo_linear_least_squares(coeffM, flatu, ic%x_internal)
end subroutine

!!> Finds the rotation matrix such that lattice'=rotation*lattice represents a diagnoal distortion of ref_lattice.
!subroutine lo_find_rotation_that_makes_strain_diagonal( ref_lattice,lattice,rotation,strain,guess )
!    !> reference lattice
!    real(flyt), dimension(3,3), intent(in) :: ref_lattice
!    !> current lattice
!    real(flyt), dimension(3,3), intent(in) :: lattice
!    !> rotation matrix
!    real(flyt), dimension(3,3), intent(out) :: rotation
!    !> strain matrix
!    real(flyt), dimension(3,3), intent(out) :: strain
!    !> guess for the symmetric strain?
!    real(flyt), dimension(3,3), intent(in), optional :: guess
!
!    integer, parameter :: nouter=1000
!    integer, parameter :: ninner=1000
!    integer, parameter :: maxctr=20
!    real(flyt), parameter :: defstep=1E-5_flyt
!    real(flyt), parameter :: errtol=1E-20_flyt
!
!    real(flyt), dimension(6,6,6) :: coeffM
!    real(flyt), dimension(6) :: bvec,eta0,eta1,grad
!    real(flyt), dimension(3,3) :: m0
!    real(flyt) :: e0,e1,f0,step
!    integer :: initer,outiter,ctr,i
!
!    ! So, what do I solve here. The idea is that any transformation of a lattice takes
!    ! a,b,c,alpha,beta,gamma to new values. I need six numbers to define that transformation
!    ! uniquely. But you often come across these transformations as 3x3-matrices, so that they
!    ! might contain an extra rotation. This routine finds the diagonal form of a general
!    ! transformation between lattice vectors, since that can always be done.
!    !
!    ! So, in practice, I write the transformed lattice vectors
!    !
!    ! a'=eta*a
!    ! b'=eta*b
!    ! c'=eta*c
!    !
!    ! And find the symmetric matrix eta that produces the correct values for
!    !
!    !  a.a b.b c.c a.b a.c b.c
!    !
!    ! when applied to the original lattice. I end up with a bunch of matrix
!    ! equations on the form X^T A X = b, that I solve with a steepest descent.
!    ! Seems to be a robust procedure.
!
!    ! Get the coefficient matrices, A in XAX=b
!    call cfvars(coeffM(:,:,1),coeffM(:,:,2),coeffM(:,:,3),coeffM(:,:,4),coeffM(:,:,5),coeffM(:,:,6),ref_lattice)
!    ! Calculate the target quantities. All combinations of scalar products of
!    ! the lattice vectors, b in XAX=b
!    bvec(1)=dot_product(lattice(:,1),lattice(:,1))
!    bvec(2)=dot_product(lattice(:,2),lattice(:,2))
!    bvec(3)=dot_product(lattice(:,3),lattice(:,3))
!    bvec(4)=dot_product(lattice(:,1),lattice(:,2))
!    bvec(5)=dot_product(lattice(:,2),lattice(:,3))
!    bvec(6)=dot_product(lattice(:,3),lattice(:,1))
!    ! Do we have a guess for the strain?
!    if ( present(guess) ) then
!        eta0(1)=guess(1,1)
!        eta0(2)=guess(2,2)
!        eta0(3)=guess(3,3)
!        eta0(4)=guess(1,2)
!        eta0(5)=guess(1,3)
!        eta0(6)=guess(2,3)
!    else
!        ! Start from something not completely crazy
!        m0=matmul(ref_lattice,lo_invert3x3matrix(lattice))
!        eta0(1)=m0(1,1)
!        eta0(2)=m0(2,2)
!        eta0(3)=m0(3,3)
!        eta0(4)=m0(1,2)
!        eta0(5)=m0(1,3)
!        eta0(6)=m0(2,3)
!    endif
!    eta1=0.0_flyt
!    grad=0.0_flyt
!
!    outloop: do outiter=1,nouter
!        ! Get the error and gradient
!        e0=0.0_flyt
!        grad=0.0_flyt
!        do i=1,6
!            f0=dot_product(eta0,matmul(coeffM(:,:,i),eta0))-bvec(i)
!            e0=e0+f0**2
!            grad=grad+matmul( coeffM(:,:,i)+transpose(coeffM(:,:,i)),eta0)*2*f0
!        enddo
!
!        ! Line search in this direction
!        ctr=0
!        step=defstep
!        inloop: do initer=1,ninner
!            ! Guess new point
!            eta1=eta0-step*grad
!            ! Get error at new point
!            e1=0.0_flyt
!            do i=1,6
!                e1=e1+( dot_product(eta1,matmul(coeffM(:,:,i),eta1))-bvec(i) )**2
!            enddo
!
!            ! See how it went:
!            if ( e1 .lt. errtol ) then
!                ! We are done
!                eta0=eta1
!                e0=e1
!                exit outloop
!            elseif ( e1 .lt. e0 ) then
!                ! Good step, keep going
!                eta0=eta1
!                e0=e1
!                step=step*1.3_flyt
!            else
!                ! Bad step, decrease and try again
!                step=step*0.1_flyt
!                ctr=ctr+1
!            endif
!
!            ! If we have decreased too many times, start over with new gradient
!            if ( ctr .ge. maxctr ) exit inloop
!        enddo inloop
!    enddo outloop
!
!    ! Convert to a strain matrix
!    strain(1,1)=eta0(1)
!    strain(2,2)=eta0(2)
!    strain(3,3)=eta0(3)
!    strain(1,2)=eta0(4)
!    strain(2,1)=eta0(4)
!    strain(1,3)=eta0(5)
!    strain(3,1)=eta0(5)
!    strain(2,3)=eta0(6)
!    strain(3,2)=eta0(6)
!    ! Get the rotation
!    m0=matmul(strain,ref_lattice)
!    rotation=matmul(m0,lo_invert3x3matrix(lattice))
!    rotation=rotation/lo_determ(rotation)**(1.0_flyt/3.0_flyt)
!
!    contains
!
!    ! These are magic matrices that define the quadratic form. Thanks Mathematica.
!    subroutine cfvars(cm1,cm2,cm3,cm4,cm5,cm6,basis)
!        real(flyt), dimension(6,6), intent(out) :: cm1,cm2,cm3,cm4,cm5,cm6
!        real(flyt), dimension(3,3), intent(in) :: basis
!        real(flyt) :: ax,ay,az,bx,by,bz,cx,cy,cz
!        ax=basis(1,1)
!        ay=basis(2,1)
!        az=basis(3,1)
!        bx=basis(1,2)
!        by=basis(2,2)
!        bz=basis(3,2)
!        cx=basis(1,3)
!        cy=basis(2,3)
!        cz=basis(3,3)
!
!        cm1(1,1)=ax**2
!        cm1(1,2)=0
!        cm1(1,3)=0
!        cm1(1,4)=2*ax*ay
!        cm1(1,5)=2*ax*az
!        cm1(1,6)=0
!        cm1(2,1)=0
!        cm1(2,2)=ay**2
!        cm1(2,3)=0
!        cm1(2,4)=2*ax*ay
!        cm1(2,5)=0
!        cm1(2,6)=2*ay*az
!        cm1(3,1)=0
!        cm1(3,2)=0
!        cm1(3,3)=az**2
!        cm1(3,4)=0
!        cm1(3,5)=2*ax*az
!        cm1(3,6)=2*ay*az
!        cm1(4,1)=0
!        cm1(4,2)=0
!        cm1(4,3)=0
!        cm1(4,4)=ax**2 + ay**2
!        cm1(4,5)=2*ay*az
!        cm1(4,6)=2*ax*az
!        cm1(5,1)=0
!        cm1(5,2)=0
!        cm1(5,3)=0
!        cm1(5,4)=0
!        cm1(5,5)=ax**2 + az**2
!        cm1(5,6)=2*ax*ay
!        cm1(6,1)=0
!        cm1(6,2)=0
!        cm1(6,3)=0
!        cm1(6,4)=0
!        cm1(6,5)=0
!        cm1(6,6)=ay**2 + az**2
!
!        cm2(1,1)=bx**2
!        cm2(1,2)=0
!        cm2(1,3)=0
!        cm2(1,4)=2*bx*by
!        cm2(1,5)=2*bx*bz
!        cm2(1,6)=0
!        cm2(2,1)=0
!        cm2(2,2)=by**2
!        cm2(2,3)=0
!        cm2(2,4)=2*bx*by
!        cm2(2,5)=0
!        cm2(2,6)=2*by*bz
!        cm2(3,1)=0
!        cm2(3,2)=0
!        cm2(3,3)=bz**2
!        cm2(3,4)=0
!        cm2(3,5)=2*bx*bz
!        cm2(3,6)=2*by*bz
!        cm2(4,1)=0
!        cm2(4,2)=0
!        cm2(4,3)=0
!        cm2(4,4)=bx**2 + by**2
!        cm2(4,5)=2*by*bz
!        cm2(4,6)=2*bx*bz
!        cm2(5,1)=0
!        cm2(5,2)=0
!        cm2(5,3)=0
!        cm2(5,4)=0
!        cm2(5,5)=bx**2 + bz**2
!        cm2(5,6)=2*bx*by
!        cm2(6,1)=0
!        cm2(6,2)=0
!        cm2(6,3)=0
!        cm2(6,4)=0
!        cm2(6,5)=0
!        cm2(6,6)=by**2 + bz**2
!
!        cm3(1,1)=cx**2
!        cm3(1,2)=0
!        cm3(1,3)=0
!        cm3(1,4)=2*cx*cy
!        cm3(1,5)=2*cx*cz
!        cm3(1,6)=0
!        cm3(2,1)=0
!        cm3(2,2)=cy**2
!        cm3(2,3)=0
!        cm3(2,4)=2*cx*cy
!        cm3(2,5)=0
!        cm3(2,6)=2*cy*cz
!        cm3(3,1)=0
!        cm3(3,2)=0
!        cm3(3,3)=cz**2
!        cm3(3,4)=0
!        cm3(3,5)=2*cx*cz
!        cm3(3,6)=2*cy*cz
!        cm3(4,1)=0
!        cm3(4,2)=0
!        cm3(4,3)=0
!        cm3(4,4)=cx**2 + cy**2
!        cm3(4,5)=2*cy*cz
!        cm3(4,6)=2*cx*cz
!        cm3(5,1)=0
!        cm3(5,2)=0
!        cm3(5,3)=0
!        cm3(5,4)=0
!        cm3(5,5)=cx**2 + cz**2
!        cm3(5,6)=2*cx*cy
!        cm3(6,1)=0
!        cm3(6,2)=0
!        cm3(6,3)=0
!        cm3(6,4)=0
!        cm3(6,5)=0
!        cm3(6,6)=cy**2 + cz**2
!
!        cm4(1,1)=ax*bx
!        cm4(1,2)=0
!        cm4(1,3)=0
!        cm4(1,4)=ay*bx + ax*by
!        cm4(1,5)=az*bx + ax*bz
!        cm4(1,6)=0
!        cm4(2,1)=0
!        cm4(2,2)=ay*by
!        cm4(2,3)=0
!        cm4(2,4)=ay*bx + ax*by
!        cm4(2,5)=0
!        cm4(2,6)=az*by + ay*bz
!        cm4(3,1)=0
!        cm4(3,2)=0
!        cm4(3,3)=az*bz
!        cm4(3,4)=0
!        cm4(3,5)=az*bx + ax*bz
!        cm4(3,6)=az*by + ay*bz
!        cm4(4,1)=0
!        cm4(4,2)=0
!        cm4(4,3)=0
!        cm4(4,4)=ax*bx + ay*by
!        cm4(4,5)=az*by + ay*bz
!        cm4(4,6)=az*bx + ax*bz
!        cm4(5,1)=0
!        cm4(5,2)=0
!        cm4(5,3)=0
!        cm4(5,4)=0
!        cm4(5,5)=ax*bx + az*bz
!        cm4(5,6)=ay*bx + ax*by
!        cm4(6,1)=0
!        cm4(6,2)=0
!        cm4(6,3)=0
!        cm4(6,4)=0
!        cm4(6,5)=0
!        cm4(6,6)=ay*by + az*bz
!
!        cm5(1,1)=bx*cx
!        cm5(1,2)=0
!        cm5(1,3)=0
!        cm5(1,4)=by*cx + bx*cy
!        cm5(1,5)=bz*cx + bx*cz
!        cm5(1,6)=0
!        cm5(2,1)=0
!        cm5(2,2)=by*cy
!        cm5(2,3)=0
!        cm5(2,4)=by*cx + bx*cy
!        cm5(2,5)=0
!        cm5(2,6)=bz*cy + by*cz
!        cm5(3,1)=0
!        cm5(3,2)=0
!        cm5(3,3)=bz*cz
!        cm5(3,4)=0
!        cm5(3,5)=bz*cx + bx*cz
!        cm5(3,6)=bz*cy + by*cz
!        cm5(4,1)=0
!        cm5(4,2)=0
!        cm5(4,3)=0
!        cm5(4,4)=bx*cx + by*cy
!        cm5(4,5)=bz*cy + by*cz
!        cm5(4,6)=bz*cx + bx*cz
!        cm5(5,1)=0
!        cm5(5,2)=0
!        cm5(5,3)=0
!        cm5(5,4)=0
!        cm5(5,5)=bx*cx + bz*cz
!        cm5(5,6)=by*cx + bx*cy
!        cm5(6,1)=0
!        cm5(6,2)=0
!        cm5(6,3)=0
!        cm5(6,4)=0
!        cm5(6,5)=0
!        cm5(6,6)=by*cy + bz*cz
!
!        cm6(1,1)=ax*cx
!        cm6(1,2)=0
!        cm6(1,3)=0
!        cm6(1,4)=ay*cx + ax*cy
!        cm6(1,5)=az*cx + ax*cz
!        cm6(1,6)=0
!        cm6(2,1)=0
!        cm6(2,2)=ay*cy
!        cm6(2,3)=0
!        cm6(2,4)=ay*cx + ax*cy
!        cm6(2,5)=0
!        cm6(2,6)=az*cy + ay*cz
!        cm6(3,1)=0
!        cm6(3,2)=0
!        cm6(3,3)=az*cz
!        cm6(3,4)=0
!        cm6(3,5)=az*cx + ax*cz
!        cm6(3,6)=az*cy + ay*cz
!        cm6(4,1)=0
!        cm6(4,2)=0
!        cm6(4,3)=0
!        cm6(4,4)=ax*cx + ay*cy
!        cm6(4,5)=az*cy + ay*cz
!        cm6(4,6)=az*cx + ax*cz
!        cm6(5,1)=0
!        cm6(5,2)=0
!        cm6(5,3)=0
!        cm6(5,4)=0
!        cm6(5,5)=ax*cx + az*cz
!        cm6(5,6)=ay*cx + ax*cy
!        cm6(6,1)=0
!        cm6(6,2)=0
!        cm6(6,3)=0
!        cm6(6,4)=0
!        cm6(6,5)=0
!        cm6(6,6)=ay*cy + az*cz
!end subroutine
!
!end subroutine
