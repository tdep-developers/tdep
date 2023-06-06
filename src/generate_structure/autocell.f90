#include "precompilerdefinitions"
module autocell
use konstanter, only: flyt, lo_tol, lo_pi
use gottochblandat, only: lo_determ
use geometryfunctions, only: lo_inscribed_sphere_in_box
use type_crystalstructure, only: lo_crystalstructure

implicit none
private
public :: return_supercellmatrix

contains

!> return a good supercell matrix
subroutine return_supercellmatrix(p, na, supercellmatrix)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: p
    !> desired number of atoms
    integer, intent(in) :: na
    !> fashionable supercell matrix
    integer, dimension(3, 3), intent(out) :: supercellmatrix

    real(flyt), dimension(3, 3) :: m0, m1, m2, basis, guessm
    real(flyt) :: det, f0, fillratio, r0, volumetolerance
    real(flyt) :: perfect_fill, sphvolf, volfactor
    integer, dimension(3, 3) :: guessmatrix
    integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, nrep
    logical :: foundsomething

    ! Some parameters. First the ideal fill-ratio of a unit sphere in a unit cube.
    ! if I achieve this, I can cancel the loop or something
    basis = p%latticevectors
    volfactor = (na*1.0_flyt)/p%na
    perfect_fill = 0.523598775598299_flyt
    sphvolf = 4.0_flyt*lo_pi/3.0_flyt

    ! Start by getting a decent guess on how to build a cubic matrix, without regard
    ! to how many atoms there are and so on.
    fillratio = 0.0_flyt ! reset the fillratio
    guessm = 0.0_flyt    ! reset the guess matrix
    ! this is probably really dangerous, but whatever
    reploop: do nrep = 1, 4
    do i1 = nrep, -nrep, -1
    do i2 = nrep, -nrep, -1
    do i3 = nrep, -nrep, -1
    do i4 = nrep, -nrep, -1
    do i5 = nrep, -nrep, -1
    do i6 = nrep, -nrep, -1
    do i7 = nrep, -nrep, -1
    do i8 = nrep, -nrep, -1
    do i9 = nrep, -nrep, -1
        m0(:, 1) = [i1, i2, i3]
        m0(:, 2) = [i4, i5, i6]
        m0(:, 3) = [i7, i8, i9]
        det = lo_determ(m0)
        if (det .lt. lo_tol) cycle
        m1 = matmul(basis, m0)               ! apply the supercell matrix
        det = abs(lo_determ(m1))            ! get the volume
        r0 = lo_inscribed_sphere_in_box(m1) ! radius of inscribed sphere
        f0 = sphvolf*r0**3/det              ! filling ratio
        if (f0 - fillratio .gt. lo_tol) then
            fillratio = f0
            guessm = m0
            if (fillratio/perfect_fill .gt. 0.96_flyt) exit reploop
        end if
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do reploop

    ! Now I have a decent guess. Adjust this matrix so that it has exactly
    ! the correct determinant
    det = lo_determ(guessm)
    m2 = guessm/(det**(1.0_flyt/3.0_flyt))
    m2 = m2*(volfactor**(1.0_flyt/3.0_flyt))
    ! for the guess, I will round this.
    guessmatrix = int(anint(m2))

    ! Now that I have an approximate supercell matrix, search a bit
    ! around the guess to find a really good one
    nrep = 2                              ! probably an ok guess
    fillratio = 0.0_flyt                  ! reset the fillratio
    m2 = 0.0_flyt                         ! reset the guess matrix
    volumetolerance = volfactor*0.2_flyt  ! initial tolerance for number of atoms
    foundsomething = .false.
    do i1 = nrep, -nrep, -1
    do i2 = nrep, -nrep, -1
    do i3 = nrep, -nrep, -1
    do i4 = nrep, -nrep, -1
    do i5 = nrep, -nrep, -1
    do i6 = nrep, -nrep, -1
    do i7 = nrep, -nrep, -1
    do i8 = nrep, -nrep, -1
    do i9 = nrep, -nrep, -1
        m0(:, 1) = [i1, i2, i3]
        m0(:, 2) = [i4, i5, i6]
        m0(:, 3) = [i7, i8, i9]
        m0 = m0 + guessmatrix
        det = lo_determ(m0)
        if (det .lt. volfactor - lo_tol) cycle
        if (abs(det - volfactor) .gt. volumetolerance) cycle
        m1 = matmul(basis, m0)               ! apply the supercell matrix
        det = abs(lo_determ(m1))            ! get the volume
        r0 = lo_inscribed_sphere_in_box(m1) ! radius of inscribed sphere
        f0 = sphvolf*r0**3/det              ! filling ratio
        if (f0 - fillratio .gt. lo_tol) then
            fillratio = f0
            m2 = m0
            foundsomething = .true.
        end if
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    ! Return the results, perhaps?
    if (foundsomething) then
        supercellmatrix = int(anint(m2))
    else
        write (*, *) 'Could not find any supercell. This is easily fixable:'
        write (*, *) 'try with a different number of desired atoms, per default only'
        write (*, *) 'supercells within 20% of this number of atoms are accepted.'
        stop
    end if

end subroutine

end module
