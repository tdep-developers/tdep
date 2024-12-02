#include "precompilerdefinitions"
module scatteringstrengths
use konstanter, only: r8, lo_twopi
use gottochblandat, only: walltime, lo_looptimer
use mpi_wrappers, only: lo_mpi_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
! unique
use phononevents, only: lo_threephononevents
implicit none

private
public :: calculate_scattering_amplitudes
contains

!> get all matrix elements
subroutine calculate_scattering_amplitudes(uc, qp, sc, dr, fct, mw)
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> integration weights
    type(lo_threephononevents), intent(inout) :: sc
    !> MPI helper
    type(lo_mpi_helper), intent(in) :: mw
    !
    real(r8), parameter :: timereport = 30.0_r8 ! report progress every 30 seconds
    !
    complex(r8), dimension(fct%na*3, 3) :: egv
    complex(r8), dimension(fct%na*3, 2) :: egviso
    complex(r8) :: c0
    real(r8), dimension(3) :: q2, q3, omega
    real(r8) :: omthres, timer, t0
    integer :: gi1, gi2, gi3
    integer :: i, j, b1, b2, b3

    timer = walltime()
    t0 = timer

    ! calculate the matrix elements
    omthres = dr%omega_min*0.5_r8
    do i = 1, sc%n_local_qpoint
        do b1 = 1, dr%n_mode
        do b2 = 1, dr%n_mode
        do b3 = 1, dr%n_mode

            ! First the plus-events
            do j = 1, sc%q(i)%plus(b1, b2, b3)%n
                ! grid-indices
                gi1 = sc%q(i)%gi1
                gi2 = sc%q(i)%plus(b1, b2, b3)%e(j)%gi2
                gi3 = sc%q(i)%plus(b1, b2, b3)%e(j)%gi3
                ! q-vectors with correct sign
                q2 = -qp%ap(gi2)%r*lo_twopi
                q3 = -qp%ap(gi3)%r*lo_twopi
                ! frequencies, eigenvectors, q-vectors
                omega(1) = dr%aq(gi1)%omega(b1)
                omega(2) = dr%aq(gi2)%omega(b2)
                omega(3) = dr%aq(gi3)%omega(b3)
                if (minval(omega) .lt. omthres) then
                    sc%q(i)%plus(b1, b2, b3)%e(j)%psisquare = 0.0_r8
                else
                    egv(:, 1) = dr%aq(gi1)%egv(:, b1)
                    egv(:, 2) = dr%aq(gi2)%egv(:, b2)
                    egv(:, 3) = dr%aq(gi3)%egv(:, b3)
                    ! and the scattering amplitude
                    c0 = fct%scatteringamplitude(omega, egv, q2, q3)
                    !sc%q(i)%plus(b1,b2,b3)%e(j)%psisquare=1E-10_r8 !abs(conjg(c0)*c0)
                    sc%q(i)%plus(b1, b2, b3)%e(j)%psisquare = abs(conjg(c0)*c0)
                    !write(*,*) abs(conjg(c0)*c0)
                    !write(*,*) abs(conjg(c0)*c0)
                end if
            end do

            ! Then the minus-events
            do j = 1, sc%q(i)%minus(b1, b2, b3)%n
                ! grid-indices
                gi1 = sc%q(i)%gi1
                gi2 = sc%q(i)%minus(b1, b2, b3)%e(j)%gi2
                gi3 = sc%q(i)%minus(b1, b2, b3)%e(j)%gi3
                ! q-vectors with correct sign
                q2 = -qp%ap(gi2)%r*lo_twopi
                q3 = -qp%ap(gi3)%r*lo_twopi
                ! frequencies, eigenvectors, q-vectors
                omega(1) = dr%aq(gi1)%omega(b1)
                omega(2) = dr%aq(gi2)%omega(b2)
                omega(3) = dr%aq(gi3)%omega(b3)
                if (minval(omega) .lt. omthres) then
                    sc%q(i)%minus(b1, b2, b3)%e(j)%psisquare = 0.0_r8
                else
                    egv(:, 1) = dr%aq(gi1)%egv(:, b1)
                    egv(:, 2) = dr%aq(gi2)%egv(:, b2)
                    egv(:, 3) = dr%aq(gi3)%egv(:, b3)
                    ! and the scattering amplitude
                    c0 = fct%scatteringamplitude(omega, egv, q2, q3)
                    !sc%q(i)%minus(b1,b2,b3)%e(j)%psisquare=1E-10_r8 ! abs(conjg(c0)*c0)
                    sc%q(i)%minus(b1, b2, b3)%e(j)%psisquare = abs(conjg(c0)*c0)
                end if
            end do
        end do
        end do
        end do

        ! And the isotope part
        if (sc%isotopescattering) then
            do b1 = 1, dr%n_mode
            do b2 = 1, dr%n_mode
                do j = 1, sc%iq(i)%band(b1, b2)%n
                    ! For q
                    gi1 = sc%iq(i)%gi1
                    egviso(:, 1) = dr%aq(gi1)%egv(:, b1)
                    ! For q'
                    gi2 = sc%iq(i)%band(b1, b2)%e(j)%gi2
                    egviso(:, 2) = dr%aq(gi2)%egv(:, b2)
                    ! Amplitude
                    sc%iq(i)%band(b1, b2)%e(j)%scatterstrength = isotope_scattering_strength(uc, egviso)
                end do
            end do
            end do
        end if

        if (mw%talk) then
        if (walltime() - t0 .gt. timereport) then
            call lo_looptimer('... calculating scatteringrates', timer, walltime(), i, sc%n_local_qpoint)
            t0 = walltime()
        end if
        end if
    end do

    if (mw%talk) then
        write (*, *) 'done matrixelements', walltime() - timer
        timer = walltime()
    end if
end subroutine

real(r8) function isotope_scattering_strength(uc, egv)
    type(lo_crystalstructure), intent(in) :: uc
    complex(r8), dimension(:, :), intent(in) :: egv
    !
    integer :: i, j
    real(r8) :: f0, f1
    complex(r8), dimension(3) :: cv0, cv1
    !
    f1 = 0.0_r8
    do i = 1, uc%na
        cv0 = egv((i - 1)*3 + 1:(i*3), 1)
        cv1 = egv((i - 1)*3 + 1:(i*3), 2)
        f0 = 0.0_r8
        do j = 1, 3
            f0 = f0 + abs(conjg(cv0(j))*cv1(j))
        end do
        f0 = f0**2
        f1 = f1 + f0*uc%isotope(i)%disorderparameter
    end do
    isotope_scattering_strength = f1
    !
end function

!!> the supposedly fast and clever matrix element thing
!pure function flat_matrixelement(om1,om2,om3,q2,q3,egv1,egv2,egv3,fc,ai,lv2,lv3) result(psi)
!        !> 1/sqrt(omega)
!        real(r8), intent(in) :: om1,om2,om3
!        !> normal eigenvectors
!        complex(r8), dimension(:), intent(in) :: egv1,egv2,egv3
!        !> normal q-vectors
!        real(r8), dimension(3), intent(in) :: q2,q3
!        !> flattened forceconstant
!        real(r8), dimension(:,:,:,:), intent(in) :: fc
!        !> indices to the atoms per forceconstant
!        integer, dimension(:,:), intent(in) :: ai
!        !> lattice vector per forceconstant
!        real(r8), dimension(:,:), intent(in) :: lv2,lv3
!        !> the actual matrix element
!        complex(r8) :: psi
!
!        ! precomputed factor to the correct unit
!        real(r8), parameter :: enhet=2.367755374064161E51_r8
!        ! some normal stuff
!        complex(r8), dimension(:,:,:,:,:,:), allocatable :: egvprod
!        complex(r8) :: expiqr,c0
!        real(r8) :: iqr,omegaprod
!        integer :: atom1,atom2,atom3,i,j,k,k1,k2,k3,t,na,nb,nf
!        !
!        nb=size(egv1,1)
!        na=nb/3
!        nf=size(ai,2)
!
!        ! rearrange the eigenvectors
!        allocate(egvprod(3,3,3,na,na,na))
!        do atom3=1,na
!        do atom2=1,na
!        do atom1=1,na
!            do k=1,3
!            do j=1,3
!            do i=1,3
!                k1=(atom1-1)*3+i
!                k2=(atom2-1)*3+j
!                k3=(atom3-1)*3+k
!                egvprod(i,j,k,atom1,atom2,atom3)=egv1(k1)*egv2(k2)*egv3(k3)
!            enddo
!            enddo
!            enddo
!        enddo
!        enddo
!        enddo
!
!        ! get the rest of the stuff
!        c0=0.0_r8
!        do t=1,nf
!            atom1=ai(1,t)
!            atom2=ai(2,t)
!            atom3=ai(3,t)
!            iqr=dot_product(q2,lv2(:,t))+dot_product(q3,lv3(:,t))
!            expiqr=dcmplx(cos(iqr),sin(iqr))
!            c0=c0+sum(fc(:,:,:,t)*egvprod(:,:,:,atom1,atom2,atom3))*expiqr
!        enddo
!        omegaprod=om1*om2*om3
!        psi=c0*enhet*omegaprod
!
!end function

!!> flatten things for speed, or something
!subroutine flatten_things(fct,uc,dr,flat_invom,flat_fc,flat_ai,flat_lv2,flat_lv3)
!    !> third order force constants
!    type(lo_forceconstant_thirdorder), intent(in) :: fct
!    !> unitcell
!    type(lo_crystalstructure), intent(in) :: uc
!    !> dispersions
!    type(lo_phonon_dispersions), intent(in) :: dr
!    !> inverse sqrt(omega)
!    real(r8), dimension(:,:), allocatable, intent(out) :: flat_invom
!    !> lattice vector 2
!    real(r8), dimension(:,:), allocatable, intent(out) :: flat_lv2
!    !> lattice vector 3
!    real(r8), dimension(:,:), allocatable, intent(out) :: flat_lv3
!    !> the force constant
!    real(r8), dimension(:,:,:,:), allocatable, intent(out) :: flat_fc
!    !> atom indices
!    integer, dimension(:,:), allocatable, intent(out) :: flat_ai
!    !
!    integer :: a1,a2,a3
!    integer :: i,j,l,n
!
!    ! flattened frequencies
!    lo_allocate(flat_invom(dr%n_mode,dr%nq_tot))
!    flat_invom=0.0_r8
!    do i=1,dr%nq_tot
!    do j=1,dr%n_mode
!        if ( dr%aq(i)%omega(j) .gt. dr%omega_min*0.5_r8 ) then
!            flat_invom(j,i)=1.0_r8/sqrt(dr%aq(i)%omega(j))
!        else
!            flat_invom(j,i)=0.0_r8
!        endif
!    enddo
!    enddo
!
!    ! and the forceconstant stuff
!    n=0
!    do a1=1,fct%na
!        n=n+fct%atom(a1)%n
!    enddo
!    allocate(flat_fc(3,3,3,n))
!    allocate(flat_lv2(3,n))
!    allocate(flat_lv3(3,n))
!    allocate(flat_ai(3,n))
!
!    l=0
!    do a1=1,fct%na
!        do i=1,fct%atom(a1)%n
!            l=l+1
!            a2=fct%atom(a1)%triplet(i)%i2
!            a3=fct%atom(a1)%triplet(i)%i3
!            flat_ai(:,l)=[a1,a2,a3]
!            flat_lv2(:,l)=fct%atom(a1)%triplet(i)%lv2
!            flat_lv3(:,l)=fct%atom(a1)%triplet(i)%lv3
!            flat_fc(:,:,:,l)=fct%atom(a1)%triplet(i)%m*uc%invsqrtmass(a1)*uc%invsqrtmass(a2)*uc%invsqrtmass(a2)
!        enddo
!    enddo
!end subroutine

end module
