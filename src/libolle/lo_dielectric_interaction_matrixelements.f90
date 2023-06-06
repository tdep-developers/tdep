submodule(lo_dielectric_interaction) lo_dielectric_interaction_matrixelements
!!
!! Routines to calculate all manner of matrix elements involving susceptibilities and dipole moments
!!
use konstanter, only: lo_freqtol,lo_twopi
implicit none
contains

!> lowest order polarizability matrix element, slow version
module function epsilon_matrixelement_firstorder(di,nuvector) result(M)
    !> interaction tensors
    class(lo_dielectric_tensor), intent(in) :: di
    !> phonon eigenvector, scaled with sqrt(1/2 m w), zero if omega=0
    complex(r8), dimension(:), intent(in) :: nuvector
    !> matrix element
    complex(r8), dimension(3,3) :: M

    real(r8), dimension(3) :: v0
    integer :: iatom,ix

    ! if no Raman activity it is zero
    if ( di%n_eps_singlet .eq. 0 ) then
        M=0.0_r8
        return
    endif

    M=0.0_r8
    do iatom=1,di%n_eps_singlet
        v0=real(nuvector( (iatom-1)*3+1:iatom*3 ),r8)
        do ix=1,3
            M=M+di%eps_singlet(iatom)%m(:,:,ix)*v0(ix)
        enddo
    enddo
end function

!> second order polarizability matrix element, slow version
module function epsilon_matrixelement_secondorder(di,q,nuvector1,nuvector2) result(M)
    !> interaction tensors
    class(lo_dielectric_tensor), intent(in) :: di
    !> q-vector (without 2*pi)
    real(r8), dimension(3) :: q
    !> phonon eigenvectors, scaled with sqrt(1/2 m w), set to zero if w=0
    complex(r8), dimension(:), intent(in) :: nuvector1,nuvector2
    !> matrix element
    complex(r8), dimension(3,3) :: M

    complex(r8), dimension(3,3) :: cv0
    real(r8) :: qdotr
    complex(r8) :: expiqr
    integer :: ipair,a1,a2,i,j,k,l,ii,jj

    M=0.0_r8
    do ipair=1,di%n_eps_pair
        qdotr=dot_product(di%eps_pair(ipair)%lv,q)*lo_twopi
        expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
        a1=di%eps_pair(ipair)%a1
        a2=di%eps_pair(ipair)%a2
        cv0=0.0_r8
        do j=1,3
        do i=1,3
            ii=(a1-1)*3+i
            jj=(a2-1)*3+j
            do l=1,3
            do k=1,3
                cv0(k,l)=cv0(k,l)+nuvector1(ii)*conjg(nuvector2(jj))*di%eps_pair(ipair)%m(k,l,i,j)
            enddo
            enddo
        enddo
        enddo
        M=M+cv0*expiqr
    enddo
end function

!> lowest order dipole matrix element, slow version
module function polarization_matrixelement_firstorder(di,nuvector) result(M)
    !> interaction tensors
    class(lo_dielectric_tensor), intent(in) :: di
    !> phonon eigenvector, scaled with sqrt(1/2 m w), zero if omega=0
    complex(r8), dimension(:), intent(in) :: nuvector
    !> matrix element
    complex(r8), dimension(3) :: M

    real(r8), dimension(3) :: v0
    integer :: iatom

    ! if no Born charges it is zero
    if ( di%n_Z_singlet .eq. 0 ) then
        M=0.0_r8
        return
    endif

    M=0.0_r8
    do iatom=1,di%n_Z_singlet
        v0=real(nuvector( (iatom-1)*3+1:iatom*3 ),r8)
        M=M+matmul(di%Z_singlet(:,:,iatom),v0)
    enddo
end function

!> second order dipole matrix element
module function polarization_matrixelement_secondorder(di,q,nuvector1,nuvector2) result(M)
    !> interaction tensors
    class(lo_dielectric_tensor), intent(in) :: di
    !> q-vector (without 2*pi)
    real(r8), dimension(3) :: q
    !> phonon eigenvectors, scaled with sqrt(1/2 m w), set to zero if w=0
    complex(r8), dimension(:), intent(in) :: nuvector1,nuvector2
    !> matrix element
    complex(r8), dimension(3) :: M

    complex(r8), dimension(3) :: cv0
    real(r8) :: qdotr
    complex(r8) :: expiqr
    integer :: ipair,a1,a2,i,j,k,ii,jj

    ! also need actual tensors
    if ( di%n_Z_pair .eq. 0 ) then
        M=0.0_r8
        return
    endif

    M=0.0_r8
    do ipair=1,di%n_Z_pair
        ! phase factor
        qdotr=dot_product(di%Z_pair(ipair)%lv,q)*lo_twopi
        expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
        a1=di%Z_pair(ipair)%a1
        a2=di%Z_pair(ipair)%a2
        ! Get the prefactor
        cv0=0.0_r8
        do i=1,3
        do j=1,3
            ii=(a1-1)*3+i
            jj=(a2-1)*3+j
            do k=1,3
                cv0(k)=cv0(k)+nuvector1(ii)*conjg(nuvector2(jj))*di%Z_pair(ipair)%m(k,i,j)
            enddo
        enddo
        enddo
        M=M+cv0*expiqr
    enddo
end function

!> third order dipole matrix element
module function polarization_matrixelement_thirdorder(di,q2,q3,nuvector1,nuvector2,nuvector3) result(M)
    !> interaction tensors
    type(lo_dielectric_tensor), intent(in) :: di
    !> q-vectors (without 2*pi)
    real(r8), dimension(3) :: q2,q3
    !> phonon eigenvectors
    complex(r8), dimension(:), intent(in) :: nuvector1,nuvector2,nuvector3
    !> matrix element
    complex(r8), dimension(3) :: M

    complex(r8), dimension(3) :: cv0
    real(r8) :: qdotr
    complex(r8) :: expiqr
    integer :: itriplet,a1,a2,a3,ix,iy,iz,imu,ialpha,ibeta,igamma

    ! also need actual tensors
    if ( di%n_Z_triplet .eq. 0 ) then
        M=0.0_r8
        return
    endif

    M=0.0_r8
    do itriplet=1,di%n_Z_pair
        qdotr=dot_product(di%Z_triplet(itriplet)%lv2,q2)*lo_twopi + dot_product(di%Z_triplet(itriplet)%lv3,q3)*lo_twopi
        expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
        a1=di%Z_triplet(itriplet)%a1
        a2=di%Z_triplet(itriplet)%a2
        a3=di%Z_triplet(itriplet)%a3
        cv0=0.0_r8
        do ix=1,3
        do iy=1,3
        do iz=1,3
            ialpha = (a1-1)*3+ix
            ibeta  = (a2-1)*3+iy
            igamma = (a3-1)*3+iz
            do imu=1,3
                cv0(imu)=cv0(imu)+&
                nuvector1(ialpha)*&
                nuvector2(ibeta)*&
                nuvector3(igamma)*&
                di%Z_triplet(itriplet)%m(imu,ix,iy,iz)
            enddo
        enddo
        enddo
        enddo
        M=M+cv0*expiqr
    enddo
end function

!> pre-transform second order Z
module subroutine pretransform_Z_secondorder(di,q,ptf)
    !> interaction tensors
    class(lo_dielectric_tensor), intent(in) :: di
    !> q-vector (without 2*pi)
    real(r8), dimension(3) :: q
    !> Fourier-transformed matrix element (3,n_mode**2)
    complex(r8), dimension(:,:), intent(out) :: ptf

    real(r8) :: qdotr
    complex(r8), dimension(3) :: cv3
    complex(r8) :: expiqr
    integer :: ipair,a1,a2,ix,iy,imu,ialpha,ibeta,ixj,n_mode

    n_mode=di%n_atom*3
    ptf=0.0_r8
    do ipair=1,di%n_Z_pair
        ! no minus?!
        qdotr=dot_product(q,di%Z_pair(ipair)%lv)*lo_twopi
        expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
        a1=di%z_pair(ipair)%a1
        a2=di%z_pair(ipair)%a2
        do ix=1,3
        do iy=1,3
            ! flatten electric field indices and fetch partial matrix element
            do imu=1,3
                cv3(imu)=di%Z_pair(ipair)%m(imu,ix,iy)
            enddo
            cv3=cv3*expiqr
            ! flatten eigenvector indices
            ialpha=(a1-1)*3+ix
            ibeta =(a2-1)*3+iy
            ixj = (ialpha-1)*n_mode + ibeta
            ! accumulate in the right place
            ptf(:,ixj)=ptf(:,ixj)+cv3
        enddo
        enddo
    enddo
end subroutine

!> pre-transform second order eps
module subroutine pretransform_epsilon_secondorder(di,q,ptf)
    !> interaction tensors
    class(lo_dielectric_tensor), intent(in) :: di
    !> q-vector (without 2*pi)
    real(r8), dimension(3) :: q
    !> Fourier-transformed matrix element (9,n_mode**2)
    complex(r8), dimension(:,:), intent(out) :: ptf

    real(r8) :: qdotr
    complex(r8), dimension(9) :: cv9
    complex(r8) :: expiqr
    integer :: ipair,a1,a2,ix,iy,imu,inu,ialpha,ibeta,ixi,ixj,n_mode

    n_mode=di%n_atom*3
    ptf=0.0_r8
    do ipair=1,di%n_eps_pair
        ! no minus?!
        qdotr=dot_product(q,di%eps_pair(ipair)%lv)*lo_twopi
        expiqr=cmplx(cos(qdotr),sin(qdotr),r8)
        a1=di%eps_pair(ipair)%a1
        a2=di%eps_pair(ipair)%a2
        do ix=1,3
        do iy=1,3
            ! flatten electric field indices and fetch partial matrix element
            do imu=1,3
            do inu=1,3
                ixi=(imu-1)*3+inu
                cv9(ixi)=di%eps_pair(ipair)%m(imu,inu,ix,iy)
            enddo
            enddo
            cv9=cv9*expiqr
            ! flatten eigenvector indices
            ialpha=(a1-1)*3+ix
            ibeta =(a2-1)*3+iy
            ixj = (ialpha-1)*n_mode + ibeta
            ! accumulate in the right place
            ptf(:,ixj)=ptf(:,ixj)+cv9
        enddo
        enddo
    enddo
end subroutine

end submodule
