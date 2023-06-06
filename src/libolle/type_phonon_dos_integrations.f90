submodule (type_phonon_dos) type_phonon_dos_integrations
use konstanter, only: lo_kappa_SI_to_au,lo_imag
use gottochblandat, only: lo_harmonic_oscillator_cv,lo_outerproduct
implicit none
contains

!> Get the spectral phonon angular momentum.
module subroutine spectral_angular_momentum(pd,uc,qp,dr,temperature,mw,mem,spec_angmom,spec_angmom_band,spec_angmom_atom)
    !> phonon dos
    class(lo_phonon_dos), intent(inout) :: pd
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> q-point mesh
    type(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> spectral
    real(r8), dimension(:,:), allocatable, intent(out) :: spec_angmom
    !> spectral per band
    real(r8), dimension(:,:,:), allocatable, intent(out) :: spec_angmom_band
    !> spectral per atom
    real(r8), dimension(:,:,:), allocatable, intent(out) :: spec_angmom_atom

    real(r8), parameter :: angmomtol=1E-30_r8
    real(r8), dimension(:,:,:), allocatable :: angmom_flat
    real(r8), dimension(9) :: angmom_total

    call mem%tick()
    ! First I have to pre-calculate the angular momentum matrix thingy.
    init: block
        complex(r8), dimension(3,3) :: Mx,My,Mz
        complex(r8), dimension(3) :: cv0,cv1,cv2
        real(r8), dimension(3) :: v0,v1,w0,w1
        real(r8), dimension(3,3) :: m0
        real(r8), dimension(:,:,:), allocatable :: kf
        real(r8) :: f0
        integer :: iq,ib,ii,i,k,o

        ! Angular momentum representation:
        Mx=0.0_r8
        My=0.0_r8
        Mz=0.0_r8
        Mx(2,3)=1
        Mx(3,2)=-1
        My(1,3)=-1
        My(3,1)=1
        Mz(1,2)=1
        Mz(2,1)=-1
        Mx=-Mx*lo_imag
        My=-My*lo_imag
        Mz=-Mz*lo_imag
        call mem%allocate(angmom_flat,[qp%n_full_point,dr%n_mode,9],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(kf,[9,qp%n_irr_point,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        angmom_flat=0.0_r8
        kf=0.0_r8
        do iq=1,qp%n_full_point
            if ( mod(iq,mw%n) .ne. mw%r ) cycle
            ii=qp%ap(iq)%irreducible_index
            do ib=1,dr%n_mode
                ! Skip acoustic
                if ( dr%aq(iq)%omega(ib) .lt. lo_freqtol ) cycle

                ! Calculate the angular momentum contribution from this mode
                cv0=0.0_r8
                do k=1,uc%na
                    cv1=dr%aq(iq)%egv( (k-1)*3+1:k*3,ib )
                    cv2=matmul(Mx,cv1)
                    cv0(1)=cv0(1)+dot_product(cv1,cv2)
                    cv2=matmul(My,cv1)
                    cv0(2)=cv0(2)+dot_product(cv1,cv2)
                    cv2=matmul(Mz,cv1)
                    cv0(3)=cv0(3)+dot_product(cv1,cv2)
                enddo
                ! Now average over the small group. Seems sensible? Yes no maybe.
                v0=0.0_r8
                w0=0.0_r8
                v1=real(cv0)
                w1=dr%aq(iq)%vel(:,ib)
                do k=1,qp%ap(iq)%n_invariant_operation
                    o=qp%ap(iq)%invariant_operation(k)
                    v0=v0+matmul(uc%sym%op(o)%m,v1)
                    w0=w0+matmul(uc%sym%op(o)%m,w1)
                enddo
                v0=v0/real(qp%ap(iq)%n_invariant_operation,r8)
                w0=w0/real(qp%ap(iq)%n_invariant_operation,r8)

                f0=lo_harmonic_oscillator_cv( temperature,dr%aq(iq)%omega(ib) )
                f0=f0/dr%aq(iq)%omega(ib)
                f0=f0*qp%ap(iq)%integration_weight
                f0=f0*(2.0_r8*dr%iq(ii)%linewidth(ib))
                v0=real(cv0,r8)*f0
                m0=lo_outerproduct(v0,w0)

                kf(1,ii,ib)=kf(1,ii,ib)+m0(1,1)
                kf(2,ii,ib)=kf(2,ii,ib)+m0(1,2)
                kf(3,ii,ib)=kf(3,ii,ib)+m0(1,3)
                kf(4,ii,ib)=kf(4,ii,ib)+m0(2,1)
                kf(5,ii,ib)=kf(5,ii,ib)+m0(2,2)
                kf(6,ii,ib)=kf(6,ii,ib)+m0(2,3)
                kf(7,ii,ib)=kf(7,ii,ib)+m0(3,1)
                kf(8,ii,ib)=kf(8,ii,ib)+m0(3,2)
                kf(9,ii,ib)=kf(9,ii,ib)+m0(3,3)
            enddo
        enddo
        call mw%allreduce('sum',kf)
        ! Now rotate this out again
        do iq=1,qp%n_full_point
            ii=qp%ap(iq)%irreducible_index
            f0=1.0_r8/qp%ip(ii)%integration_weight
            do ib=1,dr%n_mode
                angmom_flat(iq,ib,1)=kf(1,ii,ib)*f0
                angmom_flat(iq,ib,2)=kf(2,ii,ib)*f0
                angmom_flat(iq,ib,3)=kf(3,ii,ib)*f0
                angmom_flat(iq,ib,4)=kf(4,ii,ib)*f0
                angmom_flat(iq,ib,5)=kf(5,ii,ib)*f0
                angmom_flat(iq,ib,6)=kf(6,ii,ib)*f0
                angmom_flat(iq,ib,7)=kf(7,ii,ib)*f0
                angmom_flat(iq,ib,8)=kf(8,ii,ib)*f0
                angmom_flat(iq,ib,9)=kf(9,ii,ib)*f0
            enddo
        enddo
        call mem%deallocate(kf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ! And grab the total angular momentum
        do i=1,9
            angmom_total(i)=sum(angmom_flat(:,:,i))/qp%n_full_point
        enddo
        ! Remove tiny numbers
        f0=norm2(angmom_total)
        angmom_total=lo_chop(angmom_total,f0*1E-9_r8)

        ! Space for the spectral angular momentum
        allocate(spec_angmom(pd%n_dos_point,9))
        allocate(spec_angmom_band(pd%n_dos_point,dr%n_mode,9))
        allocate(spec_angmom_atom(pd%n_dos_point,uc%na,9))
        spec_angmom=0.0_r8
        spec_angmom_band=0.0_r8
        spec_angmom_atom=0.0_r8
    end block init

    select case(pd%integrationtype)
    case(1:2) ! Gaussian-type integrations
    gaussint: block
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:,:,:), allocatable :: qproj
        real(r8), dimension(:,:), allocatable :: qvals,qangmom
        real(r8) :: sigma,f0,f1,invf,foursigma,weight
        integer :: iq,lq,mode,atom,xval,ii,jj,direction,nq

        nq=0
        do iq=1,qp%n_irr_point
            if ( mod(iq,mw%n) .ne. mw%r ) cycle
            nq=nq+1
        enddo

        if ( nq .gt. 0 ) then
            call mem%allocate(qvals,[nq,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(qproj,[uc%na,nq,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(qangmom,[nq,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            qvals=0.0_r8
            qproj=0.0_r8
            qangmom=0.0_r8
            lq=0
            do iq=1,qp%n_irr_point
                if ( mod(iq,mw%n) .ne. mw%r ) cycle
                lq=lq+1
                do mode=1,dr%n_mode
                    qvals(lq,mode)=dr%iq(iq)%omega(mode)
                    do atom=1,uc%na
                        cv0=dr%iq(iq)%egv( (atom-1)*3+1:atom*3 , mode )
                        qproj(atom,lq,mode)=abs(dot_product(conjg(cv0),cv0))
                    enddo
                enddo
            enddo
        endif

        ! Do an actual integration. Should be the same integration type as before.
        directionloop: do direction=1,9
            ! skip irrelevant directions
            if ( abs(angmom_total(direction)) .lt. angmomtol ) cycle
            ! Fetch angmom values
            lq=0
            do iq=1,qp%n_irr_point
                if ( mod(iq,mw%n) .ne. mw%r ) cycle
                lq=lq+1
                ii=qp%ip(iq)%full_index
                do mode=1,dr%n_mode
                    qangmom(lq,mode)=angmom_flat(ii,mode,direction)
                enddo
            enddo

            ! Do the actual integration
            invf=(pd%n_dos_point/maxval(pd%omega))
            lq=0
            do iq=1,qp%n_irr_point
                if ( mod(iq,mw%n) .ne. mw%r ) cycle
                lq=lq+1
                weight=qp%ip(iq)%integration_weight
                do mode=1,dr%n_mode
                    f1=qvals(lq,mode) ! frequency
                    if ( f1 .lt. dr%omega_min*0.5_r8 ) cycle

                    ! Get the smearing parameter
                    select case(pd%integrationtype)
                        case(1) ! fix gaussian
                            sigma=dr%default_smearing(mode)*pd%smearing_prefactor
                        case(2) ! adaptive gaussian
                            sigma=qp%smearingparameter(dr%iq(iq)%vel(:,mode),dr%default_smearing(mode),pd%smearing_prefactor)
                    end select
                    foursigma=sigma*4

                    ! range of integration
                    ii=max( floor( (f1-foursigma)*invf ), 1)
                    jj=min( floor( (f1+foursigma)*invf ), pd%n_dos_point)
                    do xval=ii,jj
                        f0=lo_gauss(f1,pd%omega(xval),sigma)*weight
                        spec_angmom_band(xval,mode,direction)=spec_angmom_band(xval,mode,direction)+qangmom(lq,mode)*f0
                        do atom=1,uc%na
                            spec_angmom_atom(xval,atom,direction)=spec_angmom_atom(xval,atom,direction)+qangmom(lq,mode)*qproj(atom,lq,mode)*f0
                        enddo
                    enddo
                enddo
            enddo
        enddo directionloop
        ! Add together
        call mw%allreduce('sum',spec_angmom_atom)
        call mw%allreduce('sum',spec_angmom_band)
        ! And cleanup
        if ( nq .gt. 0 ) then
            call mem%deallocate(qvals,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(qproj,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(qangmom,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
        endif
    end block gaussint
    case(3) ! Tetrahedron integration
    tetint: block
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:,:,:,:), allocatable :: tetproj
        real(r8), dimension(:,:,:), allocatable :: tetvals,tetangmom
        real(r8), dimension(4) :: w,tete
        real(r8) :: sigma,mine,maxe,invf
        integer :: i,l,mode,atom,xval,ii,jj,direction
        integer :: itet,ltet,ntet

        ! Count tetrahedrons per rank
        ntet=0
        do itet=1,qp%n_irr_tet
            if ( mod(itet,mw%n) .eq. mw%r ) ntet=ntet+1
        enddo

        ! Pre-fetch data
        if ( ntet .gt. 0 ) then
            call mem%allocate(tetvals,[4,ntet,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(tetangmom,[4,ntet,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(tetproj,[4,uc%na,ntet,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            tetvals=0.0_r8
            tetangmom=0.0_r8
            tetproj=0.0_r8
            ltet=0
            do itet=1,qp%n_irr_tet
                if ( mod(itet,mw%n) .ne. mw%r ) cycle
                ltet=ltet+1
                do i=1,4
                    l=qp%it(itet)%full_index(i)
                    do mode=1,dr%n_mode
                        tetvals(i,ltet,mode)=dr%aq(l)%omega(mode)
                        do atom=1,uc%na
                            cv0=dr%aq(l)%egv( (atom-1)*3+1:atom*3 , mode )
                            tetproj(i,atom,ltet,mode)=abs(dot_product(conjg(cv0),cv0))
                        enddo
                    enddo
                enddo
            enddo
        endif

        spec_angmom_atom=0.0_r8
        spec_angmom_band=0.0_r8
        directionloop: do direction=1,9
            ! angmom values for this direction
            tetangmom=0.0_r8
            ltet=0
            do itet=1,qp%n_irr_tet
                if ( mod(itet,mw%n) .ne. mw%r ) cycle
                ltet=ltet+1
                do i=1,4
                    l=qp%it(itet)%full_index(i)
                    do mode=1,dr%n_mode
                        tetangmom(i,ltet,mode)=angmom_flat(l,mode,direction)
                    enddo
                enddo
            enddo

            invf=(pd%n_dos_point/maxval(pd%omega))
            do mode=1,dr%n_mode
                sigma=dr%default_smearing(mode)
                ltet=0
                do itet=1,qp%n_irr_tet
                    if ( mod(itet,mw%n) .ne. mw%r ) cycle
                    ltet=ltet+1
                    tete=tetvals(:,ltet,mode)
                    mine=minval(tete)
                    maxe=maxval(tete)
                    ii=max(floor( mine*invf ),1)
                    jj=min(ceiling( maxe*invf ),pd%n_dos_point)
                    do xval=ii,jj
                        w=lo_LV_tetrahedron_weights(tete,pd%omega(xval),lo_freqtol,sigma)*qp%it(itet)%integration_weight
                        spec_angmom_band(xval,mode,direction)=spec_angmom_band(xval,mode,direction)+sum(w*tetangmom(:,ltet,mode))
                        do atom=1,uc%na
                            spec_angmom_atom(xval,atom,direction)=spec_angmom_atom(xval,atom,direction)+sum(w*tetangmom(:,ltet,mode)*tetproj(:,atom,ltet,mode))
                        enddo
                    enddo
                enddo
            enddo
        enddo directionloop
        call mw%allreduce('sum',spec_angmom_atom)
        call mw%allreduce('sum',spec_angmom_band)
        if ( ntet .gt. 0 ) then
            call mem%deallocate(tetvals,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(tetangmom,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(tetproj,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
        endif
    end block tetint
    case default
        call lo_stop_gracefully(['Unknown integration type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    end select
    ! Now we don't need the flat angular momentum anymore
    call mem%deallocate(angmom_flat,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! Now make the spectrum look a little nicer.
    makenice: block
        integer, parameter :: nkernel=6
        real(r8), dimension(-nkernel:nkernel) :: kernel
        real(r8), dimension(:), allocatable :: ysm
        real(r8) :: f0,f1,f2
        integer :: direction,i,j,k,atom,mode

        ! Prepare the smearing kernel
        do i=-nkernel,nkernel
            kernel(i)=lo_gauss(0.0_r8,i*1.0_r8,1.5_r8)
        enddo
        kernel=kernel/sum(kernel)
        allocate(ysm(-nkernel:pd%n_dos_point+nkernel))
        ysm=0.0_r8

        dirloop: do direction=1,9
            if ( abs(angmom_total(direction)) .lt. angmomtol ) cycle

            ! Smear things slightly
            do atom=1,uc%na
                ysm=0.0_r8
                do i=1,pd%n_dos_point
                do j=-nkernel,nkernel
                    ysm(i+j)=ysm(i+j)+kernel(j)*spec_angmom_atom(i,atom,direction)
                enddo
                enddo
                spec_angmom_atom(:,atom,direction)=ysm(1:pd%n_dos_point)
            enddo
            do mode=1,dr%n_mode
                ysm=0.0_r8
                do i=1,pd%n_dos_point
                do j=-nkernel,nkernel
                    ysm(i+j)=ysm(i+j)+kernel(j)*spec_angmom_band(i,mode,direction)
                enddo
                enddo
                spec_angmom_band(:,mode,direction)=ysm(1:pd%n_dos_point)
            enddo

            ! Fix site degeneracies
            do i=1,uc%na
                ysm=0.0_r8
                do j=1,uc%sym%degeneracy(i)
                    k=uc%sym%degenerate_atom(j,i)
                    ysm(1:pd%n_dos_point)=ysm(1:pd%n_dos_point)+spec_angmom_atom(:,i,direction)/(1.0_r8*uc%sym%degeneracy(i))
                enddo
                do j=1,uc%sym%degeneracy(i)
                    k=uc%sym%degenerate_atom(j,i)
                    spec_angmom_atom(:,k,direction)=ysm(1:pd%n_dos_point)
                enddo
            enddo

            ! Add up to total, and normalize properly
            spec_angmom(:,direction)=0.0_r8
            do mode=1,dr%n_mode
                spec_angmom(:,direction)=spec_angmom(:,direction)+spec_angmom_band(:,mode,direction)
            enddo
            f0=lo_trapezoid_integration(pd%omega,spec_angmom(:,direction))
            if ( abs(f0) .gt. angmomtol ) then
                spec_angmom(:,direction)=spec_angmom(:,direction)*angmom_total(direction)/f0
            else
                spec_angmom(:,direction)=0.0_r8
            endif

            ! ! Make sure the projected add up to the total.
            do i=1,pd%n_dos_point
                f0=sum(spec_angmom_atom(i,:,direction))
                f1=sum(spec_angmom_band(i,:,direction))
                f2=spec_angmom(i,direction)
                ! make sure there is finite angmom
                if ( abs(f2) .gt. angmomtol ) then
                    if ( abs(f0) .gt. angmomtol ) then
                        spec_angmom_atom(i,:,direction)=spec_angmom_atom(i,:,direction)*f2/f0
                    else
                        spec_angmom_atom(i,:,direction)=0.0_r8
                    endif
                    if ( abs(f1) .gt. angmomtol ) then
                        spec_angmom_band(i,:,direction)=spec_angmom_band(i,:,direction)*f2/f1
                    else
                        spec_angmom_band(i,:,direction)=0.0_r8
                    endif
                else
                    spec_angmom(i,direction)=0.0_r8
                    spec_angmom_band(i,:,direction)=0.0_r8
                    spec_angmom_atom(i,:,direction)=0.0_r8
                endif
            enddo
        enddo dirloop
        deallocate(ysm)
    end block makenice
    call mem%tock(__FILE__,__LINE__,mw%comm)
end subroutine

!> Get the spectral thermal conductivity, very close to a DOS integration
module subroutine spectral_kappa(pd,uc,qp,dr,mw,mem,spec_kappa,spec_kappa_band,spec_kappa_atom)
    !> phonon dos
    class(lo_phonon_dos), intent(inout) :: pd
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> q-point mesh
    type(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> spectral kappa
    real(r8), dimension(:,:), allocatable, intent(out) :: spec_kappa
    !> spectral kappa per band
    real(r8), dimension(:,:,:), allocatable, intent(out) :: spec_kappa_band
    !> spectral kappa per atom
    real(r8), dimension(:,:,:), allocatable, intent(out) :: spec_kappa_atom

    real(r8), parameter :: kappatol=1E-6_r8*lo_kappa_SI_to_au
    real(r8), dimension(:,:,:), allocatable :: kappa_flat
    real(r8), dimension(9) :: kappa_total

    call mem%tick()
    ! Pre-fetch some data
    init: block
        real(r8), dimension(:,:,:), allocatable :: kf
        real(r8) :: f0
        integer :: iq,ib,ii,i

        call mem%allocate(kappa_flat,[qp%n_full_point,dr%n_mode,9],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(kf,[9,qp%n_irr_point,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        kappa_flat=0.0_r8
        kf=0.0_r8
        do iq=1,qp%n_full_point
            if ( mod(iq,mw%n) .ne. mw%r ) cycle
            ii=qp%ap(iq)%irreducible_index
            do ib=1,dr%n_mode
                kf(1,ii,ib)=kf(1,ii,ib)+dr%aq(iq)%kappa(1,1,ib)
                kf(2,ii,ib)=kf(2,ii,ib)+dr%aq(iq)%kappa(1,2,ib)
                kf(3,ii,ib)=kf(3,ii,ib)+dr%aq(iq)%kappa(1,3,ib)
                kf(4,ii,ib)=kf(4,ii,ib)+dr%aq(iq)%kappa(2,1,ib)
                kf(5,ii,ib)=kf(5,ii,ib)+dr%aq(iq)%kappa(2,2,ib)
                kf(6,ii,ib)=kf(6,ii,ib)+dr%aq(iq)%kappa(2,3,ib)
                kf(7,ii,ib)=kf(7,ii,ib)+dr%aq(iq)%kappa(3,1,ib)
                kf(8,ii,ib)=kf(8,ii,ib)+dr%aq(iq)%kappa(3,2,ib)
                kf(9,ii,ib)=kf(9,ii,ib)+dr%aq(iq)%kappa(3,3,ib)
            enddo
        enddo
        call mw%allreduce('sum',kf)
        ! Now rotate this out again
        do iq=1,qp%n_full_point
            ii=qp%ap(iq)%irreducible_index
            f0=1.0_r8/qp%ip(ii)%integration_weight
            do ib=1,dr%n_mode
                kappa_flat(iq,ib,1)=kf(1,ii,ib)*f0
                kappa_flat(iq,ib,2)=kf(2,ii,ib)*f0
                kappa_flat(iq,ib,3)=kf(3,ii,ib)*f0
                kappa_flat(iq,ib,4)=kf(4,ii,ib)*f0
                kappa_flat(iq,ib,5)=kf(5,ii,ib)*f0
                kappa_flat(iq,ib,6)=kf(6,ii,ib)*f0
                kappa_flat(iq,ib,7)=kf(7,ii,ib)*f0
                kappa_flat(iq,ib,8)=kf(8,ii,ib)*f0
                kappa_flat(iq,ib,9)=kf(9,ii,ib)*f0
            enddo
        enddo
        call mem%deallocate(kf,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        kappa_flat=lo_chop(kappa_flat/qp%n_full_point,kappatol)
        ! And grab the total Kappa
        do i=1,9
            kappa_total(i)=sum(kappa_flat(:,:,i))/qp%n_full_point
        enddo

        ! Space for the spectral kappa
        allocate(spec_kappa(pd%n_dos_point,9))
        allocate(spec_kappa_band(pd%n_dos_point,dr%n_mode,9))
        allocate(spec_kappa_atom(pd%n_dos_point,uc%na,9))
        spec_kappa=0.0_r8
        spec_kappa_band=0.0_r8
        spec_kappa_atom=0.0_r8
    end block init

    select case(pd%integrationtype)
    case(1:2) ! Gaussian-type integrations
    gaussint: block
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:,:,:), allocatable :: qproj
        real(r8), dimension(:,:), allocatable :: qvals,qkappa
        real(r8) :: sigma,f0,f1,invf,foursigma,weight
        integer :: iq,lq,mode,atom,xval,ii,jj,direction,nq

        nq=0
        do iq=1,qp%n_irr_point
            if ( mod(iq,mw%n) .ne. mw%r ) cycle
            nq=nq+1
        enddo

        if ( nq .gt. 0 ) then
            call mem%allocate(qvals,[nq,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(qproj,[uc%na,nq,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(qkappa,[nq,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            qvals=0.0_r8
            qproj=0.0_r8
            qkappa=0.0_r8
            lq=0
            do iq=1,qp%n_irr_point
                if ( mod(iq,mw%n) .ne. mw%r ) cycle
                lq=lq+1
                do mode=1,dr%n_mode
                    qvals(lq,mode)=dr%iq(iq)%omega(mode)
                    do atom=1,uc%na
                        cv0=dr%iq(iq)%egv( (atom-1)*3+1:atom*3 , mode )
                        qproj(atom,lq,mode)=abs(dot_product(conjg(cv0),cv0))
                    enddo
                enddo
            enddo
        endif

        ! Do an actual integration. Should be the same integration type as before.
        directionloop: do direction=1,9
            ! skip irrelevant directions
            if ( abs(kappa_total(direction)) .lt. kappatol ) cycle
            ! Fetch kappa values
            lq=0
            do iq=1,qp%n_irr_point
                if ( mod(iq,mw%n) .ne. mw%r ) cycle
                lq=lq+1
                ii=qp%ip(iq)%full_index
                do mode=1,dr%n_mode
                    qkappa(lq,mode)=kappa_flat(ii,mode,direction)
                enddo
            enddo

            ! Do the actual integration
            invf=(pd%n_dos_point/maxval(pd%omega))
            lq=0
            do iq=1,qp%n_irr_point
                if ( mod(iq,mw%n) .ne. mw%r ) cycle
                lq=lq+1
                weight=qp%ip(iq)%integration_weight
                do mode=1,dr%n_mode
                    f1=qvals(lq,mode) ! frequency
                    if ( f1 .lt. dr%omega_min*0.5_r8 ) cycle

                    ! Get the smearing parameter
                    select case(pd%integrationtype)
                        case(1) ! fix gaussian
                            sigma=dr%default_smearing(mode)*pd%smearing_prefactor
                        case(2) ! adaptive gaussian
                            sigma=qp%smearingparameter(dr%iq(iq)%vel(:,mode),dr%default_smearing(mode),pd%smearing_prefactor)
                    end select
                    foursigma=sigma*4

                    ! range of integration
                    ii=max( floor( (f1-foursigma)*invf ), 1)
                    jj=min( floor( (f1+foursigma)*invf ), pd%n_dos_point)
                    do xval=ii,jj
                        f0=lo_gauss(f1,pd%omega(xval),sigma)*weight
                        spec_kappa_band(xval,mode,direction)=spec_kappa_band(xval,mode,direction)+qkappa(lq,mode)*f0
                        do atom=1,uc%na
                            spec_kappa_atom(xval,atom,direction)=spec_kappa_atom(xval,atom,direction)+qkappa(lq,mode)*qproj(atom,lq,mode)*f0
                        enddo
                    enddo
                enddo
            enddo
        enddo directionloop
        ! Add together
        call mw%allreduce('sum',spec_kappa_atom)
        call mw%allreduce('sum',spec_kappa_band)
        ! And cleanup
        if ( nq .gt. 0 ) then
            call mem%deallocate(qvals,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(qproj,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(qkappa,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
        endif
    end block gaussint
    case(3) ! Tetrahedron integration
    tetint: block
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:,:,:,:), allocatable :: tetproj
        real(r8), dimension(:,:,:), allocatable :: tetvals,tetkappa
        real(r8), dimension(4) :: w,tete
        real(r8) :: sigma,mine,maxe,invf
        integer :: i,l,mode,atom,xval,ii,jj,direction
        integer :: itet,ltet,ntet

        ! Count tetrahedrons per rank
        ntet=0
        do itet=1,qp%n_irr_tet
            if ( mod(itet,mw%n) .eq. mw%r ) ntet=ntet+1
        enddo

        ! Pre-fetch data
        if ( ntet .gt. 0 ) then
            call mem%allocate(tetvals,[4,ntet,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(tetkappa,[4,ntet,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%allocate(tetproj,[4,uc%na,ntet,dr%n_mode],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            tetvals=0.0_r8
            tetkappa=0.0_r8
            tetproj=0.0_r8
            ltet=0
            do itet=1,qp%n_irr_tet
                if ( mod(itet,mw%n) .ne. mw%r ) cycle
                ltet=ltet+1
                do i=1,4
                    l=qp%it(itet)%full_index(i)
                    do mode=1,dr%n_mode
                        tetvals(i,ltet,mode)=dr%aq(l)%omega(mode)
                        do atom=1,uc%na
                            cv0=dr%aq(l)%egv( (atom-1)*3+1:atom*3 , mode )
                            tetproj(i,atom,ltet,mode)=abs(dot_product(conjg(cv0),cv0))
                        enddo
                    enddo
                enddo
            enddo
        endif

        spec_kappa_atom=0.0_r8
        spec_kappa_band=0.0_r8
        directionloop: do direction=1,9
            ! kappa values for this direction
            tetkappa=0.0_r8
            ltet=0
            do itet=1,qp%n_irr_tet
                if ( mod(itet,mw%n) .ne. mw%r ) cycle
                ltet=ltet+1
                do i=1,4
                    l=qp%it(itet)%full_index(i)
                    do mode=1,dr%n_mode
                        tetkappa(i,ltet,mode)=kappa_flat(l,mode,direction)
                    enddo
                enddo
            enddo

            invf=(pd%n_dos_point/maxval(pd%omega))
            do mode=1,dr%n_mode
                ltet=0
                sigma=dr%default_smearing(mode)
                do itet=1,qp%n_irr_tet
                    if ( mod(itet,mw%n) .ne. mw%r ) cycle
                    ltet=ltet+1
                    tete=tetvals(:,ltet,mode)
                    mine=minval(tete)
                    maxe=maxval(tete)
                    ii=max(floor( mine*invf ),1)
                    jj=min(ceiling( maxe*invf ),pd%n_dos_point)
                    do xval=ii,jj
                        w=lo_LV_tetrahedron_weights(tete,pd%omega(xval),lo_freqtol,sigma)*qp%it(itet)%integration_weight
                        spec_kappa_band(xval,mode,direction)=spec_kappa_band(xval,mode,direction)+sum(w*tetkappa(:,ltet,mode))
                        do atom=1,uc%na
                            spec_kappa_atom(xval,atom,direction)=spec_kappa_atom(xval,atom,direction)+sum(w*tetkappa(:,ltet,mode)*tetproj(:,atom,ltet,mode))
                        enddo
                    enddo
                enddo
            enddo
        enddo directionloop
        call mw%allreduce('sum',spec_kappa_atom)
        call mw%allreduce('sum',spec_kappa_band)
        if ( ntet .gt. 0 ) then
            call mem%deallocate(tetvals,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(tetkappa,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
            call mem%deallocate(tetproj,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
        endif
    end block tetint
    case default
        call lo_stop_gracefully(['Unknown integration type'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    end select

    ! Now make the spectrum look a little nicer.
    makenice: block
        integer, parameter :: nkernel=6
        real(r8), dimension(-nkernel:nkernel) :: kernel
        real(r8), dimension(:), allocatable :: ysm
        real(r8) :: f0,f1,f2
        integer :: direction,i,j,k,atom,mode

        ! Prepare the smearing kernel
        do i=-nkernel,nkernel
            kernel(i)=lo_gauss(0.0_r8,i*1.0_r8,1.5_r8)
        enddo
        kernel=kernel/sum(kernel)
        allocate(ysm(-nkernel:pd%n_dos_point+nkernel))
        ysm=0.0_r8

        ! Must be positive, I think.
        spec_kappa_atom=max(spec_kappa_atom,0.0_r8)
        spec_kappa_band=max(spec_kappa_band,0.0_r8)

        dirloop: do direction=1,9
            if ( kappa_total(direction) .lt. kappatol ) cycle

            ! Smear things slightly
            do atom=1,uc%na
                ysm=0.0_r8
                do i=1,pd%n_dos_point
                do j=-nkernel,nkernel
                    ysm(i+j)=ysm(i+j)+kernel(j)*spec_kappa_atom(i,atom,direction)
                enddo
                enddo
                spec_kappa_atom(:,atom,direction)=ysm(1:pd%n_dos_point)
            enddo
            do mode=1,dr%n_mode
                ysm=0.0_r8
                do i=1,pd%n_dos_point
                do j=-nkernel,nkernel
                    ysm(i+j)=ysm(i+j)+kernel(j)*spec_kappa_band(i,mode,direction)
                enddo
                enddo
                spec_kappa_band(:,mode,direction)=ysm(1:pd%n_dos_point)
            enddo

            ! Fix site degeneracies
            do i=1,uc%na
                ysm=0.0_r8
                do j=1,uc%sym%degeneracy(i)
                    k=uc%sym%degenerate_atom(j,i)
                    ysm(1:pd%n_dos_point)=ysm(1:pd%n_dos_point)+spec_kappa_atom(:,i,direction)/(1.0_r8*uc%sym%degeneracy(i))
                enddo
                do j=1,uc%sym%degeneracy(i)
                    k=uc%sym%degenerate_atom(j,i)
                    spec_kappa_atom(:,k,direction)=ysm(1:pd%n_dos_point)
                enddo
            enddo

            ! Has to be zero at zero
            do atom=1,uc%na
                f0=spec_kappa_atom(1,atom,direction)
                do i=1,pd%n_dos_point
                    f1=1.0_r8-(i-1.0_r8)/real(pd%n_dos_point-1,r8)
                    spec_kappa_atom(i,atom,direction)=max(spec_kappa_atom(i,atom,direction)-f1*f0,0.0_r8)
                enddo
                spec_kappa_atom(1,atom,direction)=0.0_r8
            enddo
            do mode=1,dr%n_mode
                f0=spec_kappa_band(1,mode,direction)
                do i=1,pd%n_dos_point
                    f1=1.0_r8-(i-1.0_r8)/real(pd%n_dos_point-1,r8)
                    spec_kappa_band(i,mode,direction)=max(spec_kappa_band(i,mode,direction)-f1*f0,0.0_r8)
                enddo
                spec_kappa_band(1,mode,direction)=0.0_r8
            enddo

            ! Add up to total, and normalize properly
            spec_kappa(:,direction)=0.0_r8
            do mode=1,dr%n_mode
                spec_kappa(:,direction)=spec_kappa(:,direction)+spec_kappa_band(:,mode,direction)
            enddo
            f0=lo_trapezoid_integration(pd%omega,spec_kappa(:,direction))
            if ( f0 .gt. kappatol ) then
                spec_kappa(:,direction)=spec_kappa(:,direction)*kappa_total(direction)/f0
            else
                spec_kappa(:,direction)=0.0_r8
            endif

            ! Make sure the projected add up to the total.
            do i=1,pd%n_dos_point
                f0=sum(spec_kappa_atom(i,:,direction))
                f1=sum(spec_kappa_band(i,:,direction))
                f2=spec_kappa(i,direction)
                ! make sure there is finite kappa
                if ( f2 .gt. kappatol ) then
                    if ( f0 .gt. kappatol ) then
                        spec_kappa_atom(i,:,direction)=spec_kappa_atom(i,:,direction)*f2/f0
                    else
                        spec_kappa_atom(i,:,direction)=0.0_r8
                    endif
                    if ( f1 .gt. kappatol ) then
                        spec_kappa_band(i,:,direction)=spec_kappa_band(i,:,direction)*f2/f1
                    else
                        spec_kappa_band(i,:,direction)=0.0_r8
                    endif
                else
                    spec_kappa(i,direction)=0.0_r8
                    spec_kappa_band(i,:,direction)=0.0_r8
                    spec_kappa_atom(i,:,direction)=0.0_r8
                endif
            enddo
        enddo dirloop
        deallocate(ysm)
    end block makenice
end subroutine

!> Integrate the phonon dos with the tetrahedron method
module subroutine lo_get_phonon_dos_tetrahedron(pd,qp,dr,mw,mem,verbosity)
    !> phonon dos
    type(lo_phonon_dos), intent(inout) :: pd
    !> q-point mesh
    type(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot
    integer, intent(in) :: verbosity

    complex(r8), dimension(3) :: v0
    real(r8), dimension(:,:,:), allocatable :: tetenergies
    real(r8), dimension(:,:,:), allocatable :: tetsiteproj
    real(r8), dimension(4) :: tete,w1,w2
    real(r8) :: e,f0,deltae,sigma,t0,mine,maxe,totmaxe,totmine,invf
    integer :: i,j,k,l,ii,jj,ntet,li

    ! Miniscule smearing
    deltae=(pd%omega(2)-pd%omega(1))*0.5_r8

    ! Count number of tetrahedrons on this rank.
    ntet=0
    do i=1,qp%n_irr_tet
        if ( mod(i,mw%n) .ne. mw%r ) cycle
        ntet=ntet+1
    enddo

    ! Pre-fetch things for tetrahedrons
    if ( ntet .gt. 0 ) then
        call mem%allocate(tetenergies,[4,dr%n_mode,ntet],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
        call mem%allocate(tetsiteproj,[pd%n_atom,dr%n_mode,ntet],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
    else
        call mem%allocate(tetenergies,[1,1,1],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
        call mem%allocate(tetsiteproj,[1,1,1],persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
    endif
    tetenergies=0.0_r8
    tetsiteproj=0.0_r8
    li=0
    do i=1,qp%n_irr_tet
        if ( mod(i,mw%n) .ne. mw%r ) cycle
        li=li+1
        do j=1,4
            ! energies
            l=qp%it(i)%irreducible_index(j)
            do k=1,dr%n_mode
                tetenergies(j,k,li)=dr%iq(l)%omega(k)
            enddo
            ! projections
            do ii=1,pd%n_atom
            do jj=1,dr%n_mode
                v0=dr%iq(l)%egv((ii-1)*3+1:ii*3,jj)
                tetsiteproj(ii,jj,li)=tetsiteproj(ii,jj,li)+abs(dot_product(v0,v0))
            enddo
            enddo
        enddo
    enddo

    ! ! Do the actual integration
    totmaxe=maxval(pd%omega)
    totmine=minval(pd%omega)
    invf=(pd%n_dos_point/(totmaxe-totmine))
    !foursigma=4*sigma
    t0=walltime()
    ! Start by looping over tetrahedrons, I believe is smart.
    pd%pdos_site=0.0_r8
    pd%pdos_mode=0.0_r8
    if ( verbosity .gt. 0 ) call lo_progressbar_init()
    li=0
    do i=1,qp%n_irr_tet
        if ( mod(i,mw%n) .ne. mw%r ) cycle
        li=li+1
        do j=1,dr%n_mode
            ! smearing for degenerate tetrahedrons
            sigma=dr%default_smearing(j)
            ! smallest and largest tetrahedron energies, expressed as indices in the energy-axis
            tete=tetenergies(:,j,li)
            mine=minval(tete)-totmine
            maxe=maxval(tete)-totmine
            ii=max(floor( mine*invf ),1)
            jj=min(ceiling( maxe*invf ),pd%n_dos_point)
            ! loop over energy
            do k=ii,jj
                ! integration weights
                e=pd%omega(k)
                w1=lo_LV_tetrahedron_weights(tete,e-deltae,lo_freqtol,sigma)
                w2=lo_LV_tetrahedron_weights(tete,e+deltae,lo_freqtol,sigma)
                f0=sum(w1+w2)*qp%it(i)%integration_weight*0.5_r8
                ! Increment things
                pd%pdos_mode(k,j)=pd%pdos_mode(k,j)+f0
                do l=1,pd%n_atom
                    pd%pdos_site(k,l)=pd%pdos_site(k,l)+tetsiteproj(l,j,li)*f0
                enddo
            enddo
        enddo
        if ( verbosity .gt. 0 ) then
            if ( lo_trueNtimes(li,63,ntet) ) call lo_progressbar(' ... integrating phonon DOS',li,ntet)
        endif
    enddo
    ! sum them up over ranks
    call mw%allreduce('sum',pd%pdos_mode)
    call mw%allreduce('sum',pd%pdos_site)
    ! Cleanup
    call mem%deallocate(tetenergies,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
    call mem%deallocate(tetsiteproj,persistent=.false.,scalable=.true.,file=__FILE__,line=__LINE__)
    ! finalize the message
    if ( verbosity .gt. 0 ) then
        call lo_progressbar(' ... integrating phonon DOS',ntet,ntet,walltime()-t0)
    endif
end subroutine

!> Calculate the phonon dos with Gaussian smearing, either fixed or adaptive
module subroutine lo_get_phonon_dos_gaussian(pd,qp,dr,mw,verbosity)
    !> phonon dos
    type(lo_phonon_dos), intent(inout) :: pd
    !> q-point grid
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> mpi communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity

    complex(r8), dimension(3) :: v0
    integer, parameter :: nsigma=4
    real(r8), parameter :: sigmacorr=1.000063346496191_r8
    real(r8), dimension(:), allocatable :: dy
    real(r8), dimension(pd%n_atom) :: siteweight
    real(r8) :: f0,f1,sigma,weight,t0,totmaxe,totmine,invf,foursigma
    integer :: j,k,l,ii,jj,q
    integer :: ctr

    t0=walltime()
    ! some shortand energies
    totmaxe=maxval(pd%omega)
    totmine=minval(pd%omega)
    invf=(pd%n_dos_point/(totmaxe-totmine))
    pd%pdos_site=0.0_r8
    pd%pdos_mode=0.0_r8

    allocate(dy(pd%n_dos_point))
    ! Integration
    ctr=0
    if ( verbosity .gt. 1 ) call lo_progressbar_init()
    do q=1,qp%n_irr_point
        weight=qp%ip(q)%integration_weight
        do j=1,dr%n_mode
            ctr=ctr+1
            ! Make it parallel
            if ( mod(ctr,mw%n) .ne. mw%r ) cycle

            ! Get the current frequency
            f1=dr%iq(q)%omega(j)
            if ( f1 .lt. dr%omega_min*0.5_r8 ) cycle

            ! Get the smearing parameter
            select case(pd%integrationtype)
                case(1) ! fix gaussian
                    sigma=dr%default_smearing(j)*pd%smearing_prefactor
                case(2) ! adaptive gaussian
                    sigma=qp%adaptive_sigma( qp%ip(q)%radius,dr%iq(q)%vel(:,j),dr%default_smearing(j),pd%smearing_prefactor )
                    !sigma=qp%smearingparameter(dr%iq(q)%vel(:,j),dr%default_smearing(j),pd%smearing_prefactor)
            end select
            foursigma=sigma*nsigma

            ! Fetch projection thingy
            do ii=1,pd%n_atom
                v0=dr%iq(q)%egv( (ii-1)*3+1:ii*3, j )
                siteweight(ii)=abs(dot_product(v0,conjg(v0)))
            enddo

            ! range of integration
            ii=max( floor( (f1-totmine-foursigma)*invf ), 1)
            jj=min( floor( (f1-totmine+foursigma)*invf ), pd%n_dos_point)
            ! do k=ii,jj
            !     f0=lo_gauss(f1,pd%omega(k),sigma)*weight*sigmacorr
            !     pd%pdos_mode(k,j)=pd%pdos_mode(k,j)+f0
            !     do l=1,pd%n_atom
            !         pd%pdos_site(k,l)=pd%pdos_site(k,l)+f0*siteweight(l)
            !     enddo
            ! enddo

            ! how about being really anal and make sure we are actually normalized every single time?
            dy=0.0_r8
            do k=ii,jj
                dy(k)=lo_gauss(f1,pd%omega(k),sigma)
            enddo
            pd%pdos_mode(ii:jj,j)=pd%pdos_mode(ii:jj,j)+dy(ii:jj)*weight !*sigmacorr
            do l=1,pd%n_atom
                pd%pdos_site(ii:jj,l)=pd%pdos_site(ii:jj,l)+dy(ii:jj)*siteweight(l)
            enddo

        enddo
        if ( verbosity .gt. 0 ) then
            if ( lo_trueNtimes(q,63,qp%n_irr_point) ) call lo_progressbar(' ... integrating phonon DOS',q,qp%n_irr_point)
        endif
    enddo
    ! sum them up over ranks
    deallocate(dy)
    call mw%allreduce('sum',pd%pdos_site)
    call mw%allreduce('sum',pd%pdos_mode)

    if ( verbosity .gt. 0 ) then
        call lo_progressbar(' ... integrating phonon DOS',qp%n_irr_point,qp%n_irr_point,walltime()-t0)
    endif
end subroutine

! some old thing I got from someone. Not sure about it.
real(r8) function ARCSenergyresolution(omega)
    real(r8), intent(in) :: omega
    real(r8) :: E,f0,a0,a1,a2,a3,a4

    a0 =  1.8821603859_r8
    a1 = -3.8232422875e-2_r8
    a2 =  1.7805244237e-4_r8
    a3 =  1.1792351443e-6_r8
    a4 =  1.2642974944e-8_r8
    E=omega*lo_frequency_hartree_to_meV

    f0=(a0+a1*E+a2*E**2+a3*E**3+a4*E**4)/2.35_r8
    ARCSenergyresolution=f0/lo_frequency_hartree_to_meV
end function

end submodule
