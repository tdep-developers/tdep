submodule(phonondamping) phonondamping_aux
!! Helper routines that did not fit elsewhere. Maybe deprecated.
implicit none
contains

! !> The magical thermal prefactor
! function thermal_pref(bigQ,eigenvector,uc,temperature,omega,sigma) result(pref)
!     real(r8), dimension(3), intent(in) :: bigQ
!     complex(r8), dimension(:), intent(in) :: eigenvector
!     type(lo_crystalstructure), intent(in) :: uc
!     real(r8), intent(in) :: temperature
!     real(r8), intent(in) :: omega
!     real(r8), dimension(:,:,:), intent(in) :: sigma
!     real(r8) :: pref
!
!     integer :: i
!     complex(r8) :: c0,c1,c2
!     real(r8) :: f1,xs,dw
!
!     pref=1.0_r8
!     ! Add the thermal factor
!     pref=pref*(2*lo_planck(temperature,omega)+1.0_r8)
!
!     ! Cross-product guy
!     c2=0.0_r8
!     do i=1,uc%na
!         ! Grab the cross-section?
!         xs=uc%inelastic_neutron_cross_section(i)*uc%invsqrtmass(i)
!         ! Debye-Waller factor?
!         dw=exp( -0.5_r8*dot_product(bigQ,matmul(sigma(:,:,i),bigQ))*lo_twopi**2 )
!
!         f1=dot_product(lo_twopi*bigQ,uc%rcart(:,i))
!         c0=cmplx( cos(f1), sin(f1) ,r8)
!         c1=dot_product( lo_twopi*bigQ,eigenvector( (i-1)*3+1:i*3 ) )
!         c2=c2+c0*c1*xs
!     enddo
!     pref=pref*abs(conjg(c2)*c2)
! end function
!
! subroutine spectral_tau(omega,spectralfunction,temperature,tau)
!     !> energy axis
!     real(r8), dimension(:), intent(in) :: omega
!     !> spectral function
!     real(r8), dimension(:), intent(in) :: spectralfunction
!     !> temperature
!     real(r8), intent(in) :: temperature
!     !> thermal conductivity?
!     real(r8), intent(out) :: tau
!
!     real(r8) :: n
!     integer :: ie
!
!     tau=0.0_r8
!     do ie=2,size(omega)
!         n=lo_planck(temperature,omega(ie))
!         tau=tau+n*(n+1)*(spectralfunction(ie)**2)
!     enddo
!     tau=tau*(omega(2)-omega(1))
! end subroutine

end submodule
