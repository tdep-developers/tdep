#include "precompilerdefinitions"
module type_forceconstant_alloy_secondorder
!> @todo this must be updated or retired, all of it

!    use constants
!    use gottochblandat
!    use type_crystalstructure
!    use type_forceconstant_secondorder
!    implicit none
!
!    private
!    public :: lo_forceconstant_alloy_secondorder
!    !
!    !> Forceconstants for substitutional alloys
!    type lo_forceconstant_alloy_secondorder
!        !> how many different kind of pairs
!        integer :: n
!        !> what are the pairs?
!        character(len=40), dimension(:), allocatable :: pairsymbol
!        !> the pair resolved force constants
!        type(lo_forceconstant_secondorder), dimension(:), allocatable :: fc
!        !> the effective force constants
!        type(lo_forceconstant_secondorder) :: fc_eff
!        !> how much to talk
!        integer :: verbosity
!        contains
!            procedure :: unallocate
!            procedure :: writetofile
!            procedure :: readfromfile
!            procedure :: remap
!    end type
contains

! subroutine remap(afc,fcss,uc,ss,sqs,ssorder)
!     class(lo_forceconstant_alloy_secondorder), intent(in) :: afc
!     type(lo_forceconstant_secondorder), intent(out) :: fcss
!     type(lo_crystalstructure), intent(in) :: uc,ss,sqs
!     logical, intent(in), optional :: ssorder
!     !
!     type(lo_neighbourlist) :: nn
!     type(lo_crystalstructure) :: p,ps
!     real(flyt), dimension(3) :: v1,v2,lv1,lv2,ucv1,ucv2,rv_frac
!     real(flyt) :: f0,f1
!     character(len=40), dimension(2) :: dumpar
!     character(len=6) :: dumstr
!     integer, dimension(:), allocatable :: dumind
!     integer :: i,j,k,l,a1,a2,uca,ii,jj
!     logical :: sso
!     !
!     ! Order of the atoms, very confusing
!     !
!     if ( present(ssorder) ) then
!         sso=ssorder
!     else
!         sso=.false.
!     endif
!     !
!     ! Make a copy of the alloy supercell, change the r-list to that of the SQS
!     !
!     if ( sso .eqv. .false. ) then
!         p=ss
!         ps=sqs
!     else
!         p=ss
!         ps=sqs
!         lo_allocate(dumind(p%na))
!         do i=1,p%na
!             f0=1000
!             do j=1,sqs%na
!                 f1=lo_ldist(p%r(:,i),ps%r(:,j))
!                 if ( f1 .lt. f0 ) then
!                     dumind(i)=j
!                     f0=f1
!                 endif
!             enddo
!         enddo
!         ps%r=ps%r(:,dumind)
!         ps%species=ps%species(dumind)
!     endif
!     !
!     ! This should be a good start. Now I want neighbourlists and mappings.
!     !
!     call nn%generate(uc,afc%fc_eff%cutoff*1.3_flyt,p,matchcells=.true.)
!     !
!     fcss%cutoff=afc%fc_eff%cutoff
!     fcss%na=p%na
!     lo_allocate(fcss%atom(p%na))
!     do i=1,p%na
!         uca=nn%ssatom(i)%atom_in_uc
!         fcss%atom(i)%n=afc%fc_eff%atom(uca)%n
!         lo_allocate(fcss%atom(i)%pair( fcss%atom(i)%n ))
!     enddo
!     ! start matching
!
! !    !$OMP PARALLEL DEFAULT(private) SHARED(nn,fcss,p,ps,afc)
! !    !$OMP DO
!     do i=1,p%na
!         !
!         ! Wich atom in the unit cell is it?
!         !
!         uca=nn%ssatom(i)%atom_in_uc
!         !
!         l=0
!         do j=1,afc%fc_eff%atom(uca)%n
!             do k=1,nn%ssatom(i)%n
!                 if ( norm2(afc%fc_eff%atom(uca)%pair(j)%r-nn%ssatom(i)%v(:,k)) .lt. lo_tol ) then
!                     l=l+1
!                     !
!                     ! Now the vectors are matched, just store everything
!                     !
!                     a1=i
!                     a2=nn%ssatom(i)%ind(k)
!                     fcss%atom(i)%pair(j)%i1=a1
!                     fcss%atom(i)%pair(j)%i2=a2
!                     !
!                     ! All different vectors, first the vectors in the unit cell
!                     !
!                     ucv1=p%r(:,a1)
!                     ucv2=p%r(:,a2)
!                     !
!                     ! The absolute vector in fractional coordinates
!                     !
!                     rv_frac=matmul(afc%fc_eff%atom(uca)%pair(j)%r,p%invbas)
!                     !
!                     ! The lattice vectors
!                     !
!                     lv1=0.0_flyt
!                     lv2=nint(rv_frac+ucv1-ucv2)*1.0_flyt
!                     !
!                     ! The position vectors
!                     !
!                     v1=ucv1
!                     v2=ucv2+lv2
!                     !
!                     ! Now convert everything to cartesian coordinates
!                     !
!                     fcss%atom(i)%pair(j)%v1=p%fractional_to_cartesian(v1)
!                     fcss%atom(i)%pair(j)%v2=p%fractional_to_cartesian(v2)
!                     fcss%atom(i)%pair(j)%lv1=p%fractional_to_cartesian(lv1)
!                     fcss%atom(i)%pair(j)%lv2=p%fractional_to_cartesian(lv2)
!                     fcss%atom(i)%pair(j)%r=nn%ssatom(i)%v(:,k)
!                     fcss%atom(i)%pair(j)%m=afc%fc_eff%atom(uca)%pair(j)%m
!                     !
!                     ! I should add the proper force constant
!                     !
!                     dumpar(1)=ps%atomic_symbol(ps%species(a1))
!                     dumpar(2)=ps%atomic_symbol(ps%species(a2))
!                     call qsort(dumpar)
!                     dumstr=trim(dumpar(1))//trim(dumpar(2))
!                     jj=0
!                     do ii=1,afc%n
!                         if ( dumstr .eq. afc%pairsymbol(ii) ) then
!                             jj=ii
!                             exit
!                         endif
!                     enddo
!                     if ( jj .eq. 0 ) then
!                         write(*,*) 'Error in finding correct pair for alloy forceconstant'
!                         stop
!                     endif
!                     !
!                     fcss%atom(i)%pair(j)%m=afc%fc(jj)%atom(uca)%pair(j)%m
!                     !
!                 endif
!             enddo
!         enddo
!         !
!         ! Check if it went alright
!         !
!         if ( l .ne. afc%fc_eff%atom(uca)%n ) then
!             write(*,*) 'FCREMAP: ERROR!!! FAILED MAPPING ATOM',i
!             stop
!         endif
!         !
!     enddo
! !    !$OMP END DO
! !    !$OMP END PARALLEL
!     !
!     ! Make sure the sum is 0
!     !
!     call fcss%setsumtozero()
!     !
! end subroutine
!
! subroutine writetofile(afc,uc)
!     ! dump it to file
!     class(lo_forceconstant_alloy_secondorder), intent(in) :: afc
!     type(lo_crystalstructure), intent(in) :: uc
!     !
!     integer :: i,u
!     !
!     u=open_file('out','outfile.forceconstant_alloy')
!         write(u,*) afc%n
!         do i=1,afc%n
!             write(u,*) trim(afc%pairsymbol(i))
!         enddo
!         do i=1,afc%n
!             call afc%fc(i)%writetofile(uc,enhet=u)
!         enddo
!         call afc%fc_eff%writetofile(uc,enhet=u)
!     close(u)
! end subroutine
!
! subroutine readfromfile(afc,uc)
!     ! dump it to file
!     class(lo_forceconstant_alloy_secondorder), intent(out) :: afc
!     type(lo_crystalstructure), intent(in) :: uc
!     !
!     integer :: i,u
!     !
!     u=open_file('in','infile.forceconstant_alloy')
!         read(u,*) afc%n
!        lo_allocate(afc%pairsymbol(afc%n))
!        lo_allocate(afc%fc(afc%n))
!         do i=1,afc%n
!             read(u,*) afc%pairsymbol(i)
!         enddo
!         do i=1,afc%n
!             call afc%fc(i)%readfromfile(uc,enhet=u)
!         enddo
!         call afc%fc_eff%readfromfile(uc,enhet=u)
!     close(u)
! end subroutine
!
! subroutine unallocate(afc)
!     ! deallocate the structure
!     class(lo_forceconstant_alloy_secondorder), intent(inout) :: afc
!     !
!     integer :: i
!     if ( allocated(afc%pairsymbol) )lo_deallocate(afc%pairsymbol)
!     if ( allocated(afc%fc) ) then
!         do i=1,size(afc%fc,1)
!             call afc%fc(i)%unallocate()
!         enddo
!     endif
! end subroutine

end module
