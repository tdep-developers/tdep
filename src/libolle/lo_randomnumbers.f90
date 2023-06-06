module lo_randomnumbers
!-------------------------------------------------------------------------------
!   This is a Fortran translation of the 64-bit version of
!   the Mersenne Twister pseudorandom number generator
!
!   Before using, initialize the state by using
!       call init_genrand64(seed)
!   or
!       call init_by_array64(init_key)
!
!   Translated from C-program for MT19937-64 (2004/9/29 version)
!   originally coded by Takuji Nishimura and Makoto Matsumoto
!   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
!
!   Fortran translation by RÃ©mi Piatek
!   The University of Copenhagen
!   Department of Economics
!   email: {first}.{last}@econ.ku.dk
!
!-------------------------------------------------------------------------------
!   A C-program for MT19937-64 (2004/9/29 version).
!   Coded by Takuji Nishimura and Makoto Matsumoto.
!
!   This is a 64-bit version of Mersenne Twister pseudorandom number
!   generator.
!
!   Before using, initialize the state by using init_genrand64(seed)
!   or init_by_array64(init_key, key_length).
!
!   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions
!   are met:
!
!     1. Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!
!     2. Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!
!     3. The names of its contributors may not be used to endorse or promote
!        products derived from this software without specific prior written
!        permission.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!   References:
!   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
!     ACM Transactions on Modeling and
!     Computer Simulation 10. (2000) 348--357.
!   M. Matsumoto and T. Nishimura,
!     ``Mersenne Twister: a 623-dimensionally equidistributed
!       uniform pseudorandom number generator''
!     ACM Transactions on Modeling and
!     Computer Simulation 8. (Jan. 1998) 3--30.
!
!   Any feedback is very welcome.
!   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
!-------------------------------------------------------------------------------

!
! OLLES NOTES:
! I just wrapped it into an easily passed-around container, and added some stupid things such
! as getting normal random numbers, spherically random unit vectors and stuff like that.
!
use konstanter, only: r8,i8,lo_twopi
implicit none

private
public :: lo_mersennetwister

integer(i8), parameter :: nn       = 312_i8
integer(i8), parameter :: mm       = 156_i8
integer(i8), parameter :: seed_def = 5489_i8
integer(i8), parameter :: matrix_a = -5403634167711393303_i8
integer(i8), parameter :: um       = -2147483648_i8 ! most significant 33 bits
integer(i8), parameter :: lm       =  2147483647_i8 ! least significant 31 bits

real(r8),    parameter :: pi253_1  = 1._r8/(2._r8**53 - 1._r8)
real(r8),    parameter :: pi253    = 1._r8/(2._r8**53)
real(r8),    parameter :: pi252    = 1._r8/(2._r8**52)

!> wrapper that is easily passed around.
type lo_mersennetwister
    ! Is it initialized
    logical :: initialized=.false.
    ! Magical things I don't understand, but probably should be instanced.
    integer(i8), private :: mt(nn)       ! array for the state vector
    integer,     private :: mti = nn+1   ! mti==nn+1 means mt(nn) is not initialized
    contains
        !> seed the random numbers
        procedure :: init=>init_genrand64
        !> get a random floating point thing, 0-1
        procedure :: rnd_real=>genrand64_real1
        !> get a spherically random unit vector
        procedure :: rnd_unitvector=>random_unit_vector
        !> get a random integer in the interval 1-n
        procedure :: rnd_int=>random_integer
        !> get a single gaussian random number
        procedure :: rnd_gaussian_real=>random_gaussian
        !> get a Box-Muller pair
        procedure :: rnd_boxmuller_pair=>boxmuller
        !> shuffle an integer array
        procedure :: shuffle_int_array
end type

contains

!> Randomly sorts integer array
subroutine shuffle_int_array(tw,a)
    class(lo_mersennetwister), intent(inout) :: tw
    !> the array to be shuffled
    integer, dimension(:), intent(inout) :: a

    integer :: i,j,k,n

    n=size(a,1)
    do i=n,1,-1
        j=floor(genrand64_real1(tw)*n)+1
        k=a(i)
        a(i)=a(j)
        a(j)=k
    enddo
end subroutine

!> A single random gaussian number. Will calculate an extra one that I throw away, but that's life or something.
function random_gaussian(tw,mu,sigma) result(y)
    class(lo_mersennetwister), intent(inout) :: tw
    !> mean of distribution
    real(r8), intent(in) :: mu
    !> standard deviation of distribution
    real(r8), intent(in) :: sigma
    !> Gaussian random number
    real(r8) :: y

    real(r8) :: f1,f2,z1,z2
    f1=genrand64_real3(tw)
    f2=genrand64_real1(tw)
    z1 = sqrt( -2*log(f1) )*sin(lo_twopi*f2);
    z2 = sqrt( -2*log(f1) )*cos(lo_twopi*f2);
    y=mu+z1*sigma
end function

!> Box-Muller transform, return both values
subroutine boxmuller(tw,sigma,mu,x1,x2)
    class(lo_mersennetwister), intent(inout) :: tw
    !> mean of distribution
    real(r8), intent(in) :: mu
    !> standard deviation of distribution
    real(r8), intent(in) :: sigma
    !> the box-muller pair
    real(r8), intent(out) :: x1,x2

    real(r8) :: f1,f2,z1,z2
    f1=genrand64_real3(tw)
    f2=genrand64_real1(tw)
    z1 = sqrt( -2*log(f1) )*sin(lo_twopi*f2);
    z2 = sqrt( -2*log(f1) )*cos(lo_twopi*f2);
    x1=mu+z1*sigma
    x2=mu+z2*sigma
end subroutine

!> random integer, in the interval 1-n
function random_integer(tw,n) result(i)
    class(lo_mersennetwister), intent(inout) :: tw
    integer, intent(in) :: n
    integer :: i
    i=floor(genrand64_real1(tw)*n)+1
end function

!> get a spherically random unit vector. Not sure if this is faster than using Gaussian random numbers.
function random_unit_vector(tw) result(v)
    class(lo_mersennetwister), intent(inout) :: tw
    real(r8), dimension(3) :: v

    real(r8) :: r1,r2,r3

    do
        ! using the guy that can't become zero here to make sure I never
        ! divide by zero later. Rather low odds but you never know.
        r1=(genrand64_real3(tw)-0.5_r8)*2.0_r8
        r2=(genrand64_real3(tw)-0.5_r8)*2.0_r8
        r3=(genrand64_real3(tw)-0.5_r8)*2.0_r8
        if ( r1*r1+r2*r2+r3*r3 .le. 1.0_r8 ) exit
    enddo
    v=[r1,r2,r3]/norm2([r1,r2,r3])
end function

!> Initializes mt(nn), with an optional seed. Otherwise seeded from the walltime.
subroutine init_genrand64(tw,iseed,rseed)
    !> rng container
    class(lo_mersennetwister), intent(out) :: tw
    !> integer to seed with (usually the MPI rank)
    integer, intent(in) :: iseed
    !> floating point thing to seed with (usually the Walltime)
    real(r8), intent(in) :: rseed

    real(r8) :: f0
    integer(i8) :: seed
    integer :: i

    ! Seed with the walltime? This gives me something -0.5,0.5, Kinda
    f0=rseed-floor(rseed)-0.5_r8
    ! Add an integer seed here. For example the mpi rank, to get them to
    ! start in different places.
    if ( iseed .gt. 0 ) then
        f0=mod(f0+log(real(iseed,r8)),1.0_r8)-0.5_r8
    elseif ( iseed .lt. 0 ) then
        f0=mod(f0-log(real(-iseed,r8)),1.0_r8)-0.5_r8
    endif
    ! Get it to an integer. This would usually be rather random, I hope.
    seed=int(anint(f0*(2**30)),i8)

    tw%mt(1) = seed
    do i = 1,nn-1
        tw%mt(i+1) = 6364136223846793005_i8 * ieor(tw%mt(i), ishft(tw%mt(i), -62)) + i
    end do
    tw%mti = nn

    ! And report that we are initialized
    tw%initialized=.true.
end subroutine init_genrand64

!> Main thingy, generates a random 64-bit integer on [-2^63, 2^63-1]-interval
integer(i8) function genrand64_int64(tw)
    class(lo_mersennetwister), intent(inout) :: tw
    integer(i8) :: mag01(0:1) = [0_i8, matrix_a]
    integer(i8) :: x
    integer     :: i

    if( tw%mti >= nn) then ! generate nn words at one time
        ! if init_genrand64() has not been called, a default initial seed is used
        !if(mti == nn+1) call init_genrand64(seed_def)

        do i = 1, nn-mm
            x = ior(iand(tw%mt(i),um), iand(tw%mt(i+1), lm))
            tw%mt(i) = ieor(ieor(tw%mt(i+mm), ishft(x, -1)), mag01(iand(x, 1_i8)))
        end do

        do i = nn-mm+1, nn-1
            x = ior(iand(tw%mt(i), um), iand(tw%mt(i+1), lm))
            tw%mt(i) = ieor(ieor(tw%mt(i+mm-nn), ishft(x, -1)), mag01(iand(x, 1_i8)))
        end do

        x = ior(iand(tw%mt(nn), um), iand(tw%mt(1), lm))
        tw%mt(nn) = ieor(ieor(tw%mt(mm), ishft(x, -1)), mag01(iand(x, 1_i8)))
        tw%mti = 0
    end if

    tw%mti = tw%mti + 1
    x = tw%mt(tw%mti)
    x = ieor(x, iand(ishft(x,-29), 6148914691236517205_i8))
    x = ieor(x, iand(ishft(x, 17), 8202884508482404352_i8))
    x = ieor(x, iand(ishft(x, 37),   -2270628950310912_i8))
    x = ieor(x, ishft(x, -43))
    genrand64_int64 = x
end function genrand64_int64

!> Generates a random number on [0,1]-real-interval
real(r8) function genrand64_real1(tw)
    class(lo_mersennetwister), intent(inout) :: tw
    genrand64_real1 = real(ishft(genrand64_int64(tw), -11), kind=r8) * pi253_1
end function genrand64_real1

!> Generates a random number on [0,1)-real-interval
real(r8) function genrand64_real2(tw)
    class(lo_mersennetwister), intent(inout) :: tw
    genrand64_real2 = real(ishft(genrand64_int64(tw), -11), kind=r8) * pi253
end function genrand64_real2

!> Generates a random number on (0,1)-real-interval
real(r8) function genrand64_real3(tw)
    class(lo_mersennetwister), intent(inout) :: tw
    genrand64_real3 = real(ishft(genrand64_int64(tw), -12), kind=r8)
    genrand64_real3 = (genrand64_real3 + 0.5_r8) * pi252
end function genrand64_real3

! This one I don't really understand and that makes me mad.
!  !-----------------------------------------------------------------------------
!  ! Initializes by an array with array-length
!  !   init_key is the array for initializing keys
!  subroutine init_by_array64(init_key)
!    implicit none
!    integer(i8), intent(in) :: init_key(:)
!    integer(i8), parameter  :: c1 = 3935559000370003845_i8
!    integer(i8), parameter  :: c2 = 2862933555777941757_i8
!    integer(i8) :: i, j, k, kk, key_length
!
!    call init_genrand64(19650218_i8)
!    key_length = size(init_key)
!    i = 1_i8; j = 0_i8
!    k = max(nn, key_length)
!
!    do kk = 1, k
!      mt(i+1) = ieor(mt(i+1), c1 * ieor(mt(i), ishft(mt(i), -62))) &
!                  + init_key(j+1) + j
!      i = i+1; j = j+1
!      if(i >= nn) then
!        mt(1) = mt(nn)
!        i = 1
!      end if
!      if(j >= key_length) j = 0
!    end do
!
!    do kk = 1, nn-1
!      mt(i+1) = ieor(mt(i+1), c2 * ieor(mt(i), ishft(mt(i), -62))) - i
!      i = i+1
!      if(i >= nn) then
!        mt(1) = mt(nn)
!        i = 1
!      end if
!    end do
!
!    mt(1) = ishft(1_i8, 63)  ! MSB is 1; assuring non-zero initial array
!
!  end subroutine init_by_array64

end module
