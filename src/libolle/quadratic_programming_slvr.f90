!c
!c  Copyright (C) 1995-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com>
!c
!c  This program is free software; you can redistribute it and/or modify
!c  it under the terms of the GNU General Public License as published by
!c  the Free Software Foundation; either version 2 of the License, or
!c  (at your option) any later version.
!c
!c  This program is distributed in the hope that it will be useful,
!c  but WITHOUT ANY WARRANTY; without even the implied warranty of
!c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!c  GNU General Public License for more details.
!c
!c  You should have received a copy of the GNU General Public License
!c  along with this program; if not, write to the Free Software
!c  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
!c  USA.
!c
!c  this routine uses the Goldfarb/Idnani algorithm to solve the
!c  following minimization problem:
!c
!c        minimize  -d^T x + 1/2 *  x^T D x
!c        where   A1^T x  = b1
!c                A2^T x >= b2
!c
!c  the matrix D is assumed to be positive definite.  Especially,
!c  w.l.o.g. D is assumed to be symmetric.
!c
!c  Input parameter:
!c  dmat   nxn matrix, the matrix D from above (dp)
!c         *** WILL BE DESTROYED ON EXIT ***
!c         The user has two possibilities:
!c         a) Give D (ierr=0), in this case we use routines from LINPACK
!c            to decompose D.
!c         b) To get the algorithm started we need R^-1, where D=R^TR.
!c            So if it is cheaper to calculate R^-1 in another way (D may
!c            be a band matrix) then with the general routine, the user
!c            may pass R^{-1}.  Indicated by ierr not equal to zero.
!c  dvec   nx1 vector, the vector d from above (dp)
!c         *** WILL BE DESTROYED ON EXIT ***
!c         contains on exit the solution to the initial, i.e.,
!c         unconstrained problem
!c  fddmat scalar, the leading dimension of the matrix dmat
!c  n      the dimension of dmat and dvec (int)
!c  amat   lxq matrix (dp)
!c         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!c             CHANGED SIGNES ON EXIT ***
!c  iamat  (l+1)xq matrix (int)
!c         these two matrices store the matrix A in compact form. the format
!c         is: [ A=(A1 A2)^T ]
!c           iamat(1,i) is the number of non-zero elements in column i of A
!c           iamat(k,i) for k>=2, is equal to j if the (k-1)-th non-zero
!c                      element in column i of A is A(i,j)
!c            amat(k,i) for k>=1, is equal to the k-th non-zero element
!c                      in column i of A.
!c
!c  bvec   qx1 vector, the vector of constants b in the constraints (dp)
!c         [ b = (b1^T b2^T)^T ]
!c         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!c             CHANGED SIGNES ON EXIT ***
!c  fdamat the first dimension of amat as declared in the calling program.
!c         fdamat >= n (and iamat must have fdamat+1 as first dimension)
!c  q      integer, the number of constraints.
!c  meq    integer, the number of equality constraints, 0 <= meq <= q.
!c  ierr   integer, code for the status of the matrix D:
!c            ierr =  0, we have to decompose D
!c            ierr != 0, D is already decomposed into D=R^TR and we were
!c                       given R^{-1}.
!c
!c  Output parameter:
!c  sol   nx1 the final solution (x in the notation above)
!c  lagr  qx1 the final Lagrange multipliers
!c  crval scalar, the value of the criterion at the minimum
!c  iact  qx1 vector, the constraints which are active in the final
!c        fit (int)
!c  nact  scalar, the number of constraints active in the final fit (int)
!c  iter  2x1 vector, first component gives the number of "main"
!c        iterations, the second one says how many constraints were
!c        deleted after they became active
!c  ierr  integer, error code on exit, if
!c           ierr = 0, no problems
!c           ierr = 1, the minimization problem has no solution
!c           ierr = 2, problems with decomposing D, in this case sol
!c                     contains garbage!!
!c
!c  Working space:
!c  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1
!c        where r=min(n,q)
!c
! OH: I just rewrote the routine to modern fortran
submodule (quadratic_programming) quadratic_programming_slvr

implicit none
contains

module subroutine qpgen1(dmat, dvec, fddmat, n, sol, lagr, crval, amat, iamat,bvec, q, meq, iact, nact, iter, work, ierr)
    real(flyt), dimension(:,:), intent(inout) :: dmat
    real(flyt), dimension(:), intent(inout) :: dvec
    integer, intent(in) :: fddmat
    integer, intent(in) :: n
    real(flyt), dimension(:), intent(out) :: sol
    real(flyt), dimension(:), intent(out) :: lagr
    real(flyt), intent(inout) :: crval
    real(flyt), dimension(:,:), intent(inout) :: amat
    integer, dimension(:,:), intent(in) :: iamat
    real(flyt), dimension(:), intent(inout) :: bvec
    !integer, intent(in) :: fdamat
    integer, intent(in) :: q
    integer, intent(in) :: meq
    integer, dimension(:), intent(out) :: iact
    integer, intent(inout) :: nact
    integer, dimension(2), intent(out) :: iter
    real(flyt), dimension(:), intent(inout) :: work
    integer, intent(inout) :: ierr

    real(flyt) :: temp,sum,gc,gs,nu,t1,tmpa,tmpb,tt,vsmall
    integer :: i, j, l, l1, info, it1, iwnbv, iwrm, iwrv, iwsv, iwuv, iwzv, nvl, r
    logical t1inf, t2min

      r = min(n,q)
      l = 2*n + (r*(r+5))/2 + 2*q + 1
!c
!c     code gleaned from Powell's ZQPCVX routine to determine a small
!c     number  that can be assumed to be an upper bound on the relative
!c     precision of the computer arithmetic.
!c
      vsmall = 1.0d-60
 1    vsmall = vsmall + vsmall
      tmpa = 1.0d0 + 0.1d0*vsmall
      tmpb = 1.0d0 + 0.2d0*vsmall
      if( tmpa .LE. 1.0d0 ) goto 1
      if( tmpb .LE. 1.0d0 ) goto 1
!c
!c store the initial dvec to calculate below the unconstrained minima of
!c the critical value.
!c
      do 10 i=1,n
         work(i) = dvec(i)
 10   continue
      do 11 i=n+1,l
         work(i) = 0.d0
 11   continue
      do 12 i=1,q
         iact(i) = 0
         lagr(i) = 0.d0
 12   continue
!c
!c get the initial solution
!c
      if( ierr .EQ. 0 )then
         call dpofa(dmat,size(dmat,2),info)
         if( info .NE. 0 )then
            ierr = 2
            goto 999
         endif
         call dposl(dmat,fddmat,n,dvec)
         call dpori(dmat,fddmat,n)
      else
!c
!c Matrix D is already factorized, so we have to multiply d first with
!c R^-T and then with R^-1.  R^-1 is stored in the upper half of the
!c array dmat.
!c
         do 20 j=1,n
            sol(j)  = 0.d0
            do 21 i=1,j
               sol(j) = sol(j) + dmat(i,j)*dvec(i)
 21         continue
 20      continue
         do 22 j=1,n
            dvec(j) = 0.d0
            do 23 i=j,n
               dvec(j) = dvec(j) + dmat(j,i)*sol(i)
 23         continue
 22      continue
      endif
!c
!c set lower triangular of dmat to zero, store dvec in sol and
!c calculate value of the criterion at unconstrained minima
!c
      crval = 0.d0
      do 30 j=1,n
         sol(j)  = dvec(j)
         crval   = crval + work(j)*sol(j)
         work(j) = 0.d0
         do 32 i=j+1,n
            dmat(i,j) = 0.d0
 32      continue
 30   continue
      crval = -crval/2.d0
      ierr  = 0
!c
!c calculate some constants, i.e., from which index on the different
!c quantities are stored in the work matrix
!c
      iwzv  = n
      iwrv  = iwzv + n
      iwuv  = iwrv + r
      iwrm  = iwuv + r+1
      iwsv  = iwrm + (r*(r+1))/2
      iwnbv = iwsv + q
!c
!c calculate the norm of each column of the A matrix
!c
      do 51 i=1,q
         sum = 0.d0
         do 52 j=1,iamat(1,i)
            sum = sum + amat(j,i)*amat(j,i)
 52      continue
         work(iwnbv+i) = sqrt(sum)
 51   continue
      nact = 0
      iter(1) = 0
      iter(2) = 0
 50   continue
!c
!c start a new iteration
!c
      iter(1) = iter(1)+1
!c
!c calculate all constraints and check which are still violated
!c for the equality constraints we have to check whether the normal
!c vector has to be negated (as well as bvec in that case)
!c
      l = iwsv
      do 60 i=1,q
         l = l+1
         sum = -bvec(i)
         do 61 j = 1,iamat(1,i)
            sum = sum + amat(j,i)*sol(iamat(j+1,i))
 61      continue
         if ( abs(sum) .LT. vsmall ) then
            sum = 0.0d0
         endif
         if (i .GT. meq) then
            work(l) = sum
         else
            work(l) = -abs(sum)
            if (sum .GT. 0.d0) then
               do 62 j=1,iamat(1,i)
                  amat(j,i) = -amat(j,i)
 62            continue
               bvec(i) = -bvec(i)
            endif
         endif
 60   continue
!c
!c as safeguard against rounding errors set already active constraints
!c explicitly to zero
!c
      do 70 i=1,nact
         work(iwsv+iact(i)) = 0.d0
 70   continue
!c
!c we weight each violation by the number of non-zero elements in the
!c corresponding row of A. then we choose the violated constraint which
!c has maximal absolute value, i.e., the minimum.
!c by obvious commenting and uncommenting we can choose the strategy to
!c take always the first constraint which is violated. ;-)
!c
      nvl = 0
      temp = 0.d0
      do 71 i=1,q
         if (work(iwsv+i) .LT. temp*work(iwnbv+i)) then
            nvl = i
            temp = work(iwsv+i)/work(iwnbv+i)
         endif
!c         if (work(iwsv+i) .LT. 0.d0) then
!c            nvl = i
!c            goto 72
!c         endif
 71   continue
      if (nvl .EQ. 0) then ! 72
         do 73 i=1,nact
            lagr(iact(i))=work(iwuv+i)
 73      continue
         goto 999
      endif
!c
!c calculate d=J^Tn^+ where n^+ is the normal vector of the violated
!c constraint. J is stored in dmat in this implementation!!
!c if we drop a constraint, we have to jump back here.
!c
 55   continue
      do 80 i=1,n
         sum = 0.d0
         do 81 j=1,iamat(1,nvl)
            sum = sum + dmat(iamat(j+1,nvl),i)*amat(j,nvl)
 81      continue
         work(i) = sum
 80   continue
!c
!c Now calculate z = J_2 d_2
!c
      l1 = iwzv
      do 90 i=1,n
         work(l1+i) =0.d0
 90   continue
      do 92 j=nact+1,n
         do 93 i=1,n
            work(l1+i) = work(l1+i) + dmat(i,j)*work(j)
 93      continue
 92   continue
!c
!c and r = R^{-1} d_1, check also if r has positive elements (among the
!c entries corresponding to inequalities constraints).
!c
      t1inf = .TRUE.
      do 95 i=nact,1,-1
         sum = work(i)
         l  = iwrm+(i*(i+3))/2
         l1 = l-i
         do 96 j=i+1,nact
            sum = sum - work(l)*work(iwrv+j)
            l   = l+j
 96      continue
         sum = sum / work(l1)
         work(iwrv+i) = sum
         if (iact(i) .LE. meq) goto 95
         if (sum .LE. 0.d0) goto 95
         t1inf = .FALSE.
         it1 = i
 95   continue
!c
!c if r has positive elements, find the partial step length t1, which is
!c the maximum step in dual space without violating dual feasibility.
!c it1  stores in which component t1, the min of u/r, occurs.
!c
      if ( .NOT. t1inf) then
         t1   = work(iwuv+it1)/work(iwrv+it1)
         do 100 i=1,nact
            if (iact(i) .LE. meq) goto 100
            if (work(iwrv+i) .LE. 0.d0) goto 100
            temp = work(iwuv+i)/work(iwrv+i)
            if (temp .LT. t1) then
               t1   = temp
               it1  = i
            endif
 100     continue
      endif
!c
!c test if the z vector is equal to zero
!c
      sum = 0.d0
      do 110 i=iwzv+1,iwzv+n
         sum = sum + work(i)*work(i)
 110  continue
      if (abs(sum) .LE. vsmall) then
!c
!c No step in primal space such that the new constraint becomes
!c feasible. Take step in dual space and drop a constant.
!c
         if (t1inf) then
!c
!c No step in dual space possible either, problem is not solvable
!c
            ierr = 1
            goto 999
         else
!c
!c we take a partial step in dual space and drop constraint it1,
!c that is, we drop the it1-th active constraint.
!c then we continue at step 2(a) (marked by label 55)
!c
            do 111 i=1,nact
               work(iwuv+i) = work(iwuv+i) - t1*work(iwrv+i)
 111        continue
            work(iwuv+nact+1) = work(iwuv+nact+1) + t1
            goto 700
         endif
      else
!c
!c compute full step length t2, minimum step in primal space such that
!c the constraint becomes feasible.
!c keep sum (which is z^Tn^+) to update crval below!
!c
         sum = 0.d0
         do 120 i = 1,iamat(1,nvl)
            sum = sum + work(iwzv+iamat(i+1,nvl))*amat(i,nvl)
 120     continue
         tt = -work(iwsv+nvl)/sum
         t2min = .TRUE.
         if (.NOT. t1inf) then
            if (t1 .LT. tt) then
               tt    = t1
               t2min = .FALSE.
            endif
         endif
!c
!c take step in primal and dual space
!c
         do 130 i=1,n
            sol(i) = sol(i) + tt*work(iwzv+i)
 130     continue
         crval = crval + tt*sum*(tt/2.d0 + work(iwuv+nact+1))
         do 131 i=1,nact
            work(iwuv+i) = work(iwuv+i) - tt*work(iwrv+i)
 131     continue
         work(iwuv+nact+1) = work(iwuv+nact+1) + tt
!c
!c if it was a full step, then we check wheter further constraints are
!c violated otherwise we can drop the current constraint and iterate once
!c more
         if(t2min) then
!c
!c we took a full step. Thus add constraint nvl to the list of active
!c constraints and update J and R
!c
            nact = nact + 1
            iact(nact) = nvl
!c
!c to update R we have to put the first nact-1 components of the d vector
!c into column (nact) of R
!c
            l = iwrm + ((nact-1)*nact)/2 + 1
            do 150 i=1,nact-1
               work(l) = work(i)
               l = l+1
 150        continue
!c
!c if now nact=n, then we just have to add the last element to the new
!c row of R.
!c Otherwise we use Givens transformations to turn the vector d(nact:n)
!c into a multiple of the first unit vector. That multiple goes into the
!c last element of the new row of R and J is accordingly updated by the
!c Givens transformations.
!c
            if (nact .EQ. n) then
               work(l) = work(n)
            else
               do 160 i=n,nact+1,-1
!c
!c we have to find the Givens rotation which will reduce the element
!c (l1) of d to zero.
!c if it is already zero we don't have to do anything, except of
!c decreasing l1
!c
                  if (work(i) .EQ. 0.d0) goto 160
                  gc   = max(abs(work(i-1)),abs(work(i)))
                  gs   = min(abs(work(i-1)),abs(work(i)))
                  temp = sign(gc*sqrt(1+gs*gs/(gc*gc)), work(i-1))
                  gc   = work(i-1)/temp
                  gs   = work(i)/temp
!c
!c The Givens rotation is done with the matrix (gc gs, gs -gc).
!c If gc is one, then element (i) of d is zero compared with element
!c (l1-1). Hence we don't have to do anything.
!c If gc is zero, then we just have to switch column (i) and column (i-1)
!c of J. Since we only switch columns in J, we have to be careful how we
!c update d depending on the sign of gs.
!c Otherwise we have to apply the Givens rotation to these columns.
!c The i-1 element of d has to be updated to temp.
!c
                  if (gc .EQ. 1.d0) goto 160
                  if (gc .EQ. 0.d0) then
                     work(i-1) = gs * temp
                     do 170 j=1,n
                        temp        = dmat(j,i-1)
                        dmat(j,i-1) = dmat(j,i)
                        dmat(j,i)   = temp
 170                 continue
                  else
                     work(i-1) = temp
                     nu = gs/(1.d0+gc)
                     do 180 j=1,n
                        temp        = gc*dmat(j,i-1) + gs*dmat(j,i)
                        dmat(j,i)   = nu*(dmat(j,i-1)+temp) - dmat(j,i)
                        dmat(j,i-1) = temp
 180                 continue
                  endif
 160           continue
!c
!c l is still pointing to element (nact,nact) of the matrix R.
!c So store d(nact) in R(nact,nact)
               work(l) = work(nact)
            endif
         else
!c
!c we took a partial step in dual space. Thus drop constraint it1,
!c that is, we drop the it1-th active constraint.
!c then we continue at step 2(a) (marked by label 55)
!c but since the fit changed, we have to recalculate now "how much"
!c the fit violates the chosen constraint now.
!c
            sum = -bvec(nvl)
            do 190 j = 1,iamat(1,nvl)
               sum = sum + sol(iamat(j+1,nvl))*amat(j,nvl)
 190        continue
            if( nvl .GT. meq ) then
               work(iwsv+nvl) = sum
            else
               work(iwsv+nvl) = -abs(sum)
               if( sum .GT. 0.d0) then
                  do 191 j=1,iamat(1,nvl)
                     amat(j,nvl) = -amat(j,nvl)
 191              continue
                  bvec(nvl) = -bvec(nvl)
               endif
            endif
            goto 700
         endif
      endif
      goto 50
!c
!c Drop constraint it1
!c
 700  continue
!c
!c if it1 = nact it is only necessary to update the vector u and nact
!c
      if (it1 .EQ. nact) goto 799
!c
!c After updating one row of R (column of J) we will also come back here
!c
 797  continue
!c
!c we have to find the Givens rotation which will reduce the element
!c (it1+1,it1+1) of R to zero.
!c if it is already zero we don't have to do anything except of updating
!c u, iact, and shifting column (it1+1) of R to column (it1)
!c l  will point to element (1,it1+1) of R
!c l1 will point to element (it1+1,it1+1) of R
!c
      l  = iwrm + (it1*(it1+1))/2 + 1
      l1 = l+it1
      if (work(l1) .EQ. 0.d0) goto 798
      gc   = max(abs(work(l1-1)),abs(work(l1)))
      gs   = min(abs(work(l1-1)),abs(work(l1)))
      temp = sign(gc*sqrt(1+gs*gs/(gc*gc)), work(l1-1))
      gc   = work(l1-1)/temp
      gs   = work(l1)/temp
!c
!c The Givens rotatin is done with the matrix (gc gs, gs -gc).
!c If gc is one, then element (it1+1,it1+1) of R is zero compared with
!c element (it1,it1+1). Hence we don't have to do anything.
!c if gc is zero, then we just have to switch row (it1) and row (it1+1)
!c of R and column (it1) and column (it1+1) of J. Since we swithc rows in
!c R and columns in J, we can ignore the sign of gs.
!c Otherwise we have to apply the Givens rotation to these rows/columns.
!c
      if (gc .EQ. 1.d0) goto 798
      if (gc .EQ. 0.d0) then
         do 710 i=it1+1,nact
            temp       = work(l1-1)
            work(l1-1) = work(l1)
            work(l1)   = temp
            l1 = l1+i
 710     continue
         do 711 i=1,n
            temp          = dmat(i,it1)
            dmat(i,it1)   = dmat(i,it1+1)
            dmat(i,it1+1) = temp
 711     continue
      else
         nu = gs/(1.d0+gc)
         do 720 i=it1+1,nact
            temp       = gc*work(l1-1) + gs*work(l1)
            work(l1)   = nu*(work(l1-1)+temp) - work(l1)
            work(l1-1) = temp
            l1 = l1+i
 720     continue
         do 721 i=1,n
            temp          = gc*dmat(i,it1) + gs*dmat(i,it1+1)
            dmat(i,it1+1) = nu*(dmat(i,it1)+temp) - dmat(i,it1+1)
            dmat(i,it1)   = temp
 721     continue
      endif
!c
!c shift column (it1+1) of R to column (it1) (that is, the first it1
!c elements). The posit1on of element (1,it1+1) of R was calculated above
!c and stored in l.
!c
 798  continue
      l1 = l-it1
      do 730 i=1,it1
         work(l1)=work(l)
         l  = l+1
         l1 = l1+1
 730  continue
!c
!c update vector u and iact as necessary
!c Continue with updating the matrices J and R
!c
      work(iwuv+it1) = work(iwuv+it1+1)
      iact(it1)      = iact(it1+1)
      it1 = it1+1
      if (it1 .LT. nact) goto 797
 799  work(iwuv+nact)   = work(iwuv+nact+1)
      work(iwuv+nact+1) = 0.d0
      iact(nact)        = 0
      nact = nact-1
      iter(2) = iter(2)+1
      goto 55
 999  continue
      return
end

subroutine dpori(a,lda,n)
      integer :: lda,n
      real(flyt) :: a(lda,*)
!c
!c     dpori computes the inverse of the factor of a
!c     double precision symmetric positive definite matrix
!c     using the factors computed by dpofa.
!c
!c     modification of dpodi by BaT 05/11/95
!c
!c     on entry
!c
!c        a       double precision(lda, n)
!c                the output  a  from dpofa
!c
!c        lda     integer
!c                the leading dimension of the array  a .
!c
!c        n       integer
!c                the order of the matrix  a .
!c
!c     on return
!c
!c        a       if dpofa was used to factor  a  then
!c                dpodi produces the upper half of inverse(a) .
!c                elements of  a  below the diagonal are unchanged.
!c
!c     error condition
!c
!c        a division by zero will occur if the input factor contains
!c        a zero on the diagonal and the inverse is requested.
!c        it will not occur if the subroutines are called correctly
!c        and if dpoco or dpofa has set info .eq. 0 .
!c
!c     linpack.  this version dated 08/14/78 .
!c     cleve moler, university of new mexico, argonne national lab.
!c     modified by Berwin A. Turlach 05/11/95
!c
!c     subroutines and functions
!c
!c     blas daxpy,dscal
!c     fortran mod
!c
!c     internal variables
!c
      double precision :: t
      integer :: j,k,kp1
!c
!c     compute inverse(r)
!c
      do 100 k = 1, n
         a(k,k) = 1.0d0/a(k,k)
         t = -a(k,k)
         call dscal(k-1,t,a(1,k),1)
         kp1 = k + 1
         if (n .lt. kp1) go to 90
         do 80 j = kp1, n
            t = a(k,j)
            a(k,j) = 0.0d0
            call daxpy(k,t,a(1,k),1,a(1,j),1)
 80      continue
 90      continue
 100  continue
      return
end

subroutine dposl(a,lda,n,b)
      integer :: lda,n
      real(flyt) :: a(lda,*),b(*)
!c
!c     dposl solves the double precision symmetric positive definite
!c     system a * x = b
!c     using the factors computed by dpoco or dpofa.
!c
!c     on entry
!c
!c        a       double precision(lda, n)
!c                the output from dpoco or dpofa.
!c
!c        lda     integer
!c                the leading dimension of the array  a .
!c
!c        n       integer
!c                the order of the matrix  a .
!c
!c        b       double precision(n)
!c                the right hand side vector.
!c
!c     on return
!c
!c        b       the solution vector  x .
!c
!c     error condition
!c
!c        a division by zero will occur if the input factor contains
!c        a zero on the diagonal.  technically this indicates
!c        singularity but it is usually caused by improper subroutine
!c        arguments.  it will not occur if the subroutines are called
!c        correctly and  info .eq. 0 .
!c
!c     to compute  inverse(a) * c  where  c  is a matrix
!c     with  p  columns
!c           call dpoco(a,lda,n,rcond,z,info)
!c           if (rcond is too small .or. info .ne. 0) go to ...
!c           do 10 j = 1, p
!c              call dposl(a,lda,n,c(1,j))
!c        10 continue
!c
!c     linpack.  this version dated 08/14/78 .
!c     cleve moler, university of new mexico, argonne national lab.
!c
!c     subroutines and functions
!c
!c     blas daxpy,ddot
!c
!c     internal variables
!c
      real(flyt) :: ddot,t
      integer k,kb
!c
!c     solve trans(r)*y = b
!c
      do 10 k = 1, n
         t = ddot(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k) - t)/a(k,k)
   10 continue
!c
!c     solve r*x = y
!c
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call daxpy(k-1,t,a(1,k),1,b(1),1)
   20 continue
      return
      end

      subroutine dpofa(a,n,info)
      integer n,info
      real(flyt), dimension(:,:), intent(inout) :: a
      !real(flyt) :: a(lda,1)
!c
!c     dpofa factors a double precision symmetric positive definite
!c     matrix.
!c
!c     dpofa is usually called by dpoco, but it can be called
!c     directly with a saving in time if  rcond  is not needed.
!c     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
!c
!c     on entry
!c
!c        a       double precision(lda, n)
!c                the symmetric matrix to be factored.  only the
!c                diagonal and upper triangle are used.
!c
!c        lda     integer
!c                the leading dimension of the array  a .
!c
!c        n       integer
!c                the order of the matrix  a .
!c
!c     on return
!c
!c        a       an upper triangular matrix  r  so that  a = trans(r)*r
!c                where  trans(r)  is the transpose.
!c                the strict lower triangle is unaltered.
!c                if  info .ne. 0 , the factorization is not complete.
!c
!c        info    integer
!c                = 0  for normal return.
!c                = k  signals an error condition.  the leading minor
!c                     of order  k  is not positive definite.
!c
!c     linpack.  this version dated 08/14/78 .
!c     cleve moler, university of new mexico, argonne national lab.
!c
!c     subroutines and functions
!c
!c     blas ddot
!c     fortran dsqrt
!c
!c     internal variables
!c
      real(flyt) :: t,s,ddot
      integer :: j,jm1,k
!c     begin block with ...exits to 40
!c
!c
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
!c     ......exit
            if (s .le. 0.0d0) go to 40
            a(j,j) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
end

end submodule
! Just keep a copy of the original, just in case
!!!   !c
!!!   !c  Copyright (C) 1995-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com>
!!!   !c
!!!   !c  This program is free software; you can redistribute it and/or modify
!!!   !c  it under the terms of the GNU General Public License as published by
!!!   !c  the Free Software Foundation; either version 2 of the License, or
!!!   !c  (at your option) any later version.
!!!   !c
!!!   !c  This program is distributed in the hope that it will be useful,
!!!   !c  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!   !c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!   !c  GNU General Public License for more details.
!!!   !c
!!!   !c  You should have received a copy of the GNU General Public License
!!!   !c  along with this program; if not, write to the Free Software
!!!   !c  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
!!!   !c  USA.
!!!   !c
!!!   !c  this routine uses the Goldfarb/Idnani algorithm to solve the
!!!   !c  following minimization problem:
!!!   !c
!!!   !c        minimize  -d^T x + 1/2 *  x^T D x
!!!   !c        where   A1^T x  = b1
!!!   !c                A2^T x >= b2
!!!   !c
!!!   !c  the matrix D is assumed to be positive definite.  Especially,
!!!   !c  w.l.o.g. D is assumed to be symmetric.
!!!   !c
!!!   !c  Input parameter:
!!!   !c  dmat   nxn matrix, the matrix D from above (dp)
!!!   !c         *** WILL BE DESTROYED ON EXIT ***
!!!   !c         The user has two possibilities:
!!!   !c         a) Give D (ierr=0), in this case we use routines from LINPACK
!!!   !c            to decompose D.
!!!   !c         b) To get the algorithm started we need R^-1, where D=R^TR.
!!!   !c            So if it is cheaper to calculate R^-1 in another way (D may
!!!   !c            be a band matrix) then with the general routine, the user
!!!   !c            may pass R^{-1}.  Indicated by ierr not equal to zero.
!!!   !c  dvec   nx1 vector, the vector d from above (dp)
!!!   !c         *** WILL BE DESTROYED ON EXIT ***
!!!   !c         contains on exit the solution to the initial, i.e.,
!!!   !c         unconstrained problem
!!!   !c  fddmat scalar, the leading dimension of the matrix dmat
!!!   !c  n      the dimension of dmat and dvec (int)
!!!   !c  amat   lxq matrix (dp)
!!!   !c         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!!!   !c             CHANGED SIGNES ON EXIT ***
!!!   !c  iamat  (l+1)xq matrix (int)
!!!   !c         these two matrices store the matrix A in compact form. the format
!!!   !c         is: [ A=(A1 A2)^T ]
!!!   !c           iamat(1,i) is the number of non-zero elements in column i of A
!!!   !c           iamat(k,i) for k>=2, is equal to j if the (k-1)-th non-zero
!!!   !c                      element in column i of A is A(i,j)
!!!   !c            amat(k,i) for k>=1, is equal to the k-th non-zero element
!!!   !c                      in column i of A.
!!!   !c
!!!   !c  bvec   qx1 vector, the vector of constants b in the constraints (dp)
!!!   !c         [ b = (b1^T b2^T)^T ]
!!!   !c         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!!!   !c             CHANGED SIGNES ON EXIT ***
!!!   !c  fdamat the first dimension of amat as declared in the calling program.
!!!   !c         fdamat >= n (and iamat must have fdamat+1 as first dimension)
!!!   !c  q      integer, the number of constraints.
!!!   !c  meq    integer, the number of equality constraints, 0 <= meq <= q.
!!!   !c  ierr   integer, code for the status of the matrix D:
!!!   !c            ierr =  0, we have to decompose D
!!!   !c            ierr != 0, D is already decomposed into D=R^TR and we were
!!!   !c                       given R^{-1}.
!!!   !c
!!!   !c  Output parameter:
!!!   !c  sol   nx1 the final solution (x in the notation above)
!!!   !c  lagr  qx1 the final Lagrange multipliers
!!!   !c  crval scalar, the value of the criterion at the minimum
!!!   !c  iact  qx1 vector, the constraints which are active in the final
!!!   !c        fit (int)
!!!   !c  nact  scalar, the number of constraints active in the final fit (int)
!!!   !c  iter  2x1 vector, first component gives the number of "main"
!!!   !c        iterations, the second one says how many constraints were
!!!   !c        deleted after they became active
!!!   !c  ierr  integer, error code on exit, if
!!!   !c           ierr = 0, no problems
!!!   !c           ierr = 1, the minimization problem has no solution
!!!   !c           ierr = 2, problems with decomposing D, in this case sol
!!!   !c                     contains garbage!!
!!!   !c
!!!   !c  Working space:
!!!   !c  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1
!!!   !c        where r=min(n,q)
!!!   !c
!!!   ! OH: I just rewrote the routine to modern fortran
!!!   subroutine qpgen1(dmat, dvec, fddmat, n, sol, lagr, crval, amat, iamat,bvec, fdamat, q, meq, iact, nact, iter, work, ierr)
!!!         integer :: n, i, j, l, l1
!!!         integer :: info, q, fdamat, iamat(fdamat+1,*), iact(*), iter(*), it1
!!!         integer :: ierr, nact, iwzv, iwrv, iwrm, iwsv, iwuv, nvl
!!!         integer :: r, iwnbv, meq, fddmat
!!!         real(flyt) :: dmat(fddmat,*), dvec(*),sol(*), lagr(*), bvec(*)
!!!         real(flyt) :: work(*), temp, sum, t1, tt, gc, gs, crval
!!!         real(flyt) :: nu, amat(fdamat,*), vsmall, tmpa, tmpb
!!!         logical t1inf, t2min
!!!
!!!         r = min(n,q)
!!!         l = 2*n + (r*(r+5))/2 + 2*q + 1
!!!   !c
!!!   !c     code gleaned from Powell's ZQPCVX routine to determine a small
!!!   !c     number  that can be assumed to be an upper bound on the relative
!!!   !c     precision of the computer arithmetic.
!!!   !c
!!!         vsmall = 1.0d-60
!!!    1    vsmall = vsmall + vsmall
!!!         tmpa = 1.0d0 + 0.1d0*vsmall
!!!         tmpb = 1.0d0 + 0.2d0*vsmall
!!!         if( tmpa .LE. 1.0d0 ) goto 1
!!!         if( tmpb .LE. 1.0d0 ) goto 1
!!!   !c
!!!   !c store the initial dvec to calculate below the unconstrained minima of
!!!   !c the critical value.
!!!   !c
!!!         do 10 i=1,n
!!!            work(i) = dvec(i)
!!!    10   continue
!!!         do 11 i=n+1,l
!!!            work(i) = 0.d0
!!!    11   continue
!!!         do 12 i=1,q
!!!            iact(i) = 0
!!!            lagr(i) = 0.d0
!!!    12   continue
!!!   !c
!!!   !c get the initial solution
!!!   !c
!!!         if( ierr .EQ. 0 )then
!!!            call dpofa(dmat,fddmat,n,info)
!!!            if( info .NE. 0 )then
!!!               ierr = 2
!!!               goto 999
!!!            endif
!!!            call dposl(dmat,fddmat,n,dvec)
!!!            call dpori(dmat,fddmat,n)
!!!         else
!!!   !c
!!!   !c Matrix D is already factorized, so we have to multiply d first with
!!!   !c R^-T and then with R^-1.  R^-1 is stored in the upper half of the
!!!   !c array dmat.
!!!   !c
!!!            do 20 j=1,n
!!!               sol(j)  = 0.d0
!!!               do 21 i=1,j
!!!                  sol(j) = sol(j) + dmat(i,j)*dvec(i)
!!!    21         continue
!!!    20      continue
!!!            do 22 j=1,n
!!!               dvec(j) = 0.d0
!!!               do 23 i=j,n
!!!                  dvec(j) = dvec(j) + dmat(j,i)*sol(i)
!!!    23         continue
!!!    22      continue
!!!         endif
!!!   !c
!!!   !c set lower triangular of dmat to zero, store dvec in sol and
!!!   !c calculate value of the criterion at unconstrained minima
!!!   !c
!!!         crval = 0.d0
!!!         do 30 j=1,n
!!!            sol(j)  = dvec(j)
!!!            crval   = crval + work(j)*sol(j)
!!!            work(j) = 0.d0
!!!            do 32 i=j+1,n
!!!               dmat(i,j) = 0.d0
!!!    32      continue
!!!    30   continue
!!!         crval = -crval/2.d0
!!!         ierr  = 0
!!!   !c
!!!   !c calculate some constants, i.e., from which index on the different
!!!   !c quantities are stored in the work matrix
!!!   !c
!!!         iwzv  = n
!!!         iwrv  = iwzv + n
!!!         iwuv  = iwrv + r
!!!         iwrm  = iwuv + r+1
!!!         iwsv  = iwrm + (r*(r+1))/2
!!!         iwnbv = iwsv + q
!!!   !c
!!!   !c calculate the norm of each column of the A matrix
!!!   !c
!!!         do 51 i=1,q
!!!            sum = 0.d0
!!!            do 52 j=1,iamat(1,i)
!!!               sum = sum + amat(j,i)*amat(j,i)
!!!    52      continue
!!!            work(iwnbv+i) = sqrt(sum)
!!!    51   continue
!!!         nact = 0
!!!         iter(1) = 0
!!!         iter(2) = 0
!!!    50   continue
!!!   !c
!!!   !c start a new iteration
!!!   !c
!!!         iter(1) = iter(1)+1
!!!   !c
!!!   !c calculate all constraints and check which are still violated
!!!   !c for the equality constraints we have to check whether the normal
!!!   !c vector has to be negated (as well as bvec in that case)
!!!   !c
!!!         l = iwsv
!!!         do 60 i=1,q
!!!            l = l+1
!!!            sum = -bvec(i)
!!!            do 61 j = 1,iamat(1,i)
!!!               sum = sum + amat(j,i)*sol(iamat(j+1,i))
!!!    61      continue
!!!            if ( abs(sum) .LT. vsmall ) then
!!!               sum = 0.0d0
!!!            endif
!!!            if (i .GT. meq) then
!!!               work(l) = sum
!!!            else
!!!               work(l) = -abs(sum)
!!!               if (sum .GT. 0.d0) then
!!!                  do 62 j=1,iamat(1,i)
!!!                     amat(j,i) = -amat(j,i)
!!!    62            continue
!!!                  bvec(i) = -bvec(i)
!!!               endif
!!!            endif
!!!    60   continue
!!!   !c
!!!   !c as safeguard against rounding errors set already active constraints
!!!   !c explicitly to zero
!!!   !c
!!!         do 70 i=1,nact
!!!            work(iwsv+iact(i)) = 0.d0
!!!    70   continue
!!!   !c
!!!   !c we weight each violation by the number of non-zero elements in the
!!!   !c corresponding row of A. then we choose the violated constraint which
!!!   !c has maximal absolute value, i.e., the minimum.
!!!   !c by obvious commenting and uncommenting we can choose the strategy to
!!!   !c take always the first constraint which is violated. ;-)
!!!   !c
!!!         nvl = 0
!!!         temp = 0.d0
!!!         do 71 i=1,q
!!!            if (work(iwsv+i) .LT. temp*work(iwnbv+i)) then
!!!               nvl = i
!!!               temp = work(iwsv+i)/work(iwnbv+i)
!!!            endif
!!!   !c         if (work(iwsv+i) .LT. 0.d0) then
!!!   !c            nvl = i
!!!   !c            goto 72
!!!   !c         endif
!!!    71   continue
!!!    72   if (nvl .EQ. 0) then
!!!            do 73 i=1,nact
!!!               lagr(iact(i))=work(iwuv+i)
!!!    73      continue
!!!            goto 999
!!!         endif
!!!   !c
!!!   !c calculate d=J^Tn^+ where n^+ is the normal vector of the violated
!!!   !c constraint. J is stored in dmat in this implementation!!
!!!   !c if we drop a constraint, we have to jump back here.
!!!   !c
!!!    55   continue
!!!         do 80 i=1,n
!!!            sum = 0.d0
!!!            do 81 j=1,iamat(1,nvl)
!!!               sum = sum + dmat(iamat(j+1,nvl),i)*amat(j,nvl)
!!!    81      continue
!!!            work(i) = sum
!!!    80   continue
!!!   !c
!!!   !c Now calculate z = J_2 d_2
!!!   !c
!!!         l1 = iwzv
!!!         do 90 i=1,n
!!!            work(l1+i) =0.d0
!!!    90   continue
!!!         do 92 j=nact+1,n
!!!            do 93 i=1,n
!!!               work(l1+i) = work(l1+i) + dmat(i,j)*work(j)
!!!    93      continue
!!!    92   continue
!!!   !c
!!!   !c and r = R^{-1} d_1, check also if r has positive elements (among the
!!!   !c entries corresponding to inequalities constraints).
!!!   !c
!!!         t1inf = .TRUE.
!!!         do 95 i=nact,1,-1
!!!            sum = work(i)
!!!            l  = iwrm+(i*(i+3))/2
!!!            l1 = l-i
!!!            do 96 j=i+1,nact
!!!               sum = sum - work(l)*work(iwrv+j)
!!!               l   = l+j
!!!    96      continue
!!!            sum = sum / work(l1)
!!!            work(iwrv+i) = sum
!!!            if (iact(i) .LE. meq) goto 95
!!!            if (sum .LE. 0.d0) goto 95
!!!    7       t1inf = .FALSE.
!!!            it1 = i
!!!    95   continue
!!!   !c
!!!   !c if r has positive elements, find the partial step length t1, which is
!!!   !c the maximum step in dual space without violating dual feasibility.
!!!   !c it1  stores in which component t1, the min of u/r, occurs.
!!!   !c
!!!         if ( .NOT. t1inf) then
!!!            t1   = work(iwuv+it1)/work(iwrv+it1)
!!!            do 100 i=1,nact
!!!               if (iact(i) .LE. meq) goto 100
!!!               if (work(iwrv+i) .LE. 0.d0) goto 100
!!!               temp = work(iwuv+i)/work(iwrv+i)
!!!               if (temp .LT. t1) then
!!!                  t1   = temp
!!!                  it1  = i
!!!               endif
!!!    100     continue
!!!         endif
!!!   !c
!!!   !c test if the z vector is equal to zero
!!!   !c
!!!         sum = 0.d0
!!!         do 110 i=iwzv+1,iwzv+n
!!!            sum = sum + work(i)*work(i)
!!!    110  continue
!!!         if (abs(sum) .LE. vsmall) then
!!!   !c
!!!   !c No step in primal space such that the new constraint becomes
!!!   !c feasible. Take step in dual space and drop a constant.
!!!   !c
!!!            if (t1inf) then
!!!   !c
!!!   !c No step in dual space possible either, problem is not solvable
!!!   !c
!!!               ierr = 1
!!!               goto 999
!!!            else
!!!   !c
!!!   !c we take a partial step in dual space and drop constraint it1,
!!!   !c that is, we drop the it1-th active constraint.
!!!   !c then we continue at step 2(a) (marked by label 55)
!!!   !c
!!!               do 111 i=1,nact
!!!                  work(iwuv+i) = work(iwuv+i) - t1*work(iwrv+i)
!!!    111        continue
!!!               work(iwuv+nact+1) = work(iwuv+nact+1) + t1
!!!               goto 700
!!!            endif
!!!         else
!!!   !c
!!!   !c compute full step length t2, minimum step in primal space such that
!!!   !c the constraint becomes feasible.
!!!   !c keep sum (which is z^Tn^+) to update crval below!
!!!   !c
!!!            sum = 0.d0
!!!            do 120 i = 1,iamat(1,nvl)
!!!               sum = sum + work(iwzv+iamat(i+1,nvl))*amat(i,nvl)
!!!    120     continue
!!!            tt = -work(iwsv+nvl)/sum
!!!            t2min = .TRUE.
!!!            if (.NOT. t1inf) then
!!!               if (t1 .LT. tt) then
!!!                  tt    = t1
!!!                  t2min = .FALSE.
!!!               endif
!!!            endif
!!!   !c
!!!   !c take step in primal and dual space
!!!   !c
!!!            do 130 i=1,n
!!!               sol(i) = sol(i) + tt*work(iwzv+i)
!!!    130     continue
!!!            crval = crval + tt*sum*(tt/2.d0 + work(iwuv+nact+1))
!!!            do 131 i=1,nact
!!!               work(iwuv+i) = work(iwuv+i) - tt*work(iwrv+i)
!!!    131     continue
!!!            work(iwuv+nact+1) = work(iwuv+nact+1) + tt
!!!   !c
!!!   !c if it was a full step, then we check wheter further constraints are
!!!   !c violated otherwise we can drop the current constraint and iterate once
!!!   !c more
!!!            if(t2min) then
!!!   !c
!!!   !c we took a full step. Thus add constraint nvl to the list of active
!!!   !c constraints and update J and R
!!!   !c
!!!               nact = nact + 1
!!!               iact(nact) = nvl
!!!   !c
!!!   !c to update R we have to put the first nact-1 components of the d vector
!!!   !c into column (nact) of R
!!!   !c
!!!               l = iwrm + ((nact-1)*nact)/2 + 1
!!!               do 150 i=1,nact-1
!!!                  work(l) = work(i)
!!!                  l = l+1
!!!    150        continue
!!!   !c
!!!   !c if now nact=n, then we just have to add the last element to the new
!!!   !c row of R.
!!!   !c Otherwise we use Givens transformations to turn the vector d(nact:n)
!!!   !c into a multiple of the first unit vector. That multiple goes into the
!!!   !c last element of the new row of R and J is accordingly updated by the
!!!   !c Givens transformations.
!!!   !c
!!!               if (nact .EQ. n) then
!!!                  work(l) = work(n)
!!!               else
!!!                  do 160 i=n,nact+1,-1
!!!   !c
!!!   !c we have to find the Givens rotation which will reduce the element
!!!   !c (l1) of d to zero.
!!!   !c if it is already zero we don't have to do anything, except of
!!!   !c decreasing l1
!!!   !c
!!!                     if (work(i) .EQ. 0.d0) goto 160
!!!                     gc   = max(abs(work(i-1)),abs(work(i)))
!!!                     gs   = min(abs(work(i-1)),abs(work(i)))
!!!                     temp = sign(gc*sqrt(1+gs*gs/(gc*gc)), work(i-1))
!!!                     gc   = work(i-1)/temp
!!!                     gs   = work(i)/temp
!!!   !c
!!!   !c The Givens rotation is done with the matrix (gc gs, gs -gc).
!!!   !c If gc is one, then element (i) of d is zero compared with element
!!!   !c (l1-1). Hence we don't have to do anything.
!!!   !c If gc is zero, then we just have to switch column (i) and column (i-1)
!!!   !c of J. Since we only switch columns in J, we have to be careful how we
!!!   !c update d depending on the sign of gs.
!!!   !c Otherwise we have to apply the Givens rotation to these columns.
!!!   !c The i-1 element of d has to be updated to temp.
!!!   !c
!!!                     if (gc .EQ. 1.d0) goto 160
!!!                     if (gc .EQ. 0.d0) then
!!!                        work(i-1) = gs * temp
!!!                        do 170 j=1,n
!!!                           temp        = dmat(j,i-1)
!!!                           dmat(j,i-1) = dmat(j,i)
!!!                           dmat(j,i)   = temp
!!!    170                 continue
!!!                     else
!!!                        work(i-1) = temp
!!!                        nu = gs/(1.d0+gc)
!!!                        do 180 j=1,n
!!!                           temp        = gc*dmat(j,i-1) + gs*dmat(j,i)
!!!                           dmat(j,i)   = nu*(dmat(j,i-1)+temp) - dmat(j,i)
!!!                           dmat(j,i-1) = temp
!!!    180                 continue
!!!                     endif
!!!    160           continue
!!!   !c
!!!   !c l is still pointing to element (nact,nact) of the matrix R.
!!!   !c So store d(nact) in R(nact,nact)
!!!                  work(l) = work(nact)
!!!               endif
!!!            else
!!!   !c
!!!   !c we took a partial step in dual space. Thus drop constraint it1,
!!!   !c that is, we drop the it1-th active constraint.
!!!   !c then we continue at step 2(a) (marked by label 55)
!!!   !c but since the fit changed, we have to recalculate now "how much"
!!!   !c the fit violates the chosen constraint now.
!!!   !c
!!!               sum = -bvec(nvl)
!!!               do 190 j = 1,iamat(1,nvl)
!!!                  sum = sum + sol(iamat(j+1,nvl))*amat(j,nvl)
!!!    190        continue
!!!               if( nvl .GT. meq ) then
!!!                  work(iwsv+nvl) = sum
!!!               else
!!!                  work(iwsv+nvl) = -abs(sum)
!!!                  if( sum .GT. 0.d0) then
!!!                     do 191 j=1,iamat(1,nvl)
!!!                        amat(j,nvl) = -amat(j,nvl)
!!!    191              continue
!!!                     bvec(nvl) = -bvec(nvl)
!!!                  endif
!!!               endif
!!!               goto 700
!!!            endif
!!!         endif
!!!         goto 50
!!!   !c
!!!   !c Drop constraint it1
!!!   !c
!!!    700  continue
!!!   !c
!!!   !c if it1 = nact it is only necessary to update the vector u and nact
!!!   !c
!!!         if (it1 .EQ. nact) goto 799
!!!   !c
!!!   !c After updating one row of R (column of J) we will also come back here
!!!   !c
!!!    797  continue
!!!   !c
!!!   !c we have to find the Givens rotation which will reduce the element
!!!   !c (it1+1,it1+1) of R to zero.
!!!   !c if it is already zero we don't have to do anything except of updating
!!!   !c u, iact, and shifting column (it1+1) of R to column (it1)
!!!   !c l  will point to element (1,it1+1) of R
!!!   !c l1 will point to element (it1+1,it1+1) of R
!!!   !c
!!!         l  = iwrm + (it1*(it1+1))/2 + 1
!!!         l1 = l+it1
!!!         if (work(l1) .EQ. 0.d0) goto 798
!!!         gc   = max(abs(work(l1-1)),abs(work(l1)))
!!!         gs   = min(abs(work(l1-1)),abs(work(l1)))
!!!         temp = sign(gc*sqrt(1+gs*gs/(gc*gc)), work(l1-1))
!!!         gc   = work(l1-1)/temp
!!!         gs   = work(l1)/temp
!!!   !c
!!!   !c The Givens rotatin is done with the matrix (gc gs, gs -gc).
!!!   !c If gc is one, then element (it1+1,it1+1) of R is zero compared with
!!!   !c element (it1,it1+1). Hence we don't have to do anything.
!!!   !c if gc is zero, then we just have to switch row (it1) and row (it1+1)
!!!   !c of R and column (it1) and column (it1+1) of J. Since we swithc rows in
!!!   !c R and columns in J, we can ignore the sign of gs.
!!!   !c Otherwise we have to apply the Givens rotation to these rows/columns.
!!!   !c
!!!         if (gc .EQ. 1.d0) goto 798
!!!         if (gc .EQ. 0.d0) then
!!!            do 710 i=it1+1,nact
!!!               temp       = work(l1-1)
!!!               work(l1-1) = work(l1)
!!!               work(l1)   = temp
!!!               l1 = l1+i
!!!    710     continue
!!!            do 711 i=1,n
!!!               temp          = dmat(i,it1)
!!!               dmat(i,it1)   = dmat(i,it1+1)
!!!               dmat(i,it1+1) = temp
!!!    711     continue
!!!         else
!!!            nu = gs/(1.d0+gc)
!!!            do 720 i=it1+1,nact
!!!               temp       = gc*work(l1-1) + gs*work(l1)
!!!               work(l1)   = nu*(work(l1-1)+temp) - work(l1)
!!!               work(l1-1) = temp
!!!               l1 = l1+i
!!!    720     continue
!!!            do 721 i=1,n
!!!               temp          = gc*dmat(i,it1) + gs*dmat(i,it1+1)
!!!               dmat(i,it1+1) = nu*(dmat(i,it1)+temp) - dmat(i,it1+1)
!!!               dmat(i,it1)   = temp
!!!    721     continue
!!!         endif
!!!   !c
!!!   !c shift column (it1+1) of R to column (it1) (that is, the first it1
!!!   !c elements). The posit1on of element (1,it1+1) of R was calculated above
!!!   !c and stored in l.
!!!   !c
!!!    798  continue
!!!         l1 = l-it1
!!!         do 730 i=1,it1
!!!            work(l1)=work(l)
!!!            l  = l+1
!!!            l1 = l1+1
!!!    730  continue
!!!   !c
!!!   !c update vector u and iact as necessary
!!!   !c Continue with updating the matrices J and R
!!!   !c
!!!         work(iwuv+it1) = work(iwuv+it1+1)
!!!         iact(it1)      = iact(it1+1)
!!!         it1 = it1+1
!!!         if (it1 .LT. nact) goto 797
!!!    799  work(iwuv+nact)   = work(iwuv+nact+1)
!!!         work(iwuv+nact+1) = 0.d0
!!!         iact(nact)        = 0
!!!         nact = nact-1
!!!         iter(2) = iter(2)+1
!!!         goto 55
!!!    999  continue
!!!         return
!!!   end
