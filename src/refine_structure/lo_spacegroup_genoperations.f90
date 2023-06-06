submodule(lo_spacegroup) lo_spacegroup_genoperations
!! Generates possible symmetry operations
implicit none
contains

!> generate pool of translation in a not-so-clever way
module subroutine generate_bruteforce_pool_of_translations(translations, nmax)
    !> fractional translations
    real(r8), dimension(:, :), allocatable, intent(out) :: translations
    !> largest divisor
    integer, intent(in) :: nmax

    real(r8), dimension(:, :), allocatable :: buf
    integer :: ctr, i, j, k, ii, jj, kk

    ! Count how many possible translations I can get?
    ctr = 0
    do ii = 1, nmax
    do jj = 1, nmax
    do kk = 1, nmax
        do i = 0, ii - 1
        do j = 0, jj - 1
        do k = 0, kk - 1
            ctr = ctr + 1
        end do
        end do
        end do
    end do
    end do
    end do
    allocate (buf(3, ctr))
    ctr = 0
    do ii = 1, nmax
    do jj = 1, nmax
    do kk = 1, nmax
        do i = 0, ii - 1
        do j = 0, jj - 1
        do k = 0, kk - 1
            ctr = ctr + 1
            buf(1, ctr) = real(i, r8)/real(ii, r8)
            buf(2, ctr) = real(j, r8)/real(jj, r8)
            buf(3, ctr) = real(k, r8)/real(kk, r8)
        end do
        end do
        end do
    end do
    end do
    end do
    ! Keep only the unique
    call lo_return_unique(buf, translations)
end subroutine

!> generate a pool of possible (fractional) translations
module subroutine generate_pool_of_translations(basis, position, label, pointoperations, translations, mem, verbosity)
    !> lattice vectors
    real(r8), dimension(3, 3), intent(in) :: basis
    !> positions of points (fractional coordinates)
    real(r8), dimension(:, :), intent(in) :: position
    !> label of points
    integer, dimension(:), intent(in) :: label
    !> pool of point operations
    real(r8), dimension(:, :, :), intent(in) :: pointoperations
    !> possible translations
    real(r8), dimension(:, :), allocatable, intent(out) :: translations
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), parameter :: tolerance = 1E-8_r8

    integer :: i, j, k, l, n, o
    real(r8), dimension(3, 3) :: inversebasis
    real(r8), dimension(3) :: v0, v1
    real(r8), dimension(:, :), allocatable :: r, buf
    real(r8) :: timer, t0, t1
    ! I operate under the assumption that a translation has to move one
    ! atom to another. Why do I know that? Well, I am allowed to rigidly translate the
    ! lattice within the cell without changing the symmetry. So I can put any atom at
    ! the origin. It will not move anywhere by the rotation part, so the translation
    ! part has to move it to another site with the same kind of atom. This can be repeated
    ! for all kinds of atoms. So the possible translations are then all the pair vectors in
    ! the cell, that move an atom of species A.

    timer = walltime()
    t0 = timer
    t1 = timer

    if (verbosity .gt. 0) then
        write (lo_iou, *) ''
        write (lo_iou, *) 'Generating pool of translations'
    end if

    ! get cartesian positions
    inversebasis = lo_invert3x3matrix(basis)
    n = size(position, 2)
    call mem%allocate(r, [3, n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    r = 0.0_r8
    do i = 1, n
        r(:, i) = matmul(basis, position(:, i))
    end do

    ! First rough count, I only want translations between the same species
    l = 0
    do i = 1, n
    do j = 1, n
        if (label(i) .eq. label(j)) l = l + 1
    end do
    end do
    l = l*size(pointoperations, 3)

    call mem%allocate(buf, [3, l], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf = 0.0_r8

    ! Start generating possible translations.
    l = 0
    do i = 1, n
    do j = i, n
        if (label(i) .ne. label(j)) cycle
        do o = 1, size(pointoperations, 3)
            ! possible vector
            v0 = r(:, j) - matmul(pointoperations(:, :, o), r(:, i))
            v0 = lo_clean_fractional_coordinates(matmul(inversebasis, v0))
            do k = 1, 3
                if (v0(k) .gt. tolerance) then
                    v1(k) = 1.0_r8/v0(k)
                else
                    v1(k) = 0.0_r8
                end if
            end do
            ! vector has to be a 1/n-thingy, where n is an integer
            if (sum(abs(v1 - anint(v1))) .lt. tolerance) then
                l = l + 1
                buf(:, l) = v0
            end if
        end do
    end do
    end do

    if (verbosity .gt. 0) then
        t1 = walltime()
        write (lo_iou, *) '... found ', tochar(l), ' initial trial translation (', tochar(t1 - t0), 's)'
        t0 = t1
    end if

    ! Make sure we have cleaned the positions properly
    buf(:, 1:l) = lo_chop(buf(:, 1:l), tolerance**2)
    buf(:, 1:l) = lo_clean_fractional_coordinates(buf(:, 1:l))
    buf(:, 1:l) = lo_chop(buf(:, 1:l), tolerance**2)
    buf(:, 1:l) = lo_clean_fractional_coordinates(buf(:, 1:l))
    ! And return the unique ones
    call lo_return_unique(buf, translations, tolerance)
    ! Make sure that the zero translation is first.
    do i = 1, size(translations, 2)
        if (norm2(translations(:, i)) .lt. tolerance) then
            v0 = translations(:, 1)
            translations(:, 1) = 0.0_r8
            translations(:, i) = v0
            exit
        end if
    end do

    ! And cleanup
    call mem%deallocate(r, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (verbosity .gt. 0) then
        t1 = walltime()
        write (lo_iou, *) '... got ', tochar(size(translations, 2)), ' trial translation (', tochar(t1 - timer), 's)'
    end if
end subroutine

!> generates a pool of possible (Cartesian) point operations.
module subroutine generate_pool_of_point_operations(basis, hexops, octops, hexerr, octerr, mem, verbosity)
    !> lattice vectors
    real(r8), dimension(3, 3), intent(in) :: basis
    !> list of point operations
    real(r8), dimension(3, 3, 24), intent(out) :: hexops
    real(r8), dimension(3, 3, 48), intent(out) :: octops
    !> error from point operations
    real(r8), intent(out) :: hexerr
    real(r8), intent(out) :: octerr
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), parameter :: tolerance = 1E-5_r8, tighttolerance = 1E-10_r8
    type(lo_voronoi_cell) :: voro
    real(r8), dimension(:, :), allocatable :: cluster, testvecs
    real(r8), dimension(3, 3) :: ident
    real(r8) :: timer, t0, t1

    init: block
        type(lo_distancetable) :: dt
        real(r8), dimension(:, :), allocatable :: r1
        real(r8), dimension(3, 1) :: r0
        real(r8), dimension(3) :: v0, v1
        real(r8) :: cutoff
        integer :: i, ctr

        timer = walltime()
        t0 = timer
        t1 = timer

        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'Generating pool of point operations'
        end if

        ! identity matrix
        ident = 0.0_r8
        ident(1, 1) = 1.0_r8
        ident(2, 2) = 1.0_r8
        ident(3, 3) = 1.0_r8

        ! Generate the Voronoi cell. Pick a very safe cutoff.
        cutoff = lo_bounding_sphere_of_box(basis)
        cutoff = max(cutoff, norm2(basis(:, 1)))
        cutoff = max(cutoff, norm2(basis(:, 2)))
        cutoff = max(cutoff, norm2(basis(:, 3)))
        cutoff = cutoff*5
        r0 = 0.0_r8
        call dt%generate(r0, basis, cutoff, -1)
        call voro%generate(dt%particle(1), cutoff*10, tolerance, mem)

        ! Keep a sensible cluster from the distance table to use for verification?
        cutoff = lo_bounding_sphere_of_box(basis)
        cutoff = max(cutoff, norm2(basis(:, 1)))
        cutoff = max(cutoff, norm2(basis(:, 2)))
        cutoff = max(cutoff, norm2(basis(:, 3)))
        cutoff = cutoff*1.2 + 10*tolerance
        ctr = 0
        do i = 1, dt%particle(1)%n
            if (dt%particle(1)%d(i) .lt. tolerance) cycle
            if (dt%particle(1)%d(i) .gt. cutoff) cycle
            ctr = ctr + 1
        end do
        allocate (cluster(3, ctr))
        ctr = 0
        do i = 1, dt%particle(1)%n
            if (dt%particle(1)%d(i) .lt. tolerance) cycle
            if (dt%particle(1)%d(i) .gt. cutoff) cycle
            ctr = ctr + 1
            cluster(:, ctr) = dt%particle(1)%v(:, i)
        end do

        ! Pick vectors from the Voronoi cell to use as
        ! candidate rotation axes.
        ctr = voro%n_face + voro%n_node
        allocate (r1(3, ctr))
        r0 = 0.0_r8
        ctr = 0
        do i = 1, voro%n_face
            ctr = ctr + 1
            r1(:, ctr) = voro%face(i)%neighbourvector
        end do
        do i = 1, voro%n_node
            ctr = ctr + 1
            r1(:, ctr) = voro%node(:, i)
        end do
        v1 = voro%node(:, 1) + lo_degenvector
        do i = 1, ctr
            v0 = r1(:, i)
            if (dot_product(v0, v1) .lt. 0.0_r8) v0 = -v0
            r1(:, i) = v0/norm2(v0)
        end do
        call lo_return_unique(r1, testvecs)

        ! Excellent. Have Voronoi cell.
        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... built Voronoi diagram (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block init

    ! Generate pool of hexagonal operations?
    hexvecs: block
        real(r8), dimension(3, 3, 24) :: op1
        real(r8), dimension(3, 3) :: m
        real(r8), dimension(3) :: v0, v1, v2, hv1, hv2
        real(r8) :: f0, err0, err1
        integer :: iface, jface, i

        err1 = lo_huge
        faceloop1: do iface = 1, voro%n_face
            ! Check with respect to other faces if there is an orthogonal vector?
            v0 = voro%face(iface)%neighbourvector
            v0 = v0/norm2(v0)
            faceloop2: do jface = iface + 1, voro%n_face
                v1 = voro%face(jface)%neighbourvector
                v1 = v1/norm2(v1)
                f0 = abs(dot_product(v0, v1))
                ! Skip collinear vectors
                if (abs(f0 - 1.0_r8) .lt. tolerance) cycle
                ! Get an orthogonal vector?
                if (f0 .gt. tighttolerance) then
                    v2 = lo_cross(v0, v1)
                    v2 = v2/norm2(v2)
                else
                    v2 = v1
                end if

                ! Ok, have two candidate vectors, get some operations
                call hexagonal_operations_from_two_vectors(v0, v2, tolerance, op1)
                ! Calculate error?
                err0 = cluster_transform_error(op1, cluster)
                ! Maybe we already found it?
                if (err0 .lt. tighttolerance) then
                    hv1 = v0
                    hv2 = v2
                    err1 = err0
                    exit faceloop1
                end if

                ! If not, keep the best one?
                if (err0 .lt. err1) then
                    hv1 = v0
                    hv2 = v2
                    err1 = err0
                end if
            end do faceloop2
        end do faceloop1
        hexerr = err1

        ! Put the identity first
        call hexagonal_operations_from_two_vectors(hv1, hv2, tolerance, hexops)
        do i = 1, 24
            if (norm2(hexops(:, :, i) - ident) .lt. 1E-10_r8) then
                m = hexops(:, :, 1)
                hexops(:, :, 1) = ident
                hexops(:, :, i) = m
                exit
            end if
        end do
        ! Put inversion as the second
        do i = 1, 24
            if (norm2(hexops(:, :, i) + ident) .lt. 1E-10_r8) then
                m = hexops(:, :, 1)
                hexops(:, :, 1) = -ident
                hexops(:, :, i) = m
                exit
            end if
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... found hexagonal operations (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block hexvecs

    ! Also pool of cubic operations?
    cubvecs: block
        real(r8), dimension(3, 3, 48) :: op0
        real(r8), dimension(3, 3) :: m
        real(r8), dimension(3) :: v0, v1, v2, cv1, cv2
        real(r8) :: f0, err0, err1
        integer :: ivec, jvec, i

        err1 = lo_huge
        faceloop1: do ivec = 1, size(testvecs, 2)
            ! Check with respect to other faces if there is an orthogonal vector?
            v0 = testvecs(:, ivec)
            faceloop2: do jvec = ivec + 1, size(testvecs, 2)
                v1 = testvecs(:, jvec)
                f0 = abs(dot_product(v0, v1))
                ! Skip collinear vectors
                if (abs(f0 - 1.0_r8) .lt. tolerance) cycle
                ! Get an orthogonal vector?
                if (f0 .gt. tighttolerance) then
                    v2 = lo_cross(v0, v1)
                    v2 = v2/norm2(v2)
                else
                    v2 = v1
                end if

                ! Ok, have two candidate vectors, get some operations
                call octahedral_operations_from_two_vectors(v0, v2, tolerance, op0)
                ! Calculate error?
                err0 = cluster_transform_error(op0, cluster)
                ! Maybe we already found it?
                if (err0 .lt. tighttolerance) then
                    cv1 = v0
                    cv2 = v2
                    err1 = err0
                    exit faceloop1
                end if

                ! If not, keep the best one?
                if (err0 .lt. err1) then
                    cv1 = v0
                    cv2 = v2
                    err1 = err0
                end if
            end do faceloop2
        end do faceloop1
        ! Store the error
        octerr = err1

        ! Put the identity first
        call octahedral_operations_from_two_vectors(cv1, cv2, tolerance, octops)
        do i = 1, 48
            if (norm2(octops(:, :, i) - ident) .lt. 1E-10_r8) then
                m = octops(:, :, 1)
                octops(:, :, 1) = ident
                octops(:, :, i) = m
                exit
            end if
        end do
        ! Put inversion as the second
        do i = 1, 48
            if (norm2(octops(:, :, i) + ident) .lt. 1E-10_r8) then
                m = octops(:, :, 1)
                octops(:, :, 1) = -ident
                octops(:, :, i) = m
                exit
            end if
        end do

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (lo_iou, *) '... found octahedral operations (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block cubvecs

    !@TODO add cleanup

    ! ! Build an initial pool of operations
    ! buildops: block
    !     real(r8), dimension(48) :: symerr
    !     real(r8), dimension(3,3) :: inversebasis,m,canbasis,caninvbasis
    !     real(r8) :: f0
    !     integer, dimension(48) :: ind
    !     integer, dimension(3) :: fracvals
    !     integer, dimension(3,3) :: op
    !     integer :: i,i1,i2,i3,i4,i5,i6,i7,i8,i9
    !
    !     ! Get the inverse basis
    !     inversebasis=lo_invert3x3matrix(basis)
    !
    !     ! Also the canonical basis and its inverse
    !     call delaunay_reduce_basis(basis,canbasis)
    !     caninvbasis=lo_invert3x3matrix(canbasis)
    !
    !     bflist=-lo_huge
    !     symerr=1E5_r8
    !
    !     ! Build all possible operations, all 3^9 of them, and
    !     ! keep the best ones.
    !     fracvals=[-1,0,1]
    !     il1: do i1=1,3
    !     do i2=1,3
    !     do i3=1,3
    !     do i4=1,3
    !     do i5=1,3
    !     do i6=1,3
    !     do i7=1,3
    !     do i8=1,3
    !     do i9=1,3
    !         ! possible operation in fractional coordinates
    !         op(:,1)=[fracvals(i1),fracvals(i2),fracvals(i3)]
    !         op(:,2)=[fracvals(i4),fracvals(i5),fracvals(i6)]
    !         op(:,3)=[fracvals(i7),fracvals(i8),fracvals(i9)]
    !         ! Convert it to Cartesian coordinates
    !         m=matmul(canbasis,op)
    !         m=matmul(m,caninvbasis)
    !         ! Calculate an error score?
    !         f0=norm2( matmul( m,transpose(m) )-ident )
    !         ! If this is smaller than the worst one so far, keep it.
    !         ! This is not the most efficient, I know, but I have a plan.
    !         if ( f0 .lt. symerr(48) ) then
    !             symerr(48)=f0
    !             bflist(:,:,48)=m
    !             call lo_qsort(symerr,ind)
    !             bflist=bflist(:,:,ind)
    !         endif
    !     enddo
    !     enddo
    !     enddo
    !     enddo
    !     enddo
    !     enddo
    !     enddo
    !     enddo
    !     enddo il1
    !
    !     ! Sanity check that there are no duplicates
    !     do i1=1,48
    !     do i2=i1+1,48
    !         f0=norm2(bflist(:,:,i1)-bflist(:,:,i2))
    !         if ( f0 .lt. 1E-10_r8 ) then
    !             call lo_stop_gracefully(['Found duplicate operation'],lo_exitcode_symmetry,__FILE__,__LINE__)
    !         endif
    !     enddo
    !     enddo
    !
    !     ! How many operations are valid?
    !     noperation=0
    !     do i=1,48
    !         if ( symerr(i) .lt. tolerance ) then
    !             noperation=noperation+1
    !         else
    !             exit
    !         endif
    !     enddo
    !
    !     if ( verbosity .gt. 0 ) then
    !         t1=walltime()
    !         write(lo_iou,*) '... found ',tochar(noperation),' initial trial operations (',tochar(t1-t0),'s)'
    !         t0=t1
    !     endif
    ! end block buildops
    !
    ! ! Try to align it a little to flesh out the number of operations?
    ! prepmany: block
    !     integer, dimension(48) :: angle,whatkind
    !     real(r8), dimension(:,:,:), allocatable :: buf_op
    !     real(r8), dimension(:,:), allocatable :: dr,testvecs
    !     real(r8), dimension(3,48) :: axes
    !     real(r8), dimension(3,3,48) :: op0
    !     real(r8), dimension(3,3,24) :: op1
    !     real(r8), dimension(3,3) :: m
    !     real(r8), dimension(3) :: cubvec1,cubvec2,hexvec1,hexvec2
    !     real(r8) :: f0
    !     integer :: i,j,k,l,ctr
    !
    !     ! Get some rotational axes to work with
    !     axes=0.0_r8
    !     angle=-1
    !     whatkind=-1
    !     do i=1,noperation
    !         call classify_point_operation(bflist(:,:,i),angle(i),axes(:,i),whatkind(i))
    !     enddo
    !
    !     ! Build some vectors to play with? Perhaps locate some 'cubic' vectors? I'm not
    !     ! completely sure what could be sensible candidates, but these are sensible
    !     ! guesstimates. Will likely add more in the future.
    !     ctr=6
    !     do i=1,noperation
    !         if ( norm2(axes(:,i)) .gt. tolerance ) ctr=ctr+1
    !         ctr=ctr+1
    !     enddo
    !     call mem%allocate(dr,[3,ctr],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     dr=0.0_r8
    !
    !     dr(:,1:3)=basis
    !     dr(:,4:6)=transpose( lo_invert3x3matrix(basis) )
    !
    !     ctr=6
    !     do i=1,noperation
    !         if ( norm2(axes(:,i)) .gt. tolerance ) then
    !             ctr=ctr+1
    !             dr(:,ctr)=axes(:,i)
    !         endif
    !     enddo
    !     do i=1,size(dr,2)
    !         f0=norm2(dr(:,i))
    !         if ( f0 .gt. tolerance ) then
    !             dr(:,i)=dr(:,i)/norm2(dr(:,i))
    !         else
    !             dr(:,i)=[1,0,0]
    !         endif
    !         ! Not sure if meaningful to have both positive and negative axes, get
    !         ! them pointing the same way, kind of.
    !         f0=dot_product(dr(:,i),[1.0_r8,0.1234568_r8,5.32131_r8])
    !         if ( abs(f0) .gt. tolerance ) then
    !             if ( f0 .lt. 0.0_r8 ) dr(:,i)=-dr(:,i)
    !         endif
    !     enddo
    !     call lo_return_unique(dr,testvecs)
    !     call mem%deallocate(dr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !
    !     ! Have the first set of testing vectors. Try to find a sensible match?
    !     l=0
    !     cubvec1=[1,0,0]
    !     cubvec2=[0,1,0]
    !     vl1: do i=1,size(testvecs,2)
    !     do j=i+1,size(testvecs,2)
    !         f0=dot_product(testvecs(:,i),testvecs(:,j))
    !         if ( abs(f0) .gt. tolerance ) cycle
    !
    !         call octahedral_operations_from_two_vectors(testvecs(:,i),testvecs(:,j),tolerance,op0)
    !         k=count_matching_matrices(bflist(:,:,1:noperation),op0,tolerance)
    !
    !         ! I might have a perfect match
    !         if ( k .eq. noperation ) then
    !             ! Found it! Good representation.
    !             cubvec1=testvecs(:,i)
    !             cubvec2=testvecs(:,j)
    !             l=48
    !             exit vl1
    !         endif
    !
    !         ! Maybe a decent match, keep the best one?
    !         if ( k .gt. l ) then
    !             l=k
    !             cubvec1=testvecs(:,i)
    !             cubvec2=testvecs(:,j)
    !         endif
    !     enddo
    !     enddo vl1
    !
    !     if ( verbosity .gt. 0 ) then
    !         t1=walltime()
    !         write(lo_iou,*) '... found ',tochar(l),' possible octahedral operations (',tochar(t1-t0),'s)'
    !         t0=t1
    !     endif
    !
    !     ! Do it again with the hexagonal stuff? Might as well.
    !     l=0
    !     hexvec1=[1,0,0]
    !     hexvec2=[0,1,0]
    !     vl2: do i=1,size(testvecs,2)
    !     do j=i+1,size(testvecs,2)
    !         f0=dot_product(testvecs(:,i),testvecs(:,j))
    !         if ( abs(f0) .gt. tolerance ) cycle
    !
    !         call hexagonal_operations_from_two_vectors(testvecs(:,i),testvecs(:,j),tolerance,op1)
    !         k=count_matching_matrices(bflist(:,:,1:noperation),op1,tolerance)
    !
    !         ! I might have a perfect match
    !         if ( k .eq. noperation ) then
    !             ! Found it! Good representation.
    !             l=24
    !             hexvec1=testvecs(:,i)
    !             hexvec2=testvecs(:,j)
    !             exit vl2
    !         elseif ( k .gt. l ) then
    !             l=k
    !             hexvec1=testvecs(:,i)
    !             hexvec2=testvecs(:,j)
    !         endif
    !
    !         ! Do it again but with the vectors swapped
    !         call hexagonal_operations_from_two_vectors(testvecs(:,j),testvecs(:,i),tolerance,op1)
    !         k=count_matching_matrices(bflist(:,:,1:noperation),op1,tolerance)
    !
    !         ! I might have a perfect match
    !         if ( k .eq. noperation ) then
    !             ! Found it! Good representation.
    !             l=24
    !             hexvec1=testvecs(:,j)
    !             hexvec2=testvecs(:,i)
    !             exit vl2
    !         elseif ( k .gt. l ) then
    !             l=k
    !             hexvec1=testvecs(:,j)
    !             hexvec2=testvecs(:,i)
    !         endif
    !     enddo
    !     enddo vl2
    !
    !     if ( verbosity .gt. 0 ) then
    !         t1=walltime()
    !         write(lo_iou,*) '... found ',tochar(l),' possible hexagonal operations (',tochar(t1-t0),'s)'
    !         t0=t1
    !     endif
    !
    !     ! Put everything together. Should give us a pretty complete list.
    !     ! However, this set of operations is not a proper group. Only the
    !     ! valid operations form a proper group.
    !     call mem%allocate(buf_op,[3,3,48+24+48],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !     buf_op=0.0_r8
    !     call octahedral_operations_from_two_vectors(cubvec1,cubvec2,tolerance,op0)
    !     call hexagonal_operations_from_two_vectors(hexvec1,hexvec2,tolerance,op1)
    !     l=0
    !     do i=1,48
    !         l=l+1
    !         buf_op(:,:,l)=op0(:,:,i)
    !     enddo
    !     do i=1,24
    !         l=l+1
    !         buf_op(:,:,l)=op1(:,:,i)
    !     enddo
    !     do i=1,48
    !         l=l+1
    !         buf_op(:,:,l)=bflist(:,:,i)
    !     enddo
    !     ! Decent start. Reduce to the unique?
    !     call lo_return_unique(buf_op,pointoperations,tolerance)
    !     call mem%deallocate(buf_op,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    !
    !     ! Put the identity first
    !     do i=1,size(pointoperations,3)
    !         if ( norm2(pointoperations(:,:,i)-ident) .lt. 1E-10_r8 ) then
    !             m=pointoperations(:,:,1)
    !             pointoperations(:,:,1)=ident
    !             pointoperations(:,:,i)=m
    !             exit
    !         endif
    !     enddo
    !     ! Put inversion as the second
    !     do i=1,size(pointoperations,3)
    !         if ( norm2(pointoperations(:,:,i)+ident) .lt. 1E-10_r8 ) then
    !             m=pointoperations(:,:,1)
    !             pointoperations(:,:,1)=-ident
    !             pointoperations(:,:,i)=m
    !             exit
    !         endif
    !     enddo
    !
    !     if ( verbosity .gt. 0 ) then
    !         t1=walltime()
    !         write(lo_iou,*) '... got ',tochar(size(pointoperations,3)),' possible operations (',tochar(t1-timer),'s)'
    !         t0=t1
    !     endif
    ! end block prepmany
end subroutine

!> generate the octahedral operations from two vectors.
subroutine octahedral_operations_from_two_vectors(A, B, tolerance, op)
    !> first vector
    real(r8), dimension(3), intent(in) :: A
    !> second vector
    real(r8), dimension(3), intent(in) :: B
    !> tolerance
    real(r8), intent(in) :: tolerance
    !> point operations
    real(r8), dimension(3, 3, 48) :: op

    real(r8), dimension(3, 3, 48) :: fop
    real(r8), dimension(3, 3, 5) :: generators
    real(r8), dimension(3, 3) :: m0, m1, m2
    real(r8), dimension(3) :: C, v0
    integer :: iter, ctr, i, j, k, ii

    ! First check that the vectors are orthogonal
    if (abs(dot_product(A, B)) .gt. tolerance) then
        call lo_stop_gracefully(['Generating vectors not orthogonal'], lo_exitcode_symmetry, __FILE__, __LINE__)
    end if
    if (abs(norm2(A) - 1.0_r8) .gt. tolerance) then
        call lo_stop_gracefully(['Generating vector not unit length'], lo_exitcode_symmetry, __FILE__, __LINE__)
    end if
    if (abs(norm2(B) - 1.0_r8) .gt. tolerance) then
        call lo_stop_gracefully(['Generating vector not unit length'], lo_exitcode_symmetry, __FILE__, __LINE__)
    end if

    ! Get the generators for the octahedral group, with axes defined by two vectors
    generators = 0.0_r8
    ! Octahedral:
    ! Find A orthogonal to B, and then C orthogonal to both
    C = lo_cross(A, B)

    m1(:, 1) = A
    m1(:, 2) = B
    m1(:, 3) = C
    m1 = lo_chop(m1, 1E-11_r8)
    m2 = lo_invert3x3matrix(m1)

    ! 1) Inversion
    generators(:, 1, 1) = [-1.0_r8, 0.0_r8, 0.0_r8]
    generators(:, 2, 1) = [0.0_r8, -1.0_r8, 0.0_r8]
    generators(:, 3, 1) = [0.0_r8, 0.0_r8, -1.0_r8]

    ! 2) 2-fold rotation around A
    generators(:, :, 2) = lo_rotation_matrix_from_axis_and_angle(A, lo_pi)

    ! 3) 2-fold rotation around B
    generators(:, :, 3) = lo_rotation_matrix_from_axis_and_angle(B, lo_pi)

    ! 4) 3-fold rotation around [111]
    v0 = A + B + C
    v0 = v0/norm2(v0)
    generators(:, :, 4) = lo_rotation_matrix_from_axis_and_angle(v0, lo_twopi/3.0_r8)

    ! 5) 2-fold rotation around [110]
    v0 = A + B
    v0 = v0/norm2(v0)
    generators(:, :, 5) = lo_rotation_matrix_from_axis_and_angle(v0, lo_pi)

    op = 0.0_r8
    fop = 0.0_r8
    op(:, :, 1:5) = generators(:, :, 1:5)
    do i = 1, 5
        m0 = op(:, :, i)
        m0 = matmul(matmul(m2, m0), m1)
        fop(:, :, i) = anint(m0)
    end do
    ctr = 5
    ! Fill out the rest?
    iterloop: do iter = 1, 100
        ii = ctr
        ol1: do i = 1, ii
        ol2: do j = 1, ii
            ! build combined operation?
            m0 = matmul(op(:, :, i), op(:, :, j))
            m0 = matmul(matmul(m2, m0), m1)
            m0 = anint(m0)
            ! Check if it exists
            do k = 1, ctr
                if (norm2(m0 - fop(:, :, k)) .lt. tolerance) then
                    cycle ol2
                end if
            end do
            ! Did not exist, add it and try again
            ctr = ctr + 1
            fop(:, :, ctr) = m0
            op(:, :, ctr) = matmul(matmul(m1, m0), m2)
            cycle ol1
        end do ol2
        end do ol1

        if (ctr .eq. ii) then
            ! We are done!
            exit iterloop
        end if
    end do iterloop

    ! Sanity check?
    if (ctr .ne. 48) call lo_stop_gracefully(['Did not find 48 octahedral operations'], lo_exitcode_symmetry, __FILE__, __LINE__)

    ! Clean up the operations? Probably a sensible idea.
    m1(:, 1) = A
    m1(:, 2) = B
    m1(:, 3) = C
    !m1=lo_chop(m1,1E-13_r8)
    m2 = lo_invert3x3matrix(m1)
    ! Convert to fractional?
    do i = 1, ctr
        m0 = op(:, :, i)
        m0 = matmul(matmul(m2, m0), m1)
        ! if ( norm2(m0-anint(m0)) .gt. tolerance ) then
        !     m0=matmul(m2,m1)
        !     write(*,*) m0(:,1)
        !     write(*,*) m0(:,2)
        !     write(*,*) m0(:,3)
        !     write(*,*) ''
        !     write(*,*) m1(:,1)
        !     write(*,*) m1(:,2)
        !     write(*,*) m1(:,3)
        !     write(*,*) ''
        !     write(*,*) m2(:,1)
        !     write(*,*) m2(:,2)
        !     write(*,*) m2(:,3)
        !
        !     call lo_stop_gracefully(['Failed generating octahedral group'],lo_exitcode_symmetry,__FILE__,__LINE__)
        ! endif
        m0 = anint(m0)
        op(:, :, i) = matmul(matmul(m1, m0), m2)
    end do
    op = lo_chop(op, 1E-13_r8)
end subroutine

!> generate the hexagonal operations from two vectors. A is the hexagonal axis.
subroutine hexagonal_operations_from_two_vectors(A, B, tolerance, op)
    !> first vector
    real(r8), dimension(3), intent(in) :: A
    !> second vector
    real(r8), dimension(3), intent(in) :: B
    !> tolerance
    real(r8), intent(in) :: tolerance
    !> point operations
    real(r8), dimension(3, 3, 24) :: op

    real(r8), dimension(3, 3, 3) :: generators
    real(r8), dimension(3, 3) :: m0
    integer :: iter, ctr, i, j, k, ii

    ! First check that the vectors are orthogonal
    if (abs(dot_product(A, B)) .gt. tolerance) then
        call lo_stop_gracefully(['Generating vectors not orthogonal'], lo_exitcode_symmetry, __FILE__, __LINE__)
    end if
    if (abs(norm2(A) - 1.0_r8) .gt. tolerance) then
        call lo_stop_gracefully(['Generating vector not unit length'], lo_exitcode_symmetry, __FILE__, __LINE__)
    end if
    if (abs(norm2(B) - 1.0_r8) .gt. tolerance) then
        call lo_stop_gracefully(['Generating vector not unit length'], lo_exitcode_symmetry, __FILE__, __LINE__)
    end if

    ! Get the generators for the octahedral group, with axes defined by two vectors
    generators = 0.0_r8

    ! Hexagonal:

    ! 1) Inversion
    generators(:, 1, 1) = [-1.0_r8, 0.0_r8, 0.0_r8]
    generators(:, 2, 1) = [0.0_r8, -1.0_r8, 0.0_r8]
    generators(:, 3, 1) = [0.0_r8, 0.0_r8, -1.0_r8]

    ! 2) 6-fold rotation around A
    generators(:, :, 2) = lo_rotation_matrix_from_axis_and_angle(A, lo_pi/3.0_r8)

    ! 3) 2-fold rotation around B
    generators(:, :, 3) = lo_rotation_matrix_from_axis_and_angle(B, lo_pi)

    op = 0.0_r8
    op(:, :, 1:3) = generators(:, :, 1:3)
    ctr = 3
    ! Fill out the rest?
    iterloop: do iter = 1, 100
        ii = ctr
        ol1: do i = 1, ii
        ol2: do j = 1, ii
            m0 = matmul(op(:, :, i), op(:, :, j))
            ! Check if it exists
            do k = 1, ctr
                if (norm2(m0 - op(:, :, k)) .lt. tolerance) then
                    cycle ol2
                end if
            end do
            ! Did not exist, add it and try again
            ctr = ctr + 1
            op(:, :, ctr) = m0
            cycle ol1
        end do ol2
        end do ol1

        if (ctr .eq. ii) then
            ! We are done!
            exit iterloop
        end if
    end do iterloop

    ! Sanity check?
    if (ctr .ne. 24) call lo_stop_gracefully(['Failed generating hexagonal group'], lo_exitcode_symmetry, __FILE__, __LINE__)

    ! Clean up the operations? Probably a sensible idea.
    op = lo_chop(op, 1E-13_r8)
end subroutine

!> count how many operations in set1 can be found in set2
function count_matching_matrices(set1, set2, tolerance) result(ctr)
    !> reference set
    real(r8), dimension(:, :, :), intent(in) :: set1
    !> trial set
    real(r8), dimension(:, :, :), intent(in) :: set2
    !> tolerance
    real(r8), intent(in) :: tolerance
    !> number of matches
    integer :: ctr

    integer :: i, j

    ctr = 0
    sl1: do i = 1, size(set1, 3)
    do j = 1, size(set2, 3)
        if (sum(abs(set1(:, :, i) - set2(:, :, j))) .lt. tolerance) then
            ctr = ctr + 1
            cycle sl1
        end if
    end do
    end do sl1
end function

!> cluster rotation error?
function cluster_transform_error(op, r) result(err)
    !> operations
    real(r8), dimension(:, :, :), intent(in) :: op
    !> points to transform
    real(r8), dimension(:, :), intent(in) :: r
    !> error
    real(r8) :: err

    real(r8), dimension(3, size(r, 2)) :: u
    real(r8) :: f0, f1
    integer :: i, j, n, iop

    err = 0.0_r8
    n = size(r, 2)
    do iop = 1, size(op, 3)
        ! Rotate cluster
        u = matmul(op(:, :, iop), r)
        ! Measure error
        do i = 1, n
            f1 = lo_huge
            do j = 1, n
                f0 = norm2(r(:, i) - u(:, j))
                f1 = min(f1, f0)
            end do
            err = err + f1
        end do
    end do
    err = err/real(size(r, 2)*size(op, 3), r8)
end function

end submodule
