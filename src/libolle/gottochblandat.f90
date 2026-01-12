module gottochblandat
!! A lot of helper functions, that did not logically fit somewhere else. Used to be called helpers, but it turns out it's a rather common name for a module and leads to issues, hence the weird name.
use konstanter, only: r8,i8,lo_pi,lo_twopi,lo_tol,lo_sqtol,lo_temperaturetol,lo_freqtol,lo_huge,&
                      lo_hugeint,lo_status,lo_exitcode_param,lo_exitcode_symmetry,lo_exitcode_blaslapack,&
                      lo_kb_hartree,lo_exitcode_baddim

implicit none
private
! Derived types
public :: lo_verletbox
public :: lo_array_of_array
! Interfaced things
public :: qsort
public :: tochar
public :: lo_return_unique
public :: lo_return_unique_indices
public :: lo_determ
public :: lo_flattentensor
public :: lo_outerproduct
public :: lo_transposetensor
public :: lo_chop
! Subroutines
public :: lo_complex_gram_schmidt
public :: lo_complex_hermitian_eigenvalues_eigenvectors
public :: lo_compress_equations
public :: lo_enforce_linear_constraints
public :: lo_fetch_tolerance
public :: lo_find_rotation_that_makes_strain_diagonal
public :: lo_get_axis_angles
public :: lo_identitymatrix
public :: lo_linear_least_squares
public :: lo_linspace
public :: lo_logspace
public :: lo_looptimer
public :: lo_make_coeffmatrix_tidy
public :: lo_make_complex_coeffmatrix_tidy
public :: lo_nullspace_coefficient_matrix
public :: lo_complex_nullspace_coefficient_matrix
public :: lo_real_nullspace_coefficient_matrix
public :: lo_permutations
public :: lo_points_on_sphere
public :: lo_progressbar
public :: lo_progressbar_init
public :: lo_put_function_on_new_axis
public :: lo_real_gram_schmidt
public :: lo_real_symmetric_eigenvalues_eigenvectors
public :: lo_return_tensor_transpositions
public :: lo_real_singular_value_decomposition
public :: lo_complex_singular_value_decomposition
public :: lo_stop_gracefully
public :: lo_symmetric_eigensystem_3x3matrix
public :: lo_transpositionmatrix
public :: lo_general_real_eigenvalues_eigenvectors
public :: lo_invert_real_matrix
public :: lo_real_pseudoinverse
public :: lo_triplegemm
public :: lo_untangle_one_tetrahedron
! Functions without lo_ thingy. Not good. Should fix.
public :: open_file
public :: walltime
! Functions
public :: lo_choplarge
public :: lo_classical_harmonic_oscillator_free_energy
public :: lo_clean_fractional_coordinates
public :: lo_count_lines_in_file
public :: lo_cross
public :: lo_does_file_exist
public :: lo_fermi
public :: lo_fermidirac
public :: lo_frobnorm
public :: lo_gauss
public :: lo_harmonic_oscillator_cv
public :: lo_harmonic_oscillator_entropy
public :: lo_harmonic_oscillator_free_energy
public :: lo_harmonic_oscillator_internal_energy
public :: lo_index_in_periodic_array
public :: lo_invert3x3matrix
public :: lo_sqrt3x3matrix
public :: lo_intsqrt
public :: lo_kmesh_density
public :: lo_linear_interpolation
public :: lo_linear_gradient_in_tetrahedron
public :: lo_mass_from_Z
public :: lo_lorentz
public :: lo_mean
public :: lo_negsqrt
public :: lo_planck
public :: lo_planck_deriv
public :: lo_planck_secondderiv
public :: lo_reciprocal_basis
public :: lo_rsquare
public :: lo_seconds_to_ddhhmmss
public :: lo_sgn
public :: lo_signed_tetrahedron_volume
public :: lo_sqnorm
public :: lo_stddev
public :: lo_trace
public :: lo_trapezoid_integration
public :: lo_trueNtimes
public :: lo_unflatten_2tensor
public :: lo_unflatten_3tensor
public :: lo_unflatten_4tensor
public :: lo_unsigned_tetrahedron_volume

!> Sort all kinds of different things
interface qsort
    module procedure quick_sort_real
    module procedure quick_sort_int
    module procedure quick_sortv
    module procedure quick_sortiv
    module procedure sortchar
end interface qsort

!> Take the union of different things
interface lo_return_unique
    module procedure lo_return_unique_characters
    module procedure lo_return_unique_integers
    module procedure lo_return_unique_doubles
    module procedure lo_return_unique_integer_columns
    module procedure lo_return_unique_double_columns
    module procedure lo_return_unique_double_matrices
end interface lo_return_unique

interface lo_return_unique_indices
    module procedure lo_return_unique_indices_real_vectors
    module procedure lo_return_unique_indices_integers
    module procedure lo_return_unique_indices_integer_vectors
end interface lo_return_unique_indices

!> determinant of 3x3 matrices
interface lo_determ
    module procedure lo_determ_int
    module procedure lo_determ_real
end interface lo_determ
!> convert something to characters
interface tochar
    module procedure tochar_int
    module procedure tochar_intarr
    module procedure tochar_real
    module procedure tochar_realarr
end interface
!> round and remove small numbers
interface lo_chop
    module procedure lo_chop_real
    module procedure lo_chop_complex
end interface
!> outer products
interface lo_outerproduct
    module procedure lo_real_outerproduct
    module procedure lo_complex_outerproduct
end interface
!> make tensors flat
interface lo_flattentensor
    module procedure lo_flatten_2tensor
    module procedure lo_flatten_3tensor
    module procedure lo_flatten_4tensor
end interface
interface lo_transposetensor
    module procedure lo_transpose_2tensor
    module procedure lo_transpose_3tensor
    module procedure lo_transpose_4tensor
end interface
!> rsquare and standard deviation
interface lo_rsquare
    module procedure lo_rsquare_1d
    module procedure lo_rsquare_3d
end interface
interface lo_stddev
    module procedure lo_stddev_1d
    module procedure lo_stddev_3d
end interface

!> simple array of array
type lo_array_of_array
    real(r8), dimension(:), allocatable :: buf_1d_r8
    contains
        procedure :: destroy=>destroy_array_of_arrays
end type

!> list of points in a verlet box
type lo_verletbox_box
    integer :: n=-lo_hugeint
    integer, dimension(:), allocatable :: ind
end type

!> minimal Verlet-box to generate distancelists and things like that
type lo_verletbox
    !> box divisions
    integer :: nx=-lo_hugeint,ny=-lo_hugeint,nz=-lo_hugeint
    !> lower bounds
    real(r8), dimension(3) :: rmin=lo_huge
    !> upper bounds
    real(r8), dimension(3) :: rmax=lo_huge
    !> scalefactor per dimension
    real(r8), dimension(3) :: ird=lo_huge
    !> boxes with points in them
    type(lo_verletbox_box), dimension(:,:,:), allocatable :: box
    contains
        !> stuff the particles into boxes
        procedure :: generate=>add_particles_in_boxes
        !> box-indices from a point
        procedure :: boxind=>boxind_from_coordinate
        !> locate index of a point
        procedure :: locate=>locate_index_of_point
end type

! Interfaces to submodule tensors:
interface
    module subroutine lo_return_tensor_transpositions(trm_pair,trm_triplet,trm_quartet,prm_pair,prm_triplet,prm_quartet)
        real(r8), dimension(9,9,2), intent(out), optional :: trm_pair
        real(r8), dimension(27,27,6), intent(out), optional :: trm_triplet
        real(r8), dimension(81,81,24), intent(out), optional :: trm_quartet
        integer, dimension(2,2), intent(out), optional :: prm_pair
        integer, dimension(3,6), intent(out), optional :: prm_triplet
        integer, dimension(4,24), intent(out), optional :: prm_quartet
    end subroutine
    module pure function lo_flatten_2tensor(m) result(fm)
        real(r8), dimension(3,3), intent(in) :: m
        real(r8), dimension(9) :: fm
    end function
    module pure function lo_flatten_3tensor(m) result(fm)
        real(r8), dimension(3,3,3), intent(in) :: m
        real(r8), dimension(27) :: fm
    end function
    module pure function lo_flatten_4tensor(m) result(fm)
        real(r8), dimension(3,3,3,3), intent(in) :: m
        real(r8), dimension(81) :: fm
    end function
    module pure function lo_unflatten_2tensor(fm) result(m)
        real(r8), dimension(9), intent(in) :: fm
        real(r8), dimension(3,3) :: m
    end function
    module pure function lo_unflatten_3tensor(fm) result(m)
        real(r8), dimension(27), intent(in) :: fm
        real(r8), dimension(3,3,3) :: m
    end function
    module pure function lo_unflatten_4tensor(fm) result(m)
        real(r8), dimension(81), intent(in) :: fm
        real(r8), dimension(3,3,3,3) :: m
    end function
    module pure function lo_transpose_2tensor(m,perm) result(pm)
        real(r8), dimension(3,3), intent(in) :: m
        integer, dimension(2), intent(in) :: perm
        real(r8), dimension(3,3) :: pm
    end function
    module pure function lo_transpose_3tensor(m,perm) result(pm)
        real(r8), dimension(3,3,3), intent(in) :: m
        integer, dimension(3), intent(in) :: perm
        real(r8), dimension(3,3,3) :: pm
    end function
    module pure function lo_transpose_4tensor(m,perm) result(pm)
        real(r8), dimension(3,3,3,3), intent(in) :: m
        integer, dimension(4), intent(in) :: perm
        real(r8), dimension(3,3,3,3) :: pm
    end function
end interface

! Interfaces to submodule sorting
interface
    module pure subroutine quick_sort_int(list,order)
        integer, dimension(:), intent(inout)  :: list
        integer, dimension(:), intent(inout), optional :: order
    end subroutine
    module pure subroutine quick_sortv(A,ind,tol)
        real(r8), dimension(:,:), intent(inout) :: A
        integer, dimension(:), intent(out), optional :: ind
        real(r8), intent(in) :: tol
    end subroutine
    module pure subroutine quick_sortiv(A,ind)
        integer, dimension(:,:), intent(inout) :: A
        integer, dimension(:), intent(out), optional :: ind
    end subroutine
    module pure subroutine quick_sort_real(list, order)
        real(r8), dimension(:), intent(inout)  :: list
        integer, dimension(:), intent(inout), optional :: order
    end subroutine
    module subroutine sortchar(StringArray, indexarray, CaseInsensitive)
        character (len = *), dimension (:), intent (inout) :: StringArray
        integer, dimension(:), intent(out), optional :: indexarray
        logical, intent (in), optional :: CaseInsensitive
    end subroutine
    module pure subroutine lo_return_unique_characters(a,u)
        character(len=*), dimension(:), intent(in) :: a
        character(len=len(a(1))), dimension(:), allocatable, intent(out) :: u
    end subroutine
    module pure subroutine lo_return_unique_integers(a,u,sort)
        integer, dimension(:), intent(in) :: a
        integer, dimension(:), allocatable, intent(out) :: u
        logical, intent(in), optional :: sort
    end subroutine
    module pure subroutine lo_return_unique_integer_columns(a,u)
        integer, dimension(:,:), intent(in) :: a
        integer, dimension(:,:), allocatable, intent(out) :: u
    end subroutine
    module pure subroutine lo_return_unique_double_columns(a,u,tol,ind)
        real(r8), dimension(:,:), intent(in) :: a
        real(r8), dimension(:,:), allocatable, intent(out) :: u
        real(r8), intent(in), optional :: tol
        integer, dimension(:), intent(out), optional :: ind
    end subroutine
    module pure subroutine lo_return_unique_double_matrices(a,u,tol)
        real(r8), dimension(:,:,:), intent(in) :: a
        real(r8), dimension(:,:,:), allocatable, intent(out) :: u
        real(r8), intent(in), optional :: tol
    end subroutine
    module pure subroutine lo_return_unique_doubles(a,u,tol)
        real(r8), dimension(:), intent(in) :: a
        real(r8), dimension(:), allocatable, intent(out) :: u
        real(r8), intent(in), optional :: tol
    end subroutine
    module subroutine lo_permutations(p,n)
        integer, dimension(:,:), allocatable, intent(out) :: p
        integer, intent(in) :: n
    end subroutine
    module subroutine lo_return_unique_indices_real_vectors(a,ind,redind,tol)
        real(r8), dimension(:,:), intent(in) :: a
        integer, dimension(:), allocatable, intent(out) :: ind
        integer, dimension(:), intent(out), optional :: redind
        real(r8), intent(in), optional :: tol
    end subroutine
    module subroutine lo_return_unique_indices_integers(a,ind,redind)
        integer, dimension(:), intent(in) :: a
        integer, dimension(:), allocatable, intent(out) :: ind
        integer, dimension(:), intent(out), optional :: redind
    end subroutine
    module subroutine lo_return_unique_indices_integer_vectors(a,ind,redind)
        integer, dimension(:,:), intent(in) :: a
        integer, dimension(:), allocatable, intent(out) :: ind
        integer, dimension(:), intent(out), optional :: redind
    end subroutine
end interface

! Interfaces to linalg submodule
interface
    module elemental function lo_negsqrt(x) result(y)
        real(r8), intent(in) :: x
        real(r8) :: y
    end function
    module elemental function lo_sgn(x) result(sgn)
        real(r8), intent(in) :: x
        integer :: sgn
    end function
    module pure function lo_unsigned_tetrahedron_volume(nodes) result(volume)
        real(r8), dimension(3,4), intent(in) :: nodes
        real(r8) :: volume
    end function
    module pure function lo_signed_tetrahedron_volume(nodes) result(volume)
        real(r8), dimension(3,4), intent(in) :: nodes
        real(r8) :: volume
    end function
    module pure function lo_frobnorm(m) result(nrm)
        real(r8), dimension(:,:), intent(in) :: m
        real(r8) :: nrm
    end function
    module pure function lo_trace(m) result(tr)
        real(r8), dimension(:,:), intent(in) :: m
        real(r8) :: tr
    end function
    module pure function lo_sqnorm(a) result(nrm)
        real(r8), dimension(3), intent(in) :: a
        real(r8) :: nrm
    end function
    module pure function lo_cross(b,c) result(a)
        real(r8), dimension(3), intent(in) :: b
        real(r8), dimension(3), intent(in) :: c
        real(r8), dimension(3) :: a
    end function
    module pure function lo_real_outerproduct(a,b) result(m)
        real(r8), dimension(3), intent(in) :: a
        real(r8), dimension(3), intent(in) :: b
        real(r8), dimension(3,3) :: m
    end function
    module pure function lo_complex_outerproduct(a,b) result(m)
        complex(r8), dimension(3), intent(in) :: a
        complex(r8), dimension(3), intent(in) :: b
        complex(r8), dimension(3,3) :: m
    end function
    module pure function lo_determ_real(a) result(det)
        real(r8), dimension(3,3), intent(in) :: a
        real(r8) :: det
    end function
    module pure function lo_determ_int(a) result(det)
        integer, dimension(3,3), intent(in) :: a
        integer :: det
    end function
    module pure function lo_invert3x3matrix(m) result(n)
        real(r8), dimension(3,3), intent(in) :: m
        real(r8), dimension(3,3) :: n
    end function
    module function lo_sqrt3x3matrix(matrix) result(matrix_sqrt)
        real(r8), dimension(3,3), intent(in) :: matrix
        real(r8), dimension(3,3) :: matrix_sqrt
    end function
    module pure subroutine lo_real_gram_schmidt(X)
        real(r8), dimension(:,:), intent(inout) :: X
    end subroutine
    module pure subroutine lo_complex_gram_schmidt(X)
        complex(r8), dimension(:,:), intent(inout) :: X
    end subroutine
    module pure subroutine lo_identitymatrix(m)
        real(r8), dimension(:,:), intent(out) :: m
    end subroutine
    module subroutine lo_real_singular_value_decomposition(A,S,U,V)
        real(r8), dimension(:,:), intent(in) :: A
        real(r8), dimension(:), allocatable, intent(out) :: S
        real(r8), dimension(:,:), allocatable, intent(out), optional :: U
        real(r8), dimension(:,:), allocatable, intent(out), optional :: V
    end subroutine
    module subroutine lo_complex_singular_value_decomposition(A,S,U,V)
        complex(r8), dimension(:,:), intent(in) :: A
        real(r8), dimension(:), allocatable, intent(out) :: S
        complex(r8), dimension(:,:), allocatable, intent(out), optional :: U
        complex(r8), dimension(:,:), allocatable, intent(out), optional :: V
    end subroutine
    module subroutine lo_enforce_linear_constraints(A,X,tolerance,nosquare)
        real(r8), dimension(:,:), intent(in) :: A
        real(r8), dimension(:), intent(inout) :: X
        real(r8), intent(in), optional :: tolerance
        logical, intent(in), optional :: nosquare
    end subroutine
    module pure subroutine lo_make_coeffmatrix_tidy(m,tolerance)
        real(r8), dimension(:,:), intent(inout) :: m
        real(r8), intent(in), optional :: tolerance
    end subroutine
    module pure subroutine lo_make_complex_coeffmatrix_tidy(m,tolerance)
        complex(r8), dimension(:,:), intent(inout) :: m
        real(r8), intent(in), optional :: tolerance
    end subroutine
    module subroutine lo_compress_equations(equations,n_compressed_equations,compressed_equations,trans,tolerance)
        real(r8), dimension(:,:), intent(in) :: equations
        integer, intent(out) :: n_compressed_equations
        real(r8), dimension(:,:), allocatable, intent(out) :: compressed_equations
        logical, intent(in) :: trans
        real(r8), intent(in), optional :: tolerance
    end subroutine
    module subroutine lo_nullspace_coefficient_matrix(invarM,coeffM,nvar,varind,tolerance)
        real(r8), dimension(:,:), intent(in) :: invarM
        real(r8), dimension(:,:), allocatable, intent(out) :: coeffM
        integer, intent(out) :: nvar
        integer, dimension(:), allocatable, intent(out), optional :: varind
        real(r8), intent(in), optional :: tolerance
    end subroutine
    module subroutine lo_complex_nullspace_coefficient_matrix(invarM,invariant_operations,coeff,nvar,varind,tolerance)
        complex(r8), dimension(:,:), intent(in), optional :: invarM
        complex(r8), dimension(:,:,:), intent(in), optional :: invariant_operations
        complex(r8), dimension(:,:), allocatable, intent(out) :: coeff
        integer, intent(out) :: nvar
        integer, dimension(:), allocatable, intent(out), optional :: varind
        real(r8), intent(in), optional :: tolerance
    end subroutine
    module subroutine lo_real_nullspace_coefficient_matrix(invarM,invariant_operations,coeff,nvar,varind,tolerance)
        real(r8), dimension(:,:), intent(in), optional :: invarM
        real(r8), dimension(:,:,:), intent(in), optional :: invariant_operations
        real(r8), dimension(:,:), allocatable, intent(out) :: coeff
        integer, intent(out) :: nvar
        integer, dimension(:), allocatable, intent(out), optional :: varind
        real(r8), intent(in), optional :: tolerance
    end subroutine
    module pure subroutine lo_transpositionmatrix(tm)
        real(r8), dimension(:,:), intent(out) :: tm
    end subroutine
    module subroutine lo_linear_least_squares(A,B,x,constraint_C,constraint_D,nconstraints,tolerance,subset,weights,gramified)
        real(r8), dimension(:,:), intent(in) :: A
        real(r8), dimension(:), intent(in) :: B
        real(r8), dimension(:), intent(out) :: x
        real(r8), dimension(:,:), intent(in), optional :: constraint_C
        real(r8), dimension(:), intent(in), optional :: constraint_D
        integer, intent(in), optional :: nconstraints
        real(r8), intent(in), optional :: tolerance
        integer, intent(in), dimension(:), optional :: subset
        real(r8), dimension(:), intent(in), optional :: weights
        logical, intent(in), optional :: gramified
    end subroutine
    module subroutine lo_real_symmetric_eigenvalues_eigenvectors(A,eigenvalues,eigenvectors,careful,tolerance,nzeros)
        real(r8), dimension(:,:), intent(in) :: A
        real(r8), dimension(:), intent(out) :: eigenvalues
        real(r8), dimension(:,:), intent(out), optional :: eigenvectors
        logical, intent(in), optional :: careful
        real(r8), intent(in), optional :: tolerance
        integer, intent(in), optional :: nzeros
    end subroutine
    module subroutine lo_symmetric_eigensystem_3x3matrix(matrix,eigenvalues,eigenvectors)
        real(r8), dimension(3,3), intent(in) :: matrix
        real(r8), dimension(3), intent(out) :: eigenvalues
        real(r8), dimension(3,3), intent(out) :: eigenvectors
    end subroutine
    module subroutine lo_complex_hermitian_eigenvalues_eigenvectors(A,eigenvalues,eigenvectors,careful,tolerance,nzeros)
        complex(r8), dimension(:,:), intent(in) :: A
        real(r8), dimension(:), intent(out) :: eigenvalues
        complex(r8), dimension(:,:), intent(out), optional :: eigenvectors
        logical, intent(in), optional :: careful
        real(r8), intent(in), optional :: tolerance
        integer, intent(in), optional :: nzeros
    end subroutine
    module subroutine lo_general_real_eigenvalues_eigenvectors(A,eigenvalues,vec_left,vec_right,orthogonal)
        real(r8), dimension(:,:), intent(in) :: A
        complex(r8), dimension(:), intent(out) :: eigenvalues
        real(r8), dimension(:,:), intent(out) :: vec_left,vec_right
        logical, intent(in), optional :: orthogonal
    end subroutine
    module subroutine lo_invert_real_matrix(A,iA)
       real(r8), dimension(:,:), intent(in) :: A
       real(r8), dimension(:,:), intent(out) :: iA
    end subroutine
    module subroutine lo_real_pseudoinverse(A,B,tolerance)
        real(r8), dimension(:,:), intent(in) :: A
        real(r8), dimension(:,:), intent(out) :: B
        real(r8), intent(in), optional :: tolerance
    end subroutine
    module subroutine lo_triplegemm(A,B,C,D,transa,transb,transc)
        real(r8), dimension(:,:), intent(in) :: A,B,C
        real(r8), dimension(:,:), intent(out) :: D
        character(len=1), intent(in), optional :: transa,transb,transc
    end subroutine
end interface

! Interfaces to calculus submodule
interface
    module subroutine lo_put_function_on_new_axis(x,y,xi,yi,sigma,preservenorm)
        real(r8), dimension(:), intent(in) :: x
        real(r8), dimension(:), intent(in) :: y
        real(r8), dimension(:), intent(in) :: xi
        real(r8), dimension(:), intent(out) :: yi
        real(r8), intent(in), optional :: sigma
        logical, intent(in), optional :: preservenorm
    end subroutine
    module pure function lo_linear_interpolation(x,y,xi,threshold) result(yi)
        real(r8), dimension(:), intent(in) :: x
        real(r8), dimension(:), intent(in) :: y
        real(r8), intent(in) :: xi
        real(r8) :: yi
        real(r8), intent(in), optional :: threshold
    end function
    module pure subroutine lo_points_on_sphere(r,randomseed)
        real(r8), dimension(:,:), intent(out) :: r
        real(r8), intent(in), optional :: randomseed
    end subroutine
    module pure subroutine lo_linspace(minv,maxv,x)
        real(r8), dimension(:), intent(out) :: x
        real(r8), intent(in) :: minv
        real(r8), intent(in) :: maxv
    end subroutine
    module pure subroutine lo_logspace(minv,maxv,x)
        real(r8), dimension(:), intent(out) :: x
        real(r8), intent(in) :: minv
        real(r8), intent(in) :: maxv
    end subroutine
    module pure function lo_trapezoid_integration(x,y) result(f)
        real(r8), dimension(:), intent(in) :: x
        real(r8), dimension(:), intent(in) :: y
        real(r8) :: f
    end function
    module elemental function lo_gauss(x,mu,sigma) result(g)
        real(r8), intent(in) :: x
        real(r8), intent(in) :: mu
        real(r8), intent(in) :: sigma
        real(r8) :: g
    end function
    module elemental function lo_lorentz(x,mu,sigma) result(l)
        real(r8), intent(in) :: x
        real(r8), intent(in) :: mu
        real(r8), intent(in) :: sigma
        real(r8) :: l
    end function
    module pure function lo_mean(x) result(m)
        real(r8), dimension(:), intent(in) :: x
        real(r8) :: m
    end function
    module pure function lo_stddev_1d(x) result(s)
        real(r8), dimension(:), intent(in) :: x
        real(r8) :: s
    end function
    module pure function lo_stddev_3d(x) result(s)
        real(r8), dimension(:,:,:), intent(in) :: x
        real(r8) :: s
    end function
    module pure function lo_rsquare_1d(values,predictions) result(rsq)
        real(r8), dimension(:), intent(in) :: values
        real(r8), dimension(:), intent(in) :: predictions
        real(r8) :: rsq
    end function
    module pure function lo_rsquare_3d(values,predictions) result(rsq)
        real(r8), dimension(:,:,:), intent(in) :: values
        real(r8), dimension(:,:,:), intent(in) :: predictions
        real(r8) :: rsq
    end function
    module pure function lo_linear_gradient_in_tetrahedron(corners,values,tolerance) result(gradient)
        real(r8), dimension(3,4), intent(in) :: corners
        real(r8), dimension(4), intent(in) :: values
        real(r8), intent(in) :: tolerance
        real(r8), dimension(3) :: gradient
    end function
end interface

! Interfaces to the physics submodule
interface
    module pure function lo_reciprocal_basis(a) result(b)
        real(r8), dimension(3,3), intent(in) :: a
        real(r8), dimension(3,3) :: b
    end function
    module pure function lo_kmesh_density(reciprocal_latticevectors,na,nq) result(density)
        real(r8), dimension(3,3), intent(in) :: reciprocal_latticevectors
        integer, intent(in) :: na
        integer, intent(in) :: nq
        integer, dimension(3) :: density
    end function
    module pure subroutine lo_get_axis_angles(basis,a,b,c,al,be,gm)
        real(r8), dimension(3,3), intent(in) :: basis
        real(r8), intent(out) :: a
        real(r8), intent(out) :: b
        real(r8), intent(out) :: c
        real(r8), intent(out) :: al
        real(r8), intent(out) :: be
        real(r8), intent(out) :: gm
    end subroutine
    module elemental function lo_fermidirac(x,mu,sigma) result(f)
        real(r8), intent(in) :: x
        real(r8), intent(in) :: mu
        real(r8), intent(in) :: sigma
        real(r8) :: f
    end function lo_fermidirac
    module elemental function lo_fermi(energy,efermi,temperature) result(f)
        real(r8), intent(in) :: energy
        real(r8), intent(in) :: efermi
        real(r8), intent(in) :: temperature
        real(r8) :: f
    end function lo_fermi
    module elemental function lo_planck(temperature,omega) result(n)
        real(r8), intent(in) :: temperature
        real(r8), intent(in) :: omega
        real(r8) :: n
    end function lo_planck
    module elemental function lo_planck_deriv(T,omega) result(dndt)
        real(r8), intent(in) :: T
        real(r8), intent(in) :: omega
        real(r8) :: dndt
    end function lo_planck_deriv
    module elemental function lo_planck_secondderiv(T,omega) result(ddnddt)
        real(r8), intent(in) :: T
        real(r8), intent(in) :: omega
        real(r8) :: ddnddt
    end function lo_planck_secondderiv
    module elemental function lo_harmonic_oscillator_free_energy(temp,omega) result(f)
        real(r8), intent(in) :: temp
        real(r8), intent(in) :: omega
        real(r8) :: f
    end function
    module elemental function lo_harmonic_oscillator_internal_energy(temp,omega) result(u)
        real(r8), intent(in) :: temp
        real(r8), intent(in) :: omega
        real(r8) :: u
    end function
    module elemental function lo_classical_harmonic_oscillator_free_energy(temp,omega) result(f)
        real(r8), intent(in) :: temp
        real(r8), intent(in) :: omega
        real(r8) :: f
    end function
    module elemental function lo_harmonic_oscillator_entropy(temp,omega) result(s)
        real(r8), intent(in) :: temp
        real(r8), intent(in) :: omega
        real(r8) :: s
    end function
    module elemental function lo_harmonic_oscillator_cv(temp,omega) result(cv)
        real(r8), intent(in) :: temp
        real(r8), intent(in) :: omega
        real(r8) :: cv
    end function
    module subroutine lo_untangle_one_tetrahedron(corner,energy,groupvelocity,degentol,permutation,success,error_interp_energy,error_interp_gradient)
        real(r8), dimension(3,4), intent(in) :: corner
        real(r8), dimension(:,:), intent(in) :: energy
        real(r8), dimension(:,:,:), intent(in) :: groupvelocity
        real(r8), intent(in) :: degentol
        integer, dimension(:,:), intent(out) :: permutation
        integer, intent(out) :: success
        real(r8), dimension(:), intent(out), optional :: error_interp_energy
        real(r8), dimension(:), intent(out), optional :: error_interp_gradient
    end subroutine
    module subroutine lo_find_rotation_that_makes_strain_diagonal( ref_lattice,lattice,rotation,strain,guess )
        real(r8), dimension(3,3), intent(in) :: ref_lattice
        real(r8), dimension(3,3), intent(in) :: lattice
        real(r8), dimension(3,3), intent(out) :: rotation
        real(r8), dimension(3,3), intent(out) :: strain
        real(r8), dimension(3,3), intent(in), optional :: guess
    end subroutine
    module elemental function lo_mass_from_Z(Z) result(mass)
        integer, intent(in) :: Z
        real(r8) :: mass
    end function
end interface

! Interfaces to the boxes submodule
interface
    module subroutine boxind_from_coordinate(vb,r,bi,bj,bk)
        class(lo_verletbox), intent(in) :: vb
        real(r8), dimension(3), intent(in) :: r
        integer, intent(out) :: bi,bj,bk
    end subroutine
    module function locate_index_of_point(vb,points,r,singlebox) result(ind)
        class(lo_verletbox), intent(in) :: vb
        real(r8), dimension(:,:), intent(in) :: points
        real(r8), dimension(3), intent(in) :: r
        logical, intent(in), optional :: singlebox
        integer :: ind
    end function
    module subroutine add_particles_in_boxes(vb,r,ndim)
        class(lo_verletbox), intent(out) :: vb
        real(r8), dimension(:,:), intent(in) :: r
        integer, dimension(3), intent(in) :: ndim
    end subroutine
end interface

contains

!> Wallclock time. Wrapper around system_clock, I think.
function walltime() result(time)
    !> elapsed time
    real(r8) :: time

    logical, save :: firstcall=.true.
    integer(i8), save :: countrate=-1,countmax=-1
    integer(i8) :: i

    if ( firstcall ) then
        call system_clock(count_rate=countrate)
        call system_clock(count_max=countmax)
    endif

    call system_clock(i)
    time=real(i,r8)/real(countrate,r8)
end function

!> helper function that returns true nrep times when i cycles over nloop values. I use it for progressbars, so that if it is a very long and fast loop printing the progress bar does not become a bottleneck.
pure function lo_trueNtimes(i,nrep,nloop) result(tf)
    !> loop counter
    integer, intent(in) :: i
    !> number of times to say true
    integer, intent(in) :: nrep
    !> total number of loops
    integer, intent(in) :: nloop
    logical :: tf

    integer :: j
    ! trivial case
    if ( i .eq. nloop ) then
        tf=.false.
        return
    endif
    if ( nloop .le. nrep ) then
        tf=.true.
        return
    endif
    ! normal case
    j=max(nloop/nrep,1)
    if ( mod(i,j) .eq. 0 ) then
        tf=.true.
    else
        tf=.false.
    endif
end function

!> convert one int to a character
pure function tochar_int(i,padding) result(s)
    !> integer to convert
    integer, intent(in) :: i
    !> pad the integer? Positive number means zer-padding, negative means pad with whitespace
    integer, intent(in), optional :: padding
    !> resulting string
    character(len=:), allocatable :: s

    character(len=range(i)) :: tmp
    character(len=range(i)+5) :: ttmp
    integer :: j,k

    if ( present(padding) ) then
        write(tmp,'(i0)') i ! get i to a string
        ! build a string that contains the padding
        do j=1,range(i)+5
            if ( padding > 0 ) then
                ttmp(j:j)='0'
            else
                ttmp(j:j)=' '
            endif
        enddo
        ! wrap it together to a nice string
        k=len(trim(adjustl(tmp))) ! how many digits where relevant
        s=ttmp(1:abs(padding)-k)//trim(adjustl(tmp))
    else
        ! much easier
        write(tmp,'(i0)') i
        s=trim(adjustl(tmp))
    endif
end function

!> convert array of integers to a character
pure function tochar_intarr(i,padding) result(s)
    integer, dimension(:), intent(in) :: i
    integer, intent(in), optional :: padding
    character(len=:), allocatable :: s

    character(len=size(i,1)*100) :: tmp
    integer :: j
    tmp=''
    do j=1,size(i,1)
        if ( present(padding) ) then
            tmp=trim(tmp)//' '//tochar_int(i(j),padding)//' '
        else
            tmp=trim(tmp)//' '//tochar_int(i(j))//' '
        endif
    enddo
    s=trim(tmp)
end function

!> convert a real to a character
pure function tochar_real(f,ndecimals,frmt) result(s)
    !> number to convert
    real(r8), intent(in) :: f
    !> how many decimals
    integer, intent(in), optional :: ndecimals
    !> maybe format instead
    character(len=*), intent(in), optional :: frmt
    !> resulting string
    character(len=:), allocatable :: s

    integer :: ndec
    character(len=100) :: tmp

    ! how many decimal places
    if ( present(ndecimals) ) then
        ndec=ndecimals
    else
        ndec=5
    endif

    ! convert to string
    if ( present(frmt) ) then
        write(tmp,trim(frmt)) f
        s=trim(tmp)
    else
        write(tmp,"(F30."//tochar_int(ndec)//")") f
        s=trim(adjustl(tmp))
    endif
end function

!> convert an array of reals to a character
pure function tochar_realarr(f,ndecimals,frmt) result(s)
    !> number to convert
    real(r8), dimension(:), intent(in) :: f
    !> how many decimals
    integer, intent(in), optional :: ndecimals
    !> maybe format instead
    character(len=*), intent(in), optional :: frmt
    !> resulting string
    character(len=:), allocatable :: s

    integer :: j,ndec
    character(len=size(f,1)*50) :: tmp

    ! how many decimal places
    if ( present(ndecimals) ) then
        ndec=ndecimals
    else
        ndec=5
    endif

    tmp=''
    do j=1,size(f,1)
        if ( present(frmt) ) then
            tmp=trim(tmp)//tochar_real(f(j),frmt=frmt)
        else
            tmp=trim(tmp)//' '//tochar_real(f(j),ndec)//' '
        endif
    enddo
    if ( present(frmt) ) then
        s=trim(tmp)
    else
        s=trim(adjustl(tmp))
    endif
end function

!> convert seconds to hh:mm:ss
pure function lo_seconds_to_ddhhmmss(t) result(ddhhmmss)
    !> time, in seconds
    real(r8), intent(in) :: t
    !> formatted string
    character(len=:), allocatable :: ddhhmmss

    integer :: total,days,hours,minutes,seconds
    total=ceiling(abs(t)+lo_tol) ! at least a second
    days=total/24/3600
    total=total-days*24*3600
    hours=total/3600
    total=total-hours*3600
    minutes=total/60
    seconds=total-minutes*60
    if ( days .gt. 0 ) then
        ddhhmmss=tochar(days,2)//':'//tochar(hours,2)//':'//tochar(minutes,2)//':'//tochar(seconds,2)
    else
        ddhhmmss=tochar(hours,2)//':'//tochar(minutes,2)//':'//tochar(seconds,2)
    endif
end function

!> Set all the possible tolerances from a global tolerance
pure subroutine lo_fetch_tolerance(tol,basis,realspace_cart_tol,realspace_fractional_tol,relative_tol,reciprocal_cart_tol,reciprocal_fractional_tol,temperature_tol,radian_tol,degree_tol,freq_tol,squared)
    !> global tolerance
    real(r8), intent(in) :: tol
    !> metric for cartesian distances, realspace
    real(r8), dimension(3,3), intent(in), optional :: basis
    !> this tolerance in different metrics
    real(r8), intent(out), optional :: realspace_cart_tol
    real(r8), intent(out), optional :: realspace_fractional_tol
    real(r8), intent(out), optional :: relative_tol
    real(r8), intent(out), optional :: reciprocal_cart_tol
    real(r8), intent(out), optional :: reciprocal_fractional_tol
    real(r8), intent(out), optional :: temperature_tol
    real(r8), intent(out), optional :: radian_tol
    real(r8), intent(out), optional :: degree_tol
    real(r8), intent(out), optional :: freq_tol
    logical, intent(in), optional :: squared

    real(r8) :: f0
    real(r8), dimension(3,3) :: recbasis
    logical :: sq

    if ( present(squared) ) then
        sq=squared
    else
        sq=.false.
    endif

    if ( present(realspace_cart_tol) ) then
        ! Realspace distances, in A
        realspace_cart_tol=tol
        if ( sq ) realspace_cart_tol=realspace_cart_tol**2
    endif

    if ( present(realspace_fractional_tol) ) then
        ! Realspace distances in fractional coordinates
        f0=( norm2(basis(:,1))+norm2(basis(:,2))+norm2(basis(:,3)) )/3.0_r8
        realspace_fractional_tol=tol/f0
        if ( sq ) realspace_fractional_tol=realspace_fractional_tol**2
    endif

    if ( present(relative_tol) ) then
        ! For relative quantities
        relative_tol=tol*1E-2_r8
        if ( sq ) relative_tol=relative_tol**2
    endif

    if ( present(reciprocal_cart_tol) ) then
        ! in reciprocal space
        reciprocal_cart_tol=tol/abs(lo_determ(basis))
        if ( sq ) reciprocal_cart_tol=reciprocal_cart_tol**2
    endif

    if ( present(reciprocal_fractional_tol) ) then
        ! in reciprocal fractional space
        recbasis=lo_reciprocal_basis(basis)
        f0=( norm2(recbasis(:,1))+norm2(recbasis(:,2))+norm2(recbasis(:,3)) )/3.0_r8
        reciprocal_fractional_tol=tol/f0
        if ( sq ) reciprocal_fractional_tol=reciprocal_fractional_tol**2
    endif

    if ( present(temperature_tol) ) then
        temperature_tol=tol*1E2_r8
        if ( sq ) temperature_tol=temperature_tol**2
    endif

    if ( present(radian_tol) ) then
        radian_tol=tol*10.0_r8*180.0_r8/lo_pi
        if ( sq ) radian_tol=radian_tol**2
    endif

    if ( present(degree_tol) ) then
        degree_tol=tol*10.0_r8
        if ( sq ) degree_tol=degree_tol**2
    endif

    if ( present(freq_tol) ) then
        freq_tol=tol*1E-4_r8
        if ( sq ) freq_tol=freq_tol**2
    endif
end subroutine

!> Find an available unit to write to @todo Retire in favour of the intrinsic
integer function unit_number(filename)
    !> The filename
    character(len=*),intent(in) :: filename

    integer :: unit, i
    logical :: file_open
    i=100
    inquire(file=filename,number=unit,opened=file_open)
    if ( file_open .eqv. .true. ) then
        unit_number=unit
    else
        do
            inquire(unit=i,opened=file_open)
            if ( file_open .eqv. .false. ) then
                unit_number=i
                exit
            endif
            i=i+1
            if ( i == 1000 ) then
                write(*,*) 'No available units between 100 and 1000. In my opinion, too many open files.'
                stop
            endif
        enddo
    endif
end function unit_number

!> Opens a file, returns the unit @todo This should also be fixed to the intrinsinc or something.
integer function open_file(inut,filename)
    !> The file
    character(len=*),intent(in) :: filename
    !> 'in' or 'out'
    character(len=*),intent(in) :: inut
    !
    integer :: unit, status
    !
    status=0
    !
    unit=unit_number(filename)
    select case(trim(inut))
        case('in')
            open(unit=unit, file=filename, status='old', action='read',iostat=status)
            open_file=unit
        case('out')
            open(unit=unit, file=filename, status='replace', action='write',iostat=status)
            open_file=unit
        case default
            call lo_stop_gracefully(["Please open the file with either 'in' or 'out', not something in between"],lo_exitcode_param)
    end select
    !
    if ( status .ne. 0 ) then
        write(*,*) 'Could not open ',filename
        stop
    endif
end function open_file

!> Counts the number of lines in a file
function lo_count_lines_in_file(filename,unsafe) result(nlines)
    !> the filename
    character(len=*), intent(in) :: filename
    !> unsafe way, does a system call to `wc -l`
    logical, intent(in), optional :: unsafe
    !> number of lines
    integer :: nlines

    integer :: n,u
    character(len=1000):: cmd
    logical :: us

    ! the unsafe way that call wc, or the safe fortran way
    if ( present(unsafe) ) then
        us=unsafe
    else
        us=.false.
    endif

    ! Slightly unsafe, but fast way. Relies on 'wc' being on the path, which I have no way of verifying.
    if ( us ) then
        cmd = "wc -l "//trim(filename)//" > nlines.txt"
        call execute_command_line(trim(cmd))
        u=open_file('in','nlines.txt')
            read(u,*) n
        close(u)
        cmd = 'rm nlines.txt'
        call execute_command_line(trim(cmd))
    else
        ! A slightly nicer way, without the system calls, but way slower.
        n=0
        u=open_file('in',trim(filename))
            do
                read(u,*,END=10)
                n=n+1
            enddo
        10 close(u)
    endif
    nlines=n
end function

!> Test if a file exist. Just wrapped it since I never remember the syntax.
function lo_does_file_exist(filename) result(isthere)
    !> filename to check
    character(len=*), intent(in) :: filename
    !> does it exist
    logical :: isthere

    inquire(file=trim(filename), exist=isthere)
end function

!> Print loop timings
subroutine lo_looptimer(msg,t_start,t_curr,i,n)
    !> message
    character(len=*), intent(in) :: msg
    !> starting time for the loop
    real(r8), intent(in) :: t_start
    !> current time
    real(r8), intent(in) :: t_curr
    !> counter
    integer, intent(in) :: i
    !> total count
    integer, intent(in) :: n

    real(r8) :: t_total,t_elapsed,t_remaining,percent_done
    ! elapsed time
    t_elapsed=t_curr-t_start
    ! estimate total and remaining time
    if ( i .gt. 1 ) then
        t_total=t_elapsed*(n*1.0_r8-1.0_r8)/(i*1.0_r8-1.0_r8)
        t_remaining=t_total-t_elapsed
    else
        t_total=0.0_r8
        t_remaining=0.0_r8
    endif
    percent_done=100*(i*1.0_r8-1.0_r8)/(n*1.0_r8-1.0_r8)

    ! print this
    write(*,"(1X,A)") trim(adjustl(msg))//', remaining time: '//lo_seconds_to_ddhhmmss(t_remaining)//&
               ', estimated total time: '//lo_seconds_to_ddhhmmss(t_total)//'  ('//tochar(percent_done,1)//'%)'
end subroutine

!> Initializes the progress bar
subroutine lo_progressbar_init()
    ! gfortran and ifort has slightly different ideas how things should be done.
#if ifortprogressbar
    write(*,*) ''
#else
    ! do nothing
#endif
end subroutine lo_progressbar_init

!> Prints a progressbar to stdout, with both some kind of bar thing and percentages. Looks pretty and is oddly satisfying. If you have a really fast loop, don't use this since it will become a bottleneck. Also don't call this too often, since some things, such as SLURM/PBS output files, don't really like this. Use it like this:
!>
!> ```fortran
!> call lo_progressbar_init()
!> do i=1,n
!>     ! do stuff
!>     call lo_progressbar('doing stuff',i,n)
!> enddo
!> ```
subroutine lo_progressbar(message,j,totn,elapsed_time)
    !> Current index
    integer, intent(in) :: j
    !> Max index
    integer, intent(in) :: totn
    !> What is being calculated
    character(len=*), intent(in) :: message
    !> Optionally, how long time it has taken
    real(r8), intent(in), optional :: elapsed_time
    !
    real(r8) :: f0
    integer :: i,k
    character(len=50) :: bar
    character(len=40) :: msg
    character(len=15) :: time
    !
    bar="     % |                                        |"
    msg="                                        "
    ! write the percentage in the bar
    if ( totn .gt. 1 ) then
        f0=100.0_r8*(j-1.0_r8)/(totn*1.0_r8-1.0_r8)
    else
        f0=100.0_r8
    endif
    write(unit=bar(1:5),fmt="(F5.1)") f0
    ! fill in the bar
    k=ceiling(40*f0/100)
    do i = 1, k
        bar(8+i:8+i)="="
    enddo
    ! add the message
    do i=1,min(len_trim(message),40)
        msg(i:i)=message(i:i)
    enddo
    ! maybe the elapsed time?
    if ( present(elapsed_time) ) then
        time=' '//tochar(elapsed_time)//'s'
    else
        time=''
    endif
    ! Print the progress bar. Different compilers like it differently.
#if ifortprogressbar
    write(unit=6,fmt="(1X,a1,a1,a40,a50)") '+',char(13),msg,bar,trim(time)
#elif gfortranprogressbar
    if ( j .lt. totn ) then
        write(*,'(1X,a,a,a,a)',advance='no') char(13),msg,bar,trim(time)
    else
        write(*,'(1X,a,a,a,a)') char(13),msg,bar,trim(time)
    endif
#elif clusterprogressbar
    if ( j .eq. totn ) then
        write(*,'(1X,a,a,a)') msg,bar,trim(time)
    endif
#else
    ! Boring version, just dump line after line.
    if ( totn .gt. 100 ) then
        if ( mod(j,totn/100) .eq. 0 ) write(*,*) msg,real(f0),'%'
    else
        write(*,*) msg,real(f0),'%'
    endif
#endif
end subroutine lo_progressbar

!> Returns an index in a periodic array of period d. For example: if I have an array of 5 periodic values, and I want to access element 7, this function will return 2.
pure function lo_index_in_periodic_array(i,d) result(j)
    !> index
    integer, intent(in) :: i
    !> period
    integer, intent(in) :: d
    !> periodicity-adjusted index
    integer :: j

    j=i
    if ( j .gt. d ) then
        do while ( j .gt. d )
            j=j-d
        enddo
    elseif ( j .lt. 1 ) then
        do while ( j .lt. 1 )
            j=j+d
        enddo
    else
        j=i
    endif
end function lo_index_in_periodic_array

!> Clean up fractional coordinates in a way that I like, so that they become 0-1, inclusive 0 not inclusive 1. Some quick testing suggested this was faster than mod if I am close to 0-1, which almost always is the case. Also means that you should never use this for large numbers.
elemental function lo_clean_fractional_coordinates(x,tol) result(y)
    !> unclean coordinate
    real(r8), intent(in) :: x
    !> tolerance
    real(r8), intent(in), optional :: tol
    !> the clean coordinate
    real(r8) :: y

    real(r8) :: dl
    if ( present(tol) ) then
        dl=tol
    else
        dl=lo_sqtol
    endif
    y=x
    do
        if ( abs(y) .lt. dl ) y=0.0_r8
        if ( abs(y-1.0_r8) .lt. dl ) y=0.0_r8
        if ( y .gt. 1.0_r8 ) y=y-1.0_r8
        if ( y .lt. 0.0_r8 ) y=y+1.0_r8
        if ( y .gt. -dl .and. y .lt. 1.0_r8 ) exit
    enddo
end function

!> Chop number to a short list of well-defined numbers, within a tolerance. Contains all n/m with m from 1 to 10 and n<m, and sqrt(3)/2. No problem adding more if that is needed. Only tests between -1 and 1
elemental function lo_chop_real(x,tol) result(y)
    !> number to be chopped
    real(r8), intent(in) :: x
    !> tolerance
    real(r8), intent(in) :: tol
    !> the chopped number
    real(r8) :: y

    integer :: i
    real(r8), dimension(69), parameter :: welldefined=[&
          0.0000000000000000_r8,&
          0.1000000000000000_r8,&
         -0.1000000000000000_r8,&
          0.1111111111111111_r8,&
         -0.1111111111111111_r8,&
          0.1250000000000000_r8,&
         -0.1250000000000000_r8,&
          0.1428571428571428_r8,&
         -0.1428571428571428_r8,&
          0.1666666666666667_r8,&
         -0.1666666666666667_r8,&
          0.2000000000000000_r8,&
         -0.2000000000000000_r8,&
          0.2222222222222222_r8,&
         -0.2222222222222222_r8,&
          0.2500000000000000_r8,&
         -0.2500000000000000_r8,&
          0.2857142857142857_r8,&
         -0.2857142857142857_r8,&
          0.3000000000000000_r8,&
         -0.3000000000000000_r8,&
          0.3333333333333333_r8,&
         -0.3333333333333333_r8,&
          0.3750000000000000_r8,&
         -0.3750000000000000_r8,&
          0.4000000000000000_r8,&
         -0.4000000000000000_r8,&
          0.4285714285714285_r8,&
         -0.4285714285714285_r8,&
          0.4444444444444444_r8,&
         -0.4444444444444444_r8,&
          0.5000000000000000_r8,&
         -0.5000000000000000_r8,&
          0.5555555555555556_r8,&
         -0.5555555555555556_r8,&
          0.5714285714285714_r8,&
         -0.5714285714285714_r8,&
          0.6000000000000000_r8,&
         -0.6000000000000000_r8,&
          0.6250000000000000_r8,&
         -0.6250000000000000_r8,&
          0.6666666666666667_r8,&
         -0.6666666666666667_r8,&
          0.7000000000000000_r8,&
         -0.7000000000000000_r8,&
          0.7071067811865475_r8,&
         -0.7071067811865475_r8,&
          0.7142857142857143_r8,&
         -0.7142857142857143_r8,&
          0.7500000000000000_r8,&
         -0.7500000000000000_r8,&
          0.7777777777777778_r8,&
         -0.7777777777777778_r8,&
          0.8000000000000000_r8,&
         -0.8000000000000000_r8,&
          0.8333333333333334_r8,&
         -0.8333333333333334_r8,&
          0.8571428571428571_r8,&
         -0.8571428571428571_r8,&
          0.8660254037844386_r8,&
         -0.8660254037844386_r8,&
          0.8750000000000000_r8,&
         -0.8750000000000000_r8,&
          0.8888888888888888_r8,&
         -0.8888888888888888_r8,&
          0.9000000000000000_r8,&
         -0.9000000000000000_r8,&
          1.0000000000000000_r8,&
         -1.0000000000000000_r8]

    ! my list of well defined values. Could very well grow a little bit.
    ! starts with numbers that are common in symmetry operations
    y=x
    do i=1,69
        if ( abs(y-welldefined(i)) .lt. tol ) then
            y=welldefined(i)
            return
        endif
    enddo
end function

!> cutoff at small numbers for complex numbers
elemental function lo_chop_complex(x,tol) result(y)
    !> number to be chopped
    complex(r8), intent(in) :: x
    !> tolerance
    real(r8), intent(in) :: tol
    !> chopped number
    complex(r8) :: y

    real(r8) :: re,im
    re=real(x)
    im=aimag(x)
    if ( abs(re) .lt. tol ) re=0.0_r8
    if ( abs(im) .lt. tol ) im=0.0_r8
    y=cmplx(re,im,r8)
end function

!> Slower chopping function, identifies integers and common sqrt things, and multiples of those. Not that fast perhaps, but quite useful.
elemental function lo_choplarge(x,tol) result(y)
    !> number to be chopped
    real(r8), intent(in) :: x
    !> tolerance
    real(r8), intent(in) :: tol
    !> the chopped number
    real(r8) :: y

    integer :: i
    real(r8), dimension(51), parameter :: welldefined=[&
-10.000000000000000000_r8, &
-9.000000000000000000_r8, &
-8.660254037844385522_r8, &
-8.000000000000000000_r8, &
-7.794228634059947147_r8, &
-7.000000000000000000_r8, &
-6.928203230275508773_r8, &
-6.062177826491070398_r8, &
-6.000000000000000000_r8, &
-5.196152422706632024_r8, &
-5.000000000000000000_r8, &
-4.500000000000000000_r8, &
-4.330127018922192761_r8, &
-4.000000000000000000_r8, &
-3.500000000000000000_r8, &
-3.464101615137754386_r8, &
-3.000000000000000000_r8, &
-2.598076211353316012_r8, &
-2.500000000000000000_r8, &
-2.000000000000000000_r8, &
-1.732050807568877193_r8, &
-1.500000000000000000_r8, &
-1.000000000000000000_r8, &
-0.866025403784438597_r8, &
-0.500000000000000000_r8, &
 0.000000000000000000_r8, &
 0.500000000000000000_r8, &
 0.866025403784438597_r8, &
 1.000000000000000000_r8, &
 1.500000000000000000_r8, &
 1.732050807568877193_r8, &
 2.000000000000000000_r8, &
 2.500000000000000000_r8, &
 2.598076211353316012_r8, &
 3.000000000000000000_r8, &
 3.464101615137754386_r8, &
 3.500000000000000000_r8, &
 4.000000000000000000_r8, &
 4.330127018922192761_r8, &
 4.500000000000000000_r8, &
 5.000000000000000000_r8, &
 5.196152422706632024_r8, &
 6.000000000000000000_r8, &
 6.062177826491070398_r8, &
 6.928203230275508773_r8, &
 7.000000000000000000_r8, &
 7.794228634059947147_r8, &
 8.000000000000000000_r8, &
 8.660254037844385522_r8, &
 9.000000000000000000_r8, &
10.000000000000000000_r8]

    ! my list of well defined values. Could very well grow a little bit.
    ! starts with numbers that are common in symmetry operations
    y=x
    do i=1,51
        if ( abs(y-welldefined(i)) .lt. tol ) then
            y=welldefined(i)
            return
        endif
    enddo
end function

!> sqrt of an integer. Only works for those that actually have integer square roots
elemental function lo_intsqrt(i) result(j)
    !> integer to take root of
    integer, intent(in) :: i
    !> square root of i
    integer :: j

    ! list of pre-calculated squared integers. Can fill it out once I need more
    integer, dimension(50), parameter :: defi=[1,4,9,16,25,36,49,64,81,100,121,144,169,&
        196,225,256,289,324,361,400,441,484,529,576,625,676,729,784,841,900,961,1024,1089,&
        1156,1225,1296,1369,1444,1521,1600,1681,1764,1849,1936,2025,2116,2209,2304,2401,2500]
    integer :: l
    ! just to make sure it breaks if you send in a bad integer
    j=-lo_hugeint
    ! compare with pre-calculated squares
    do l=1,50
        if ( i .eq. defi(l) ) then
            j=l
            return
        endif
    enddo
end function

!> kill the program gracefully, with an exit code and message. Optionally say what line in what file that called it.
subroutine lo_stop_gracefully(msg,exitcode,filename,line)
    !> message
    character(len=*), dimension(:), intent(in) :: msg
    !> what exit code to give?
    integer, intent(in) :: exitcode
    !> which file?
    character(len=*), intent(in), optional :: filename
    !> which line?
    integer, intent(in), optional :: line

    integer :: i

    write(*,*) ''
    write(*,*) 'ERROR'
    select case(exitcode)
    case(1)
        write(*,*) 'exit code 1: unexpected dimensions'
    case(2)
        write(*,*) 'exit code 1: blas/lapack returned nonzero exitcode'
    case(3)
        write(*,*) 'exit code 3: unphysical value detected'
    case(4)
        write(*,*) 'exit code 4: symmetry error'
    case(5)
        write(*,*) 'exit code 5: bad parameters sent to routine'
    case(6)
        write(*,*) 'exit code 6: I/O error'
    case(7)
        write(*,*) 'exit code 7: MPI error'
    end select
    write(*,*) ''
    do i=1,size(msg)
        write(*,*) trim(msg(i))
    enddo
    write(*,*) ''
    if ( present(filename) ) write(*,*) '    occurs in file: ',filename
    if ( present(line) )     write(*,*) '    occurs on line: ',tochar(line)

    ! Seems I can not use a variable as the errorcode with Ifort. No, worries,
    ! got a brute force solution:
    select case(exitcode)
    case(1)
        error stop 1
    case(2)
        error stop 2
    case(3)
        error stop 3
    case(4)
        error stop 4
    case(5)
        error stop 5
    case(6)
        error stop 6
    case(7)
        error stop 7
    case default
        error stop
    end select
end subroutine

subroutine destroy_array_of_arrays(self)
    class(lo_array_of_array), intent(inout) :: self

    if ( allocated(self%buf_1d_r8) ) deallocate(self%buf_1d_r8)
end subroutine

end module
