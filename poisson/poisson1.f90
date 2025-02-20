! Poisson problem on a unit square with Dirichlet boundary conditions
program poisson1
 use kind_defs_m
 use mesh_m
 use meshgen_m
 use problem_m
 use system_m
 use hsl_ma57_m
 use poisson_elements_m
 use io_utils_m
 use math_defs_m
 implicit none

! constants
 integer, parameter :: &
   uintpl = 8,         & ! scalar interpolation
   gauss = 3,          & ! 3x3 integration of quads
   gaussb = 3,         & ! 3 point integration of boundary elements
   nx=20,              & ! number of elements in x
   ny=20,              & ! number of elements in y
   funcnr=1              ! function number for the right-hand side
 real(dp) :: &
   alpha = 1._dp         ! diffusion coefficient
 real(dp) :: beta = 0.5_dp ! rhs parameter

! definitions
 type(meshgen_options_t) :: meshgen_options
 type(mesh_t) :: mesh
 type(input_probdef_t) :: input_probdef
 type(problem_t) :: problem
 type(sysmatrix_t) :: sysmatrix
 type(sysvector_t) :: sol, rhsd
 type(coefficients_t) :: coefficients
 type(subscript_t) :: subscript
 character(len=30) :: filename
 namelist /comppar/ alpha, beta
 read ( unit=*, nml=comppar )
! fill coefficients
 call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )
 coefficients%i = 0
 coefficients%i(1) = uintpl
 coefficients%i(10:12) = (/ gauss, gaussb, funcnr /)
 coefficients%r(1) = alpha
 coefficients%r(2:) = 0
 coefficients%func => func
! create mesh
 meshgen_options%elshape = 6
 meshgen_options%nx = nx
 meshgen_options%ny = ny
 call quadrilateral2d ( mesh, meshgen_options )
 call fill_mesh_parts ( mesh )
! problem definition
 call create_input_probdef ( mesh, input_probdef )
 input_probdef%elementdof(1)%a = (/1,1,1,1,1,1,1,1,1/)
 call define_essential ( mesh, input_probdef, curve1=1, curve2=4 )
 call problem_definition ( input_probdef, mesh, problem )
! create system vectors (solution and right-hand side)
 call create_sysvector ( problem, sol )
 call create_sysvector ( problem, rhsd )
 call create_subscript ( mesh, problem, subscript)
! fill solution vector with essential boundary conditions
 call fill_sysvector ( mesh, problem, sol, &
   curve1=1, curve2=4, value=0._dp )
! create system matrix
 call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )
 call create_sysmatrix_data ( sysmatrix )
! build (assemble) matrix and vector from elements
 call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
   coefficients=coefficients )
 call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )
 call solve_system_ma57 ( sysmatrix, rhsd, sol )
 print *, sol%u(subscript%s)
 write(filename, '(a,f6.2,a)') 'sol', alpha, '.vtk'
 call write_scalar_vtk(mesh, problem, filename, dataname = 'u',sysvector=sol)
! delete all data including all allocated memory
 call delete ( problem )
 call delete ( input_probdef )
 call delete ( mesh )
 call delete ( sol, rhsd )
 call delete ( sysmatrix )
 call delete ( coefficients )
contains
function func ( nr, x )
 integer, intent(in) :: nr
 real(dp), intent(in), dimension(:) :: x
 real(dp) :: func
 real(dp) :: mu_x1 = 0.25, mu_y1 = 0.75
 real(dp) :: mu_x2 = 0.75, mu_y2 = 0.25
 real(dp) :: sigma = 0.1
 select case(nr)
   case(1)
     func = beta*(exp(-((x(1) - mu_x1)**2 / (2.0 * sigma**2) + &
                       (x(2) - mu_y1)**2 / (2.0 * sigma**2))))
   case default
     write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
     stop
 end select
end function func
end program poisson1