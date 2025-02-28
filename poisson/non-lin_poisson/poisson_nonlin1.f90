! non-linear Poisson problem on a unit square with Dirichlet 
! boundary conditions
! diffusion coefficient assumed to be 
!     D = 1 + delta*u^2 + gamma*|\nabla u|
! where delta and gamma are constants

program poisson_nonlin1

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use hsl_ma41_m
  use poisson_elements_m
  use poisson_functions_m
  use io_utils_m

  implicit none


! constants

  integer :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx = 20,            & ! number of elements in x
    ny = 20,            & ! number of elements in y
    nPicard = 10,       & ! number of Picard iterations before Newton-Raphson
    itermax = 100         ! maximum interations

  real(dp) :: &
    delta = 3._dp,    &   ! delta coefficient
    gamma = 0.4_dp,   &   ! gamma coefficient
    eps = 1.e-9_dp        ! accuracy in Picard iteration


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(sysvector_t), target :: soln
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(subscript_t) :: resu

  integer :: iter
  real(dp) :: udiff
  logical :: newton_raphson
  character(len=999) :: outputfile = 'sol.vtk'

  ! namelist for input of variables; read from standard input

  namelist /comppar/ outputfile, delta, gamma, nx, ny, eps, nPicard, itermax

  read ( unit=*, nml=comppar )

  !  pass to global variables  

  ldelta = delta
  lgamma = gamma

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:11) = (/ gauss, gaussb /)

  coefficients%r = 0 

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

  call create_subscript (mesh, problem, resu, &
    essentialpart=.false.)

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln )
  call create_sysvector ( problem, rhsd )

! create oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => soln

! initial guess
! NOTE: do not use a uniform field for nPicard = 0, since
! in the Newton-Raphson iteration there is a 1 / |\nabla u|

  soln%u = 0._dp

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, soln, curve1=3, value=1._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.false. )

  call create_sysmatrix_data ( sysmatrix )

  sol%u = 0._dp ! initialize

  ! iteration

  iter = 0

  iterate: do 

    iter = iter + 1

    newton_raphson = iter > nPicard

  ! build (assemble) matrix and vector from elements 

    if ( .not. newton_raphson  ) then
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=poisson_elem_picard, coefficients=coefficients, &
        oldvectors=oldvectors )
      call fill_sysvector ( mesh, problem, sol, curve1=3, value=1._dp )

    else
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=poisson_elem_newton_raphson, coefficients=coefficients, &
        oldvectors=oldvectors )
      call fill_sysvector ( mesh, problem, sol, curve1=3, value=0._dp )
    end if

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma41 ( sysmatrix, rhsd, sol )

    if ( newton_raphson ) then
        sol%u = sol%u + soln%u
        print *,'residual = ',maxval(abs(rhsd%u(resu%s)))
    end if

    udiff = maxval(abs(sol%u-soln%u))

    print *,'iter = ',iter, 'udiff = ',udiff

    if ( udiff < eps ) exit iterate

    if ( iter >= itermax ) then
      write(*,'(a,i0)') ' too many iterations: ', itermax
      exit iterate
      ! stop
    end if

    call copy ( sol, soln )

  end do iterate

! plot solution

  call write_scalar_vtk (mesh, problem, trim(adjustl(outputfile)), &
   'u', sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson_nonlin1
