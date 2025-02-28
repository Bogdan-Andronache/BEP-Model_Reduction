module poisson_functions_m

  use tfem_m
  use poisson_elements_m
  use math_defs_m

  implicit none

  real(dp) :: ldelta, lgamma

contains

  function alpha_func ( u, gradu_magn )

    real(dp), intent(in) :: u, gradu_magn
    real(dp) :: alpha_func

    alpha_func = 1 + ldelta*u**2._dp + lgamma*gradu_magn
  
  end function alpha_func

  function dalpha_du_func ( u, gradu_magn )

    real(dp), intent(in) :: u, gradu_magn
    real(dp) :: dalpha_du_func

    dalpha_du_func = 2._dp*ldelta*u
  
  end function dalpha_du_func
  
  function dalpha_dmagngradu_func ( u, gradu_magn )

    real(dp), intent(in) :: u, gradu_magn
    real(dp) :: dalpha_dmagngradu_func

    dalpha_dmagngradu_func = lgamma
  
  end function dalpha_dmagngradu_func
 
  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    func = 0

    select case(nr)
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func 

  function vfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    vfunc = 0
                                                                                
    select case(nr)
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select
                                                                                
  end function vfunc

! Internal element routine for the Poisson Equation
! ( with alpha given by a function).

  subroutine poisson_elem_picard ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )
  
    use poisson_globals_m
  
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec
  
  
    integer :: i, j, ip, funcnr
    ! real(dp) :: alpha

    real(dp), dimension(:), allocatable, save :: alphag, u_iter, ug
    real(dp), dimension(:,:), allocatable, save :: gradu
  
    if ( first ) then
  
  !     first element in this group
  
  !     set globals
  
      call set_globals_poisson ( mesh, coefficients, elgrp )
  
      allocate ( wg(ninti), fg(ninti), detF(ninti), work(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( u_iter(ndf), gradu(ninti,ndim), alphag(ninti), ug(ninti) )

  !     set Gauss integration and shape function
  
      call set_Gauss_integration ( gauss, xig, wg )
  
      call set_shape_function ( shapefunc, xig, phi, dphi )
  
    end if
  
    call get_coordinates ( mesh, elgrp, elem, x )
  
    call isoparametric_deformation ( x, dphi, F, Finv, detF )
  
    call isoparametric_coordinates ( x, phi, xg )
  
    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if
  
    call shape_derivative ( dphi, Finv, dphidx )
  
   ! get the non-linear diffusion coefficient

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, &
    elem, u_iter, layer=layer )

    ug = matmul ( phi, u_iter )
  
    do i = 1, ndim
      gradu(:,i) = matmul ( dphidx(:,:,i), u_iter )
    end do

    do ip = 1, ninti
      alphag(ip) = alpha_func(ug(ip),sqrt(sum(gradu(ip,:)**2)))
    end do
  
    if ( vector ) then
  
      funcnr = coefficients%i(12)
  
      if ( funcnr > 0 ) then
  
        do ip = 1, ninti
          fg(ip) = coefficients%func ( funcnr, xg(ip,:) )
        end do
  
        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do
  
      else
  
        elemvec = 0
  
      end if
  
    end if
  
    if ( matrix ) then
  
      ! alpha = coefficients%r(1)
  
      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( alphag(ip) * dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do
  
    end if 
  
    if ( last ) then
  
  !     last element in this group
   
      deallocate ( Finv, dphidx )
      deallocate ( dphi, F )
      deallocate ( xg )
      deallocate ( xig, phi, x, work )
      deallocate ( wg, fg, detF )
      deallocate ( u_iter, gradu, alphag, ug )

    end if
  
  end subroutine poisson_elem_picard

! Internal element routine for the Poisson Equation
! ( with alpha given by a function).
  
  subroutine poisson_elem_newton_raphson ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )
  
    use poisson_globals_m
  
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec
  
  
    integer :: i, j, ip, funcnr
    ! real(dp) :: alpha

    real(dp), dimension(:), allocatable, save :: alphag, u_iter, ug, &
      dalpha_du_funcg, dalpha_dmagngradu_funcg, gradu_magn
    real(dp), dimension(:,:), allocatable, save :: gradu
  
    if ( first ) then
  
  !     first element in this group
  
  !     set globals
  
      call set_globals_poisson ( mesh, coefficients, elgrp )
  
      allocate ( wg(ninti), fg(ninti), detF(ninti), work(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( u_iter(ndf), gradu(ninti,ndim), alphag(ninti), ug(ninti) )
      allocate ( dalpha_du_funcg(ninti), dalpha_dmagngradu_funcg(ninti) )
      allocate ( gradu_magn(ninti) )

  !     set Gauss integration and shape function
  
      call set_Gauss_integration ( gauss, xig, wg )
  
      call set_shape_function ( shapefunc, xig, phi, dphi )
  
    end if
  
    call get_coordinates ( mesh, elgrp, elem, x )
  
    call isoparametric_deformation ( x, dphi, F, Finv, detF )
  
    call isoparametric_coordinates ( x, phi, xg )
  
    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if
  
    call shape_derivative ( dphi, Finv, dphidx )
  
   ! get the non-linear diffusion coefficient

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, &
    elem, u_iter, layer=layer )

    ug = matmul ( phi, u_iter )
  
    do i = 1, ndim
      gradu(:,i) = matmul ( dphidx(:,:,i), u_iter )
    end do

    gradu_magn = norm2 ( gradu, dim=2)

    do ip = 1, ninti
      alphag(ip) = alpha_func(ug(ip),gradu_magn(ip))
      dalpha_du_funcg(ip) = dalpha_du_func(ug(ip),gradu_magn(ip))
      dalpha_dmagngradu_funcg(ip) = dalpha_dmagngradu_func(ug(ip),gradu_magn(ip))
    end do
  
    if ( vector ) then
  
      funcnr = coefficients%i(12)
  
      if ( funcnr > 0 ) then
  
        do ip = 1, ninti
          fg(ip) = coefficients%func ( funcnr, xg(ip,:) )
        end do
  
        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do
  
      else
  
        elemvec = 0
  
      end if
  
      do i = 1, ndf
        work = alphag * sum ( dphidx(:,i,:) * gradu, dim=2 )
        elemvec(i) = elemvec(i) - sum ( work * detF * wg )
      end do

    end if
  
    if ( matrix ) then
  
      ! alpha = coefficients%r(1)
  
      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( alphag(ip) * dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do
  
      do i = 1, ndf
        do j = 1, ndf
            work = sum ( dphidx(:,i,:) * gradu, dim=2 )
            elemmat(i,j) = elemmat(i,j) + &
                sum ( work * dalpha_du_funcg * phi(:,j) * detF * wg )
        end do
      end do

      do i = 1, ndf
        do j = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,i,:) * gradu(ip,:) ) * &
                       sum ( dalpha_dmagngradu_funcg(ip) * &
                            gradu(ip,:) / gradu_magn(ip) * dphidx(ip,j,:)  )
          end do
          elemmat(i,j) = elemmat(i,j) + sum ( work * detF * wg )
          ! elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if 
  
    if ( last ) then
  
  !     last element in this group
   
      deallocate ( Finv, dphidx )
      deallocate ( dphi, F )
      deallocate ( xg )
      deallocate ( xig, phi, x, work )
      deallocate ( wg, fg, detF )
      deallocate ( u_iter, gradu, alphag, ug )
      deallocate ( dalpha_du_funcg, dalpha_dmagngradu_funcg )
      deallocate ( gradu_magn )
    end if
  
  end subroutine poisson_elem_newton_raphson
                                                                                
end module poisson_functions_m
 


