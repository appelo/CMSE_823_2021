module type_defs
  integer, parameter:: sp = kind(1.0),&
    dp = selected_real_kind(2*precision(1.0_sp)),&
    qp = selected_real_kind(2*precision(1.0_dp))
end module type_defs

module problem_setup
  use type_defs
  implicit none
  integer,  parameter :: Nx = 20
  logical, parameter :: dirichlet_bc = .true.
  logical, parameter :: use_direct = .false. 
  real(dp), parameter :: TOL = 1.0e-12_dp
end module problem_setup

module arrs
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: u,b,x
end module arrs

module afuns
  
contains
  subroutine apply_1D_laplacian_D(au,u,n)
    use type_defs
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out) :: au(n)
    real(dp), intent(in)  ::  u(n)
    integer :: i
    Au(1) = (u(2) - 2.0_dp*u(1)         )
    Au(n) = (     - 2.0_dp*u(n) + u(n-1))
    do i = 2,n-1
     Au(i)= (u(i+1) - 2.0_dp*u(i) + u(i-1))
    end do
  end subroutine apply_1D_laplacian_D

end module afuns

module iterative_solver_D
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: r_isd,x_isd,Ar_isd
  real(dp) :: TOL_ISD
  integer :: iterations
contains
  subroutine set_tol_isd(tol)
    implicit none
    real(dp) :: tol
    TOL_ISD = tol
  end subroutine set_tol_isd

  subroutine get_iter_isd(iter)
    implicit none
    integer :: iter
    iter = iterations
  end subroutine get_iter_isd
  
  subroutine allocate_isd(n)
    implicit none
    integer, intent(in) :: n
    allocate(r_isd(n),x_isd(n),Ar_isd(n))
    iterations  = 0
  end subroutine allocate_isd

  subroutine deallocate_isd
    implicit none
    deallocate(r_isd,x_isd,Ar_isd)
  end subroutine deallocate_isd
  
  subroutine steep_descent_d(x,b,n)
    use type_defs
    use afuns
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout)  :: x(n)
    real(dp), intent(inout)  :: b(n)
    real(dp) :: res_norm2,res_norm20,alpha
    integer :: iter    
    x_isd = 0.0_dp
    r_isd = b
    res_norm20 = sum(r_isd**2)
    res_norm2 = res_norm20
    iter = 0
    do while ((res_norm2/res_norm20 .gt. TOL_ISD**2) &
         .and. (iter .lt. 1000000)) 
       call apply_1D_laplacian_D(Ar_isd,r_isd,n)
       alpha = res_norm2 / sum(r_isd*Ar_isd)
       x_isd = x_isd + alpha*r_isd
       r_isd = r_isd - alpha*Ar_isd       
       res_norm2 = sum(r_isd**2)
       iter = iter + 1
       ! write(*,*) iter, sqrt(res_norm2) 
    end do
    x = x_isd
    iterations = iter
  end subroutine steep_descent_d
end module iterative_solver_D

program ins
  use type_defs
  use problem_setup
  use arrs
  use afuns
  use iterative_solver_D
  implicit none
  ! This program solves u_xx = f 
  ! on the domain [x] \in [0,1] with Dirichlet BC
  ! hx = 1/Nx
  real(dp) :: hx
  integer :: i,n_iter,N_sys,info  
  real(dp), allocatable, dimension(:,:) :: A
  real(dp), allocatable, dimension(:) :: action_of_A,u_inner
  integer, allocatable, dimension(:) ::  ipiv
  ! Set up the grid
  hx = 1.0_dp/real(Nx,dp)
  allocate(x(0:nx))
  allocate(u(0:nx))
  do i = 0,nx
   x(i) = real(i,dp)*hx
  end do

  if(dirichlet_bc) then
     ! For Dirichlet BC we don't solve for the boundary conditions. 
     N_sys = nx-1 ! system size
  else
  end if
  
  allocate(b(N_sys))
  
  if (use_direct) then
     allocate(A(N_sys,N_sys),&
          action_of_A(N_sys),&
          u_inner(N_sys),&
          ipiv(N_sys))
     if(dirichlet_bc) then
        do i = 1,N_sys
           u_inner = 0.0d0
           u_inner(i) = 1.0d0
           call apply_1D_laplacian_D(action_of_A,u_inner,N_sys)
           A(:,i) = action_of_A
        end do
     else
     end if
  end if

  ! Set up the problem.
  ! We cook up a problem with a known solution.
  ! Say that u(x) = exp(-x), then we must have that
  ! f = u_xx  = exp(-x),
  ! u(x = 0) = exp(0), u(x = 1) = exp(-1)
  !
  ! We move the hx^2 from the denominator over to the right hand side of the
  ! system of equations.

  if (dirichlet_bc) then
     do i = 1,nx-1
        b(i) = hx*hx*exp(-x(i))
     end do
     ! We must also account for the boundary conditions
     b(1) = b(1) - exp(-x(0))
     b(nx-1) = b(nx-1) - exp(-x(nx))
     ! Here we either solve with a direct method or with steepest descent
     if (use_direct) then
        CALL DGETRF(N_sys,N_sys,A,N_sys,ipiv,INFO)
        CALL DGETRS('N',N_sys,1,A,N_sys,IPIV,b,N_sys,INFO)
        u(1:nx-1) = b 
        n_iter = 0
     else
        call allocate_isd(N_sys)
        call set_tol_isd(tol)
        call steep_descent_d(u(1:nx),b,N_sys)
        call get_iter_isd(n_iter)
        call deallocate_isd
     end if
     u(0) = exp(-x(0))
     u(nx) = exp(-x(nx))
     write(*,*) tol, n_iter, maxval(abs(u - exp(-x)))
  else
  end if
    
end program ins
