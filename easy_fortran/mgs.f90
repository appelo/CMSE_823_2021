module type_defs
  integer, parameter:: sp = kind(1.0),&
       dp = selected_real_kind(2*precision(1.0_sp)),&
       qp = selected_real_kind(4*precision(1.0_sp))
end module type_defs

program mgs
  use type_defs
  implicit none
  integer :: n,i,j,k
  real(kind = dp), allocatable, dimension(:,:) :: A,q,r,QTQ,V
  real(dp) :: length_v, time1, time2


  do n = 5,500,10

     allocate(A(n,n),Q(n,n),R(1:n,1:n),v(n,n),QTQ(n,n))

     A = 0.0_dp
     Q = 0.0_dp
     R = 0.0_dp
     ! Create the typical FD matrix
     do i = 1,n,1
        A(i,i) = -2.0_dp
     end do
     do i = 1,n-1
        A(i,i+1) = 1.0_dp
        A(i+1,i) = 1.0_dp
     end do

     call cpu_time(time1)

     V = A
     do i = 1,n
        R(i,i) = sqrt(sum(V(:,i)**2))
        Q(:,i) = V(:,i) / R(i,i)
        do j = i+1,n
           R(i,j) = sum(Q(:,i)*V(:,j))
           V(:,j) = V(:,j) - R(i,j)*Q(:,i)
        end do
     end do
     call cpu_time(time2)

     QTQ = matmul(transpose(Q),Q)
     do k = 1,n
        QTQ(k,k) = QTQ(k,k) - 1.0_dp
     end do
     write(*,*) n, time2-time1, maxval(abs(A-matmul(Q,R))), maxval(abs(QTQ))
     !  call print_mat(A,1,n,1,n,'A.txt')
     !  call print_mat(Q,1,n,1,n,'Q.txt')
     !  call print_mat(R,1,n,1,n,'R.txt')

     deallocate(A,Q,R,v,qtq)

  end do
end program mgs

subroutine print_mat(u,nx1,nx2,ny1,ny2,str)
  use type_defs
  implicit none
  integer, intent(in) :: nx1,nx2,ny1,ny2
  real(dp), intent(in) :: u(nx1:nx2,ny1:ny2)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do i=nx1,nx2,1
     do j=ny1,ny2,1
        if(abs(u(i,j)) .lt. 1e-40) then
           write(2,fmt='(E24.16)',advance='no') 0.d0
        else
           write(2,fmt='(E24.16)',advance='no') u(i,j)
        end if
     end do
     write(2,'()')
  end do
  close(2)
end subroutine print_mat
