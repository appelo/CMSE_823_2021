module type_defs
  integer, parameter:: sp = kind(1.0),&
       dp = selected_real_kind(2*precision(1.0_sp)),&
       qp = selected_real_kind(4*precision(1.0_sp))
end module type_defs

program cgs
  use type_defs
  implicit none
  integer :: n,i,j,k
  real(kind = dp), allocatable, dimension(:,:) :: A,q,r,QTQ
  real(kind = dp), allocatable, dimension(:) :: funny_vector,v
  real(dp) :: length_v, time1, time2


  do n = 5,500,10

     allocate(A(n,n),Q(n,n),R(1:n,1:n),funny_vector(-23:8),v(n),QTQ(n,n))

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

     do j = 1,n
        ! v_j = a_j
        do k = 1,n
           v(k) = A(k,j)
        end do
        ! v = A(:,j)
        ! v = A(1:n,j)
        do i = 1,j-1
           ! R_ij = dot(q_i,A_j)
           R(i,j) = 0.0_dp
           do k = 1,n
              R(i,j) = R(i,j) + Q(k,i)*A(k,j)
           end do

           ! Subtract old cols. v_j = v_j - r_ij*qi
           do k = 1,n
              v(k) = v(k) - R(i,j)*Q(k,i)
           end do
           ! v = v - R(i,j)*Q(:,i)
        end do
        ! compute length of v
        length_v = 0.0_dp
        do k = 1,n
           length_v = length_v + v(k)**2
        end do
        length_v = sqrt(length_v)
        ! length_v = sqrt(sum(v**2))
        ! scale v_j to get Q_j
        R(i,j) = length_v
        do k = 1,n
           Q(k,j) = v(k)/length_v
        end do
        ! Q(:,j) = v/length_v
     end do

     call cpu_time(time2)

     QTQ = matmul(transpose(Q),Q)
     do k = 1,n
        QTQ(k,k) = QTQ(k,k) -1.0_dp
     end do
     write(*,*) n, time2-time1, maxval(abs(A-matmul(Q,R))), maxval(abs(QTQ))
     !  call print_mat(A,1,n,1,n,'A.txt')
     !  call print_mat(Q,1,n,1,n,'Q.txt')
     !  call print_mat(R,1,n,1,n,'R.txt')

     deallocate(A,Q,R,funny_vector,v,qtq)

  end do
end program cgs

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
