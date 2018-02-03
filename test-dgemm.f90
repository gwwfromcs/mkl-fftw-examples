program test_dgemm

  implicit none

  integer, parameter :: m = 10
  integer, parameter :: k = 5
  integer, parameter :: n = 7 
  integer :: i , j

  real(8) :: a(m,k), b(k, n), c(m, n), d(k,k), e(m,m), alpha, beta

  alpha = 1.0
  beta = 0.0

  write(*,'(" a = ")')
  do i = 1, m
     do j = 1, k
        a(i,j) = real(i)/real(j)
        write(*,'(f12.3)',advance='no') a(i,j)
     enddo
     write(*,'()')
  enddo

  write(*,'(" b = ")')
  do i = 1, k
     do j = 1, n
        b(i,j) = real(i)/real(j)
        write(*,'(f12.3)',advance='no') b(i,j)
     enddo
     write(*,'()')
  enddo

  call dgemm('n','n',m,n,k,alpha,a,m,b,k,beta,c,m)

  write(*,'(" c = a * b")')
  do i = 1, m 
     do j = 1, n
        write(*,'(f12.3)',advance='no') c(i,j)
     enddo
     write(*,'()')
  enddo

  call dgemm('t','n',k,k,m,alpha,a,m,a,m,beta,d,k)

  write(*,'(" d = a^T * a ")')
  do i = 1, k
     do j = 1, k
        write(*,'(f12.3)',advance='no') d(i,j)
     enddo
     write(*,'()')
  enddo
       
  call dgemm('n','t',m,m,k,alpha,a,m,a,m,beta,e,m)

  write(*,'(" e = a * a^T")')
  do i = 1, m
     do j = 1, m
        write(*,'(f12.3)',advance='no') e(i,j)
     enddo
     write(*,'()')
  enddo

end program test_dgemm
