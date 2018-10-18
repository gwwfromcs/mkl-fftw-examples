! On utexas-ICES computers:
!   module load mkl
!   ifort -mkl test-dgemm.f90 -o dgemmtest
! In this example, we call dgemm to calculate matrix multiplications from 
!   submatrices of x(:,:) and y(:,:)
program testdgemm

integer, parameter :: dp = kind(1.0d0)
real(dp) :: x(6,3), y(3,2), z(2,2)
integer :: i, j

write(*, *) " x = [ "
do i = 1, 6
  do j = 1, 3
    x(i,j) = i-j
  enddo
  write(*, '(3f8.3)') (x(i,j),j=1,3)
enddo
write(*, *) " ] "

write(*, *) " y = [ "
do i = 1, 3
  do j = 1, 2
    y(i,j) = i+j
  enddo
  write(*, '(2f8.3)') (y(i,j),j=1,2)
enddo
write(*, *) " ] "

call dgemm('n', 'n', 2, 2, 3, 1.0d0, x(1,1), 6, y(1,1), 3, 0.0d0, z(1,1), 2)

write(*, *) " z = x(1:2,1:3)*y(1:3,1:2) = [ "
do i = 1, 2
  write(*, '(2f8.3)') (z(i,j),j=1,2)
enddo
write(*, *) " ] "

call dgemm('n', 'n', 2, 2, 3, 1.0d0, x(3,1), 6, y(1,1), 3, 0.0d0, z(1,1), 2)

write(*, *) " z = x(3:4,1:3)*y(1:3,1:2) = [ "
do i = 1, 2
  write(*, '(2f8.3)') (z(i,j),j=1,2)
enddo
write(*, *) " ] "

call dgemm('n', 'n', 2, 2, 3, 1.0d0, x(5,1), 6, y(1,1), 3, 0.0d0, z(1,1), 2)

write(*, *) " z = x(5:6,1:3)*y(1:3,1:2) = [ "
do i = 1, 2
  write(*, '(2f8.3)') (z(i,j),j=1,2)
enddo
write(*, *) " ] "

call dgemm('n', 'n', 2, 2, 2, 1.0d0, x(1,2), 6, y(2,1), 3, 0.0d0, z(1,1), 2) 

write(*, *) " z = x(1:2,2:3)*y(2:3,1:2) = [ "
do i = 1, 2
  write(*, '(2f8.3)') (z(i,j),j=1,2)
enddo
write(*, *) " ] "

call dgemm('t', 'n', 2, 2, 2, 1.0d0, x(1,2), 6, y(2,1), 3, 0.0d0, z(1,1), 2) 

write(*, *) " z = x(1:2,2:3)*y(2:3,1:2)^T = [ "
do i = 1, 2
  write(*, '(2f8.3)') (z(i,j),j=1,2)
enddo
write(*, *) " ] "

end program
