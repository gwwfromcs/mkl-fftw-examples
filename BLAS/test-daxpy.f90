! On utexas-ICES computers:
!  module load mkl
!  ifort -mkl test-daxpy.f90 -o daxpytest
program testdaxpy

  integer, parameter :: dp = kind(1.0d0)
  real(dp), dimension(4)  :: x
  real(dp), dimension(12) :: y
  integer :: i

  x = (/ 1.0d0, 1.5d0, 2.2d0, -8.0d0 /)
  y = 0d0

  call daxpy(4,1.0d0,x(1),1,y(1),3)  ! Note: that you must use 1.0d0 here
  call daxpy(4,1.5d0,x(1),1,y(2),3)
  call daxpy(4,2.0d0,x(1),1,y(3),3)

  write(*,'((f7.3))') (x(i), i = 1,4)
  write(*,'((f7.3))') (y(i), i = 1,12)

end program
