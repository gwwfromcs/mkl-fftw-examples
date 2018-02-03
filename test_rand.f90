program test_rand

 use ifport
 implicit none
 integer, parameter :: rseed = 86456
 real :: a

 call srand(rseed)
 a = rand(0) 
 !print *, rand(), rand(), rand(), rand()
 !print *, rand(seed), rand(), rand(), rand()
 print *, a

end program test_rand
