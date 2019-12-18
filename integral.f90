program quadr
    use func
    use integrator
    implicit none
    
    real(8) :: a, b
    integer(4) :: n, i
    character*30 :: math_function
    
    open(2, file='func.f90')
    do i = 1, 10
        read(2,'(A)') math_function
    enddo
    close(2)
    
    write(*,*) 'Please, enter left integration limit'
    read(*,*) a
    write(*,*) 'Please, enter right integration limit'
    read(*,*) b
    write(*,*) 'Please, enter number of intervals'
    read(*,*) n
    write(*,*) 'Integrated function is', trim(math_function(12:30))
    write(*,*) 'integral value is', gauss_integral(n, a, b, f)
    
end program quadr

