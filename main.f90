program mesh
	use gauss
    use func
    implicit none
    	
!	-----------------------------------------------------
!	|	Program integrates 1D function with        		|
!	|			Gauss quadratures     					|
!	|													|
!	|	Just wirte wanted order after ./prog			|
!	|													|
!	----------------------------------------------------
    real(16), allocatable :: coef(:), roots(:)
    real(16) :: a, b
    integer(4) :: n, i
    character*10 :: file_name
    character*2 :: node_num
    
    call getarg(1, node_num)
    read(node_num,*) n
    
    if(n == 0) then
        write(*,*) 'Incorrect input! n=0'
        read(*,*) n
    endif
    
    write (file_name,"('quad',i0,'.dat')") n
    allocate(coef(n), roots(n))
    
    roots = legendre_roots(n)
    coef = quadrature_coefficents(roots)
    
    open(1, file=trim(file_name))
    do i = 1, n
        write(1,*) coef(i), roots(i)
    enddo    
    close(1)
    
    deallocate(coef, roots)

end program mesh
