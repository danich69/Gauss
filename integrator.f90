module integrator

use func

contains

    function gauss_integral(n, a, b, f) result(integral)            ! integration function
        
            interface

                function f(x)
                    real(8) :: x, f
                end function f
                
            end interface
                
            real(8), dimension(n) :: coeff, t
            real(8) :: a, b, integral, hfsum, hfsub
            integer(4) :: n, i
            character*10 :: file_name
            
            write (file_name,"('quad',i0,'.dat')") n
            open(1, file=trim(file_name))
            do i = 1, n
                read(1,*) coeff(i), t(i)
            enddo    
            close(1)
            
            integral = 0
            hfsub = (b-a)/2
            hfsum = (a+b)/2
            
            do i = 1, n
            
                integral = integral + coeff(i)*f(hfsub*t(i)+hfsum)
            
            enddo
            
            integral = integral*hfsub
        
    end function gauss_integral

end module integrator
