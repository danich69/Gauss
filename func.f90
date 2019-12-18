module func
    implicit none
    
    contains
    
        function f(x)
        
            real(8) :: x, f
            
            f = sqrt(x)
        
        end function f
    
end module func
