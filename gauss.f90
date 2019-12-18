module gauss
    use linear_system
    
    implicit none
  
    contains
        
        function quadrature_coefficents(roots) result(A)                ! function which gives roots' weights for integration
        
            real(16), dimension(:) :: roots
            real(16), dimension(size(roots),size(roots)) :: roots_matrix
            real(16), dimension(size(roots)) :: A, coefficent
            integer(4) :: n, k
            
            n = size(roots)
            coefficent = 0
            
            do k = 0, n-1, 2        
                coefficent(k+1) = 2.0_16/(k+1)        
            enddo
            
            do k = 0, n-1        
                roots_matrix(k+1,:) = roots**k        
            enddo
            
            call Lin_sys_solver(roots_matrix,coefficent,A,n,'ModGauss')
        
        end function quadrature_coefficents
        
        function legendre_roots(n) result(t)                            ! function which calculates Legendre polynom roots with Bernulli iterative method
        
            real(16), dimension(n) :: t
            real(16), dimension(n+1) :: poly
            real(16), allocatable :: new_poly(:)
            real(16) :: current_root
            integer(4) :: n, k, i
            
            poly = legendre_poly(n)
            t = 0
            
            if(mod(n,2) == 0) then
            
                allocate(new_poly(n/2))
                new_poly = squeeze(poly)
                
            else
            
                allocate(new_poly((n+1)/2))
                poly = cshift(poly,1)
                new_poly = squeeze(poly)
                
            endif
            
            k = size(new_poly)
            
            do i = 1, k-1
            
                current_root = bernulli(new_poly(1:k-i+1)) 
                t(i) = -sqrt(current_root)
                t(n-i+1) = sqrt(current_root)
            
                new_poly(1:k-i) = horner(new_poly(1:k-i+1), [-current_root,1.0_16])
                new_poly(k-i+1) = 0
            
            enddo
            
            deallocate(new_poly)
            
        end function legendre_roots
        
        function legendre_poly(n) result(T)                             ! function which calculates Legendre polynom coeficents
        
            real(16), dimension(0:n, 0:n) :: P
            real(16), dimension(0:n) :: T
            integer :: n, k
            
            P = 0
            P(0,0) = 1
            P(1,1) = 1
            
            do k = 2, n        
                P(k,:) = real(2*k-1.0_16)/k*cshift(P(k-1,:),-1)-real(k-1)/k*P(k-2,:)        
            enddo
            
            T = P(n,:)
        
        end function legendre_poly
        
        function squeeze(v)                                             ! this one replaces x^2 with t in Legendre polynom
        
            real(16), dimension(:) :: v
            real(16), dimension((size(v)+1)/2) :: squeeze
            integer :: n, i
            
            n = size(v)
            
            squeeze(1) = v(1)
            
            do i = 1, (n+1)/2-1        
                squeeze(i+1) = v(i*2+1)        
            enddo
            
        end function squeeze
        
        function horner(a, b) result(res)                               ! this one implimends polynom division
        
            real(16), dimension(:) :: a, b
            real(16), allocatable :: c(:,:), d(:), res(:)
            integer :: n, k, i
            
            n = size(a)
            k = size(b)
            allocate(c(n-k+1,n), d(n), res(n-k+1))
            
            d = 0
            d(n-k+1:n) = b(1:k)
            c = 0
            
            c(1,:) = a-d*a(n)/d(n)
            res(n-k+1) = a(n)/d(n)
            
            do i = 2, n-k+1
            
                d = cshift(d,1)
                c(i,:) = c(i-1,:)-d*c(i-1,n-i+1)/d(n-i+1)
                res(n-k-i+2) = c(i-1,n-i+1)/d(n-i+1)
            
            enddo
            
        end function horner
        
        function bernulli(a) result(root)                               ! Bernulli method 
        
            real(16), dimension(:) :: a
            real(16), dimension(size(a,dim=1)) :: y
            real(16) :: eps, root
            integer :: n, i, j
            
            n = size(a)
            eps = 1e-8
            
            if (n == 2) then
            
                root = -a(1)/a(2)
                return
                
            else if (n == 1) then
                
                write(*,*) 'Warning! Corrupted data: 0 degree poly at&
                & Bernulli method!'
                
            endif
            
            do while (a(n) == 0)
            
                n = n-1
            
            enddo
            
            y = 0
            call random_number(y)
            do while(abs(y(n)/y(n-1)-y(n-1)/y(n-2)) >= eps)

                y = cshift(y,1)
                y(n) = 0
                y(n) = -dot_product(y(1:n-1),a(1:n-1))/a(n)
                
            enddo
            
            root = y(n)/y(n-1)
            
        end function bernulli
 
end module gauss
