module LinearAlgebraInterfaces
    
    use MathAndProgConstants
    use BaseMeansToolbox
    use LinearAlgebraLowLevel
    
    contains
    
    ! About function, open for F2PY wrapping
    function linearalgebrainterfaces_about() result(i)
            integer :: i
            i = 0
            write(*,*) 'LinearAlgebraInterfaces: Version 1.0.'
            write(*,*) 'Author: Maksim Shirobokov.'
            write(*,*) 'Date: 19.01.2022.'
    end function linearalgebrainterfaces_about
    
    function linsolve(A,b) result(x)
    
        real(dp), dimension(:,:) :: A
        real(dp), dimension(:)   :: b
        real(dp), dimension(size(A,1)) :: x
        
        integer, dimension(size(b)) :: ipiv
        integer :: info
        real(dp), dimension(size(b),1) :: bm
        
        bm(:,1) = b

        ! Before
        !call getrf(A,ipiv)
        !call getrs(A,ipiv,bm)

        ! After
        call dgetrf(size(A,1), size(A,2), A, size(A,1), ipiv, info)
        call dgetrs('N', size(A,2), 1, A, size(A,1), ipiv, bm, size(b), info)
        
        x = bm(:,1)
    
    end function linsolve
    
    function dotmm(M1,M2) result(M)
    
        real(dp), dimension(:,:) :: M1
        real(dp), dimension(:,:) :: M2
        real(dp), dimension(size(M1,1),size(M2,2)) :: M
        
        M = 0.0D0
        call dgemm('N', 'N', size(M1,1), size(M2,2), size(M1,2), 1.0D0, M1, size(M1,1), M2, size(M2,1), 1.0D0, M, size(M,1))
    
    end function dotmm
    function dotvm(v,M) result(w)
    
        real(dp), dimension(:) :: v
        real(dp), dimension(:,:) :: M
        real(dp), dimension(size(M,2)) :: w
        
        w = reshape( dotmm (reshape(v,(/1,size(v,1)/)) , M) , (/size(M,2)/) )
    
    end function dotvm
    function dotmv(M,v) result(w)
    
        real(dp), dimension(:) :: v
        real(dp), dimension(:,:) :: M
        real(dp), dimension(size(M,1)) :: w
        
        w = reshape( dotmm( M , reshape(v,(/size(v,1),1/))) , (/size(M,1)/) )
    
    end function dotmv
    function doto(v,w) result(M)
    
        real(dp), dimension(:) :: v
        real(dp), dimension(:) :: w
        real(dp), dimension(size(v,1),size(w,1)) :: M
        
        M = 0
        call dgemm('N', 'T', size(v,1), size(w,1), 1, 1.0D0, v, size(v,1), w, size(w,1), 1.0D0, M, size(M,1))
        
    end function doto
    
end module LinearAlgebraInterfaces