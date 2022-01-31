module mymodule

    use VecMatInterfaces

contains

function dotmm() result(x)

    implicit none

    real(8) :: A(2,2)
    real(8) :: M2(2,2)
    real(8) :: M(2,2)
    real(8) :: b(2), bm(2,1), x(2)
    integer :: ipiv(2)
    integer :: info

    A(1,:) = (/ 1.0D0, 2.0D0/)
    A(2,:) = (/ 3.0D0, 4.0D0/)
    b = (/ 1.0D0, 3.0D0 /)
    bm(:,1) = b

    call dgetrf(size(A,1), size(A,2), A, size(A,1), ipiv, info)
    call dgetrs('N', size(A,2), 1, A, size(A,1), ipiv, bm, size(b), info)

    x = bm(:,1)

end function dotmm

end module mymodule