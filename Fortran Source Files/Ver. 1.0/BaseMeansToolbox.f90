module BaseMeansToolbox

    use MathAndProgConstants

    contains

        ! Maximum norm
        function norminf(y) result(output)
            real(dp), dimension(:) :: y
            real(dp)               :: output
            output = maxval(abs(y))
        end function

        ! Second norm
        function norm(y) result(output)
            real(dp), dimension(:) :: y
            real(dp)               :: output
            output = sqrt(sum(y*y))
        end function

        ! Zeros matrix
        function zeros(m,n) result(output)
            integer(8) :: m, n
            real(dp), dimension(m,n) :: output
            output = 0.0D0
        end function

        ! Ones matrix
        function ones(m,n) result(output)
            integer(8) :: m, n
            real(dp), dimension(m,n) :: output
            output = 1.0D0
        end function

        ! Identity matrix
        function eye(n) result (output)
            integer :: n, i
            real(dp), dimension(n,n) :: output
            output = 0.0D0
            forall (i=1:n) output(i,i) = 1.0D0
        end function

        ! Cross product
        function cross(a,b) result(c)

            real(dp) :: a(3), b(3), c(3)

            c(1) = a(2)*b(3) - a(3)*b(2)
            c(2) = a(3)*b(1) - a(1)*b(3)
            c(3) = a(1)*b(2) - a(2)*b(1)

        end function

        ! MATLAB diff function for vectors
        function diff(v) result(w)

            real(dp) :: v(:)
            real(dp) :: w(size(v)-1)
            integer :: k, n

            n = size(v)
            do k = 2,n
                w(k-1) = v(k) - v(k-1)
            end do

        end function

        ! Concatenate matrices by columns
        function ConCol(A,B) result(C)
            real(dp), dimension(:,:)  :: A, B
            real(dp), dimension(size(A,1),size(A,2)+size(B,2)) :: C
            C(:,1:size(A,2)) = A
            C(:,size(A,2)+1:size(A,2)+size(B,2)) = B
        end function

        ! Concatenate matrices by rows
        function ConRow(A,B) result(C)
            real(dp) :: A(:,:), B(:,:)
            real(dp) :: C(size(A,1)+size(B,1),size(A,2))
            C(1:size(A,1),:) = A
            C(size(A,1)+1:size(A,1)+size(B,1),:) = B
        end function

        ! MATLAB colon operator and linspace function
        function ColonOperator(i1,i2,step) result(v)
            integer i1, i2, step
            integer i
            integer :: v(i2-i1+1)
            do i = 1,i2-i1+1
                v(i) = i1 + step * (i - 1)
            end do
        end function
        function linspace(a,b,n) result(v)

            integer :: n
            real(dp) :: a, b, rho
            real(dp) :: v(n)
            integer :: k

            rho = (b - a)/(n-1)
            do k = 1,n
                v(k) = a + (k-1) * rho
            end do

        end function

        !!! Python returns -1 if true and 0 if false.
        function double2logical(d) result(l)
            implicit none
            real(dp) :: d
            logical :: l
            if (int(d+1.0D-10) .eq. 1) then
                l = .true.
            else
                l = .false.
            end if
        end function

        function basemeanstoolbox_about() result(i)
            integer :: i
            i = 0
            write(*,*) 'This is BaseMeansToolbox.'
        end function basemeanstoolbox_about

end module
