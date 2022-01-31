module OdeToolbox
    
    use BaseMeansToolbox
    use MathAndProgConstants
    
    type OdeSol
        real(dp), allocatable :: T(:)
        real(dp), allocatable :: Y(:,:)
    end type

    type OdeOptions
        real(dp) :: NumberOfIntervals = 20000
        real(dp) :: AbsTol = 1.0D-10
        real(dp) :: RelTol = 1.0D-10
        logical :: OnlyFinalState = .true.
        integer :: chunk = 200000
    end type
    
    contains
    
        ! About function, opened for F2PY wrapping
        function odetoolbox_about() result(i)
            integer :: i
            i = 0
            write(*,*) 'OdeToolbox: Version 1.0.'
            write(*,*) 'Author: Maksim Shirobokov.'
            write(*,*) 'Date: 19.01.2022.'
        end function odetoolbox_about
    
        ! The classical Runge-Kutta method, 4th order
        function ode4(fcn,tspan,y0,options) result(sol)
        
        implicit none
        
            interface
                function fcn(t,y) result(dydt)
                    use MathAndProgConstants
                    real(dp),               intent(in)   :: t
                    real(dp), dimension(:), intent(in)   :: y
                    real(dp), dimension(size(y))         :: dydt
                end function fcn
            end interface
            
            real(dp), dimension(2),intent(in) :: tspan
            real(dp), dimension(:),intent(in) :: y0
            real(dp), dimension(size(y0))     :: yend, k1, k2, k3, k4
            real(dp) :: t0, tf, t, h
            integer :: N, k
            
            type(OdeOptions) :: options
            type(OdeSol)     :: sol
            
            N = options%NumberOfIntervals
            
            if (options%OnlyFinalState) then
                allocate(sol%T(2));          sol%T = 0.0D0
                allocate(sol%Y(size(y0),2)); sol%Y = 0.0D0
            else
                allocate(sol%T(N+1));          sol%T = 0.0D0
                allocate(sol%Y(size(y0),N+1)); sol%Y = 0.0D0
            end if
            
            t0 = tspan(1)
            tf = tspan(2)
            
            sol%T(1) = t0
            sol%Y(:,1) = y0
            
            h = (tf-t0)/N
            
            t    = t0
            yend = y0
            do k = 1,N
                   
                k1 = fcn(t,         yend)
                k2 = fcn(t+h/2.0D0, yend + h * k1/2.0D0)
                k3 = fcn(t+h/2.0D0, yend + h * k2/2.0D0)
                k4 = fcn(t+h,       yend + h * k3)
                
                t = t0 + h*k
                yend = yend + h*(k1 + 2.0D0*k2 + 2.0D0*k3 + k4)/6.0D0
                
                if (.NOT. options%OnlyFinalState) then
                    sol%T(k+1) = t
                    sol%Y(:,k+1) = yend
                end if
                
            end do
            
            if (options%OnlyFinalState) then
                sol%T(2) = t
                sol%Y(:,2) = yend
            end if
        
        end function ode4
        
        ! Runge-Kutta method, 8th order
        function ode8(fcn,tspan,y0,options) result(sol)
        
            implicit none
        
            interface
                function fcn(t,y) result(dydt)
                    use MathAndProgConstants
                    real(dp),               intent(in)   :: t
                    real(dp), dimension(:), intent(in)   :: y
                    real(dp), dimension(size(y))         :: dydt
                end function fcn
            end interface
            
            real(dp), dimension(2),intent(in) :: tspan
            real(dp), dimension(:),intent(in) :: y0
            real(dp), dimension(size(y0))     :: yend, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12
            real(dp) :: t0, tf, t, h
            integer :: N, k
            
            type(OdeOptions) :: options
            type(OdeSol)     :: sol
            
            N = options%NumberOfIntervals
            
            if (options%OnlyFinalState) then
                allocate(sol%T(2));          sol%T = 0.0D0
                allocate(sol%Y(size(y0),2)); sol%Y = 0.0D0
            else
                allocate(sol%T(N+1));          sol%T = 0.0D0
                allocate(sol%Y(size(y0),N+1)); sol%Y = 0.0D0
            end if
            
            t0 = tspan(1)
            tf = tspan(2)
            
            sol%T(1) = t0
            sol%Y(:,1) = y0
            
            h = (tf-t0)/N
            
            t    = t0
            yend = y0
            do k = 1,N
                
                k1  = fcn(t,yend);
                k2  = fcn(t+1.0D0/9.0D0*h,  yend+1.0D0/9.0D0 *h*k1);
                k3  = fcn(t+1.0D0/6.0D0*h,  yend+1.0D0/2.40D1*h*(k1+3.0D0*k2));
                k4  = fcn(t+1.0D0/4.0D0*h,  yend+1.0D0/1.60D1*h*(k1+3.0D0*k3));
                k5  = fcn(t+1.0D0/1.0D1*h,  yend+1.0D0/5.00D2*h*(2.9D1*k1+3.3D1*k3-1.2D1 *k4));
                k6  = fcn(t+1.0D0/6.0D0*h,  yend+1.0D0/9.72D2*h*(3.3D1*k1+4.0D0*k4+1.25D2*k5));
                k7  = fcn(t+1.0D0/2.0D0*h,  yend+1.0D0/3.60D1*h*(-2.1D1*k1+7.6D1*k4+1.25D2*k5-1.62D2*k6));
                k8  = fcn(t+2.0D0/3.0D0*h,  yend+1.0D0/2.43D2*h*(-3.0D1*k1-3.2D1*k4+1.25D2*k5+9.9D1 *k7));
                k9  = fcn(t+1.0D0/3.0D0*h,  yend+1.0D0/3.24D2*h*(1.175D3*k1-3.456D3*k4-6.25D3*k5+8.424D3*k6+2.42D2*k7-2.7D1*k8));
                k10 = fcn(t+5.0D0/6.0D0*h,  yend+1.0D0/3.24D2*h*(2.93D2*k1-8.52D2*k4-1.375D3*k5+1.836D3*k6-1.18D2*k7+1.62D2*k8+3.24D2*k9));
                k11 = fcn(t+5.0D0/6.0D0*h,  yend+1.0D0/1.62D3*h*(1.303D3*k1-4.260D3*k4-6.875D3*k5+9.990D3*k6+1.030D3*k7+1.62D2*k10));
                k12 = fcn(t+h,              yend+1.0D0/4.428D3*h*(-8.595D3*k1+3.0720D4*k4+4.8750D4*k5-6.6096D4*k6+3.78D2*k7-7.29D2*k8-1.944D3*k9-1.296D3*k10+3.240D3*k11));
    
                t = t0 + h*k
                yend = yend + h*(4.1D1/8.40D2*k1+9.0D0/3.5D1*k6+3.4D1/1.05D2*k7+9.0D0/2.80D2*k8+9.0D0/2.80D2*k9+3.0D0/7.0D1*k10+3.0D0/1.4D1*k11+4.1D1/8.40D2*k12);
                
                if (.NOT. options%OnlyFinalState) then
                    sol%T(k+1) = t
                    sol%Y(:,k+1) = yend
                end if
                
            end do
            
            if (options%OnlyFinalState) then
                sol%T(2) = t
                sol%Y(:,2) = yend
            end if
        
        end function

        ! The Dormand-Prince 5(4) method, 5th order
        function ode45(fcn,tspan,y0,options) result(sol)
        
            implicit none

            interface
                function fcn(t,y) result(dydt)
                    use MathAndProgConstants
                    real(dp),               intent(in)   :: t
                    real(dp), dimension(:), intent(in)   :: y
                    real(dp), dimension(size(y))         :: dydt
                end function fcn
            end interface

            real(dp), dimension(2),intent(in) :: tspan
            real(dp), dimension(:),intent(in) :: y0
            type(OdeOptions),  intent(in) :: options
            type(OdeSol)                  :: sol

            logical :: done, nofail
            integer :: chunk, nout
            real(dp) :: t, t0, tf, tnew, pow, h, rh, hmin, hmax, rtol, threshold, LocalErrEstimate, temp
            real(dp), dimension(6)          :: c, hc
            real(dp), dimension(7)          :: E
            real(dp), dimension(7,6)        :: AT, hAT
            real(dp), dimension(size(y0))   :: y, ynew, f0
            real(dp), dimension(size(y0),7) :: K

            ! Parsing, allocation
            t0 = tspan(1)
            tf = tspan(2)
            hmax = 0.1D0*(tf-t0)
            if (options%OnlyFinalState) then
                allocate(sol%T(2));          sol%T = 0.0D0
                allocate(sol%Y(size(y0),2)); sol%Y = 0.0D0
            else
                chunk = options%chunk
                allocate(sol%T(chunk));          sol%T = 0.0D0
                allocate(sol%Y(size(y0),chunk)); sol%Y = 0.0D0
                nout = 1
            end if
            sol%T(1)   = t0
            sol%Y(:,1) = y0
            f0 = fcn(t0,y0)
            rtol = options%RelTol
            threshold = options%AbsTol/options%RelTol
            K = 0.0D0

            ! Current node
            t = t0
            y = y0

            pow = 1.0D0/5.0D0
            c = [ 1.0D0/5.0D0, 3.0D0/10.0D0, 4.0D0/5.0D0, 8.0D0/9.0D0, 1.0D0, 1.0D0 ]
            AT(1,:) = [1.0D0/5.0D0,  3.0D0/40.0D0,    44.0D0/45.0D0,   19372.0D0/6561.0D0,    9017.0D0/3168.0D0,       35.0D0/384.0D0]
            AT(2,:) = [0.0D0,        9.0D0/40.0D0,   -56.0D0/15.0D0,  -25360.0D0/2187.0D0,    -355.0D0/33.0D0,          0.0D0]
            AT(3,:) = [0.0D0,        0.0D0,           32.0D0/9.0D0,    64448.0D0/6561.0D0,   46732.0D0/5247.0D0,      500.0D0/1113.0D0]
            AT(4,:) = [0.0D0,        0.0D0,            0.0D0,           -212.0D0/729.0D0,       49.0D0/176.0D0,       125.0D0/192.0D0]
            AT(5,:) = [0.0D0,        0.0D0,            0.0D0,              0.0D0,            -5103.0D0/18656.0D0,   -2187.0D0/6784.0D0]
            AT(6,:) = [0.0D0,        0.0D0,            0.0D0,              0.0D0,                0.0D0,                11.0D0/84.0D0]
            AT(7,:) = [0.0D0,        0.0D0,            0.0D0,              0.0D0,                0.0D0,                 0.0D0]
            E = [ 71.0D0/57600.0D0, 0.0D0, -71.0D0/16695.0D0, 71.0D0/1920.0D0, -17253.0D0/339200.0D0, 22.0D0/525.0D0, -1.0D0/40.0D0 ]
            hmin = 16.0D0*epsilon(t)

            ! Initial step
            h = hmax
            rh = norminf(f0/max(abs(y),threshold))/(0.8D0*rtol**pow)
            if (h*rh > 1.0D0) then
                h = 1.0D0/rh
            end if
            h = max(h,hmin)
            K(:,1) = f0

            ! Main loop
            done = .false.
            do while (.not. done)

                ! Step correction
                hmin = 16.0D0*epsilon(t)
                h = min(hmax, max(hmin, h)); !write(*,'(e13.6)') h
                if (1.1D0*h >= abs(tf - t)) then
                    h = tf - t
                    done = .true.
                end if

                ! Loop for advancing a step
                nofail = .true.
                do while (.true.)

                    hc  = h*c
                    hAT = h*AT
                    K(:,2) = fcn(t+hc(1),y+matmul(K,hAT(:,1)));
                    K(:,3) = fcn(t+hc(2),y+matmul(K,hAT(:,2)));
                    K(:,4) = fcn(t+hc(3),y+matmul(K,hAT(:,3)));
                    K(:,5) = fcn(t+hc(4),y+matmul(K,hAT(:,4)));
                    K(:,6) = fcn(t+hc(5),y+matmul(K,hAT(:,5)));

                    tnew = t + h;
                    if (done) then
                        tnew = tf
                    end if
                    h = tnew - t

                    ynew = y + matmul(K,hAT(:,6))

                    K(:,7) = fcn(tnew,ynew)

                    ! Error estimation
                    LocalErrEstimate = h * norminf(matmul(K,E) / max(max(abs(y),abs(ynew)),threshold))
                    if (LocalErrEstimate > rtol) then

                        if (h <= hmin) then
                            h = hmin
                        end if

                        if (nofail) then   ! soft correction
                            nofail = .false.
                            h = max(hmin,h*max(0.1D0,0.8D0*(rtol/LocalErrEstimate)**pow))
                        else
                            h = max(hmin,0.5D0*h)  ! rough correction
                        end if
                        done = .false.

                    else ! Successfull step

                        exit

                    end if

                end do ! Loop for advancing a step

                if (.not. options%OnlyFinalState) then
                    nout = nout + 1
                    if (nout > size(sol%T)) then
                        sol%T = [ sol%T, zeros(1,chunk) ]
                        sol%Y = ConCol(sol%Y,zeros(size(y0),chunk))
                    end if
                    sol%T(nout) = tnew
                    sol%Y(:,nout) = ynew
                end if

                if (done) then
                    exit
                end if

                if (nofail) then
                    temp = 1.25D0*(LocalErrEstimate/rtol)**pow
                    if (temp > 0.2D0) then
                        h = h/temp
                    else
                        h = h*5.0D0
                    end if
                end if

                t = tnew
                y = ynew
                K(:,1) = K(:,7) ! FSAL property

            end do ! Main loop

            if (.not. options%OnlyFinalState) then
                sol%T = sol%T(1:nout)
                sol%Y = sol%Y(:,1:nout)
            else
                sol%T(2) = tnew
                sol%Y(:,2) = ynew
            end if

        end function
        
        ! The Adams variable-step variable-order method, maximum 13 order
        function ode113(fcn,tspan,y0,options) result(sol)
        
            implicit none

            interface
                function fcn(t,y) result(dydt)
                    use MathAndProgConstants
                    real(dp),               intent(in)   :: t
                    real(dp), dimension(:), intent(in)   :: y
                    real(dp), dimension(size(y))         :: dydt
                end function fcn
            end interface
            
            ! Input and output variables
            real(dp), dimension(2),intent(in) :: tspan
            real(dp), dimension(:),intent(in) :: y0
            type(OdeOptions),  intent(in)    :: options
            type(OdeSol)                     :: sol
            
            ! real(8) variables
            real(dp) :: t0, tf, hmax, hmin, absh, h, hlast, t, tlast, gstar(13), rh, invwt(size(y0)), temp1, temp2, atol, rtol, threshold, err
            real(dp) :: phi(size(y0),14), psi(12), alpha(12), beta(12), sig(13), g(13), temp3, erk, erkm1, erkm2, reduce, erkp1, temp4(1)
            real(dp), dimension(size(y0)) :: f0, y, yp, yout, p, phikp1
            real(dp), allocatable :: v(:), w(:)
            
            ! integer variables
            integer :: neq, maxk, two(13), k, ns, klast, failed, i, iq, j, knew, kold, nout, chunk, nfailed
            integer, allocatable :: Kb(:)
            
            ! logical variables
            logical phase1, done
            
            neq = size(y0)
            
            t0 = tspan(1)
            tf = tspan(2)
            
            if (options%OnlyFinalState) then
                allocate(sol%T(2));          sol%T = 0.0D0
                allocate(sol%Y(neq,2));      sol%Y = 0.0D0
            else
                chunk = options%chunk
                allocate(sol%T(chunk));     sol%T = 0.0D0
                allocate(sol%Y(neq,chunk)); sol%Y = 0.0D0
                nout = 1
            end if
            sol%T(1)   = t0
            sol%Y(:,1) = y0
            
            f0 = fcn(t0,y0)
            
            atol = options%AbsTol
            rtol = options%RelTol
            threshold = atol / rtol
            
            hmax = 0.1D0*(tf-t0)
            
            t = t0
            y = y0
            yp = f0
            yout = y
            
            ! Initialize method parameters
            maxk = 12
            two = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192 ]
            gstar = [ 0.500000D0,  0.083300D0, 0.041700D0, 0.026400D0, &
                      0.018800D0,  0.014300D0, 0.011400D0, 0.009360D0, &
                      0.007890D0,  0.006790D0, 0.005920D0, 0.005240D0, 0.004680D0 ]
        
            hmin = 16.0D0*epsilon(t)
            absh = min(hmax, tf-t0)
            rh = norminf(yp / max(abs(y),threshold)) / (0.250D0 * sqrt(rtol))
            if (absh * rh > 1.0D0) then
	            absh = 1.0D0 / rh
            end if
            absh = max(absh, hmin)
            
            ! Initialize
            k = 1
            allocate(Kb(1))
            Kb(1) = 1
            phi = zeros(neq,14)
            phi(:,1) = yp
            psi = 0.0D0
            alpha = 0.0D0
            beta = 0.0D0
            sig = 0.0D0
            sig(1) = 1.0D0
            allocate(w(1))
            allocate(v(1))
            g = 0.0D0
            g(1) = 1.0D0
            g(2) = 0.5D0
            
            ns    = 0
            nfailed = 0
            hlast = 0.0D0
            klast = 0
            phase1 = .true.
            
            ! The main loop
            done = .false.
            do while (.not. done)
                
                hmin = 16.0D0*epsilon(t)
                absh = min(hmax, max(hmin, absh))
                h = absh
                
                ! Stretch the step if within 10% of tf-t
                if (1.1D0*absh >= abs(tf - t)) then
                    h = tf - t
                    absh = abs(h)
                    done = .true.
                end if
                
                ! Loop for advancing one step
                failed = 0
                invwt = 1.0D0 / max(abs(y),threshold)
                
                do while (.TRUE.)

                    if (h .ne. hlast) then
                        ns = 0
                    end if
                    
                    if (ns .le. klast) then
                        ns = ns + 1
                    end if
                
                    if (k .ge. ns) then
                        beta(ns) = 1.0D0
                        alpha(ns) = 1.0D0 / ns
                        temp1 = h * ns
                        sig(ns+1) = 1.0D0
                        do i = ns+1,k
                            temp2 = psi(i-1)
                            psi(i-1) = temp1
                            temp1 = temp2 + h
                            beta(i) = beta(i-1) * psi(i-1) / temp2
                            alpha(i) = h / temp1
                            sig(i+1) = i * alpha(i) * sig(i)
                        end do
                        psi(k) = temp1

                        ! Compute coefficients g
                        if (ns .eq. 1) then
                            v = 1.0D0 / (Kb * (Kb + 1.0D0))
                            w = v
                        else
                            ! If order was raised, update diagonal part of v
                            if (k > klast) then
                                temp4(1) = 1.0D0 / (k * (k+1.0D0))
                                v = [v, temp4]
                                do j = 1,ns-2
                                    v(k-j) = v(k-j) - alpha(j+1) * v(k-j+1);
                                end do
                            end if
                            ! Update v and set w
                            do iq = 1,k+1-ns
                                v(iq) = v(iq) - alpha(ns) * v(iq+1)
                                w(iq) = v(iq)
                            end do
                            g(ns+1) = w(1)
                        end if

                        ! Compute g in the work vector w
                        do i = ns+2,k+1
                            do iq = 1,k+2-i
                                w(iq) = w(iq) - alpha(i-1) * w(iq+1)
                            end do
                            g(i) = w(1)
                        end do
                    end if 

                    ! Change phi to phi star
                    do i = ns+1,k
                        do j = 1,neq
                            phi(j,i) = phi(j,i) * beta(i)
                        end do
                    end do

                    ! Predict solution and differences.
                    phi(:,k+2) = phi(:,k+1)
                    phi(:,k+1) = 0.0D0
                    p = 0.0D0
                    do i = k,1,-1
                        p = p + g(i) * phi(:,i)
                        phi(:,i) = phi(:,i) + phi(:,i+1)
                    end do
    
                    p = y + h * p
                    tlast = t
                    t = tlast + h
                    if (done) then
                        t = tf
                    end if

                    yp = fcn(t,p)

                    ! Estimate errors at orders k, k-1, k-2.
                    phikp1 = yp - phi(:,1)
    
                    temp3 = norminf(phikp1 * invwt)
                    err = absh * (g(k) - g(k+1)) * temp3;
                    erk = absh * sig(k+1) * gstar(k) * temp3;
                    if (k .ge. 2) then
                        erkm1 = absh * sig(k) * gstar(k-1) * norminf((phi(:,k)+phikp1) * invwt)
                    else
                        erkm1 = 0.0D0
                    end if
                    if (k .ge. 3) then
                        erkm2 = absh * sig(k-1) * gstar(k-2) * norminf((phi(:,k-1)+phikp1) * invwt)
                    else
                        erkm2 = 0.0D0
                    end if
    
                    ! Test if order should be lowered
                    knew = k;
                    if ((k .eq. 2) .and. (erkm1 .le. 0.50D0*erk)) then
                        knew = k - 1
                    end if
                    if ((k > 2) .and. (max(erkm1,erkm2) .le. erk)) then
                        knew = k - 1
                    end if
   
                    ! Test if step successful
                    if (err > rtol) then     
                        
                        nfailed = nfailed + 1
                        
                        if (absh .le. hmin) then
                            write(*,*) "ODE113: Unable to meet tolerance"
                            return
                        end if
      
                        ! Restore t, phi, and psi
                        phase1 = .false.
                        t = tlast
                        do i = 1,k
                            phi(:,i) = (phi(:,i) - phi(:,i+1)) / beta(i)
                        end do
                        do i = 2,k
                            psi(i-1) = psi(i) - h
                        end do

                        failed = failed + 1
                        reduce = 0.50D0
                        if (failed .eq. 3) then
                            knew = 1
                        elseif (failed > 3) then
                            reduce = min(0.5D0, sqrt(0.5D0*rtol/erk))
                        end if
                        absh = max(reduce * absh, hmin)
                        h = absh
                        k = knew
                        Kb = ColonOperator(1,k,1)
                        done = .false.
                    else ! Successful step
                        exit
                    end if
                end do
                
                klast = k
                hlast = h

                ! Correct and evaluate
                y = p + h * g(k+1) * phikp1
                yp = fcn(t,y)
  
                ! Update differences for next step
                phi(:,k+1) = yp - phi(:,1)
                phi(:,k+2) = phi(:,k+1) - phi(:,k+2)
                do i = 1,k
                    phi(:,i) = phi(:,i) + phi(:,k+1)
                end do

                if ((knew .eq. k-1) .OR. (k .eq. maxk)) then
                    phase1 = .false.
                end if

                ! Select a new order.
                kold = k
                if (phase1) then             ! Always raise the order in phase1
                    k = k + 1
                elseif (knew .eq. k-1) then  ! Already decided to lower the order
                    k = k - 1
                    erk = erkm1
                elseif (k+1 <= ns) then          ! Estimate error at higher order
                    erkp1 = absh * gstar(k+1) * norminf(phi(:,k+2) * invwt)
                    if (k .eq. 1) then
                        if (erkp1 < 0.50D0*erk) then
                            k = k + 1
                            erk = erkp1
                        end if
                    else
                        if (erkm1 .le. min(erk,erkp1)) then
                            k = k - 1
                            erk = erkm1
                        elseif ((k < maxk) .and. (erkp1 < erk)) then
                            k = k + 1
                            erk = erkp1
                        end if
                    end if
                end if
                if (k .ne. kold) then
                    Kb = ColonOperator(1,k,1)
                end if
    
                yout = y
                
                if (.not. options%OnlyFinalState) then
                    nout = nout + 1
                    if (nout > size(sol%T)) then
                        sol%T = [ sol%T, zeros(1,chunk) ]
                        sol%Y = ConCol(sol%Y,zeros(neq,chunk))
                    end if
                    sol%T(nout) = t
                    sol%Y(:,nout) = yout
                end if
  
                if (done) then
                    exit
                end if
  
                ! Select a new step size.
                if (phase1) then
                    absh = 2.00D0 * absh
                elseif (0.50D0*rtol .GE. erk*two(k+1)) then
                    absh = 2.00D0 * absh    
                elseif (0.50D0*rtol < erk) then
                    reduce = (0.50D0 * rtol / erk)**(1.00D0 / (k+1.00D0))
                    absh = absh * max(0.50D0, min(0.90D0, reduce))
                end if
                
            end do ! while ~done
            
            deallocate(Kb)
            deallocate(v)
            deallocate(w)
            
            if (.not. options%OnlyFinalState) then
                sol%T = sol%T(1:nout)
                sol%Y = sol%Y(:,1:nout)
            else
                sol%T(2) = t
                sol%Y(:,2) = yout
            end if

        end function ode113

        ! The Dormand-Prince 8(7) method, 8th order
        function ode87(fcn,tspan,y0,options) result(sol)
        
            implicit none

            interface
                function fcn(t,y) result(dydt)
                    use MathAndProgConstants
                    real(dp),               intent(in)   :: t
                    real(dp), dimension(:), intent(in)   :: y
                    real(dp), dimension(size(y))         :: dydt
                end function fcn
            end interface
            
            ! Input and output variables
            real(dp), dimension(2),intent(in) :: tspan
            real(dp), dimension(:),intent(in) :: y0
            type(OdeOptions),  intent(in)    :: options
            type(OdeSol)                     :: sol
        
            ! real(8) variables
            real(dp) :: rtol, pow, t, t0, tfinal, y(size(y0)), f(size(y0),13), threshold
            real(dp) :: c(12), A(13,12), b_8(13), b_7(13), h, hmin
            real(dp) :: yn8(size(y0)), yn7(size(y0)), err, temp
            
            ! integer variables
            integer :: nout, neq, chunk, j
            
            ! logical variables
            logical :: nofailed, done

            ! ��������� �����������
            rtol = options%RelTol
            threshold = options%AbsTol/rtol
            pow  = 1.0D0/8.0D0
            neq  = size(y0)

            t0     = tspan(1)
            tfinal = tspan(2)
            y      = y0
            t      = t0
            f      = 0.0D0

            ! ������������ ������ (P. J. Prince and J. R. Dormand, High order embedded Runge-Kutta formulae, 1981)
            c = [ 1.0D0/1.8D1, 1.0D0/1.2D1, 1.0D0/8.0D0, 5.0D0/1.6D1, 3.0D0/8.0D0, 5.9D1/4.0D2, 9.3D1/2.0D2, 5.490023248D0/9.719169821D0, 1.3D1/2.0D1, 1.201146811D0/1.299019798D0, 1.0D0, 1.0D0];

            A(1,:) = [ 1.0D0/1.8D1, 1.0D0/4.8D1, 1.0D0/3.2D1, 5.0D0/1.6D1, 3.0D0/8.0D1, 2.9443841D0/6.14563906D1, 1.6016141D0/9.46692911D1, 3.9632708D0/5.73591083D1, 2.46121993D0/1.340847787D1, -1.028468189D1/8.46180014D0, 1.85892177D0/7.18116043D0, 4.03863854D0/4.91063109D0 ]
            A(2,:) = [ 0.0D0, 1.0D0/1.6D1, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 ]
            A(3,:) = [ 0.0D0, 0.0D0, 3.0D0/3.2D1, -7.5D0/6.4D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
            A(4,:) = [ 0.0D0, 0.0D0, 0.0D0, 7.5D0/6.4D0, 3.D0/1.6D1, 7.7736538D0/6.92538347D1, 6.1564180D0/1.58732637D1, -4.33636366D0/6.83701615D0, -3.7695042795D0/1.5268766246D0, 8.478235783D1/5.08512852D0, -3.185094517D1/6.67107341D0, -5.068492393D1/4.34740067D0 ]
            A(5,:) = [ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 3.0D0/2.0D1, -2.8693883D0/1.125D2, 2.2789713D0/6.33445777D1, -4.21739975D0/2.616292301D1,             -3.09121744D0/1.061227803D1, 1.311729495D0/1.432422823D0, -4.77755414D0/1.098053517D1, -4.11421997D0/5.43043805D0 ]
            A(6,:) = [ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 2.3124283D0/1.8D2, 5.45815736D0/2.771057229D1, 1.00302831D0/7.23423059D0, -1.2992083D0/4.90766935D1, -1.0304129995D1/1.701304382D0, -7.03635378D0/2.30739211D0, 6.52783627D0/9.14296604D0 ]
            A(7,:) = [ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -1.80193667D0/1.043307555D1, 7.90204164D0/8.39813087D0, 6.005943493D0/2.108947869D0, -4.8777925059D1/3.047939560D0, 5.731566787D0/1.027545527D0, 1.1173962825D2/9.25320556D0 ]
            A(8,:) = [ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 8.00635310D0/3.783071287D1, 3.93006217D0/1.396673457D1, 1.5336726248D1/1.032824649D0, 5.232866602D1/8.50066563D0, -1.3158990841D1/6.184727034D0 ]
            A(9,:) = [ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.23872331D0/1.001029789D1, -4.5442868181D1/3.398467696D0, -4.093664535D1/8.08688257D0, 3.936647629D0/1.978049680D0 ]
            A(10,:)= [ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 3.065993473D1/5.97172653D0, 3.962137247D0/1.805957418D0, -1.60528059D0/6.85178525D0 ]
            A(11,:)= [ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 6.5686358D0/4.87910083D1, 2.48638103D0/1.413531060D1 ]
            A(12,:)=   0.0D0
            A(13,:)=   0.0D0
            
            b_8 = [ 1.4005451D0/3.35480064D1, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -5.9238493D0/1.068277825D2, 1.81606767D0/7.58867731D0,   5.61292985D0/7.97845732D0,   -1.041891430D0/1.371343529D0,  7.60417239D0/1.151165299D1, 1.18820643D0/7.51138087D0, -5.28747749D0/2.220607170D1,  1.0D0/4.0D0 ]
            b_7 = [ 1.3451932D0/4.55176623D1, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -8.08719846D0/9.76000145D0, 1.757004468D0/5.645159321D0, 6.56045339D0/2.65891186D0,   -3.867574721D0/1.518517206D0,  4.65885868D0/3.22736535D0,  5.3011238D0/6.67516719D1,  2.0D0/4.5D1, 0.0D0]

            ! ��������� ���
            h = min((tfinal - t0)/50.0D0, 0.1D0)

            ! ������������
            if (options%OnlyFinalState) then
                allocate(sol%T(2));          sol%T = 0.0D0
                allocate(sol%Y(neq,2));      sol%Y = 0.0D0
            else
                chunk = options%chunk
                allocate(sol%T(chunk));     sol%T = 0.0D0
                allocate(sol%Y(neq,chunk)); sol%Y = 0.0D0
                nout = 1
            end if
            sol%T(1)   = t0
            sol%Y(:,1) = y0

            ! ������� ����
            nofailed = .TRUE.
            done = .FALSE.
            do while (.NOT. done)
    
                hmin = 16.0D0*epsilon(t)
	            if ((t + h) > tfinal) then
                    h = tfinal - t
                    done = .TRUE.
	            end if

                ! ���������� ������
	            f(:,1) = fcn(t,y)
	            do j = 1,12
                    f(:,j+1) = fcn(t+c(j)*h,y+h*matmul(f,A(:,j)))
	            end do

                ! ��� �������
                yn8 = y + h*matmul(f,b_8)
	            yn7 = y + h*matmul(f,b_7)

                ! ������ ������
                err = norminf((yn7-yn8)/max(max(abs(y),abs(yn8)),threshold))
                if ((nofailed) .AND. (err .LE. rtol)) then
        
                    ! ��� ������� � ��������
                    t = t + h
                    y = yn8
                    
                    if (.not. options%OnlyFinalState) then
                        nout = nout + 1
                        if (nout > size(sol%T)) then
                            sol%T = [ sol%T, zeros(1,chunk) ]
                            sol%Y = ConCol(sol%Y,zeros(size(y0),chunk))
                        end if
                        sol%T(nout) = t
                        sol%Y(:,nout) = y
                    end if
        
                    ! ����������� ���, �� �� ������� ������
                    temp = 1.25D0*(err/rtol)**pow
                    if (temp > 1.0D0/6.0D0) then
                        h = h / temp
                    else
                        h = 6.0D0*h
                    end if
                    
                    nofailed = .TRUE.
        
                else if ((.NOT. nofailed) .AND. (err .LE. rtol)) then
        
                    ! ��� ������� � ��������
                    t = t + h
                    y = yn8
                    
                    if (.not. options%OnlyFinalState) then
                        nout = nout + 1
                        if (nout > size(sol%T)) then
                            sol%T = [ sol%T, zeros(1,chunk) ]
                            sol%Y = ConCol(sol%Y,zeros(size(y0),chunk))
                        end if
                        sol%T(nout) = t
                        sol%Y(:,nout) = y
                    end if

                    nofailed = .TRUE.
        
                else if ((nofailed) .AND. (err > rtol)) then
        
                    ! ������� ��������� ������� ������, ������ �������� ���
                    nofailed = .FALSE.
                    h = max(hmin, h * max(0.1D0, 7.5D-1*(rtol/err)**pow))

                else if ((.NOT. nofailed) .AND. (err > rtol)) then
        
                    ! ����� ��������� ���
                    h = max(hmin, 0.5D0 * h)
        
                end if
                
            end do
        
            if (.not. options%OnlyFinalState) then
                sol%T = sol%T(1:nout)
                sol%Y = sol%Y(:,1:nout)
            else
                sol%T(2) = t
                sol%Y(:,2) = y
            end if
        
        end function ode87
        
end module