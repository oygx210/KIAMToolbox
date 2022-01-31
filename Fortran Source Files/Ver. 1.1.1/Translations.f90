module Translations
    
    use MathAndProgConstants
    use LinearAlgebraInterfaces
    
    contains
    
        ! About function, opened for F2PY wrapping
        function translations_about() result(i)
            integer :: i
            i = 0
            write(*,*) 'Translations: Version 1.0.'
            write(*,*) 'Author: Maksim Shirobokov.'
            write(*,*) 'Date: 19.01.2022.'
        end function translations_about
    
        ! Vars transltaions, opened for F2PY wrapping.
        subroutine kee2rv(ee, mu, grad_req, rv, drv)
        
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp), intent(in) :: ee(6), mu
            logical(dp), intent(in) :: grad_req
            real(dp), intent(out) :: rv(6), drv(6,6)
            real(dp) :: h, ex, ey, ix, iy, L, p, ksi
            real(dp) :: cosL, sinL, K, f(3), g(3), ix2, iy2, ixiy
            real(dp) :: dh(6), dex(6), dey(6), dix(6), diy(6), dL(6), dpp(6), dksi(6), dK(6)
            real(dp) :: T1(3,6), T2(3,6), df(3,6), dg(3,6)

            h  = ee(1)
            ex = ee(2)
            ey = ee(3)
            ix = ee(4)
            iy = ee(5)
            L  = ee(6)
            
            cosL = cos(L)
            sinL = sin(L)
            ix2 = ix**2
            iy2 = iy**2
            ixiy = ix*iy

            p = h**2*mu
            ksi = 1 + ex*cosL + ey*sinL
    
            K = h*mu/(1+ix**2+iy**2)

            f = K * [ 1.0D0+ix2-iy2, 2.0D0*ixiy, -2.0D0*iy ]
            g = K * [ 2.0D0*ixiy, 1.0D0-ix2+iy2,  2.0D0*ix ]

            rv = [ h*(f*cosL + g*sinL)/ksi, (ex+cosL)/p*g - (ey+sinL)/p*f ]
            
            drv = 0.0D0
            if (grad_req) then
        
                dh =  [1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
                dex = [0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
                dey = [0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0]
                dix = [0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0]
                diy = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0]
                dL =  [0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0]

                dpp = 2.0D0 * h * dh * mu
                dksi = cosL*dex - ex*sinL*dL + sinL*dey + ey*cosL*dL

                dK = mu/(1+ix**2+iy**2)*dh - h*mu/(1+ix**2+iy**2)**2*(2*ix*dix+2*iy*diy)
        
                T1(1,:) = 2*ix*dix - 2*iy*diy
                T1(2,:) = 2*iy*dix + 2*ix*diy
                T1(3,:) = -2*diy
                T2(1,:) = T1(2,:)
                T2(2,:) = -T1(1,:)
                T2(3,:) = 2*dix
                df = doto([1+ix**2-iy**2, 2*ix*iy, -2*iy], dK) + K * T1
                dg = doto([2*ix*iy, 1-ix**2+iy**2,  2*ix], dK) + K * T2

                drv(1:3,:) = doto((f*cosL + g*sinL)/ksi, dh) + h*(cosL*df-sinL*doto(f,dL)+sinL*dg+cosL*doto(g,dL))/ksi - h*doto((f*cosL + g*sinL)/ksi**2, dksi)
                drv(4:6,:) = doto(g,dex-sinL*dL)/p - (ex+cosL)/p**2*doto(g,dpp) + (ex+cosL)/p*dg - doto(f,dey+cosL*dL)/p + (ey+sinL)/p**2*doto(f,dpp) - (ey+sinL)/p*df

            end if
    
        end subroutine kee2rv
        subroutine krv2ee(rv, mu, grad_req, ee, dee)
        
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp), intent(in) :: rv(6), mu
            logical, intent(in) :: grad_req
            real(dp), intent(out) :: ee(6), dee(6,6)
            real(dp) :: rvect(3), vvect(3), cvect(3)
            real(dp) :: r, c, Avect(3), f(3), g(3)
            real(dp) :: L1, L2, h, ex, ey, ix, iy, L
            real(dp) :: drvect(3,6), dvvect(3,6), dcvect(3,6), dr(6), dc(6)
            real(dp) :: dT1(3,6), dT2(3,6), dAvect(3,6), df(3,6), dg(3,6)
            real(dp) :: dL1(6), dL2(6)

            rvect = rv(1:3)
            vvect = rv(4:6)
            cvect = cross(rvect,vvect)

            r = norm(rvect)
            c = norm(cvect)

            Avect = cross(vvect,cvect) - mu*rvect/r
            
            f(1) = c-cvect(1)**2/(c+cvect(3))
            f(2) = -cvect(1)*cvect(2)/(c+cvect(3))
            f(3) = -cvect(1)
            
            g(1) = f(2)
            g(2) = c-cvect(2)**2/(c+cvect(3))
            g(3) = -cvect(2)
    
            L1 = sum(rvect*g)
            L2 = sum(rvect*f)

            h  = c/mu
            ex = sum(Avect*f)/c/mu
            ey = sum(Avect*g)/c/mu
            ix = -cvect(2)/(c+cvect(3))
            iy =  cvect(1)/(c+cvect(3))
            L = atan2(L1,L2)

            ee = [ h, ex, ey, ix, iy, L ]
            if (ee(6) < 0) then
                ee(6) = ee(6) + 2*pi
            end if
            
            dee = 0.0D0
            if (grad_req) then
                
                drvect(1,:) = [1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
                drvect(2,:) = [0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
                drvect(3,:) = [0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0]
                
                dvvect(1,:) = [0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0]
                dvvect(2,:) = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0]
                dvvect(3,:) = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0]
                
                dcvect(:,1:3) = -kcrossProdMatrix(vvect)
                dcvect(:,4:6) =  kcrossProdMatrix(rvect)
        
                dr = dotvm(rvect, drvect)/r
                dc = dotvm(cvect, dcvect)/c
        
                dT1 = -mu*drvect/r + mu*doto(rvect,dr)/r**2
                dT2 = dotmm(kcrossProdMatrix(vvect), dcvect) - dotmm(kcrossProdMatrix(cvect), dvvect)
                dAvect = dT1 + dT2
        
                df(1,:) = dc - 2*cvect(1)*dcvect(1,:)/(c+cvect(3)) + cvect(1)**2/(c+cvect(3))**2*(dc + dcvect(3,:))
                df(2,:) = -dcvect(1,:)*cvect(2)/(c+cvect(3)) - cvect(1)*dcvect(2,:)/(c+cvect(3)) + cvect(1)*cvect(2)/(c+cvect(3))**2*(dc+dcvect(3,:))
                df(3,:) = -dcvect(1,:)
        
                dg(1,:) = df(2,:)
                dg(2,:) = dc - 2*cvect(2)*dcvect(2,:)/(c+cvect(3)) + cvect(2)**2/(c+cvect(3))**2*(dc+dcvect(3,:))
                dg(3,:) = -dcvect(2,:)
        
                dee(1,:) = dc/mu
                dee(2,:) = dotvm(f, dAvect)/c/mu + dotvm(Avect, df)/c/mu - sum(Avect * f)/c**2*dc/mu
                dee(3,:) = dotvm(g, dAvect)/c/mu + dotvm(Avect, dg)/c/mu - sum(Avect * g)/c**2*dc/mu
                dee(4,:) = -dcvect(2,:)/(c+cvect(3)) + cvect(2)/(c+cvect(3))**2*(dc+dcvect(3,:))
                dee(5,:) = dcvect(1,:)/(c+cvect(3)) - cvect(1)/(c+cvect(3))**2*(dc+dcvect(3,:))
        
                dL1 = dotvm(g, drvect) + dotvm(rvect, dg)
                dL2 = dotvm(f, drvect) + dotvm(rvect, df)
                
                dee(6,:) = 1.0D0/(1.0D0+(L1/L2)**2)*(dL1/L2 - L1/L2**2*dL2)
        
            end if
    
        end subroutine krv2ee
        
        ! System tranlations, opened for F2PY wrapping.
        subroutine kscrs2pa(XSCRS,JD,XPA,dXPA)
        
            use Ephemeris
            use LinearAlgebraInterfaces
            
            implicit none
            
            real(dp), intent(in) :: XSCRS(3), JD
            real(dp), intent(out) :: XPA(3), dXPA(3,3)
            real(dp) :: angles(3), cos1, cos2, cos3, sin1, sin2, sin3

            angles = MoonLibration(JD)

            cos1 = cos(angles(1))
            cos2 = cos(angles(2))
            cos3 = cos(angles(3))
            sin1 = sin(angles(1))
            sin2 = sin(angles(2))
            sin3 = sin(angles(3))

            dXPA(1,:) = (/ cos1*cos3-sin1*cos2*sin3,   sin1*cos3+cos1*cos2*sin3,  sin2*sin3/)
            dXPA(2,:) = (/-cos1*sin3-sin1*cos2*cos3,  -sin1*sin3+cos1*cos2*cos3,  sin2*cos3/)
            dXPA(3,:) = (/ sin1*sin2, -cos1*sin2, cos2 /)

            XPA = dotmv(dXPA, XSCRS)

        end subroutine kscrs2pa
        subroutine kitrs2gcrs(XITRS,JD,XGCRS,dXGCRS)
        
            use LinearAlgebraInterfaces

            implicit none
            
            real(dp), intent(in) :: XITRS(3), JD
            real(dp), intent(out) :: XGCRS(3), dXGCRS(3,3)
            real(dp) :: Om, t, X, Y, JD_UTC, Q(3,3), R(3,3)
            real(dp), parameter :: leapsec = 37.0D0

            XGCRS = 0.0D0
            dXGCRS = 0.0D0

            t = (JD - 2451545.0D0)/36525.0D0
    
            Om = 125.04455501D0 + ( -6962890.5431D0*t + 7.47220D0*t**2 + 0.007702D0*t**3 - 0.00005939D0*t**4)/3600.0D0
            Om = Om*pi/180.0D0
    
            X = 2004.191898D0*t - 6.84431844D0*sin(Om)
            Y = -22.4072747D0*t**2 + 9.20523626D0*cos(Om)
    
            JD_UTC = JD - (32.184D0 + leapsec)/86400.0D0
    
            Q = getQ(X,Y)
            R = getR(JD_UTC)
            dXGCRS = dotmm(Q,R)
            XGCRS = dotmv(dXGCRS,XITRS)

        end subroutine kitrs2gcrs
        subroutine kgcrs2itrs(XGCRS,JD,XITRS,dXITRS)
        
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp), intent(in) :: XGCRS(3), JD
            real(dp), intent(out) :: XITRS(3), dXITRS(3,3)
            real(dp) :: XGCRS_EYE(3), dXGCRS_EYE(3,3)
            
            call kitrs2gcrs([1.0D0, 0.0D0, 0.0D0], JD, XGCRS_EYE, dXGCRS_EYE)
            dXITRS = transpose(dXGCRS_EYE)
            
            XITRS = dotmv(dXITRS,XGCRS)

        end subroutine kgcrs2itrs
        
        function getQ(X,Y) result(Q)
        
            implicit none
            
            real(dp) :: X, Y, X2, Y2, XY, a, Q(3,3)

            ! Conversion to radians
            X = X/3600.0D0*pi/180.0D0
            Y = Y/3600.0D0*pi/180.0D0

            X2 = X**2
            Y2 = Y**2
            XY = X*Y

            a = 0.5D0 + (X2 + Y2)/8.0D0
            
            Q(1,:) = (/ 1.0D0-a*X2, -a*XY, X /)
            Q(2,:) = (/ -a*XY, 1.0D0-a*Y2, Y /)
            Q(3,:) = (/ -X, -Y, 1.0D0-a*(X2+Y2) /)

        end function getQ
        function getR(JD_UTC) result(R)
        
            implicit none
            
            real(dp) :: JD_UTC, T_u, ERA, c, s, R(3,3)

            T_u = JD_UTC - 2451545.0D0

            ERA = 2*pi*( T_u - floor(T_u) + 0.7790572732640D0 + 0.00273781191135448D0*T_u )

            c = cos(ERA)
            s = sin(ERA)
            
            R(1,:) = (/ c, -s,  0.0D0 /)
            R(2,:) = (/ s,  c,  0.0D0 /)
            R(3,:) = (/ 0.0D0, 0.0D0, 1.0D0 /)

        end function getR

end module