module Translations
    
    use MathAndProgConstants
    use LinearAlgebraInterfaces
    
    type TranslationStructure
        real(dp) :: NewVec(3)
        real(dp) :: dNewVec(3,3)
    end type
    
    contains
    
        ! Vars transltaions.
        function kee2rv(ee,mu) result(rv)
        
            implicit none
            
            real(dp) :: ee(6), mu, rv(6)
            real(dp) :: h, ex, ey, ix, iy, L, p, ksi
            real(dp) :: cosL, sinL, K, f(3), g(3), ix2, iy2, ixiy

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

            f = K * (/ 1.0D0+ix2-iy2, 2.0D0*ixiy, -2.0D0*iy /)
            g = K * (/ 2.0D0*ixiy, 1.0D0-ix2+iy2,  2.0D0*ix /)

            rv = (/ h*(f*cosL + g*sinL)/ksi, (ex+cosL)/p*g - (ey+sinL)/p*f /)
    
        end function kee2rv
        function krv2ee(rv,mu) result(ee)
        
            use BaseMeansToolbox
        
            implicit none
            
            real(dp) :: rv(6), mu, ee(6)
            real(dp) :: rvect(3), vvect(3), cvect(3)
            real(dp) :: r, c, Avect(3), f(3), g(3)
            real(dp) :: L1, L2, h, ex, ey, ix, iy, L

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

            ee = (/ h, ex, ey, ix, iy, L /)
            if (ee(6) < 0) then
                ee(6) = ee(6) + 2*pi
            end if
    
        end function krv2ee
        
        ! System tranlations.
        function kscrs2pa(SCRS_vec,JD) result(kscrs2pa_output)
        
            use Ephemeris
            use LinearAlgebraInterfaces
            
            implicit none
            
            real(dp) :: angles(3), cos1, cos2, cos3, sin1, sin2, sin3
            real(dp) :: jd, SCRS_vec(3), PA_vec(3), M(3,3)
            type(TranslationStructure) :: kscrs2pa_output

            angles = MoonLibration(JD)

            cos1 = cos(angles(1))
            cos2 = cos(angles(2))
            cos3 = cos(angles(3))
            sin1 = sin(angles(1))
            sin2 = sin(angles(2))
            sin3 = sin(angles(3))

            M(1,:) = (/ cos1*cos3-sin1*cos2*sin3,   sin1*cos3+cos1*cos2*sin3,  sin2*sin3/)
            M(2,:) = (/-cos1*sin3-sin1*cos2*cos3,  -sin1*sin3+cos1*cos2*cos3,  sin2*cos3/)
            M(3,:) = (/ sin1*sin2, -cos1*sin2, cos2 /)

            PA_vec = dotmv(M, SCRS_vec)
            
            kscrs2pa_output.NewVec = PA_vec
            kscrs2pa_output.dNewVec = M

            ! � ������� [PA_vec,invM] �� �����
            ! ����� �� ����� [Pa_vec, M] !!!

        end function kscrs2pa
        function kitrs2gcrs(XITRS,JD) result(kitrs2gcrs_output)
        
            use LinearAlgebraInterfaces

            implicit none
            
            real(dp) :: XITRS(3), XGCRS(3), dXGCRS(3,3), JD, Om, t, X, Y, JD_UTC, Q(3,3), R(3,3)
            real(dp), parameter :: leapsec = 37.0D0
            type(TranslationStructure) :: kitrs2gcrs_output

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
            
            kitrs2gcrs_output.NewVec = XGCRS
            kitrs2gcrs_output.dNewVec = dXGCRS

        end function kitrs2gcrs
        function kgcrs2itrs(XGCRS,JD) result(kgcrs2itrs_output)
        
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp) :: XGCRS(3), JD, dXITRS(3,3), XITRS(3)
            type(TranslationStructure) :: kitrs2gcrs_output, kgcrs2itrs_output
            
            kitrs2gcrs_output = kitrs2gcrs((/1.0D0,0.0D0,0.0D0/),JD)
            dXITRS = transpose(kitrs2gcrs_output.dNewVec)
            XITRS = dotmv(dXITRS,XGCRS)
            
            kgcrs2itrs_output.NewVec = XITRS
            kgcrs2itrs_output.dNewVec = dXITRS

        end function kgcrs2itrs
        
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

        function translations_about() result(i)
            integer :: i
            i = 0
        end function translations_about

end module