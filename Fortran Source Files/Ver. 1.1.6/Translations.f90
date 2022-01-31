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
    
        ! Vars translations, opened for F2PY wrapping.
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
        subroutine koe2rv(oe, mu, grad_req, rv, drv)
        
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp), intent(in) :: oe(6), mu
            logical, intent(in) :: grad_req
            real(dp), intent(out) :: rv(6), drv(6,6)
            real(dp) :: a, e, inc, Om, w, th, p, r
            real(dp) :: n1(3), n2(3), n3(3), er(3), en(3)
            real(dp) :: dadx(6), dedx(6), dincdx(6), dOmdx(6)
            real(dp) :: dwdx(6), dthdx(6), dpdx(6), drdx(6)
            real(dp) :: dn1dx(3,6), dn2dx(3,6), dn3dx(3,6)
            real(dp) :: derdx(3,6), dendx(3,6)
            real(dp) :: cosOm, sinOm, costh, sinth, cosinc, sininc
            
            a   = OE(1)
            e   = OE(2)
            inc = OE(3)
            Om  = OE(4)
            w   = OE(5)
            th  = OE(6)
    
            p = a*(1-e**2)
            
            cosOm = cos(Om)
            sinOm = sin(Om)
            costh = cos(th)
            sinth = sin(th)
            cosinc = cos(inc)
            sininc = sin(inc)
            
            r = p/(1+e*costh)
    
            n1 = [cosOm, sinOm, 0.0D0]
            n3 = [sinOm*sininc, -cosOm*sininc, cosinc]
            n2 = cross(n3,n1)
    
            er = n1*cos(w+th) + n2*sin(w+th)
            en = cross(n3,er)
    
            rv(1:3) = r * er
            rv(4:6) = sqrt(mu/p) * (e*sinth*er + (1+e*costh)*en)
    
            drv = 0.0D0
            if (grad_req) then
        
                dadx   = [1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
                dedx   = [0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
                dincdx = [0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0]
                dOmdx  = [0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0]
                dwdx   = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0]
                dthdx  = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0]

                dpdx = (1-e**2)*dadx - a*2*e*dedx
                drdx = dpdx/(1+e*costh) - p/(1+e*costh)**2*(dedx*costh - e*sinth*dthdx)
                
                dn1dx(1,:) = -sinOm*dOmdx
                dn1dx(2,:) = cosOm*dOmdx
                dn1dx(3,:) = 0.0D0
                
                dn3dx(1,:) = cosOm*sininc*dOmdx + cosinc*sinOm*dincdx
                dn3dx(2,:) = sinOm*sininc*dOmdx - cosinc*cosOm*dincdx
                dn3dx(3,:) = -sininc*dincdx

                dn2dx = dotmm(kcrossProdMatrix(n3), dn1dx) - dotmm(kcrossProdMatrix(n1), dn3dx)
    
                derdx = dn1dx*cos(w+th) - sin(w+th)*doto(n1, dwdx+dthdx) + dn2dx*sin(w+th) + cos(w+th)*doto(n2, dwdx+dthdx)
                dendx = dotmm(kcrossProdMatrix(n3), derdx) - dotmm(kcrossProdMatrix(er), dn3dx)

                drv(1:3,:) = doto(er, drdx) + r*derdx
                drv(4:6,:) = doto(e*sinth*er + (1.0D0+e*costh)*en, dpdx)/2.0D0/sqrt(mu/p)*(-mu/p**2) + &
                    sqrt(mu/p)*( sinth*doto(er, dedx) + e*costh*doto(er, dthdx) + e*sinth*derdx + &
                    (1.0D0+e*costh)*dendx + doto(en, dedx*costh-e*sinth*dthdx) )
                
            end if
        
        end subroutine koe2rv
        subroutine krv2oe(rv, mu, grad_req, oe, doe)
        
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp), intent(in) :: rv(6), mu
            logical, intent(in) :: grad_req
            real(dp), intent(out) :: oe(6), doe(6,6)
            real(dp) :: rvect(3), vvect(3)
            real(dp) :: r, v, cvect(3), c, h, p, e, a
            real(dp) :: inc, Om, T1, T2, th, T3, T4, w
            real(dp) :: drvect(3,6), dvvect(3,6), dr(6), dv(6)
            real(dp) :: dcvect(3,6), dc(6), dh(6), dpp(6), de(6)
            real(dp) :: da(6), dinc(6), dOm(6), dT1(6), dT2(6), dth(6)
            real(dp) :: dT3(6), dT4(6), dw(6)
            
            rvect = rv(1:3)
            r = sqrt(sum(rvect*rvect))
            
            vvect = rv(4:6)
            v = sqrt(sum(vvect*vvect))
    
            cvect = cross(rvect,vvect)
            c = sqrt(sum(cvect*cvect))
    
            h = v**2 - 2.0D0*mu/r
            p = c**2/mu
    
            e = sqrt(max(0.0D0, 1.0D0 + h*p/mu))
            a = -mu/h
    
            inc = acos(min(1.0D0, max(-1.0D0, cvect(3)/c)))
            Om  = atan2(cvect(1), -cvect(2))
    
            T1 = sum(rvect*vvect)/c
            T2 = 1.0D0 - r/p
            th  = atan2(T1,T2)
    
            T3 = rvect(3)
            T4 = rvect(2)*cvect(1)/c - rvect(1)*cvect(2)/c
            w  = atan2(T3,T4) - th
    
            oe(1) = a
            oe(2) = e
            oe(3) = inc
            oe(4) = Om
            oe(5) = w
            oe(6) = th
            
            doe = 0.0D0
            if (grad_req) then
                
                drvect(1,:) = [1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
                drvect(2,:) = [0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0] 
                drvect(3,:) = [0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0] 
                
                dr = dotvm(rvect, drvect)/r
                
                dvvect(1,:) = [0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0]
                dvvect(2,:) = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0] 
                dvvect(3,:) = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0] 

                dv = dotvm(vvect, dvvect)/v

                dcvect = dotmm(kcrossProdMatrix(rvect), dvvect) - dotmm(kcrossProdMatrix(vvect), drvect)
        
                dc = dotvm(cvect, dcvect)/c

                dh = 2.0D0*v*dv + 2.0D0*mu/r**2*dr
                dpp = 2.0D0*c/mu*dc
                de = (dh*p/mu + h*dpp/mu)/2.0D0/e
                da = mu/h**2*dh
                dinc = -(dcvect(3,:)/c - cvect(3)/c**2*dc)/abs(sin(inc))
                dOm = -1.0D0/(1+cvect(1)**2/cvect(2)**2)*(dcvect(1,:)/cvect(2)-cvect(1)/cvect(2)**2*dcvect(2,:))

                dT1 = dotvm(rvect, dvvect)/c + dotvm(vvect, drvect)/c - sum(rvect*vvect)/c**2*dc
                dT2 = -dr/p + r/p**2*dpp
                dth = 1.0D0/(1.0D0+T1**2/T2**2)*(dT1/T2 - T1/T2**2*dT2)

                dT3 = drvect(3,:)
                dT4 = drvect(2,:)*cvect(1)/c + rvect(2)*dcvect(1,:)/c - rvect(2)*cvect(1)/c**2*dc - &
                      drvect(1,:)*cvect(2)/c - rvect(1)*dcvect(2,:)/c + rvect(1)*cvect(2)/c**2*dc
                dw = 1.0D0/(1+T3**2/T4**2)*(dT3/T4 - T3/T4**2*dT4) - dth

                doe(1,:) = da
                doe(2,:) = de
                doe(3,:) = dinc
                doe(4,:) = dOm
                doe(5,:) = dw
                doe(6,:) = dth
                
            end if
        
        end subroutine krv2oe
        subroutine kcart2sphere(x, q)
		
			use BaseMeansToolbox
		
			implicit none
			
			real(dp), intent(in) :: x(:,:)
			real(dp), intent(out) :: q(3,size(x,2))
			real(dp) :: r(size(x,2))
			
			r = vecnorm(x)
			q(1,:) = r
			q(2,:) = atan2(x(2,:), x(1,:))
			q(3,:) = acos(min(1.0D0,max(-1.0D0,x(3,:)/r)))
		
		end subroutine kcart2sphere
        subroutine ksphere2cart(q, x)
		
			implicit none
			
			real(dp), intent(in) :: q(:,:)
			real(dp), intent(out) :: x(3,size(q,2))
			
			x(1,:) = q(1,:)*cos(q(2,:))*sin(q(3,:))
			x(2,:) = q(1,:)*sin(q(2,:))*sin(q(3,:))
			x(3,:) = q(1,:)*cos(q(3,:))
			
		end subroutine ksphere2cart
		subroutine kcart2latlon(x, latlon)
		
			use BaseMeansToolbox
		
			implicit none
			
			real(dp), intent(in) :: x(:,:)
			real(dp), intent(out) :: latlon(2,size(x,2))
			integer :: i
			real(dp) :: nx(size(x,1),size(x,2))
			
			forall (i=1:size(x,2)) nx(:,i) = x(:,i) / sqrt(sum(x(:,i)*x(:,i)))
			latlon(1,:) = pi/2.0D0 - acos(min(1.0D0,max(-1.0D0,nx(3,:))))
			latlon(2,:) = atan2(nx(2,:), nx(1,:))
		
		end subroutine kcart2latlon
		subroutine klatlon2cart(latlon, x)
		
			implicit none
			
			real(dp), intent(in) :: latlon(:,:)
			real(dp), intent(out) :: x(3,size(latlon,2))
			
			x(1,:) = cos(latlon(1,:))*cos(latlon(2,:))
			x(2,:) = cos(latlon(1,:))*sin(latlon(2,:))
			x(3,:) = sin(latlon(1,:))
		
		end subroutine klatlon2cart
		
		! System translations, opened for F2PY wrapping.
		subroutine kine2rotEph(XINE, JD, FirstBody, SecondBody, DistUnit, VelUnit, XROT)
		
			use BaseMeansToolbox
			use Ephemeris
		
			implicit none
			
			real(dp), intent(in) :: XINE(:,:), JD(size(XINE,2))
			character(*), intent(in) :: FirstBody, SecondBody
			real(dp), intent(in) :: DistUnit, VelUnit
			real(dp), intent(out) :: XROT(size(XINE,1),size(XINE,2))
			integer :: ntimes, i
			real(dp) :: xS(6), RS(3), VS(3), CS(3), Mine2rot(6,6), Omega(3,3), C(3,3)
			real(dp) :: dtheta(size(XROT,2))
		
			!  Translate phase states from inertial CS to the rotating CS (w.r.t. the 2nd body).
			!
			!  XINE       -- (6xN) phase states in the inertial CS
			!  JD         -- (1xN,Nx1) julian dates at which states are identified (in days)
			!  FirstBody  -- name of the first body in three-body system
			!  SecondBody -- name of the second body in the three-body system
			!  DistUnit   -- unit of distance ( in km )
			!  VelUnit    -- unit of velocity ( in km/s )
			!
			!  XROT       -- (6xN) phase states in rotating frame with center at the 2nd body
			!
			! Author: Maksim Shirobokov
			! Date: 30.01.2018
			! UPD: 19.08.2020 -- JDs = JD0 + tspanDays !
			! UPD: 26.01.2022 -- simple analitics for the inverse matrix
			Omega = 0.0D0
			Mine2rot = 0.0D0
			ntimes = size(XINE,2)
			do i = 1, ntimes
			
				xS = PlanetState(JD(i), SecondBody, FirstBody)
				RS = xS(1:3)/DistUnit
				VS = xS(4:6)/VelUnit
				CS = cross(RS, VS)
				
				C(1,:) = RS/sqrt(sum(RS*RS))
				C(3,:) = CS/sqrt(sum(CS*CS))
				C(2,:) = cross(C(3,:), C(1,:))
				
				dtheta(i) = sqrt(sum(CS*CS))/sum(RS*RS)
				
				Omega(1,:) = dtheta(i)*C(2,:)
				Omega(2,:) = -dtheta(i)*C(1,:)
				
				Mine2rot(1:3,1:3) = C
				Mine2rot(4:6,1:3) = Omega
				Mine2rot(4:6,4:6) = C
				
				XROT(:,i) = dotmv(Mine2rot, XINE(:,i))
				
			end do
		
		end subroutine kine2rotEph
		subroutine krot2ineEph(XROT, JD, FirstBody, SecondBody, DistUnit, VelUnit, XINE)
		
			use BaseMeansToolbox
			use Ephemeris
		
			implicit none
			
			real(dp), intent(in) :: XROT(:,:), JD(size(XROT,2))
			character(*), intent(in) :: FirstBody, SecondBody
			real(dp), intent(in) :: DistUnit, VelUnit
			real(dp), intent(out) :: XINE(size(XROT,1),size(XROT,2))
			integer :: ntimes, i
			real(dp) :: xS(6), RS(3), VS(3), CS(3), Mrot2ine(6,6), Omega(3,3), C(3,3)
			real(dp) :: dtheta(size(XROT,2))
		
			!  Translate phase states from rotating CS to the inertial CS (w.r.t. the 2nd body).
			!
			!  XROT       -- (6xN) phase states in the inertial CS
			!  JD         -- (1xN,Nx1) julian dates at which states are identified (in days)
			!  FirstBody  -- name of the first body in three-body system
			!  SecondBody -- name of the second body in the three-body system
			!  DistUnit   -- unit of distance ( in km )
			!  VelUnit    -- unit of velocity ( in km/s )
			!
			!  XINE       -- (6xN) phase states in rotating frame with center at the 2nd body
			!
			! Author: Maksim Shirobokov
			! Date: 30.01.2018
			! UPD: 19.08.2020 -- JDs = JD0 + tspanDays !
			Omega = 0.0D0
			Mrot2ine = 0.0D0
			ntimes = size(XROT,2)
			do i = 1, ntimes
			
				xS = PlanetState(JD(i), SecondBody, FirstBody)
				RS = xS(1:3)/DistUnit
				VS = xS(4:6)/VelUnit
				CS = cross(RS, VS)
				
				C(:,1) = RS/sqrt(sum(RS*RS))
				C(:,3) = CS/sqrt(sum(CS*CS))
				C(:,2) = cross(C(:,3), C(:,1))
				
				dtheta(i) = sqrt(sum(CS*CS))/sum(RS*RS)
				
				Omega(:,1) = dtheta(i)*C(:,2)
				Omega(:,2) = -dtheta(i)*C(:,1)
				
				Mrot2ine(1:3,1:3) = C
				Mrot2ine(4:6,1:3) = Omega
				Mrot2ine(4:6,4:6) = C
				
				XINE(:,i) = dotmv(Mrot2ine, XROT(:,i))
				
			end do
		
		end subroutine krot2ineEph
		subroutine kine2rot(XINE, T, T0, XROT)
		
			implicit none
			
			real(dp), intent(in) :: XINE(:,:), T(size(XINE,2)), T0(:)
			real(dp), intent(out) :: XROT(size(XINE,1),size(XINE,2))
			integer :: ntimes, ntimest0, i
			real(dp) :: A(3,3), dA(3,3), c, s
			
			! t0 - момент(ы) времени, когда вращ. и невращ. с/к совпадают

			! Author: Maksim Shirobokov
			! Date: Mar 24, 2021
			
			ntimest0 = size(T0,1)
			ntimes = size(XINE,2)
			A = 0.0D0
			A(3,3) = 1.0D0
			dA = 0.0D0
			do i = 1, ntimes
			
				if (ntimest0 .eq. 1) then
					c = cos(T(i)-T0(1))
					s = sin(T(i)-T0(1))
				else
					c = cos(T(i)-T0(i))
					s = sin(T(i)-T0(i))
				end if
				
				A(1,1:2) = [ c, s ]
				A(2,1:2) = [-s, c ]
				
				dA(1,1:2) = [ -s,  c ]
				dA(2,1:2) = [ -c, -s ]
    
				XROT(1:3,i) = dotmv(A, XINE(1:3,i))
				XROT(4:6,i) = dotmv(A, XINE(4:6,i)) + dotmv(dA, XINE(1:3,i))
			
			end do
		
		end subroutine kine2rot
		subroutine krot2ine(XROT, T, T0, XINE)
		
			implicit none
			
			real(dp), intent(in) :: XROT(:,:), T(size(XROT,2)), T0(:)
			real(dp), intent(out) :: XINE(size(XROT,1),size(XROT,2))
			integer :: ntimes, ntimest0, i
			real(dp) :: A(3,3), dA(3,3), c, s
			
			! t0 - момент(ы) времени, когда вращ. и невращ. с/к совпадают

			! Author: Maksim Shirobokov
			! Date: Mar 24, 2021
			
			ntimest0 = size(T0,1)
			ntimes = size(XROT,2)
			A = 0.0D0
			A(3,3) = 1.0D0
			dA = 0.0D0
			do i = 1, ntimes
			
				if (ntimest0 .eq. 1) then
					c = cos(T(i)-T0(1))
					s = sin(T(i)-T0(1))
				else
					c = cos(T(i)-T0(i))
					s = sin(T(i)-T0(i))
				end if
				
				A(1,1:2) = [ c, -s ]
				A(2,1:2) = [ s,  c ]
				
				dA(1,1:2) = [ -s, -c ]
				dA(2,1:2) = [  c, -s ]
    
				XINE(1:3,i) = dotmv(A, XROT(1:3,i))
				XINE(4:6,i) = dotmv(A, XROT(4:6,i)) + dotmv(dA, XROT(1:3,i))
			
			end do
		
		end subroutine krot2ine
        subroutine kscrs2pa(XSCRS, JD, XPA, dXPA)
		
        
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
        subroutine kscrs2mer(XSCRS, JD, XMER, dXMER)
		
			use Ephemeris
			use LinearAlgebraInterfaces
		
			implicit none
			
			real(dp), intent(in) :: XSCRS(:,:), JD(size(XSCRS,2))
			real(dp), intent(out) :: XMER(size(XSCRS,1),size(XSCRS,2))
			real(dp), intent(out) :: dXMER(size(XSCRS,1),size(XSCRS,1),size(XSCRS,2))
			integer :: dim, ntimes, i
			real(dp) :: ang(3), c1, c2, c3, s1, s2, s3
			real(dp) :: angles(size(XSCRS,2),3)
			real(dp) :: cos1(size(XSCRS,2)), cos2(size(XSCRS,2)), cos3(size(XSCRS,2))
			real(dp) :: sin1(size(XSCRS,2)), sin2(size(XSCRS,2)), sin3(size(XSCRS,2))
			real(dp) :: e_1(size(XSCRS,2),3), e_2(size(XSCRS,2),3), e_3(size(XSCRS,2),3)
			real(dp) :: M(3,3), E(3,3), dotME(3,3)
			
			! ----------------------------------------------------------
			! XSCRS is a 3xN or 6xN array of vectors in the SCRS axes
			! ----------------------------------------------------------

			dim = size(XSCRS,1)
			ntimes = size(XSCRS,2)

			! -------------------------------------------------------------------------------
			! The MER axes are defined differently for different ephemeris models
			! A 1-2-3 series of rotations transforms the MER axes to the lunar principal axes
			! -------------------------------------------------------------------------------
			
			! We use 430 model
			ang = [0.285D0, 78.580D0, 67.573D0]/3600.0D0*pi/180.0D0
			
			c1 = cos(ang(1))
			c2 = cos(ang(2))
			c3 = cos(ang(3))
			s1 = sin(ang(1))
			s2 = sin(ang(2))
			s3 = sin(ang(3))

			! 3-1-3 Euler angles describe librations of the lunar principal axes wrt the SCRS axes
			do i = 1,ntimes
                angles(i,:) = MoonLibration(JD(i))
            end do

			cos1 = cos(angles(:,1))
			cos2 = cos(angles(:,2))
			cos3 = cos(angles(:,3))
			sin1 = sin(angles(:,1))
			sin2 = sin(angles(:,2))
			sin3 = sin(angles(:,3))
			
			e_1(:,1) = cos1*cos3 - sin1*cos2*sin3
			e_1(:,2) = sin1*cos3 + cos1*cos2*sin3
			e_1(:,3) = sin2*sin3
			
			e_2(:,1) = -cos1*sin3 - sin1*cos2*cos3
			e_2(:,2) = -sin1*sin3 + cos1*cos2*cos3
			e_2(:,3) = sin2*cos3
			
			e_3(:,1) = sin1*sin2
			e_3(:,2) = -cos1*sin2
			e_3(:,3) = cos2
			
			M(1,:) = [c2*c3, -c2*s3, s2]
			M(2,:) = [s1*s2*c3+c1*s3, -s1*s2*s3+c1*c3, -s1*c2]
			M(3,:) = [-c1*s2*c3+s1*s3,  c1*s2*s3+s1*c3,  c1*c2]
 
			dXMER = 0.0D0
			XMER = 0.0D0
			do i = 1,ntimes
				E(1,:) = e_1(i,:)
				E(2,:) = e_2(i,:)
				E(3,:) = e_3(i,:)
				dotME = dotmm(M, E)
				if (dim .eq. 3) then
					dXMER(:,:,i) = dotME
				else if (dim .eq. 6) then
					dXMER(1:3,1:3,i) = dotME
					dXMER(4:6,4:6,i) = dotME
				end if
				XMER(:,i) = dotmv(dXMER(:,:,i), XSCRS(:,i))
			end do
		
		end subroutine kscrs2mer
		subroutine kmer2scrs(XMER, JD, XSCRS, dXSCRS)
		
			use Ephemeris
			use LinearAlgebraInterfaces
		
			implicit none
			
			real(dp), intent(in) :: XMER(:,:), JD(size(XMER,2))
			real(dp), intent(out) :: XSCRS(size(XMER,1),size(XMER,2))
			real(dp), intent(out) :: dXSCRS(size(XMER,1),size(XMER,1),size(XMER,2))
			integer :: dim, ntimes, i
			real(dp) :: c1, c2, c3, s1, s2, s3, angles(size(XMER,2),3), ang(3)
			real(dp) :: cos1(size(XMER,2)), cos2(size(XMER,2)), cos3(size(XMER,2))
			real(dp) :: sin1(size(XMER,2)), sin2(size(XMER,2)), sin3(size(XMER,2))
			real(dp) :: e_1(3,size(XMER,2)), e_2(3,size(XMER,2)), e_3(3,size(XMER,2))
			real(dp) :: M(3,3), E(3,3), dotEM(3,3)
			
			! --------------------------------------------------------
			! MER_vec is a 3xN or 6xN array of vectors in the MER axes
			! --------------------------------------------------------

			dim = size(XMER,1)
			ntimes = size(XMER,2)

			! -------------------------------------------------------------------------------
			! The MER axes are defined differently for different ephemeris models
			! A 1-2-3 series of rotations transforms the MER axes to the lunar principal axes
			! -------------------------------------------------------------------------------
			
			! We use 430 model.
			ang = [0.285D0, 78.580D0, 67.573D0]/3600.D0*pi/180.0D0

			c1 = cos(ang(1))
			c2 = cos(ang(2))
			c3 = cos(ang(3))
			s1 = sin(ang(1))
			s2 = sin(ang(2))
			s3 = sin(ang(3))

			! 3-1-3 Euler angles describe librations of the lunar principal axes wrt the SCRS axes
			do i = 1,ntimes
				angles(i,:) = MoonLibration(JD(i))
			end do

			cos1 = cos(angles(:,1))
			cos2 = cos(angles(:,2))
			cos3 = cos(angles(:,3))
			sin1 = sin(angles(:,1))
			sin2 = sin(angles(:,2))
			sin3 = sin(angles(:,3))
			
			e_1(1,:) = cos1*cos3 - sin1*cos2*sin3
			e_1(2,:) = sin1*cos3 + cos1*cos2*sin3
			e_1(3,:) = sin2*sin3
			
			e_2(1,:) = -cos1*sin3 - sin1*cos2*cos3
			e_2(2,:) = -sin1*sin3 + cos1*cos2*cos3
			e_2(3,:) = sin2*cos3
			
			e_3(1,:) = sin1*sin2
			e_3(2,:) = -cos1*sin2
			e_3(3,:) = cos2
			
			M(1,:) = [ c2*c3, s1*s2*c3+c1*s3, -c1*s2*c3+s1*s3 ]
			M(2,:) = [-c2*s3, -s1*s2*s3+c1*c3, c1*s2*s3+s1*c3 ]
			M(3,:) = [s2, -s1*c2, c1*c2]
 
			dXSCRS = 0.0D0
			XSCRS = 0.0D0
			do i = 1,ntimes
				E(:,1) = e_1(:,i)
				E(:,2) = e_2(:,i)
				E(:,3) = e_3(:,i)
				dotEM = dotmm(E, M)
				if (dim .eq. 3) then
					dXSCRS(:,:,i) = dotEM
				else if (dim == 6) then
					dXSCRS(1:3,1:3,i) = dotEM
					dXSCRS(4:6,4:6,i) = dotEM
				end if
				XSCRS(:,i) = dotmv(dXSCRS(:,:,i), XMER(:,i))
			end do
		
		end subroutine kmer2scrs
		subroutine kgcrs2scrs(XGCRS_DIM, JD, XSCRS_DIM)
		
			use Ephemeris
		
			implicit none
			
			real(dp), intent(in) :: XGCRS_DIM(:,:), JD(size(XGCRS_DIM,2))
			real(dp), intent(out) :: XSCRS_DIM(size(XGCRS_DIM,1),size(XGCRS_DIM,2))
			real(dp) :: xMoon(6)
			integer :: i
			
			do i = 1,size(XGCRS_DIM,2)
				xMoon = PlanetState(JD(i),'Moon','Earth')
				XSCRS_DIM(:,i) = XGCRS_DIM(:,i) - xMoon
			end do		
		
		end subroutine kgcrs2scrs
		subroutine kscrs2gcrs(XSCRS_DIM, JD, XGCRS_DIM)
		
			use Ephemeris
			
			implicit none
			
			real(dp), intent(in) :: XSCRS_DIM(:,:), JD(size(XSCRS_DIM,2))
			real(dp), intent(out) :: XGCRS_DIM(size(XSCRS_DIM,1),size(XSCRS_DIM,2))
			real(dp) :: xEarth(6)
			integer :: i
			
			do i = 1,size(XGCRS_DIM,2)
				xEarth = PlanetState(JD(i),'Earth','Moon')
				XGCRS_DIM(:,i) = XSCRS_DIM(:,i) - xEarth
			end do	
		
		end subroutine kscrs2gcrs
		subroutine kitrs2gcrs(XITRS, JD, XGCRS, dXGCRS)
        
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
        subroutine kgcrs2itrs(XGCRS, JD, XITRS, dXITRS)
		
        
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp), intent(in) :: XGCRS(3), JD
            real(dp), intent(out) :: XITRS(3), dXITRS(3,3)
            real(dp) :: XGCRS_EYE(3), dXGCRS_EYE(3,3)
            
            call kitrs2gcrs([1.0D0, 0.0D0, 0.0D0], JD, XGCRS_EYE, dXGCRS_EYE)
            dXITRS = transpose(dXGCRS_EYE)
            
            XITRS = dotmv(dXITRS,XGCRS)

        end subroutine kgcrs2itrs
		subroutine kscrs2sors(XSCRS, JD, XSORS, dXSORS)
		
			use BaseMeansToolbox
			use Ephemeris
		
			implicit none
			
			real(dp), intent(in) :: XSCRS(:,:), JD(size(XSCRS,2))
			real(dp), intent(out) :: XSORS(size(XSCRS,1),size(XSCRS,2))
			real(dp), intent(out) :: dXSORS(size(XSCRS,1),size(XSCRS,1),size(XSCRS,2))
			integer :: dim, ntimes, i
			real(dp) :: xEarth(6), RE(3), VE(3), CE(3)
			real(dp) :: E1(3), E2(3), E3(3)
		
			! ----------------------------------------------------------
			! XSCRS is a 3xN or 6xN array of vectors in the SCRS axes
			! ----------------------------------------------------------

			dim = size(XSCRS,1)
			ntimes = size(XSCRS,2)
			
			do i = 1,ntimes
			
				xEarth = PlanetState(JD(i),'Earth','Moon')
				RE = xEarth(1:3)/sqrt(sum(xEarth(1:3)*xEarth(1:3)))
				VE = xEarth(4:6)/sqrt(sum(xEarth(4:6)*xEarth(4:6)))
				CE = cross(RE,VE)
				
				E1 = RE
				E3 = CE/sqrt(sum(CE*CE))
				E2 = cross(E3, E1)
				
				dXSORS(1,1:3,i) = E1
				dXSORS(2,1:3,i) = E2
				dXSORS(3,1:3,i) = E3
				if (dim .eq. 6) then
					dXSORS(4,4:6,i) = E1
					dXSORS(5,4:6,i) = E2
					dXSORS(6,4:6,i) = E3
				end if
				
				XSORS(:,i) = dotmv(dXSORS(:,:,i), XSCRS(:,i))
				
			end do
		
		end subroutine kscrs2sors
		subroutine ksors2scrs(XSORS, JD, XSCRS, dXSCRS)
		
			use BaseMeansToolbox
			use Ephemeris
		
			implicit none
			
			real(dp), intent(in) :: XSORS(:,:), JD(size(XSORS,2))
			real(dp), intent(out) :: XSCRS(size(XSORS,1),size(XSORS,2))
			real(dp), intent(out) :: dXSCRS(size(XSORS,1),size(XSORS,1),size(XSORS,2))
			integer :: dim, ntimes, i
			real(dp) :: xEarth(6), RE(3), VE(3), CE(3)
			real(dp) :: E1(3), E2(3), E3(3)
		
			! ----------------------------------------------------------
			! XSORS is a 3xN or 6xN array of vectors in the SORS axes
			! ----------------------------------------------------------

			dim = size(XSCRS,1)
			ntimes = size(XSCRS,2)
			
			do i = 1,ntimes
			
				xEarth = PlanetState(JD(i),'Earth','Moon')
				RE = xEarth(1:3)/sqrt(sum(xEarth(1:3)*xEarth(1:3)))
				VE = xEarth(4:6)/sqrt(sum(xEarth(4:6)*xEarth(4:6)))
				CE = cross(RE,VE)
				
				E1 = RE
				E3 = CE/sqrt(sum(CE*CE))
				E2 = cross(E3, E1)
				
				dXSCRS(1:3,1,i) = E1
				dXSCRS(1:3,2,i) = E2
				dXSCRS(1:3,3,i) = E3
				if (dim .eq. 6) then
					dXSCRS(4:6,4,i) = E1
					dXSCRS(4:6,5,i) = E2
					dXSCRS(4:6,6,i) = E3
				end if
				
				XSCRS(:,i) = dotmv(dXSCRS(:,:,i), XSORS(:,i))
				
			end do
			
		
		end subroutine ksors2scrs
		
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