module PropagationModule
    
    use OdeToolbox
    use EquationsModule
    use Ephemeris
    use BaseMeansToolbox
    use MathAndProgConstants
    
    contains
    
        ! Opened for F2PY wrapping.
        subroutine propagate_nbp(central_body, tspan, x0, sources, dat, stm, vars, neq, T, Y)
    
            implicit none
            
            character(*), intent(in) :: central_body, vars
            real(dp), intent(in) :: tspan(:), x0(6)
            integer(2), intent(in) :: sources(14), stm
            real(dp), intent(in) :: dat(4)
            real(dp), allocatable :: s0(:)
            integer(2) :: neq
            type(OdeOptions) :: options
            type(OdeSol)     :: sol
            real(dp), intent(out) :: T(size(tspan)), Y(neq, size(tspan))
    
            if (central_body .eq. 'Moon') then
                DistUnit = MoonDistUnit
                VelUnit = MoonVelUnit
                TimeUnit = MoonTimeUnit
                AccUnit = MoonAccUnit
                RSun = MoonRSun
                REarth = MoonREarth
                RMoon = MoonRMoon
            else
                DistUnit = EarthDistUnit
                VelUnit = EarthVelUnit
                TimeUnit = EarthTimeUnit
                AccUnit = EarthAccUnit
                RSun = EarthRSun
                REarth = EarthREarth
                RMoon = EarthRMoon
            end if
    
            sun       = int2logical(sources(1))
            mercury   = int2logical(sources(2))
            venus     = int2logical(sources(3))
            earth     = int2logical(sources(4))
            moon      = int2logical(sources(5))
            mars      = int2logical(sources(6))
            jupiter   = int2logical(sources(7))
            saturn    = int2logical(sources(8))
            uranus    = int2logical(sources(9))
            neptune   = int2logical(sources(10))
            srp       = int2logical(sources(11))
            cmplxmoon = int2logical(sources(12))
            atm       = int2logical(sources(13))
            j2        = int2logical(sources(14))
            
            JD_Zero = dat(1)
            area    = dat(2)
            mass    = dat(3)
            order   = int(dat(4))
    
            stm_required = stm
    
            ! State the propagation problem.
            if (stm_required) then
                s0 = [ x0, reshape(eye(6), [36,1]) ]
            else
                s0 = x0
            end if
    
            if (central_body .eq. 'Moon') then
                if (vars .eq. 'rv') then
                    call ode113(knbp_rv_Moon, tspan, s0, options, T, Y)
                else
                    call ode113(knbp_ee_Moon, tspan, s0, options, T, Y)
                end if
            else
                if (vars .eq. 'rv') then
                    call ode113(knbp_rv_Earth, tspan, s0, options, T, Y)
                else
                    call ode113(knbp_ee_Earth, tspan, s0, options, T, Y)
                end if
            end if
    
        end subroutine propagate_nbp
        subroutine propagate_r2bp(tspan, x0, T, Y)
        
            implicit none
            
            real(dp), intent(in) :: tspan(:), x0(6)
            real(dp), intent(out) :: T(size(tspan)), Y(6, size(tspan))
            type(OdeOptions) :: options
            
            call ode113(kr2bp, tspan, x0, options, T, Y)

        end subroutine propagate_r2bp
        subroutine propagate_cr3bp(central_body, tspan, x0, mu, stm_required, neq, T, Y)
        
            implicit none
            
            character(*), intent(in) :: central_body
            real(dp), intent(in) :: tspan(:), x0(6), mu
            logical, intent(in) :: stm_required
            integer, intent(in) :: neq
            real(dp), intent(out) :: T(size(tspan)), Y(neq, size(tspan))
            real(dp), allocatable :: s0(:)
            type(OdeOptions) :: options
            
            MassParameter = mu
            
            ! State the propagation problem.
            if (stm_required) then
                s0 = [ x0, reshape(eye(6), [36,1]) ]
            else
                s0 = x0
            end if
            
            if (central_body == 'First') then
                call ode113(kcr3bp_fb, tspan, s0, options, T, Y)
            else if (central_body == 'Secondary') then
                call ode113(kcr3bp_sb, tspan, s0, options, T, Y)
            else if (central_body == 'Center') then
                call ode113(kcr3bp, tspan, s0, options, T, Y)
            else
                write(*,*) 'Unknown body.'
                return
            end if
        
        end subroutine propagate_cr3bp
        subroutine propagate_br4bp(central_body, tspan, x0, mu, GM4b, a4b, theta0, stm_required, neq, T, Y)
        
            implicit none
            
            character(*), intent(in) :: central_body
            real(dp), intent(in) :: tspan(:), x0(6), mu, GM4b, a4b, theta0
            logical, intent(in) :: stm_required
            integer, intent(in) :: neq
            real(dp), intent(out) :: T(size(tspan)), Y(neq, size(tspan))
            real(dp), allocatable :: s0(:)
            type(OdeOptions) :: options
            
            MassParameter = mu
	        GravParameterFoursBody = GM4b
	        DistanceToFoursBody = a4b
	        InitialSynodicPhase = theta0
            
            ! State the propagation problem.
            if (stm_required) then
                s0 = [ x0, reshape(eye(6), [36,1]) ]
            else
                s0 = x0
            end if
            
            if (central_body == 'First') then
                call ode113(kbr4bp_fb, tspan, s0, options, T, Y)
            else if (central_body == 'Secondary') then
                call ode113(kbr4bp_sb, tspan, s0, options, T, Y)
            else if (central_body == 'Center') then
                call ode113(kbr4bp, tspan, s0, options, T, Y)
            else
                write(*,*) 'Unknown body.'
                return
            end if
        
        end subroutine propagate_br4bp
        
        
end module PropagationModule