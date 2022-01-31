module PropagationModule
    
    use OdeToolbox
    use EquationsModule
    use Ephemeris
    use BaseMeansToolbox
    use MathAndProgConstants
    
    contains
    
        ! Opened for F2PY wrapping.
        subroutine PropagateTrajectory(central_body, t0, tf, x0, sources, dat, stm, vars, neq, T, Y)
    
            implicit none
            
            character(*), intent(in) :: central_body, vars
            real(dp), intent(in) :: t0, tf, x0(6)
            integer(2), intent(in) :: sources(14), stm
            real(dp), intent(in) :: dat(4)
            real(dp), allocatable :: s0(:)
            integer(2) :: neq
            type(OdeOptions) :: options
            type(OdeSol)     :: sol
            real(dp), intent(out) :: T(2), Y(neq,2)
    
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
                s0 = (/ x0, reshape(eye(6),(/36,1/)) /)
            else
                s0 = x0
            end if
    
            if (central_body .eq. 'Moon') then
                if (vars .eq. 'rv') then
                    sol = ode113(knbp_rv_Moon,(/t0,tf/),s0,options)
                else
                    sol = ode113(knbp_ee_Moon,(/t0,tf/),s0,options)
                end if
            else
                if (vars .eq. 'rv') then
                    sol = ode113(knbp_rv_Earth,(/t0,tf/),s0,options)
                else
                    sol = ode113(knbp_ee_Earth,(/t0,tf/),s0,options)
                end if
            end if
    
            T = sol.T
            Y = sol.Y
    
    end subroutine PropagateTrajectory
    
end module PropagationModule