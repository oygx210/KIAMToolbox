module EquationsModule
    
    use MathAndProgConstants
    
    ! Moon units.
    real(dp), parameter :: MoonMu = 4.902800145956161D+03
    real(dp), parameter :: MoonDistUnit = 1.737400000000000D+03
    real(dp), parameter :: MoonVelUnit = 1.679856509416075D+00
    real(dp), parameter :: MoonTimeUnit = 1.197054402181422D-02
    real(dp), parameter :: MoonAccUnit = 1.624218885822240D+00
    real(dp), parameter :: MoonRSun = 4.004259237941752D+02
    real(dp), parameter :: MoonREarth = 3.671081270864510D+00
    real(dp), parameter :: MoonRMoon = 1.0D0
    
    ! Earth units.
    real(dp), parameter :: EarthMu = 3.986004356000000D+05
    real(dp), parameter :: EarthDistUnit = 6.371008400000000D+03
    real(dp), parameter :: EarthVelUnit = 7.909787126714006D+00
    real(dp), parameter :: EarthTimeUnit = 9.322440916154166D-03
    real(dp), parameter :: EarthAccUnit = 9.820224438870715D+00
    real(dp), parameter :: EarthRSun = 1.091977841372804D+02
    real(dp), parameter :: EarthREarth = 1.001118849568618D+00
    real(dp), parameter :: EarthRMoon = 2.727040824494911D-01

    logical :: atm = .false.
    logical :: j2 = .false.
    logical :: srp = .false.
    logical :: sun = .false.
    logical :: mercury = .false.
    logical :: venus = .false.
    logical :: earth = .false.
    logical :: mars = .false.
    logical :: jupiter = .false.
    logical :: saturn = .false.
    logical :: uranus = .false.
    logical :: neptune = .false.
    logical :: moon = .false.
    logical :: cmplxmoon = .false.

    real(dp) :: JD_Zero
    integer :: order
    real(dp) :: area
    real(dp) :: mass

    real(dp) :: DistUnit
    real(dp) :: VelUnit
    real(dp) :: TimeUnit
    real(dp) :: AccUnit
    real(dp) :: RSun
    real(dp) :: REarth
    real(dp) :: RMoon

    type ForceStructure
        real(dp) :: Force(3)
        real(dp) :: dFdx(3,6)
        real(dp) :: dFdr(3,3)
        real(dp) :: dFdv(3,3)
    end type
    type IntensityStructure
        real(dp) :: inten
        real(dp) :: didr(3)
        character(1) :: whichEclipse
    end type
    
    logical :: stm_required
    
    contains
    
        ! About function, open for F2PY wrapping
        function equationsmodule_about() result(i)
            integer :: i
            i = 0
            write(*,*) 'EquationsModule: Version 1.0.'
            write(*,*) 'Author: Maksim Shirobokov.'
            write(*,*) 'Date: 19.01.2022.'
        end function equationsmodule_about
        
        ! Open for F2PY wrapping
        function knbp_rv_Moon(t,x) result(dxdt)
        
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp),               intent(in)   :: t
            real(dp), dimension(:), intent(in)   :: x
            real(dp), dimension(size(x))     :: dxdt
            real(dp) :: r, jd
            real(dp) :: rvect(3), vvect(3), Fcent(3)
            real(dp) :: Phi(6,6), dPhidt(6,6), dfdr(3,3), dfdx(6,6)
            real(dp) :: muEarth, muSun, muMercury, muVenus, muMars
            real(dp) :: muJupiter, muSaturn, muUranus, muNeptune
            type(ForceStructure) :: SunGravForce, MercuryGravForce
            type(ForceStructure) :: VenusGravForce, EarthGravForce
            type(ForceStructure) :: MarsGravForce, JupiterGravForce
            type(ForceStructure) :: SaturnGravForce, UranusGravForce
            type(ForceStructure) :: NeptuneGravForce, SRPForce
            type(ForceStructure) :: FCentCmplx
            real(dp), dimension(3) :: FSun, FMercury, FVenus, FEarth, FMars
            real(dp), dimension(3) :: FJupiter, FSaturn, FUranus, FNeptune, FSRP
            real(dp), dimension(3,3) :: dFSun, dFMercury, dFVenus, dFEarth, dFMars
            real(dp), dimension(3,3) :: dFJupiter, dFSaturn, dFUranus, dFNeptune, dFSRP, dFcent
            
            ! Moon system of units.
            ! mu = 1.0

            rvect = x(1:3)
            vvect = x(4:6)
            r = norm(rvect)
            if (stm_required) then
                Phi = reshape(x(7:42),(/6,6/))
            end if

            if (JD_Zero .gt. 0) then
                jd = t*TimeUnit + JD_Zero
            end if

            Fcent = -rvect/r**3
            if (.not. stm_required) then
                dxdt = (/ vvect,  Fcent /)
            else
                dxdt = (/ vvect, Fcent, zeros(36,1) /)
            end if
            
            if (stm_required) then
                dfdr = -eye(3)/r**3 + 3.0D0*doto(rvect,rvect)/r**5
            end if
            
            if (CmplxMoon) then
                FCentCmplx = cmplx_moon(rvect,jd,stm_required)
                Fcent = FCentCmplx.Force
                dFcent = FCentCmplx.dFdr
                dxdt(4:6) = Fcent
                if (stm_required) then
                    dfdr = dFcent
                end if
            end if
            
            
            
            if (SRP) then
                SRPForce = kSRPForce(jd,rvect,'Moon',stm_required)
                FSRP = SRPForce.Force
                dFSRP = SRPForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FSRP
                if (stm_required) then
                    dfdr = dfdr + dFSRP
                end if
            end if
            
            if (Earth) then
                muEarth = 8.130056778446628D+01
                EarthGravForce = kgravForce(muEarth,jd,rvect,'Moon','Earth',stm_required)
                FEarth = EarthGravForce.Force
                dFEarth = EarthGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FEarth
                if (stm_required) then
                    dfdr = dfdr + dFEarth
                end if
            end if
            
            if (Sun) then
                muSun = 2.706870280026843D+07
                SunGravForce = kgravForce(muSun,jd,rvect,'Moon','Sun',stm_required)
                FSun = SunGravForce.Force
                dFSun = SunGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FSun
                if (stm_required) then
                    dfdr = dfdr + dFSun
                end if
            end if

            if (Mercury) then
                muMercury = 4.493774951900596D+00
                MercuryGravForce = kgravForce(muMercury,jd,rvect,'Moon','Mercury',stm_required)
                FMercury = MercuryGravForce.Force
                dFMercury = MercuryGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FMercury
                if (stm_required) then
                    dfdr = dfdr + dFMercury
                end if
            end if

            if (Venus) then
                muVenus = 6.625980689754914D+01
                VenusGravForce = kgravForce(muVenus,jd,rvect,'Moon','Venus',stm_required)
                FVenus = VenusGravForce.Force
                dFVenus = VenusGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FVenus
                if (stm_required) then
                    dfdr = dfdr + dFVenus
                end if
            end if

            if (Mars) then
                muMars = 8.735492767886337D+00
                MarsGravForce = kgravForce(muMars,jd,rvect,'Moon','Mars',stm_required)
                FMars = MarsGravForce.Force
                dFMars = MarsGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FMars
                if (stm_required) then
                    dfdr = dfdr + dFMars
                end if
            end if

            if (Jupiter) then
                muJupiter = 2.584497813153080D+04
                JupiterGravForce = kgravForce(muJupiter,jd,rvect,'Moon','Jupiter',stm_required)
                FJupiter = JupiterGravForce.Force
                dFJupiter = JupiterGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FJupiter
                if (stm_required) then
                    dfdr = dfdr + dFJupiter
                end if
            end if

            if (Saturn) then
                muSaturn = 7.738554238506189D+03
                SaturnGravForce = kgravForce(muSaturn,jd,rvect,'Moon','Saturn',stm_required)
                FSaturn = SaturnGravForce.Force
                dFSaturn = SaturnGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FSaturn
                if (stm_required) then
                    dfdr = dfdr + dFSaturn
                end if
            end if

            if (Uranus) then
                muUranus = 1.181885623629258D+03
                UranusGravForce = kgravForce(muUranus,jd,rvect,'Moon','Uranus',stm_required)
                FUranus = UranusGravForce.Force
                dFUranus = UranusGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FUranus
                if (stm_required) then
                    dfdr = dfdr + dFUranus
                end if
            end if

            if (Neptune) then
                muNeptune = 1.394412747421909D+03
                NeptuneGravForce = kgravForce(muNeptune,jd,rvect,'Moon','Neptune',stm_required)
                FNeptune = NeptuneGravForce.Force
                dFNeptune = NeptuneGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FNeptune
                if (stm_required) then
                    dfdr = dfdr + dFNeptune
                end if
            end if

            if (stm_required) then
                dfdx = 0.0D0
                dfdx(1:3,4:6) = eye(3)
                dfdx(4:6,1:3) = dfdr
                dPhidt = dotmm(dfdx, Phi)
                dxdt(7:42) = reshape(dPhidt,(/36/))
            end if
        
        end function knbp_rv_Moon
        
        ! Open for F2PY wrapping
        function knbp_rv_Earth(t,x) result(dxdt)
        
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp),               intent(in)   :: t
            real(dp), dimension(:), intent(in)   :: x
            real(dp), dimension(size(x))     :: dxdt
            real(dp) :: r, jd
            real(dp) :: rvect(3), vvect(3), Fcent(3)
            real(dp) :: Phi(6,6), dPhidt(6,6), dfdx(6,6), dFatm(3,6)
            real(dp) :: muMoon, muSun, muMercury, muVenus, muMars
            real(dp) :: muJupiter, muSaturn, muUranus, muNeptune
            type(ForceStructure) :: SunGravForce, MercuryGravForce
            type(ForceStructure) :: VenusGravForce, MoonGravForce
            type(ForceStructure) :: MarsGravForce, JupiterGravForce
            type(ForceStructure) :: SaturnGravForce, UranusGravForce
            type(ForceStructure) :: NeptuneGravForce, SRPForce, AtmForce, J2Force
            real(dp), dimension(3) :: FSun, FMercury, FVenus, FMoon, FMars, Fatm, FJ2
            real(dp), dimension(3) :: FJupiter, FSaturn, FUranus, FNeptune, FSRP
            real(dp), dimension(3,3) :: dFSun, dFMercury, dFVenus, dFMoon, dFMars, dFJ2
            real(dp), dimension(3,3) :: dFJupiter, dFSaturn, dFUranus, dFNeptune, dFSRP, dFcentdr
            
            ! Earth system of units.
            ! mu = 1.0

            rvect = x(1:3)
            vvect = x(4:6)
            r = norm(rvect)
            if (stm_required) then
                Phi = reshape(x(7:42),(/6,6/))
            end if

            if (JD_Zero .gt. 0) then
                jd = t*TimeUnit + JD_Zero
            end if

            Fcent = -rvect/r**3
            if (.not. stm_required) then
                dxdt = (/ vvect,  Fcent /)
            else
                dxdt = (/ vvect, Fcent, zeros(36,1) /)
            end if
            
            dfdx = 0.0D0
            
            if (stm_required) then
                dfdx(4:6,1:3) = -eye(3)/r**3 + 3.0D0*doto(rvect,rvect)/r**5
            end if
            
            if (atm) then
                AtmForce = kEarthAtm(x(1:6),stm_required)
                Fatm = AtmForce.Force
                dFatm = AtmForce.dFdx
                dxdt(4:6) = dxdt(4:6) + Fatm
                if (stm_required) then
                    dfdx(4:6,:) = dfdx(4:6,:) + dFatm
                end if
            end if
            
            if (J2) then
                J2Force = kEarthJ2(jd,rvect,stm_required)
                FJ2 = J2Force.Force
                dFJ2 = J2Force.dFdr
                dxdt(4:6) = dxdt(4:6) + FJ2
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFJ2
                end if
            end if
            
            if (SRP) then
                SRPForce = kSRPForce(jd,rvect,'Earth',stm_required)
                FSRP = SRPForce.Force
                dFSRP = SRPForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FSRP
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFSRP
                end if
            end if
            
            if (Moon) then
                muMoon = 1.230003710000000D-02
                MoonGravForce = kgravForce(muMoon,jd,rvect,'Earth','Moon',stm_required)
                FMoon = MoonGravForce.Force
                dFMoon = MoonGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FMoon
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFMoon
                end if
            end if
            
            if (Sun) then
                muSun = 3.329460486921756D+05
                SunGravForce = kgravForce(muSun,jd,rvect,'Earth','Sun',stm_required)
                FSun = SunGravForce.Force
                dFSun = SunGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FSun
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFSun
                end if
            end if

            if (Mercury) then
                muMercury = 5.527359862742804D-02
                MercuryGravForce = kgravForce(muMercury,jd,rvect,'Earth','Mercury',stm_required)
                FMercury = MercuryGravForce.Force
                dFMercury = MercuryGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FMercury
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFMercury
                end if
            end if

            if (Venus) then
                muVenus = 8.149980830786903D-01
                VenusGravForce = kgravForce(muVenus,jd,rvect,'Earth','Venus',stm_required)
                FVenus = VenusGravForce.Force
                dFVenus = VenusGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FVenus
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFVenus
                end if
            end if

            if (Mars) then
                muMars = 1.074468851317836D-01
                MarsGravForce = kgravForce(muMars,jd,rvect,'Earth','Mars',stm_required)
                FMars = MarsGravForce.Force
                dFMars = MarsGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FMars
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFMars
                end if
            end if

            if (Jupiter) then
                muJupiter = 3.178941898665174D+02
                JupiterGravForce = kgravForce(muJupiter,jd,rvect,'Earth','Jupiter',stm_required)
                FJupiter = JupiterGravForce.Force
                dFJupiter = JupiterGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FJupiter
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFJupiter
                end if
            end if

            if (Saturn) then
                muSaturn = 9.518450423398838D+01
                SaturnGravForce = kgravForce(muSaturn,jd,rvect,'Earth','Saturn',stm_required)
                FSaturn = SaturnGravForce.Force
                dFSaturn = SaturnGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FSaturn
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFSaturn
                end if
            end if

            if (Uranus) then
                muUranus = 1.453723701859651D+01
                UranusGravForce = kgravForce(muUranus,jd,rvect,'Earth','Uranus',stm_required)
                FUranus = UranusGravForce.Force
                dFUranus = UranusGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FUranus
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFUranus
                end if
            end if

            if (Neptune) then
                muNeptune = 1.715132852600241D+01
                NeptuneGravForce = kgravForce(muNeptune,jd,rvect,'Earth','Neptune',stm_required)
                FNeptune = NeptuneGravForce.Force
                dFNeptune = NeptuneGravForce.dFdr
                dxdt(4:6) = dxdt(4:6) + FNeptune
                if (stm_required) then
                    dfdx(4:6,1:3) = dfdx(4:6,1:3) + dFNeptune
                end if
            end if

            if (stm_required) then
                dfdx(1:3,4:6) = eye(3)
                dPhidt = dotmm(dfdx, Phi)
                dxdt(7:42) = reshape(dPhidt,(/36/))
            end if
        
        end function knbp_rv_Earth    
        
        ! Open for F2PY wrapping
        function knbp_ee_Moon(t,x) result(dxdt)
        
            use Translations
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            ! No gradients.
        
            ! Moon system of units.
            ! mu = 1.0
            
            implicit none
            
            real(dp),               intent(in)   :: t
            real(dp), dimension(:), intent(in)   :: x
            real(dp), dimension(size(x))         :: dxdt
            real(dp) :: rv(6), rvect(3), vvect(3), cvect(3)
            real(dp) :: drvdt(6), Fcent(3), F(3)
            real(dp) :: h, ex, ey, ix, iy, L, E1(3), E2(3), E3(3)
            real(dp) :: trans_mat(3,3), Forb(3)
            real(dp) :: sinL, cosL, phi, eta, sigma, M(6,3)

            call kee2rv(x, 1.0D0, rv)
            rvect = rv(1:3)
            vvect = rv(4:6)
            cvect = cross(rvect,vvect)

            drvdt = knbp_rv_Moon(t,rv)
            Fcent = -rvect/norm(rvect)**3
            F = drvdt(4:6) - Fcent

            h = x(1)
            ex = x(2)
            ey = x(3)
            ix = x(4)
            iy = x(5)
            L = x(6)

            E1 = rvect/norm(rvect)
            E3 = cvect/norm(cvect)
            E2 = cross(E3,E1)
            
            trans_mat(1,:) = E1
            trans_mat(2,:) = E2
            trans_mat(3,:) = E3
            
            Forb = dotmv(trans_mat,F)

            sinL = sin(L)
            cosL = cos(L)

            phi = 1 + ix**2 + iy**2
            eta = ix*sinL - iy*cosL
            sigma = 1 + ex*cosL + ey*sinL

            M(1,:) = (/ 0.0D0, h**2/sigma, 0.0D0 /)
            M(2,:) = h* (/ sinL, (1.0D0+1.0D0/sigma)*cosL+ex/sigma, -ey*eta/sigma /)
            M(3,:) = h* (/-cosL, (1.0D0+1.0D0/sigma)*sinL+ey/sigma,  ex*eta/sigma /)
            M(4,:) = (/ 0.0D0, 0.0D0, h*phi*cosL/2.0D0/sigma /)
            M(5,:) = (/ 0.0D0, 0.0D0, h*phi*sinL/2.0D0/sigma /)
            M(6,:) = (/ 0.0D0, 0.0D0, h*eta/sigma /)

            dxdt = dotmv(M,Forb) + (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, sigma**2/h**3 /)
        
        end function knbp_ee_Moon
        
        ! Open for F2PY wrapping
        function knbp_ee_Earth(t,x) result(dxdt)
        
            use Translations
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            ! No gradients.
        
            ! Earth system of units.
            ! mu = 1.0
            
            implicit none
            
            real(dp),               intent(in)   :: t
            real(dp), dimension(:), intent(in)   :: x
            real(dp), dimension(size(x))     :: dxdt
            real(dp) :: rv(6), rvect(3), vvect(3), cvect(3)
            real(dp) :: drvdt(6), Fcent(3), F(3)
            real(dp) :: h, ex, ey, ix, iy, L, E1(3), E2(3), E3(3)
            real(dp) :: trans_mat(3,3), Forb(3)
            real(dp) :: sinL, cosL, phi, eta, sigma, M(6,3)

            call kee2rv(x, 1.0D0, rv)
            rvect = rv(1:3)
            vvect = rv(4:6)
            cvect = cross(rvect,vvect)

            drvdt = knbp_rv_Earth(t,rv)
            Fcent = -rvect/norm(rvect)**3
            F = drvdt(4:6) - Fcent

            h = x(1)
            ex = x(2)
            ey = x(3)
            ix = x(4)
            iy = x(5)
            L = x(6)

            E1 = rvect/norm(rvect)
            E3 = cvect/norm(cvect)
            E2 = cross(E3,E1)
            
            trans_mat(1,:) = E1
            trans_mat(2,:) = E2
            trans_mat(3,:) = E3
            
            Forb = dotmv(trans_mat,F)

            sinL = sin(L)
            cosL = cos(L)

            phi = 1 + ix**2 + iy**2
            eta = ix*sinL - iy*cosL
            sigma = 1 + ex*cosL + ey*sinL

            M(1,:) = (/ 0.0D0, h**2/sigma, 0.0D0 /)
            M(2,:) = h* (/ sinL, (1.0D0+1.0D0/sigma)*cosL+ex/sigma, -ey*eta/sigma /)
            M(3,:) = h* (/-cosL, (1.0D0+1.0D0/sigma)*sinL+ey/sigma,  ex*eta/sigma /)
            M(4,:) = (/ 0.0D0, 0.0D0, h*phi*cosL/2.0D0/sigma /)
            M(5,:) = (/ 0.0D0, 0.0D0, h*phi*sinL/2.0D0/sigma /)
            M(6,:) = (/ 0.0D0, 0.0D0, h*eta/sigma /)

            dxdt = dotmv(M,Forb) + (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, sigma**2/h**3 /)
        
        end function knbp_ee_Earth
        
        function cmplx_moon(rvect,jd,grad_req) result(FCentCmplx)
        
            use Translations
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp) :: rvect(3), jd
            type(ForceStructure) :: FCentCmplx, GravityMoonForce
            real(dp) :: rvectPA(3), M_SCRS_2_PA(3,3), Fcent_dim_pa(3), dFcent_dim_pa_d_rvectPA(3,3)
            real(dp) :: Fcent(3), dFcent_d_rvectPA(3,3), dFcent(3,3), Fcent_pa(3), dFcent_pa_d_rvectPA(3,3)
            logical :: grad_req
        
            call kscrs2pa(rvect*DistUnit, jd, rvectPA, M_SCRS_2_PA)
            GravityMoonForce = kGravityMoonGradient50(rvectPA,order,order,grad_req)
            Fcent_dim_pa = GravityMoonForce.Force
            Fcent_pa = Fcent_dim_pa/AccUnit
            Fcent = dotmv(transpose(M_SCRS_2_PA),Fcent_pa)
            
            dFcent = 0.0D0
            if (grad_req) then
                dFcent_dim_pa_d_rvectPA = GravityMoonForce.dFdr
                dFcent_pa_d_rvectPA = dFcent_dim_pa_d_rvectPA/AccUnit
                dFcent_d_rvectPA = dotmm(transpose(M_SCRS_2_PA),dFcent_pa_d_rvectPA)
                dFcent = dotmm(dFcent_d_rvectPA,M_SCRS_2_PA) * DistUnit
            end if
            
            FCentCmplx.Force = Fcent
            FCentCmplx.dFdr = dFcent
            
        end function cmplx_moon
        function kgravForce(muP,jd,rsc,CenterPlanet,TargetPlanet,grad_req) result(GravForce)
        
            use Ephemeris
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp) :: muP, jd, rsc(3), ps(6), rpl(3), dr(3), drnorm
            real(dp) :: drc(3), drcnorm, sc2pl(3), c2pl(3), F(3), dFdr(3,3)
            character(*) :: CenterPlanet, TargetPlanet
            logical :: grad_req
            type(ForceStructure) :: GravForce
            
            ps = PlanetState(jd,TargetPlanet,CenterPlanet)
            rpl = ps(1:3)/DistUnit

            dr = rsc - rpl;
            drnorm = norm(dr)

            drc = - rpl
            drcnorm = norm(drc)

            sc2pl = -muP*dr/drnorm**3
            c2pl = -muP*drc/drcnorm**3

            F = sc2pl - c2pl

            dFdr = 0.0D0
            if (grad_req) then
                dFdr = -muP*eye(3)/drnorm**3 + 3*muP*doto(dr,dr)/drnorm**5;
            end if
            
            GravForce.Force = F
            GravForce.dFdr = dFdr
        
        end function kgravForce
        function kSRPForce(jd,rsc,CentralBody,grad_req) result(SRPForce)
        
            use BaseMeansToolbox
            use Ephemeris
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp) :: sSunEarth(6), sSunMoon(6), rSunEarth(3), rSunMoon(3)
            real(dp) :: rsc(3), rMoonSc(3), rSunSc(3), rMoonEarth(3), rMoonSun(3), rEarthSc(3), F(3), dFdr(3,3)
            real(dp) :: jd, PAdm, Coeff, drnormSun, inten, didr(3)
            logical :: grad_req
            character(*) :: CentralBody
            type(ForceStructure) :: SRPForce
            type(IntensityStructure) :: Intensity

            sSunEarth = planetState(jd,'Earth','Sun')
            sSunMoon = planetState(jd,'Moon','Sun')
            
            rSunEarth = sSunEarth(1:3)/DistUnit
            rSunMoon = sSunMoon(1:3)/DistUnit
    
            if (CentralBody == 'Moon') then
                rMoonSc = rsc
                rSunSc = rSunMoon + rMoonSc
                rMoonEarth = rSunEarth - rSunMoon
                rMoonSun = - rSunMoon
            elseif (CentralBody == 'Earth') then
                rEarthSc = rsc
                rMoonEarth = rSunEarth - rSunMoon
                rMoonSc = rMoonEarth + rEarthSc
                rMoonSun = - rSunMoon
                rSunSc = rSunMoon + rMoonSc
            end if
    
            drnormSun = norm(rSunSc)

            ! Area to mass ratio (0.005--0.02 for standard s/c without solar sail or big panels)
            PAdm = area/mass; ! m^2/kg
            Coeff = 4.56D-06*PAdm/AccUnit*(1.495978707000000D+08/DistUnit)**2
            ! N/m^2 * m^2/kg / (m/s^2)

            Intensity = kintensity(rMoonSc,rMoonSun,rMoonEarth,grad_req)
            inten = Intensity.inten
            didr = Intensity.didr
            
            F = inten*Coeff*rSunSc/drnormSun**3

            dFdr = 0.0D0
            if (grad_req) then
                dFdr = inten*(Coeff*eye(3)/drnormSun**3-3*Coeff*doto(rSunSc,rSunSc)/drnormSun**5) + Coeff*doto(rSunSc,didr)/drnormSun**3
            end if
            
            SRPForce.Force = F
            SRPForce.dFdr = dFdr

        end function kSRPForce
        function kintensity(rsc,rvectSun,rvectEarth,grad_req) result(Intensity)
        
            use BaseMeansToolbox
            use LinearAlgebraInterfaces

            implicit none
            
            real(dp) :: rsc(3), rvectSun(3), rvectEarth(3)
            real(dp) :: value
            logical :: grad_req
            type(IntensityStructure) :: Intensity
            real(dp) :: CM, CE, DS, DM, DE, rhoS, rhoM, rhoE, cosTM, thetaM, cosTE, thetaE
            real(dp) :: cosS, sinS, cosM, sinM, cosE, sinE, sinTM, sinTE
            real(dp) :: cos1, cos2, cos3, acos1, acos2, acos3, m, n
            real(dp) :: dacos1, dacos2, dacos3
            real(dp), dimension(3) :: didr, drhoSdr, nE, nM, drhoMdr, drhoEdr, n1, n2M, n2E
            real(dp), dimension(3) :: dthetaEdr, dthetaMdr
            real(dp), dimension(3) :: dcSdr, dsSdr, dcMdr, dsMdr, dcEdr, dsEdr, dsTMdr, dsTEdr
            real(dp), dimension(3) :: dacos1dr, dacos2dr, dacos3dr, dmdr, dndr, dcTMdr, dcTEdr
            real(dp), dimension(3,3) :: dn1dr, dn2Mdr, dn2Edr
            character(1) :: whichEclipse
            logical :: PartialEclipseByTheMoon, TotalEclipseByTheMoon, AnnularEclipseByTheMoon
            logical :: PartialEclipseByTheEarth, TotalEclipseByTheEarth, AnnularEclipseByTheEarth
            

            if (.not. grad_req) then
                didr = 0.0D0
            end if

            CM = RMoon/(RSun - RMoon)
            CE = REarth/(RSun - REarth)

            ! ���������� �� �� �� ������ ������, ������ ���� � ������ �����
            DS = norm(rvectSun - rsc)
            DM = norm(rsc)
            DE = norm(rvectEarth - rsc)

            ! ������� ������� ������ ������
            if (RSun/DS < 1) then
                rhoS = asin(RSun/DS)
                if (grad_req) then
                    drhoSdr  = 1/sqrt(1-RSun**2/DS**2) * RSun *   1  / DS**3 * (rvectSun - rsc)
                end if
            else
                ! Spacecraft is under the surface of the Sun.
                value = 1.0D0
                if (grad_req) then
                    didr = 0.0D0
                end if
                whichEclipse = 'L'
                nE = 0.0D0
                nM = 0.0D0
                Intensity.inten = value
                Intensity.didr = didr
                Intensity.whichEclipse = whichEclipse
                return
            end if

            ! ������� ������� ������ ����
            if (RMoon/DM < 1) then
                rhoM = asin(RMoon/DM)
                if (grad_req) then
                    drhoMdr = 1/sqrt(1-RMoon**2/DM**2) * RMoon * (-1) / DM**3 * rsc
                end if
            else
                ! Spacecraft is under the surface of the Moon.
                value = 0.0D0
                if (grad_req) then
                    didr = 0.0D0
                end if
                whichEclipse = 'M'
                nE = 0.0D0
                nM = 0.0D0
                Intensity.inten = value
                Intensity.didr = didr
                Intensity.whichEclipse = whichEclipse
                return
            end if

            ! ������� ������� ������ �����
            if (REarth/DE < 1) then
                rhoE = asin(REarth/DE)
                if (grad_req) then
                    drhoEdr = 1/sqrt(1-REarth**2/DE**2) * REarth * 1 / DE**3 * (rvectEarth - rsc)
                end if
            else
                ! Spacecraft is under the surface of the Earth.
                value = 0.0D0
                if (grad_req) then
                    didr = 0.0D0
                end if
                whichEclipse = 'E'
                nE = 0.0D0
                nM = 0.0D0
                Intensity.inten = value
                Intensity.didr = didr
                Intensity.whichEclipse = whichEclipse
                return
            end if

            ! ������� ���������� ����� �������� ������ � ����
            n1 = (rvectSun - rsc)/DS
            n2M =      - rsc /DM
            cosTM = max(-1.0,min(1.0,sum(n1*n2M)))
            thetaM = acos(cosTM)

            ! ������� ���������� ����� �������� ������ � �����
            n2E = (rvectEarth - rsc)/DE
            cosTE = max(-1.0,min(1.0,sum(n1*n2E)))
            thetaE = acos(cosTE)

            if (grad_req) then
                dn1dr = -eye(3)/DS + (doto(rvectSun-rsc,rvectSun-rsc))/DS**3
                dn2Mdr = -eye(3)/DM + (doto(rsc,rsc))/DM**3
                dn2Edr = -eye(3)/DE + (doto(rvectEarth-rsc,rvectEarth-rsc))/DE**3
                if (abs(cosTM) < 1) then
                    dthetaMdr = -1/sqrt(1-cosTM**2) * (dotvm(n2M,dn1dr) + dotvm(n1,dn2Mdr))
                else
                    dthetaMdr = 0.0D0
                end if
                if (abs(cosTE) < 1) then
                    dthetaEdr = -1/sqrt(1-cosTE**2) * (dotvm(n2E,dn1dr) + dotvm(n1,dn2Edr))
                else
                    dthetaEdr = 0.0D0
                end if
            end if

            cosS = cos(rhoS)
            sinS = sin(rhoS)
            cosM = cos(rhoM)
            sinM = sin(rhoM)
            cosE = cos(rhoE)
            sinE = sin(rhoE)
            sinTM = sin(thetaM)
            sinTE = sin(thetaE)

            if (grad_req) then
                dcSdr = -sinS * drhoSdr
                dsSdr =  cosS * drhoSdr
                dcMdr = -sinM * drhoMdr
                dsMdr =  cosM * drhoMdr
                dcEdr = -sinE * drhoEdr
                dsEdr =  cosE * drhoEdr
                dcTMdr = -sinTM * dthetaMdr
                dsTMdr =  cosTM * dthetaMdr
                dcTEdr = -sinTE * dthetaEdr
                dsTEdr =  cosTE * dthetaEdr
            end if

            PartialEclipseByTheMoon = (DS > norm(rvectSun)) .and. (rhoS + rhoM > thetaM) .and. (thetaM > abs(rhoM - rhoS))
            TotalEclipseByTheMoon   = (DS > norm(rvectSun)) .and. (DS < (1 + CM)*norm(rvectSun)) .and. (rhoM - rhoS > thetaM)
            AnnularEclipseByTheMoon = (DS > (1 + CM)*norm(rvectSun)) .and. (rhoS - rhoM > thetaM)

            PartialEclipseByTheEarth = (DS > norm(rvectSun-rvectEarth)) .and. (rhoS + rhoE > thetaE) .and. (thetaE > abs(rhoE - rhoS))
            TotalEclipseByTheEarth   = (DS > norm(rvectSun-rvectEarth)) .and. (DS < (1 + CE)*norm(rvectSun-rvectEarth)) .and. (rhoE - rhoS > thetaE)
            AnnularEclipseByTheEarth = (DS > (1 + CE)*norm(rvectSun-rvectEarth)) .and. (rhoS - rhoE > thetaE)

            ! MATLAB CODE
            ! if any([PartialEclipseByTheMoon,TotalEclipseByTheMoon,AnnularEclipseByTheMoon]) && ...
            !    any([PartialEclipseByTheEarth,TotalEclipseByTheEarth,AnnularEclipseByTheEarth])
            !     fprintf('EclipseByTheMoon AND EclipseByTheEarth!');
            ! end

            ! �������� � �������������, ��� �������� ���� ��������� �����, ���� ������
            if (PartialEclipseByTheMoon) then
    
                !MATLAB CODE
                !fprintf('%6.6E  %6.6E  %6.6E PartialEclipseByTheMoon\n',asind(RSun/DS),asind(RMoon/DM),asind(REarth/DE));
    
                cos1 = max(-1.0,min(1.0,(cosM-cosS*cosTM)/sinS/sinTM))
                cos2 = max(-1.0,min(1.0,(cosS-cosM*cosTM)/sinM/sinTM))
                cos3 = max(-1.0,min(1.0,(cosTM-cosS*cosM)/sinS/sinM))

                acos1 = acos(cos1)
                acos2 = acos(cos2)
                acos3 = acos(cos3)
    
                m = pi - cosS * acos1 - cosM * acos2 - acos3
                n = 1 - cosS
    
                value = 1 - m/n/pi
    
                if (grad_req) then
        
                    if (abs(cos1) < 1) then
                        dacos1 = -1/sqrt(1-cos1**2)
                        dacos1dr = dacos1 * (dcMdr/sinS/sinTM - cosM/sinS**2/sinTM*dsSdr - cosM/sinTM**2/sinS*dsTMdr + 1/sinS**2*cotan(thetaM)*drhoSdr + 1/sinTM**2*cotan(rhoS)*dthetaMdr)
                    else
                        dacos1dr = 0.0D0
                    end if
        
                    if (abs(cos2) < 1) then
                        dacos2 = -1/sqrt(1-cos2**2)
                        dacos2dr = dacos2 * (dcSdr/sinM/sinTM - cosS/sinM**2/sinTM*dsMdr - cosS/sinTM**2/sinM*dsTMdr + 1/sinM**2*cotan(thetaM)*drhoMdr + 1/sinTM**2*cotan(rhoM)*dthetaMdr)
                    else
                        dacos2dr = 0.0D0
                    end if
        
                    if (abs(cos3) < 1) then
                        dacos3 = -1/sqrt(1-cos3**2)
                        dacos3dr = dacos3 * (dcTMdr/sinS/sinM - cosTM/sinS**2/sinM*dsSdr - cosTM/sinM**2/sinS*dsMdr  + 1/sinS**2*cotan(rhoM)  *drhoSdr + 1/sinM**2 *cotan(rhoS)*drhoMdr)
                    else
                        dacos3dr = 0.0D0
                    end if

                    dmdr = - dcSdr*acos1 - cosS*dacos1dr - dcMdr*acos2 - cosM*dacos2dr - dacos3dr
                    dndr = - dcSdr

                    didr  = -(dmdr*n  - dndr*m)/n**2/pi
    
                end if

            elseif (TotalEclipseByTheMoon) then
                
                !MATLAB CODE
                !fprintf('%6.6E  %6.6E  %6.6E TotalEclipseByTheMoon\n',asind(RSun/DS),asind(RMoon/DM),asind(REarth/DE));

                value = 0.0D0
    
                if (grad_req) then
                    didr  = 0.0D0
                end if
    
            elseif (AnnularEclipseByTheMoon) then
                
                !MATLAB CODE
                !fprintf('%6.6E  %6.6E  %6.6E AnnularEclipseByTheMoon\n',asind(RSun/DS),asind(RMoon/DM),asind(REarth/DE));
    
                value = 1 - (1-cosM)/(1-cosS)
    
                if (grad_req) then
                    didr =  - (-dcMdr* (1-cosS) + dcSdr* (1-cosM))/(1-cosS)**2
                end if
    
            elseif (PartialEclipseByTheEarth) then
                
                !MATLAB CODE
                !fprintf('%6.6E  %6.6E  %6.6E PartialEclipseByTheEarth\n',asind(RSun/DS),asind(RMoon/DM),asind(REarth/DE));
    
                cos1 = max(-1.0,min(1.0,(cosE-cosS*cosTE)/sinS/sinTE))
                cos2 = max(-1.0,min(1.0,(cosS-cosE*cosTE)/sinE/sinTE))
                cos3 = max(-1.0,min(1.0,(cosTE-cosS*cosE)/sinS/sinE))
    
                acos1 = acos(cos1)
                acos2 = acos(cos2)
                acos3 = acos(cos3)
    
                m = pi - cosS * acos1 - cosE * acos2 - acos3
                n = 1 - cosS
    
                value = 1 - m/n/pi
    
                if (grad_req) then
    
                    if (abs(cos1) < 1) then
                        dacos1 = -1/sqrt(1-cos1**2)
                        dacos1dr = dacos1 * (dcEdr/sinS/sinTE - cosE/sinS**2/sinTE*dsSdr - cosE/sinTE**2/sinS*dsTEdr + 1/sinS**2*cotan(thetaE)*drhoSdr + 1/sinTE**2*cotan(rhoS)*dthetaEdr)
                    else
                        dacos1dr = 0.0D0
                    end if
        
                    if (abs(cos2) < 1) then
                        dacos2 = -1/sqrt(1-cos2**2)
                        dacos2dr = dacos2 * (dcSdr/sinE/sinTE - cosS/sinE**2/sinTE*dsEdr - cosS/sinTE**2/sinE*dsTEdr + 1/sinE**2*cotan(thetaE)*drhoEdr + 1/sinTE**2*cotan(rhoE)*dthetaEdr)
                    else
                        dacos2dr = 0.0D0
                    end if
        
                    if (abs(cos3) < 1) then
                        dacos3 = -1/sqrt(1-cos3**2)
                        dacos3dr = dacos3 * (dcTEdr/sinS/sinE - cosTE/sinS**2/sinE*dsSdr - cosTE/sinE**2/sinS*dsEdr  + 1/sinS**2*cotan(rhoE)  *drhoSdr + 1/sinE**2 *cotan(rhoS)*drhoEdr)
                    else
                        dacos3dr = 0.0D0
                    end if

                    dmdr = - dcSdr*acos1 - cosS*dacos1dr - dcEdr*acos2 - cosE*dacos2dr - dacos3dr
                    dndr = - dcSdr

                    didr  = -(dmdr*n  - dndr*m)/n**2/pi
    
                end if
    
            elseif (TotalEclipseByTheEarth) then
                
                !MATLAB CODE
                !fprintf('%6.6E  %6.6E  %6.6E TotalEclipseByTheEarth\n',asind(RSun/DS),asind(RMoon/DM),asind(REarth/DE));

                value = 0.0D0
    
                if (grad_req) then
                    didr = 0.0D0
                end if
    
            elseif (AnnularEclipseByTheEarth) then
    
                !MATLAB CODE
                !fprintf('%6.6E  %6.6E  %6.6E AnnularEclipseByTheEarth\n',asind(RSun/DS),asind(RMoon/DM),asind(REarth/DE));
    
                value = 1 - (1-cosE)/(1-cosS)
    
                if (grad_req) then
                    didr =  - (-dcEdr* (1-cosS) + dcSdr* (1-cosE))/(1-cosS)**2
                end if
    
            else
                
                ! No eclipse
    
                value = 1.0D0
    
                if (grad_req) then
                    didr  = 0.0D0
                end if
    
            end if

            if (PartialEclipseByTheEarth .or. TotalEclipseByTheEarth .or. AnnularEclipseByTheEarth) then
                whichEclipse = 'E'
            elseif (PartialEclipseByTheMoon .or. TotalEclipseByTheMoon .or. AnnularEclipseByTheMoon) then
                whichEclipse = 'M'
            else
                whichEclipse = 'L'
            end if
            
            Intensity.inten = value
            Intensity.didr = didr
            Intensity.whichEclipse = whichEclipse

        end function kintensity
        function kGravityMoonGradient50(r_vec,N,M,grad_req) result(ComplexMoonForce)
        
            use GravityMoonCoefficients50
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp) :: ReferenceRadius, mu, a0
            integer :: N, M, L, ii, ji
            real(dp) :: i, j, temp_number
            real(dp) :: A(N+1,min(M+2,N+1)), B(N,min(M+2,N)), d(M+1)
            real(dp) :: F(N-1,M+1), G(N-1,M), H(N-1,M+1)
            real(dp) :: r_vec(3), rnorm, dr(3)
            real(dp) :: V(N+2,M+2), W(N+2,M+2)
            real(dp) :: dV(N+2,M+2,3), dW(N+2,M+2,3)
            real(dp) :: Qa(3), dQa(3,3), rel_harmonics(3,N-1), drh(3,N-1,3)
            real(dp) :: C(N*(N+3)/2,6)
            real(dp) :: QR, QX, QY, QZ, dQR(3), dQX(3), dQY(3), dQZ(3)
            type(ForceStructure) :: ComplexMoonForce
            real(dp) :: a_grav(3), a_cent(3), da_grav(3,3)
            logical :: grad_req
    
            ReferenceRadius = 1.7380000000000000D+03 ! Reference radius, [km]
            mu = 4.9028001224453001D+03 ! Gravitational parameter, [km^3/s^2]

            a0 = mu/ReferenceRadius**2*1000           ! Standard gravitational acceleration, [m/s^2]
            L = N*(N+3)/2               ! We need to read only the first L rows of coefficients.

            ! Scanning the first L rows of coefficients and extracting n, m, C_nm, and S_nm values...
            C = GVCoeffs(1:L,:)

            A = 0.0D0
            B = 0.0D0
            d = 0.0D0
            F = 0.0D0
            G = 0.0D0
            H = 0.0D0

            do i = 1,N+1 ! n -> i
    
                do j = 1,i ! m -> j
                    
                    if (j .le. min(M+2,N+1)) then
                        A(i,j) = sqrt( (2*i+1) * (2*i-1) / (i+j-1) / (i-j+1) )
                    end if
        
                    if (i > j .and. j <= min(M+2,N)) then
                        B(i-1,j) = A(i,j) / A(i-1,j)
                    end if
        
                end do
    
                if (i .le. M+1) then
                    if (i .eq. 1) then
                        temp_number = 1.0D0
                    else
                        temp_number = 0.0D0
                    end if
                    d(i) = sqrt( (2*i+1) / (2-temp_number) / i )
                end if
    
            end do

            do i = 1,N-1 ! n -> i
    
                do j = 1,i+2 ! m -> j
                    
                    if (j .eq. 1) then
                        temp_number = 1.0D0
                    else
                        temp_number = 0.0D0
                    end if
        
                    if (j .le. M+1) then
                        F(i,j) = sqrt( (2*i+3) * (i+j+2) * (i+j+1) / 2 / (2-temp_number) / (2*i+5) )
                        H(i,j) = sqrt( (2*i+3) * (i+j+1) * (i-j+3) / (2*i+5) )
                    end if
        
                    if (j .le. M) then
                        G(i,j) = sqrt( (2*i+3) * (i-j+3) * (i-j+2) / 2 / (2-temp_number) / (2*i+5) )
                    end if
        
                end do
    
            end do
    
            rnorm = sqrt(sum(r_vec*r_vec))
            
            QR = ReferenceRadius/rnorm
            dr = r_vec/rnorm
            Qx = dr(1)
            Qy = dr(2)
            Qz = dr(3)
            V = 0.0D0
            W = 0.0D0
            V(1,1) = QR

            if (grad_req) then
                dQR = -ReferenceRadius/rnorm**2*dr
                dQx = (/1,0,0/)/rnorm - r_vec(1)/rnorm**2*dr
                dQy = (/0,1,0/)/rnorm - r_vec(2)/rnorm**2*dr
                dQz = (/0,0,1/)/rnorm - r_vec(3)/rnorm**2*dr
                dV = 0.0D0
                dW = 0.0D0
                dV(1,1,1:3) = dQR
            end if

            do j = 1,M+1
    
                ! Calculating the sectorial harmonics...
                V(j+1,j+1) = d(j) * QR * ( Qx*V(j,j) - Qy*W(j,j) )
                W(j+1,j+1) = d(j) * QR * ( Qx*W(j,j) + Qy*V(j,j) )
                
                ! Calculating the semi-sectorial harmonics...
                V(j+1,j) = A(j,j) * QR * Qz * V(j,j)
                W(j+1,j) = A(j,j) * QR * Qz * W(j,j)
                
                if (grad_req) then
    
                    dV(j+1,j+1,:) = d(j) * dQR * ( Qx* V(j,j)   -  Qy* W(j,j) ) + &
                                    d(j) *  QR * (dQx* V(j,j)   - dQy* W(j,j) ) + &
                                    d(j) *  QR * ( Qx*dV(j,j,:) -  Qy*dW(j,j,:))
                    dW(j+1,j+1,:) = d(j) * dQR * ( Qx* W(j,j)   +  Qy* V(j,j) ) + &
                                    d(j) *  QR * (dQx* W(j,j)   + dQy* V(j,j) ) + &
                                    d(j) *  QR * ( Qx*dW(j,j,:) +  Qy*dV(j,j,:))

                    dV(j+1,j,:) = A(j,j) * dQR *  Qz *  V(j,j) + &
                                  A(j,j) *  QR * dQz *  V(j,j) + &
                                  A(j,j) *  QR *  Qz * dV(j,j,:)
                    dW(j+1,j,:) = A(j,j) * dQR *  Qz *  W(j,j) + &
                                  A(j,j) *  QR * dQz *  W(j,j) + &
                                  A(j,j) *  QR *  Qz * dW(j,j,:)
                    
                end if
    
            end do

            ! Calculating all the other harmonics...
            do ii = 1,N
    
                do ji = 1,min(M+2,ii)
        
                    V(ii+2,ji) = A(ii+1,ji) * QR * Qz * V(ii+1,ji) - B(ii,ji) * QR**2 * V(ii,ji)
                    W(ii+2,ji) = A(ii+1,ji) * QR * Qz * W(ii+1,ji) - B(ii,ji) * QR**2 * W(ii,ji)
        
                    if (grad_req) then
                        dV(ii+2,ji,:) = A(ii+1,ji) * dQR *  Qz *  V(ii+1,ji)   - B(ii,ji) * 2*QR*dQR *  V(ii,ji) + &
                                        A(ii+1,ji) *  QR * dQz *  V(ii+1,ji)   - B(ii,ji) * QR**2     * dV(ii,ji,:) + &
                                        A(ii+1,ji) *  QR *  Qz * dV(ii+1,ji,:)
                        dW(ii+2,ji,:) = A(ii+1,ji) * dQR *  Qz *  W(ii+1,ji)   - B(ii,ji) * 2*QR*dQR *  W(ii,ji) + &
                                        A(ii+1,ji) *  QR * dQz *  W(ii+1,ji)   - B(ii,ji) * QR**2     * dW(ii,ji,:) + &
                                        A(ii+1,ji) *  QR *  Qz * dW(ii+1,ji,:)
                        
                    end if
        
                end do
    
            end do

            ! Perturbing non-dimensional acceleration
            Qa = 0.0D0
            rel_harmonics = 0.0D0
            
            if (grad_req) then
                dQa = 0.0D0
                drh = 0.0D0
            end if

            do i = 1,N-1
    
                ! Calculating the zonal harmonic acceleration...
                rel_harmonics(1,i) = F(i,1) * C(i*(i+3)/2+1,3) * V(i+3,2)
                rel_harmonics(2,i) = F(i,1) * C(i*(i+3)/2+1,3) * W(i+3,2)
                rel_harmonics(3,i) = H(i,1) * C(i*(i+3)/2+1,3) * V(i+3,1)
    
                if (grad_req) then
                    drh(1,i,:) = F(i,1) * C(i*(i+3)/2+1,3) * dV(i+3,2,:)
                    drh(2,i,:) = F(i,1) * C(i*(i+3)/2+1,3) * dW(i+3,2,:)
                    drh(3,i,:) = H(i,1) * C(i*(i+3)/2+1,3) * dV(i+3,1,:)
                end if
    
                ! Calculating the acceleration due to the tesseral and sectorial harmonics...
                if (M .gt. 0) then
                    do j = 1,M
                        rel_harmonics(1,i) = rel_harmonics(1,i) + F(i,j+1) * ( C(i*(i+3)/2+j+1,3) * V(i+3,j+2) + C(i*(i+3)/2+j+1,4) * W(i+3,j+2) ) - G(i,j) * ( C(i*(i+3)/2+j+1,3) * V(i+3,j) + C(i*(i+3)/2+j+1,4) * W(i+3,j) )
                        rel_harmonics(2,i) = rel_harmonics(2,i) + F(i,j+1) * ( C(i*(i+3)/2+j+1,3) * W(i+3,j+2) - C(i*(i+3)/2+j+1,4) * V(i+3,j+2) ) + G(i,j) * ( C(i*(i+3)/2+j+1,3) * W(i+3,j) - C(i*(i+3)/2+j+1,4) * V(i+3,j) )
                        rel_harmonics(3,i) = rel_harmonics(3,i) + H(i,j+1) * ( C(i*(i+3)/2+j+1,3) * V(i+3,j+1) + C(i*(i+3)/2+j+1,4) * W(i+3,j+1) )
                        if (grad_req) then
                            drh(1,i,:) = drh(1,i,:) + F(i,j+1) * ( C(i*(i+3)/2+j+1,3) * dV(i+3,j+2,:) + C(i*(i+3)/2+j+1,4) * dW(i+3,j+2,:) ) - G(i,j) * ( C(i*(i+3)/2+j+1,3) * dV(i+3,j,:) + C(i*(i+3)/2+j+1,4) * dW(i+3,j,:) )
                            drh(2,i,:) = drh(2,i,:) + F(i,j+1) * ( C(i*(i+3)/2+j+1,3) * dW(i+3,j+2,:) - C(i*(i+3)/2+j+1,4) * dV(i+3,j+2,:) ) + G(i,j) * ( C(i*(i+3)/2+j+1,3) * dW(i+3,j,:) - C(i*(i+3)/2+j+1,4) * dV(i+3,j,:) )
                            drh(3,i,:) = drh(3,i,:) + H(i,j+1) * ( C(i*(i+3)/2+j+1,3) * dV(i+3,j+1,:) + C(i*(i+3)/2+j+1,4) * dW(i+3,j+1,:) )
                        end if
                    end do
                    
                    Qa = Qa + rel_harmonics(:,i)
                    
                    if (grad_req) then
                        dQa = dQa + drh(:,i,:)
                    end if
        
                end if
    
            end do

            rel_harmonics = -rel_harmonics
            a_grav = -a0 * (QR**2*r_vec/rnorm + Qa)
            a_cent = -a0 * QR**2*r_vec/rnorm
            
            da_grav = 0.0D0
            if (grad_req) then
                da_grav = -a0 * (2*QR*doto(r_vec,dQR)/rnorm + QR**2*eye(3)/rnorm - QR**2*doto(r_vec,dr)/rnorm**2 + dQa)
            end if
            
            ComplexMoonForce.Force = a_grav
            ComplexMoonForce.dFdr = da_grav
            
        end function kGravityMoonGradient50
        function kEarthAtm(state,grad_req) result(AtmForce)
        
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp) :: state(6), rvect(3), vvect(3), r, height
            real(dp) :: h0, ro_l, ro_h, cx, v, power, coeff, Fm, F(3), dF(3,6)
            real(dp) :: dheight(6), dvvect(3,6), dv(6), dpower(6), dFm(6)
            logical :: grad_req
            type(ForceStructure) :: AtmForce
        
            ! Earth units.

            rvect = state(1:3)
            vvect = state(4:6)

            r = norm(rvect)

            height = (r - 1.0D0)*6371.0D0

            if (height > 900) then
    
                h0   = 900.0D0
                ro_l = 1.0D0
                ro_h = 0.0D0
    
            elseif (height > 800 .and. height <= 900) then
    
                h0 = 800.0D0
                ro_l = 1.57D-14
                ro_h = 6.59D-15
    
            elseif (height > 700 .and. height <= 800) then
    
                h0 = 700.0D0
                ro_l = 4.80D-14
                ro_h = 1.57D-14
    
            elseif (height > 600 .and. height <= 700) then
    
                h0 = 600.0D0
                ro_l = 1.80D-13
                ro_h = 4.80D-14
    
            elseif (height > 500 .and. height <= 600) then
    
                h0 = 500.0D0
                ro_l = 7.85D-13
                ro_h = 1.80D-13
    
            elseif (height > 400 .and. height <= 500) then
    
                h0 = 400.0D0
                ro_l = 3.96D-12
                ro_h = 7.85D-13
    
            elseif (height > 300 .and. height <= 400) then
    
                h0 = 300.0D0
                ro_l = 2.52D-11
                ro_h = 3.96D-12
    
            elseif (height > 200 .and. height <= 300) then
    
                h0 = 200.0D0
                ro_l = 2.84D-10
                ro_h = 2.52D-11
    
            else ! height <= 200
    
                h0 = 100.0D0
                ro_l = 5.73D-07
                ro_h = 2.84D-10
    
            end if

            cx = 2.2D0

            v = norm(vvect)
            power = (ro_h/ro_l)**((height-h0)/100.0D0)
            coeff = 6371000.D0*(cx/2)*(area/mass)
            Fm = coeff*v**2*ro_l*power
            F = -Fm*vvect/v

            dF = 0.0D0
            if (grad_req) then
                if (height > 900.0D0) then
                    AtmForce.Force = F
                    AtmForce.dFdx = dF
                    return
                end if
                dheight = (/ 6371*rvect/r, 0.0D0, 0.0D0, 0.0D0 /)
                dvvect(1,:) = (/ 0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0 /)
                dvvect(2,:) = (/ 0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0 /)
                dvvect(3,:) = (/ 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0 /)
                
                dv = (/ 0.0D0, 0.0D0, 0.0D0, vvect/v /)
                dpower = power*log(ro_h/ro_l)*dheight/100
                dFm = coeff*2*v*dv*ro_l*power + coeff*v**2*ro_l*dpower
                dF = -doto(vvect,dFm)/v - Fm*(dvvect/v - doto(vvect,dv)/v**2)
            end if
            
            AtmForce.Force = F
            AtmForce.dFdx = dF
            
        end function kEarthAtm
        function kEarthJ2(jd,rvect,grad_req) result(J2Force)
        
            use Translations
            use BaseMeansToolbox
            use LinearAlgebraInterfaces
        
            implicit none
            
            real(dp) :: jd, rvect(3), rvect_ITRS(3), drvect_ITRS(3,3), r, x, y, z, J2c
            real(dp) :: RE, T1, T2, T3, T4, F_ITRS(3), F(3), dF_dF_ITRS(3,3), T5(3,3)
            real(dp) :: dF_ITRS(3,3), dF(3,3)
            real(dp) :: dr(3), dx(3), dy(3), dz(3), dT1(3), VT(3), dT2(3), dT3(3), dT4(3)
            logical :: grad_req
            type(ForceStructure) :: J2Force

            call kgcrs2itrs(rvect, jd, rvect_ITRS, drvect_ITRS)

            r = norm(rvect_ITRS)
            x = rvect_ITRS(1)
            y = rvect_ITRS(2)
            z = rvect_ITRS(3)
            J2c = 1.0826D-03
            RE = 6378.13D0/DistUnit

            T1 = - 3*J2c*RE/2/r**5
            T2 = (1-5*z**2/r**2)*x
            T3 = (1-5*z**2/r**2)*y
            T4 = (3-5*z**2/r**2)*z

            F_ITRS = T1 * (/ T2, T3, T4 /)

            call kitrs2gcrs(F_ITRS, jd, F, dF_dF_ITRS)

            if (grad_req) then
                dr = dotvm(rvect_ITRS,drvect_ITRS)/r
                dx = drvect_ITRS(1,:)
                dy = drvect_ITRS(2,:)
                dz = drvect_ITRS(3,:)
                dT1 = -5*T1/r * dr
                VT = -z*dz/r**2 + z**2/r**3*dr
                dT2 = (1-5*z**2/r**2)*dx + 10*x*VT
                dT3 = (1-5*z**2/r**2)*dy + 10*y*VT
                dT4 = (3-5*z**2/r**2)*dz + 10*z*VT
                T5(1,:) = dT2
                T5(2,:) = dT3
                T5(3,:) = dT4
                dF_ITRS = doto( (/ T2, T3, T4 /) , dT1 ) + T1 * T5
                dF = dotmm(dF_dF_ITRS,dF_ITRS)
            end if
            
            J2Force.Force = F
            J2Force.dFdr = dF

        end function kEarthJ2

end module
