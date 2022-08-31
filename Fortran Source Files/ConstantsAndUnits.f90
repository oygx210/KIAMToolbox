module ConstantsAndUnits

use MathAndProgConstants

! -------------------------------------------------------------
! AstroConst  Contains the IAU 2009/2012 astronomical constants
!
! N.B. Latest IAU and NIST current best estimates are taken into account.
! The orbital elements of all the planets relate to the mean ecliptic and
! equinox of the J2000.0 epoch. The lunar semimajor axis value is derived
! from the Lunar Laser Ranging Experiment (LLRE).
!
! uniConst        Structure with the universal constants
! star            Structure with the parameters of the Sun
! planet          Structure with the parameters of the Solar System planets
! moon            Structure with the parameters of the planetary satellites
! smallBody       Structure with the parameters of the Solar System small bodies
!
! Authors: Sergey Trofimov (Matlab), Maksim Shirobokov (Fortran)
! Dates: Feb 03, 2020, Jan 29, 2022

! Universal constants
real(dp), parameter :: uniConst_SoL = 299792.458D0	! speed of light, km/s
real(dp), parameter :: uniConst_AU  = 149597870.7D0	! astronomical unit, km
real(dp), parameter :: uniConst_G   = 6.67408D-20	! constant of gravitation, km^3/kg/s^2
real(dp), parameter :: uniConst_g0  = 9.80665D0		! standard acceleration due to gravity, m/s^2
real(dp), parameter :: uniConst_RAD = 180.0D0/pi    ! degrees in 1 radian

! Sun
real(dp), parameter :: Sun_GM         = 1.3271244004D+11	! km^3/s^2  (TDB)
real(dp), parameter :: Sun_MeanRadius = 695700.0D0          ! km

! Earth
character(3), parameter :: Earth_OrbitsAround  = 'Sun'
real(dp), parameter :: Earth_GM            = 398600.4356D0                            ! km^3/s^2  (TDB)
real(dp), parameter :: Earth_MeanRadius    = 6371.0084D0                              ! km
real(dp), parameter :: Earth_EquatorRadius = 6378.1366D0                              ! km
real(dp), parameter :: Earth_SemimajorAxis = 1.00000261D0*uniConst_AU                 ! km        (J2000.0)
real(dp), parameter :: Earth_Obliquity     = 23.0D0+(26.0D0*60.0D0+21.406D0)/3600.0D0 ! deg       (J2000.0)
real(dp), parameter :: Earth_dObliquitydt  = -46.836769D0                             ! arcsec/cy (TDB)

! Moon
character(5), parameter :: Moon_OrbitsAround  = 'Earth'
real(dp), parameter :: Moon_GM            = 0.0123000371D0*Earth_GM      ! km^3/s^2
real(dp), parameter :: Moon_MeanRadius    = 1737.4D0                     ! km
real(dp), parameter :: Moon_EquatorRadius = 1738.4D0                     ! km
real(dp), parameter :: Moon_SemimajorAxis = 384402.0D0                   ! km        (LLRE)

! Mercury
character(3), parameter :: Mercury_OrbitsAround  = 'Sun'
real(dp), parameter :: Mercury_GM            = Sun_GM/6.0236D+06              ! km^3/s^2
real(dp), parameter :: Mercury_MeanRadius    = 2439.4D0                       ! km
real(dp), parameter :: Mercury_EquatorRadius = 2440.53D0                      ! km
real(dp), parameter :: Mercury_SemimajorAxis = 0.38709927D0*uniConst_AU       ! km        (J2000.0)

! Venus
character(3), parameter :: Venus_OrbitsAround  = 'Sun'
real(dp), parameter :: Venus_GM            = Sun_GM/4.0852372e+05             ! km^3/s^2
real(dp), parameter :: Venus_MeanRadius    = 6051.8D0                         ! km
real(dp), parameter :: Venus_EquatorRadius = 6051.8D0                         ! km
real(dp), parameter :: Venus_SemimajorAxis = 0.72333566D0*uniConst_AU         ! km        (J2000.0)

! Mars
character(3), parameter :: Mars_OrbitsAround  = 'Sun'
real(dp), parameter :: Mars_GM            = Sun_GM/3.09870359D+06             ! km^3/s^2
real(dp), parameter :: Mars_MeanRadius    = 3389.50D0                         ! km
real(dp), parameter :: Mars_EquatorRadius = 3396.19D0                         ! km
real(dp), parameter :: Mars_SemimajorAxis = 1.52371034D0*uniConst_AU          ! km        (J2000.0)

! Jupiter
character(3), parameter :: Jupiter_OrbitsAround  = 'Sun'
real(dp), parameter :: Jupiter_GM            = Sun_GM/1.047348644D+03         ! km^3/s^2
real(dp), parameter :: Jupiter_MeanRadius    = 69911.0D0                      ! km
real(dp), parameter :: Jupiter_EquatorRadius = 71492.0D0                      ! km
real(dp), parameter :: Jupiter_SemimajorAxis = 5.20288700D0                   ! km        (J2000.0)

! Saturn
character(3), parameter :: Saturn_OrbitsAround  = 'Sun'
real(dp), parameter :: Saturn_GM            = Sun_GM/3.4979018D+03            ! km^3/s^2
real(dp), parameter :: Saturn_MeanRadius    = 58232.0D0                       ! km
real(dp), parameter :: Saturn_EquatorRadius = 60268.0D0                       ! km
real(dp), parameter :: Saturn_SemimajorAxis = 9.53667594D0                    ! km        (J2000.0)

! Uranus
character(3), parameter :: Uranus_OrbitsAround  = 'Sun'
real(dp), parameter :: Uranus_GM            = Sun_GM/2.290298D+04             ! km^3/s^2
real(dp), parameter :: Uranus_MeanRadius    = 25362.0D0                       ! km
real(dp), parameter :: Uranus_EquatorRadius = 25559.0D0                       ! km
real(dp), parameter :: Uranus_SemimajorAxis = 19.18916464D0                   ! km        (J2000.0)

! Neptune
character(3), parameter :: Neptune_OrbitsAround  = 'Sun'
real(dp), parameter :: Neptune_GM            = Sun_GM/1.941226D+04            ! km^3/s^2
real(dp), parameter :: Neptune_MeanRadius    = 24622.0D0                      ! km
real(dp), parameter :: Neptune_EquatorRadius = 24764.0D0                      ! km
real(dp), parameter :: Neptune_SemimajorAxis = 30.06992276D0                  ! km        (J2000.0)

! Pluto
character(3), parameter :: Pluto_OrbitsAround  = 'Sun'
real(dp), parameter :: Pluto_GM            = Sun_GM/1.36566D+08               ! km^3/s^2
real(dp), parameter :: Pluto_MeanRadius    = 1188.3D0                         ! km
real(dp), parameter :: Pluto_EquatorRadius = 1188.3D0                         ! km
real(dp), parameter :: Pluto_SemimajorAxis = 39.48211675D0                    ! km        (J2000.0)

	contains
	
		subroutine kunits_onebody(PrimaryBody, GM, DistUnit, VelUnit, TimeUnit, AccUnit)
		
			implicit none
			
			character(*), intent(in) :: PrimaryBody
			real(dp), intent(out) :: GM, DistUnit, VelUnit, TimeUnit, AccUnit
			
			if (PrimaryBody .eq. 'moon') then
				GM = Moon_GM 
				DistUnit = Moon_MeanRadius
			else if (PrimaryBody .eq. 'earth') then
				GM = Earth_GM 
				DistUnit = Earth_MeanRadius
			else if (PrimaryBody .eq. 'sun') then
				GM = Sun_GM
				DistUnit = Sun_MeanRadius
			else if (PrimaryBody .eq. 'mercury') then
				GM = Mercury_GM
				DistUnit = Mercury_MeanRadius
			else if (PrimaryBody .eq. 'venus') then
				GM = Venus_GM
				DistUnit = Venus_MeanRadius
			else if (PrimaryBody .eq. 'mars') then
				GM = Mars_GM
				DistUnit = Mars_MeanRadius
			else if (PrimaryBody .eq. 'jupiter') then
				GM = Jupiter_GM
				DistUnit = Jupiter_MeanRadius
			else if (PrimaryBody .eq. 'saturn') then
				GM = Saturn_GM
				DistUnit = Saturn_MeanRadius
			else if (PrimaryBody .eq. 'uranus') then
				GM = Uranus_GM
				DistUnit = Uranus_MeanRadius
			else if (PrimaryBody .eq. 'neptune') then
				GM = Neptune_GM
				DistUnit = Neptune_MeanRadius
			else if (PrimaryBody .eq. 'pluto') then
				GM = Pluto_GM
				DistUnit = Pluto_MeanRadius
			end if
				
			VelUnit = sqrt(GM/DistUnit)               ! km/s
			TimeUnit = DistUnit/VelUnit/86400.0D0     ! days
			AccUnit  = VelUnit**2/DistUnit*1000.0D0   ! m/s^2
		
		end subroutine kunits_onebody

		subroutine kunits_twobody(PrimaryBody, SecondaryBody, GM, mu, DistUnit, VelUnit, TimeUnit, AccUnit)
		
			implicit none
			
			character(*), intent(in) :: PrimaryBody, SecondaryBody
			real(dp), intent(out) :: GM, mu, DistUnit, VelUnit, TimeUnit, AccUnit
			real(dp) :: GMPB, GMSB, SMASB
			
			if ((PrimaryBody .eq. 'earth') .and. (SecondaryBody .eq. 'moon')) then
				GMPB = Earth_GM
				GMSB = Moon_GM
				SMASB = Moon_SemimajorAxis
			else if ((PrimaryBody .eq. 'sun') .and. (SecondaryBody .eq. 'earth')) then
				GMPB = Sun_GM
				GMSB = Earth_GM + Moon_GM
				SMASB = Earth_SemimajorAxis
			else if ((PrimaryBody .eq. 'sun') .and. (SecondaryBody .eq. 'venus')) then
				GMPB = Sun_GM
				GMSB = Venus_GM
				SMASB = Venus_SemimajorAxis
			else if ((PrimaryBody .eq. 'sun') .and. (SecondaryBody .eq. 'mars')) then
				GMPB = Sun_GM
				GMSB = Mars_GM
				SMASB = Mars_SemimajorAxis
			else if ((PrimaryBody .eq. 'sun') .and. (SecondaryBody .eq. 'jupiter')) then
				GMPB = Sun_GM
				GMSB = Jupiter_GM
				SMASB = Jupiter_SemimajorAxis
			else if ((PrimaryBody .eq. 'sun') .and. (SecondaryBody .eq. 'saturn')) then
				GMPB = Sun_GM
				GMSB = Saturn_GM
				SMASB = Saturn_SemimajorAxis
			end if
			
			GM = GMPB + GMSB
			mu = GMSB / GM
			DistUnit = SMASB
            VelUnit  = sqrt(GM/DistUnit)             ! km/s
            TimeUnit = DistUnit/VelUnit/86400.0D0    ! days
            AccUnit  = VelUnit**2/DistUnit*1000.0D0  ! m/s^2
		
		end subroutine kunits_twobody

end module ConstantsAndUnits