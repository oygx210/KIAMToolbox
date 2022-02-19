module VisibilityModule

	use MathAndProgConstants
    
    contains
    
        subroutine kIsvisible(r_sat, lat_deg, long_deg, R, threshold_deg, vis_status, elev_deg, azim_deg)

			! Spacecraft visibility function for an observer located on the spherical body surface


			! -----------------------------------------------------
			! r_sat is a 3xN array of spacecraft position vectors
			! R is the body radius (in the same units as r_sat)
			! -----------------------------------------------------
			! -90 deg <= lat_deg <= 90 deg, -180 deg <= long_deg <= 180 deg
			! -----------------------------------------------------
			
			implicit none
			
			real(dp), intent(in) :: r_sat(:,:), lat_deg, long_deg, R, threshold_deg
			real(dp), intent(out) :: elev_deg(size(r_sat,2)), azim_deg(size(r_sat,2))
			integer, intent(out) :: vis_status(size(r_sat,2))
			integer :: N, i
			real(dp) :: d_sat(size(r_sat,2)), cent(size(r_sat,2)), lat, long, r_obs(3)
			real(dp) :: lat_sat(size(r_sat,2)), long_sat(size(r_sat,2))

			! Number of spacecraft position vectors
			N = size(r_sat, 2)

			! Distance from the body center to the spacecraft
			d_sat = 0.0D0
			do i = 1,N
				d_sat(i) = sqrt(sum(r_sat(:,i)*r_sat(:,i)))
				if (d_sat(i) <= R) then
					write(*,*) 'Spacecraft altitude is below zero!'
					return
				end if
			end do

			! Conversion to radians
			lat = lat_deg*pi/180.0D0
			long = long_deg*pi/180.0D0

			! Observer unit vector
			r_obs = [cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)]

			! Central angle
			cent = 0.0D0

			! Elevation angle (in degrees)
			elev_deg = 0.0D0

			! Subsatellite point latitude
			lat_sat = 0.0D0

			! Subsatellite point longitude
			long_sat = 0.0D0

			! Geodesic azimuth angle (in degrees)
			azim_deg = 0.0D0

			! Status is 1 if a spacecraft is visible, and 0 otherwise
			vis_status = 0.0D0

			do i = 1, N
    
				lat_sat(i) = asin( r_sat(3,i) / d_sat(i) )
				long_sat(i) = atan2( r_sat(2,i), r_sat(1,i) )
    
				cent(i) = acos( sum(r_sat(:,i)*r_obs) / d_sat(i) )
				elev_deg(i) = atan2( cos(cent(i)) - R/d_sat(i), sin(cent(i)) ) * 180.0D0/pi
				azim_deg(i) = atan2( sin(long_sat(i)-long)*cos(lat_sat(i)), cos(lat)*sin(lat_sat(i)) - sin(lat)*cos(lat_sat(i))*cos(long_sat(i)-long) ) * 180.0D0/pi
    
				if (azim_deg(i) < 0.0D0) then
					azim_deg(i) = 360.0D0 + azim_deg(i)
				end if
    
				if (elev_deg(i) >= threshold_deg) then
					vis_status(i) = 1
				end if
    
			end do

		end subroutine kIsvisible
    
    
end module VisibilityModule