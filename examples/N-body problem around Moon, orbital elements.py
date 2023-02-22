import numpy as np
from kiam_astro import kiam
from kiam_astro.trajectory import Trajectory

# Restricted n-body problem, moon, orbital elements
t0 = 0.0
s0 = np.array([3.0, 0.01, kiam.deg2rad(80), 0.0, 0.0, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory(s0, t0, jd0, 'oe', 'mer', 'moon')
tr.set_model('rv', 'nbp', 'moon', ['earth', 'sun', 'cmplxmoon'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 10  # order of the Moon's gravity field
tr.propagate(10000, 10000)  # (time of flight, number of points)
tr.change_system('mer')  # change back to MER system
tr.change_vars('oe')
tr.show('e')  # show the eccentricity
