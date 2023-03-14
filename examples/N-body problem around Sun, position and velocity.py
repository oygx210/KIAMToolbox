import numpy
from numpy import pi
from kiam_astro import kiam
from kiam_astro import trajectory

units = kiam.units('sun')

t0 = 0.0
s0 = numpy.array([1, 0, 0, 0, 1, 0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = trajectory.Trajectory(s0, t0, jd0, 'rv', 'hcrs', 'sun')
tr.set_model('rv', 'nbp', 'sun', ['earth', 'jupiter', 'saturn'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.propagate(2*pi, 100)  # (time of flight, number of points)
tr.show('xy')
