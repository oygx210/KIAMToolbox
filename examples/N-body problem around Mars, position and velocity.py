import numpy
from numpy import sqrt, pi
from kiam_astro import kiam
from kiam_astro import trajectory

units = kiam.units('mars')

t0 = 0.0
s0 = numpy.array([1.0, 0, 0, 0, sqrt(1/1.0), 0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = trajectory.Trajectory(s0, t0, jd0, 'rv', 'mcrs', 'mars')
tr.set_model('rv', 'nbp', 'mars', ['jupiter', 'saturn'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.propagate(10*2*pi, 10000)  # (time of flight, number of points)
tr.show('xy')
