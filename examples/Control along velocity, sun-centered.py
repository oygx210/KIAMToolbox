import numpy
from numpy.linalg import norm
from kiam_astro import kiam
from kiam_astro.engine import SPT50
from kiam_astro.trajectory import Trajectory
from numpy import pi

units = kiam.units('sun')

engine = SPT50()
engine.force /= units['AccUnit']  # kg * (nondim. acceler.)
engine.specific_impulse /= units['TimeUnit'] * 24 * 3600  # (nondim. time)

def control(t, x):
    force_vector = x[3:6]/norm(x[3:6]) * engine.force
    specific_impulse = engine.specific_impulse
    return force_vector, specific_impulse

t0 = 0.0
jd0 = kiam.juliandate(2023, 1, 1, 0, 0, 0)
s0 = numpy.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 100.0])
tr = Trajectory(s0, t0, jd0, 'rvm', 'hcrs', 'sun')
tr.set_model('rvm', 'nbp', 'sun', ['jupiter'])
tr.model['data']['jd_zero'] = jd0
tr.model['data']['area'] = 2.0
tr.model['data']['mass'] = s0[6]
tr.model['control'] = control
tr.propagate(2*pi/12, 100)
tr.show('xy')
# See also tr.control_history, tr.specific_impulse_history
