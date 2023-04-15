import numpy as np
from kiam_astro import kiam
from kiam_astro.trajectory import Trajectory

t0 = 0.0
s0 = np.array([2.0, 0.0, 0.0, 0.0, 1/np.sqrt(2.0), 0.3])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory(s0, t0, jd0, 'rv', 'scrs', 'moon')
tr.set_model('rv', 'nbp', 'moon', ['earth', 'sun', 'cmplxmoon'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 10  # order of the Moon's gravity field
tr.propagate(24*np.pi, 20000)  # (time of flight, number of points)
fig = tr.show('3d', draw=False)  # show the trajectory in 3d
fig = kiam.set_axis_equal(fig)
fig.show()