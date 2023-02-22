import numpy as np
from kiam_astro import kiam
from kiam_astro.trajectory import Trajectory

t0 = 0.0
s0 = [2.0, 0.0, 0.0, 0.0, 1/np.sqrt(2.0), 0.0]
s0.extend(list(kiam.eye2vec(6)))
s0 = np.array(s0)
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory(s0, t0, jd0, 'rv_stm', 'gcrs', 'earth')
tr.set_model('rv_stm', 'nbp', 'earth', ['moon', 'sun'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 0  # order of the Moon's gravity field
tr.propagate(6*np.pi, 20000)  # (time of flight, number of points)
tr.show('3d')  # show the trajectory in 3d
