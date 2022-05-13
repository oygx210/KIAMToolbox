import numpy as np
import kiam
import Trajectory
from numpy import sqrt

t0 = 0.0
s0 = np.array([2.0, 0.0, 0.0, 0.0, 1/sqrt(2.0), 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'scrs', 'moon')
tr.set_model('rv', 'nbp', 'moon', [])
tr.model['data']['jd_zero'] = jd0
tr.model['data']['mass'] = 100.0
tr.model['data']['area'] = 2.0
tr.model['data']['order'] = 1
tr.propagate(2*np.pi, 200000)
tr.change_system('mer')

tr.show('3d')
