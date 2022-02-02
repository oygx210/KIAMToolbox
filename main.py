import Trajectory
import numpy as np
import kiam

t0 = 0.0
x0 = np.array([1, 0, 0, 0, 1, 0])
jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gcrs', 'earth')
tr.set_model('rv', 'nbp', 'earth', ['Sun', 'Moon'])
tr.vars_transform('rv', 'ee')

print(tr)
