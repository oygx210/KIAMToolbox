import numpy as np
from kiam_astro import kiam
from kiam_astro.trajectory import Trajectory

# Circular restricted three-body problem
t0 = 0.0
s0 = np.array([0.5, 0.0, 0.0, 0.0, 0.1, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory(s0, t0, jd0, 'rv', 'rot_fb', 'earth_moon')
tr.set_model('rv', 'cr3bp_fb', 'earth_moon', [])
tr.model['data']['t0'] = 0.0  # time for which rot and ine systems coincide
tr.propagate(2*np.pi, 1000)  # time of flight, number of points
fig = tr.show('xy', draw=False)  # show the trajectory in 3d
fig = kiam.set_axis_equal(fig)  # set the axis to be equal
fig.show()
