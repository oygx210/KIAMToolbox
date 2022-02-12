import Trajectory
import numpy as np
import kiam

t0 = 0.0
x0 = np.array([1.5, 0.0, 0.0, 0.0, 1.0/np.sqrt(1.5), 0.1])
jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'scrs', 'moon')
tr.set_model('rv', 'nbp', 'moon', ['Sun'])
tr.model.data['jd_zero'] = jd0
tr.propagate(4*np.pi, 10000)
vis_status, elev_deg, azim_deg = kiam.is_visible(tr.states[0:3, :], 0.0, 0.0, 1.0, 0.0)
print(1)
