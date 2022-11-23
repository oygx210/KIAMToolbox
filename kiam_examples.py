import kiam
import numpy as np
import Trajectory

# Juliandate (year, month, day, hour, minute, second)
jd = kiam.juliandate(2022, 2, 15, 0, 0, 0)

# (x, y, z, vx, vy, vz) -> (a, e, i, Omega, omega, theta)
oe = kiam.rv2oe([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 1.0)

# Moon units
ku = kiam.units('Moon')

# two-body problem right hand of equations (mu = 1.0)
dxdt = kiam.r2bp(0.0, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0])

# n-body problem, moon
t0 = 0.0
s0 = np.array([2.0, 0.0, 0.0, 0.0, 1/np.sqrt(2.0), 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'scrs', 'moon')
tr.set_model('rv', 'nbp', 'moon', ['earth', 'sun'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 1  # order of the Moon's gravity field
tr.propagate(6*np.pi, 20000)  # (time of flight, number of points)
# tr.show('3d')  # show the trajectory in 3d

# n-body problem, earth, variational equations
t0 = 0.0
s0 = [2.0, 0.0, 0.0, 0.0, 1/np.sqrt(2.0), 0.0]
s0.extend(list(kiam.eyevec(6)))
s0 = np.array(s0)
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv_stm', 'gcrs', 'earth')
tr.set_model('rv_stm', 'nbp', 'earth', ['moon', 'sun'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 0  # order of the Moon's gravity field
tr.propagate(6*np.pi, 20000)  # (time of flight, number of points)
tr.show('3d')  # show the trajectory in 3d
