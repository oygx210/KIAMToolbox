import kiam
import numpy as np
import Trajectory
import os

# os.environ['BROWSER'] = 'firefox'

# Juliandate (year, month, day, hour, minute, second)
jd = kiam.juliandate(2022, 4, 30, 0, 0, 0)

# (x, y, z, vx, vy, vz) -> (a, e, i, Omega, omega, theta)
oe = kiam.rv2oe(np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]), 1.0)

# Moon units
ku = kiam.units('Moon')

# Two-body problem right hand of equations (mu = 1.0)
dxdt = kiam.r2bp(0.0, np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0]))

# State of the Moon wrt to the Earth
x = kiam.planet_state(kiam.juliandate(2022, 9, 1, 0, 0, 0), 'Earth', 'Moon')
print(x)

# Restricted two-body problem, earth
t0 = 0.0
s0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'gcrs', 'earth')
tr.set_model('rv', 'r2bp', 'earth', [])
tr.propagate(2*np.pi, 1000)  # (time of flight, number of points)
fig = tr.show('xy', draw=False)  # show the trajectory in 3d
fig = kiam.axis_equal(fig)  # set the axis to be equal
fig.show()

# Hill's problem
t0 = 0.0
s0 = np.array([-0.5, 0.0, 0.0, 0.0, 1.5, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'rot_sb', 'earth_moon')
tr.set_model('rv', 'hill', 'earth_moon', [])
tr.model['data']['t0'] = 0.0  # time for which rot and ine systems coincide
tr.propagate(2*np.pi, 1000)  # time of flight, number of points
fig = tr.show('xy', draw=False)  # show the trajectory in 3d
fig = kiam.axis_equal(fig)  # set the axis to be equal
fig.show()

# Circular restricted three-body problem
t0 = 0.0
s0 = np.array([0.5, 0.0, 0.0, 0.0, 0.1, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'rot_fb', 'earth_moon')
tr.set_model('rv', 'cr3bp_fb', 'earth_moon', [])
tr.model['data']['t0'] = 0.0  # time for which rot and ine systems coincide
tr.propagate(2*np.pi, 1000)  # time of flight, number of points
fig = tr.show('xy', draw=False)  # show the trajectory in 3d
fig = kiam.axis_equal(fig)  # set the axis to be equal
fig.show()

# Non-pertubed two-body problem, earth
t0 = 0.0
s0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'gcrs', 'earth')
tr.set_model('rv', 'nbp', 'earth', [])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 1  # order of the Moon's gravity field
tr.propagate(2*np.pi, 1000)  # (time of flight, number of points)
tr.show('xy')  # show the trajectory in 3d

# Restricted n-body problem, moon, orbital elements
t0 = 0.0
s0 = np.array([3.0, 0.01, kiam.deg2rad(80), 0.0, 0.0, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'oe', 'mer', 'moon')
tr.set_model('rv', 'nbp', 'moon', ['earth', 'sun', 'cmplxmoon'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 10  # order of the Moon's gravity field
tr.propagate(10000, 10000)  # (time of flight, number of points)
tr.change_system('mer')  # change back to MER system
tr.change_vars('oe')
tr.show('e')  # show the eccentricity

# Pertubed n-body problem, moon
t0 = 0.0
s0 = np.array([2.0, 0.0, 0.0, 0.0, 1/np.sqrt(2.0), 0.3])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'scrs', 'moon')
tr.set_model('rv', 'nbp', 'moon', ['earth', 'sun'])
tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
tr.model['data']['mass'] = 100.0  # spacecraft mass, kg
tr.model['data']['area'] = 2.0  # spacecraft area, m^2
tr.model['data']['order'] = 1  # order of the Moon's gravity field
tr.propagate(24*np.pi, 20000)  # (time of flight, number of points)
fig = tr.show('3d')  # show the trajectory in 3d
fig = kiam.axis_equal(fig)
fig.show()

# Pertubed n-body problem, earth, variational equations
t0 = 0.0
s0 = [2.0, 0.0, 0.0, 0.0, 1/np.sqrt(2.0), 0.0]
s0.extend(list(kiam.eye2vec(6)))
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
