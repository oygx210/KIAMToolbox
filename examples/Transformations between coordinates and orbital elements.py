import numpy as np
from kiam_astro import kiam

# Cartesian coordinates (x, y, z, vx, vy, vz)
rv = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

# (x, y, z, vx, vy, vz) -> (a, e, i, Omega, omega, theta)
oe = kiam.rv2oe(rv, mu=1.0)

# (x, y, z, vx, vy, vz) -> (h, ex, ey, ix, iy, L)
ee = kiam.rv2ee(rv, mu=1.0)

# Orbital elements (a, e, i, Omega, omega, theta)
oe = np.array([1.1, 0.01, np.pi, 0.0, 0.0, 0.0])

# (a, e, i, Omega, omega, theta) -> (x, y, z, vx, vy, vz)
rv = kiam.oe2rv(oe, mu=1.0)




