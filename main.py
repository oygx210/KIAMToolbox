import FKIAMToolbox as fkt
import numpy as np


print('Running test_nbp_equations...')

fkt.equationsmodule.stm_required = 0

fkt.equationsmodule.atm = 1
fkt.equationsmodule.j2 = 1
fkt.equationsmodule.srp = 1
fkt.equationsmodule.sun = 0
fkt.equationsmodule.mercury = 0
fkt.equationsmodule.venus = 0
fkt.equationsmodule.earth = 0
fkt.equationsmodule.mars = 0
fkt.equationsmodule.jupiter = 0
fkt.equationsmodule.saturn = 0
fkt.equationsmodule.uranus = 0
fkt.equationsmodule.neptune = 0
fkt.equationsmodule.moon = 0
fkt.equationsmodule.cmplxmoon = 0

fkt.equationsmodule.jd_zero = 2459599.5
fkt.equationsmodule.order = 50
fkt.equationsmodule.area = 1
fkt.equationsmodule.mass = 100

fkt.equationsmodule.distunit = fkt.equationsmodule.earthdistunit
fkt.equationsmodule.velunit = fkt.equationsmodule.earthvelunit
fkt.equationsmodule.timeunit = fkt.equationsmodule.earthtimeunit
fkt.equationsmodule.accunit = fkt.equationsmodule.earthaccunit
fkt.equationsmodule.rsun = fkt.equationsmodule.earthrsun
fkt.equationsmodule.rearth = fkt.equationsmodule.earthrearth
fkt.equationsmodule.rmoon = fkt.equationsmodule.earthrmoon

print('Running knbp_rv_earth...')
result_rv_earth = fkt.equationsmodule.knbp_rv_earth(0, np.array([1.1, 0, 0, 0, 0.9, 0]))
print('Done knbp_rv_earth.')
result_rv_true = np.array([0.,  0.9,  0., -8.275566130556367e-01,
                           -1.0601077609733587e-07, -4.751264387875831e-06])
