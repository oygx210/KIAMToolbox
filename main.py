import Trajectory
import numpy as np
import kiam
import networkx as nx
import matplotlib.pyplot as plt
import math

t0 = 0.0
x0 = np.array([1.5, 0, 0, 0, 1/math.sqrt(1.5), 0.1])
jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gcrs', 'earth')
tr.set_model('rv', 'nbp', 'earth', ['Sun', 'Moon'])
tr.model.data['jd_zero'] = jd0
tr.propagate(4*math.pi, 10000)
tr.change_vars('ee')
tr.show('ex')
