import numpy as np
import kiam
import Trajectory
from Model import Model

t0 = 0.0
s0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
jd0 = kiam.juliandate(2022, 4, 30, 0, 0, 0)
tr = Trajectory.Trajectory(s0, t0, jd0, 'rv', 'scrs', 'moon')
tr.set_model('rv', 'nbp', 'moon', [])
db = {'asdf': 123, 'tr': tr}

print(db)

kiam.save(db, 'trajectory')

db2 = kiam.load('trajectory')

print(db2)