import numpy as np
import kiam

r = np.linspace(0, 1, 100)
phi = np.linspace(0, 2*np.pi, 100)

kiam.polarplot(r, phi, 'r--', rmax=1)