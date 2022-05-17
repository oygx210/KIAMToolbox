import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import kiam

x = np.linspace(0, 1, 100)
y1 = np.exp(x)
y2 = np.exp(x**2)

ax = kiam.plot(x, y1)
_, ax = kiam.plot(x, y2, ax, show=True)
