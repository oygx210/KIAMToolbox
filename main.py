import Trajectory
import numpy as np
import kiam
import networkx as nx
import matplotlib.pyplot as plt
import math

t0 = 0.0
x0 = np.array([1.5, 0, 0, 0, 1/math.sqrt(1.5), 0])
jd0 = kiam.juliandate(2022, 2, 2, 0, 0, 0)
tr = Trajectory.Trajectory(x0, t0, jd0, 'rv', 'gcrs', 'earth')
tr.set_model('rv', 'nbp', 'earth', ['Sun', 'Moon'])
tr.model.data['jd_zero'] = jd0
tr.propagate(2*math.pi)


tr.change_units('moon')
tr.vars_transform('rv', 'ee')
tr.vars_transform('ee', 'rv')
tr.allocate_units_graph()
G = tr.units_graph
ax = plt.plot()
pos = nx.kamada_kawai_layout(G)
# pos = nx.spectral_layout(G)
nx.draw(G, pos=pos)
nx.draw_networkx_labels(G, pos=pos)
plt.show()

G = nx.DiGraph()
G.add_edge(('a', 'qwer'), 2, weight=0.1)
G.add_edge(2, 3, weight=0.2)
G.add_edge(2, 4, weight=0.3)
G.add_edge(3, 5, weight=0.3)
G.add_edge(1, 5, weight=10)
ax = plt.plot()
pos = nx.shell_layout(G)
nx.draw(G, pos=pos)
nx.draw_networkx_labels(G, pos=pos)
nx.draw_networkx_edge_labels(G, pos=pos)
plt.show()
print(nx.shortest_path(G, 1, 5))



#print(tr)
