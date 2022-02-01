import Model as md

model = md.Model('rv', 'nbp', 'moon', ['Moon', 'Sun', 'SRP'], 2459599.5)
model.data['mass'] = 100.0
model.data['order'] = 10
model.data['area'] = 1
dfdt = model.eqs(0, [10, 0, 0, 0, 1, 0])
print(dfdt)
