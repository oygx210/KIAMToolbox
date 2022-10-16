import kiam

# Juliandate
jd = kiam.juliandate(2022, 2, 15, 0, 0, 0)

# rv 2 oe
oe = kiam.rv2oe([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 1.0)

# Moon units
ku = kiam.units('Moon')

# r2bp right hand of equations (mu = 1.0)
dxdt = kiam.r2bp(0.0, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0])

# n-body problem, moon
sources = {}  # order of perturbations is important!
sources['sun'] = False
sources['mercury'] = False
sources['venus'] = False
sources['earth'] = True
sources['moon'] = False
sources['mars'] = False
sources['jupiter'] = False
sources['saturn'] = False
sources['uranus'] = False
sources['neptune'] = False
sources['srp'] = False
sources['cmplxmoon'] = True
sources['atm'] = False
sources['j2'] = False
data = {}
data['jd_zero'] = kiam.juliandate(2022, 2, 15, 0, 0, 0)
data['area'] = 2.0
data['mass'] = 100.0
data['order'] = 50
units_data = {}
ku = kiam.units('Moon')
_, star, planet, moon, _ = kiam.astro_const()
units_data['DistUnit'] = ku['DistUnit']
units_data['VelUnit'] = ku['VelUnit']
units_data['TimeUnit'] = ku['TimeUnit']
units_data['AccUnit'] = ku['AccUnit']
units_data['RSun'] = star['Sun']['MeanRadius'] / ku['DistUnit']  # 695700 km
units_data['REarth'] = planet['Earth']['EquatorRadius'] / ku['DistUnit']  # 6378.1366 km
units_data['RMoon'] = moon['Moon']['MeanRadius'] / ku['DistUnit']  # 1737.4 km
dxdt = kiam.nbp_rv_moon(0.0, [1.5, 0.0, 0.0, 0.0, 0.9, 0.0], False, sources, data, units_data)
print(dxdt)
