from kiam_astro import kiam

# Moon's state with respect to the Earth (ITRS)
jd = kiam.juliandate(2022, 9, 1, 0, 0, 0)  # TDB
x = kiam.planet_state(jd, center='Earth', target='Moon')
r_Moon = x[0:3]  # Moon's position (km)
v_Moon = x[3:6]  # Moon's velocity (km/s)

print(r_Moon)
print(v_Moon)
