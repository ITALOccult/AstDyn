#!/usr/bin/env python3
"""
Check Earth position on 2025-01-12 (MJD 60693.484150)
"""

from astropy.time import Time
from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
import numpy as np

# Observation time
mjd_utc = 60693.484150
t = Time(mjd_utc, format='mjd', scale='utc')

print(f"Date: {t.iso}")
print(f"MJD:  {mjd_utc}")

# Use builtin ephemeris
with solar_system_ephemeris.set('builtin'):
    earth_pos = get_body_barycentric('earth', t)
    
print(f"\nEarth position (ICRS barycentric):")
print(f"  X = {earth_pos.x.value:.6f} AU")
print(f"  Y = {earth_pos.y.value:.6f} AU")
print(f"  Z = {earth_pos.z.value:.6f} AU")

# Our code gives (ecliptic): [-0.466, 0.866, -0.00005]
# Convert to check
print(f"\nOur code (ecliptic): [-0.466, 0.866, -0.00005]")

# Obliquity
epsilon = 23.439291 * np.pi / 180.0
c = np.cos(epsilon)
s = np.sin(epsilon)

# Rotation from equatorial to ecliptic: R_x(-epsilon)
R_eq_to_ecl = np.array([
    [1.0,  0.0, 0.0],
    [0.0,    c,  -s],
    [0.0,    s,   c]
])

earth_eq = np.array([earth_pos.x.value, earth_pos.y.value, earth_pos.z.value])
earth_ecl = R_eq_to_ecl @ earth_eq

print(f"Astropy (converted to ecliptic):")
print(f"  X = {earth_ecl[0]:.6f} AU")
print(f"  Y = {earth_ecl[1]:.6f} AU")
print(f"  Z = {earth_ecl[2]:.6f} AU")
