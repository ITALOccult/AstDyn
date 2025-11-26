#!/usr/bin/env python3
"""
Debug coordinate transformations
Following EXACTLY the OrbFit algorithm
"""

import numpy as np

# Constants
epsilon = 23.439291 * np.pi / 180.0  # Obliquity

# Rotation matrix from ECLIPTIC to EQUATORIAL
# This is R_x(+epsilon)
def roteceq():
    c = np.cos(epsilon)
    s = np.sin(epsilon)
    return np.array([
        [1.0,  0.0, 0.0],
        [0.0,    c,   s],
        [0.0,   -s,   c]
    ])

# Test data from our C++ output
# Asteroid position (ecliptic): [-1.150, -2.653, -0.157] AU
# Earth position (ecliptic): [-0.466, 0.866, -0.00005] AU
# Observatory offset (ecliptic): [0, 0, 0] (geocenter for now)

xast = np.array([-1.150, -2.653, -0.157])
xea = np.array([-0.466, 0.866, -0.00005])

print("╔════════════════════════════════════════════════════════╗")
print("║ OrbFit Algorithm Step-by-Step                         ║")
print("╚════════════════════════════════════════════════════════╝\n")

print("INPUT:")
print(f"  Asteroid position (ecliptic): {xast}")
print(f"  Earth position (ecliptic):    {xea}")

# Step 1: Geocentric vector in ECLIPTIC
d_ecl = xast - xea
print(f"\nStep 1: Geocentric vector (ecliptic)")
print(f"  d = xast - xea = {d_ecl}")

# Step 2: Topocentric correction (skip for geocenter)
print(f"\nStep 2: Topocentric correction (geocenter, skip)")

# Step 3: Distance
dis0 = np.linalg.norm(d_ecl)
print(f"\nStep 3: Distance")
print(f"  dis0 = |d| = {dis0:.6f} AU")

# Step 4: Rotation to EQUATORIAL
R = roteceq()
d_eq = R @ d_ecl
print(f"\nStep 4: Rotation to EQUATORIAL")
print(f"  Rotation matrix R_x(+{epsilon * 180/np.pi:.6f}°):")
print(f"    {R[0]}")
print(f"    {R[1]}")
print(f"    {R[2]}")
print(f"  d_eq = R @ d_ecl = {d_eq}")

# Step 5: Normalize
d_norm = d_eq / dis0
print(f"\nStep 5: Normalize")
print(f"  d_norm = d_eq / dis0 = {d_norm}")

# Step 6: Compute RA and Dec
alpha = np.arctan2(d_eq[1], d_eq[0])
if alpha < 0:
    alpha += 2 * np.pi
    
delta = np.arcsin(d_eq[2] / dis0)

print(f"\nStep 6: Compute RA and Dec")
print(f"  RA  = atan2(y, x) = {alpha:.8f} rad = {alpha * 180/np.pi:.6f} deg")
print(f"  Dec = asin(z/r)   = {delta:.8f} rad = {delta * 180/np.pi:.6f} deg")

# Expected values from observation
obs_ra_deg = 20.205458
obs_dec_deg = 11.528833
print(f"\n╔════════════════════════════════════════════════════════╗")
print(f"║ COMPARISON                                             ║")
print(f"╚════════════════════════════════════════════════════════╝")
print(f"  Observed:  RA = {obs_ra_deg:.6f}°, Dec = {obs_dec_deg:.6f}°")
print(f"  Computed:  RA = {alpha * 180/np.pi:.6f}°, Dec = {delta * 180/np.pi:.6f}°")
print(f"  Residual:  ΔRA = {obs_ra_deg - alpha * 180/np.pi:.6f}°, ΔDec = {obs_dec_deg - delta * 180/np.pi:.6f}°")

# Check also if we invert the rotation
print(f"\n╔════════════════════════════════════════════════════════╗")
print(f"║ TEST: What if we use WRONG rotation R_x(-ε)?          ║")
print(f"╚════════════════════════════════════════════════════════╝")
c = np.cos(-epsilon)
s = np.sin(-epsilon)
R_wrong = np.array([
    [1.0,  0.0, 0.0],
    [0.0,    c,   s],
    [0.0,   -s,   c]
])
d_eq_wrong = R_wrong @ d_ecl
alpha_wrong = np.arctan2(d_eq_wrong[1], d_eq_wrong[0])
if alpha_wrong < 0:
    alpha_wrong += 2 * np.pi
delta_wrong = np.arcsin(d_eq_wrong[2] / dis0)

print(f"  RA  = {alpha_wrong * 180/np.pi:.6f}° (wrong rotation)")
print(f"  Dec = {delta_wrong * 180/np.pi:.6f}° (wrong rotation)")
