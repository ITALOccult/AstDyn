#!/usr/bin/env python3
"""
Debug: verifica che la rotazione equatoriale->eclittica sia corretta
"""

import numpy as np

# Obliquità
epsilon = 23.439291 * np.pi / 180.0

# Matrice da EQUATORIALE a ECLITTICO: R_x(-epsilon)
def equatorial_to_ecliptic():
    c = np.cos(-epsilon)
    s = np.sin(-epsilon)
    return np.array([
        [1.0,  0.0, 0.0],
        [0.0,    c,   s],
        [0.0,   -s,   c]
    ])

# Test: punto sull'equatore celeste (Dec=0, RA=0) dovrebbe rimanere sull'eclittica
eq_point = np.array([1.0, 0.0, 0.0])  # RA=0°, Dec=0°
ecl_point = equatorial_to_ecliptic() @ eq_point

print("Test 1: Punto sull'equatore celeste")
print(f"  Equatoriale: {eq_point}")
print(f"  Eclittico:   {ecl_point}")
print(f"  Atteso:      [1, 0, 0] (l'asse X è comune)")

# Test 2: Polo Nord celeste dovrebbe andare sul piano eclittico inclinato
eq_north = np.array([0.0, 0.0, 1.0])  # Dec=90°
ecl_north = equatorial_to_ecliptic() @ eq_north

print(f"\nTest 2: Polo Nord celeste")
print(f"  Equatoriale: {eq_north}")
print(f"  Eclittico:   {ecl_north}")
print(f"  Latitudine eclittica: {np.arcsin(ecl_north[2]) * 180/np.pi:.2f}°")
print(f"  Atteso: ~66.56° (90° - 23.44°)")

# Test 3: verifica che la matrice sia ortonormale
R = equatorial_to_ecliptic()
RT = R.T
identity = R @ RT
print(f"\nTest 3: Matrice ortonormale?")
print(f"  R @ R^T =")
print(f"    {identity[0]}")
print(f"    {identity[1]}")
print(f"    {identity[2]}")
print(f"  Determinante: {np.linalg.det(R):.10f} (deve essere 1.0)")
