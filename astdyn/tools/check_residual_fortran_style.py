#!/usr/bin/env python3
"""
Calcolo residui in stile OrbFit/Fortran
Basato sulla logica standard di calcolo residui astrometrici
"""
import numpy as np

# Dati dal debug test
print("="*60)
print("CALCOLO RESIDUI - STILE FORTRAN/ORBFIT")
print("="*60)

# 1. POSIZIONI IN COORDINATE ECLITTICHE J2000
asteroid_ecl = np.array([-1.14994460, -2.65293284, -0.15727587])  # AU
earth_ecl = np.array([-0.46622061, 0.86639855, -0.00004927])      # AU (from PlanetaryEphemeris)

print("\n1. POSIZIONI HELIOCENTRICHE (ECLIPTIC J2000):")
print(f"   Asteroid: {asteroid_ecl}")
print(f"   Earth:    {earth_ecl}")

# 2. OFFSET GEOCENTRICO OSSERVATORIO
# Dal codice: rho_cos_phi, rho_sin_phi, longitude, LST
# L'osservatorio D29 ha coordinate che producono un offset geocentrico
# Il codice calcola: obs_geocentric_equatorial = [rho_cos_phi * cos(LST), rho_cos_phi * sin(LST), rho_sin_phi]
# Poi converte in eclittico

# Dal debug, la posizione finale dell'osservatore era: [-0.46622061, 0.86639855, -0.00004927]
# Questo include Terra + offset geocentrico

observer_ecl = earth_ecl  # Per ora usiamo solo geocentro per semplicità
print(f"   Observer: {observer_ecl} (geocenter)")

# 3. VETTORE TOPOCENTRICO (ECLIPTIC)
rho_ecl = asteroid_ecl - observer_ecl
print(f"\n2. VETTORE TOPOCENTRICO (ECLIPTIC J2000):")
print(f"   ρ = {rho_ecl} AU")
print(f"   |ρ| = {np.linalg.norm(rho_ecl):.6f} AU")

# 4. CONVERSIONE ECLIPTIC -> EQUATORIAL
# ATTENZIONE: Questa è la parte critica!
# Ecliptic to Equatorial: Ruota di +ε attorno all'asse X
epsilon = 23.439291 * np.pi / 180.0

# Matrice di rotazione: R_x(+ε)
R_ecl_to_eq = np.array([
    [1, 0, 0],
    [0, np.cos(epsilon), np.sin(epsilon)],
    [0, -np.sin(epsilon), np.cos(epsilon)]
])

rho_eq = R_ecl_to_eq @ rho_ecl
direction = rho_eq / np.linalg.norm(rho_eq)

print(f"\n3. VETTORE TOPOCENTRICO (EQUATORIAL J2000):")
print(f"   ρ = {rho_eq} AU")
print(f"   direction = {direction}")

# 5. CALCOLO RA/DEC
ra = np.arctan2(direction[1], direction[0])
if ra < 0:
    ra += 2 * np.pi
dec = np.arcsin(direction[2])

ra_deg = ra * 180 / np.pi
dec_deg = dec * 180 / np.pi

print(f"\n4. RA/DEC CALCOLATE:")
print(f"   RA  = {ra_deg:.6f}° = {ra_deg/15:.6f}h")
print(f"   Dec = {dec_deg:.6f}°")

# 6. OSSERVAZIONI
obs_ra_deg = 20.205458
obs_dec_deg = 11.528833

print(f"\n5. RA/DEC OSSERVATE:")
print(f"   RA  = {obs_ra_deg:.6f}°")
print(f"   Dec = {obs_dec_deg:.6f}°")

# 7. RESIDUI
residual_ra = (obs_ra_deg - ra_deg) * 3600  # arcsec
residual_dec = (obs_dec_deg - dec_deg) * 3600  # arcsec

print(f"\n6. RESIDUI (O-C):")
print(f"   RA:  {residual_ra:.1f} arcsec")
print(f"   Dec: {residual_dec:.1f} arcsec")
print(f"   RMS: {np.sqrt(residual_ra**2 + residual_dec**2):.1f} arcsec")

print("\n" + "="*60)
if abs(residual_ra) > 100000:
    print("❌ ERRORE ENORME! C'è un bug fondamentale")
else:
    print("✓ Residui ragionevoli")
print("="*60)
