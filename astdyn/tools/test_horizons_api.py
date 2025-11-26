#!/usr/bin/env python3
"""Test JPL Horizons API to see actual response format"""

import requests

url = 'https://ssd.jpl.nasa.gov/api/horizons.api'

# Test 1: Orbital elements
print("="*70)
print("TEST 1: Orbital Elements")
print("="*70)

params = {
    'format': 'text',
    'COMMAND': "'203'",
    'OBJ_DATA': 'YES',
    'MAKE_EPHEM': 'YES',
    'EPHEM_TYPE': 'ELEMENTS',
    'CENTER': '@10',  # Sun
    'START_TIME': "'2027-01-01'",
    'STOP_TIME': "'2027-01-01'",
    'STEP_SIZE': '1d',
    'REF_PLANE': 'ECLIPTIC',
    'REF_SYSTEM': 'J2000',
    'OUT_UNITS': 'AU-D',
}

response = requests.get(url, params=params)
print(f"Status: {response.status_code}\n")

if response.status_code == 200:
    text = response.text
    
    # Find the elements section (between $$SOE and $$EOE)
    in_data = False
    print("DATA SECTION:")
    for line in text.split('\n'):
        if '$$SOE' in line:
            in_data = True
            continue
        if '$$EOE' in line:
            break
        if in_data:
            print(line)
    
    print("\n" + "="*70)
    print("FULL RESPONSE (first 2000 chars):")
    print("="*70)
    print(text[:2000])
else:
    print(f"Error: {response.status_code}")
    print(response.text[:500])

# Test 2: State vectors
print("\n" + "="*70)
print("TEST 2: State Vectors")
print("="*70)

params2 = {
    'format': 'text',
    'COMMAND': "'203'",
    'OBJ_DATA': 'NO',
    'MAKE_EPHEM': 'YES',
    'EPHEM_TYPE': 'VECTORS',
    'CENTER': '@10',  # Heliocentric
    'START_TIME': "'2027-01-01'",
    'STOP_TIME': "'2027-01-01'",
    'STEP_SIZE': '1d',
    'VEC_TABLE': '2',
    'REF_PLANE': 'ECLIPTIC',
    'REF_SYSTEM': 'J2000',
    'OUT_UNITS': 'AU-D',
    'VEC_CORR': 'NONE',
}

response2 = requests.get(url, params=params2)
print(f"Status: {response2.status_code}\n")

if response2.status_code == 200:
    text2 = response2.text
    
    in_data = False
    print("VECTOR DATA:")
    for line in text2.split('\n'):
        if '$$SOE' in line:
            in_data = True
            continue
        if '$$EOE' in line:
            break
        if in_data:
            print(line)
