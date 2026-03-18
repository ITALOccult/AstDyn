import requests
import json
import numpy as np

def fetch_horizons(command, center, jd_start, jd_end):
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    params = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": f"'{center}'",
        "START_TIME": f"'JD{jd_start}'",
        "STOP_TIME": f"'JD{jd_end}'",
        "STEP_SIZE": "'1d'",
        "CSV_FORMAT": "'NO'",
        "REF_PLANE": "'ECLIPTIC'",
        "REF_SYSTEM": "'J2000'",
        "OUT_UNITS": "'KM-S'",
        "VEC_TABLE": "'3'"
    }
    
    response = requests.get(url, params=params)
    data = response.json()
    if 'result' not in data:
        print(f"Error fetching {command}: {data}")
        return None
    
    result = data['result']
    # Look for $$SOE and $$EOE
    start = result.find("$$SOE")
    end = result.find("$$EOE")
    if start == -1 or end == -1:
        print(f"Data not found for {command}")
        return None
    
    block = result[start+5:end].strip()
    return block

def parse_block(block):
    lines = block.split('\n')
    # Each block is usually:
    # JD...
    # X = ... Y = ... Z = ...
    # VX= ... VY= ... VZ= ...
    # LT= ... RG= ... RR= ...
    
    # Let's extract X, Y, Z, VX, VY, VZ
    res = {}
    for line in lines:
        if "=" in line:
            parts = line.split("=")
            for i in range(len(parts)-1):
                key = parts[i].strip().split()[-1]
                val_str = parts[i+1].strip().split()[0]
                try:
                    res[key] = float(val_str)
                except:
                    pass
    return res

jd = 2461164.5 # 2026-05-04 00:00 TDB (from haumea_mapper.cpp)

# 1. Haumea Barycenter (w.r.t Sun)
haumea_bary = fetch_horizons("20136108", "500@10", jd, jd+0.1)
h_bary_data = parse_block(haumea_bary) if haumea_bary else None

# 2. Hi'iaka (w.r.t Haumea Barycenter)
hiiaka = fetch_horizons("120136108", "500@20136108", jd, jd+0.1)
hi_data = parse_block(hiiaka) if hiiaka else None

# 3. Namaka (w.r.t Haumea Barycenter)
namaka = fetch_horizons("220136108", "500@20136108", jd, jd+0.1)
nm_data = parse_block(namaka) if namaka else None

print("\n=== HAUMEA SYSTEM ELEMENTS (JD 2461164.5 Ecliptic J2000) ===\n")
if h_bary_data:
    print(f"Haumea Bary (Heliocentric):")
    print(f"  X  = {h_bary_data.get('X'):+18.15e} km, Y  = {h_bary_data.get('Y'):+18.15e} km, Z  = {h_bary_data.get('Z'):+18.15e} km")
    print(f"  VX = {h_bary_data.get('VX'):+18.15e} km/s, VY = {h_bary_data.get('VY'):+18.15e} km/s, VZ = {h_bary_data.get('VZ'):+18.15e} km/s")

if hi_data:
    print(f"Hi'iaka (Barycentric):")
    print(f"  X  = {hi_data.get('X'):+18.15e} km, Y  = {hi_data.get('Y'):+18.15e} km, Z  = {hi_data.get('Z'):+18.15e} km")
    print(f"  VX = {hi_data.get('VX'):+18.15e} km/s, VY = {hi_data.get('VY'):+18.15e} km/s, VZ = {hi_data.get('VZ'):+18.15e} km/s")

if nm_data:
    print(f"Namaka (Barycentric):")
    print(f"  X  = {nm_data.get('X'):+18.15e} km, Y  = {nm_data.get('Y'):+18.15e} km, Z  = {nm_data.get('Z'):+18.15e} km")
    print(f"  VX = {nm_data.get('VX'):+18.15e} km/s, VY = {nm_data.get('VY'):+18.15e} km/s, VZ = {nm_data.get('VZ'):+18.15e} km/s")
