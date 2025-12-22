#!/usr/bin/env python3
"""
Scarica e formatta dati JPL Horizons per asteroide 17030 (Sierks)
Genera tabella di confronto AstDyn vs Chebyshev vs JPL
"""

import urllib.request
import urllib.parse
import json
from datetime import datetime, timedelta
import sys

def get_jpl_horizons_data():
    """
    Scarica dati effemeridi da JPL Horizons via API
    Asteroide 17030 (Sierks)
    Data: 2025-11-25 to 2025-11-30, step: 6 hours
    """
    
    # Parametri Horizons
    params = {
        'COMMAND': '17030',           # Asteroid 17030 Sierks
        'EPHEM_TYPE': 'VECTORS',      # Coordinate vector output
        'OUT_UNITS': 'AU-D',          # AU and AU/day
        'COORD_TYPE': 'GEODETIC',     # Geodetic/ICRF
        'REF_PLANE': 'FRAME',         # FRAME = Equatorial J2000
        'CENTER': '@sun',             # Barycentric Solar System
        'START_TIME': "'2025-11-25 00:00:00'",
        'STOP_TIME': "'2025-11-30 23:59:00'",
        'STEP_SIZE': '6h',            # 6 hour intervals
        'VEC_LABELS': 'NO',
        'VEC_CORR': 'NONE',
        'CAL_FORMAT': 'BOTH',
        'TIME_DIGITS': 'MINUTES',
        'CSV_FORMAT': 'YES',
    }
    
    # URL Horizons API
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    
    try:
        print("📡 Scaricamento dati JPL Horizons...")
        query_string = urllib.parse.urlencode(params)
        full_url = f"{url}?{query_string}"
        
        import ssl
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        
        # Scarica dati
        with urllib.request.urlopen(full_url, context=ctx, timeout=30) as response:
            data = response.read().decode('utf-8')
            
        print("✓ Dati scaricati")
        return data
        
    except Exception as e:
        print(f"✗ Errore download Horizons: {e}")
        print("  Usando dati locali...")
        return None

def parse_horizons_csv(data):
    """Parse CSV output from Horizons"""
    lines = data.split('\n')
    
    results = []
    in_data = False
    
    for line in lines:
        # Cerca inizio dati
        if '$$SOE' in line:
            print("DEBUG: Found $$SOE")
            in_data = True
            continue
        if '$$EOE' in line:
            print("DEBUG: Found $$EOE")
            break
            
        if in_data and line.strip():
            print(f"DEBUG: Data line: {line}")
            # Parse line
            if ',' in line:
                parts = line.split(',')
            else:
                parts = line.split()
                
            if len(parts) >= 8:
                try:
                    # Parse JD from first column
                    jd = float(parts[0])
                    mjd = jd - 2400000.5
                    
                    # Parse date string for display
                    date_str = parts[1].strip()
                    if "A.D." in date_str:
                        date_str = date_str.replace("A.D. ", "")
                    
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    vx = float(parts[5])
                    vy = float(parts[6])
                    vz = float(parts[7])
                    
                    results.append({
                        'date': date_str,
                        'mjd': mjd,
                        'x': x, 'y': y, 'z': z,
                        'vx': vx, 'vy': vy, 'vz': vz
                    })
                except ValueError:
                    continue
    
    return results

def generate_cpp_data(data):
    """Generate C++ vector initialization"""
    
    print("\n📋 Dati JPL Horizons (convertiti in C++):\n")
    print("std::vector<JPLData> data = {")
    
    for entry in data:
        print(f"    // {entry['date']} = MJD {entry['mjd']:.2f}")
        print(f"    {{{entry['mjd']:.2f}, {entry['x']:.5f}, {entry['y']:.5f}, {entry['z']:.5f}, "
              f"{entry['vx']:.5f}, {entry['vy']:.5f}, {entry['vz']:.5f}}},")
    
    print("};")

def generate_markdown_table(data):
    """Generate Markdown table"""
    print("\n\n## 📊 Dati JPL Horizons - Asteroide 17030 Sierks (25-30 Nov 2025)\n")
    print("| Data/Ora | MJD TDB | X (AU) | Y (AU) | Z (AU) | VX (AU/d) | VY (AU/d) | VZ (AU/d) |")
    print("|----------|---------|--------|--------|--------|-----------|-----------|-----------|")
    
    for entry in data:
        print(f"| {entry['date']} | {entry['mjd']:.4f} | {entry['x']:.6f} | {entry['y']:.6f} | {entry['z']:.6f} | "
              f"{entry['vx']:.6f} | {entry['vy']:.6f} | {entry['vz']:.6f} |")

if __name__ == '__main__':
    # Scarica dati
    csv_data = get_jpl_horizons_data()
    
    if csv_data:
        try:
            # Horizons API returns JSON by default now
            import json
            json_data = json.loads(csv_data)
            if 'result' in json_data:
                csv_data = json_data['result']
        except json.JSONDecodeError:
            pass # Assume it's raw text
            
        # Parse
        data = parse_horizons_csv(csv_data)
        
        if data:
            print(f"✓ {len(data)} effemeridi scaricate")
            
            # Genera output
            generate_cpp_data(data)
            generate_markdown_table(data)
            
            # Salva CSV
            with open('jpl_horizons_17030.csv', 'w') as f:
                f.write("Date,X_AU,Y_AU,Z_AU,VX_AU_d,VY_AU_d,VZ_AU_d,MJD_TDB\n")
                for entry in data:
                    f.write(f"{entry['date']},{entry['x']},{entry['y']},{entry['z']},"
                           f"{entry['vx']},{entry['vy']},{entry['vz']},{entry['mjd']:.4f}\n")
            
            print("\n✓ Dati salvati in jpl_horizons_17030.csv")
