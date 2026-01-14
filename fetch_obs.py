import requests
import json

def get_horizons_obs(obj_id, epoch_mjd):
    epoch_jd = epoch_mjd + 2400000.5
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    params = {
        'format': 'json',
        'COMMAND': f"'{obj_id}'",
        'OBJ_DATA': 'NO',
        'MAKE_EPHEM': 'YES',
        'EPHEM_TYPE': 'OBSERVER',
        'CENTER': '500',  # Earth Geocentric
        'QUANTITIES': '1', # RA/Dec only
        'TLIST': f"{epoch_jd}",
    }
    
    response = requests.get(url, params=params)
    try:
        data = response.json()
        if 'result' in data:
            lines = data['result'].split('\n')
            for line in lines:
                if "2026-Jan" in line:
                    print(f"[{epoch}] found: {line.strip()}")
        else:
            print(f"Error: {data}")
    except Exception as e:
        print(f"Error: {e}")

# Asteroid 34713
# Asteroid 249 (Ilse)
# Epochs: 2026-Jan-09, 10, 11 -> MJD 61049.0, 61050.0, 61051.0
epochs = [61049.0, 61050.0, 61051.0]

results = []

print("Fetching data...")
for epoch in epochs:
    print(f"--- Fecthing for Epoch {epoch} ---")
    get_horizons_obs('249', epoch)
    # Note: The function prints the result but doesn't return parsed RA/Dec efficiently for the table without major refactoring.
    # I will rely on the printed output which should now be short enough if I didn't truncate it, or I can parse it here.
    # Actually, simpler: just let it run. The truncation was due to long previous output. 
    # I will modify the function to print the RAW LINE with "R.A." only.
