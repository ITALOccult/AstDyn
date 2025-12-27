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
    except Exception as e:
        print(f"Error parsing JSON: {e}")
        print(f"Raw response: {response.text}")
        return
    
    if 'result' in data:
        print(data['result'])
    else:
        print(f"Error in API response: {data}")

# Asteroid 34713
get_horizons_obs('34713', 61050.0)
