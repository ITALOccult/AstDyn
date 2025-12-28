import requests
import json

def get_horizons_elements(obj_id, epoch_mjd):
    epoch_jd = epoch_mjd + 2400000.5
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    params = {
        'format': 'json',
        'COMMAND': f"'{obj_id}'",
        'OBJ_DATA': 'YES',
        'MAKE_EPHEM': 'YES',
        'EPHEM_TYPE': 'ELEMENTS',
        'CENTER': '500@10',  # Sun center
        'REF_SYSTEM': 'ICRF',
        'TLIST': f"{epoch_jd}",
    }
    
    response = requests.get(url, params=params)
    data = response.json()
    
    if 'result' in data:
        print(data['result'])
    else:
        print("Error fetching data")

# Asteroid 249 (Ilse)
get_horizons_elements('249', 61050.5)
