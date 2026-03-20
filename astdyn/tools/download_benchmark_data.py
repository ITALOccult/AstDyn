import requests
import os
import json
import datetime

def parse_mpc_date(line):
    try:
        # Col 16-32: YYYY MM DD.ddddd
        date_str = line[15:32].strip()
        parts = date_str.split()
        yr = int(parts[0])
        mo = int(parts[1])
        day_frac = float(parts[2])
        day = int(day_frac)
        frac = day_frac - day
        return datetime.datetime(yr, mo, day) + datetime.timedelta(days=frac)
    except:
        return datetime.datetime(1900, 1, 1)

def download_obs(designation, output_path):
    print(f"Downloading observations for {designation}...")
    url = "https://data.minorplanetcenter.net/api/get-obs"
    payload = {"desigs": [designation], "output_format": ["OBS80"]}
    try:
        response = requests.request("GET", url, json=payload, headers={"Content-Type": "application/json"})
        if response.status_code == 200:
            data = response.json()
            if isinstance(data, list) and len(data) > 0: content = data[0].get('OBS80', '')
            elif isinstance(data, dict): content = data.get('OBS80', '')
            else: content = response.text
            # Cleanup non-OBS80 lines (like JSON baggage if any)
            clean_lines = [l for l in content.splitlines() if len(l.strip()) >= 80]
            with open(output_path, 'w') as f: f.write('\n'.join(clean_lines))
            return True
        return False
    except: return False

def extract_subset(full_file, iod_file):
    if not os.path.exists(full_file): return
    with open(full_file, 'r', encoding='latin-1') as f:
        lines = [l for l in f if len(l.strip()) >= 80]
    if len(lines) < 3: return

    # Sort lines by date
    lines.sort(key=parse_mpc_date)

    iod_lines = []
    # Find a 1-10 day arc
    for i in range(len(lines) - 2):
        t1 = parse_mpc_date(lines[i])
        for j in range(i + 1, len(lines) - 1):
            t2 = parse_mpc_date(lines[j])
            if (t2 - t1).total_seconds() < 86400 * 0.5: continue
            for k in range(j + 1, len(lines)):
                t3 = parse_mpc_date(lines[k])
                dt = (t3 - t1).total_seconds() / 86400.0
                if 1.0 < dt < 10.0:
                    iod_lines = [lines[i], lines[j], lines[k]]
                    break
            if iod_lines: break
        if iod_lines: break
    
    if not iod_lines: iod_lines = lines[:3]

    with open(iod_file, 'w') as f:
        for l in iod_lines: f.write(l + '\n')
    print(f"Extracted to {iod_file} (Arc: {(parse_mpc_date(iod_lines[2]) - parse_mpc_date(iod_lines[0])).total_seconds()/86400.0:.2f} d)")

data_dir = "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data"
os.makedirs(data_dir, exist_ok=True)
targets = [{"desig": "99942", "name": "apophis"}, {"desig": "4", "name": "vesta"}, {"desig": "3200", "name": "phaethon"}]
for t in targets:
    full_path = os.path.join(data_dir, f"{t['name']}_mpc_full.txt")
    iod_path = os.path.join(data_dir, f"{t['name']}_iod_3obs.txt")
    if download_obs(t['desig'], full_path): extract_subset(full_path, iod_path)
