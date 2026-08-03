#!/usr/bin/env python3
"""
verifica_corde.py — confronta una predizione con le corde effettivamente
osservate.

Il confronto fra due predizioni dice quanto due calcoli differiscono; il
confronto con una corda dice quale dei due aveva ragione. La corda registra
dove l'ombra e' passata davvero, e l'istante in cui e' passata sopra un punto
noto: la distanza fra la traccia predetta e quel punto e' l'errore vero, non
una stima.

Da ciascuna corda si ricavano due quantita' indipendenti:

  - l'ISTANTE centrale, che si confronta direttamente con quello predetto per
    quel sito: uno scarto in tempo dice se la predizione anticipa o ritarda;

  - la DURATA, che rapportata alla durata massima da' la lunghezza della corda
    e quindi la distanza del sito dall'asse dell'ombra:
        L = D * t_corda / t_max,      d = sqrt((D/2)^2 - (L/2)^2)
    Il segno di quella distanza non e' determinato da una corda sola, ma con
    due corde su lati diversi si vincola la posizione della traccia.

Uso:
    verifica_corde.py predizione.xml --evento ID --corda "nome,lat,lon,alt,t_centro,durata" ...

Latitudine e longitudine in gradi decimali, quota in metri, istante in UT
decimale, durata in secondi.
"""

import argparse
import math
import sys
import xml.etree.ElementTree as ET

A_TERRA_KM = 6378.137
F_TERRA = 1.0 / 298.257223563
E2 = F_TERRA * (2.0 - F_TERRA)
DEG = math.pi / 180.0


def leggi_predizione(percorso, evento_id=None):
    radice = ET.parse(percorso).getroot()
    for ev in radice.findall(".//Event"):
        def campi(tag):
            n = ev.find(tag)
            return [c.strip() for c in n.text.split(",")] if n is not None and n.text else []
        idd = campi("ID")
        if evento_id and (not idd or idd[0] != evento_id):
            continue
        el, ea, ob = campi("Elements"), campi("Earth"), campi("Object")
        return {
            "sorgente": el[0], "durata_max": float(el[1]),
            "anno": int(el[2]), "mese": int(el[3]), "giorno": int(el[4]),
            "ut": float(el[5]),
            "x": float(el[6]), "y": float(el[7]),
            "dx": float(el[8]), "dy": float(el[9]),
            "d2x": float(el[10]), "d2y": float(el[11]),
            "d3x": float(el[12]), "d3y": float(el[13]),
            "sub_lon": float(ea[0]), "sub_lat": float(ea[1]),
            "oggetto": ob[1], "diametro": float(ob[3]),
            "id": idd[0] if idd else "?",
        }
    return None


def asse(p, t_ore):
    """Posizione dell'asse dell'ombra sul piano fondamentale, in raggi terrestri."""
    t = t_ore
    xi = p["x"] + p["dx"]*t + p["d2x"]*t*t + p["d3x"]*t**3
    eta = p["y"] + p["dy"]*t + p["d2y"]*t*t + p["d3y"]*t**3
    return xi, eta


def posizione_osservatore(lat_deg, lon_deg, alt_m, p, t_ore):
    """Coordinate dell'osservatore sul piano fondamentale, in raggi terrestri.

    Il piano e' perpendicolare alla direzione della stella; la declinazione
    sub-stellare da' l'inclinazione dell'asse, e la longitudine sub-stellare
    ruota con la Terra di 15 gradi l'ora.
    """
    lat = lat_deg * DEG
    # latitudine geocentrica dall'ellissoide
    lat_gc = math.atan((1.0 - E2) * math.tan(lat))
    r = A_TERRA_KM / math.sqrt(1.0 - E2 * math.sin(lat)**2)
    rho = ((r + alt_m/1000.0) * math.cos(lat)) / A_TERRA_KM
    z = ((r * (1.0 - E2) + alt_m/1000.0) * math.sin(lat)) / A_TERRA_KM

    # angolo orario rispetto al punto sub-stellare, che ruota con la Terra
    h = (lon_deg - (p["sub_lon"] - 15.0 * t_ore)) * DEG
    ds = p["sub_lat"] * DEG

    # proiezione sugli assi del piano fondamentale
    xi = rho * math.sin(h)
    eta = z * math.cos(ds) - rho * math.sin(ds) * math.cos(h)
    return xi, eta


def distanza_dall_asse(p, lat, lon, alt, t_ore):
    """Distanza dell'osservatore dall'asse dell'ombra, in km, all'istante dato."""
    ax, ay = asse(p, t_ore)
    ox, oy = posizione_osservatore(lat, lon, alt, p, t_ore)
    return math.hypot(ox - ax, oy - ay) * A_TERRA_KM


def istante_minimo(p, lat, lon, alt, mezza_finestra=0.5):
    """Istante in cui l'osservatore e' piu' vicino all'asse: e' quando vede
    l'occultazione. Ricerca per sezione aurea."""
    a, b = -mezza_finestra, mezza_finestra
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - gr*(b-a), a + gr*(b-a)
    for _ in range(80):
        if distanza_dall_asse(p, lat, lon, alt, c) < distanza_dall_asse(p, lat, lon, alt, d):
            b = d
        else:
            a = c
        c, d = b - gr*(b-a), a + gr*(b-a)
    t = 0.5*(a+b)
    return t, distanza_dall_asse(p, lat, lon, alt, t)


def da_dms(testo):
    """'43 28 21.6' oppure '43.4727' -> gradi decimali."""
    parti = testo.replace(":", " ").split()
    if len(parti) == 1:
        return float(parti[0])
    segno = -1.0 if parti[0].strip().startswith("-") else 1.0
    g = abs(float(parti[0]))
    m = float(parti[1]) if len(parti) > 1 else 0.0
    s = float(parti[2]) if len(parti) > 2 else 0.0
    return segno * (g + m/60.0 + s/3600.0)


def main():
    ap = argparse.ArgumentParser(
        description="Confronta una predizione con le corde osservate")
    ap.add_argument("predizione", help="file occelmnt")
    ap.add_argument("--evento", default=None, help="ID dell'evento")
    ap.add_argument("--corda", action="append", required=True,
                    metavar="NOME,LAT,LON,ALT,T_CENTRO,DURATA",
                    help="lat/lon in gradi decimali o 'g m s'; t_centro in UT "
                         "'hh mm ss.s'; durata in secondi")
    args = ap.parse_args()

    p = leggi_predizione(args.predizione, args.evento)
    if not p:
        sys.exit("evento non trovato nella predizione")

    print(f"=== ({p['oggetto']}) — {p['giorno']:02d}/{p['mese']:02d}/{p['anno']} ===")
    print(f"predizione   {p['sorgente']}")
    print(f"diametro     {p['diametro']:.3f} km")
    print(f"durata max   {p['durata_max']:.3f} s   (sulla linea centrale)")
    print(f"UT centrale  {p['ut']:.6f} h = {int(p['ut']):02d}:"
          f"{int(p['ut']%1*60):02d}:{(p['ut']%1*60)%1*60:05.2f}\n")

    print(f"{'osservatore':<18}{'d predetta':>12}{'d da corda':>12}"
          f"{'t predetto':>13}{'t osservato':>13}{'scarto':>10}")
    print("-" * 78)

    risultati = []
    for c in args.corda:
        campi = [x.strip() for x in c.split(",")]
        nome = campi[0]
        lat, lon = da_dms(campi[1]), da_dms(campi[2])
        alt = float(campi[3])
        t_oss_h = da_dms(campi[4])          # UT in ore
        durata = float(campi[5])

        # dalla predizione: quando e a che distanza dall'asse
        t_pred_rel, d_pred = istante_minimo(p, lat, lon, alt)
        t_pred_h = p["ut"] + t_pred_rel

        # dalla corda: la lunghezza da' la distanza dall'asse
        L = p["diametro"] * durata / p["durata_max"]
        semi = p["diametro"] / 2.0
        d_oss = math.sqrt(max(0.0, semi*semi - (L/2.0)**2))

        dt_s = (t_oss_h - t_pred_h) * 3600.0
        risultati.append((nome, d_pred, d_oss, dt_s))

        def hms(h):
            hh = int(h); mm = int((h-hh)*60); ss = ((h-hh)*60-mm)*60
            return f"{hh:02d}:{mm:02d}:{ss:05.2f}"

        print(f"{nome:<18}{d_pred:9.1f} km{d_oss:9.1f} km"
              f"{hms(t_pred_h):>13}{hms(t_oss_h):>13}{dt_s:+9.1f} s")

    print()
    if len(risultati) >= 2:
        print("Con due corde la posizione dell'ombra e' vincolata: la differenza")
        print("fra distanza predetta e distanza dalla corda, se dello stesso segno")
        print("per entrambe, misura lo spostamento della traccia.")
        scarti = [dp - do for _, dp, do, _ in risultati]
        print(f"scarti: {', '.join(f'{s:+.1f} km' for s in scarti)}")
    dt_medio = sum(r[3] for r in risultati) / len(risultati)
    print(f"\nscarto medio sull'istante: {dt_medio:+.1f} s")
    if abs(dt_medio) > 1.0:
        v_ombra = p["diametro"] / p["durata_max"]
        print(f"  a {v_ombra:.1f} km/s corrisponde a {abs(dt_medio)*v_ombra:.0f} km "
              f"lungo la traccia")


if __name__ == "__main__":
    main()
