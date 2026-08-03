#!/usr/bin/env python3
"""
importa_diametri.py — porta i diametri misurati nel catalogo locale.

Il diametro decide la larghezza dell'ombra e la durata massima dell'evento: dice
all'osservatore quanto e' larga la fascia da cui vedra' l'occultazione e per
quanti secondi. Finora lo stimavamo da H con albedo fissa 0.10, e su
(316) Goberta questo dava 31.9 km contro i 47.9 misurati — un fattore 1.7, cioe'
una durata prevista di 2.27 s invece di 3.40. La corda registrata da F. Garcia
il 13 luglio 2026 dura 3.4 s: il diametro stimato rendeva quell'osservazione
impossibile.

Fonti, in ordine di preferenza:

  1. JPL SBDB — 139 582 oggetti con diametro, albedo e INCERTEZZA, da IRAS,
     WISE/NEOWISE, Akari, radar e occultazioni. E' la compilazione piu' ampia e
     aggiornata disponibile pubblicamente.
       https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=spkid,full_name,diameter,albedo,diameter_sigma,class&sb-kind=a&sb-cdata=%7B%22AND%22%3A%5B%22diameter%7CDF%22%5D%7D

  2. astorb.dat del Lowell Observatory, colonne 60-64: diametro IRAS. Piu'
     povero (2 140 oggetti) ma indipendente, utile come ricaduta.
       https://ftp.lowell.edu/pub/elgb/astorb.dat.gz

  3. Stima da H con albedo 0.10, per tutto il resto.

Il campo `diameter_source` registra la provenienza, cosi' chi legge sa se il
diametro e' misurato o dedotto, e l'export puo' dichiarare un'incertezza
sensata invece di scrivere zero.

Uso:
    importa_diametri.py --jpl diametri.json [--astorb astorb.dat] [--dry-run]
"""

import argparse
import json
import math
import os
import re
import sqlite3
import sys

ALBEDO_PREDEFINITA = 0.10


def diametro_da_H(h, albedo=ALBEDO_PREDEFINITA):
    """Relazione standard: D = 1329 / sqrt(p) * 10^(-H/5), in km."""
    return 1329.0 / math.sqrt(albedo) * 10.0 ** (-h / 5.0)


def leggi_jpl(percorso):
    """Numero -> (diametro, albedo, sigma, classe). Solo asteroidi numerati."""
    with open(percorso) as f:
        d = json.load(f)
    campi = {n: i for i, n in enumerate(d["fields"])}
    fuori = {}
    for r in d["data"]:
        # spkid 20000001 = asteroide numerato 1
        spkid = r[campi["spkid"]]
        try:
            spkid = int(spkid)
        except (TypeError, ValueError):
            continue
        if not (20000000 < spkid < 21000000):
            continue          # non numerato, o cometa
        numero = spkid - 20000000

        def val(nome):
            v = r[campi[nome]] if nome in campi else None
            if v in (None, ""):
                return None
            try:
                return float(v)
            except (TypeError, ValueError):
                return None

        diam = val("diameter")
        if not diam or diam <= 0:
            continue

        # "     1 Ceres (A801 AA)" -> "Ceres";  "  4769 Castalia (1989 PB)" -> "Castalia"
        # Gli asteroidi senza nome proprio danno "(2015 BK290)": nome vuoto.
        nome = None
        grezzo = r[campi["full_name"]] if "full_name" in campi else ""
        m = re.match(r"\s*\d+\s+([^(]+)", grezzo or "")
        if m:
            nome = m.group(1).strip() or None

        fuori[numero] = (diam, val("albedo"), val("diameter_sigma"),
                         r[campi["class"]] if "class" in campi else None,
                         nome)
    return fuori


def leggi_astorb(percorso):
    """Numero -> diametro IRAS in km. Colonne 60-64 del record."""
    fuori = {}
    with open(percorso, "r", errors="replace") as f:
        for riga in f:
            if len(riga) < 70:
                continue
            try:
                numero = int(riga[0:6])
            except ValueError:
                continue
            t = riga[59:65].strip()
            if not t:
                continue
            try:
                d = float(t)
            except ValueError:
                continue
            if d > 0:
                fuori[numero] = d
    return fuori


def main():
    ap = argparse.ArgumentParser(description="Importa i diametri nel catalogo locale")
    ap.add_argument("--jpl", default="diametri.json",
                    help="JSON dalla query SBDB di JPL")
    ap.add_argument("--astorb", default=os.path.expanduser("~/tmp/astorb/astorb.dat"))
    ap.add_argument("--db", default=os.path.expanduser("~/.ioccultcalc/database/allnum.db"))
    ap.add_argument("--dry-run", action="store_true", help="non scrive, riporta soltanto")
    args = ap.parse_args()

    jpl = {}
    if os.path.exists(args.jpl):
        print(f"lettura di {args.jpl} ...")
        jpl = leggi_jpl(args.jpl)
        print(f"  {len(jpl)} asteroidi numerati con diametro misurato")
    else:
        print(f"JPL non trovato ({args.jpl}), si prosegue senza")

    astorb = {}
    if os.path.exists(args.astorb):
        print(f"lettura di {args.astorb} ...")
        astorb = leggi_astorb(args.astorb)
        print(f"  {len(astorb)} con diametro IRAS")

    if not os.path.exists(args.db):
        sys.exit(f"catalogo non trovato: {args.db}")

    conn = sqlite3.connect(args.db)
    cur = conn.cursor()

    esistenti = {r[1] for r in cur.execute("PRAGMA table_info(allnum_asteroids)")}
    for nome, tipo in (("diameter_km", "REAL"),
                       ("diameter_sigma_km", "REAL"),
                       ("diameter_source", "TEXT"),
                       ("albedo", "REAL"),
                       ("taxonomy", "TEXT"),
                       ("orbit_class", "TEXT")):
        if nome not in esistenti:
            if args.dry_run:
                print(f"  [prova] aggiungerei la colonna {nome}")
            else:
                cur.execute(f"ALTER TABLE allnum_asteroids ADD COLUMN {nome} {tipo}")
                print(f"  aggiunta la colonna {nome}")

    righe = list(cur.execute("SELECT number, H, name FROM allnum_asteroids"))
    print(f"\n{len(righe)} asteroidi nel catalogo locale")

    conteggio = {"JPL": 0, "IRAS": 0, "H+albedo 0.10": 0, "senza dati": 0}
    nomi_nuovi = 0
    aggiornamenti = []

    for numero, h, nome_db in righe:
        nome = None
        if numero in jpl:
            diam, albedo, sigma, classe, nome = jpl[numero]
            fonte = "JPL"
        elif numero in astorb:
            diam, sigma, classe = astorb[numero], None, None
            albedo = ((1329.0 / diam) ** 2 * 10.0 ** (-2.0 * h / 5.0)) if h else None
            fonte = "IRAS"
        elif h is not None:
            diam, albedo, sigma, classe = diametro_da_H(h), ALBEDO_PREDEFINITA, None, None
            fonte = "H+albedo 0.10"
        else:
            conteggio["senza dati"] += 1
            continue
        conteggio[fonte] += 1

        # Il nome nel catalogo locale ha la forma "(316) Goberta"; JPL lo da'
        # pulito. Si aggiorna solo se il nostro manca o porta il prefisso.
        nome_finale = nome_db
        if nome and (not nome_db or nome_db.startswith("(")):
            nome_finale = nome
            nomi_nuovi += 1

        aggiornamenti.append((round(diam, 4),
                              round(sigma, 4) if sigma else None,
                              fonte,
                              round(albedo, 4) if albedo else None,
                              classe, nome_finale, numero))

    print("\nprovenienza dei diametri:")
    tot = sum(conteggio.values())
    for k, v in conteggio.items():
        print(f"  {k:16s} {v:8d}  {100.0*v/tot:5.1f}%")

    if args.dry_run:
        print("\n[prova] nessuna scrittura")
        conn.close()
        return

    cur.executemany(
        "UPDATE allnum_asteroids SET diameter_km=?, diameter_sigma_km=?, "
        "diameter_source=?, albedo=?, orbit_class=?, name=? WHERE number=?",
        aggiornamenti)
    conn.commit()
    print(f"\n{len(aggiornamenti)} record aggiornati, {nomi_nuovi} nomi ripuliti")

    print("\nverifica:")
    for n, atteso in ((316, "Occult4 usa 54.6"),
                      (1216, "Askania, riferimento 10.66"),
                      (4769, "Castalia, riferimento 1.40"),
                      (65803, "Didymos, 0.78")):
        r = cur.execute("SELECT name, diameter_km, diameter_sigma_km, diameter_source, "
                        "albedo, orbit_class FROM allnum_asteroids WHERE number=?",
                        (n,)).fetchone()
        if r:
            nome, d, sg, fonte, alb, cls = r
            print(f"  ({n}) {nome or '':10s} {d or 0:9.3f} +/- {sg or 0:6.3f} km  "
                  f"{fonte or '?':14s} albedo {alb or 0:.3f}  {cls or '':5s} [{atteso}]")

    conn.close()


if __name__ == "__main__":
    main()
