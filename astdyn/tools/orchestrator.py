#!/usr/bin/env python3
"""
orchestrator.py - Orchestratore di campagne ITALOccult (Fase 1)

Sta DAVANTI a ioccultcalc. Legge una config di campagna, prepara gli elementi
orbitali (.eq1) in astdyn_base/elements/, invoca ioccultcalc e raccoglie i
risultati in astdyn_base/results/.

Filosofia: Python e' un corriere di file grezzi + orchestrazione. La conversione
KEP->EQU e' una riparametrizzazione algebrica NELLO STESSO frame (ECLIPTIC_J2000),
validata contro l'oracolo BK290 al livello del rumore aritmetico. Python NON
cambia mai frame: verifica che il frame della sorgente sia ECLIPTIC_J2000 e si
ferma se non lo e'.

Fonti (config: source):
  - db            : elementi dal DB locale allnum.db (default; richiede 'database')
  - user          : gli .eq1 sono gia' in astdyn_base/elements/, non li tocchiamo
  - jpl-horizons  : non prepara .eq1; ioccultcalc usa Horizons (comportamento originale)
  - astdys-neodys : predisposto (Fase 3, download online); per ora si comporta come 'db'

Config (YAML o JSON):
  astdyn_base: /percorso/campagna      # OBBLIGATORIA (nessun default)
  source: db                           # db | user | jpl-horizons | astdys-neodys
  database: ~/.ioccultcalc/database/allnum.db   # richiesto se source usa il DB
  objects:
    asteroids: "1,4,704,820987"        # lista esplicita (Fase 1)
  time:
    start: "2026-07-27"
    duration_days: 2
  mag_limit: 12
  ephemeris_file: ~/.ioccultcalc/ephemerides/de440s.bsp
  ioccultcalc: build/tools/bin/ioccultcalc   # path eseguibile (default: cerca in PATH)
"""

import argparse
import copy
import json
import math
import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Conversione KEP -> EQU (stesso frame ECLIPTIC_J2000; validata su BK290)
# ---------------------------------------------------------------------------
def kep_to_equ(a, e, i_deg, node_deg, peri_deg, M_deg):
    """Kepleriani -> equinoziali AstDyS. Angoli in gradi in ingresso.
    Ritorna (a, h, k, p, q, lambda) con a in AU, lambda in gradi."""
    i = math.radians(i_deg)
    node = math.radians(node_deg)
    peri = math.radians(peri_deg)
    M = math.radians(M_deg)
    varpi = node + peri                       # longitudine del perielio
    h = e * math.sin(varpi)
    k = e * math.cos(varpi)
    thi = math.tan(i / 2.0)
    p = thi * math.sin(node)
    q = thi * math.cos(node)
    lam = math.degrees(M + varpi) % 360.0     # mean longitude, gradi
    return a, h, k, p, q, lam


# ---------------------------------------------------------------------------
# Scrittura .eq1 (formato OEF2.0, come AstDyS; senza covarianza)
# ---------------------------------------------------------------------------
def write_eq1(path, number, a, h, k, p, q, lam, epoch_mjd, H=None, G=None):
    """Scrive un .eq1 nel formato che read_eq1 digerisce.
    Solo elementi (niente COV/NOR): il DB non ha covarianza."""
    lines = []
    lines.append("format  = 'OEF2.0'       ! file format")
    lines.append("rectype = '1L'           ! record type (1L/ML)")
    lines.append("refsys  = ECLM J2000     ! default reference system")
    lines.append("END_OF_HEADER")
    lines.append(str(number))
    lines.append("! Equinoctial elements: a, e*sin(LP), e*cos(LP), "
                 "tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.")
    # riga EQU: a in notazione E, gli altri in decimale, lambda in gradi
    lines.append(" EQU   {:.16E}  {:.15f}   {:.15f}   {:.15f}   {:.15f} {:.13f}".format(
        a, h, k, p, q, lam))
    lines.append(" MJD     {:.9f} TDT".format(epoch_mjd))
    if H is not None:
        g_val = G if G is not None else 0.15
        lines.append(" MAG  {:.3f}  {:.3f}".format(H, g_val))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Sorgente DB: estrae kepleriani, verifica frame, converte, scrive .eq1
# ---------------------------------------------------------------------------
def prepare_from_db(db_path, ids, elements_dir):
    """Per ogni id: query DB -> verifica frame -> KEP->EQU -> scrive .eq1.
    Ritorna (preparati, mancanti)."""
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    cur = con.cursor()
    preparati, mancanti = [], []
    for aid in ids:
        row = cur.execute(
            "SELECT number, epoch_mjd, a, e, i, node_long, peri_arg, M, H, G, "
            "frame_type, element_type FROM allnum_asteroids WHERE number=?",
            (int(aid),),
        ).fetchone()
        if row is None:
            mancanti.append(aid)
            continue
        # GUARDIA FRAME: Python non converte mai frame. Se non e'
        # ECLIPTIC_J2000, ci fermiamo (l'.eq1 dichiara ECLM J2000).
        frame = (row["frame_type"] or "").upper()
        if frame not in ("ECLIPTIC_J2000", "ECLM_J2000", "ECLMJ2000"):
            con.close()
            raise SystemExit(
                f"ERRORE: asteroide {aid} ha frame '{row['frame_type']}', "
                f"atteso ECLIPTIC_J2000. Interrotto per non scrivere un .eq1 "
                f"nel frame sbagliato."
            )
        a, h, k, p, q, lam = kep_to_equ(
            row["a"], row["e"], row["i"],
            row["node_long"], row["peri_arg"], row["M"],
        )
        out = Path(elements_dir) / f"{aid}.eq1"
        write_eq1(out, row["number"] if row["number"] is not None else aid,
                  a, h, k, p, q, lam, row["epoch_mjd"], row["H"], row["G"])
        preparati.append(aid)
    con.close()
    return preparati, mancanti


# ---------------------------------------------------------------------------
# Config loader (YAML o JSON)
# ---------------------------------------------------------------------------
def load_config(path):
    text = Path(path).read_text()
    if path.endswith((".yaml", ".yml")):
        import yaml
        return yaml.safe_load(text)
    return json.loads(text)


def expand_path(p):
    return os.path.expanduser(os.path.expandvars(p)) if p else p


def parse_id_list(spec):
    """'1,4,704' o '1-5' -> lista di stringhe. (Fase 1: liste e range semplici.)"""
    ids = []
    for tok in str(spec).split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok and not tok.startswith("-"):
            lo, hi = tok.split("-", 1)
            ids.extend(str(n) for n in range(int(lo), int(hi) + 1))
        else:
            ids.append(tok)
    return ids


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Orchestratore campagne ITALOccult (Fase 1)")
    ap.add_argument("config", help="config di campagna (YAML o JSON)")
    ap.add_argument("--dry-run", action="store_true",
                    help="prepara gli .eq1 e la config del motore, ma non invoca ioccultcalc")
    args = ap.parse_args()

    cfg = load_config(args.config)

    # astdyn_base OBBLIGATORIA
    base = cfg.get("astdyn_base")
    if not base:
        sys.exit("ERRORE: 'astdyn_base' e' obbligatoria nella config.")
    base = Path(expand_path(base))

    source = cfg.get("source", "db")
    objects = cfg.get("objects", {})
    spec = objects.get("asteroids", "")
    if not spec:
        sys.exit("ERRORE: 'objects.asteroids' vuoto (Fase 1: lista esplicita).")
    ids = parse_id_list(spec)

    # struttura cartelle
    elements_dir = base / "elements"
    results = base / "results"
    for d in (elements_dir, results / "xml", results / "map", results / "log"):
        d.mkdir(parents=True, exist_ok=True)

    print(f"[orchestrator] campagna: {base}")
    print(f"[orchestrator] sorgente: {source}")
    print(f"[orchestrator] oggetti : {len(ids)} ({spec})")

    # --- preparazione elementi secondo la fonte ---
    use_elements_dir = True
    if source in ("db", "astdys", "astdys-neodys"):
        # astdys-neodys: per ora dal DB (che gia' contiene il merge). Fase 3 fara'
        # il download online fresco.
        db = cfg.get("database")
        if not db:
            sys.exit("ERRORE: source='{}' richiede 'database' (path allnum.db).".format(source))
        db = expand_path(db)
        if not os.path.exists(db):
            sys.exit(f"ERRORE: DB non trovato: {db}")
        prep, missing = prepare_from_db(db, ids, elements_dir)
        print(f"[orchestrator] elementi preparati: {len(prep)}")
        if missing:
            print(f"[orchestrator] NON nel DB ({len(missing)}): {','.join(missing)}")
    elif source == "user":
        # gli .eq1 sono gia' in elements/: verifichiamo che ci siano
        prep = [aid for aid in ids if (elements_dir / f"{aid}.eq1").exists()]
        missing = [aid for aid in ids if aid not in prep]
        print(f"[orchestrator] .eq1 utente trovati: {len(prep)}")
        if missing:
            print(f"[orchestrator] .eq1 MANCANTI ({len(missing)}): {','.join(missing)}")
    elif source == "jpl-horizons":
        # non prepariamo nulla: ioccultcalc usera' Horizons (no elements_dir)
        use_elements_dir = False
        prep = ids
        print("[orchestrator] fonte Horizons: elementi presi dal motore (online).")
    else:
        sys.exit(f"ERRORE: sorgente sconosciuta '{source}'. "
                 f"Usa: db | user | jpl-horizons | astdys-neodys")

    if not prep and source != "jpl-horizons":
        sys.exit("ERRORE: nessun elemento preparato. Niente da processare.")

    # --- config per il motore: PARTIAMO DALLA CONFIG DELL'UTENTE e la
    # arricchiamo. Un solo file di verita': l'utente scrive UN file, qui
    # aggiungiamo solo le chiavi che il motore ricava dalla campagna
    # (elements_dir, out-dir, ...). ioccultcalc ignora le chiavi che non
    # conosce (astdyn_base, source, database), quindi possono restare.
    engine_cfg = copy.deepcopy(cfg)

    # la lista effettivamente preparata (i corpi trovati), come stringa
    engine_cfg.setdefault("objects", {})
    engine_cfg["objects"]["asteroids"] = ",".join(prep) if prep else spec
    if use_elements_dir:
        engine_cfg["objects"]["elements_dir"] = str(elements_dir)
    if cfg.get("ephemeris_file"):
        engine_cfg["ephemeris_file"] = expand_path(cfg["ephemeris_file"])

    # Output: dirigi XML e mappa nelle cartelle results/. Sovrascrivibili da
    # config (blocco output: {xml, svg, write_empty}).
    out_cfg = cfg.get("output", {}) or {}
    engine_cfg["out-dir"] = str(results)
    engine_cfg["xml-output"] = out_cfg.get("xml", "xml/occultations.xml")
    if out_cfg.get("svg", True):
        svg_name = out_cfg["svg"] if isinstance(out_cfg.get("svg"), str) else "map/worldmap.svg"
        engine_cfg["svg-output"] = svg_name
    engine_cfg["output"] = {"write_empty": bool(out_cfg.get("write_empty", True))}

    # scriviamo la config arricchita accanto alla base (per ispezione/riproducibilita')
    engine_cfg_path = base / "engine_config.json"
    engine_cfg_path.write_text(json.dumps(engine_cfg, indent=2))
    print(f"[orchestrator] config motore (arricchita): {engine_cfg_path}")

    if args.dry_run:
        print("[orchestrator] --dry-run: preparazione completata, motore NON invocato.")
        return

    # --- invocazione ioccultcalc ---
    exe = expand_path(cfg.get("ioccultcalc", "ioccultcalc"))
    cmd = [exe, "--conf", str(engine_cfg_path)]
    if cfg.get("mag_limit") is not None:
        cmd += ["--mag", str(cfg["mag_limit"])]

    env = dict(os.environ)
    if cfg.get("ephemeris_file"):
        env["ASTDYN_EPHEMERIS_PATH"] = expand_path(cfg["ephemeris_file"])

    log_path = results / "log" / "campaign.log"
    print(f"[orchestrator] invoco: {' '.join(cmd)}")
    print(f"[orchestrator] log: {log_path}")

    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True)
        logf.write(proc.stdout)

    # riepilogo: righe salienti dal log
    salienti = [ln for ln in proc.stdout.splitlines()
                if any(w in ln for w in ("Trovate", "ellisse", "caricato", "SKIP"))]
    print("[orchestrator] --- risultati ---")
    for ln in salienti[:40]:
        print("  " + ln)
    print(f"[orchestrator] fine (exit {proc.returncode}). Risultati in {results}/")


if __name__ == "__main__":
    main()
