#!/usr/bin/env python3
"""
update_asteroids.py — aggiorna il database degli asteroidi (allnum.db) con un
merge selettivo per fonte, preservando tutto cio' che le fonti non coprono.

Politica:
  - NEODyS (neodys.cat)  -> i NEO: hanno la priorita', sono le orbite curate.
  - AstDyS (allnum.cat)  -> tutti gli altri numerati non gia' presi da NEODyS.
  - Qualsiasi corpo gia' nel DB e non coperto da nessuna delle due fonti
    (es. i non-numerati) viene LASCIATO INTATTO. Nessun DELETE, mai.
  - I corpi presenti vengono aggiornati; i numerati nuovi vengono aggiunti.

Sicurezza: lavora su una COPIA del DB in un file temporaneo, applica gli
upsert, verifica, e solo se tutto torna sostituisce l'originale (swap atomico,
con backup .bak). Un download interrotto non tocca il DB in uso.

Formato dei cataloghi (OEF, single-line, KEP), da wro1lr di OrbFit:
  'name'  epoch(MJD)  a  e  i  node  peri  M   H  G   [flag]
  - a in AU; i, node, peri, M in gradi; elementi in notazione e25.16.
  - H < -100  => magnitudine assente.

Uso tipico (sul sistema con accesso a newton.spacedys.com):
    python3 update_asteroids.py                    # scarica entrambe, aggiorna
    python3 update_asteroids.py --neodys-cat F     # usa file NEODyS locale
    python3 update_asteroids.py --astdys-cat F     # usa file AstDyS locale
    python3 update_asteroids.py --dry-run          # non sostituisce, solo report
"""

import argparse
import os
import shutil
import sqlite3
import sys
import tempfile
import time
import urllib.request
from datetime import datetime, timezone

NEODYS_URL = "https://newton.spacedys.com/~neodys2/neodys.cat"
ASTDYS_URL = "https://newton.spacedys.com/~astdys2/catalogs/allnum.cat"
DEFAULT_DB = os.path.expanduser("~/.ioccultcalc/database/allnum.db")


def log(msg):
    print(f"[update_asteroids] {msg}", flush=True)


def download(url, dest):
    log(f"download {url}")
    start = time.time()
    last = [start]

    def hook(blocks, bsize, total):
        now = time.time()
        if now - last[0] >= 1.0 or (total > 0 and blocks * bsize >= total):
            mb = blocks * bsize / 1e6
            pct = f" ({100.0*blocks*bsize/total:5.1f}%)" if total > 0 else ""
            log(f"  {mb:7.1f} MB{pct}")
            last[0] = now

    urllib.request.urlretrieve(url, dest, reporthook=hook)
    log(f"  scaricato in {time.time()-start:.1f}s ({os.path.getsize(dest)/1e6:.1f} MB)")


def parse_oef(path):
    """
    Generatore di record dal catalogo OEF single-line KEP.
    Produce dict con chiavi: number|None, designation, epoch, a, e, i,
    node, peri, M, H|None, G|None.
    """
    in_body = False
    with open(path, "r", errors="replace") as f:
        for line in f:
            if not in_body:
                # header finisce a END_OF_HEADER; alcune varianti non ce l'hanno
                if line.strip() == "END_OF_HEADER":
                    in_body = True
                    continue
                # se la riga inizia con apice, siamo gia' nel corpo
                if not line.lstrip().startswith("'"):
                    continue
                in_body = True
            s = line.strip()
            if not s or s.startswith("!") or s[0] != "'":
                continue
            end = s.find("'", 1)
            if end < 0:
                continue
            name = s[1:end].strip()
            rest = s[end + 1:].split()
            if len(rest) < 7:
                continue
            try:
                epoch = float(rest[0])
                a = float(rest[1]); e = float(rest[2]); inc = float(rest[3])
                node = float(rest[4]); peri = float(rest[5]); M = float(rest[6])
            except ValueError:
                continue
            H = None
            G = None
            if len(rest) > 7:
                try:
                    h = float(rest[7])
                    H = h if h > -100.0 else None   # H < -100 = assente
                except ValueError:
                    pass
            if len(rest) > 8:
                try:
                    G = float(rest[8])
                except ValueError:
                    pass
            number = int(name) if name.isdigit() else None
            yield {
                "number": number, "designation": name, "epoch": epoch,
                "a": a, "e": e, "i": inc, "node": node, "peri": peri,
                "M": M, "H": H, "G": G,
            }


def key_of(rec):
    """Chiave d'identita': il numero se numerato, altrimenti la designazione."""
    return ("N", rec["number"]) if rec["number"] is not None else ("D", rec["designation"])


def upsert(cur, rec):
    cur.execute(
        """INSERT INTO allnum_asteroids
           (number, designation, name, epoch_mjd, a, e, i, node_long, peri_arg, M, H, G)
           VALUES (:number,:designation,:designation,:epoch,:a,:e,:i,:node,:peri,:M,:H,:G)
           ON CONFLICT(number) DO UPDATE SET
             designation=excluded.designation, epoch_mjd=excluded.epoch_mjd,
             a=excluded.a, e=excluded.e, i=excluded.i, node_long=excluded.node_long,
             peri_arg=excluded.peri_arg, M=excluded.M, H=excluded.H, G=excluded.G,
             updated_at=CURRENT_TIMESTAMP""",
        rec,
    )


def main():
    ap = argparse.ArgumentParser(description="Aggiorna allnum.db (NEODyS + AstDyS, preserva il resto)")
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--neodys-url", default=NEODYS_URL)
    ap.add_argument("--astdys-url", default=ASTDYS_URL)
    ap.add_argument("--neodys-cat", default=None, help="file NEODyS locale (salta download)")
    ap.add_argument("--astdys-cat", default=None, help="file AstDyS locale (salta download)")
    ap.add_argument("--dry-run", action="store_true", help="non sostituisce il DB, solo report")
    ap.add_argument("--keep-cat", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.db):
        log(f"ERRORE: DB non trovato: {args.db}")
        log("Questo tool aggiorna un DB esistente preservandone il contenuto.")
        return 1

    tmpdir = tempfile.mkdtemp(prefix="update_ast_")
    tmp_db = os.path.join(tmpdir, "allnum.db")
    try:
        # 1. copia del DB esistente: cio' che non tocchiamo resta com'e'
        shutil.copy2(args.db, tmp_db)
        con = sqlite3.connect(tmp_db)
        cur = con.cursor()
        before = cur.execute("SELECT COUNT(*) FROM allnum_asteroids").fetchone()[0]
        log(f"DB di partenza: {before} record (copia di lavoro)")

        # 2. NEODyS
        neodys_path = args.neodys_cat
        if not neodys_path:
            neodys_path = os.path.join(tmpdir, "neodys.cat")
            download(args.neodys_url, neodys_path)
        neo_keys = set()
        n_neo = 0
        for rec in parse_oef(neodys_path):
            upsert(cur, rec)
            neo_keys.add(key_of(rec))
            n_neo += 1
        con.commit()
        log(f"NEODyS: {n_neo} NEO applicati")

        # 3. AstDyS, saltando cio' che NEODyS ha gia' coperto
        astdys_path = args.astdys_cat
        if not astdys_path:
            astdys_path = os.path.join(tmpdir, "allnum.cat")
            download(args.astdys_url, astdys_path)
        n_ast = 0
        n_skip = 0
        for rec in parse_oef(astdys_path):
            if key_of(rec) in neo_keys:
                n_skip += 1
                continue
            upsert(cur, rec)
            n_ast += 1
        con.commit()
        log(f"AstDyS: {n_ast} applicati, {n_skip} saltati (gia' NEO)")

        after = cur.execute("SELECT COUNT(*) FROM allnum_asteroids").fetchone()[0]
        # aggiorna metadata
        cur.execute(
            """INSERT INTO allnum_metadata
               (download_date, total_records, file_url, file_format)
               VALUES (?,?,?,?)""",
            (datetime.now(timezone.utc).isoformat(), after,
             f"NEODyS:{args.neodys_url} + AstDyS:{args.astdys_url}", "OEF-KEP"),
        )
        con.commit()

        # 4. verifica
        def h_of(num):
            r = cur.execute("SELECT H FROM allnum_asteroids WHERE number=?", (num,)).fetchone()
            return r[0] if r else None
        ceres = h_of(1)
        checks_ok = (ceres is not None and 3.0 <= ceres <= 4.0)
        con.close()

        log(f"risultato: {before} -> {after} record  (delta {after-before:+d})")
        log(f"verifica Ceres H={ceres}: {'OK' if checks_ok else 'FALLITO'}")

        if after < before:
            log(f"ABORT: il DB si e' ridotto ({before}->{after}); non sostituito.")
            return 1
        if not checks_ok:
            log("ABORT: verifica corpi noti fallita; non sostituito.")
            return 1

        if args.dry_run:
            log(f"DRY-RUN: DB di lavoro pronto in {tmp_db}, originale NON toccato.")
            return 0

        # 5. swap atomico con backup
        backup = args.db + ".bak"
        shutil.copy2(args.db, backup)
        shutil.move(tmp_db, args.db)
        log(f"OK: {args.db} aggiornato ({after} record). Backup: {backup}")
        return 0

    except Exception as ex:
        log(f"ERRORE: {ex}")
        return 1
    finally:
        if not args.keep_cat:
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
