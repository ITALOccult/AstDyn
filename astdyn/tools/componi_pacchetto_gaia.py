#!/usr/bin/env python3
"""
componi_pacchetto_gaia.py — estrae un sottoinsieme distribuibile del catalogo
stellare, limitato in magnitudine.

Il catalogo completo (gaia_dr3_occult_pro.db) contiene 77.952.398 stelle fino a
magnitudine G 16 e pesa 14.7 GB: troppo per una distribuzione ordinaria. La
distribuzione cumulativa e':

    mag <= 14    16.844.164 stelle   (22%)
    mag <= 15    36.909.365 stelle   (47%)
    mag <= 16    77.952.398 stelle  (100%)

Per le occultazioni asteroidali il limite pratico e' 15-16: oltre, il calo di
luce e' difficilmente misurabile con strumentazione amatoriale.

Lo script copia la tabella `stars` filtrata e RICOSTRUISCE l'indice spaziale
R-tree, che e' indispensabile alle ricerche per cono.

Uso:
    python3 componi_pacchetto_gaia.py --mag 15 --out ~/pacchetti/gaia_mag15.db
    python3 componi_pacchetto_gaia.py --mag 14 --stima     # solo stima, non crea
"""

import argparse
import os
import sqlite3
import sys
import time

SORGENTE = os.path.expanduser("~/.catalog/crossreference/gaia_dr3_occult_pro.db")


def umano(byte):
    for unita in ("B", "KB", "MB", "GB"):
        if byte < 1024:
            return f"{byte:.1f} {unita}"
        byte /= 1024
    return f"{byte:.1f} TB"


def main():
    ap = argparse.ArgumentParser(description="Estrae un sottoinsieme del catalogo stellare")
    ap.add_argument("--sorgente", default=SORGENTE)
    ap.add_argument("--mag", type=float, required=True, help="magnitudine G limite")
    ap.add_argument("--out", help="database di destinazione")
    ap.add_argument("--stima", action="store_true", help="solo stima, non crea nulla")
    args = ap.parse_args()

    if not os.path.exists(args.sorgente):
        sys.exit(f"catalogo non trovato: {args.sorgente}")

    src = sqlite3.connect(f"file:{args.sorgente}?mode=ro", uri=True)
    tot_byte = os.path.getsize(args.sorgente)
    tot_righe = src.execute("SELECT COUNT(*) FROM stars").fetchone()[0]
    sel_righe = src.execute("SELECT COUNT(*) FROM stars WHERE mag <= ?", (args.mag,)).fetchone()[0]
    frazione = sel_righe / tot_righe if tot_righe else 0

    print(f"catalogo   : {args.sorgente}")
    print(f"             {tot_righe:,} stelle, {umano(tot_byte)}")
    print(f"selezione  : mag <= {args.mag}")
    print(f"             {sel_righe:,} stelle ({frazione*100:.1f}%)")
    print(f"stima esito: ~{umano(tot_byte * frazione)} (l'indice R-tree e' proporzionale)")

    if args.stima:
        return
    if not args.out:
        sys.exit("serve --out (oppure usare --stima)")

    dest_path = os.path.expanduser(args.out)
    os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)
    if os.path.exists(dest_path):
        sys.exit(f"esiste gia': {dest_path}")

    # spazio disponibile: servono il risultato piu' altrettanto di temporanei
    st = os.statvfs(os.path.dirname(dest_path) or ".")
    libero = st.f_bavail * st.f_frsize
    servono = tot_byte * frazione * 2.2
    print(f"\nspazio     : {umano(libero)} disponibili, ~{umano(servono)} necessari")
    if libero < servono:
        sys.exit("spazio insufficiente: liberare disco o scegliere una magnitudine minore")

    tmp_path = dest_path + ".parziale"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)

    print(f"\ncreo {tmp_path} ...")
    t0 = time.time()
    dst = sqlite3.connect(tmp_path)
    dst.execute("PRAGMA journal_mode=OFF")
    dst.execute("PRAGMA synchronous=OFF")

    # struttura identica all'originale
    dst.executescript("""
        CREATE TABLE stars (
            sid INTEGER PRIMARY KEY,
            ra REAL, dec REAL,
            pmra REAL, pmdec REAL,
            plx REAL, mag REAL,
            ruwe REAL,
            name TEXT, bayer TEXT,
            flam INTEGER, const TEXT,
            hd INTEGER, hip INTEGER,
            sao INTEGER
        );
    """)

    dst.execute("ATTACH DATABASE ? AS src", (args.sorgente,))

    # Copia a blocchi, per poter mostrare l'avanzamento: una INSERT..SELECT unica
    # non da' alcun segnale di vita per decine di minuti.
    print(f"  copio {sel_righe:,} stelle...")
    BLOCCO = 500_000
    copiate = 0
    ultimo_sid = -1
    while True:
        righe = src.execute(
            "SELECT * FROM stars WHERE mag <= ? AND sid > ? ORDER BY sid LIMIT ?",
            (args.mag, ultimo_sid, BLOCCO)).fetchall()
        if not righe:
            break
        dst.executemany(
            "INSERT INTO stars VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", righe)
        dst.commit()
        ultimo_sid = righe[-1][0]
        copiate += len(righe)

        pct = 100.0 * copiate / sel_righe if sel_righe else 100.0
        trascorso = time.time() - t0
        stima = trascorso / copiate * (sel_righe - copiate) if copiate else 0
        print(f"\r  {copiate:,} / {sel_righe:,}  ({pct:5.1f}%)"
              f"  trascorso {trascorso/60:.1f} min"
              f"  mancano ~{stima/60:.1f} min   ", end="", flush=True)
    print()
    print(f"  copiate {copiate:,} righe in {(time.time()-t0)/60:.1f} min")

    print("  ricostruisco l'indice spaziale...")
    t_rt = time.time()
    dst.execute("CREATE VIRTUAL TABLE stars_spatial USING rtree(id, min_ra, max_ra, min_dec, max_dec)")
    # Anche qui a blocchi: l'indice R-tree su decine di milioni di voci e' lento.
    inserite = 0
    ultimo = -1
    while True:
        righe = dst.execute(
            "SELECT sid, ra, ra, dec, dec FROM stars WHERE sid > ? ORDER BY sid LIMIT ?",
            (ultimo, BLOCCO)).fetchall()
        if not righe:
            break
        dst.executemany("INSERT INTO stars_spatial VALUES (?,?,?,?,?)", righe)
        dst.commit()
        ultimo = righe[-1][0]
        inserite += len(righe)
        pct = 100.0 * inserite / copiate if copiate else 100.0
        print(f"\r  indice: {inserite:,} / {copiate:,}  ({pct:5.1f}%)"
              f"  {(time.time()-t_rt)/60:.1f} min   ", end="", flush=True)
    print()

    print("  indici secondari...")
    dst.execute("CREATE INDEX idx_mag ON stars(mag)")
    dst.execute("CREATE INDEX idx_name ON stars(name) WHERE name IS NOT NULL")
    dst.execute("CREATE INDEX idx_hip ON stars(hip) WHERE hip IS NOT NULL")
    dst.commit()
    dst.execute("DETACH DATABASE src")

    # verifica prima di considerarlo valido
    n_stelle = dst.execute("SELECT COUNT(*) FROM stars").fetchone()[0]
    n_spaziale = dst.execute("SELECT COUNT(*) FROM stars_spatial").fetchone()[0]
    dst.close()
    src.close()

    print(f"\nverifica   : {n_stelle:,} stelle, {n_spaziale:,} voci nell'indice spaziale")
    if n_stelle != sel_righe or n_spaziale != n_stelle:
        print("ATTENZIONE: i conteggi non tornano, il file resta come .parziale")
        return

    os.rename(tmp_path, dest_path)
    print(f"fatto      : {dest_path}  {umano(os.path.getsize(dest_path))}  in {time.time()-t0:.0f} s")
    print("\nPer la distribuzione, comprimere con:")
    print(f"    zstd -19 -T0 {dest_path}      (oppure xz -9)")
    print("Le colonne dei cataloghi classici sono nulle per quasi tutte le stelle,")
    print("quindi la compressione dovrebbe rendere parecchio.")


if __name__ == "__main__":
    main()
