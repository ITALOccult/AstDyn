#!/usr/bin/env python3
"""
ioccultcalc_setup.py — installa e aggiorna i dati esterni di ioccultcalc.

Il sistema dipende da dati che non stanno nel repository: codici osservatorio,
effemeridi planetarie, catalogo asteroidi, catalogo stellare. Senza di essi
funziona lo stesso, ma peggio — e non sempre lo dice. L'esempio misurato: senza
il catalogo degli osservatori le osservazioni vengono ridotte dal geocentro e
il fit orbitale peggiora di un fattore sei (RMS 1.59" invece di 0.26").

Politica, ripresa da update_asteroids.py:
  - si scarica in un file temporaneo, si verifica, e solo allora si installa;
  - un download interrotto non tocca i dati in uso;
  - nulla viene cancellato: gli aggiornamenti sostituiscono, con backup .bak;
  - se un dato manca, lo si dice chiaramente. Mai degradare in silenzio.

Uso:
    ioccultcalc_setup.py                  # --check: cosa c'e' e cosa manca
    ioccultcalc_setup.py --obscodes       # codici osservatorio MPC
    ioccultcalc_setup.py --ephemerides    # effemeridi planetarie (31 MB)
    ioccultcalc_setup.py --perturbers     # perturbatori asteroidali (645 MB)
    ioccultcalc_setup.py --catalog        # catalogo asteroidi (delega)
    ioccultcalc_setup.py --essential      # il minimo per un sistema funzionante
    ioccultcalc_setup.py --force ...      # riscarica anche se presente
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
import urllib.request

BASE = os.path.expanduser("~/.ioccultcalc")
UA = "ioccultcalc-setup/1.0 (+https://github.com/ITALOccult/AstDyn)"

# Sorgenti pubbliche stabili
URL_OBSCODES = "https://www.minorplanetcenter.net/iau/lists/ObsCodes.html"
URL_DE440S = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp"
URL_SB441 = "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp"

# Manifesto dei pacchetti distribuiti da noi (catalogo stellare, effemeridi di
# satelliti). GitHub non ospita file di queste dimensioni, quindi l'indirizzo e'
# configurabile: si cambia il manifesto senza toccare il codice.
URL_MANIFESTO = os.environ.get("IOCCULTCALC_PACKAGES_URL", "")


def dice(msg):
    print(msg, flush=True)


def umano(n):
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} TB"


def eta_file(path):
    if not os.path.exists(path):
        return None
    giorni = (time.time() - os.path.getmtime(path)) / 86400
    return giorni


def scarica(url, dest_tmp, atteso_min=0):
    """Scarica in dest_tmp mostrando l'avanzamento. Ritorna True se plausibile."""
    dice(f"    da {url}")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r, open(dest_tmp, "wb") as f:
            totale = int(r.headers.get("Content-Length", 0))
            fatto = 0
            ultimo = 0.0
            while True:
                blocco = r.read(1 << 20)
                if not blocco:
                    break
                f.write(blocco)
                fatto += len(blocco)
                ora = time.time()
                if ora - ultimo > 1.0:
                    if totale:
                        print(f"\r    {umano(fatto)} / {umano(totale)} "
                              f"({100*fatto/totale:.0f}%)", end="", flush=True)
                    else:
                        print(f"\r    {umano(fatto)}", end="", flush=True)
                    ultimo = ora
        print()
    except Exception as e:
        dice(f"    ERRORE: {e}")
        return False

    dim = os.path.getsize(dest_tmp)
    if dim < atteso_min:
        dice(f"    ERRORE: file troppo piccolo ({umano(dim)}), atteso almeno {umano(atteso_min)}")
        dice("    Probabilmente e' una pagina di errore, non il dato richiesto.")
        return False
    return True


def installa(tmp, dest):
    """Sostituzione atomica con backup."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        shutil.move(dest, dest + ".bak")
    shutil.move(tmp, dest)
    dice(f"    installato: {dest}  ({umano(os.path.getsize(dest))})")


# ---------------------------------------------------------------- osservatori

DEST_OBSCODES = os.path.join(BASE, "observatories", "ObsCodes.txt")


def stato_obscodes():
    if not os.path.exists(DEST_OBSCODES):
        return False, ("ASSENTE — le osservazioni saranno ridotte dal GEOCENTRO, "
                       "perdendo la parallasse topocentrica (~4\" a 2.2 AU). "
                       "Misurato: RMS del fit 1.59\" invece di 0.26\".")
    n = sum(1 for _ in open(DEST_OBSCODES, errors="ignore"))
    g = eta_file(DEST_OBSCODES)
    return True, f"{n} righe, aggiornato {g:.0f} giorni fa"


def fai_obscodes(force):
    presente, info = stato_obscodes()
    if presente and not force:
        dice(f"  codici osservatorio: gia' presenti ({info}) — usare --force per riscaricare")
        return True
    dice("  codici osservatorio (MPC)")
    tmp_html = DEST_OBSCODES + ".html.tmp"
    if not scarica(URL_OBSCODES, tmp_html, atteso_min=50_000):
        return False

    # La pagina e' HTML con il contenuto a colonne fisse dentro <pre>.
    import re
    testo = open(tmp_html, errors="ignore").read()
    righe = [re.sub(r"<[^>]*>", "", r) for r in testo.splitlines()]
    righe = [r for r in righe if r.strip()]
    os.remove(tmp_html)

    if len(righe) < 2000:
        dice(f"    ERRORE: solo {len(righe)} righe dopo la ripulitura, attese ~2700")
        return False
    contenuto = "\n".join(righe) + "\n"
    noti = ("500", "691", "F51", "G96")
    mancanti = [c for c in noti if not any(r.startswith(c) for r in righe)]
    if mancanti:
        dice(f"    ERRORE: codici noti assenti ({', '.join(mancanti)}): formato inatteso")
        return False

    tmp = DEST_OBSCODES + ".tmp"
    open(tmp, "w").write(contenuto)
    dice(f"    verificato: {len(righe)} righe, codici noti presenti")
    installa(tmp, DEST_OBSCODES)
    return True


# ---------------------------------------------------------------- effemeridi

DEST_DE440S = os.path.join(BASE, "ephemerides", "de440s.bsp")
DEST_SB441 = os.path.join(BASE, "ephemerides", "sb441-n16.bsp")


def stato_ephem():
    if not os.path.exists(DEST_DE440S):
        return False, "ASSENTE — senza effemeridi planetarie non si propaga nulla."
    return True, f"{umano(os.path.getsize(DEST_DE440S))} (de440s, copre 1849-2150)"


def stato_perturbers():
    if not os.path.exists(DEST_SB441):
        return False, ("assente — i perturbatori asteroidali non verranno applicati. "
                       "AstDyS li include: su archi lunghi l'effetto e' dell'ordine dell'arcsec.")
    return True, umano(os.path.getsize(DEST_SB441))


def fai_ephemerides(force):
    if os.path.exists(DEST_DE440S) and not force:
        dice(f"  effemeridi planetarie: gia' presenti — usare --force per riscaricare")
        return True
    dice("  effemeridi planetarie (JPL NAIF, ~31 MB)")
    tmp = DEST_DE440S + ".tmp"
    if not scarica(URL_DE440S, tmp, atteso_min=20_000_000):
        return False
    if open(tmp, "rb").read(8)[:4] not in (b"DAF/", b"NAIF"):
        dice("    ERRORE: non sembra un kernel SPK")
        return False
    installa(tmp, DEST_DE440S)
    dice("    NOTA: de441.bsp (3 GB) copre un intervallo molto piu' ampio ma e' "
         "sconsigliato su macchine con poca memoria.")
    return True


def fai_perturbers(force):
    if os.path.exists(DEST_SB441) and not force:
        dice("  perturbatori asteroidali: gia' presenti — usare --force per riscaricare")
        return True
    dice("  perturbatori asteroidali (JPL, 645 MB — richiede tempo)")
    tmp = DEST_SB441 + ".tmp"
    if not scarica(URL_SB441, tmp, atteso_min=500_000_000):
        return False
    installa(tmp, DEST_SB441)
    return True


# ------------------------------------------------------------------ catalogo

DEST_DB = os.path.join(BASE, "database", "allnum.db")


def stato_catalogo():
    if not os.path.exists(DEST_DB):
        return False, "ASSENTE — nessun catalogo di asteroidi su cui lavorare."
    g = eta_file(DEST_DB)
    return True, f"{umano(os.path.getsize(DEST_DB))}, aggiornato {g:.0f} giorni fa"


def fai_catalogo(force):
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "update_asteroids.py")
    if not os.path.exists(script):
        dice(f"  catalogo asteroidi: script non trovato ({script})")
        return False
    dice("  catalogo asteroidi (delega a update_asteroids.py)")
    return subprocess.call([sys.executable, script]) == 0


# ------------------------------------------------------------------- stellare

CATALOGO_STELLE = os.path.expanduser("~/.catalog/crossreference/gaia_dr3_occult_pro.db")


def stato_stelle():
    if not os.path.exists(CATALOGO_STELLE):
        return False, ("ASSENTE — senza catalogo stellare non si trovano occultazioni. "
                       "Non e' scaricabile da sorgenti pubbliche: va richiesto il "
                       "pacchetto (vedi --packages).")
    return True, umano(os.path.getsize(CATALOGO_STELLE))


def fai_packages(force):
    if not URL_MANIFESTO:
        dice("  pacchetti: nessun manifesto configurato.")
        dice("    Impostare IOCCULTCALC_PACKAGES_URL con l'indirizzo del manifesto.")
        dice("    Serve per i dati che non hanno sorgente pubblica: catalogo stellare,")
        dice("    effemeridi di satelliti per i sistemi binari.")
        return False
    dice(f"  pacchetti: manifesto da {URL_MANIFESTO}")
    dice("    (scarico dei pacchetti non ancora implementato)")
    return False


# --------------------------------------------------------------------- check

def controlla():
    voci = [
        ("codici osservatorio", stato_obscodes, "--obscodes", True),
        ("effemeridi planetarie", stato_ephem, "--ephemerides", True),
        ("catalogo asteroidi", stato_catalogo, "--catalog", True),
        ("catalogo stellare", stato_stelle, "--packages", True),
        ("perturbatori asteroidali", stato_perturbers, "--perturbers", False),
    ]
    dice(f"Dati in {BASE}\n")
    mancano_essenziali = []
    for nome, fn, opzione, essenziale in voci:
        ok, info = fn()
        segno = "  ok  " if ok else ("MANCA" if essenziale else " --  ")
        dice(f"[{segno}] {nome}")
        dice(f"         {info}")
        if not ok:
            dice(f"         installare con: {opzione}")
            if essenziale:
                mancano_essenziali.append(nome)
        dice("")
    if mancano_essenziali:
        dice(f"Mancano dati essenziali: {', '.join(mancano_essenziali)}")
        dice("Per installare il minimo necessario:  ioccultcalc_setup.py --essential")
    else:
        dice("Tutti i dati essenziali sono presenti.")


def main():
    ap = argparse.ArgumentParser(description="Installa e aggiorna i dati esterni di ioccultcalc")
    ap.add_argument("--check", action="store_true", help="verifica cosa c'e' e cosa manca (default)")
    ap.add_argument("--obscodes", action="store_true", help="codici osservatorio MPC")
    ap.add_argument("--ephemerides", action="store_true", help="effemeridi planetarie (31 MB)")
    ap.add_argument("--perturbers", action="store_true", help="perturbatori asteroidali (645 MB)")
    ap.add_argument("--catalog", action="store_true", help="catalogo asteroidi")
    ap.add_argument("--packages", action="store_true", help="pacchetti distribuiti (catalogo stellare)")
    ap.add_argument("--essential", action="store_true", help="il minimo per un sistema funzionante")
    ap.add_argument("--force", action="store_true", help="riscarica anche se presente")
    args = ap.parse_args()

    azioni = any([args.obscodes, args.ephemerides, args.perturbers,
                  args.catalog, args.packages, args.essential])
    if args.check or not azioni:
        controlla()
        return

    esiti = []
    if args.essential:
        dice("Installazione del minimo necessario:\n")
        esiti.append(fai_obscodes(args.force))
        esiti.append(fai_ephemerides(args.force))
        esiti.append(fai_catalogo(args.force))
    else:
        if args.obscodes:    esiti.append(fai_obscodes(args.force))
        if args.ephemerides: esiti.append(fai_ephemerides(args.force))
        if args.perturbers:  esiti.append(fai_perturbers(args.force))
        if args.catalog:     esiti.append(fai_catalogo(args.force))
        if args.packages:    esiti.append(fai_packages(args.force))

    dice("")
    if all(esiti):
        dice("Fatto.")
    else:
        dice("Alcune operazioni non sono riuscite: rivedere i messaggi sopra.")
        sys.exit(1)


if __name__ == "__main__":
    main()
