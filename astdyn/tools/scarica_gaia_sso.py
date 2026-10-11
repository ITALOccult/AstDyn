#!/usr/bin/env python3
"""
scarica_gaia_sso.py — scarica l'astrometria Gaia FPR di un asteroide.

Le osservazioni Gaia (tabella gaiafpr.sso_observation) non sono distribuite
dall'MPC ne' da AstDyS, ma JPL le usa: senza, l'orbita di (3830) restava a
70-85 km da quella JPL sul piano dell'ombra. ioccultcalc le legge da
<observations_dir>/<numero>_gaia_fpr.csv, accanto al file .rwo.

Uso:
    scarica_gaia_sso.py NUMERO CARTELLA_OSSERVAZIONI
"""

import argparse
import pathlib
import sys
import urllib.parse
import urllib.request

TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"


def scarica(numero, cartella):
    query = f"SELECT * FROM gaiafpr.sso_observation WHERE number_mp={int(numero)} ORDER BY epoch_utc"
    dati = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query})
    with urllib.request.urlopen(TAP_SYNC, data=dati.encode(), timeout=300) as risposta:
        testo = risposta.read().decode()
    righe = testo.count("\n") - 1
    destinazione = pathlib.Path(cartella).expanduser() / f"{int(numero)}_gaia_fpr.csv"
    if righe <= 0:
        print(f"({numero}): nessuna osservazione Gaia FPR", file=sys.stderr)
        return None
    destinazione.write_text(testo)
    print(f"({numero}): {righe} osservazioni Gaia FPR in {destinazione}")
    return destinazione


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("numero", type=int)
    ap.add_argument("cartella")
    args = ap.parse_args()
    return 0 if scarica(args.numero, args.cartella) else 1


if __name__ == "__main__":
    sys.exit(main())
