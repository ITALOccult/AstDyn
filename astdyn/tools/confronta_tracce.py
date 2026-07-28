#!/usr/bin/env python3
"""
confronta_tracce.py — disegna su una stessa mappa le tracce d'ombra di piu'
predizioni della stessa occultazione.

Legge file in formato occelmnt (l'XML che Occult4 e OccultWatcher si scambiano)
e ne ricostruisce il percorso dagli elementi besseliani, senza propagare nulla:
tutta l'informazione geometrica e' gia' nel file.

Serve a confrontare la propria predizione con quelle pubblicate, e a mostrare
quanto si discostano sul terreno — che e' cio' che conta per chi deve decidere
dove andare a osservare.

Uso:
    confronta_tracce.py nostra.xml riferimento.xml -o confronto.png
    confronta_tracce.py *.xml --evento 20260815_4720.9 -o askania.png

Ogni file puo' contenere piu' eventi: con --evento se ne sceglie uno per ID,
altrimenti si prende il primo che tutti i file hanno in comune.
"""

import argparse
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Terra: WGS84
A_TERRA_KM = 6378.137
F_TERRA = 1.0 / 298.257223563
E2 = F_TERRA * (2.0 - F_TERRA)


@dataclass
class Predizione:
    """Una predizione letta da un file occelmnt."""
    file: str
    sorgente: str
    evento_id: str
    # elementi besseliani, raggi terrestri e ore
    x: float; y: float
    dx: float; dy: float
    d2x: float; d2y: float
    d3x: float; d3y: float
    durata_s: float
    anno: int; mese: int; giorno: int
    ut_ca_h: float
    # orientamento della Terra
    substellar_lon: float
    substellar_lat: float
    subsolar_lon: float
    subsolar_lat: float
    # oggetto e stella
    oggetto: str
    diametro_km: float
    stella: str
    # incertezza
    err_maggiore_as: float = 0.0
    err_minore_as: float = 0.0
    err_pa_deg: float = 0.0
    larghezze_percorso: float = 0.0

    def posizione_asse(self, t_ore):
        """Asse dell'ombra sul piano fondamentale, in raggi terrestri.

        Il polinomio del formato include gia' i fattoriali:
            x(t) = x + dx t + d2x t^2 + d3x t^3
        """
        t = np.asarray(t_ore)
        xi = self.x + self.dx * t + self.d2x * t**2 + self.d3x * t**3
        eta = self.y + self.dy * t + self.d2y * t**2 + self.d3y * t**3
        return xi, eta


def _numeri(testo):
    return [c.strip() for c in testo.split(",")]


def leggi_occelmnt(percorso):
    """Legge tutte le predizioni contenute in un file occelmnt."""
    try:
        radice = ET.parse(percorso).getroot()
    except ET.ParseError as e:
        print(f"  {percorso}: XML non valido ({e})", file=sys.stderr)
        return []

    fuori = []
    for ev in radice.findall(".//Event"):
        def campo(tag, idx, default=0.0):
            nodo = ev.find(tag)
            if nodo is None or nodo.text is None:
                return default
            parti = _numeri(nodo.text)
            if idx >= len(parti) or parti[idx] == "":
                return default
            try:
                return float(parti[idx])
            except ValueError:
                return parti[idx]

        def testo(tag, idx, default=""):
            nodo = ev.find(tag)
            if nodo is None or nodo.text is None:
                return default
            parti = _numeri(nodo.text)
            return parti[idx] if idx < len(parti) else default

        p = Predizione(
            file=Path(percorso).name,
            sorgente=testo("Elements", 0, "?"),
            evento_id=testo("ID", 0, "?"),
            durata_s=campo("Elements", 1),
            anno=int(campo("Elements", 2)),
            mese=int(campo("Elements", 3)),
            giorno=int(campo("Elements", 4)),
            ut_ca_h=campo("Elements", 5),
            x=campo("Elements", 6), y=campo("Elements", 7),
            dx=campo("Elements", 8), dy=campo("Elements", 9),
            d2x=campo("Elements", 10), d2y=campo("Elements", 11),
            d3x=campo("Elements", 12), d3y=campo("Elements", 13),
            substellar_lon=campo("Earth", 0),
            substellar_lat=campo("Earth", 1),
            subsolar_lon=campo("Earth", 2),
            subsolar_lat=campo("Earth", 3),
            stella=testo("Star", 0, "?"),
            oggetto=testo("Object", 1, "?"),
            diametro_km=campo("Object", 3),
            larghezze_percorso=campo("Errors", 0),
            err_maggiore_as=campo("Errors", 1),
            err_minore_as=campo("Errors", 2),
            err_pa_deg=campo("Errors", 3),
        )
        fuori.append(p)
    return fuori


def proietta(pred, t_ore, offset_perp_km=0.0):
    """Dall'asse dell'ombra alla superficie: latitudine e longitudine.

    Il piano fondamentale e' perpendicolare alla direzione della stella e passa
    per il centro della Terra. Il punto sub-stellare, dato nel file, e' dove la
    stella sta allo zenit all'istante di massimo avvicinamento; ruotando con la
    Terra si ottiene la longitudine agli altri istanti.
    """
    xi, eta = pred.posizione_asse(t_ore)

    # spostamento perpendicolare alla traccia (bordi dell'ombra)
    if offset_perp_km != 0.0:
        vx = pred.dx + 2*pred.d2x*np.asarray(t_ore) + 3*pred.d3x*np.asarray(t_ore)**2
        vy = pred.dy + 2*pred.d2y*np.asarray(t_ore) + 3*pred.d3y*np.asarray(t_ore)**2
        n = np.hypot(vx, vy)
        n = np.where(n > 0, n, 1.0)
        d = offset_perp_km / A_TERRA_KM
        xi = xi + d * (-vy / n)
        eta = eta + d * (vx / n)

    # intersezione dell'asse (parallelo alla direzione della stella) con
    # l'ellissoide. La declinazione sub-stellare da' l'inclinazione dell'asse.
    d_rad = np.radians(pred.substellar_lat)
    sd, cd = np.sin(d_rad), np.cos(d_rad)

    # coordinate nel sistema del piano fondamentale, schiacciate per l'ellissoide
    eta1 = eta / np.sqrt(1.0 - E2)
    rho2 = 1.0 - xi**2 - eta1**2
    # rho2 <= 0: l'asse dell'ombra passa ACCANTO alla Terra, non c'e' traccia da
    # disegnare. Accade per l'intero evento quando l'occultazione e'
    # geometricamente mancata: ioccultcalc le esporta comunque, perche'
    # l'incertezza potrebbe far cadere l'ombra sulla superficie.
    valido = rho2 > 0.0
    zeta1 = np.sqrt(np.where(valido, rho2, 1.0))

    # dal piano fondamentale al sistema geocentrico ruotato
    zeta = zeta1 * np.sqrt(1.0 - E2)
    z = eta * cd + zeta * sd          # verso il polo
    yq = -eta * sd + zeta * cd        # verso il punto sub-stellare

    lat = np.degrees(np.arctan2(z / np.sqrt(1.0 - E2), np.hypot(xi, yq)))
    # longitudine: angolo orario dal punto sub-stellare, che ruota con la Terra
    lon_locale = np.degrees(np.arctan2(xi, yq))
    lon = pred.substellar_lon + lon_locale - 15.0 * np.asarray(t_ore)
    lon = ((lon + 180.0) % 360.0) - 180.0

    lat = np.where(valido, lat, np.nan)
    lon = np.where(valido, lon, np.nan)
    return lat, lon


def altezza_sole(pred, lat, lon, t_ore):
    """Altezza del Sole sull'orizzonte in ciascun punto della traccia, in gradi.

    Il punto sub-solare — dove il Sole e' allo zenit — e' dato nel file. La
    distanza angolare da quel punto e' la distanza zenitale del Sole, quindi
    l'altezza e' 90 gradi meno quella distanza.

    Serve perche' una traccia lunga migliaia di chilometri attraversa spesso il
    terminatore: sapere che l'evento e' "di notte" al centro non dice nulla su
    dove si possa effettivamente osservare.
    """
    # il punto sub-solare ruota con la Terra: 15 gradi all'ora verso ovest
    slon = pred.subsolar_lon - 15.0 * np.asarray(t_ore)
    slat = pred.subsolar_lat

    la, lo = np.radians(lat), np.radians(lon)
    sa, so = np.radians(slat), np.radians(slon)
    cos_z = (np.sin(la) * np.sin(sa) +
             np.cos(la) * np.cos(sa) * np.cos(lo - so))
    return np.degrees(np.arcsin(np.clip(cos_z, -1.0, 1.0)))


def marcatori(ax, pred, t_ore, lat, lon, colore, passo_min=10.0):
    """Segna gli istanti lungo la traccia: e' cio' che dice all'osservatore a che
    ora l'ombra passa sopra di lui."""
    ore_ca = pred.ut_ca_h
    # istanti tondi rispetto all'ora, non rispetto al massimo avvicinamento
    t_min = t_ore.min() * 60.0
    t_max = t_ore.max() * 60.0
    primo = np.ceil((ore_ca * 60.0 + t_min) / passo_min) * passo_min
    segnati = 0
    for minuti_assoluti in np.arange(primo, ore_ca * 60.0 + t_max, passo_min):
        dt_ore = (minuti_assoluti - ore_ca * 60.0) / 60.0
        i = int(np.argmin(np.abs(t_ore - dt_ore)))
        if np.isnan(lat[i]) or np.isnan(lon[i]):
            continue
        ax.plot(lon[i], lat[i], "o", color=colore, markersize=3.4,
                markeredgecolor="white", markeredgewidth=0.5,
                transform=ccrs.PlateCarree(), zorder=7)
        h = int(minuti_assoluti // 60) % 24
        m = int(minuti_assoluti % 60)
        ax.text(lon[i], lat[i], f" {h:02d}:{m:02d}", fontsize=6.5,
                color=colore, transform=ccrs.PlateCarree(),
                ha="left", va="bottom", zorder=8,
                bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                          edgecolor="none", alpha=0.65))
        segnati += 1
    return segnati


def main():
    ap = argparse.ArgumentParser(
        description="Confronta su mappa le tracce d'ombra di piu' predizioni")
    ap.add_argument("file", nargs="+", help="file occelmnt (.xml)")
    ap.add_argument("-o", "--out", default="confronto.png")
    ap.add_argument("--evento", default=None,
                    help="ID dell'evento; se assente si usa il primo comune")
    ap.add_argument("--ore", type=float, default=1.5,
                    help="semiampiezza della finestra temporale (default 1.5 h)")
    ap.add_argument("--sigma", action="store_true",
                    help="disegna il corridoio di incertezza a 1 sigma (tratteggiato)")
    ap.add_argument("--orari", action="store_true",
                    help="segna gli istanti lungo la traccia")
    ap.add_argument("--passo-min", type=float, default=10.0, dest="passo_min",
                    help="intervallo fra i marcatori orari, minuti (default 10)")
    ap.add_argument("--zoom", nargs=4, type=float, metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
                    help="inquadratura esplicita invece di quella automatica")
    ap.add_argument("--zoom-su", nargs=3, type=float, metavar=("LON", "LAT", "RAGGIO_KM"),
                    dest="zoom_su",
                    help="inquadra attorno a un punto, con raggio in km")
    ap.add_argument("--sole-max", type=float, default=-12.0, dest="sole_max",
                    help="altezza del Sole sotto la quale si considera osservabile "
                         "(default -12, crepuscolo nautico)")
    ap.add_argument("--terminatore", action="store_true",
                    help="ombreggia la parte di Terra in luce")
    ap.add_argument("--titolo", default=None)
    args = ap.parse_args()

    # lettura
    per_file = {}
    for f in args.file:
        pp = leggi_occelmnt(f)
        if pp:
            per_file[f] = pp
            print(f"{Path(f).name}: {len(pp)} eventi")
        else:
            print(f"{Path(f).name}: nessun evento", file=sys.stderr)

    if not per_file:
        sys.exit("nessuna predizione da disegnare")

    # scelta dell'evento
    if args.evento:
        scelto = args.evento
    else:
        insiemi = [set(p.evento_id for p in pp) for pp in per_file.values()]
        comuni = set.intersection(*insiemi) if len(insiemi) > 1 else insiemi[0]
        if not comuni:
            print("\nI file non hanno eventi in comune. ID disponibili:",
                  file=sys.stderr)
            for f, pp in per_file.items():
                print(f"  {Path(f).name}: " +
                      ", ".join(sorted(p.evento_id for p in pp)[:6]), file=sys.stderr)
            sys.exit(1)
        scelto = sorted(comuni)[0]

    predizioni = []
    for pp in per_file.values():
        for p in pp:
            if p.evento_id == scelto:
                predizioni.append(p)

    if not predizioni:
        sys.exit(f"evento {scelto} non trovato")

    rif = predizioni[0]
    print(f"\nevento {scelto}: ({rif.oggetto}) occulta {rif.stella} "
          f"il {rif.giorno}/{rif.mese}/{rif.anno}")
    print(f"{len(predizioni)} predizioni a confronto\n")

    # ---- disegno ----
    fig = plt.figure(figsize=(13, 7.5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#f2efe9")
    ax.add_feature(cfeature.OCEAN, facecolor="#dbe9f4")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor="#999999")
    ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.5)

    colori = ["#c1272d", "#0066b3", "#2e8b45", "#e08a00", "#7a3b9e"]
    t = np.linspace(-args.ore, args.ore, 1200)

    if args.terminatore:
        # zona in luce all'istante di massimo avvicinamento
        gl = np.linspace(-180, 180, 361)
        gb = np.linspace(-90, 90, 181)
        LO, LA = np.meshgrid(gl, gb)
        h = altezza_sole(rif, LA, LO, 0.0)
        ax.contourf(LO, LA, h, levels=[args.sole_max, 90],
                    colors=["#ffd97a"], alpha=0.28,
                    transform=ccrs.PlateCarree(), zorder=1)

    lat_tutte, lon_tutte = [], []
    for i, p in enumerate(predizioni):
        col = colori[i % len(colori)]
        lat, lon = proietta(p, t)

        # spezza dove la traccia esce dalla Terra o supera l'antimeridiano
        salto = np.abs(np.diff(lon)) > 180.0
        lon_p = lon.copy()
        lon_p[:-1][salto] = np.nan

        etichetta = f"{p.sorgente[:34]}"
        if p.larghezze_percorso:
            etichetta += f"  (err {p.larghezze_percorso:.2f} pw)"

        # La traccia si spezza al terminatore: tinta piena dove il Sole e' sotto
        # l'orizzonte, tratteggio dove e' sopra. Una traccia lunga migliaia di
        # chilometri attraversa spesso il confine, e la parte diurna non e'
        # osservabile con strumentazione amatoriale.
        h_sole = altezza_sole(p, lat, lon_p, t)
        notte = h_sole < args.sole_max
        lon_notte = np.where(notte, lon_p, np.nan)
        lon_giorno = np.where(~notte, lon_p, np.nan)

        ax.plot(lon_notte, lat, color=col, linewidth=1.9,
                transform=ccrs.Geodetic(), label=etichetta, zorder=5)
        if np.any(~notte & ~np.isnan(lat)):
            ax.plot(lon_giorno, lat, color=col, linewidth=1.0, linestyle=":",
                    alpha=0.65, transform=ccrs.Geodetic(), zorder=5)

        frazione = float(np.sum(notte & ~np.isnan(lat))) / max(1, int(np.sum(~np.isnan(lat))))
        if frazione < 0.999:
            print(f"       osservabile (Sole sotto {args.sole_max:g} gradi) "
                  f"sul {frazione*100:.0f}% della traccia")

        # Bordi dell'ombra: dove l'occultazione e' effettivamente visibile se la
        # predizione e' esatta.
        if p.diametro_km > 0:
            for segno in (-1, 1):
                lb, lnb = proietta(p, t, segno * p.diametro_km / 2.0)
                s2 = np.abs(np.diff(lnb)) > 180.0
                lnb = lnb.copy(); lnb[:-1][s2] = np.nan
                ax.plot(lnb, lb, color=col, linewidth=0.7, alpha=0.6,
                        transform=ccrs.Geodetic(), zorder=4)

        # Corridoio a 1 sigma: entro questa fascia l'ombra puo' realmente cadere.
        # L'incertezza cross-track e' `path_widths` volte la larghezza dell'ombra,
        # che e' la definizione di quel campo nel formato occelmnt.
        if args.sigma and p.larghezze_percorso > 0 and p.diametro_km > 0:
            sigma_km = p.larghezze_percorso * p.diametro_km
            for segno in (-1, 1):
                ls, lns = proietta(p, t, segno * (p.diametro_km / 2.0 + sigma_km))
                s3 = np.abs(np.diff(lns)) > 180.0
                lns = lns.copy(); lns[:-1][s3] = np.nan
                ax.plot(lns, ls, color=col, linewidth=0.8, linestyle="--",
                        alpha=0.45, transform=ccrs.Geodetic(), zorder=3)
            print(f"       corridoio 1 sigma: +/-{sigma_km:.0f} km "
                  f"({p.larghezze_percorso:.2f} larghezze d'ombra)")

        if args.orari:
            n = marcatori(ax, p, t, lat, lon_p, col, args.passo_min)
            if n == 0:
                print("       nessun marcatore orario sulla parte visibile")

        buoni = ~np.isnan(lat)
        if not buoni.any():
            r_min = float(np.min(np.hypot(*p.posizione_asse(t))))
            print(f"       ombra fuori dalla Terra: minimo {r_min:.3f} raggi "
                  f"terrestri ({(r_min - 1) * 6378:.0f} km dalla superficie)")
            continue
        lat_tutte.append(lat[buoni]); lon_tutte.append(lon[buoni])
        print(f"  {p.file:28s} {p.sorgente[:30]:32s} "
              f"diam {p.diametro_km:6.2f} km  durata {p.durata_s:5.2f} s")

    # inquadratura sulla zona interessata
    if args.zoom:
        ax.set_extent(args.zoom, crs=ccrs.PlateCarree())
    elif args.zoom_su:
        clon, clat, raggio = args.zoom_su
        dlat = raggio / 111.0
        dlon = raggio / (111.0 * max(0.15, np.cos(np.radians(clat))))
        ax.set_extent([clon - dlon, clon + dlon, clat - dlat, clat + dlat],
                      crs=ccrs.PlateCarree())
        ax.plot(clon, clat, "*", color="black", markersize=11,
                markeredgecolor="white", markeredgewidth=0.8,
                transform=ccrs.PlateCarree(), zorder=9)
    elif lat_tutte:
        la = np.concatenate(lat_tutte); lo = np.concatenate(lon_tutte)
        mlat = max(6.0, (la.max() - la.min()) * 0.18)
        mlon = max(8.0, (lo.max() - lo.min()) * 0.18)
        ax.set_extent([lo.min() - mlon, lo.max() + mlon,
                       max(-89, la.min() - mlat), min(89, la.max() + mlat)],
                      crs=ccrs.PlateCarree())

    # Scarto fra le predizioni: e' il numero che conta per chi deve decidere dove
    # posizionarsi. Misurato perpendicolarmente alla traccia di riferimento.
    if len(predizioni) > 1:
        print()
        rifer = predizioni[0]
        lat_r, lon_r = proietta(rifer, t)
        for altra in predizioni[1:]:
            lat_a, lon_a = proietta(altra, t)
            buoni = ~(np.isnan(lat_r) | np.isnan(lat_a))
            if not buoni.any():
                continue
            # Distanza PERPENDICOLARE alla traccia di riferimento: e' di quanto
            # un osservatore dovrebbe spostarsi, non la distanza fra punti
            # contemporanei — che sulle code, dove le tracce divergono in
            # direzione, sovrastima molto.
            La = np.radians(lat_a[buoni]); Lo = np.radians(lon_a[buoni])
            Lr = np.radians(lat_r[buoni]); Or_ = np.radians(lon_r[buoni])
            d_km = np.empty(len(La))
            for k in range(len(La)):
                # distanza dal punto k della traccia A a TUTTA la traccia di
                # riferimento: il minimo e' lo scarto trasversale
                dlat = Lr - La[k]
                dlon = (Or_ - Lo[k]) * np.cos((Lr + La[k]) / 2.0)
                d_km[k] = A_TERRA_KM * np.min(np.hypot(dlat, dlon))
            # Scarto AL CENTRO dell'occultazione: e' il punto in cui la
            # geometria e' definita senza ambiguita', e il numero che dice di
            # quanto le due predizioni collocano diversamente l'ombra.
            i_ca = int(np.argmin(np.abs(t)))
            if not (np.isnan(lat_r[i_ca]) or np.isnan(lat_a[i_ca])):
                dl = np.radians(lat_a[i_ca] - lat_r[i_ca])
                do = np.radians(lon_a[i_ca] - lon_r[i_ca]) * np.cos(
                    np.radians((lat_a[i_ca] + lat_r[i_ca]) / 2.0))
                d_ca = A_TERRA_KM * np.hypot(dl, do)
            else:
                d_ca = float("nan")

            print(f"scarto {rifer.sorgente[:22]} - {altra.sorgente[:22]}:")
            print(f"        al massimo avvicinamento: {d_ca:.0f} km")
            print(f"        trasversale lungo la traccia: medio {d_km.mean():.0f} km, "
                  f"massimo {d_km.max():.0f} km")
            if rifer.diametro_km > 0:
                print(f"        ({d_km.mean() / rifer.diametro_km:.1f} larghezze d'ombra)")

    titolo = args.titolo or (
        f"({rif.oggetto}) occulta {rif.stella}\n"
        f"{rif.giorno:02d}/{rif.mese:02d}/{rif.anno}  "
        f"{int(rif.ut_ca_h):02d}:{int((rif.ut_ca_h % 1) * 60):02d} UT")
    ax.set_title(titolo, fontsize=12, pad=12)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(args.out, dpi=180, bbox_inches="tight")
    print(f"\n{args.out}")


if __name__ == "__main__":
    main()
