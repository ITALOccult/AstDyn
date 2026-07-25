# TODO — Tool di aggiornamento dati esterni

Un unico comando che aggiorna tutti i dati esterni da cui dipende il sistema,
con cache locale e verifica di integrita'. Nasce dalla scoperta (2026-07-22) che
il database osservatori conteneva solo 16 codici hardcoded.

## Dati da gestire

### 1. Codici osservatorio MPC  [PRIORITA' ALTA — bug attivo]
- Sorgente: https://www.minorplanetcenter.net/iau/lists/ObsCodes.html
  (pagina HTML, contenuto in <pre>: va ripulita dai tag)
- Destinazione: ~/.ioccultcalc/observatories/ObsCodes.txt
- Formato: colonne fisse — codice 1-3, longitudine est 5-13,
  rho*cos(phi') 14-21, rho*sin(phi') 22-30, nome 31-80.
  E' esattamente il formato letto da ObservatoryDatabase::loadFromMPCFile().
- Aggiornamento: raro (nuovi codici qualche volta l'anno), ma indispensabile.

### 2. Elementi orbitali AstDyS
- Gia' implementato nell'orchestrator (Fase 3): download .eq1 per asteroide,
  con cache e throttle di cortesia.
- Da integrare nel tool unico, con opzione di refresh forzato.

### 3. Osservazioni astrometriche (.rwo)
- Gia' implementato (Fase 3), scaricate per i positivi.

### 4. Effemeridi planetarie
- de440s.bsp (31 MB, 1849-2150) e sb441-n16.bsp (asteroidi perturbatori).
- Sorgente: JPL NAIF. Cambiano di rado; verificare presenza e dimensione.

### 5. Tabelle di debiasing dei cataloghi stellari  [DA VALUTARE]
- AstDyS/OrbFit applicano correzioni di bias di catalogo (Farnocchia et al.).
  Nel .rwo sono gia' tabulate per osservazione (colonne bias RA/Dec).
- Verificare se ci serve la tabella completa o se basta leggere il .rwo.

### 6. Catalogo asteroidi (allnum.db)
- Gia' esiste update_asteroids.py. Da assorbire nel tool unico.

## Requisiti del tool
- Comando unico (es. `ioccultcalc-update` o sottocomando dell'orchestrator).
- Aggiornamento selettivo: --obscodes, --ephemerides, --catalog, --all.
- Cache con timestamp: non riscaricare se recente, salvo --force.
- User-Agent di cortesia e throttle, come gia' fatto per AstDyS.
- Verifica di integrita' dopo il download (righe attese, formato, dimensione).
- Fallback esplicito: se un dato manca, dirlo, NON degradare in silenzio.

## Bug collegato da correggere nella libreria  [PRIORITA' ALTA]
ObservatoryDatabase:
- loadDefaultObservatories() contiene solo 16 codici hardcoded e non risulta
  chiamata da alcun costruttore: il database puo' essere VUOTO.
- ResidualCalculator::get_observer_position() in caso di codice non trovato fa
  `return earth_center_vec;` — ricade sul GEOCENTRO senza alcun avviso,
  perdendo la parallasse topocentrica.
- Effetto misurato: residui medi 6.13 arcsec su Eros (osservato da stazioni
  quasi tutte assenti dalla lista) contro 1.41 su BK290 (stazioni presenti).
  Riferimento AstDyS: 0.42 e 0.17 arcsec.
- Correzione: caricare ObsCodes.txt completo all'inizializzazione; rendere
  esplicito (warning o errore) il fallback al geocentro.
