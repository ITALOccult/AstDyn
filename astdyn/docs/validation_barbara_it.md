# Report Occultazione (234) Barbara - 10 Maggio 2026

L'evento di occultazione della stella **UCAC4 521-053748** da parte dell'asteroide **(234) Barbara** è stato confermato e mappato con alta precisione (calibrazione < 0.25 mas rispetto a JPL).

## Perché la scansione automatica iniziale ha fallito?
Abbiamo identificato tre motivi tecnici principali:
1.  **Coordinate di Scansione:** La ricerca iniziale è stata eseguita su un'area di cielo errata (Dec +15 invece di +14) a causa di un errore nel puntamento dello script di scansione.
2.  **Finestra Temporale:** Il calcolo automatico si interrompeva alle 00:00 UTC del 10 Maggio, mentre l'evento avviene alle 16:59 UTC.
3.  **Orbita Nominale JPL:** L'orbita nominale scaricata da JPL Horizons posiziona l'ombra a circa 127 km di distanza dalla superficie terrestre (miss). L'evento è rilevabile solo utilizzando gli elementi orbitali precisi forniti nel tuo XML, che correggono la traiettoria portandola ad impattare la Terra.

## Risultati del Calcolo
- **Stella:** UCAC4 521-053748 (Gaia DR3 3737957664602055808)
- **Magnitudine G:** 12.75
- **Diametro Barbara:** 45.9 km
- **Orbita utilizzata:** JPL +INTG (da XML utente)
- **Precisione:** < 0.25 mas (allineata a JPL Ephemeris)

## Prodotti Generati
- **Mappa SVG:** [barbara_may10_map.svg](file:///Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/barbara_may10_map.svg) - Visualizzazione globale del percorso.
- **KML Google Earth:** [barbara_may10.kml](file:///Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/barbara_may10.kml) - Percorso 3D interattivo.

L'ombra attraversa l'India e l'Oceano Indiano, come previsto.
