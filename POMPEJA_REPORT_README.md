# Report Confronto Pompeja - Istruzioni

## File Disponibili

- **POMPEJA_COMPARISON_REPORT.md** (21 KB) - Report completo in formato Markdown
- **POMPEJA_COMPARISON_REPORT.html** (29 KB) - Report formattato in HTML
- **POMPEJA_COMPARISON_REPORT.pdf** - PDF da generare manualmente

## Come Generare il PDF

### Metodo 1: Conversione Manuale (Raccomandato)

1. Apri `POMPEJA_COMPARISON_REPORT.html` nel browser (doppio click)
2. Premi `Cmd + P` (Mac) o `Ctrl + P` (Windows/Linux)
3. Seleziona "Salva come PDF" come destinazione
4. Clicca "Salva" e scegli dove salvare il file

### Metodo 2: Script Automatico

```bash
./convert_to_pdf.sh
```

Lo script proverà automaticamente diversi tool (wkhtmltopdf, weasyprint) e in caso di fallimento aprirà l'HTML nel browser.

### Metodo 3: Installare Tool di Conversione

**macOS (con Homebrew):**
```bash
brew install wkhtmltopdf
# oppure
pip3 install weasyprint
```

**Linux:**
```bash
sudo apt-get install wkhtmltopdf
# oppure
pip3 install weasyprint
```

Poi esegui:
```bash
wkhtmltopdf --page-size A4 POMPEJA_COMPARISON_REPORT.html POMPEJA_COMPARISON_REPORT.pdf
```

## Contenuti del Report

Il report contiene un'analisi completa del confronto tra:
- **OrbFit** (sistema di riferimento)
- **AstDyn** (sistema in validazione)
- **JPL Horizons** (riferimento indipendente)

### Sezioni Principali

1. **Executive Summary** - Risultati principali
2. **Dati di Input** - Elementi orbitali e osservazioni
3. **Configurazione AstDyn** - Parametri e modello dinamico
4. **Risultati AstDyn** - Convergenza e elementi fittati
5. **Confronto OrbFit vs AstDyn** - Differenze dettagliate
6. **Confronto JPL Horizons** - Validazione indipendente
7. **Vettori di Stato** - Confronto completo in coordinate
8. **Verifica Bug Fix** - Documentazione correzione matrice rotazione
9. **Conclusioni** - Validazione e raccomandazioni

### Risultati Chiave

✅ **AstDyn converge correttamente**: 8 iterazioni, RMS = 0.658 arcsec  
✅ **Accordo eccellente con OrbFit**: Δa = 578 km, Δe = 0.0006  
✅ **Bug fix verificato**: Miglioramento 105,000× nei residui  

## Visualizzazione del Report

Il report HTML è ottimizzato per:
- ✅ Stampa su carta A4
- ✅ Visualizzazione su schermo
- ✅ Conversione in PDF
- ✅ Tabelle formattate
- ✅ Codice con syntax highlighting
- ✅ Simboli matematici

## Note

- Il report è generato automaticamente dal test `test_pompeja_diffcorr_simple`
- I dati provengono da AstDyS-2 (newton.spacedys.com)
- Il bug della matrice di rotazione è stato completamente risolto
- Tutti i test passano ✅

## Contatti

Per domande o problemi:
- Repository: https://github.com/manvalan/ITALOccultLibrary
- Test: astdyn/tests/test_pompeja_diffcorr_simple.cpp

---

**Generato il:** 25 Novembre 2025  
**Software:** AstDyn v1.0.0 (ITALOccult Library)
