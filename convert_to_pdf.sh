#!/bin/bash
# Script per convertire il report HTML in PDF
# Uso: ./convert_to_pdf.sh

echo "📄 Conversione Report Pompeja in PDF"
echo "===================================="
echo ""

HTML_FILE="POMPEJA_COMPARISON_REPORT.html"
PDF_FILE="POMPEJA_COMPARISON_REPORT.pdf"

if [ ! -f "$HTML_FILE" ]; then
    echo "❌ Errore: $HTML_FILE non trovato"
    echo "   Esegui prima: python3 generate_pdf.py"
    exit 1
fi

echo "✓ File HTML trovato: $HTML_FILE"
echo ""

# Prova diversi metodi di conversione

# Metodo 1: wkhtmltopdf (se installato)
if command -v wkhtmltopdf &> /dev/null; then
    echo "📝 Usando wkhtmltopdf..."
    wkhtmltopdf \
        --enable-local-file-access \
        --page-size A4 \
        --margin-top 20mm \
        --margin-bottom 20mm \
        --margin-left 15mm \
        --margin-right 15mm \
        "$HTML_FILE" "$PDF_FILE"
    
    if [ $? -eq 0 ]; then
        echo "✅ PDF generato: $PDF_FILE"
        echo "   Dimensione: $(du -h "$PDF_FILE" | cut -f1)"
        open "$PDF_FILE"
        exit 0
    fi
fi

# Metodo 2: weasyprint (se installato)
if command -v weasyprint &> /dev/null; then
    echo "📝 Usando weasyprint..."
    weasyprint "$HTML_FILE" "$PDF_FILE"
    
    if [ $? -eq 0 ]; then
        echo "✅ PDF generato: $PDF_FILE"
        echo "   Dimensione: $(du -h "$PDF_FILE" | cut -f1)"
        open "$PDF_FILE"
        exit 0
    fi
fi

# Metodo 3: Apertura manuale in browser
echo "⚠️  Tool automatici non disponibili"
echo ""
echo "📖 Conversione manuale:"
echo "   1. Apri il file HTML nel browser"
echo "   2. Premi Cmd+P (o File > Stampa)"
echo "   3. Seleziona 'Salva come PDF' nel menu a tendina"
echo "   4. Salva come: $PDF_FILE"
echo ""

open "$HTML_FILE"

echo "✓ File HTML aperto nel browser"
echo "  Segui le istruzioni sopra per completare la conversione"

exit 1
