#!/bin/bash

# Script per testare la compilazione di ogni capitolo singolarmente

BUILD_DIR="../../build/manual"
LOG_DIR="/tmp/chapter_tests"
mkdir -p "$LOG_DIR"

# Estrai il preambolo dal main_it.tex (fino a \begin{document})
create_test_doc() {
    local chapter_file="$1"
    local test_file="$2"
    
    # Copia preambolo
    sed -n '1,/\\begin{document}/p' main_it.tex > "$test_file"
    
    # Aggiungi frontmatter minimale
    echo "\\frontmatter" >> "$test_file"
    echo "\\tableofcontents" >> "$test_file"
    echo "\\mainmatter" >> "$test_file"
    
    # Includi il capitolo
    echo "\\include{$chapter_file}" >> "$test_file"
    
    # Chiudi documento
    echo "\\end{document}" >> "$test_file"
}

# Testa un singolo capitolo
test_chapter() {
    local chapter_num="$1"
    local chapter_name="$2"
    local chapter_file="${chapter_num}_${chapter_name}"
    
    echo "=========================================="
    echo "Testing Chapter $chapter_num: $chapter_name"
    echo "=========================================="
    
    # Crea documento di test
    local test_file="test_ch${chapter_num}.tex"
    create_test_doc "$chapter_file" "$test_file"
    
    # Compila
    pdflatex -interaction=nonstopmode -output-directory="$BUILD_DIR" "$test_file" \
        > "$LOG_DIR/ch${chapter_num}.log" 2>&1
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo "✓ Chapter $chapter_num: SUCCESS"
        rm -f "$test_file"
        return 0
    else
        echo "✗ Chapter $chapter_num: FAILED"
        echo "  Log saved to: $LOG_DIR/ch${chapter_num}.log"
        
        # Estrai errori principali
        echo "  Main errors:"
        grep "^!" "$LOG_DIR/ch${chapter_num}.log" | head -5 | sed 's/^/    /'
        
        # Mostra la linea con l'errore
        echo "  Context:"
        grep -A 2 "^l\.[0-9]" "$LOG_DIR/ch${chapter_num}.log" | head -6 | sed 's/^/    /'
        
        return 1
    fi
}

# Array di capitoli da testare
declare -a chapters=(
    "01:introduzione"
    "02:sistemi_tempo"
    "03:sistemi_coordinate"
    "04:sistemi_riferimento"
    "05:elementi_orbitali"
    "06:problema_due_corpi"
    "07:perturbazioni"
    "08:integrazione_numerica"
    "09:propagazione_orbite"
    "10:matrice_transizione"
    "11:effemeridi"
    "12:osservazioni"
    "13:orbita_iniziale"
    "14:correzione_differenziale"
    "15:residui"
)

echo "Starting chapter-by-chapter compilation test..."
echo "Logs will be saved to: $LOG_DIR"
echo ""

failed_chapters=()
successful_chapters=()

for chapter in "${chapters[@]}"; do
    IFS=':' read -r num name <<< "$chapter"
    
    if test_chapter "$num" "$name"; then
        successful_chapters+=("$num")
    else
        failed_chapters+=("$num")
        echo ""
        echo "Press Enter to continue to next chapter, or Ctrl+C to stop..."
        read -r
    fi
    echo ""
done

# Riepilogo
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Successful chapters: ${#successful_chapters[@]}"
echo "Failed chapters: ${#failed_chapters[@]}"
echo ""

if [ ${#failed_chapters[@]} -gt 0 ]; then
    echo "Failed: ${failed_chapters[*]}"
    exit 1
else
    echo "All chapters compiled successfully!"
    exit 0
fi
