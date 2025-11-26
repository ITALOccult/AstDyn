# Makefile for AstDyn Scientific Manual
# Requires: pdflatex, bibtex (TeX Live or similar)

# Main document
MAIN_EN = main
MAIN_IT = main_it

# Output directory
BUILD_DIR = ../../build/manual

# LaTeX compiler
LATEX = pdflatex
BIBTEX = bibtex

# Compiler flags
LATEX_FLAGS = -interaction=nonstopmode -halt-on-error -file-line-error

.PHONY: all clean en it

all: en

# English version
en: $(BUILD_DIR)/$(MAIN_EN).pdf

$(BUILD_DIR)/$(MAIN_EN).pdf: $(MAIN_EN).tex $(wildcard *.tex)
	@echo "Building English manual..."
	@mkdir -p $(BUILD_DIR)
	$(LATEX) $(LATEX_FLAGS) -output-directory=$(BUILD_DIR) $(MAIN_EN).tex
	@# Run twice for references
	$(LATEX) $(LATEX_FLAGS) -output-directory=$(BUILD_DIR) $(MAIN_EN).tex
	@echo "English manual built: $(BUILD_DIR)/$(MAIN_EN).pdf"

# Italian version (to be implemented)
it: $(BUILD_DIR)/$(MAIN_IT).pdf

$(BUILD_DIR)/$(MAIN_IT).pdf: $(MAIN_IT).tex
	@echo "Building Italian manual..."
	@mkdir -p $(BUILD_DIR)
	$(LATEX) $(LATEX_FLAGS) -output-directory=$(BUILD_DIR) $(MAIN_IT).tex
	$(LATEX) $(LATEX_FLAGS) -output-directory=$(BUILD_DIR) $(MAIN_IT).tex
	@echo "Italian manual built: $(BUILD_DIR)/$(MAIN_IT).pdf"

# Clean build artifacts
clean:
	@echo "Cleaning build files..."
	rm -rf $(BUILD_DIR)/*.aux $(BUILD_DIR)/*.log $(BUILD_DIR)/*.toc
	rm -rf $(BUILD_DIR)/*.lof $(BUILD_DIR)/*.lot $(BUILD_DIR)/*.out
	rm -rf $(BUILD_DIR)/*.bbl $(BUILD_DIR)/*.blg
	@echo "Clean complete."

# Clean everything including PDFs
distclean: clean
	rm -f $(BUILD_DIR)/*.pdf
	@echo "All build files removed."

# Quick build (single pass, for testing)
quick:
	@mkdir -p $(BUILD_DIR)
	$(LATEX) $(LATEX_FLAGS) -output-directory=$(BUILD_DIR) $(MAIN_EN).tex

# Help
help:
	@echo "AstDyn Manual Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all        - Build English manual (default)"
	@echo "  en         - Build English manual"
	@echo "  it         - Build Italian manual"
	@echo "  quick      - Quick build (single pass, for testing)"
	@echo "  clean      - Remove auxiliary files"
	@echo "  distclean  - Remove all generated files including PDFs"
	@echo "  help       - Show this help message"
