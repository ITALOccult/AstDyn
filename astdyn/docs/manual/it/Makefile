# Makefile per il manuale italiano di AstDyn

MAIN = main_it
BUILD_DIR = ../../build/manual
TEX = pdflatex
BIBTEX = bibtex

# Target principale
all: $(BUILD_DIR)/$(MAIN).pdf

# Compilazione completa (con bibliografia)
$(BUILD_DIR)/$(MAIN).pdf: $(MAIN).tex
	$(TEX) -output-directory=$(BUILD_DIR) $(MAIN).tex
	-$(BIBTEX) $(BUILD_DIR)/$(MAIN)
	$(TEX) -output-directory=$(BUILD_DIR) $(MAIN).tex
	$(TEX) -output-directory=$(BUILD_DIR) $(MAIN).tex

# Compilazione veloce (senza bibliografia)
quick:
	$(TEX) -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(MAIN).tex

# Pulizia file ausiliari
clean:
	rm -f $(BUILD_DIR)/$(MAIN).aux $(BUILD_DIR)/$(MAIN).log $(BUILD_DIR)/$(MAIN).out
	rm -f $(BUILD_DIR)/$(MAIN).toc $(BUILD_DIR)/$(MAIN).lof $(BUILD_DIR)/$(MAIN).lot
	rm -f $(BUILD_DIR)/$(MAIN).bbl $(BUILD_DIR)/$(MAIN).blg

# Pulizia completa (incluso PDF)
distclean: clean
	rm -f $(BUILD_DIR)/$(MAIN).pdf

.PHONY: all quick clean distclean
