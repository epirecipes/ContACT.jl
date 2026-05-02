.PHONY: build mdgen html pdf clean all

PANDOC_FLAGS = --standalone --toc --mathjax \
    --css=style.css \
    --metadata title="Formal Verification of ContACT.jl" \
    --metadata author="Simon Frost" \
    --highlight-style=tango

MD_FILES = \
    ContactCat.md Coarsening.md Symmetrisation.md \
    Composition.md Stratification.md

all: build mdgen html

# Build the Lean 4 project (verifies all proofs)
build:
	lake build

# Generate markdown from Lean source via mdgen
mdgen:
	rm -rf .mdgen_tmp
	mkdir -p .mdgen_tmp
	lake exe mdgen ContACTProofs .mdgen_tmp
	cat $(addprefix .mdgen_tmp/, $(MD_FILES)) > ContACTProofs.md
	rm -rf .mdgen_tmp
	@echo "Generated ContACTProofs.md"

# Generate HTML from markdown
html: ContACTProofs.md
	pandoc ContACTProofs.md -o ContACTProofs.html $(PANDOC_FLAGS)
	@echo "Generated ContACTProofs.html"

# Generate PDF from markdown
pdf: ContACTProofs.md
	pandoc ContACTProofs.md -o ContACTProofs.pdf \
	    --pdf-engine=lualatex \
	    --toc \
	    --metadata title="Formal Verification of ContACT.jl" \
	    --metadata author="Simon Frost" \
	    --highlight-style=tango
	@echo "Generated ContACTProofs.pdf"

clean:
	rm -f ContACTProofs.html ContACTProofs.pdf ContACTProofs.md
	rm -rf .mdgen_tmp
