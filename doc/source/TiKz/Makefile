OUT_DIR = build

TEX_CMD = latexmk \
	--pdf \
	--output-directory=$(OUT_DIR) \
	--halt-on-error \
	--interaction=nonstopmode

PDF_2_PNG = convert \
	-density 600 \
	-transparent none

../pyBindings/figures/sector.png: sector.tex | build
	$(TEX_CMD) sector.tex
	$(PDF_2_PNG) build/sector.pdf $@

build:
	@mkdir -p build

clean:
	rm -rf *.pdf build

