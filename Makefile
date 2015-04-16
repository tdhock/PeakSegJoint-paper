HOCKING-PeakSegJoint-paper.pdf: HOCKING-PeakSegJoint-paper.tex
	rm -f *.aux *.bbl
	pdflatex HOCKING-PeakSegJoint-paper
	bibtex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
HOCKING-PeakSegJoint-slides.pdf: HOCKING-PeakSegJoint-slides.tex figure-profiles.tex table-H3K36me3.tex table-H3K4me3.tex
	rm -f *.aux *.bbl
	pdflatex HOCKING-PeakSegJoint-slides
figure-profiles.tex: figure-profiles.R
	R --no-save < $<
table-H3K36me3.tex: table-H3K36me3.R
	R --no-save < $<
table-H3K4me3.tex: table-H3K4me3.R
	R --no-save < $<