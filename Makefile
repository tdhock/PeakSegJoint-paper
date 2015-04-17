HOCKING-PeakSegJoint-paper.pdf: HOCKING-PeakSegJoint-paper.tex figure-PeakSegJoint.png refs.bib
	rm -f *.aux *.bbl
	pdflatex HOCKING-PeakSegJoint-paper
	bibtex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
HOCKING-PeakSegJoint-slides.pdf: HOCKING-PeakSegJoint-slides.tex figure-profiles.tex table-H3K36me3.tex table-H3K4me3.tex figure-bin-factor.pdf
	rm -f *.aux *.bbl
	pdflatex HOCKING-PeakSegJoint-slides
figure-profiles.tex: figure-profiles.R
	R --no-save < $<
table-H3K36me3.tex: table-H3K36me3.R
	R --no-save < $<
table-H3K4me3.tex: table-H3K4me3.R
	R --no-save < $<
figure-PeakSegJoint.png: figure-PeakSegJoint.R
	R --no-save < $<
figure-bin-factor.pdf: figure-bin-factor.R 
	R --no-save < $<