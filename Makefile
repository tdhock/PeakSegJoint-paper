HOCKING-PeakSegJoint-slides.pdf: HOCKING-PeakSegJoint-slides.tex figure-profiles.tex table-timings.tex
	rm -f *.aux *.bbl
	pdflatex HOCKING-PeakSegJoint-slides
figure-profiles.tex: figure-profiles.R
	R --no-save < $<
table-timings.tex: table-timings.R
	R --no-save < $<