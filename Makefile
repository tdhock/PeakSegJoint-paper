HOCKING-PeakSegJoint-paper.pdf: HOCKING-PeakSegJoint-paper.tex figure-PeakSegJoint.png refs.bib 
	rm -f *.aux *.bbl
	pdflatex HOCKING-PeakSegJoint-paper
	bibtex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
HOCKING-PeakSegJoint-slides.pdf: HOCKING-PeakSegJoint-slides.tex figure-profiles.tex table-H3K36me3.tex table-H3K4me3.tex figure-bin-factor.pdf table-nrsf.tex table-H3K27ac.tex figure-heuristic-loss.pdf
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
table-nrsf.tex: table-nrsf.R 
	R --no-save < $<
H3K27ac_TDH_control.RData: H3K27ac_TDH_control.R H3K27ac_TDH_control.txt
	R --no-save < $<
table-H3K27ac.tex: table-H3K27ac.R H3K27ac_TDH_control.RData
	R --no-save < $<
figure-heuristic-loss.pdf: figure-heuristic-loss.R
	R --no-save < $<