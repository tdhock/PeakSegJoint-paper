HOCKING-PeakSegJoint-paper.pdf: HOCKING-PeakSegJoint-paper.tex figure-PeakSegJoint.png refs.bib figure-test-error-dots.pdf figure-timings.tex figure-heuristic-algo.pdf figure-good-bad.pdf
	rm -f *.aux *.bbl
	pdflatex HOCKING-PeakSegJoint-paper
	bibtex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
	pdflatex HOCKING-PeakSegJoint-paper
HOCKING-PeakSegJoint-slides.pdf: HOCKING-PeakSegJoint-slides.tex figure-profiles.tex table-H3K36me3.tex table-H3K4me3.tex figure-bin-factor.pdf table-nrsf.tex table-H3K27ac.tex figure-heuristic-loss.pdf figure-weighted-error.pdf figure-label-problem-size.pdf figure-lasso-path.pdf figure-scatter-cheating-step1.pdf
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
chunk.problems.RData: chunk.problems.R
	R --no-save < $<
figure-weighted-error.pdf: figure-weighted-error.R chunk.problems.RData
	R --no-save < $<
histone.sets.RData: histone.sets.R
	R --no-save < $<
train.sets.RData: train.sets.R histone.sets.RData chunk.problems.RData
	R --no-save < $<
selected.by.set.RData: selected.by.set.R 
	R --no-save < $<
figure-label-problem-size.pdf: figure-label-problem-size.R selected.by.set.RData
	R --no-save < $<
figure-lasso-path.pdf: figure-lasso-path.R
	R --no-save < $<
step1.error.RData: step1.error.R selected.by.set.RData
	R --no-save < $<
figure-test-error-dots.pdf: figure-test-error-dots.R step2.error.RData 
	R --no-save < $<
cheating.error.RData: cheating.error.R selected.by.set.RData
	R --no-save < $<
figure-scatter-cheating-step1.pdf: figure-scatter-cheating-step1.R cheating.error.RData selected.by.set.RData
	R --no-save < $<
figure-timings.tex: figure-timings.R timings.RData
	R --no-save < $<
timings.RData: timings.R
	R --no-save < $<
figure-heuristic-algo.pdf: figure-heuristic-algo.R
	R --no-save < $<
figure-good-bad.pdf: figure-good-bad.R
	R --no-save < $<
