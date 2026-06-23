.PHONY: all clean data task1 task2 task3 report

all: data task1 task2 task3 report

data:
	python3 scripts/01_fetch_sdss.py

task1:
	python3 scripts/02_bpt_whan.py

task2:
	python3 scripts/03_fetch_wise.py
	python3 scripts/04_wise_diagram.py

task3:
	python3 scripts/05_filter_sii.py
	python3 scripts/06_density_mass.py

report:
	cd report && pdflatex report.tex && bibtex report && pdflatex report.tex && pdflatex report.tex

clean:
	rm -f data/*.fits output/*.pdf tables/*.tex
	rm -f report/report.pdf report/report.aux report/report.log report/report.bbl report/report.blg report/report.out
