# Requires 
#  - python3 with matplotlib, numpy, scipy
#  - dmfortfactor
#  - pandoc
#

test:
	python3 validation.py
	python3 plots.py
	pandoc validation.md -o validation.pdf

clean:
	rm -r *.dat log.* perr.* validation.md data
