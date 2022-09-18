# Requires 
#  - python3 with matplotlib, numpy, scipy
#  - dmscatter
#  - pandoc
#

test:
	tar -xf data.tar.gz
	python3 validation.py
	python3 plots.py
	pandoc validation.md -o validation.pdf

clean:
	rm -r *.dat log.* perr.* validation.md data
