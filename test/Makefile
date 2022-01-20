# Requires 
#  - python3 with matplotlib, numpy, scipy
#  - dmfortfactor
#  - pandoc
#

pdf:
	python3 plots.py
	pandoc validation.md -o validation.pdf

unpack:
	tar -xf data.tar.gz

test:
	tar -xf data.tar.gz
	python3 validation.py
	python3 plots.py
	pandoc validation.md -o validation.pdf

debug:
	rm *.dat
	python3 validation.py

clean:
	rm -r *.dat log.* perr.* validation.md data
