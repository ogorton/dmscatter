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
	$(unpack)
	python3 validation.py
	$(pdf)

debug:
	rm *.dat
	python3 validation.py
	$(pdf)

clean:
	rm -r *.dat log.* perr.* validation.md data
