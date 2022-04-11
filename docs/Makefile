manual:
	pandoc $(DOCS) -o manual.pdf --citeproc

html: 
	pandoc $(DOCS) -o index.html --citeproc --mathjax --standalone --toc

DOCS =\
	  source/intro.md\
	  source/python.md\
	  source/usage.md\
	  source/nuclear.md\
	  source/details.md\
	  source/keywords.md\
	  source/references.md
