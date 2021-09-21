.DEFAULT_GOAL := default

COMP = gfortran
exec=dmfortfactor.x
bindir = ~/bin/

F77_CFLAGS = -c -O2
F90_CFLAGS = -c -O2
F90_LFLAGS=

openmp_cflags = "-O2 -fopenmp -c"
openmp_lflags = "-fopenmp"
debug_cflags = "-O2 -fbounds-check -Wall -c"

src = $(wildcard *.f90)
obj = $(src:.f90=.o)

#####################################
# STANDARD OBJECT-SOURCE DEPENDENCY #
#####################################
%.o : %.f90
	$(COMP) $(F90_CFLAGS) $(<F)

##################################
# SUBROUTINE/MODULE DEPENDENCIES #
##################################
bessel.o : wignerfunctions.o norm.o

mjl.o : bessel.o norm.o

phi.o : bessel.o norm.o

sigma.o : bessel.o norm.o

operme.o : sigma.o phi.o mjl.o

densities.o : orbitals.o

dmresponse.o : dmresponselib.o

nucresponse.o : orbitals.o

controlfiler.o : orbitals.o

transition.o : dmresponse.o nucresponse.o

crosssection.o : transition.o 

eventrate.o : gausquad.o gaus.o crosssection.o velocitydist.o 

spectra.o : eventrate.o

totaleventrate.o : eventrate.o gaus.o

dmfortfactor.o : spectra.o wignerfunctions.o modules.o nucresponse.o

################
# USER RECIPES #
################
default: modules.o $(obj)
	$(COMP) -o $(exec) $(F90_LFLAGS) $^

install:
	$(MAKE)
	mv $(exec) $(bindir)

install-openmp:
	$(MAKE) F90_CFLAGS=$(openmp_cflags) F90_LFLAGS=$(openmp_lflags)
	mv $(exec) $(bindir)

make new:
	make clean 
	$(MAKE)

dmfortfactor:
	$(MAKE)

serial: 
	$(MAKE)

openmp:
	$(MAKE) F90_CFLAGS=$(openmp_cflags) F90_LFLAGS=$(openmp_lflags)

debug:
	$(MAKE) F90_CFLAGS=$(debug_cflags) F77_CFLAGS=$(debug_cflags)

clean:
	touch make.o make.x make.mod
	rm *.o *.x *.mod
