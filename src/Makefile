.DEFAULT_GOAL := default

COMP = gfortran
exec=dmfortfactor
bindir = ~/bin/

F77_CFLAGS = -c -Ofast
F90_CFLAGS = -c -Ofast
F90_LFLAGS=

openmp_cflags = "-Ofast -fopenmp -c"
openmp_lflags = "-fopenmp"
debug_cflags = "-Ofast -fbounds-check -Wall -c"

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

phi.o : bessel.o norm.o sigma.o

sigma.o : bessel.o norm.o

operme.o : sigma.o phi.o mjl.o

densities.o : orbitals.o

dmresponse.o : dmresponselib.o

nucresponse.o : orbitals.o mjl.o phi.o sigma.o densities.o

nucresponse_spectra.o : orbitals.o spectra.o

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
