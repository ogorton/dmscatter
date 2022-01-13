.DEFAULT_GOAL := default

#COMP = <compiler>
COMP = gfortran
exec=dmfortfactor
bindir = ~/bin/

F77_CFLAGS = -c -Ofast
F90_CFLAGS = -c -Ofast
F90_LFLAGS=

openmp_cflags = "-Ofast -fopenmp -c"
openmp_lflags = "-fopenmp"
debug_cflags = "-Ofast -fbounds-check -Wall -c"

src = $(wildcard *.f*)
obj = $(src:.f90=.o)

#####################################
# STANDARD OBJECT-SOURCE DEPENDENCY #
#####################################
%.o : %.f90
	$(COMP) $(F90_CFLAGS) $(<F)

%.o : %.f
	$(COMP) $(F77_CFLAGS) $(<F)

##################################
# SUBROUTINE/MODULE DEPENDENCIES #
##################################

setparameters.o : orbitals.o
bessel.o : wigner.o norm.o
mjl.o : bessel.o norm.o
phi.o : bessel.o norm.o sigma.o
sigma.o : bessel.o norm.o
operme.o : sigma.o phi.o mjl.o
densities.o : orbitals.o
dmresponse.o : dmresponselib.o densities.o
nucresponse.o : orbitals.o mjl.o phi.o sigma.o densities.o
nucresponse_spectra.o : orbitals.o spectra.o
controlfile.o : modules.o orbitals.o quadrature.o
transition.o : dmresponse.o nucresponse.o
crosssection.o : transition.o 
eventrate.o : quadrature.o crosssection.o velocitydist.o 
spectra.o : eventrate.o
totaleventrate.o : eventrate.o quadrature.o
dmfortfactor.o : spectra.o wigner.o modules.o nucresponse.o controlfile.o 

################
# USER RECIPES #
################
default: $(obj)
	$(COMP) -o $(exec) $(F90_LFLAGS) $^

install:
	mv $(exec) $(bindir)

install-openmp:
	$(MAKE) F90_CFLAGS=$(openmp_cflags) F90_LFLAGS=$(openmp_lflags)
	make install

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

tidy:
	mv *.o obj
	mv *.mod mod

clean:
	touch make.o make.x make.mod
	rm *.o *.x *.mod dmfortfactor
