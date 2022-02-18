.DEFAULT_GOAL := all

SRC_DIR = ../src
BIN_DIR = ../bin

EXE = $(BIN_DIR)/dmfortfactor
SRC = $(wildcard $(SRC_DIR)/*.f*)
OBJ = $(SRC:$(SRC_DIR)/%.f90=%.o)

FC = gfortran

F77_CFLAGS = -Ofast 
F90_CFLAGS = -Ofast
F90_LFLAGS=

openmp_cflags = "-Ofast -fopenmp"
openmp_lflags = "-fopenmp"
debug_cflags = "-Ofast -fbounds-check -Wall"

#####################################
# STANDARD OBJECT-SOURCE DEPENDENCY #
#####################################
$(EXE): $(OBJ) | $(BIN_DIR)
	$(FC) $(F90_CFLAGS) $^ -o $@

%.o : $(SRC_DIR)/%.f90 
	$(FC) $(F90_CFLAGS) -c $< -o $@

%.o : $(SRC_DIR)/%.f 
	$(FC) $(F77_CFLAGS) -c $< -o $@

$(BIN_DIR):
	mkdir -p $@

##################################
# SUBROUTINE/MODULE DEPENDENCIES #
##################################

setparameters.o : orbitals.o keywords.o
bessel.o : wigner.o norm.o
mjl.o : bessel.o norm.o
phi.o : bessel.o norm.o sigma.o
sigma.o : bessel.o norm.o
operme.o : sigma.o phi.o mjl.o
densities.o : orbitals.o
dmresponse.o : dmresponselib.o densities.o
nucresponse.o : orbitals.o mjl.o phi.o sigma.o densities.o
controlfile.o : modules.o orbitals.o quadrature.o
transition.o : dmresponse.o nucresponse.o
crosssection.o : transition.o 
eventrate.o : quadrature.o crosssection.o velocitydist.o main.o
spectra.o : eventrate.o 
totaleventrate.o : eventrate.o quadrature.o settings.o
dmfortfactor.o : spectra.o wigner.o nucresponse.o keywords.o 
velocitydist.o: constants.o	

################
# USER RECIPES #
################

all: $(EXE)

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
	touch make.o make.mod obj
	rm *.o *.mod
	rm -r obj
	rm -r ../bin
