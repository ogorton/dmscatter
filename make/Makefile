#############################
# Makefile for DMFortFactor #
#############################
.DEFAULT_GOAL := all

SRC_DIR = ../src
BIN_DIR = ../bin

EXE = $(BIN_DIR)/dmfortfactor
SRC = $(wildcard $(SRC_DIR)/*.f*)
OBJ = $(SRC:$(SRC_DIR)/%.f90=%.o)

# To use a different Fortran compiler, edit the following line:
FC = gfortran

F90_CFLAGS =-Ofast -fbounds-check
F90_LFLAGS =

OPENMP_CFLAGS = "-Ofast -fopenmp"
OPENMP_LFLAGS = "-fopenmp"

#####################################
# STANDARD OBJECT-SOURCE DEPENDENCY #
#####################################
$(EXE): $(OBJ) | $(BIN_DIR)
	$(FC) $(F90_CFLAGS) $^ -o $@

%.o : $(SRC_DIR)/%.f90 
	$(FC) $(F90_CFLAGS) -c $< -o $@

$(BIN_DIR):
	mkdir -p $@

##################################
# SUBROUTINE/MODULE DEPENDENCIES #
##################################
bessel.o:           wigner.o
coefficients.o:     main.o constants.o 	
constants.o:        
crosssection.o:     main.o constants.o transition.o
densities.o:        main.o orbitals.o
distributions.o:    constants.o 	
dmfortfactor.o:     spectra.o wigner.o nucresponse.o keywords.o totaleventrate.o parameters.o
dmresponse.o:       dmresponselib.o densities.o
dmresponselib.o:    constants.o
eventrate.o:        main.o quadrature.o constants.o crosssection.o distributions.o
keywords.o:         main.o constants.o settings.o quadrature.o wigner.o orbitals.o coefficients.o
main.o:             types.o
mjl.o:              norm.o constants.o wigner.o bessel.o
norm.o:     
nucinputs.o:        main.o settings.o orbitals.o densities.o
nucresponse.o:      main.o orbitals.o wigner.o densities.o mjl.o norm.o sigma.o phi.o
orbitals.o:         
parameters.o:       main.o constants.o quadrature.o orbitals.o settings.o keywords.o
phi.o:              norm.o bessel.o wigner.o sigma.o
quadrature.o:       
settings.o:         
sigma.o:            norm.o bessel.o wigner.o
spectra.o:          main.o constants.o crosssection.o transition.o eventrate.o quadrature.o settings.o orbitals.o nucresponse.o densities.o
totaleventrate.o:   main.o constants.o eventrate.o settings.o quadrature.o
transition.o:       main.o constants.o orbitals.o dmresponse.o nucresponse.o
types.o:                
wigner.o:           
zhangjin.o:   

################
# USER RECIPES #
################
all: $(EXE)

dmfortfactor:
	$(MAKE)

openmp:
	$(MAKE) F90_CFLAGS=$(OPENMP_CFLAGS) F90_LFLAGS=$(OPENMP_FLAGS)

clean:
	touch make.o make.mod
	rm *.o *.mod

cleaner:
	touch make.o make.mod $(BIN_DIR)
	rm *.o *.mod
	rm -r $(BIN_DIR)

help:
	@echo "Makefile for DMFortFactor"
	@echo "---"
	@echo "Here are the compiler options:"
	@echo "    make all -- (default)"
	@echo "    make dmfortfactor -- same as make all"
	@echo "    make openmp -- compile with OpenMP parallelization"
	@echo "    make clean -- remove object (.o) and module (.mod) files"
	@echo "    make cleaner -- remove object (.o) files, module (.mod) files, and bin directory"
