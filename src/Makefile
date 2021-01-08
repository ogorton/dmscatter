COMP =  gfortran -pg
OPTIONS = -O2 -Wall -fopenmp 

OBJ1 = \
norm.o\
modules.o\
setparameters.o\
setcoeffs.o\
controlfiler.o\
getsps.o\
densities.o\
readdensities.o\
bessel.o\
summand.o\
phi.o\
mjl.o\
sigma.o\
zhangjin.o\
operme.o\
dmresponsefun.o\
dmresponse_lib.o\
dmresponse_iso.o\
dmresponse_pn.o\
nucresponse.o\
test_nucresponse.o\
transition.o\
diffcrosssection.o\
velocitydist.o\
quadrature.o\
eventrate.o\
eventrate_spectra.o\
cg.o\
setup_nuclearinputs.o\
darkmattermain.o

%.o : %.f90
	$(COMP) -c $(OPTIONS) $(<F)
.f.o:;  ($(COMP) -c $(OPTIONS) $(<F)) 

darkmatter.x: $(OBJ1) 
	$(COMP) -o darkmatter.x  $(OPTIONS) $(OBJ1) 
	cp darkmatter.x ~/bin/

clean:
	rm *.o *.x *.mod
