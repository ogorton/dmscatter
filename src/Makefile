COMP = gfortran

F77_CFLAGS = -c -O2
F90_CFLAGS = -c -O2
F90_LFLAGS=

openmp_cflags = "-O2 -fopenmp -c"
openmp_lflags = "-fopenmp"
debug_cflags = "-O2 -fbounds-check -Wall -c"

exec=dmf90factor.x

OBJ1 = \
fileobject.o\
intlib.o\
sj2i.o\
progressmessage.o\
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
deventrate.o\
spectra.o\
totaleventrate.o\
main.o\
cg.o\
nucinputs.o\
dmf90factor.o

%.o : %.f90
	$(COMP) $(F90_CFLAGS) $(<F)

%.o : %.f
	$(COMP) $(F77_CFLAGS) $(<F)

default: $(OBJ1)
	$(COMP) -o $(exec) $(F90_LFLAGS) $(OBJ1) 

dmf90factor:
	$(MAKE)

serial: $(OBJ1) 
	$(MAKE)

openmp:
	$(MAKE) F90_CFLAGS=$(openmp_cflags) F90_LFLAGS=$(openmp_lflags)

debug:
	$(MAKE) F90_CFLAGS=$(debug_cflags) F77_CFLAGS=$(debug_cflags)

clean:
	rm *.o *.x *.mod
