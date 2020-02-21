COMP =  gfortran
OPTIONS = -O2 -fbounds-check   # -r8 

OBJ1 = \
modules.o\
setcoeffs.o\
getsps.o\
densities.o\
read.o\
bessel.o\
summand.o\
phi.o\
mjl.o\
sigma.o\
zhangjin.o\
operme.o\
dmresponse.o\
nucresponse.o\
test_nucresponse.o\
transition.o\
cg.o\
darkmattermain.o

%.o : %.f90
	$(COMP) -c $(OPTIONS) $(<F)
.f.o:;  ($(COMP) -c $(OPTIONS) $(<F)) 

darkmatter.x: $(OBJ1) 
	$(COMP) -o darkmatter.x  $(OPTIONS) $(OBJ1) 

clean:
	rm *.o *.x *.mod
