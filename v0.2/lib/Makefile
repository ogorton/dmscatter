COMP =  gfortran
OPTIONS = -O2    # -r8 

OBJ1 = \
modules.o\
getsps.o\
densities.o\
read.o\
darkmattermain.o\
bessel.o\
summand.o\
phi.o\
mjl.o\
sigma.o\
zhangjin.o\
operme.o\
dmresponse.f90\
cg.o

%.o : %.f90
	$(COMP) -c $(OPTIONS) $(<F)
.f.o:;  ($(COMP) -c $(OPTIONS) $(<F)) 

darkmatter.x: $(OBJ1) 
	$(COMP) -o darkmatter.x -lblas -llapack $(OPTIONS) $(OBJ1) 

clean:
	rm *.o *.x
