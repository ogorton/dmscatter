COMP =  gfortran
OPTIONS = -O2    # -r8 

OBJ1 = \
darkmattermain.o\
ResponseNuclear.o\
cg.o

%.o : %.f90
	$(COMP) -c $(OPTIONS) $(<F)
.f.o:;  ($(COMP) -c $(OPTIONS) $(<F)) 

darkmatter.x: $(OBJ1) 
	$(COMP) -o darkmatter.x -lblas -llapack $(OPTIONS) $(OBJ1) 

clean:
	rm *.o *.x
