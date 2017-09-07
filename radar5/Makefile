FC = f90

FCFLAGS = -c

DEBUGFLAGS = -O

LIB = /lib

.f.o:
	$(FC) $(FCFLAGS) $(DEBUGFLAGS) $*.f

OBJ1 =    radar5.o dc_decdel.o decsol.o contr5.o 
OBJ2 =    dr-hepat.o     

prog :  $(OBJ1) $(OBJ2) 
	$(FC) $(DEBUGFLAGS) -o driver $(OBJ1) $(OBJ2)
