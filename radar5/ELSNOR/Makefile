FC = f90

FCFLAGS = -c

DEBUGFLAGS = -O

LIB = /lib /usr/lib

.f.o:
	$(FC) $(FCFLAGS) $(DEBUGFLAGS) $*.f

# OBJ1 =    ../radar5-ver2.o ../dc_decdel.o ../decsol.o ../contr5.o
OBJ1 =    ../radar5.o ../dc_decdel.o ../decsol.o ../contr5.o
OBJ2 =    dr-exen.o     

prog :  $(OBJ2) 
	$(FC) $(DEBUGFLAGS) -o driver $(OBJ1) $(OBJ2)
