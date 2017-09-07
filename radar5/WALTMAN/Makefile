FC = f90

FCFLAGS = -c

DEBUGFLAGS = -O

LIB = /lib /usr/lib

.f.o:
	$(FC) $(FCFLAGS) $(DEBUGFLAGS) $*.f

OBJ2 =    ../radar5.o ../dc_decdel.o ../decsol.o ../contr5.o 
OBJ1 = dr-waltman.o     

prog :  $(OBJ1) 
	$(FC) $(DEBUGFLAGS) -o driver $(OBJ1) $(OBJ2)
