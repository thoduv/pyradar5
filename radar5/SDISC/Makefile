FC = f90

FCFLAGS = -c

DEBUGFLAGS = -g

LIB = /lib /usr/lib

.f.o:
	$(FC) $(FCFLAGS) $(DEBUGFLAGS) $*.f


OBJ1 =    ../radar5.o ../dc_decdel.o ../decsol.o ../contr5.o 
OBJ2 =    dr-sdisc.o         

prog :  $(OBJ2) 
	$(FC) $(DEBUGFLAGS) -o driver $(OBJ1) $(OBJ2)
