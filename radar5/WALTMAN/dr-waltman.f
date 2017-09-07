C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

C --- PARAMETERS FOR RADAR5 (FULL JACOBIAN)
        REAL*4 TARRAY(2)
        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=6
        INTEGER, PARAMETER :: NRDENS=3
        INTEGER, PARAMETER :: NGRID=2
        INTEGER, PARAMETER :: NLAGS=2
C
        INTEGER, PARAMETER :: NJACL=9
        INTEGER, PARAMETER :: MXST=50000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        INTEGER, PARAMETER :: NPAR=15
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        REAL(kind=DP), dimension(ND) :: RTOL
        REAL(kind=DP), dimension(ND) :: ATOL
        REAL(kind=DP), dimension(ND) :: YT
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+1) :: IPAST
        REAL(kind=DP), dimension(NPAR) :: RPAR
        INTEGER, dimension(22) :: ISTAT
        EXTERNAL  FCN,PHI,ARGLAG,JACLAG,SOLOUT

C ------ FILE DE DONNEES ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10

  
C --- PARAMETER IN THE DIFFERENTIAL EQUATION
        
	WRITE(6,*) 'TOLL = '
        READ(*,*) TOLL
C ---   alpha
        RPAR(1)= 1.8D0   
C ---   beta
        RPAR(2)= 20.D0
C ---   gamma
        RPAR(3)= 0.002D0
C ---   r
        RPAR(4)= 5.0D4
C ---   s
        RPAR(5)= 1.0D5
C ---   m1
        RPAR(6)= 0.3D-13
C ---   m2
        RPAR(7)= 1.1D-7
C ---   x0
        RPAR(8)= 0.5D-5
C ---   y0
        RPAR(9)= 1.0D-15
C ---   w0,z0
        RPAR(10)=0.D0
C ---   t0
        RPAR(11)=35.D0
C ---   t1
        RPAR(12)=197.D0

C -----------------------------------------------------------------------
C
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE STANDARD JACOBIAN NUMERICALLY  
        IJAC=0
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM 
        IMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES 
        X=0.D0     
        Y(1)=RPAR(8) 
        Y(2)=RPAR(9) 
        Y(3)=RPAR(10) 
        Y(4)=RPAR(10) 
        Y(5)=0.0D0
        Y(6)=0.0D0
C       Consistent with initial function
C --- ENDPOINT OF INTEGRATION
        XEND=300.0D0
        XEND=300.0D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE: DIFFER FOR DIFFERENT COMP.
        ITOL=1
        DO I=1,4
         RTOL(I)=TOLL  
         ATOL(I)=RTOL(I)*1.0D-12
        END DO
        RTOL(5)=TOLL
        RTOL(6)=TOLL
        ATOL(5)=RTOL(5)
        ATOL(6)=RTOL(6)
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.0D0  
C --- PARAMETERS FOR ERROR ESTIMATION
        WORK(10)=1.D0
        IWORK(11)=0
C --- PARAMETER FOR BREAKING POINT SEARCH
        WORK(11)=1.0D0
C --- WORKSPACE FOR PAST
        IWORK(12)=MXST
C --- COMPONENTS WHICH USE RETARDED ARGUMENTS
        IWORK(15)=NRDENS
        IPAST(1)=1
        IPAST(2)=2
        IPAST(3)=3
C ---  SET THE PRESCRIBED GRID-POINTS
C       If you don't want to set the grid-points
        IWORK(13)=NGRID
        IF (NGRID.EQ.1) THEN
         GRID(1)=RPAR(11)
        ELSE IF (NGRID.EQ.2) THEN
         GRID(1)=RPAR(11)
         GRID(2)=RPAR(12)
        END IF
C ---  CONTROL OF NEWTON ITERATION (FULL ENABLED)
        IWORK(14)=2

C _____________________________________________________________________
C --- CALL OF THE SUBROUTINE RADAR5   
        CALL DTIME(TARRAY)
        CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  DUMMY,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  IMAS,SOLOUT,IOUT,
     &                  WORK,IWORK,RPAR,IPAR,IDID,
     &                  GRID,IPAST,DUMMY,MLMAS,MUMAS)
        CALL DTIME(TARRAY)

C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,90) X,Y(1),Y(2),Y(3),Y(4)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=IWORK(J)
         END DO
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
        WRITE (6,92) ISTAT(13)
 90     FORMAT(1X,'X =',F8.2,'    Y =',4E18.10)
 91     FORMAT(' fcn=',I7,' jac=',I6,' step=',I6,
     &        ' accpt=',I6,' rejct=',I6,' dec=',I6,
     &        ' sol=',I7)
 92     FORMAT(' full its =',I7)
        WRITE(6,*) 'SOLUTION IS TABULATED IN FILES: sol.out & cont.out'

        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,HSOL,Y,CONT,LRC,N,
     &                     RPAR,IPAR,IRTRN)
C ----- PRINTS THE DISCRETE OUTPUT AND THE CONTINUOUS OUTPUT
C       AT EQUIDISTANT OUTPUT-POINTS
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), PARAMETER :: XSTEP=0.01D0
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LRC) :: CONT
        REAL(kind=DP), dimension(1) :: RPAR
        EXTERNAL PHI
C       XOUT IS USED FOR THE DENSE OUTPUT
        COMMON /INTERN/XOUT

        WRITE (9,99) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
C
        IF (NR.EQ.1) THEN
           WRITE (10,99) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
C ---   CONTINUOUS OUTPUT FOR RADAU5
              WRITE (10,99) XOUT,
     &         CONTR5(1,N,XOUT,CONT,X,HSOL),
     &         CONTR5(2,N,XOUT,CONT,X,HSOL),
     &         CONTR5(3,N,XOUT,CONT,X,HSOL), 
     &         CONTR5(4,N,XOUT,CONT,X,HSOL),
     &         CONTR5(5,N,XOUT,CONT,X,HSOL),
     &         CONTR5(6,N,XOUT,CONT,X,HSOL)
              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F12.6,'    Y =',6E18.10)
        RETURN
        END
C
        FUNCTION ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: Y
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        INTEGER, dimension(1) :: IPAR
        REAL(kind=DP), dimension(1) :: RPAR
        EXTERNAL PHI
        SELECT CASE (IL)
        CASE (1)
          ARGLAG=Y(5)*FHS(X-RPAR(11))
        CASE (2) 
          ARGLAG=Y(6)*FHS(X-RPAR(12))
        END SELECT
        RETURN
        END
C
        SUBROUTINE FCN(N,X,Y,F,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(N) :: F
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(1) :: RPAR
        EXTERNAL PHI

C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       SECOND DELAY
        CALL LAGR5(2,X,Y,ARGLAG,PAST,THETA2,IPOS2,RPAR,IPAR,
     &             PHI,IPAST,NRDS)

C	  WRITE(6,*) IPOS1,IPOS2,THETA1,THETA2

        Y1L1=YLAGR5(1,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y2L1=YLAGR5(2,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y3L1=YLAGR5(3,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y1L2=YLAGR5(1,THETA2,IPOS2,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y2L2=YLAGR5(2,THETA2,IPOS2,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y3L2=YLAGR5(3,THETA2,IPOS2,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        F1L=F1(Y1L1,Y2L1,Y3L1)
        F2L=F2(Y2L2,Y3L2)

        G1=FHS(X-RPAR(11))
        G2=FHS(X-RPAR(12))
        F(1)= -RPAR(4)*Y(1)*Y(2) - RPAR(5)*Y(1)*Y(4) 
        F(2)= -RPAR(4)*Y(1)*Y(2) + RPAR(1)*RPAR(4)*Y1L1*Y2L1*G1
        F(3)=  RPAR(4)*Y(1)*Y(2)
        F(4)= -RPAR(5)*Y(1)*Y(4) - RPAR(3)*Y(4) + RPAR(2)*RPAR(4)*
     &         Y1L2*Y2L2*G2
        IF (X.LT.RPAR(11)) THEN
	   F(5)= 0.D0
	  ELSE
	   F(5)=  F1(Y(1),Y(2),Y(3))/F1L*G1
        END IF
        IF (X.LT.RPAR(12)) THEN
	   F(6)= 0.D0
        ELSE
	   F(6)=  F2(Y(2),Y(3))/F2L*G2
        END IF

        RETURN
        END
C
        FUNCTION FHS(T)
        IMPLICIT REAL*8 (A-H,O-Z)

 1      IF (T.LT.0.D0) THEN
          FHS=0.D0
        ELSE
          FHS=1.D0
        END IF
        RETURN
        END
C
        FUNCTION F1(X,Y,W)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        F1=X*Y+W
        RETURN
        END
C
        FUNCTION F1X(X,Y,W)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        F1X=Y
        RETURN
        END
C
        FUNCTION F1Y(X,Y,W)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        F1Y=X
        RETURN
        END
C
        FUNCTION F1W(X,Y,W)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        F1W=1.D0
        RETURN
        END
C
        FUNCTION DFL1(I,X,Y,W,XL,YL,WL)
        IMPLICIT REAL*8 (A-H,O-Z)

        SELECT CASE (I)
        CASE (1)
          DFL1= -(F1(X,Y,W)*F1X(XL,YL,WL))/(F1(XL,YL,WL)**2)
        CASE (2)              
          DFL1= -(F1(X,Y,W)*F1Y(XL,YL,WL))/(F1(XL,YL,WL)**2)
        CASE (3)              
          DFL1= -(F1(X,Y,W)*F1W(XL,YL,WL))/(F1(XL,YL,WL)**2)
        END SELECT
        RETURN
        END
C
        FUNCTION F2(Y,W)
        IMPLICIT REAL*8 (A-H,O-Z)
        CONST=0.D0

        F2=CONST+Y+W
        RETURN
        END
C
        FUNCTION F2Y(Y,W)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        F2Y=1.D0
        RETURN
        END
C
        FUNCTION F2W(Y,W)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        F2W=1.D0
        RETURN
        END
C
        FUNCTION DFL2(I,Y,W,YL,WL)
        IMPLICIT REAL*8 (A-H,O-Z)

        SELECT CASE (I)
        CASE (1)
          DFL2= -(F2(Y,W)*F2Y(YL,WL))/(F2(YL,WL)**2)
        CASE (2)             
          DFL2= -(F2(Y,W)*F2W(YL,WL))/(F2(YL,WL)**2)
        END SELECT
        RETURN
        END
C
        SUBROUTINE JACLAG(N,X,Y,DFYL,ARGLAG,PHI,IVE,IVC,IVL,
     &                    RPAR,IPAR,PAST,IPAST,NRDS)
C ----- JACOBIAN OF DELAY TERMS IN THE EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(1) :: DFYL
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(1) :: RPAR
        INTEGER, dimension(1) :: IVE,IVC,IVL
        EXTERNAL PHI

C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       SECOND DELAY
        CALL LAGR5(2,X,Y,ARGLAG,PAST,THETA2,IPOS2,RPAR,IPAR,
     &             PHI,IPAST,NRDS)

        Y1L1=YLAGR5(1,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y2L1=YLAGR5(2,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y3L1=YLAGR5(3,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y1L2=YLAGR5(1,THETA2,IPOS2,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y2L2=YLAGR5(2,THETA2,IPOS2,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y3L2=YLAGR5(3,THETA2,IPOS2,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        G1=FHS(X-RPAR(11))
        G2=FHS(X-RPAR(12))
        IVL(1)=1
        IVE(1)=2
        IVC(1)=1
        DFYL(1)=RPAR(1)*RPAR(4)*Y2L1*G1
        IVL(2)=1
        IVE(2)=2
        IVC(2)=2
        DFYL(2)=RPAR(1)*RPAR(4)*Y1L1*G1

        IVL(3)=2
        IVE(3)=4
        IVC(3)=1
        DFYL(3)=RPAR(2)*RPAR(4)*Y2L2*G2
        IVL(4)=2
        IVE(4)=4
        IVC(4)=2
        DFYL(4)=RPAR(2)*RPAR(4)*Y1L2*G2

        IVL(5)=1
        IVE(5)=5
        IVC(5)=1
        DFYL(5)=DFL1(1,Y(1),Y(2),Y(3),Y1L1,Y2L1,Y3L1)*G1
        IVL(6)=1
        IVE(6)=5
        IVC(6)=2
        DFYL(6)=DFL1(2,Y(1),Y(2),Y(3),Y1L1,Y2L1,Y3L1)*G1
        IVL(7)=1
        IVE(7)=5
        IVC(7)=3
        DFYL(7)=DFL1(3,Y(1),Y(2),Y(3),Y1L1,Y2L1,Y3L1)*G1

        IVL(8)=2
        IVE(8)=6
        IVC(8)=2
        DFYL(8)=DFL2(1,Y(2),Y(3),Y2L2,Y3L2)*G2
        IVL(9)=2
        IVE(9)=6
        IVC(9)=3
        DFYL(9)=DFL2(2,Y(2),Y(3),Y2L2,Y3L2)*G2

	RETURN
        END

C
        DOUBLE PRECISION FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR
        SELECT CASE (I)
        CASE (1)
            PHI=RPAR(8)
        CASE (2)             
            PHI=RPAR(9)
        CASE (3)             
            PHI=RPAR(10)
        CASE (4)             
            PHI=RPAR(10)
        CASE (5)             
            PHI=0.D0
        CASE (6)             
            PHI=0.D0
        END SELECT
        RETURN
        END

C
