C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=2
        INTEGER, PARAMETER :: NRDENS=2
        INTEGER, PARAMETER :: NGRID=0
        INTEGER, PARAMETER :: NLAGS=1
        INTEGER, PARAMETER :: NJACL=2
        INTEGER, PARAMETER :: MXST=10000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        REAL, dimension(2) :: TARRAY
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        REAL(kind=DP), dimension(1) :: RTOL
        REAL(kind=DP), dimension(1) :: ATOL
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+1) :: IPAST
        INTEGER, dimension(22) :: ISTAT
        EXTERNAL  FCN,PHI,ARGLAG,JACLAG,QFUN,SOLOUT
  
C ------ FILE TO OPEN ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10

        WRITE(*,*) 'C = '
	READ(*,*) PAR
        RPAR=PAR

C -----------------------------------------------------------------------
C
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE STANDARD JACOBIAN NUMERICALLY 
        IJAC=0
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
        IMAS=1
        MLMAS=N
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES 
        X=0.0D0
        Y(1)= SIN(X)
        Y(2)= COS(X)
C       Consistent with initial function
C --- ENDPOINT OF INTEGRATION
        XEND=3.1415926535897932385D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.D-8
        ATOL=RTOL
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.0D0  
C --- WORKSPACE FOR PAST 
        IWORK(12)=MXST
C --- BOTH COMPONENTS USE RETARDED ARGUMENT
        IWORK(15)=NRDENS
        IPAST(1)=1
        IPAST(2)=2
C --- CONTROL OF NEWTON ITERATION
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
     &                  GRID,IPAST,QFUN,MLMAS,MUMAS)
        CALL DTIME(TARRAY)

C ---
        WRITE (6,*) 'Time = ',TARRAY(1)
C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,90) X,Y(1),Y(2)
C --- PRINT STATISTICS
         DO J=13,20
            ISTAT(J)=IWORK(J)
         END DO
        WRITE(6,*)' ***** TOL=',RTOL,' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
        WRITE (6,92) ISTAT(13)
 90     FORMAT(1X,'X =',F8.2,'    Y =',2E18.10)
 91     FORMAT(' fcn=',I7,' jac=',I6,' step=',I6,
     &        ' accpt=',I6,' rejct=',I6,' dec=',I6,
     &        ' sol=',I7)
 92     FORMAT(' full Newt. its =',I7)
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

        WRITE (9,99) X,Y(1),Y(2)
C
        IF (NR.EQ.1) THEN
           WRITE (10,99) X,Y(1),Y(2)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (10,99) XOUT,CONTR5(1,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(2,N,XOUT,CONT,X,HSOL)
              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F12.8,'    Y =',2E18.10)
        RETURN
        END
C
        FUNCTION ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: Y
	  REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        INTEGER, dimension(1) :: IPAR
        REAL(kind=DP), dimension(1) :: RPAR

C       ARGLAG=MIN(X,X*Y(1)**2)
        ARGLAG=X*Y(1)**2
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

        P=RPAR(1)

        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA,IPOS,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y1L1=YLAGR5(1,THETA,IPOS,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)
        Y2L1=YLAGR5(2,THETA,IPOS,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        F(1)= Y(2) 
        F(2)=-Y(2)+COS(X)*(1.D0+Y1L1)+ P*Y(1)*Y2L1 + 
     &        (1.D0-P)*SIN(X)*COS(X*SIN(X)*SIN(X)) -
     &        SIN(X+X*SIN(X)*SIN(X))

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

        P=RPAR(1)
        
        IVL(1)=1
        IVE(1)=2
        IVC(1)=1
        IVL(2)=1
        IVE(2)=2
        IVC(2)=2
        DFYL(1)=COS(X)
        DFYL(2)=P*Y(1)

        RETURN
        END

C
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR
        SELECT CASE (I)
        CASE (1)
            PHI=0.0D0        
        CASE (2)
            PHI=1.0D0 
        END SELECT
        RETURN
        END

        SUBROUTINE QFUN(N,Q,LQ,RPAR,IPAR)
C --- MATRIX "M" FOR THE TEST PROBLEM
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR
        REAL(kind=DP), dimension(LQ,N) :: Q

        Q(1,1)=1.D0
        Q(1,2)=0.D0
        Q(2,1)=0.D0
        Q(2,2)=0.D0
        RETURN
        END
