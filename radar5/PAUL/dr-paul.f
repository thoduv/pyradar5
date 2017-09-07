C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=1
        INTEGER, PARAMETER :: NRDENS=1
        INTEGER, PARAMETER :: NGRID=0
        INTEGER, PARAMETER :: NLAGS=1
        INTEGER, PARAMETER :: NJACL=1
        INTEGER, PARAMETER :: MXST=5000
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

        WRITE(*,*) 'INSERT TOL '
        READ(*,*) TOL
C -----------------------------------------------------------------------
C
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE STANDARD JACOBIAN ANALYTICALLY 
        IJAC=0
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
        MLMAS=N
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES 
        X=2.0D0
        Y(1)=1.0D0
C       Consistent with initial function
C --- ENDPOINT OF INTEGRATION
        XEND=5.5D0
C       XEND=2.D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=TOL   
        ATOL=RTOL
C --- INITIAL STEP SIZE
        H=1.0D-2
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.0D0
C --- ERROR CONTROL
        IWORK(11)=1
C --- STEPSIZE RATIOS
C       WORK(8)=0.5D0
C       WORK(9)=2.0D0
C --- WORKSPACE FOR PAST 
        IWORK(12)=MXST
C --- BOTH COMPONENTS USE RETARDED ARGUMENT
        IWORK(15)=NRDENS
        IPAST(1)=1
C --- CONTROL OF NEWTON ITERATION
        IWORK(3)=10
        IWORK(14)=1
C --- GRID
        IWORK(13)=NGRID

C _____________________________________________________________________
C --- CALL OF THE SUBROUTINE RADAR5   
        CALL DTIME(TARRAY)
        CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JFCN,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  IMAS,SOLOUT,IOUT,
     &                  WORK,IWORK,RPAR,IPAR,IDID,
     &                  GRID,IPAST,QFUN,MLMAS,MUMAS)
        CALL DTIME(TARRAY)

C ---
        WRITE (6,*) 'Time = ',TARRAY(1)
        YT=-2.D0*LOG(EXP(-2.D0)*(1.D0-XEND+LOG(4.D0*EXP(4.D0))))
C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,90) X,Y(1),ABS(Y(1)-YT)
C --- PRINT STATISTICS
         DO J=13,20
            ISTAT(J)=IWORK(J)
         END DO
        WRITE(6,*)' ***** TOL=',RTOL,' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
        WRITE (6,92) ISTAT(13)
 90     FORMAT(1X,'X =',F8.2,'    Y =',2E23.15)
 91     FORMAT(' fcn=',I7,' jac=',I6,' step=',I6,
     &        ' accpt=',I6,' rejct=',I6,' dec=',I6,
     &        ' sol=',I7)
 92     FORMAT(' full Newt. its =',I7)
        WRITE(6,*) 'SOLUTION IS TABULATED IN FILES: sol.out & cont.out'

        CLOSE(9)
        CLOSE(10)

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
        REAL(kind=DP), PARAMETER :: XSTEP=0.002D0
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LRC) :: CONT
        REAL(kind=DP), dimension(1) :: RPAR
C       EXTERNAL PHI
C       XOUT IS USED FOR THE DENSE OUTPUT
        COMMON /INTERN/XOUT

        WRITE (9,99) X,Y(1)

        IF (NR.EQ.1) THEN
           WRITE (10,99) X,Y(1)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (10,99) XOUT,CONTR5(1,N,XOUT,CONT,X,HSOL)
              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F12.8,'    Y =',2E18.10)
        RETURN
        END
C
        DOUBLE PRECISION FUNCTION ARGLAG(IL,X,Y,RPAR,IPAR,
     &                                   PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: Y
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        INTEGER, dimension(1) :: IPAR
        REAL(kind=DP), dimension(1) :: RPAR

        ARGLAG=Y(1)

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
        REAL RN
        EXTERNAL PHI

        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA,IPOS,RPAR,IPAR,
     &                 PHI,IPAST,NRDS)
        Y1L1=YLAGR5(1,THETA,IPOS,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        F(1)= Y1L1
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

        IVL(1)=1
        IVE(1)=1
        IVC(1)=1
        DFYL(1)=1.D0 

        RETURN
        END

C
CD      SUBROUTINE PHI(I,X,VAL,RPAR,IPAR)
        DOUBLE PRECISION FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR 

        PHI=0.5D0

        RETURN 
        END 

        SUBROUTINE QFUN(N,Q,LQ,RPAR,IPAR)
C --- MATRIX "M" FOR THE TEST PROBLEM
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR
        REAL(kind=DP), dimension(LQ,N) :: Q

        Q(1,1)=1.D0
        RETURN
        END

