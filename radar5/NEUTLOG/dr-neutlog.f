C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

        REAL*4 TARRAY(2)
        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=3
        INTEGER, PARAMETER :: NRDENS=2
        INTEGER, PARAMETER :: NGRID=0
        INTEGER, PARAMETER :: NLAGS=1
        INTEGER, PARAMETER :: NJACL=4
        INTEGER, PARAMETER :: MXST=4000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+2) :: IPAST
        REAL(kind=DP), dimension(3) :: RPAR
        REAL(kind=DP), dimension(1) :: Y1
        REAL(kind=DP), dimension(1) :: Y2
        REAL(kind=DP), dimension(1) :: Y3
        REAL(kind=DP), dimension(3) :: RTOL
        REAL(kind=DP), dimension(3) :: ATOL
        INTEGER, dimension(23) :: ISTAT
        EXTERNAL  FCN,PHI,ARGLAG,JFCN,JACLAG,QFUN,SOLOUT
  
C ------ FILE TO OPEN ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10

        TOL=1.D-9

C       PARAMETERS IN THE DIFFERENTIAL EQUATION
C ---   rho
        RPAR(1)=2.9D0
C ---   tau
        RPAR(2)=21.D0/50.D0
C ---   alpha
        RPAR(3)=0.1D0
        
C -----------------------------------------------------------------------
C
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
        IMAS=2
	  MLMAS=N
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES 
        X=0.0D0
        Y(1)=33.D-2
        Y(2)=222.D-2
        Y(3)=-0.1D0
C       Consistent with initial function
C --- DELAY
        TAU=RPAR(2)
C --- ENDPOINT OF INTEGRATION
        XEND=30.0D0
C       XEND=6.D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=1
C       TOL=1.D-6
        TFAC=1.D-3
C       TFAC=1.D0
        RTOL(1)=TOL
        RTOL(2)=TOL
        RTOL(3)=TOL   
        ATOL(1)=RTOL(1)
        ATOL(2)=RTOL(2)
        ATOL(3)=RTOL(3)*TFAC 
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.0D0  
C       IWORK(3)=50
C --- MAX NUMBER OF STEPS
        IWORK(2)=100000
C --- ERROR CONTROL
        IWORK(8)=1
        IWORK(11)=0
C     BREAKING POINTS SEARCH
        WORK(11)=2.D0
C --- WORKSPACE FOR PAST 
        IWORK(12)=MXST
C --- THE SECOND COMPONENT USES RETARDED ARGUMENT
        IWORK(15)=NRDENS
        IPAST(1)=1
        IPAST(2)=3
C --- SET THE PRESCRIBED GRID-POINTS
       DO I=1,NGRID
        GRID(I)=I*TAU
       END DO
C --- WORKSPACE FOR GRID 
       IWORK(13)=NGRID
C --- CONTROL OF NEWTON ITERATION
       IWORK(14)=1
C --- NEUTRAL PROBLEM
C --- DERIVATIVE COMPONENTS
       IWORK(16)=1
       IPAST(NRDENS+1)=1
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

C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,*) X,Y(1),Y(2),Y(3)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=IWORK(J)
         END DO
        WRITE(6,*)' ***** TOL=',TOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
        WRITE (6,92) ISTAT(23)

 90     FORMAT(1X,'X =',F8.2,'    Y =',3E20.12)
 91     FORMAT(' fcn=',I5,' jac=',I5,' step=',I5,
     &        ' accpt=',I5,' rejct=',I5,' dec=',I5,
     &        ' sol=',I5)
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
        REAL(kind=DP), PARAMETER :: XSTEP=0.001D0
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LRC) :: CONT
        REAL(kind=DP), dimension(1) :: RPAR
        EXTERNAL PHI
C       XOUT IS USED FOR THE DENSE OUTPUT
        COMMON /INTERN/XOUT

        WRITE (9,99) X,Y(1),Y(2),Y(3)

        IF (NR.EQ.1) THEN
           WRITE (10,99) X,Y(1),Y(2),Y(3)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
C ---   CONTINUOUS OUTPUT FOR RADAU5
              WRITE (10,99) XOUT,
     &         CONTR5(1,N,XOUT,CONT,X,HSOL),
     &         CONTR5(2,N,XOUT,CONT,X,HSOL),
     &         CONTR5(3,N,XOUT,CONT,X,HSOL)
              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF
  99     FORMAT(1X,'X =',F12.6,'    Y =',3F12.6)
C 99     FORMAT(1X,'X =',F12.6,'    Y =',3E23.15)

        RETURN
        END
C
        FUNCTION ARGLAG(IL,X,Y,RPAR,IPAR,
     &                  PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: Y
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        INTEGER, dimension(1) :: IPAR
        REAL(kind=DP), dimension(1) :: RPAR

        ARGLAG=X-RPAR(2)   

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
        
        rho=RPAR(1)
        tau=RPAR(2)
        alpha=RPAR(3)

        Y12=Y(1)*Y(1)
        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA,IPOS,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y1L1=YLAGR5(1,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS)
        Y3L1=YLAGR5(3,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS)

        F1=Y(1)*(1.D0-Y1L1-rho*Y3L1) - Y(2)*Y12/(1.D0+Y12)
        F(1)=F1  
        F(2)=Y(2)*(Y12/(1.D0+Y12)-alpha)
        F(3)=F1 - Y(3)

        RETURN
        END
C
        SUBROUTINE JFCN(N,X,Y,DFY,LDFY,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
C ----- STANDARD JACOBIAN OF THE EQUATION
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LDFY,N) :: DFY
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(1) :: RPAR
        EXTERNAL PHI

        rho=RPAR(1)
        tau=RPAR(2)
        alpha=RPAR(3)

        Y12=Y(1)*Y(1)
        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA,IPOS,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y1L1=YLAGR5(1,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS)
        Y3L1=YLAGR5(3,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS)

        DFY(1,1)=1.D0-Y1L1+(2.D0*Y12*Y(1)*Y(2))/((1.D0+Y12)*(1.D0+Y12))
     1          -(2.D0*Y(1)*Y(2))/(1.D0+Y12)-rho*Y3L1
        DFY(1,2)=-Y12/(1.D0+Y12)
        DFY(1,3)=0.D0
        DFY(2,1)=2.D0*Y(1)*Y(2)/((1.D0+Y12)**2) 
        DFY(2,2)=Y12/(1.D0+Y12)-alpha 
        DFY(2,3)=0.D0
        DFY(3,1)=DFY(1,1)
        DFY(3,2)=DFY(1,2) 
        DFY(3,3)=-1.D0 

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
                 
        rho=RPAR(1)
        tau=RPAR(2)
        alpha=RPAR(3)

        IVL(1)=1
        IVE(1)=1
        IVC(1)=1
        IVL(2)=1
        IVE(2)=1
        IVC(2)=3
        DFYL(1)=-Y(1)
        DFYL(2)=-rho*Y(1)
        IVL(3)=1
        IVE(3)=3
        IVC(3)=1
        IVL(4)=1
        IVE(4)=3
        IVC(4)=3
        DFYL(3)=DFYL(2)
        DFYL(4)=DFYL(1)

        RETURN
        END

C
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR
        SELECT CASE (I)
        CASE (1) 
            PHI= 33.D-2-X/10.D0
        CASE (2) 
            PHI= 222.D-2+X/10.D0
        CASE (3)
            PHI=-0.1D0
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
        Q(1,3)=0.D0
        Q(2,1)=0.D0
        Q(2,2)=1.D0
        Q(2,3)=0.D0
        Q(3,1)=0.D0
        Q(3,2)=0.D0
        Q(3,3)=0.D0
        RETURN
        END
