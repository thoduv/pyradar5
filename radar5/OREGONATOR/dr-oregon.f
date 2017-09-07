C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

        REAL*4 TARRAY(2)
        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=2
        INTEGER, PARAMETER :: NRDENS=1
        INTEGER, PARAMETER :: NGRID=1
        INTEGER, PARAMETER :: NLAGS=1
        INTEGER, PARAMETER :: NJACL=2
        INTEGER, PARAMETER :: MXST=4000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+1) :: IPAST
        REAL(kind=DP), dimension(10) :: RPAR
        INTEGER, dimension(23) :: ISTAT
        EXTERNAL  FCN,PHI,ARGLAG,JFCN,JACLAG,SOLOUT
  
C ------ FILE TO OPEN ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10

C       PARAMETERS IN THE DIFFERENTIAL EQUATION
C ---   kM1
        RPAR(1)=1.34  
C ---   kM2
        RPAR(2)=1.6D9
C ---   kM3
        RPAR(3)=8.0D3
C ---   kM4
        RPAR(4)=4.0D7
C ---   kM5
        RPAR(5)=1.D0
C ---   fr  
        RPAR(6)=1.D0
C ---   A  
        RPAR(7)=6.D-2
C ---   B  
        RPAR(8)=6.D-2
C ---   Tau
        RPAR(9)=15.D-2

C -----------------------------------------------------------------------
C
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES 
        X=0.0D0
        Y(1)= 1.D-10
        Y(2)= 1.D-5
C       Consistent with initial function
C --- DELAY
        TAU=RPAR(9)
C --- ENDPOINT OF INTEGRATION
        XEND=100.5D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.D-9
        ATOL=RTOL*1.D-9
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.0D0  
C --- MAX NUMBER OF STEPS
        IWORK(2)=100000
C --- WORKSPACE FOR PAST 
        IWORK(12)=MXST
C --- THE SECOND COMPONENT USES RETARDED ARGUMENT
        IWORK(15)=NRDENS
        IPAST(1)=2
C --- SET THE PRESCRIBED GRID-POINTS
       DO I=1,NGRID
        GRID(I)=I*TAU
       END DO
C --- WORKSPACE FOR GRID 
       IWORK(13)=NGRID
C --- CONTROL OF NEWTON ITERATION
       IWORK(14)=1

C _____________________________________________________________________
C --- CALL OF THE SUBROUTINE RADAR5   
        CALL DTIME(TARRAY)
        CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JFCN,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  IMAS,SOLOUT,IOUT,
     &                  WORK,IWORK,RPAR,IPAR,IDID,
     &                  GRID,IPAST,DUMMY,MLMAS,MUMAS)
        CALL DTIME(TARRAY)

C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,90) X,Y(1),Y(2)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=IWORK(J)
         END DO
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
        WRITE (6,92) ISTAT(23)
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

        ARGLAG=X-RPAR(9)   

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
        
        kM1=RPAR(1)
        kM2=RPAR(2)
        kM3=RPAR(3)
        kM4=RPAR(4)
        kM5=RPAR(5)
        fr =RPAR(6)
        A  =RPAR(7)
        B  =RPAR(8)

        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y2L1=YLAGR5(2,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        F(1)=  kM1*A*Y(2) - kM2*Y(1)*Y2L1 + kM3*B*Y(1)-2.D0*kM4*Y(1)**2  
        F(2)= -kM1*A*Y(2) - kM2*Y(1)*Y2L1 + fr*kM3*B*Y(1)
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

        kM1=RPAR(1)
        kM2=RPAR(2)
        kM3=RPAR(3)
        kM4=RPAR(4)
        kM5=RPAR(5)
        fr =RPAR(6)
        A  =RPAR(7)
        B  =RPAR(8)

        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
        Y2L1=YLAGR5(2,THETA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

C       STANDARD JACOBIAN MATRIX 2x2
        DFY(1,1)= -kM2*Y2L1 + kM3*B - 4.D0*kM4*Y(1)
        DFY(1,2)=  kM1*A 
        DFY(2,1)= -kM2*Y2L1 + fr*kM3*B
        DFY(2,2)= -kM1*A
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
                 
        kM1=RPAR(1)
        kM2=RPAR(2)
        A  =RPAR(7)
        B  =RPAR(8)

C ---   ``B CASE''
        IVL(1)=1
        IVE(1)=1
        IVC(1)=2
        IVL(2)=1
        IVE(2)=2
        IVC(2)=2
        DFYL(1)=-kM2*Y(1)
        DFYL(2)=-kM2*Y(1)

        RETURN
        END

C
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR
        SELECT CASE (I)
C       CASE (1) 
C           PHI= 1.D-10       
C           This component is not required
        CASE (2) 
            PHI= 1.D-5
        END SELECT
        RETURN
        END
