C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=10
        INTEGER, PARAMETER :: NRDENS=5
        INTEGER, PARAMETER :: NGRID=3
        INTEGER, PARAMETER :: NLAGS=5
        INTEGER, PARAMETER :: NJACL=13
        INTEGER, PARAMETER :: MXST=10000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+1) :: IPAST
        REAL(kind=DP), dimension(45) :: RPAR
        INTEGER, dimension(22) :: ISTAT
        EXTERNAL  FCN,PHI,ARGLAG,JFCN,JACLAG,SOLOUT

C ------ FILE TO OPEN ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10

C ---   PARAMETER IN THE DIFFERENTIAL EQUATION
C ---   alpha1
        RPAR(1)= 83.0D0
C ---   alpha2
        RPAR(2)= 5.0D0
C ---   alpha3
        RPAR(3)= 6.6D14
C ---   alpha4
        RPAR(4)= 3.0D11
C ---   alpha5
        RPAR(5)= 0.4D0

C ---   alpha6
        RPAR(6)= 2.5D7
C ---   alpha7
        RPAR(7)= 0.5D-12
C ---   alpha8
        RPAR(8)= 2.3D9
C ---   alpha9
        RPAR(9)= 0.052D0
C ---   alpha10
        RPAR(10)= 0.15D0

C ---   alpha11
        RPAR(11)= 9.4D9
C ---   alpha12
        RPAR(12)= 1.0D-15
C ---   alpha13
        RPAR(13)= 1.2D0
C ---   alpha14
        RPAR(14)= 2.7D16
C ---   alpha15
        RPAR(15)= 2.0D0

C ---   alpha16
        RPAR(16)= 5.3D27
C ---   alpha17
        RPAR(17)= 1.0D0
C ---   alpha18
        RPAR(18)= 1.0D-18
C ---   alpha19
        RPAR(19)= 2.7D16
C ---   alpha20
        RPAR(20)= 2.0D0

C ---   alpha21
        RPAR(21)= 8.0D28
C ---   alpha22
        RPAR(22)= 1.0D0
C ---   alpha23
        RPAR(23)= 1.0D-19
C ---   alpha24
        RPAR(24)= 5.3D33
C ---   alpha25
        RPAR(25)= 16.0D0

C ---   alpha26
        RPAR(26)= 1.6D14
C ---   alpha27
        RPAR(27)= 0.4D0 
C ---   alpha28
        RPAR(28)= 1.0D-18
C ---   alpha29
        RPAR(29)= 8.0D32
C ---   alpha30
        RPAR(30)= 16.0D0

C ---   alpha31
        RPAR(31)= 0.1D0
C ---   alpha32
        RPAR(32)= 1.0D-18
C ---   alpha33
        RPAR(33)= 1.7D30
C ---   alpha34
        RPAR(34)= 3.0D0
C ---   alpha35
        RPAR(35)= 0.4D0

C ---   alpha36
        RPAR(36)= 4.3D-22
C ---   alpha37
        RPAR(37)= 0.85D7
C ---   alpha38
        RPAR(38)= 8.6D11
C ---   alpha39
        RPAR(39)= 0.043D0

C  ---  Tau1
        RPAR(41)= 0.6D0
C  ---  Tau2
        RPAR(42)= 0.6D0
C  ---  Tau3
        RPAR(43)= 2.0D0
C  ---  Tau4
        RPAR(44)= 2.0D0
C  ---  Tau5
        RPAR(45)= 3.0D0

C -----------------------------------------------------------------------
C
C ---   DIMENSION OF THE SYSTEM
        N=ND
C ---   COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C ---   JACOBIAN IS A FULL MATRIX
        MLJAC=N
C ---   DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
        MLMAS=N

C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES 
        X=0.0D0
        Y(1)=2.9D-16
        Y(2)=0.0D0
        Y(3)=0.0D0
        Y(4)=0.0D0
        Y(5)=RPAR(18) 
        Y(6)=RPAR(23) 
        Y(7)=RPAR(28) 
        Y(8)=RPAR(32)
        Y(9)=RPAR(36)
        Y(10)=RPAR(37)*RPAR(36)/RPAR(39)
C       Consistent with initial function
C --- ENDPOINT OF INTEGRATION
        XEND=130.0D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.D-11
        ATOL=RTOL*1.0D-20
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
        ISTAT(I)=0
  10    WORK(I)=0.0D0  
C --- WORKSPACE FOR PAST 
        IWORK(12)=MXST
C --- COMPONENTS 4,5,6,7,8 USE RETARDED ARGUMENTS
        IWORK(15)=NRDENS
        IPAST(1)=4
        IPAST(2)=5
        IPAST(3)=6
        IPAST(4)=7
        IPAST(5)=8
C ----  
C --- SET THE PRESCRIBED GRID-POINTS
        GRID(1)=RPAR(41)
        GRID(2)=RPAR(43)
        GRID(3)=RPAR(45)      
C       WRITE(6,*) (GRID(I),I=1,IWORK(13))
C _____________________________________________________________________
C --- CALL OF THE SUBROUTINE RADAR5   
        CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JFCN,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  IMAS,SOLOUT,IOUT,
     &                  WORK,IWORK,RPAR,IPAR,IDID,
     &                  GRID,IPAST,DUMMY,MLMAS,MUMAS)

C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,90) X,Y(1),Y(2),Y(3),Y(4),Y(5)
        WRITE (6,90) X,Y(6),Y(7),Y(8),Y(9),Y(10)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=IWORK(J)
         END DO
        WRITE(6,*)' ***** TOL=',RTOL,' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
        WRITE (6,92) ISTAT(13)
 90     FORMAT(1X,'X =',F8.2,'    Y =',5E18.10)
 91     FORMAT(' fcn=',I7,' jac=',I6,' step=',I6,
     &        ' accpt=',I6,' rejct=',I6,' dec=',I6,
     &        ' sol=',I7)
 92     FORMAT(' full Newt. its = ',I7)
        WRITE(6,*) 'SOLUTION IS TABULATED IN FILES: sol.out & cont.out'

        STOP
        END

C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,HSOL,Y,CONT,LRC,N,
     &                     RPAR,IPAR,IRTRN)
C ----- PRINTS THE DISCRETE OUTPUT AND THE CONTINUOUS OUTPUT
C       OF THE FIRST TWO COMPONENTS AT EQUIDISTANT OUTPUT-POINTS
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), PARAMETER :: XSTEP=0.01D0
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LRC) :: CONT
        REAL(kind=DP), dimension(1) :: RPAR
        EXTERNAL PHI
C       XOUT IS USED FOR THE DENSE OUTPUT
        COMMON /INTERN/XOUT

C
        IF (NR.EQ.1) THEN
           WRITE (10,99) X,Y(1),Y(3)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (10,99) XOUT,CONTR5(1,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(3,N,XOUT,CONT,X,HSOL)
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

        SELECT CASE (IL) 
        CASE (1) 
                 ARGLAG=X-RPAR(41) 
        CASE (2) 
                 ARGLAG=X-RPAR(42) 
        CASE (3) 
                 ARGLAG=X-RPAR(43) 
        CASE (4) 
                 ARGLAG=X-RPAR(44) 
        CASE (5) 
                 ARGLAG=X-RPAR(45) 
        END SELECT
        RETURN
        END

        SUBROUTINE FCN(N,X,Y,F,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(N) :: F
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(1) :: RPAR
        EXTERNAL PHI,ARGLAG
        XI(Z)=1.0D0-Z/RPAR(7)

C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       SECOND DELAY
        CALL LAGR5(2,X,Y,ARGLAG,PAST,THETA2,IPOS2,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       THIRD DELAY
        CALL LAGR5(3,X,Y,ARGLAG,PAST,THETA3,IPOS3,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       FOURTH DELAY
        CALL LAGR5(4,X,Y,ARGLAG,PAST,THETA4,IPOS4,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       FIFTH DELAY
        CALL LAGR5(5,X,Y,ARGLAG,PAST,THETA5,IPOS5,RPAR,IPAR,
     &             PHI,IPAST,NRDS)

        Y4L1=YLAGR5(4,THETA1,IPOS1,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L2=YLAGR5(4,THETA2,IPOS2,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L3=YLAGR5(4,THETA3,IPOS3,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L4=YLAGR5(4,THETA4,IPOS4,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L5=YLAGR5(4,THETA5,IPOS5,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y5L1=YLAGR5(5,THETA1,IPOS1,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y5L3=YLAGR5(5,THETA3,IPOS3,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y6L2=YLAGR5(6,THETA2,IPOS2,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y6L4=YLAGR5(6,THETA4,IPOS4,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y6L5=YLAGR5(6,THETA5,IPOS5,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y7L3=YLAGR5(7,THETA3,IPOS3,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y8L4=YLAGR5(8,THETA4,IPOS4,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y8L5=YLAGR5(8,THETA5,IPOS5,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)

        F(1)= RPAR(1)*Y(2) + RPAR(2)*RPAR(3)*Y(2)*Y(7) - 
     &        RPAR(4)*Y(1)*Y(10) - RPAR(5)*Y(1) - 
     &        RPAR(6)*Y(1)*(RPAR(7)-Y(2)-Y(3)) 
        F(2)= RPAR(8)*Y(1)*(RPAR(7)-Y(2)-Y(3)) - RPAR(3)*Y(2)*Y(7) -
     &        RPAR(9)*Y(2)
        F(3)= RPAR(3)*Y(2)*Y(7) + RPAR(9)*Y(2) - RPAR(10)*Y(3)
        F(4)= RPAR(11)*RPAR(12)*Y(1) - RPAR(13)*Y(4)
        F(5)= RPAR(14)*(XI(Y(3))*RPAR(15)*Y4L1*Y5L1-Y(4)*Y(5)) -
     &        RPAR(16)*Y(4)*Y(5)*Y(7) + RPAR(17)*(RPAR(18)-Y(5))
        F(6)= RPAR(19)*(XI(Y(3))*RPAR(20)*Y4L2*Y6L2-Y(4)*Y(6)) -
     &        RPAR(21)*Y(4)*Y(6)*Y(8) + RPAR(22)*(RPAR(23)-Y(6))
        F(7)= RPAR(24)*(XI(Y(3))*RPAR(25)*Y4L3*Y5L3*Y7L3 - 
     &        Y(4)*Y(5)*Y(7)) - RPAR(26)*Y(2)*Y(7) + 
     &        RPAR(27)*(RPAR(28)-Y(7))
        F(8)= RPAR(29)*(XI(Y(3))*RPAR(30)*Y4L4*Y6L4*Y8L4 -
     &        Y(4)*Y(6)*Y(8)) + RPAR(31)*(RPAR(32)-Y(8))
        F(9)= RPAR(33)*XI(Y(3))*RPAR(34)*Y4L5*Y6L5*Y8L5 +
     &        RPAR(35)*(RPAR(36)-Y(9))
        F(10)= RPAR(37)*Y(9) - RPAR(38)*Y(10)*Y(1) - RPAR(39)*Y(10)
        
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
        EXTERNAL PHI,ARGLAG
        XI(Z)=1.0D0-Z/RPAR(7)

C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       SECOND DELAY
        CALL LAGR5(2,X,Y,ARGLAG,PAST,THETA2,IPOS2,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       THIRD DELAY
        CALL LAGR5(3,X,Y,ARGLAG,PAST,THETA3,IPOS3,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       FOURTH DELAY
        CALL LAGR5(4,X,Y,ARGLAG,PAST,THETA4,IPOS4,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       FIFTH DELAY
        CALL LAGR5(5,X,Y,ARGLAG,PAST,THETA5,IPOS5,RPAR,IPAR,
     &             PHI,IPAST,NRDS)

        Y4L1=YLAGR5(4,THETA1,IPOS1,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L2=YLAGR5(4,THETA2,IPOS2,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L3=YLAGR5(4,THETA3,IPOS3,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L4=YLAGR5(4,THETA4,IPOS4,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y4L5=YLAGR5(4,THETA5,IPOS5,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y5L1=YLAGR5(5,THETA1,IPOS1,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y5L3=YLAGR5(5,THETA3,IPOS3,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y6L2=YLAGR5(6,THETA2,IPOS2,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y6L4=YLAGR5(6,THETA4,IPOS4,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y6L5=YLAGR5(6,THETA5,IPOS5,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y7L3=YLAGR5(7,THETA3,IPOS3,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y8L4=YLAGR5(8,THETA4,IPOS4,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)
        Y8L5=YLAGR5(8,THETA5,IPOS5,PHI,RPAR,IPAR,
     &            PAST,IPAST,NRDS)

C       STANDARD JACOBIAN MATRIX 10x10
        DFY(1,1)=-RPAR(4)*Y(10) - RPAR(5) - RPAR(6)*(RPAR(7)-Y(2)-Y(3))
        DFY(1,2)= RPAR(1) + RPAR(2)*RPAR(3)*Y(7) + RPAR(6)*Y(1)
        DFY(1,3)= RPAR(6)*Y(1)
        DFY(1,4)=0.D0
        DFY(1,5)=0.D0
        DFY(1,6)=0.D0
        DFY(1,7)= RPAR(2)*RPAR(3)*Y(2)
        DFY(1,8)=0.D0
        DFY(1,9)=0.D0
        DFY(1,10)=-RPAR(4)*Y(1)
        DFY(2,1)= RPAR(8)*(RPAR(7)-Y(2)-Y(3))
        DFY(2,2)=-RPAR(8)*Y(1) - RPAR(3)*Y(7) - RPAR(9)
        DFY(2,3)=-RPAR(8)*Y(1) 
        DFY(2,4)=0.D0
        DFY(2,5)=0.D0
        DFY(2,6)=0.D0
        DFY(2,7)=-RPAR(3)*Y(2)
        DFY(2,8)=0.D0
        DFY(2,9)=0.D0
        DFY(2,10)=0.D0
        DFY(3,1)=0.D0
        DFY(3,2)= RPAR(3)*Y(7) + RPAR(9) 
        DFY(3,3)=-RPAR(10)
        DFY(3,4)=0.D0
        DFY(3,5)=0.D0
        DFY(3,6)=0.D0
        DFY(3,7)= RPAR(3)*Y(2)
        DFY(3,8)=0.D0
        DFY(3,9)=0.D0
        DFY(3,10)=0.D0
        DFY(4,1)= RPAR(11)*RPAR(12)
        DFY(4,2)=0.D0
        DFY(4,3)=0.D0
        DFY(4,4)=-RPAR(13)
        DFY(4,5)=0.D0
        DFY(4,6)=0.D0
        DFY(4,7)=0.D0
        DFY(4,8)=0.D0
        DFY(4,9)=0.D0
        DFY(4,10)=0.D0
        DFY(5,1)=0.D0
        DFY(5,2)=0.D0
        DFY(5,3)=-RPAR(14)*RPAR(15)*Y4L1*Y5L1/RPAR(7)
        DFY(5,4)=-RPAR(14)*Y(5) - RPAR(16)*Y(5)*Y(7)
        DFY(5,5)=-RPAR(14)*Y(4) - RPAR(16)*Y(4)*Y(7) - RPAR(17)
        DFY(5,6)=0.D0
        DFY(5,7)=-RPAR(16)*Y(4)*Y(5)
        DFY(5,8)=0.D0
        DFY(5,9)=0.D0
        DFY(5,10)=0.D0
        DFY(6,1)=0.D0
        DFY(6,2)=0.D0
        DFY(6,3)=-RPAR(19)*RPAR(20)*Y4L2*Y6L2/RPAR(7)
        DFY(6,4)=-RPAR(19)*Y(6) - RPAR(21)*Y(6)*Y(8)
        DFY(6,5)=0.D0
        DFY(6,6)=-RPAR(19)*Y(4) - RPAR(21)*Y(4)*Y(8) - RPAR(22)
        DFY(6,7)=0.D0
        DFY(6,8)=-RPAR(21)*Y(4)*Y(6)
        DFY(6,9)=0.D0
        DFY(6,10)=0.D0
        DFY(7,1)=0.D0
        DFY(7,2)=-RPAR(26)*Y(7)
        DFY(7,3)=-RPAR(24)*RPAR(25)*Y4L3*Y5L3*Y7L3/RPAR(7)
        DFY(7,4)=-RPAR(24)*Y(5)*Y(7)
        DFY(7,5)=-RPAR(24)*Y(4)*Y(7) 
        DFY(7,6)=0.D0
        DFY(7,7)=-RPAR(24)*Y(4)*Y(5) - RPAR(26)*Y(2) - RPAR(27)
        DFY(7,8)=0.D0
        DFY(7,9)=0.D0
        DFY(7,10)=0.D0
        DFY(8,1)=0.D0
        DFY(8,2)=0.D0
        DFY(8,3)=-RPAR(29)*Y4L4*Y6L4*Y8L4/RPAR(7)
        DFY(8,4)=-RPAR(29)*Y(6)*Y(8)
        DFY(8,5)=0.D0
        DFY(8,6)=-RPAR(29)*Y(4)*Y(8)
        DFY(8,7)=0.D0
        DFY(8,8)=-RPAR(29)*Y(4)*Y(6) - RPAR(31)
        DFY(8,9)=0.D0
        DFY(8,10)=0.D0
        DFY(9,1)=0.D0
        DFY(9,2)=0.D0
        DFY(9,3)=-RPAR(33)*RPAR(34)*Y4L5*Y6L5*Y8L5/RPAR(7)
        DFY(9,4)=0.D0
        DFY(9,5)=0.D0
        DFY(9,6)=0.D0
        DFY(9,7)=0.D0
        DFY(9,8)=0.D0
        DFY(9,9)=-RPAR(35)
        DFY(9,10)=0.D0
        DFY(10,1)=-RPAR(38)*Y(10)
        DFY(10,2)=0.D0
        DFY(10,3)=0.D0
        DFY(10,4)=0.D0
        DFY(10,5)=0.D0
        DFY(10,6)=0.D0
        DFY(10,7)=0.D0
        DFY(10,8)=0.D0
        DFY(10,9)= RPAR(37)
        DFY(10,10)=-RPAR(38)*Y(1) - RPAR(39)
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
        EXTERNAL PHI,ARGLAG
        XI(Z)=1.0D0-Z/RPAR(7)

C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,THETA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       SECOND DELAY
        CALL LAGR5(2,X,Y,ARGLAG,PAST,THETA2,IPOS2,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       THIRD DELAY
        CALL LAGR5(3,X,Y,ARGLAG,PAST,THETA3,IPOS3,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       FOURTH DELAY
        CALL LAGR5(4,X,Y,ARGLAG,PAST,THETA4,IPOS4,RPAR,IPAR,
     &             PHI,IPAST,NRDS)
C       FIFTH DELAY
        CALL LAGR5(5,X,Y,ARGLAG,PAST,THETA5,IPOS5,RPAR,IPAR,
     &             PHI,IPAST,NRDS)

          IVL(1)=1
          IVE(1)=5        
          IVC(1)=4
          IVL(2)=1
          IVE(2)=5        
          IVC(2)=5
          Y4L1=YLAGR5(4,THETA1,IPOS1,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y5L1=YLAGR5(5,THETA1,IPOS1,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          DFYL(1)=RPAR(14)*RPAR(15)*XI(Y(3))*Y5L1
          DFYL(2)=RPAR(14)*RPAR(15)*XI(Y(3))*Y4L1

          IVL(3)=2
          IVE(3)=6        
          IVC(3)=4
          IVL(4)=2
          IVE(4)=6        
          IVC(4)=6
          Y4L2=YLAGR5(4,THETA2,IPOS2,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y6L2=YLAGR5(6,THETA2,IPOS2,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          DFYL(3)=RPAR(19)*RPAR(20)*XI(Y(3))*Y6L2
          DFYL(4)=RPAR(19)*RPAR(20)*XI(Y(3))*Y4L2

          IVL(5)=3
          IVE(5)=7        
          IVC(5)=4
          IVL(6)=3
          IVE(6)=7        
          IVC(6)=5
          IVL(7)=3
          IVE(7)=7
          IVC(7)=7
          Y4L3=YLAGR5(4,THETA3,IPOS3,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y5L3=YLAGR5(5,THETA3,IPOS3,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y7L3=YLAGR5(7,THETA3,IPOS3,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          DFYL(5)=RPAR(24)*RPAR(25)*XI(Y(3))*Y5L3*Y7L3
          DFYL(6)=RPAR(24)*RPAR(25)*XI(Y(3))*Y4L3*Y7L3
          DFYL(7)=RPAR(24)*RPAR(25)*XI(Y(3))*Y4L3*Y5L3

          IVL(8)=4
          IVE(8)=8
          IVC(8)=4
          IVL(9)=4
          IVE(9)=8
          IVC(9)=6
          IVL(10)=4
          IVE(10)=8
          IVC(10)=8
          Y4L4=YLAGR5(4,THETA4,IPOS4,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y6L4=YLAGR5(6,THETA4,IPOS4,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y8L4=YLAGR5(8,THETA4,IPOS4,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          DFYL(8)=RPAR(29)*RPAR(30)*XI(Y(3))*Y6L4*Y8L4
          DFYL(9)=RPAR(29)*RPAR(30)*XI(Y(3))*Y4L4*Y8L4
          DFYL(10)=RPAR(29)*RPAR(30)*XI(Y(3))*Y4L4*Y6L4

          IVL(11)=5
          IVE(11)=9        
          IVC(11)=4
          IVL(12)=5
          IVE(12)=9        
          IVC(12)=6
          IVL(13)=5
          IVE(13)=9
          IVC(13)=8
          Y4L5=YLAGR5(4,THETA5,IPOS5,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y6L5=YLAGR5(6,THETA5,IPOS5,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          Y8L5=YLAGR5(8,THETA5,IPOS5,PHI,RPAR,IPAR,
     &                PAST,IPAST,NRDS)
          DFYL(11)=RPAR(33)*RPAR(34)*XI(Y(3))*Y6L5*Y8L5
          DFYL(12)=RPAR(33)*RPAR(34)*XI(Y(3))*Y4L5*Y8L5
          DFYL(13)=RPAR(33)*RPAR(34)*XI(Y(3))*Y4L5*Y6L5

        RETURN
        END
C
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR

        SELECT CASE (I)
        CASE (1)
            PHI=2.9D-16
        CASE (2)
            PHI=0.0D0
        CASE (3)
            PHI=0.0D0
        CASE (4)
            PHI=0.0D0
        CASE (5)
            PHI=RPAR(18)
        CASE (6)
            PHI=RPAR(23)
        CASE (7)
            PHI=RPAR(28)
        CASE (8)
            PHI=RPAR(32)
        CASE (9)
            PHI=0.0D0
        CASE (10)
            PHI=0.0D0
        END SELECT
        RETURN
        END

