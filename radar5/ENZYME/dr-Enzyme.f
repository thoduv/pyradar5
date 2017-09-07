C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=4
        INTEGER, PARAMETER :: NRDENS=1
        INTEGER, PARAMETER :: NGRID=11
        INTEGER, PARAMETER :: NLAGS=1
        INTEGER, PARAMETER :: NJACL=2
        INTEGER, PARAMETER :: MXST=1000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        INTEGER, PARAMETER :: NPAR=2
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+1) :: IPAST
        REAL(kind=DP), dimension(NPAR) :: RPAR
        INTEGER, dimension(22) :: ISTAT
        EXTERNAL  FCN,PHI,ARGLAG,JFCN,JACLAG,SOLOUT
  
C ------ FILE TO OPEN ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10

C --- PARAMETERS IN THE DIFFERENTIAL EQUATION
C ---   Tau
        RPAR(1)= 4.0D0
C ---   K1
        RPAR(2)= 0.0005D0
C -----------------------------------------------------------------------
C
C ---   DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXLPICIT FORM
        IMAS=0
        MLMAS=N

C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES 
        X=0.0D0
        Y(1)=6.0D1
        Y(2)=1.0D1
        Y(3)=1.0D1
        Y(4)=2.0D1
C       Consistent with initial function
C --- ENDPOINT OF INTEGRATION
        XEND=160.0D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.D-6
        ATOL=RTOL*1.0D0
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.0D0  

C --- WORKSPACE FOR PAST 
        IWORK(12)=MXST
C --- THE FOURTH COMPONENT USES RETARDED ARGUMENTS
        IWORK(15)=NRDENS
        IPAST(1)=4
C ---  SET THE PRESCRIBED GRID-POINTS
        DO I=1,NGRID
          GRID(I)=RPAR(1)*I
        END DO
C --- WORKSPACE FOR GRID
        IWORK(13)=NGRID
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
        WRITE (6,90) X,Y(1),Y(2),Y(3),Y(4)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=IWORK(J)
         END DO
        WRITE(6,*)' ***** TOL=',RTOL,' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
 90     FORMAT(1X,'X =',F8.2,'    Y =',4E18.10)
 91     FORMAT(' fcn=',I7,' jac=',I6,' step=',I6,
     &        ' accpt=',I6,' rejct=',I6,' dec=',I6,
     &        ' sol=',I7)
        WRITE(6,*) 'SOLUTION IS TABULATED IN FILES: sol.out & cont.out'

        STOP
        END


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

        WRITE (9,99) X,Y(1),Y(2),Y(3),Y(4)
C
        IF (NR.EQ.1) THEN
           WRITE (10,99) X,Y(1),Y(2),Y(3),Y(4)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (10,99) XOUT,CONTR5(1,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(2,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(3,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(4,N,XOUT,CONT,X,HSOL)
              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F12.8,'    Y =',4E18.10)
        RETURN
        END
C
        FUNCTION ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: Y
        REAL(kind=DP), dimension(1) :: PAST
	  INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(1) :: RPAR

        ARGLAG=X-RPAR(1)

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

C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,ALPHA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)

        Y4L1=YLAGR5(4,ALPHA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        U=1.D0/(1.D0 + RPAR(2)*Y4L1**3)

        F(1)= 10.5D0 - Y(1)*U
        F(2)= Y(1)*U - Y(2)
        F(3)= Y(2) - Y(3)
        F(4)= Y(3) - 0.5D0*Y(4)
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

C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,ALPHA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAST,NRDS)

        Y4L1=YLAGR5(4,ALPHA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        U=1.D0/(1.D0 + RPAR(2)*Y4L1**3)

C       STANDARD JACOBIAN MATRIX 4x4
        DFY(1,1)=-U
        DFY(1,2)= 0.D0
        DFY(1,3)= 0.D0
        DFY(1,4)= 0.D0
        DFY(2,1)= U
        DFY(2,2)=-1.D0
        DFY(2,3)= 0.D0
        DFY(2,4)= 0.D0
        DFY(3,1)= 0.D0
        DFY(3,2)= 1.D0
        DFY(3,3)=-1.D0
        DFY(3,4)= 0.D0
        DFY(4,1)= 0.D0
        DFY(4,2)= 0.D0
        DFY(4,3)= 1.D0
        DFY(4,4)=-0.5D0
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
        CALL LAGR5(1,X,Y,ARGLAG,PAST,ALPHA1,IPOS1,RPAR,IPAR,
     &             PHI,IPAT,NRDS)

        Y4L1=YLAGR5(4,ALPHA1,IPOS1,PHI,RPAR,IPAR,
     &              PAST,IPAST,NRDS)

        UD=(3.D0*RPAR(2)*Y(1)*Y4L1**2)/((1.D0+RPAR(2)*Y4L1**3)**2)

        IVL(1)=1
        IVE(1)=1
        IVC(1)=4
        IVL(2)=1
        IVE(2)=2
        IVC(2)=4
        DFYL(1)= UD
        DFYL(2)=-UD

        RETURN
        END
C
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(1) :: RPAR

        SELECT CASE (I)
        CASE (1)
            PHI=6.0D1
        CASE (2) 
            PHI=1.0D1
        CASE (3) 
            PHI=1.0D1
        CASE (4)
            PHI=2.0D1
        END SELECT
        RETURN
        END

