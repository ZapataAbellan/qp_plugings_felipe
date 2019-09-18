
       SUBROUTINE gauleg(x1,x2,x,w,n)
       IMPLICIT NONE
       INTEGER n
       REAL(8):: x1,x2,x(n),w(n)
       real(8):: EPS
       PARAMETER (EPS=3.d-14)
       INTEGER i,j,m
       real(8)::p1,p2,p3,pp,xl,xm,z,z1,pi
       pi = 4.d0*datan(1.d0)
       m=(n+1)/2
       xm=0.5d0*(x2+x1)
       xl=0.5d0*(x2-x1)
       do 12 i=1,m
          z=cos(pi*(i-.25d0)/(n+.5d0))
1         continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
             p3=p2
             p2=p1
             p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11           continue
             pp=n*(z*p1-p2)/(z*z-1.d0)
             z1=z
             z=z1-p1/pp
             if(abs(z-z1).gt.EPS)goto 1
             x(i)=xm-xl*z
             x(n+1-i)=xm+xl*z
             w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
             w(n+1-i)=w(i)
12           continue
             return
       END SUBROUTINE gauleg


         SUBROUTINE INTERV(XT,LXT,X,LEFT,MFLAG)
      !COMPUTES LEFT= MAX(1:XT(I) .LT. XT(LXT) .AND. XT(I) .LE. X)

      !INPUTS:
      !XT.......A REAL SEQUENCE OF LENGTH LXT, ASSUMED TO BE NONDECREASING
      !LXT......NUMBER OF TERMS IN THE SEQUENCE XT
      !X........THE POINT WHOSE LOCATION WITH RESPECT TO THE SEQUENCE XT IS
      !         TO BE DETERMINED.
         
      !OUTPUTS:
      !LEFT, MFLAG....BOTH INTEGERS, WHOSE VALUE IS:
      !
      !1    -1    IF             X .LT. XT(1)
      !I     0    IF  XT(I) .LE. X .LT. XT(I+1)
      !I     0    IF  XT(I) .LT. X .EQ. XT(I+1) .EQ. XT(LXT)
      !I     0    IF  XT(I) .LT.        XT(I+1) .EQ. XT(LXT) .LT. X
         
      !IN PARTICULAR, MFLAG=0 IS THE USUAL CASE. MFLAG .NE. 0 INDICATES THAT
      !X LIES OUTSIDE THE CLOSED INTERVAL XT(1).LE.Y.LE.XT(LXT).
      !THE ASYMMETRIC TREATEMENT OF THE INTERVALES IS DUE TO THE DECISION TO MAKE ALL 
      !PP-FUNCTIONS CONTINIUS FROM THE RIGTH, BUT, BY RETURNING MFLAG = 0 EVEN IF 
      !X=XT(LXT),THERE IS THE OPTION OF HAVING THE COMPUTED PP-FUNCTION CONTINUOUS FROM 
      !THE KEFT AT XT(LXT).
       IMPLICIT NONE
       INTEGER LEFT,LXT,MFLAG,IHI,ILO,ISTEP,MIDDLE
       REAL(8):: X,XT(LXT)
       DATA ILO /1/
       SAVE ILO

       IHI = ILO + 1

       IF (IHI.LT.LXT)        GOTO 20
       IF (X.GE.XT(LXT))      GOTO 110
       IF (LXT.LE.1)          GOTO 90
       ILO = LXT - 1
       IHI = LXT

20     IF (X.GE.XT(IHI))      GOTO 40
       IF (X.GE.XT(ILO))      GOTO 100

       ISTEP = 1
31     IHI = ILO
       ILO = IHI - ISTEP
       IF (ILO.LE.1)          GOTO 35
       IF (X.GE.XT(ILO))      GOTO 50
       ISTEP = ISTEP*2

35     ILO = 1                
       IF (X.LT.XT(1))        GOTO 90
                              GOTO 50

40     ISTEP = 1
41     ILO = IHI 
       IHI = ILO + ISTEP      
       IF (IHI.GE.LXT)        GOTO 45
       IF (X.LT.XT(IHI))      GOTO 50
       ISTEP = ISTEP*2
                              GOTO 41
45     IF (X.GE.XT(LXT))      GOTO 110
       IHI = LXT

50     MIDDLE = (ILO + IHI)/2
       IF (MIDDLE.EQ.ILO)           GOTO 100
      !IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
       IF(X.LT.XT(MIDDLE))          GOTO 53
       ILO = MIDDLE
                                    GOTO 50

53     IHI = MIDDLE
                                    GOTO 50
       
90     MFLAG = -1
       LEFT = 1 
                                    RETURN

100    MFLAG = 0 
       LEFT = ILO
                                    RETURN

110    MFLAG = 1 
       LEFT = LXT                               
       IF (X.EQ.XT(LXT)) MFLAG = 0
       LEFT = LXT                                         
111    IF(LEFT.EQ.1)                RETURN
       LEFT = LEFT - 1
       IF (XT(LEFT).LT.XT(LXT))     RETURN 
                                    GOTO 111

       END SUBROUTINE INTERV


      SUBROUTINE BSPLVD(T,NT,K,X,LEFT,A,DBIATX,NDERIV)

      !CALCULATEAS VALUE AND DERIVATIVES OF ALL B-SPLINES WHICH DO NOT VANISH
      !AT X

      !INPUT.................................................................
      !T : THE KNOT SEQUENCE
      !NT: THE NUMBER OF KNOT POINTS
      !K : THE ORDER OF THE B-SPLINES TO BE EVALUATED
      !X : THE POINT TO BE EVALUATED
      !LEFT : AN INTEGER INDICATING THE LEFT ENDPOINT OF THE INTERVAL
      !OF INTEREST. THE K B-SPLINES WHOSE SUPPORT CONTAINS THE INTERVAL.
      !NDERIV : AN INTEGER INDICATING THAT VALUES OF B-SPLINES AND THEIR
      !DERIVATIVES UP TO BUT NOT INCLUDING THE NDERIVATIVES ARE ASKED FOR.
      
      !WORK AREA.............................................................
      ! A : AN ARRAY OF ORDER(K,K) TO CONTAINT B-COEF.S OF THE DERIVATIVES OF A 
      ! CERTAIN ORDER OF THE K B-SPLINES OF INTEREST.
      
      !OUTPUT................................................................
      !DBIATX : AN ARRAY OF ORDER (K,NDERIV). ITS ENTRY CONTAINS VALUE OF 
      !(M-1)ST DERIVATIVE OF (LEFT-K+I)-TH B-SPLINE OF ORDER K FOR KNOT 
      !SEQUENCE T, I=1,...,K, M=1,....,NDERIV. 

       IMPLICIT NONE 

       INTEGER::NT,K,LEFT,NDERIV,I,IDERIV,IL,J,JLOW,JP1MID,KP1,KP1MM
       INTEGER::LDUMMY,M,MHIGH
       REAL(8):: A(K,K),DBIATX(K,NDERIV),T(NT),X,FACTOR,FKP1MM,SUM

       
       MHIGH = MAX(MIN(NDERIV,K),1)
      
       KP1 = K+1

       CALL BSPLVB(T,NT,KP1-MHIGH,1,X,LEFT,DBIATX)

       IF (MHIGH.EQ.1)                       GOTO 99

       IDERIV = MHIGH 

       DO 15 M=2,MHIGH
          JP1MID = 1
          DO 11 J=IDERIV,K
             DBIATX(J,IDERIV) = DBIATX(JP1MID,1)
11           JP1MID = JP1MID + 1
          IDERIV = IDERIV - 1
          CALL BSPLVB(T,NT,KP1-IDERIV,2,X,LEFT,DBIATX)
15        CONTINUE

       JLOW = 1
       DO 20 I=1,K
          DO 19 J=JLOW,K
19            A(J,I) = 0.D0
          JLOW = I
20        A(I,I) = 1.D0
          
       DO 40 M=2,MHIGH
          KP1MM = KP1 - M
          FKP1MM = FLOAT(KP1MM)
          IL = LEFT
          I = K

          DO 25 LDUMMY=1,KP1MM
             FACTOR = FKP1MM/(T(IL+KP1MM)-T(IL))
             DO 24 J=1,I
24              A(I,J)=(A(I,J)-A(I-1,J))*FACTOR
             IL = IL - 1
25           I = I - 1

         DO 40 I=1,K
            SUM = 0.D0
            JLOW = MAX(I,M)
            DO 35 J=JLOW,K
35             SUM = A(J,I)*DBIATX(J,M) + SUM
40          DBIATX(I,M) = SUM
99                                        RETURN
            
       END SUBROUTINE BSPLVD

       
       SUBROUTINE BSPLVB(T,NT,JHIGH,INDEX,X,LEFT,BIATX)
       IMPLICIT NONE 
      !THIS SUBROUTINE CALCULATES THE VALUE OF ALL POSSIBLE NONZERO 
      !B-SPLINES AT X OF ORDER 
      !
      !          JOUT = MAX(JHIGH,(J+1)*(INDEX-1))
      !
      !WITH KNOT SEQUENCE T.
      !
      !INPUTS:
      !T..........KNOT SEQUENCE, OF LENGTH LEFT+JOUT, ASSUMED TO BE NONDECREASING
      !
      !NT.........NUMBER OF KNOT POINTS
      !
      !INDEX......INTEGERS WHICH DETERMINE THE ORDER JOUT=(JHIGH,(J+1)*(INDEX-1))
      !           OF B-SPLINE VALUE ARE NEEDED.  
      !           THE RESTRICTION JOUT.LE.JMAX (=20) IS IMPOSED ARBITRARIALY BY THE 
      !           DIMENSION STATEMENT FOR DELTAL AND DELTAR, BUT IS NOWHERE CHECKED FOR.
      !
      !X..........THE POINT AT WHICH THE B-SPLINES ARE TO BE EVALUATED.
      !
      !LEFT.......AN INTEGER CHOSEN (USUALLY) SO THAT: T(LEFT).LE.X.LE.T(LEFT+1)
      

      !OUTPUT:
      !BIATX......ARRAY OF LENGTH JOUT, WITH BIATX(I) CONTAINING THE VALUE AT X OF THE 
      !           POLYNOMIAL OF ORDER JOUT WHICH AGREES WITH THE B-SPLINE
      !           B(LEFT-JOUT+I,JOUT,T) ON THE INTERVAL (T(LEFT)...T(LEFT+1))
      !
      !

      !METHOD: THE RECURRENCE RELATION....
      !
      !                 X - T(I)                T(I+J+1) - X
      ! B(I,J+1)(X) = ------------ B(I,J)(X) + --------------- B(I+1,J)(X)
      !                T(I+J)-T(I)             T(I+J+1)-T(I+1)
      !

       INTEGER INDEX,JHIGH,LEFT,I,J,JMAX,JP1,NT
       PARAMETER (JMAX=20)
     !  REAL(8)::BIATX(JHIGH),T(LEFT+JHIGH),X,DELTAL(JMAX),DELTAR(JMAX),SAVED,TERM
       REAL(8)::BIATX(JHIGH),T(NT),X,DELTAL(JMAX),DELTAR(JMAX),SAVED,TERM
     ! DIMENSION BIATX(JOUT),T(LEFT+JOUT)
       DATA J/1/
       SAVE J,DELTAL,DELTAR

                          GOTO (10,20), INDEX
10     J = 1
       BIATX(1) = 1.d0
       IF (J .GE. JHIGH)  GOTO 99

20     JP1 = J + 1

       DELTAR(J) = T(LEFT+J) - X
       DELTAL(J) = X - T(LEFT+1-J)

       SAVED = 0.

       DO 26 I=1,J
          TERM = BIATX(I)/(DELTAR(I) + DELTAL(JP1-I))
          BIATX(I) = SAVED + DELTAR(I)*TERM
26        SAVED = DELTAL(JP1-I)*TERM

       BIATX(JP1)= SAVED
       J = JP1
       IF (J .LT. JHIGH)  GOTO 20
       
99     RETURN

       END SUBROUTINE BSPLVB
